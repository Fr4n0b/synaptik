"""
LLM-based structured extraction from full text.

Design decisions that matter for quality:

1. We extract from Methods + Results + figure legends jointly, because the
   experimental parameters (DIV, protocol, markers) live in Methods while the
   outcome lives in Results. Extracting from either alone produces records
   that cannot be compared across papers.

2. The model is asked to emit a strict JSON array and to quote the supporting
   span for every record. The quote is validated against the source text; a
   record whose quote does not appear in the source is dropped as a
   hallucination. This is the single most effective quality control in the
   pipeline.

3. Negative results are explicitly solicited. The prompt asks for "no change"
   findings by name, because they are the scarcest and most valuable evidence
   and models otherwise skip them as uninteresting.

4. The model is instructed never to infer a number that is not stated. A
   missing effect size is a null, not a guess.
"""

from __future__ import annotations

import json
import logging
import os
import re
from dataclasses import dataclass
from typing import Any, Callable

from ..schema import (
    AssayClass, CellModel, CellModelClass, Direction, EvidenceType,
    ExperimentRecord, NormalizedTerm, Outcome, Perturbation, PerturbationType,
    Provenance, SourceMetadata, SplicingDetail, SplicingEventType,
)

log = logging.getLogger(__name__)

EXTRACTOR_VERSION = "llm-v1.2"

SYSTEM_PROMPT = """You extract structured experimental records from neuroscience \
and stem-cell papers. You are precise and conservative.

An "experiment record" is one measured outcome: a specific biological system, \
subjected to a specific perturbation (or none), assayed by a specific method, \
producing a specific measured change.

RULES, in priority order:

1. NEVER invent a number. If the paper does not state an effect size, sample \
size, or p-value, use null. Do not estimate from a figure description.

2. EVERY record must include a `quote` field containing a VERBATIM span of \
15-200 characters copied exactly from the provided text that supports the \
outcome. If you cannot find a verbatim supporting span, do not emit the record.

3. REPORT NULL RESULTS. If the paper states that something was measured and \
did NOT change, emit a record with direction "no_change". These are as \
important as positive findings. Look for phrases like "no significant \
difference", "unchanged", "comparable to controls".

4. Distinguish what was MEASURED from what the authors CONCLUDE. Set \
evidence_type to "direct_measurement" only when the outcome was assayed in \
that experiment. Use "inferred" for the authors' interpretation, "rescue" when \
a treatment reverses a phenotype, "correlation" for observed associations.

5. Use the exact enum values given. If nothing fits, use the "other" value.

6. Do not emit records for background statements about other papers' findings.
Only extract experiments performed in THIS paper.

Output a JSON array. No prose, no markdown fences. If no extractable \
experiments exist, output [].
"""

RECORD_TEMPLATE = """{
  "cell_model_class": "one of: ipsc_neuron_ngn2, ipsc_neuron_small_molecule, ipsc_npc, ipsc_organoid, primary_rodent_neuron, immortalized_line, in_vivo_rodent, post_mortem_human, other",
  "cell_model_description": "verbatim-ish description of the system",
  "species": "human|mouse|rat|other|null",
  "differentiation_protocol": "e.g. 'doxycycline-inducible NGN2' or null",
  "days_in_vitro": "integer or null",
  "coculture": "e.g. 'mouse glia' or null",
  "genetic_background": "cell line ID if stated, else null",
  "is_isogenic_pair": "true/false/null",
  "n_donor_lines": "integer or null (distinct donor/patient iPSC lines per genotype)",
  "n_clones": "integer or null (distinct iPSC clones per genotype)",

  "perturbation_type": "one of: mutation_patient, mutation_engineered, knockout, knockdown, overexpression, small_molecule, antisense_oligonucleotide, isoform_switch, none, other",
  "perturbation_target": "gene/protein symbol, or null",
  "perturbation_agent": "compound or reagent name, or null",
  "variant_hgvs": "HGVS string if stated, else null",
  "dose": "string or null",
  "duration": "string or null",
  "control_description": "what it was compared to",

  "assay_description": "the method used, as described",
  "markers_used": ["antibody or reporter names"],

  "phenotype": "what was measured, e.g. 'SYN1 puncta density'",
  "direction": "one of: increase, decrease, no_change, biphasic, unspecified",
  "effect_size": "number or null",
  "effect_size_unit": "percent_change|fold|cohens_d|absolute|null",
  "p_value": "number or null",
  "n_biological_replicates": "integer or null",
  "n_cells": "integer or null",
  "statistical_test": "string or null",

  "evidence_type": "one of: direct_measurement, rescue, correlation, inferred, review_assertion",
  "pathways_mentioned": ["pathway names the authors link to this result"],
  "genes_mentioned": ["gene symbols relevant to this result"],

  "splicing": {
    "event_type": "SE|A5SS|A3SS|MXE|RI|AFE|ALE|APA|unspecified",
    "gene": "gene symbol or null",
    "delta_psi": "number between -1 and 1, or null",
    "detection_method": "rMATS|MAJIQ|RT-PCR|other|null",
    "regulator": "RBP implicated, or null"
  } or null,

  "section": "which section this came from",
  "quote": "VERBATIM span from the text, 15-200 chars",
  "confidence": "0.0-1.0, your confidence this record is accurate"
}"""


@dataclass
class ExtractionResult:
    records: list[ExperimentRecord]
    n_proposed: int
    n_quote_failed: int
    n_schema_failed: int
    raw_response: str = ""

    @property
    def hallucination_rate(self) -> float:
        return self.n_quote_failed / self.n_proposed if self.n_proposed else 0.0


class LLMExtractor:
    """Calls an LLM to extract records, then validates every one.

    `call_fn` is injected so this works with any provider or a local model.
    It takes (system_prompt, user_prompt) and returns the model's text.
    """

    def __init__(
        self,
        call_fn: Callable[[str, str], str] | None = None,
        model: str = "claude-sonnet-4-6",
        max_chars_per_chunk: int = 24000,
        min_quote_len: int = 15,
    ):
        self.call_fn = call_fn or self._default_anthropic_call
        self.model = model
        self.max_chars = max_chars_per_chunk
        self.min_quote_len = min_quote_len

    # ---------------- provider ----------------

    def _default_anthropic_call(self, system: str, user: str) -> str:
        """Default provider. Requires ANTHROPIC_API_KEY in the environment."""
        try:
            import anthropic
        except ImportError as e:
            raise RuntimeError(
                "anthropic package not installed. Run: pip install anthropic\n"
                "Or pass your own call_fn to LLMExtractor."
            ) from e
        key = os.environ.get("ANTHROPIC_API_KEY")
        if not key:
            raise RuntimeError("ANTHROPIC_API_KEY not set in environment.")
        client = anthropic.Anthropic(api_key=key)
        resp = client.messages.create(
            model=self.model,
            max_tokens=8000,
            system=system,
            messages=[{"role": "user", "content": user}],
        )
        return "".join(
            b.text for b in resp.content if getattr(b, "type", "") == "text"
        )

    # ---------------- prompt assembly ----------------

    def build_prompt(self, meta: SourceMetadata, sections: dict[str, str]) -> str:
        parts = [
            f"PAPER: {meta.title}",
            f"JOURNAL: {meta.journal or 'unknown'} ({meta.year or 'n.d.'})",
            "",
        ]
        for name in ("abstract", "methods", "results", "figures"):
            txt = sections.get(name, "")
            if txt:
                parts.append(f"--- {name.upper()} ---")
                parts.append(txt)
                parts.append("")

        body = "\n".join(parts)
        if len(body) > self.max_chars:
            body = self._truncate_preserving_sections(sections, meta)

        return (
            body
            + "\n\nExtract all experiment records from the text above. "
            + "Each record must follow this exact structure:\n\n"
            + RECORD_TEMPLATE
            + "\n\nOutput ONLY a JSON array of such objects."
        )

    def _truncate_preserving_sections(
        self, sections: dict[str, str], meta: SourceMetadata
    ) -> str:
        """Budget characters across sections rather than truncating the tail.

        Naive truncation drops the Results section entirely on long papers,
        which is exactly the content we need. We allocate a fixed share to each
        section so every one survives.
        """
        budget = {
            "abstract": 0.10, "methods": 0.30, "results": 0.45, "figures": 0.15,
        }
        parts = [
            f"PAPER: {meta.title}",
            f"JOURNAL: {meta.journal or 'unknown'} ({meta.year or 'n.d.'})",
            "",
        ]
        for name, frac in budget.items():
            txt = sections.get(name, "")
            if not txt:
                continue
            allowed = int(self.max_chars * frac)
            if len(txt) > allowed:
                txt = txt[:allowed] + " [...truncated...]"
            parts.append(f"--- {name.upper()} ---")
            parts.append(txt)
            parts.append("")
        return "\n".join(parts)

    # ---------------- extraction ----------------

    def extract(
        self, meta: SourceMetadata, sections: dict[str, str]
    ) -> ExtractionResult:
        source_text = "\n".join(sections.values())
        prompt = self.build_prompt(meta, sections)

        try:
            raw = self.call_fn(SYSTEM_PROMPT, prompt)
        except Exception as e:
            log.error("LLM call failed for %s: %s", meta.source_id, e)
            return ExtractionResult([], 0, 0, 0, "")

        proposals = self._parse_json_array(raw)
        records, n_quote_fail, n_schema_fail = [], 0, 0

        for p in proposals:
            if not self._validate_quote(p.get("quote"), source_text):
                n_quote_fail += 1
                continue
            try:
                rec = self._to_record(p, meta)
                records.append(rec)
            except Exception as e:
                log.debug("Schema build failed for %s: %s", meta.source_id, e)
                n_schema_fail += 1

        return ExtractionResult(
            records, len(proposals), n_quote_fail, n_schema_fail, raw
        )

    # ---------------- validation ----------------

    @staticmethod
    def _parse_json_array(raw: str) -> list[dict]:
        """Tolerant JSON extraction. Models sometimes add fences or prose."""
        if not raw:
            return []
        txt = raw.strip()
        txt = re.sub(r"^```(?:json)?\s*", "", txt)
        txt = re.sub(r"\s*```$", "", txt)
        try:
            data = json.loads(txt)
            return data if isinstance(data, list) else [data]
        except json.JSONDecodeError:
            pass
        start, end = txt.find("["), txt.rfind("]")
        if start != -1 and end > start:
            try:
                data = json.loads(txt[start : end + 1])
                return data if isinstance(data, list) else [data]
            except json.JSONDecodeError:
                pass
        log.warning("Could not parse LLM output as JSON (%d chars)", len(raw))
        return []

    def _validate_quote(self, quote: str | None, source_text: str) -> bool:
        """Verify the quote actually appears in the source. Anti-hallucination.

        Exact match is too brittle (whitespace, unicode), so we normalize both
        sides. We still require a genuine substring match -- fuzzy matching here
        would defeat the purpose of the check.
        """
        if not quote or len(quote) < self.min_quote_len:
            return False
        norm = lambda s: re.sub(r"\s+", " ", s.lower()).strip()
        return norm(quote) in norm(source_text)

    # ---------------- record construction ----------------

    def _to_record(self, p: dict, meta: SourceMetadata) -> ExperimentRecord:
        from ..normalize.matcher import classify_assay

        assay_raw = p.get("assay_description") or ""
        assay_cls, _ = classify_assay(assay_raw)
        if assay_cls == AssayClass.OTHER:
            assay_cls = self._enum(AssayClass, p.get("assay_class"), AssayClass.OTHER)

        cm = CellModel(
            model_class=self._enum(
                CellModelClass, p.get("cell_model_class"), CellModelClass.OTHER
            ),
            raw_description=p.get("cell_model_description") or "",
            species=p.get("species"),
            differentiation_protocol=p.get("differentiation_protocol"),
            days_in_vitro=self._int(p.get("days_in_vitro")),
            coculture=p.get("coculture"),
            genetic_background=p.get("genetic_background"),
            is_isogenic_pair=self._bool(p.get("is_isogenic_pair")),
            n_donor_lines=self._int(p.get("n_donor_lines")),
            n_clones=self._int(p.get("n_clones")),
        )

        target = p.get("perturbation_target")
        pert = Perturbation(
            ptype=self._enum(
                PerturbationType, p.get("perturbation_type"), PerturbationType.OTHER
            ),
            target=NormalizedTerm(raw=target) if target else None,
            agent=p.get("perturbation_agent"),
            variant_hgvs=p.get("variant_hgvs"),
            dose=p.get("dose"),
            duration=p.get("duration"),
            control_description=p.get("control_description"),
        )

        markers = p.get("markers_used") or []
        if isinstance(markers, str):
            markers = [markers]

        outcome = Outcome(
            phenotype=NormalizedTerm(raw=p.get("phenotype") or "unspecified"),
            direction=self._enum(Direction, p.get("direction"), Direction.UNSPECIFIED),
            effect_size=self._float(p.get("effect_size")),
            effect_size_unit=p.get("effect_size_unit"),
            p_value=self._float(p.get("p_value")),
            n_biological_replicates=self._int(p.get("n_biological_replicates")),
            n_cells=self._int(p.get("n_cells")),
            statistical_test=p.get("statistical_test"),
            markers_used=[str(m) for m in markers],
            is_quantified=self._float(p.get("effect_size")) is not None,
        )

        splicing = None
        sp = p.get("splicing")
        if isinstance(sp, dict) and any(sp.get(k) for k in
                                        ("event_type", "gene", "delta_psi")):
            gene = sp.get("gene")
            reg = sp.get("regulator")
            splicing = SplicingDetail(
                event_type=self._enum(
                    SplicingEventType, sp.get("event_type"),
                    SplicingEventType.UNSPECIFIED
                ),
                gene=NormalizedTerm(raw=gene) if gene else None,
                delta_psi=self._clamp(self._float(sp.get("delta_psi")), -1.0, 1.0),
                detection_method=sp.get("detection_method"),
                regulator=NormalizedTerm(raw=reg) if reg else None,
            )

        pathways = [NormalizedTerm(raw=x) for x in (p.get("pathways_mentioned") or [])
                    if isinstance(x, str) and x.strip()]
        genes = [NormalizedTerm(raw=x) for x in (p.get("genes_mentioned") or [])
                 if isinstance(x, str) and x.strip()]

        return ExperimentRecord(
            source=meta,
            cell_model=cm,
            perturbation=pert,
            outcome=outcome,
            assay_class=assay_cls,
            assay_raw=assay_raw,
            evidence_type=self._enum(
                EvidenceType, p.get("evidence_type"), EvidenceType.DIRECT_MEASUREMENT
            ),
            pathways=pathways,
            genes=genes,
            splicing=splicing,
            provenance=Provenance(
                section=p.get("section"),
                quote=p.get("quote"),
                extractor_version=EXTRACTOR_VERSION,
                extraction_confidence=self._clamp(
                    self._float(p.get("confidence")) or 0.5, 0.0, 1.0
                ),
            ),
        )

    # ---------------- coercion helpers ----------------

    @staticmethod
    def _enum(enum_cls, value, default):
        if value is None:
            return default
        v = str(value).strip().lower().replace(" ", "_").replace("-", "_")
        for m in enum_cls:
            if m.value.lower() == v or m.name.lower() == v:
                return m
        return default

    @staticmethod
    def _int(v) -> int | None:
        if v is None or isinstance(v, bool):
            return None
        try:
            return int(float(str(v).strip().rstrip("+")))
        except (ValueError, TypeError):
            m = re.search(r"-?\d+", str(v))
            return int(m.group()) if m else None

    @staticmethod
    def _float(v) -> float | None:
        if v is None or isinstance(v, bool):
            return None
        s = str(v).strip()
        if s.lower() in ("null", "none", "na", "n/a", ""):
            return None
        # Handle "<0.001", "~34", "34%"
        s = s.lstrip("<>~≈").rstrip("%")
        try:
            return float(s)
        except ValueError:
            m = re.search(r"-?\d*\.?\d+(?:[eE][-+]?\d+)?", s)
            return float(m.group()) if m else None

    @staticmethod
    def _bool(v) -> bool | None:
        if isinstance(v, bool):
            return v
        if v is None:
            return None
        s = str(v).strip().lower()
        if s in ("true", "yes", "1"):
            return True
        if s in ("false", "no", "0"):
            return False
        return None

    @staticmethod
    def _clamp(v: float | None, lo: float, hi: float) -> float | None:
        return None if v is None else max(lo, min(hi, v))
