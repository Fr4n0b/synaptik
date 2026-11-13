"""Full-text JATS XML experimental extractor for Europe PMC OA articles.

Extracts sentence-grounded ExperimentRecords from Methods, Results, Figures,
and Abstract sections, recovering experimental-rigor metadata that rarely
appears in abstracts:
  - n_donor_lines (distinct patient/donor iPSC lines)
  - n_clones (independent clones per line)
  - is_isogenic_pair (CRISPR-corrected isogenic controls)
  - days_in_vitro (DIV / days of neuronal maturation)
  - cell model class (NGN2, small-molecule/dual-SMAD, organoid, NPC, primary)
  - null / negative results ("no significant difference", "did not alter", etc.)
"""

from __future__ import annotations

import re

from ..normalize.matcher import PhenotypeMatcher, classify_assay
from ..schema import (
    CellModel,
    CellModelClass,
    Direction,
    EvidenceType,
    ExperimentRecord,
    NormalizedTerm,
    Outcome,
    Perturbation,
    PerturbationType,
    Provenance,
    SourceMetadata,
)
from .heuristic import (
    _DECREASE_RE,
    _GENE_RE,
    _INCREASE_RE,
    _NULL_RE,
    _PATHWAY_HINTS,
    _infer_model_class,
)

_WORD_TO_INT = {
    "one": 1, "two": 2, "three": 3, "four": 4, "five": 5,
    "six": 6, "seven": 7, "eight": 8, "nine": 9, "ten": 10,
}

_DONOR_PATTERNS = [
    re.compile(
        r"\b(\d+|two|three|four|five|six|seven|eight|nine|ten)\s+"
        r"(?:independent\s+|unrelated\s+|distinct\s+|different\s+)?"
        r"(?:patient|donor|als|ad|pd|hd|ftd|sma|mutant|diseased|affected|subject|individual)s?\s*"
        r"(?:ipsc\s*)?(?:cell\s*)?(?:lines?|donors?|patients?|subjects?)\b",
        re.I,
    ),
    re.compile(
        r"\bfrom\s+(\d+|two|three|four|five|six|seven|eight|nine|ten)\s+"
        r"(?:unrelated\s+|distinct\s+|independent\s+)?"
        r"(?:patients?|donors?|individuals?|subjects?)\b",
        re.I,
    ),
]

_CLONE_PATTERNS = [
    re.compile(
        r"\b(\d+|two|three|four|five|six)\s+"
        r"(?:independent\s+|distinct\s+)?(?:ipsc\s+)?clones?\b",
        re.I,
    ),
    re.compile(r"\bclone(?:s)?\s*(?:#|no\.?\s*)?(\d+)\s*(?:and|,)\s*(?:#|no\.?\s*)?(\d+)\b", re.I),
]

_ISOGENIC_POS_RE = re.compile(
    r"\b(?:isogenic\s+(?:controls?|corrected|lines?|pairs?|wt|wild[- ]type|revertants?)|"
    r"crispr[- /]cas9[- ]corrected\s+isogenic|gene[- ]corrected\s+isogenic)\b",
    re.I,
)
_ISOGENIC_NEG_RE = re.compile(r"\bnon[- ]isogenic\b", re.I)

_DIV_PATTERNS = [
    re.compile(r"\bDIV\s*[-:]?\s*(\d{1,3})\b", re.I),
    re.compile(r"\b(\d{1,3})\s*days?\s+in\s+vitro\b", re.I),
    re.compile(r"\bday\s+(\d{1,3})\s+(?:of\s+)?(?:differentiation|maturation|post[- ]induction)\b", re.I),
]

_NULL_EXTENDED_RE = re.compile(
    r"\b(no significant (?:change|difference|effect|alteration)|unchanged|"
    r"did not (?:alter|affect|change|rescue|significantly)|comparable to (?:control|wt|wild[- ]type)|"
    r"failed to (?:rescue|alter|affect)|no difference between)\b",
    re.I,
)

_ENGINEERED_RE = re.compile(r"\b(crispr|knock[- ]?in|isogenic|engineered)\b", re.I)
_KNOCKDOWN_RE = re.compile(r"\b(sirna|shrna|aso|antisense|knockdown|knocked down|deplet\w*|crispri)\b", re.I)
_KNOCKOUT_RE = re.compile(r"\b(knockout|knock[- ]?out|ko\b|null\s+allele|crispr[- ]ko)\b", re.I)
_OVEREXPR_RE = re.compile(r"\b(overexpress\w*|ectopic|transgen\w*|lentiviral\s+expression)\b", re.I)
_DRUG_RE = re.compile(r"\b(treated with|compound|inhibitor|agonist|antagonist|small molecule|\d+\s*[nuµ]m)\b", re.I)


def _parse_num(token: str) -> int | None:
    tok = token.strip().lower()
    if tok.isdigit():
        val = int(tok)
        return val if 1 <= val <= 100 else None
    return _WORD_TO_INT.get(tok)


def extract_experimental_design_metadata(text: str) -> dict:
    """Extract iPSC rigor metadata (donors, clones, isogenic, DIV, model) from text."""
    if not text:
        return {
            "n_donor_lines": None,
            "n_clones": None,
            "is_isogenic_pair": None,
            "days_in_vitro": None,
            "model_class": CellModelClass.OTHER,
        }

    donors: list[int] = []
    for pat in _DONOR_PATTERNS:
        for m in pat.finditer(text):
            v = _parse_num(m.group(1))
            if v is not None:
                donors.append(v)

    clones: list[int] = []
    for m in _CLONE_PATTERNS[0].finditer(text):
        v = _parse_num(m.group(1))
        if v is not None:
            clones.append(v)
    if not clones and _CLONE_PATTERNS[1].search(text):
        clones.append(2)

    isogenic: bool | None = None
    if _ISOGENIC_POS_RE.search(text):
        isogenic = True
    elif _ISOGENIC_NEG_RE.search(text):
        isogenic = False

    divs: list[int] = []
    for pat in _DIV_PATTERNS:
        for m in pat.finditer(text):
            val = int(m.group(1))
            if 1 <= val <= 365:
                divs.append(val)

    return {
        "n_donor_lines": max(donors) if donors else None,
        "n_clones": max(clones) if clones else None,
        "is_isogenic_pair": isogenic,
        "days_in_vitro": divs[0] if divs else None,
        "model_class": _infer_model_class(text),
    }


def _infer_perturbation_type(sent: str) -> PerturbationType:
    if _KNOCKOUT_RE.search(sent):
        return PerturbationType.KNOCKOUT
    if _KNOCKDOWN_RE.search(sent):
        return PerturbationType.KNOCKDOWN
    if _OVEREXPR_RE.search(sent):
        return PerturbationType.OVEREXPRESSION
    if _DRUG_RE.search(sent):
        return PerturbationType.SMALL_MOLECULE
    if _ENGINEERED_RE.search(sent):
        return PerturbationType.MUTATION_ENGINEERED
    return PerturbationType.MUTATION_PATIENT


def extract_records_from_fulltext(
    source: SourceMetadata,
    sections: dict[str, str],
    matcher: PhenotypeMatcher | None = None,
) -> list[ExperimentRecord]:
    """Extract ExperimentRecords from Europe PMC JATS full-text section dict."""
    if not sections:
        return []

    pm = matcher or PhenotypeMatcher()
    methods_txt = sections.get("methods", "")
    results_txt = sections.get("results", "")
    abstract_txt = sections.get("abstract", "")
    figures_txt = sections.get("figures", "")

    combined_design_txt = "\n".join(p for p in (methods_txt, abstract_txt, results_txt) if p)
    design = extract_experimental_design_metadata(combined_design_txt)

    # Prefer neuronal cell-model classification if any section specifies NGN2/iPSC
    model_class = design["model_class"]
    if model_class == CellModelClass.OTHER:
        model_class = _infer_model_class(combined_design_txt)

    doc_genes = _GENE_RE.findall(combined_design_txt)
    default_gene = (
        "TARDBP"
        if doc_genes and doc_genes[0] == "TDP-43"
        else (doc_genes[0] if doc_genes else None)
    )

    # Scan Results first (richest experimental claims), then Figures, then Abstract
    scan_sections = [
        ("results", results_txt),
        ("figures", figures_txt),
        ("abstract", abstract_txt),
    ]

    records: list[ExperimentRecord] = []
    seen_keys: set[tuple[str, str, str]] = set()

    for sec_name, sec_text in scan_sections:
        if not sec_text or not sec_text.strip():
            continue
        sentences = re.split(r"(?<=[.!?])\s+", sec_text.strip())
        for idx, sent in enumerate(sentences):
            if len(sent) < 25:
                continue

            if _NULL_EXTENDED_RE.search(sent) or _NULL_RE.search(sent):
                direction = Direction.NO_CHANGE
            elif _DECREASE_RE.search(sent):
                direction = Direction.DECREASE
            elif _INCREASE_RE.search(sent):
                direction = Direction.INCREASE
            else:
                continue

            pheno = pm.match(sent)
            if not pheno:
                continue

            sent_genes = _GENE_RE.findall(sent)
            gene = (
                "TARDBP"
                if sent_genes and sent_genes[0] == "TDP-43"
                else (sent_genes[0] if sent_genes else default_gene)
            )
            if not gene:
                continue

            pathway_term: NormalizedTerm | None = None
            for pat, label, go_id in _PATHWAY_HINTS:
                if pat.search(sent):
                    pathway_term = NormalizedTerm(
                        raw=label,
                        ontology_id=go_id,
                        ontology_label=label,
                        mapping_confidence=0.80,
                        mapping_method="fulltext_keyword",
                    )
                    break

            key = (gene, pheno.key(), direction.value)
            if key in seen_keys:
                continue
            seen_keys.add(key)

            # Check sentence-local DIV override if stated
            sent_div_m = _DIV_PATTERNS[0].search(sent) or _DIV_PATTERNS[1].search(sent)
            div = int(sent_div_m.group(1)) if sent_div_m else design["days_in_vitro"]

            assay = classify_assay(sent)
            ptype = _infer_perturbation_type(sent)
            rec_id = f"{source.source_id}:{sec_name}:{idx}"

            records.append(
                ExperimentRecord(
                    record_id=rec_id,
                    source=source,
                    cell_model=CellModel(
                        model_class=model_class,
                        raw_description=sent[:180],
                        days_in_vitro=div,
                        is_isogenic_pair=design["is_isogenic_pair"],
                        n_donor_lines=design["n_donor_lines"],
                        n_clones=design["n_clones"],
                    ),
                    perturbation=Perturbation(
                        ptype=ptype,
                        target=NormalizedTerm(
                            raw=gene,
                            ontology_label=gene,
                            mapping_confidence=0.90,
                            mapping_method="gene_symbol",
                        ),
                    ),
                    assay=assay,
                    outcome=Outcome(
                        phenotype=pheno,
                        pathway=pathway_term,
                        direction=direction,
                        evidence_type=EvidenceType.DIRECT_MEASUREMENT,
                    ),
                    provenance=Provenance(
                        section=sec_name,
                        quote=sent[:200],
                        extractor_version="fulltext-jats-v1",
                        extraction_confidence=0.80 if sec_name == "results" else 0.72,
                    ),
                )
            )

    return records
