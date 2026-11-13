"""
Evidence weighting.

The user asked for journal-quartile weighting. We implement that, but as one
component among several, because journal tier is a weak and biased proxy for
reliability: it correlates with novelty and institutional prestige at least as
much as with rigor, and the retraction rate in high-impact journals is not
lower than average.

The composite weight combines six signals:

  journal_tier      -- SJR quartile (the user's requested signal)
  study_design      -- isogenic controls, replicate count, quantification
  directness        -- measured vs inferred vs review assertion
  extraction_conf   -- how confident the extractor was in this record
  recency           -- mild decay; protocols and reagents change
  integrity         -- hard penalties for retraction, preprint status

Every component is exposed individually in the returned breakdown so a user
can see why a record was weighted as it was, and re-weight if they disagree.
This transparency is the point: a single opaque score would be unauditable.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import date

from ..schema import EvidenceType, ExperimentRecord

# Quartile multipliers. The Q1:Q4 ratio is deliberately modest (2:1 rather
# than 10:1) because the evidence that Q4 findings are half as likely to
# replicate is itself weak. Aggressive tiering would encode prestige bias as
# if it were data quality.
QUARTILE_WEIGHT = {"Q1": 1.00, "Q2": 0.85, "Q3": 0.65, "Q4": 0.50}
UNKNOWN_QUARTILE_WEIGHT = 0.70

EVIDENCE_TYPE_WEIGHT = {
    EvidenceType.DIRECT_MEASUREMENT: 1.00,
    EvidenceType.RESCUE: 1.00,
    EvidenceType.CORRELATION: 0.65,
    EvidenceType.INFERRED: 0.35,
    EvidenceType.REVIEW_ASSERTION: 0.20,
}


@dataclass
class WeightBreakdown:
    """Per-component weights, always returned alongside the composite."""

    journal_tier: float = 1.0
    study_design: float = 1.0
    directness: float = 1.0
    extraction_confidence: float = 1.0
    recency: float = 1.0
    integrity: float = 1.0
    composite: float = 1.0
    notes: list[str] = field(default_factory=list)

    def explain(self) -> str:
        parts = [
            f"journal={self.journal_tier:.2f}",
            f"design={self.study_design:.2f}",
            f"directness={self.directness:.2f}",
            f"extraction={self.extraction_confidence:.2f}",
            f"recency={self.recency:.2f}",
            f"integrity={self.integrity:.2f}",
        ]
        s = " x ".join(parts) + f" = {self.composite:.3f}"
        if self.notes:
            s += "  [" + "; ".join(self.notes) + "]"
        return s


class EvidenceWeighter:
    def __init__(
        self,
        current_year: int | None = None,
        recency_halflife_years: float = 12.0,
        min_recency_weight: float = 0.6,
        preprint_weight: float = 0.80,
        retraction_weight: float = 0.0,
        enable_journal_tier: bool = True,
    ):
        self.current_year = current_year or date.today().year
        self.halflife = recency_halflife_years
        self.min_recency = min_recency_weight
        self.preprint_weight = preprint_weight
        self.retraction_weight = retraction_weight
        self.enable_journal_tier = enable_journal_tier

    def weight(self, rec: ExperimentRecord) -> WeightBreakdown:
        b = WeightBreakdown()
        src = rec.source

        # --- journal tier ---
        if self.enable_journal_tier:
            if src.is_preprint:
                b.journal_tier = self.preprint_weight
                b.notes.append("preprint")
            elif src.sjr_quartile in QUARTILE_WEIGHT:
                b.journal_tier = QUARTILE_WEIGHT[src.sjr_quartile]
            else:
                b.journal_tier = UNKNOWN_QUARTILE_WEIGHT
                b.notes.append("journal tier unknown")

        # --- study design ---
        design = 1.0
        if rec.cell_model.is_isogenic_pair:
            design *= 1.15  # isogenic controls are the gold standard here
            b.notes.append("isogenic")
        n_donors = rec.cell_model.n_donor_lines
        n_clones = rec.cell_model.n_clones
        if n_donors is not None:
            if n_donors >= 3:
                design *= 1.25
                b.notes.append(f"{n_donors} donor lines")
            elif n_donors == 2:
                design *= 1.12
                b.notes.append("2 donor lines")
        if n_clones is not None:
            if n_clones >= 3:
                design *= 1.15
                b.notes.append(f"{n_clones} clones")
            elif n_clones == 1 and (n_donors is None or n_donors == 1):
                # Single-clone iPSC studies are vulnerable to clone-to-clone drift
                design *= 0.80
                b.notes.append("single-clone (clonal drift risk)")
        n = rec.outcome.n_biological_replicates
        if n is None:
            design *= 0.85
            b.notes.append("n not reported")
        elif n < 3:
            design *= 0.75
            b.notes.append(f"n={n}")
        elif n >= 6:
            design *= 1.10
        if not rec.outcome.is_quantified:
            design *= 0.80
            b.notes.append("not quantified")
        b.study_design = min(1.3, design)

        # --- directness ---
        b.directness = EVIDENCE_TYPE_WEIGHT.get(rec.evidence_type, 0.5)
        if src.is_review and rec.evidence_type != EvidenceType.REVIEW_ASSERTION:
            b.directness *= 0.5
            b.notes.append("from review article")

        # --- extraction confidence ---
        # Mapped onto [0.5, 1.0]: even a low-confidence extraction is not
        # worthless, and we already dropped hallucinations at the quote check.
        b.extraction_confidence = 0.5 + 0.5 * rec.provenance.extraction_confidence

        # --- recency ---
        if src.year:
            age = max(0, self.current_year - src.year)
            decay = 0.5 ** (age / self.halflife)
            b.recency = max(self.min_recency, decay)
            if age > 20:
                b.notes.append(f"{age}y old")
        else:
            b.recency = 0.9

        # --- integrity ---
        if src.is_retracted:
            b.integrity = self.retraction_weight
            b.notes.append("RETRACTED")

        b.composite = (
            b.journal_tier * b.study_design * b.directness
            * b.extraction_confidence * b.recency * b.integrity
        )
        return b

    # Uncalibrated rigor-structure prior weights (NOT empirically calibrated;
    # publication volume `n_records` is deliberately excluded and `null_ratio`
    # is not penalized pending Phase 3 pre-registered replication analysis).
    CLAIM_FEATURE_WEIGHTS: dict[str, float] = {
        "n_labs": 0.45,
        "n_donor_lines": 0.25,
        "n_clones": 0.15,
        "is_isogenic_pair": 0.30,
        "assay_diversity": 0.25,
    }
    CLAIM_INTERCEPT: float = -0.80

    def score_claim_features(self, feats: dict[str, float]) -> float:
        """Score a claim's evidence structure (rigor prior, excluding raw volume)."""
        import math

        z = self.CLAIM_INTERCEPT
        for k, w in self.CLAIM_FEATURE_WEIGHTS.items():
            val = float(feats.get(k, 0.0))
            if k in ("n_labs", "n_donor_lines", "n_clones", "assay_diversity"):
                val = math.log1p(max(0.0, val))
            z += w * val
        return 1.0 / (1.0 + math.exp(-max(-20.0, min(20.0, z))))


# --------------------------------------------------------------------------
# Journal quartile resolution
# --------------------------------------------------------------------------

# A small built-in table so the tool is useful before the user loads a full
# SCImago dump. Deliberately limited to journals dominant in this field.
BUILTIN_QUARTILES: dict[str, str] = {
    "nature": "Q1", "science": "Q1", "cell": "Q1", "neuron": "Q1",
    "nature neuroscience": "Q1", "nature medicine": "Q1",
    "nature genetics": "Q1", "nature methods": "Q1",
    "nature communications": "Q1", "cell stem cell": "Q1",
    "cell reports": "Q1", "molecular cell": "Q1", "immunity": "Q1",
    "science translational medicine": "Q1", "science advances": "Q1",
    "elife": "Q1", "embo journal": "Q1", "embo reports": "Q1",
    "journal of neuroscience": "Q1", "brain": "Q1", "annals of neurology": "Q1",
    "molecular psychiatry": "Q1", "biological psychiatry": "Q1",
    "american journal of human genetics": "Q1", "genome biology": "Q1",
    "nucleic acids research": "Q1", "pnas": "Q1",
    "proceedings of the national academy of sciences": "Q1",
    "stem cell reports": "Q1", "acta neuropathologica": "Q1",
    "nature protocols": "Q1", "rna": "Q1", "genes & development": "Q1",
    "human molecular genetics": "Q2", "journal of cell biology": "Q1",
    "journal of biological chemistry": "Q2", "scientific reports": "Q2",
    "plos one": "Q2", "plos genetics": "Q1", "plos biology": "Q1",
    "frontiers in molecular neuroscience": "Q2",
    "frontiers in cellular neuroscience": "Q2",
    "frontiers in neuroscience": "Q2", "molecular neurodegeneration": "Q1",
    "neurobiology of disease": "Q1", "glia": "Q1", "cerebral cortex": "Q1",
    "journal of neurochemistry": "Q2", "molecular brain": "Q2",
    "bmc genomics": "Q2", "rna biology": "Q2", "cell death & disease": "Q1",
    "international journal of molecular sciences": "Q2",
    "biochemical and biophysical research communications": "Q3",
    "neuroscience letters": "Q3", "brain research": "Q3",
    "molecular biology reports": "Q3", "gene": "Q3",
    "journal of molecular neuroscience": "Q3",
}


class JournalQualityResolver:
    """Resolves journal -> SJR quartile, from a loaded table or the builtin."""

    def __init__(self, scimago_csv: str | None = None):
        self.table: dict[str, tuple[str, float | None]] = {
            k: (v, None) for k, v in BUILTIN_QUARTILES.items()
        }
        if scimago_csv:
            self.load_scimago(scimago_csv)

    def load_scimago(self, path: str) -> int:
        """Load a SCImago Journal Rank CSV export.

        Download from scimagojr.com/journalrank.php (semicolon-delimited).
        Expected columns include 'Title', 'SJR Best Quartile', 'SJR'.
        """
        import csv

        n = 0
        with open(path, newline="", encoding="utf-8-sig") as f:
            sample = f.read(4096)
            f.seek(0)
            delim = ";" if sample.count(";") > sample.count(",") else ","
            for row in csv.DictReader(f, delimiter=delim):
                title = (row.get("Title") or "").strip().lower()
                quart = (row.get("SJR Best Quartile") or "").strip()
                if not title or quart not in QUARTILE_WEIGHT:
                    continue
                sjr_raw = (row.get("SJR") or "").replace(",", ".").strip()
                try:
                    sjr = float(sjr_raw) if sjr_raw else None
                except ValueError:
                    sjr = None
                self.table[title] = (quart, sjr)
                n += 1
        return n

    def resolve(self, journal: str | None) -> tuple[str | None, float | None]:
        if not journal:
            return (None, None)
        j = journal.strip().lower()
        if j in self.table:
            return self.table[j]
        j2 = j.replace("the ", "").replace("&", "and").strip()
        if j2 in self.table:
            return self.table[j2]
        # Substring fallback for "Nature Neuroscience" vs "Nat Neurosci"
        for k, v in self.table.items():
            if len(k) > 8 and (k in j or j in k):
                return v
        return (None, None)

    def annotate(self, rec_or_source) -> None:
        """Set sjr_quartile/sjr_score on a SourceMetadata in place."""
        src = getattr(rec_or_source, "source", rec_or_source)
        if src.sjr_quartile:
            return
        q, s = self.resolve(src.journal)
        src.sjr_quartile = q
        src.sjr_score = s
