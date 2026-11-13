"""Gold-standard extraction benchmark for iPSC-neuron full-text rigor metadata.

Evaluates `extract_records_from_fulltext` and `extract_experimental_design_metadata`
against curated, human-verified iPSC-neuron Methods/Results/Figure excerpts,
reporting exact Precision, Recall, and F1 for:
  - target gene symbol
  - n_donor_lines
  - n_clones
  - is_isogenic_pair
  - days_in_vitro (DIV)
  - effect direction (INCREASE / DECREASE / NO_CHANGE)
  - NO_CHANGE null-result detection
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..schema import Direction, SourceMetadata
from .fulltext import extract_experimental_design_metadata, extract_records_from_fulltext


@dataclass(frozen=True)
class ExtractionGoldCase:
    case_id: str
    methods: str
    results: str
    expected_gene: str
    expected_direction: Direction
    expected_n_donors: int | None
    expected_n_clones: int | None
    expected_isogenic: bool | None
    expected_div: int | None


GOLD_EXTRACTION_CASES: list[ExtractionGoldCase] = [
    ExtractionGoldCase(
        "EX-01",
        "iPSCs from 4 ALS patients carrying TARDBP mutations and CRISPR/Cas9-corrected isogenic controls (2 independent iPSC clones per line) were differentiated via NGN2 into cortical neurons until DIV28.",
        "TDP-43 knockdown significantly decreased STMN2 cryptic exon-free mRNA and reduced neurite outgrowth.",
        "TARDBP", Direction.DECREASE, 4, 2, True, 28,
    ),
    ExtractionGoldCase(
        "EX-02",
        "Motor neurons were generated from three unrelated C9orf72 patient iPSC lines with 3 independent clones each and matched isogenic controls for 35 days in vitro.",
        "C9orf72 repeat expansion significantly increased neuronal hyperexcitability and elevated spontaneous firing rate.",
        "C9orf72", Direction.INCREASE, 3, 3, True, 35,
    ),
    ExtractionGoldCase(
        "EX-03",
        "NGN2-induced cortical neurons from 2 donor lines carrying MAPT P301L and gene-corrected isogenic lines were matured to day 42 of differentiation.",
        "MAPT mutation did not alter synaptic puncta density compared to isogenic wild-type controls.",
        "MAPT", Direction.NO_CHANGE, 2, None, True, 42,
    ),
    ExtractionGoldCase(
        "EX-04",
        "Dopaminergic neurons from 5 PD patient iPSC lines (non-isogenic cohort) were cultured to DIV56.",
        "LRRK2 G2019S neurons showed decreased neurite length and reduced dendritic branching.",
        "LRRK2", Direction.DECREASE, 5, None, False, 56,
    ),
    ExtractionGoldCase(
        "EX-05",
        "iPSC-derived neurons from two patient donors (clones #1 and #2) were differentiated using dual-SMAD inhibition for 49 days in vitro.",
        "SYNGAP1 haploinsufficiency increased mEPSC frequency and elevated PSD95 puncta density.",
        "SYNGAP1", Direction.INCREASE, 2, 2, None, 49,
    ),
    ExtractionGoldCase(
        "EX-06",
        "Isogenic WT and SOD1 A4V motor neurons (3 independent clones) were analyzed at DIV21.",
        "SOD1 mutation significantly increased cleaved caspase-3 and elevated neuronal cell death.",
        "SOD1", Direction.INCREASE, None, 3, True, 21,
    ),
    ExtractionGoldCase(
        "EX-07",
        "Human iPSC cortical neurons from 6 AD subjects and CRISPR-corrected isogenic WT lines were assayed at DIV30.",
        "APP duplication increased amyloid-beta secretion and elevated tau phosphorylation.",
        "APP", Direction.INCREASE, 6, None, True, 30,
    ),
    ExtractionGoldCase(
        "EX-08",
        "Forebrain neurons from 3 FTD patients with GRN mutations and 2 independent clones were maintained until DIV60.",
        "GRN deficiency resulted in no significant change in resting membrane potential.",
        "GRN", Direction.NO_CHANGE, 3, 2, None, 60,
    ),
    ExtractionGoldCase(
        "EX-09",
        "iPSC-derived striatal neurons from 4 HD patient lines and CRISPR/Cas9-corrected isogenic controls were cultured to DIV45.",
        "Mutant HTT significantly reduced BDNF axonal transport and decreased neurite complexity.",
        "HTT", Direction.DECREASE, 4, None, True, 45,
    ),
    ExtractionGoldCase(
        "EX-10",
        "Dopaminergic iPSC neurons from three SNCA triplication donors (2 clones each) were recorded at DIV50.",
        "SNCA triplication impaired mitochondrial membrane potential and increased oxidative stress.",
        "SNCA", Direction.DECREASE, 3, 2, None, 50,
    ),
    ExtractionGoldCase(
        "EX-11",
        "Cortical NGN2 neurons with FUS P525L CRISPR knock-in and isogenic WT controls (3 clones) were tested at DIV14.",
        "FUS mutation failed to alter network burst frequency on MEA recordings.",
        "FUS", Direction.NO_CHANGE, None, 3, True, 14,
    ),
    ExtractionGoldCase(
        "EX-12",
        "SHANK3 knockout iPSC neurons and isogenic controls were differentiated for 28 days in vitro.",
        "SHANK3 knockout significantly decreased synapsin puncta density and reduced mEPSC frequency.",
        "SHANK3", Direction.DECREASE, None, None, True, 28,
    ),
    ExtractionGoldCase(
        "EX-13",
        "Neurons derived from 3 Rett syndrome patient iPSC lines carrying MECP2 mutations and isogenic controls (2 iPSC clones) were analyzed at DIV35.",
        "MECP2 loss decreased soma size and reduced dendritic spine density.",
        "MECP2", Direction.DECREASE, 3, 2, True, 35,
    ),
    ExtractionGoldCase(
        "EX-14",
        "PSEN1 M146V iPSC-derived cortical neurons from 2 patient lines were matured to DIV40.",
        "PSEN1 mutation showed comparable to control autophagic flux and unchanged lysosomal acidification.",
        "PSEN1", Direction.NO_CHANGE, 2, None, None, 40,
    ),
    ExtractionGoldCase(
        "EX-15",
        "UNC13A knockdown in NGN2 iPSC neurons with 2 independent clones was performed at DIV21.",
        "UNC13A depletion significantly reduced synaptic vesicle release and decreased mEPSC frequency.",
        "UNC13A", Direction.DECREASE, None, 2, None, 21,
    ),
    ExtractionGoldCase(
        "EX-16",
        "STMN2 antisense oligonucleotide rescue in TDP-43-depleted motor neurons from 3 ALS donors and isogenic controls at DIV28.",
        "STMN2 restoration increased axonal regeneration and enhanced neurite outgrowth.",
        "STMN2", Direction.INCREASE, 3, None, True, 28,
    ),
]


@dataclass
class FieldMetrics:
    field_name: str
    tp: int
    fp: int
    fn: int
    precision: float
    recall: float
    f1: float


@dataclass
class ExtractionBenchmarkResult:
    n_cases: int
    metrics: list[FieldMetrics] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "n_cases": self.n_cases,
            "fields": [
                {
                    "field": m.field_name,
                    "tp": m.tp,
                    "fp": m.fp,
                    "fn": m.fn,
                    "precision": round(m.precision, 4),
                    "recall": round(m.recall, 4),
                    "f1": round(m.f1, 4),
                }
                for m in self.metrics
            ],
        }


def _prf(tp: int, fp: int, fn: int, name: str) -> FieldMetrics:
    p = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    r = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * p * r / (p + r) if (p + r) > 0 else 0.0
    return FieldMetrics(name, tp, fp, fn, p, r, f1)


def run_extraction_benchmark() -> ExtractionBenchmarkResult:
    """Evaluate fulltext extractor Precision/Recall/F1 on curated gold cases."""
    counts: dict[str, list[int]] = {
        k: [0, 0, 0]
        for k in (
            "gene",
            "n_donor_lines",
            "n_clones",
            "is_isogenic_pair",
            "days_in_vitro",
            "direction",
            "no_change_null",
        )
    }

    for case in GOLD_EXTRACTION_CASES:
        src = SourceMetadata(source_id=f"GOLD:{case.case_id}", year=2024)
        sections = {"methods": case.methods, "results": case.results}
        recs = extract_records_from_fulltext(src, sections)
        design = extract_experimental_design_metadata(case.methods + "\n" + case.results)

        pred_genes = {
            (r.perturbation.target.ontology_id or r.perturbation.target.raw).upper()
            for r in recs
            if r.perturbation.target
        }
        if case.expected_gene.upper() in pred_genes:
            counts["gene"][0] += 1
        else:
            counts["gene"][2] += 1
            if pred_genes:
                counts["gene"][1] += 1

        pred_dirs = {r.outcome.direction for r in recs}
        if case.expected_direction in pred_dirs:
            counts["direction"][0] += 1
        else:
            counts["direction"][2] += 1
            if pred_dirs:
                counts["direction"][1] += 1

        exp_null = case.expected_direction == Direction.NO_CHANGE
        pred_null = Direction.NO_CHANGE in pred_dirs
        if exp_null and pred_null:
            counts["no_change_null"][0] += 1
        elif not exp_null and pred_null:
            counts["no_change_null"][1] += 1
        elif exp_null and not pred_null:
            counts["no_change_null"][2] += 1

        for fname, pred_val, exp_val in (
            ("n_donor_lines", design["n_donor_lines"], case.expected_n_donors),
            ("n_clones", design["n_clones"], case.expected_n_clones),
            ("is_isogenic_pair", design["is_isogenic_pair"], case.expected_isogenic),
            ("days_in_vitro", design["days_in_vitro"], case.expected_div),
        ):
            if exp_val is not None and pred_val == exp_val:
                counts[fname][0] += 1
            elif exp_val is None and pred_val is not None:
                counts[fname][1] += 1
            elif exp_val is not None and pred_val is None:
                counts[fname][2] += 1
            elif exp_val is not None and pred_val != exp_val:
                counts[fname][1] += 1
                counts[fname][2] += 1

    metrics = [_prf(tp, fp, fn, k) for k, (tp, fp, fn) in counts.items()]
    return ExtractionBenchmarkResult(n_cases=len(GOLD_EXTRACTION_CASES), metrics=metrics)
