"""Synthetic smoke tests for iPSC-neuron full-text rigor extraction.

IMPORTANT: These 16 cases are SYNTHETIC test sentences with zero PMIDs/DOIs,
written to verify regex/parser behavior. They are NOT a real-paper gold
standard and must not be cited as an empirical F1 benchmark. Real-paper
extraction recall/precision requires Phase 1 independent human annotation.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass

from synaptik.extract.fulltext import (
    extract_experimental_design_metadata,
    extract_records_from_fulltext,
)
from synaptik.schema import Direction, SourceMetadata


@dataclass(frozen=True)
class SyntheticExtractionCase:
    case_id: str
    methods: str
    results: str
    expected_gene: str
    expected_direction: Direction
    expected_n_donors: int | None
    expected_n_clones: int | None
    expected_isogenic: bool | None
    expected_div: int | None


SYNTHETIC_EXTRACTION_CASES: list[SyntheticExtractionCase] = [
    SyntheticExtractionCase(
        "SYNTH-01",
        "iPSCs from 4 ALS patients carrying TARDBP mutations and CRISPR/Cas9-corrected isogenic controls (2 independent iPSC clones per line) were differentiated via NGN2 into cortical neurons until DIV28.",
        "TDP-43 knockdown significantly decreased STMN2 cryptic exon-free mRNA and reduced neurite outgrowth.",
        "TARDBP", Direction.DECREASE, 4, 2, True, 28,
    ),
    SyntheticExtractionCase(
        "SYNTH-02",
        "Motor neurons were generated from three unrelated C9orf72 patient iPSC lines with 3 independent clones each and matched isogenic controls for 35 days in vitro.",
        "C9orf72 repeat expansion significantly increased neuronal hyperexcitability and elevated spontaneous firing rate.",
        "C9orf72", Direction.INCREASE, 3, 3, True, 35,
    ),
    SyntheticExtractionCase(
        "SYNTH-03",
        "NGN2-induced cortical neurons from 2 donor lines carrying MAPT P301L and gene-corrected isogenic lines were matured to day 42 of differentiation.",
        "MAPT mutation did not alter synaptic puncta density compared to isogenic wild-type controls.",
        "MAPT", Direction.NO_CHANGE, 2, None, True, 42,
    ),
    SyntheticExtractionCase(
        "SYNTH-04",
        "Dopaminergic neurons from 5 PD patient iPSC lines (non-isogenic cohort) were cultured to DIV56.",
        "LRRK2 G2019S neurons showed decreased neurite length and reduced dendritic branching.",
        "LRRK2", Direction.DECREASE, 5, None, False, 56,
    ),
    SyntheticExtractionCase(
        "SYNTH-05",
        "iPSC-derived neurons from two patient donors (clones #1 and #2) were differentiated using dual-SMAD inhibition for 49 days in vitro.",
        "SYNGAP1 haploinsufficiency increased mEPSC frequency and elevated PSD95 puncta density.",
        "SYNGAP1", Direction.INCREASE, 2, 2, None, 49,
    ),
    SyntheticExtractionCase(
        "SYNTH-06",
        "Isogenic WT and SOD1 A4V motor neurons (3 independent clones) were analyzed at DIV21.",
        "SOD1 mutation significantly increased cleaved caspase-3 and elevated neuronal cell death.",
        "SOD1", Direction.INCREASE, None, 3, True, 21,
    ),
    SyntheticExtractionCase(
        "SYNTH-07",
        "Forebrain neurons from 3 FTD patients with GRN mutations and 2 independent clones were maintained until DIV60.",
        "GRN deficiency resulted in no significant change in resting membrane potential.",
        "GRN", Direction.NO_CHANGE, 3, 2, None, 60,
    ),
    SyntheticExtractionCase(
        "SYNTH-08",
        "Dopaminergic iPSC neurons from three SNCA triplication donors (2 clones each) were recorded at DIV50.",
        "SNCA triplication impaired mitochondrial membrane potential and increased oxidative stress.",
        "SNCA", Direction.DECREASE, 3, 2, None, 50,
    ),
    SyntheticExtractionCase(
        "SYNTH-09",
        "Cortical NGN2 neurons with FUS P525L CRISPR knock-in and isogenic WT controls (3 clones) were tested at DIV14.",
        "FUS mutation failed to alter network burst frequency on MEA recordings.",
        "FUS", Direction.NO_CHANGE, None, 3, True, 14,
    ),
    SyntheticExtractionCase(
        "SYNTH-10",
        "SHANK3 knockout iPSC neurons and isogenic controls were differentiated for 28 days in vitro.",
        "SHANK3 knockout significantly decreased synapsin puncta density and reduced mEPSC frequency.",
        "SHANK3", Direction.DECREASE, None, None, True, 28,
    ),
    SyntheticExtractionCase(
        "SYNTH-11",
        "Neurons derived from 3 Rett syndrome patient iPSC lines carrying MECP2 mutations and isogenic controls (2 iPSC clones) were analyzed at DIV35.",
        "MECP2 loss decreased soma size and reduced dendritic spine density.",
        "MECP2", Direction.DECREASE, 3, 2, True, 35,
    ),
    SyntheticExtractionCase(
        "SYNTH-12",
        "PSEN1 M146V iPSC-derived cortical neurons from 2 patient lines were matured to DIV40.",
        "PSEN1 mutation showed comparable to control autophagic flux and unchanged lysosomal acidification.",
        "PSEN1", Direction.NO_CHANGE, 2, None, None, 40,
    ),
    SyntheticExtractionCase(
        "SYNTH-13",
        "UNC13A knockdown in NGN2 iPSC neurons with 2 independent clones was performed at DIV21.",
        "UNC13A depletion significantly reduced synaptic vesicle release and decreased mEPSC frequency.",
        "UNC13A", Direction.DECREASE, None, 2, None, 21,
    ),
    SyntheticExtractionCase(
        "SYNTH-14",
        "STMN2 antisense oligonucleotide rescue in TDP-43-depleted motor neurons from 3 ALS donors and isogenic controls at DIV28.",
        "STMN2 restoration increased axonal regeneration and enhanced neurite outgrowth.",
        "STMN2", Direction.INCREASE, 3, None, True, 28,
    ),
]


def main() -> int:
    failures: list[str] = []
    for case in SYNTHETIC_EXTRACTION_CASES:
        src = SourceMetadata(source_id=f"SYNTH:{case.case_id}", year=2024)
        recs = extract_records_from_fulltext(src, {"methods": case.methods, "results": case.results})
        design = extract_experimental_design_metadata(case.methods + "\n" + case.results)

        pred_genes = {
            (r.perturbation.target.ontology_id or r.perturbation.target.raw).upper()
            for r in recs
            if r.perturbation.target
        }
        if case.expected_gene.upper() not in pred_genes:
            failures.append(f"{case.case_id}: expected gene {case.expected_gene}, got {pred_genes}")
        pred_dirs = {r.outcome.direction for r in recs}
        if case.expected_direction not in pred_dirs:
            failures.append(f"{case.case_id}: expected direction {case.expected_direction}, got {pred_dirs}")
        for fname, got, exp in (
            ("n_donor_lines", design["n_donor_lines"], case.expected_n_donors),
            ("n_clones", design["n_clones"], case.expected_n_clones),
            ("is_isogenic_pair", design["is_isogenic_pair"], case.expected_isogenic),
            ("days_in_vitro", design["days_in_vitro"], case.expected_div),
        ):
            if got != exp:
                failures.append(f"{case.case_id}: {fname} expected {exp}, got {got}")

    if failures:
        for f in failures:
            print("FAIL:", f)
        return 1
    print(f"All {len(SYNTHETIC_EXTRACTION_CASES)} synthetic extraction smoke cases passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
