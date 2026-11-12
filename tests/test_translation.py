"""
Regression tests for Dish-to-Brain Human Post-Mortem Translational Concordance.

Checks cell-type marker mapping, post-mortem disease concordance against
NYGC ALS and SEA-AD snRNA-seq benchmarks, developmental maturation gap detection,
and discordant direction penalty scoring.

Run: python3 -m tests.test_translation
"""

from __future__ import annotations

import sys

from synaptik.schema import Direction
from synaptik.translate.atlas import (
    BrainCellType, DiseaseAtlas, HumanTranslationalConcordance,
    PostMortemStatus, evaluate_translational_concordance,
)

FAILURES: list[str] = []


def check(cond: bool, msg: str) -> None:
    print(f"  {'PASS' if cond else 'FAIL'}  {msg}")
    if not cond:
        FAILURES.append(msg)


def test_als_concordance() -> None:
    print("\n[translation: ALS post-mortem concordance]")
    entries = [
        ("STMN2", Direction.INCREASE, 0.42),   # Cryptic exon inclusion increase
        ("UNC13A", Direction.INCREASE, 0.31),  # Cryptic donor activation
        ("KALRN", Direction.DECREASE, -0.28),  # Downregulation
        ("TARDBP", Direction.DECREASE, -0.50), # Nuclear depletion
    ]

    res = evaluate_translational_concordance(entries, disease="als")

    check(res.n_genes_evaluated == 4, f"evaluated 4 genes (got {res.n_genes_evaluated})")
    check(res.n_concordant == 4, f"all 4 hallmark ALS genes are concordant (got {res.n_concordant})")
    check(res.n_discordant == 0, "zero discordant genes on true positive set")
    check(res.concordance_score_pct == 100.0, f"concordance score is 100% (got {res.concordance_score_pct}%)")
    check("HIGH TRANSLATIONAL FIDELITY" in res.translational_verdict, "high fidelity verdict rendered")

    calls = {c.gene: c for c in res.gene_calls}
    check("Betz cell" in calls["STMN2"].primary_cell_type or "Motor Neuron" in calls["STMN2"].primary_cell_type,
          "STMN2 mapped to human motor neurons / Betz cells")
    check(calls["STMN2"].is_cryptic_exon, "STMN2 recognized as cryptic exon pathology")
    check("Brown" in calls["UNC13A"].evidence_source or "Nature" in calls["UNC13A"].evidence_source,
          "UNC13A evidence cites landmark Nature study")


def test_developmental_maturation_gap() -> None:
    print("\n[translation: developmental maturation gap]")
    entries = [
        ("NEFH", Direction.NO_CHANGE, 0.0),    # Adult neurofilament heavy
        ("GRIN2A", Direction.NO_CHANGE, 0.0),  # Adult NMDA subunit
        ("TUBB3", Direction.NO_CHANGE, 0.0),   # Early pan-neuronal
    ]

    res = evaluate_translational_concordance(entries, disease="als")
    calls = {c.gene: c for c in res.gene_calls}

    check(calls["NEFH"].status == PostMortemStatus.DEVELOPMENTAL_GAP,
          "NEFH recognized as adult developmental maturation gap in early iPSCs")
    check(calls["GRIN2A"].status == PostMortemStatus.DEVELOPMENTAL_GAP,
          "GRIN2A recognized as adult developmental maturation gap")
    check(res.n_developmental_gap == 2, f"counted 2 developmental gaps (got {res.n_developmental_gap})")
    check(any("DIV 45-60" in r for r in res.recommendations),
          "prolonged culture recommendation provided for developmental gaps")


def test_direction_discordance() -> None:
    print("\n[translation: direction discordance]")
    # Claiming TARDBP is increased in ALS motor cortex directly conflicts with post-mortem pathology
    entries = [
        ("TARDBP", Direction.INCREASE, 0.60),
        ("KALRN", Direction.INCREASE, 0.40),
    ]

    res = evaluate_translational_concordance(entries, disease="als")
    check(res.n_discordant == 2, f"opposing directions flagged as discordant (got {res.n_discordant})")
    check("MIXED TRANSLATIONAL FIDELITY" in res.translational_verdict or "contradiction" in res.translational_verdict.lower(),
          "discordant finding triggers caution verdict")
    check(any("Address discordant" in r for r in res.recommendations),
          "recommendation to investigate in vitro culture confounds")


def test_cell_type_and_other_diseases() -> None:
    print("\n[translation: cell-type specificity and other diseases]")
    # FTD test
    ftd_entries = [("GRN", Direction.DECREASE, -0.4), ("MAPT", Direction.INCREASE, 0.3)]
    ftd_res = evaluate_translational_concordance(ftd_entries, disease="ftd")
    check(ftd_res.n_concordant >= 1, "GRN down-regulation concordant in FTD")

    # AD test
    ad_entries = [("TREM2", Direction.INCREASE, 0.5), ("APOE", Direction.INCREASE, 0.4)]
    ad_res = evaluate_translational_concordance(ad_entries, disease="ad")
    check(ad_res.n_concordant >= 2, "TREM2 and APOE concordant in AD post-mortem atlas")

    # Glial cell type mapping
    entries = [("GFAP", Direction.INCREASE, 0.5), ("AIF1", Direction.INCREASE, 0.5)]
    res = evaluate_translational_concordance(entries, disease="als")
    calls = {c.gene: c for c in res.gene_calls}
    check(calls["GFAP"].primary_cell_type == BrainCellType.ASTROCYTE.value, "GFAP mapped to Astrocyte")
    check(calls["AIF1"].primary_cell_type == BrainCellType.MICROGLIA.value, "AIF1 mapped to Microglia")


def main() -> int:
    print("=" * 70)
    print("SYNAPTIK post-mortem translational concordance test suite")
    print("=" * 70)
    test_als_concordance()
    test_developmental_maturation_gap()
    test_direction_discordance()
    test_cell_type_and_other_diseases()
    print("\n" + "=" * 70)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S):")
        for f in FAILURES:
            print("  -", f)
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
