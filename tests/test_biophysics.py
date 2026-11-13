"""
3D RNA Biophysics & Steric Accessibility Test Suite.

Run: python -m tests.test_biophysics
"""

from __future__ import annotations

import sys

from synaptik.biophysics import (
    AccessibilityTier,
    compute_delta_g_open,
    compute_effective_binding,
    compute_steric_displacement,
    compute_target_accessibility,
    fold_mfe_pure_python,
    get_hairpin_energy,
    get_stacking_energy,
    get_terminal_penalty,
    get_tetraloop_bonus,
)

FAILURES: list[str] = []


def check(cond: bool, msg: str) -> None:
    print(f"  {'PASS' if cond else 'FAIL'}  {msg}")
    if not cond:
        FAILURES.append(msg)


def test_parameters() -> None:
    print("\n[biophysics: Turner 2004 parameters]")
    check(abs(get_stacking_energy("GC", "CG") - (-3.42)) < 0.1, "GC/CG stacking is -3.42 kcal/mol")
    check(abs(get_stacking_energy("AA", "UU") - (-0.93)) < 0.1, "AA/UU stacking is -0.93 kcal/mol")
    check(abs(get_terminal_penalty("A", "U") - 0.45) < 0.05, "terminal AU penalty is +0.45 kcal/mol")
    check(get_terminal_penalty("G", "C") == 0.0, "terminal GC penalty is 0.0 kcal/mol")
    check(abs(get_hairpin_energy(4) - 5.60) < 0.1, "4-nt hairpin initiation is 5.60 kcal/mol")
    check(get_hairpin_energy(2) > 100.0, "2-nt hairpin is sterically forbidden")
    check(abs(get_tetraloop_bonus("GAAA") - (-3.00)) < 0.1, "GAAA tetraloop bonus is -3.00 kcal/mol")


def test_folding_and_accessibility() -> None:
    print("\n[biophysics: Zuker MFE folding & Delta G_open]")
    mfe, db = fold_mfe_pure_python("GCGCUUUUGCGC")
    check(mfe < -3.0 and db == "((((....))))", f"GCGCUUUUGCGC folds into stem-loop ({mfe} kcal/mol, {db})")

    mfe_a, db_a = fold_mfe_pure_python("AAAAAAAAAAAAAAAAAAAA")
    check(mfe_a == 0.0 and db_a == "." * 20, "poly-A remains completely unstructured")

    dg_loop = compute_delta_g_open("GCGCUUUUGCGC", 4, 4)
    check(dg_loop == 0.0, "single-stranded loop has Delta G_open = 0.0")

    dg_stem = compute_delta_g_open("GCGCUUUUGCGC", 0, 4)
    check(abs(dg_stem - (-mfe)) < 0.01, f"5' stem arm incurs full unfolding penalty ({dg_stem} kcal/mol)")

    acc = compute_target_accessibility("GCGCUUUUGCGC", 4, 8)
    check(acc.accessibility_tier is AccessibilityTier.HIGHLY_ACCESSIBLE, "loop classified as highly accessible")


def test_effective_affinity_and_competition() -> None:
    print("\n[biophysics: effective affinity & spliceosomal displacement]")
    eff = compute_effective_binding(duplex_delta_g=-22.0, delta_g_open=3.6)
    check(abs(eff.delta_g_effective - (-18.4)) < 1e-4, "Delta G_effective = -22.0 + 3.6 = -18.4 kcal/mol")
    check(eff.is_favorable and eff.kd_effective_nm < 1e-3, "effective binding remains sub-picomolar")

    disp = compute_steric_displacement(oligo_kd_nm=0.10, competitor_name="u1_snrnp", target_gene="STMN2")
    check(abs(disp.displacement_ratio - 250.0) < 0.1, "U1 snRNP displacement ratio is 250x")
    check(disp.dominant_binder == "oligo" and disp.fractional_occupancy > 0.98, "ASO dominates U1 snRNP")

    weak = compute_steric_displacement(oligo_kd_nm=100.0, competitor_name="u1_snrnp")
    check(weak.dominant_binder == "competitor", "weak oligo is outcompeted by endogenous factor")


def main() -> int:
    print("=" * 70)
    print("SYNAPTIK 3D RNA Biophysics Test Suite")
    print("=" * 70)
    test_parameters()
    test_folding_and_accessibility()
    test_effective_affinity_and_competition()
    print("\n" + "=" * 70)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S):")
        for f in FAILURES:
            print("  -", f)
        return 1
    print("All biophysics checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
