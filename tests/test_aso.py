"""
Regression tests for In Silico Antisense Oligonucleotide (ASO) Engineering.

Checks nearest-neighbor thermodynamics, biophysical risk filters (G-quadruplex,
poly-T/U, hairpins), neural off-target screening, and clinical drug architectures
(Nusinersen-style splice-switching and Tofersen-style 5-10-5 gapmers).

Run: python3 -m tests.test_aso
"""

from __future__ import annotations

import sys

from synaptik.drug.aso import (
    ASOChemistry, ASODesigner, ASOMode,
    _gc_content, _reverse_complement_rna,
    calculate_molecular_weight, compute_thermodynamics,
    format_vendor_spec, predict_hairpin_delta_g,
    predict_self_dimer_score, screen_off_targets,
)

FAILURES: list[str] = []


def check(cond: bool, msg: str) -> None:
    print(f"  {'PASS' if cond else 'FAIL'}  {msg}")
    if not cond:
        FAILURES.append(msg)


def test_thermodynamics() -> None:
    print("\n[aso: thermodynamics]")
    # Nusinersen 18-mer: TCACTTTCATAATGCTGG
    # Target RNA: CCAGCAUUAUGAAAGUGA
    oligo = "TCACTTTCATAATGCTGG"
    target = "CCAGCAUUAUGAAAGUGA"

    dh, ds, dg37, tm = compute_thermodynamics(oligo, target)

    check(-24.0 <= dg37 <= -18.0, f"Nusinersen dG_37 in therapeutic range (-21.4 kcal/mol expected, got {dg37})")
    check(60.0 <= tm <= 70.0, f"Nusinersen Tm with 2'-MOE in clinical range 60-70C (got {tm}C)")
    check(dh < -120.0 and ds < -350.0, "enthalpy and entropy are negative and favorable for duplex formation")

    # High GC vs Low GC
    gc_high = "GCGCGCGCGCGCGCGCGC"
    gc_low = "ATATATATATATATATAT"
    _, _, _, tm_high = compute_thermodynamics(gc_high, "C" * 18)
    _, _, _, tm_low = compute_thermodynamics(gc_low, "A" * 18)
    check(tm_high > tm_low + 20.0, f"high GC duplex Tm ({tm_high}C) far exceeds low GC duplex ({tm_low}C)")


def test_sequence_manipulation() -> None:
    print("\n[aso: sequence manipulation]")
    # RNA target 5' -> 3': AUGCAU
    # Complement 3' -> 5': UACGUA (TACGTA)
    # Reverse complement 5' -> 3': AUGCAU (ATGCAT in DNA)
    rc = _reverse_complement_rna("AUGCAU")
    check(rc == "ATGCAT", f"palindromic sequence reverse-complements correctly (got {rc})")

    rc2 = _reverse_complement_rna("AAAACCCC")
    check(rc2 == "GGGGTTTT", f"directionality 5'->3' preserved as reverse complement (got {rc2})")

    check(_gc_content("GGCCAA") == 66.7, "GC content rounds to 1 decimal place")
    check(_gc_content("AAAAAA") == 0.0, "zero GC content computed correctly")


def test_biophysical_filters() -> None:
    print("\n[aso: biophysical liability filters]")
    designer = ASODesigner()

    # G-quadruplex: target with CCCC produces oligo with GGGG
    target_c = "ACGTCCCCACGT"
    res_quad = designer.design(target_c, mode=ASOMode.SPLICE_SWITCH, window_sizes=(10,))
    check(any(c.has_quadruplex_risk for c in res_quad.candidates),
          "candidate containing GGGG flagged for G-quadruplex neurotoxicity risk")

    # Poly-T / Poly-U: target with AAAA produces oligo with TTTT
    target_a = "ACGTAAAAACGT"
    res_u = designer.design(target_a, mode=ASOMode.SPLICE_SWITCH, window_sizes=(10,))
    check(any(c.has_poly_u_risk for c in res_u.candidates),
          "candidate containing TTTT flagged for premature termination risk")

    # Self-complementary hairpin
    # Stem: CGCG, Loop: TTTT, Stem: CGCG -> CGCGTTTTCGCG
    hairpin_seq = "CGCGTTTTCGCG"
    dg_hp = predict_hairpin_delta_g(hairpin_seq)
    check(dg_hp <= -3.0, f"strong hairpin predicted with negative free energy (got {dg_hp} kcal/mol)")

    # Dimer score
    dimer_seq = "ACGTACGTACGT"
    dimer_s = predict_self_dimer_score(dimer_seq)
    check(dimer_s >= 4.0, f"self-dimerizing palindromic sequence scores high dimer risk (got {dimer_s})")


def test_off_target_screening() -> None:
    print("\n[aso: off-target screening]")
    # An oligo identical to ACTB housekeeping gene must trigger off-target alert
    actb_sub = "CCACCATGTACCCAGGCATT"  # 20 nt from ACTB
    # An oligo complementary to ACTB
    oligo_actb = _reverse_complement_rna(actb_sub)

    off_0, off_1, min_dist = screen_off_targets(oligo_actb, target_gene="OTHER_GENE")
    check(off_0 >= 1, "exact match in core neural transcriptome detected as 0-mismatch off-target")

    # Clean oligo from synthetic random-like sequence
    clean_oligo = "CTGAAAGCGAGCATGTCAGT"
    off_0_clean, off_1_clean, min_dist_clean = screen_off_targets(clean_oligo, target_gene="TEST")
    check(off_0_clean == 0, "clean synthetic oligo has zero exact off-targets in neural transcriptome")
    check(min_dist_clean >= 4, f"clean oligo maintains high mismatch distance to neural transcripts (got {min_dist_clean})")


def test_vendor_and_architectures() -> None:
    print("\n[aso: vendor formatting & clinical architectures]")
    seq = "TCACTTTCATAATGCTGG"  # 18-mer

    # Splice-switching uniform 2'-MOE with PS
    moe_spec = format_vendor_spec(seq, ASOMode.SPLICE_SWITCH, ASOChemistry.MOE_PS)
    check(moe_spec.startswith("mU*mC*") or moe_spec.startswith("mT*mC*"),
          "splice-switching ASO formatted with uniform 2'-MOE and phosphorothioate linkages")

    # Gapmer 5-10-5
    gapmer_20 = "CATTACTTTCCTTCTGCTCG"
    gapmer_spec = format_vendor_spec(gapmer_20, ASOMode.GAPMER, ASOChemistry.GAPMER_5_10_5)
    check("+d" in gapmer_spec and gapmer_spec.startswith("mC*mA*mU*mU*mA*"),
          "gapmer formatted with 5' MOE wing, central +d DNA gap, and 3' MOE wing")

    # PMO
    pmo_spec = format_vendor_spec(seq, ASOMode.SPLICE_SWITCH, ASOChemistry.PMO)
    check(pmo_spec.startswith("5'-PMO-") and pmo_spec.endswith("-3'"),
          "PMO formatted with standard morpholino notation")

    # Molecular weight
    mw_moe = calculate_molecular_weight(seq, ASOChemistry.MOE_PS)
    check(6800.0 <= mw_moe <= 7500.0, f"18-mer 2'-MOE molecular weight in accurate range (got {mw_moe} Da)")


def test_real_world_targets() -> None:
    print("\n[aso: real-world disease target validation]")
    designer = ASODesigner()

    # 1. STMN2 cryptic exon (ALS/FTD)
    stmn2_seq = "TTCTGTAACAGTTTTTGTTTTTTCAGGTGAGTATTCAGGT"
    stmn2_res = designer.design(stmn2_seq, mode=ASOMode.SPLICE_SWITCH, target_gene="STMN2")

    check(stmn2_res.n_windows_evaluated > 30, "tiling generates comprehensive candidate window series")
    check(stmn2_res.n_passed_filters > 0, "filter pipeline retains high-quality candidate oligos")
    lead_stmn2 = stmn2_res.candidates[0]
    check(lead_stmn2.composite_score >= 85.0, f"lead STMN2 candidate achieves high composite score (got {lead_stmn2.composite_score})")
    check(lead_stmn2.melting_temp_c >= 58.0, f"lead candidate has strong hybridization Tm (got {lead_stmn2.melting_temp_c}C)")
    check("CAGGTGAG" in lead_stmn2.target_sequence_5_to_3 or "CAGGUGAG" in lead_stmn2.target_sequence_5_to_3,
          "lead candidate targets the STMN2 cryptic donor splice junction")

    # 2. SOD1 gapmer knockdown (ALS)
    sod1_seq = "ATGGCGACGAAGGCCGTGTGCGTGCTGAAGGGCGACGGCCCAGTGCAGGGCATCATCAATTTCGAGCAGAAGGAAAGTAATGGACC"
    sod1_res = designer.design(sod1_seq, mode=ASOMode.GAPMER, target_gene="SOD1")
    lead_sod1 = sod1_res.candidates[0]
    check(lead_sod1.chemistry == ASOChemistry.GAPMER_5_10_5, "gapmer mode selects 5-10-5 architecture")
    check(lead_sod1.composite_score >= 85.0, f"lead SOD1 gapmer achieves high composite score (got {lead_sod1.composite_score})")
    check(lead_sod1.off_target_0_mismatch == 0, "zero 0-mismatch off-targets in neural transcriptome")


def test_aso_tiling_benchmark() -> None:
    print("\n[aso: splice-regulatory context & published tiling benchmark]")
    from synaptik.drug.aso import score_splice_regulatory_context
    from synaptik.drug.benchmark_aso import run_aso_tiling_benchmark

    bonus, hits = score_splice_regulatory_context("UUUUUUUUCAGGUGAGUAUU")
    check(bonus >= 4.0, f"detects PPT + splice junction + silencer (bonus={bonus})")
    check("polypyrimidine tract" in hits and "splice junction consensus" in hits, f"reports motifs: {hits}")

    report = run_aso_tiling_benchmark()
    check(report.n_loci == 4, f"evaluates 4 published loci (got {report.n_loci})")
    by_gene = {o.gene: o for o in report.outcomes}
    check(by_gene["SMN2"].best_rank_with_biophysics == 1, f"SMN2 ISS-N1 ranked #1 (got {by_gene['SMN2'].best_rank_with_biophysics})")
    check(by_gene["STMN2"].best_rank_with_biophysics == 1, f"STMN2 cryptic exon 2a ranked #1 (got {by_gene['STMN2'].best_rank_with_biophysics})")
    check(report.mrr_biophysics >= 0.5, f"biophysics MRR >= 0.50 (got {report.mrr_biophysics})")


def test_calibrated_empirical_knockdown() -> None:
    print("\n[aso: ASO-Atlas calibrated empirical knockdown]")
    from synaptik.drug.aso_calibrated import CALIBRATION_METADATA, predict_empirical_knockdown_pct

    check(CALIBRATION_METADATA["n_train"] >= 131000, "records >=131k non-neuro training pairs")
    check(CALIBRATION_METADATA["heldout_neuro_auroc_70pct"] > 0.70, "held-out neuro AUROC > 0.70")

    designer = ASODesigner()
    target = "AUGGCUAGCUAGCUAGCUAGCUGACUGACUGACUGACUGACUGACU"
    gap_res = designer.design(target, mode=ASOMode.GAPMER, window_sizes=(20,))
    check(len(gap_res.candidates) > 0, "produces gapmer candidates")
    kd = gap_res.candidates[0].empirical_knockdown_pct
    check(kd is not None and 0.0 <= kd <= 100.0, f"gapmer gets bounded empirical_knockdown_pct (got {kd})")

    sso_res = designer.design(target, mode=ASOMode.SPLICE_SWITCH, window_sizes=(18,))
    check(len(sso_res.candidates) > 0 and sso_res.candidates[0].empirical_knockdown_pct is None,
          "splice-switching ASO leaves empirical_knockdown_pct=None")

    kd_direct = predict_empirical_knockdown_pct("TCACTTTCAACCACATGGTG", -22.5, 62.0, -1.2, 2.0, 0)
    check(0.0 <= kd_direct <= 100.0, f"direct predictor bounded (got {kd_direct})")


def main() -> int:
    print("=" * 70)
    print("SYNAPTIK antisense oligonucleotide (ASO) engineering test suite")
    print("=" * 70)
    test_thermodynamics()
    test_sequence_manipulation()
    test_biophysical_filters()
    test_off_target_screening()
    test_vendor_and_architectures()
    test_real_world_targets()
    test_aso_tiling_benchmark()
    test_calibrated_empirical_knockdown()
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
