"""Tests for Autonomous N-of-1 Rare Disease & Translational Pipeline (`synaptik nof1`)."""

from synaptik.nof1.pipeline import (
    ActionabilityTier,
    Nof1PatientCase,
    Nof1Pipeline,
    nof1_dossier_to_dict,
    nof1_dossier_to_markdown,
    run_nof1_pipeline,
)
from synaptik.workflows import run_nof1_workflow


def check(desc, cond):
    if cond:
        print(f"  PASS  {desc}")
    else:
        print(f"  FAIL  {desc}")
        raise AssertionError(desc)


class MockArgs:
    def __init__(self, gene="STMN2", patient="PATIENT-TEST", disease="als", json=False, dossier=None):
        self.gene = gene
        self.patient = patient
        self.disease = disease
        self.variant = "cryptic_splice_variant"
        self.sequence = None
        self.nfl = 65.0
        self.div = 28
        self.json = json
        self.dossier = dossier


def run_tests():
    print("=" * 70)
    print("SYNAPTIK autonomous N-of-1 translational pipeline test suite")
    print("=" * 70)

    # ------------------------------------------------------------------
    # 1. Pipeline execution on validated cryptic splice target (STMN2 in ALS)
    # ------------------------------------------------------------------
    print("\n[n-of-1: STMN2 cryptic exon in ALS]")
    case_stmn2 = Nof1PatientCase(
        patient_id="PATIENT-ALS-STMN2",
        target_gene="STMN2",
        disease="als",
        variant_id="cryptic_exon_2a",
        baseline_nfl_pg_ml=68.0,
    )
    dossier = run_nof1_pipeline(case_stmn2)

    check("actionability tier is Tier 1 (IND-Ready)",
          dossier.actionability_tier == ActionabilityTier.TIER_1_CLINICALLY_ACTIONABLE)
    check("stage 1 identifies splice-switching modality",
          "Splice-Switching" in dossier.stage1_splice_and_mechanism["therapeutic_modality"])
    check("stage 2 confirms human post-mortem concordance",
          dossier.stage2_human_concordance["concordance_status"] == "concordant")
    check("stage 2 maps vulnerable human cell type",
          "Motor Neuron" in dossier.stage2_human_concordance["primary_human_cell_type"] or
          "Betz" in dossier.stage2_human_concordance["primary_human_cell_type"])
    check("stage 3 lead ASO Tm in physiological window (55-75 C)",
          55.0 <= dossier.stage3_therapeutic_aso["melting_temp_c"] <= 75.0)
    check("stage 3 lead ASO cleared of biophysical liabilities",
          dossier.stage3_therapeutic_aso["biophysical_liabilities_cleared"] is True)
    check("stage 3 lead ASO has zero 0-mismatch neural off-targets",
          dossier.stage3_therapeutic_aso["off_target_0_mismatch"] == 0)
    check("stage 3 predicted CNS half-life is ~135 days",
          dossier.stage3_therapeutic_aso["predicted_cns_half_life_days"] == 135)
    check("stage 3 IDT format string provided",
          dossier.stage3_therapeutic_aso["vendor_idt_format"].startswith("m"))
    check("stage 4 composite rescue score > 75",
          dossier.stage4_counterfactual_simulation["composite_rescue_score"] >= 75.0)
    check("stage 4 counterfactual simulation clinically viable",
          dossier.stage4_counterfactual_simulation["is_clinically_viable"] is True)
    check("stage 4 cosine alignment positive",
          dossier.stage4_counterfactual_simulation["cosine_rescue_alignment"] > 0.4)
    check("stage 4 collateral toxicity is low (< 0.25)",
          dossier.stage4_counterfactual_simulation["collateral_toxicity_index"] < 0.25)
    check("stage 5 projects significant 48-week NfL reduction (> 30%)",
          dossier.stage5_biomarker_trajectory["maximal_nfl_reduction_pct"] >= 30.0)
    check("stage 5 achieves 60-80% clinical trial sample size reduction",
          60 <= dossier.stage5_biomarker_trajectory["sample_size_reduction_pct"] <= 80)
    check("stage 6 cryptographic SHA-256 fingerprint generated (64 hex)",
          len(dossier.stage6_preclinical_protocol["protocol_sha256"]) == 64)
    check("stage 6 blinding verification token generated",
          dossier.stage6_preclinical_protocol["blinding_token"].startswith("SYN-"))
    check("citation metadata embeds author Francesco Nobile with ORCID",
          dossier.citation["author"] == "Francesco Nobile" and "0009-0001-9477-3332" in dossier.citation["orcid"])

    # ------------------------------------------------------------------
    # 2. Pipeline execution on knockdown target (ATXN2 in ALS)
    # ------------------------------------------------------------------
    print("\n[n-of-1: ATXN2 gapmer in ALS]")
    case_atxn2 = Nof1PatientCase(
        patient_id="PATIENT-ALS-ATXN2",
        target_gene="ATXN2",
        disease="als",
        baseline_nfl_pg_ml=85.0,
    )
    dossier_atxn2 = run_nof1_pipeline(case_atxn2)
    check("ATXN2 routed to RNase H1 gapmer modality",
          "Gapmer" in dossier_atxn2.stage1_splice_and_mechanism["therapeutic_modality"])
    check("ATXN2 intervention is knockdown",
          dossier_atxn2.stage4_counterfactual_simulation["intervention"] == "knockdown")

    # ------------------------------------------------------------------
    # 3. Pipeline execution on developmental maturation gap gene (NEFH)
    # ------------------------------------------------------------------
    print("\n[n-of-1: developmental gap detection]")
    case_dev = Nof1PatientCase(
        patient_id="PATIENT-DEV-TEST",
        target_gene="NEFH",
        disease="als",
    )
    dossier_dev = run_nof1_pipeline(case_dev)
    check("NEFH detected as developmental gap marker",
          dossier_dev.stage2_human_concordance["concordance_status"] == "developmental_gap")
    check("developmental gap assigned Tier 3 Exploratory",
          dossier_dev.actionability_tier == ActionabilityTier.TIER_3_PRECLINICAL_EXPLORATORY)

    # ------------------------------------------------------------------
    # 4. Dossier serialization and Markdown generation
    # ------------------------------------------------------------------
    print("\n[n-of-1: serialization and markdown export]")
    d_dict = nof1_dossier_to_dict(dossier)
    check("dossier_to_dict returns dict with all 6 stages",
          all(f"stage{i}" in str(d_dict.keys()) for i in range(1, 7)))

    md = nof1_dossier_to_markdown(dossier)
    check("markdown dossier contains IND title", "Translational Preclinical IND Dossier" in md)
    check("markdown dossier contains lead oligo sequence", dossier.stage3_therapeutic_aso["lead_oligo_sequence"] in md)
    check("markdown dossier contains vendor IDT format", dossier.stage3_therapeutic_aso["vendor_idt_format"] in md)
    check("markdown dossier contains protocol SHA-256 digest", dossier.stage6_preclinical_protocol["protocol_sha256"] in md)
    check("markdown dossier contains Francesco Nobile citation and ORCID", "Nobile" in md and "0009-0001-9477-3332" in md)

    # ------------------------------------------------------------------
    # 5. CLI workflow and handler smoke test
    # ------------------------------------------------------------------
    print("\n[n-of-1: cli handler]")
    from synaptik.cli.commands.nof1 import cmd_nof1
    ret_default = cmd_nof1(MockArgs())
    check("cmd_nof1 default returns 0", ret_default == 0)

    ret_json = cmd_nof1(MockArgs(json=True))
    check("cmd_nof1 --json returns 0", ret_json == 0)

    print("\n" + "=" * 70)
    print("ALL N-OF-1 PIPELINE TESTS PASSED")
    print("=" * 70)


if __name__ == "__main__":
    run_tests()
