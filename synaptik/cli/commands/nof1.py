"""`synaptik nof1`: Autonomous N-of-1 Rare Disease & Translational Therapeutics Pipeline."""

from __future__ import annotations

import json
from pathlib import Path
from ..main import Fmt, header
from ...workflows import run_nof1_workflow


def cmd_nof1(args) -> int:
    gene = args.gene.strip().upper()
    patient = getattr(args, "patient", None) or "PATIENT-N1"
    disease = getattr(args, "disease", "als")
    variant = getattr(args, "variant", "cryptic_splice_variant")
    seq = getattr(args, "sequence", None)
    nfl = getattr(args, "nfl", 65.0)
    div = getattr(args, "div", 28)

    dossier_dict, dossier_md = run_nof1_workflow(
        patient_id=patient,
        target_gene=gene,
        disease=disease,
        variant_id=variant,
        target_sequence=seq,
        baseline_nfl=nfl,
        div=div,
    )

    if getattr(args, "json", False):
        print(json.dumps(dossier_dict, indent=2))
        return 0

    if getattr(args, "dossier", None):
        out_path = Path(args.dossier)
        try:
            out_path.write_text(dossier_md, encoding="utf-8")
            print(Fmt.green(f"Exported IND-enabling preclinical dossier to {out_path}"))
        except OSError as e:
            print(Fmt.red(f"Failed to export dossier: {e}"))

    tier_str = dossier_dict["actionability_tier"]
    if "Tier 1" in tier_str:
        colored_tier = Fmt.green(tier_str)
    elif "Tier 2" in tier_str or "Tier 3" in tier_str:
        colored_tier = Fmt.yellow(tier_str)
    else:
        colored_tier = Fmt.red(tier_str)

    s1 = dossier_dict["stage1_splice_and_mechanism"]
    s2 = dossier_dict["stage2_human_concordance"]
    s3 = dossier_dict["stage3_therapeutic_aso"]
    s4 = dossier_dict["stage4_counterfactual_simulation"]
    s5 = dossier_dict["stage5_biomarker_trajectory"]
    s6 = dossier_dict["stage6_preclinical_protocol"]

    header(f"Autonomous N-of-1 Translational Pipeline: {patient} ({gene})")
    print(f"Target Gene / Disease Context : {gene} | {s1['disease_context']}")
    print(f"Translational Actionability   : {colored_tier}")
    print(f"Therapeutic Modality Routed   : {s1['therapeutic_modality']}")
    print(f"\nExecutive Verdict:\n  {dossier_dict['executive_verdict']}\n")

    print("Stage 1-2: Genomic Mechanism & Human Brain snRNA-seq Concordance")
    print("=" * 78)
    print(f"  Molecular Mechanism Class   : {s1['mechanism_class']}")
    print(f"  Human In Vivo Concordance   : {s2['concordance_status'].upper()} (Score: {s2['concordance_score_pct']:.0f}%)")
    print(f"  Vulnerable Brain Cell Type  : {s2['primary_human_cell_type']}")
    print(f"  Single-Cell Evidence Source : {s2['evidence_source']}")

    print("\nStage 3: Lead In Silico Therapeutic ASO Candidate")
    print("=" * 78)
    print(f"  Lead Oligo (5'->3')         : {s3['lead_oligo_sequence']} ({s3['chemistry']})")
    print(f"  Duplex Thermodynamics       : Tm = {s3['melting_temp_c']} C | dG_37 = {s3['gibbs_free_energy_37c']} kcal/mol | GC = {s3['gc_content_pct']}%")
    print(f"  Biophysical Liabilities     : {'CLEARED' if s3['biophysical_liabilities_cleared'] else 'WARNING'}")
    print(f"  Neural Off-Target Profile   : 0-mismatch: {s3['off_target_0_mismatch']}, 1-mismatch: {s3['off_target_1_mismatch']}")
    print(f"  Predicted Human CSF t_1/2   : ~{s3['predicted_cns_half_life_days']} days (intrathecal)")
    print(f"  Vendor Order Specification  : {s3['vendor_idt_format']}")

    print("\nStage 4: Causal Counterfactual Simulation (Transcriptomic Rescue)")
    print("=" * 78)
    print(f"  Vector Rescue Efficiency    : {s4['rescue_efficiency_pct']:.1f}%")
    print(f"  Cosine Directional Alignment: {s4['cosine_rescue_alignment']:.3f} (-1 to +1)")
    print(f"  Collateral Toxicity Index   : {s4['collateral_toxicity_index']:.2f} (0=clean, 1=lethal)")
    print(f"  Composite Rescue Score      : {s4['composite_rescue_score']:.1f} / 100 ({'VIABLE' if s4['is_clinically_viable'] else 'NON-VIABLE'})")

    print("\nStage 5: Longitudinal Fluid Biomarker Trajectory (Plasma NfL)")
    print("=" * 78)
    print(f"  Baseline Plasma NfL         : {s5['baseline_nfl_pg_ml']} pg/mL ({s5['patient_progression_tier']} Progressor)")
    print(f"  Projected In Vivo Trajectory: Week 12: {s5['predicted_week12_nfl']:.1f} pg/mL | Week 24: {s5['predicted_week24_nfl']:.1f} pg/mL | Week 48: {s5['predicted_week48_nfl']:.1f} pg/mL")
    print(f"  Maximal NfL Suppression     : -{s5['maximal_nfl_reduction_pct']:.1f}% ({s5['clinical_surrogate_milestone']})")
    print(f"  Trial Cohort Optimization   : {s5['sample_size_reduction_pct']}% patient reduction per arm via stratification")

    print("\nStage 6: Pre-Registered Preclinical Study Protocol")
    print("=" * 78)
    print(f"  Cryptographic SHA-256 Digest: {s6['protocol_sha256']}")
    print(f"  Blinding Verification Token : {s6['blinding_token']}")
    print(f"  Required Paired Sample Size : {s6['recommended_pairs']} matched isogenic pairs (assumed rho = {s6['inter_batch_correlation_rho']:.2f})")
    print(f"  Plate Evaporation Barrier   : {'VERIFIED (PBS perimeter wells)' if s6['evaporation_barrier_verified'] else 'MISSING'}")

    print("\n" + "-" * 78)
    print("Cite: Francesco Nobile (ORCID: 0009-0001-9477-3332); run 'synaptik cite'")
    return 0


def register(sub, common) -> None:
    s = sub.add_parser(
        "nof1", parents=[common],
        help="autonomous N-of-1 rare disease & translational therapeutics pipeline",
    )
    s.add_argument("--gene", required=True, metavar="GENE",
                   help="target gene symbol (e.g. STMN2, UNC13A, ATXN2, SOD1, SCN2A)")
    s.add_argument("--patient", default="PATIENT-N1",
                   help="patient or cohort identifier (default: PATIENT-N1)")
    s.add_argument("--disease", choices=["als", "ftd", "ad"], default="als",
                   help="disease indication (default: als)")
    s.add_argument("--variant", default="cryptic_splice_variant",
                   help="variant description or HGVS string (default: cryptic_splice_variant)")
    s.add_argument("--sequence", metavar="SEQ",
                   help="target RNA sequence window (auto-resolved from curated gold standard if omitted)")
    s.add_argument("--nfl", type=float, default=65.0,
                   help="baseline plasma NfL in pg/mL (default: 65.0)")
    s.add_argument("--div", type=int, default=28,
                   help="in vitro neuronal maturity DIV (default: 28)")
    s.add_argument("--dossier", metavar="PATH",
                   help="file path to export full IND-enabling preclinical dossier in Markdown")
    s.add_argument("--json", action="store_true", help="output complete dossier as JSON")
    s.set_defaults(func=cmd_nof1)
