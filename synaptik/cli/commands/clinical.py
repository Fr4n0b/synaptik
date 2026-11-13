"""
`synaptik clinical`: First-in-Human (FIH) clinical protocol & adaptive dose escalation.

Principal Investigator & Author:
Francesco Nobile (ORCID: 0009-0001-9477-3332)
Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School;
University of Pavia.
"""

from __future__ import annotations

import json
from pathlib import Path

from ...clinical.allometry import calculate_allometric_dose
from ...clinical.escalation import (
    BOINEscalationEngine,
    CohortOutcome,
    MTPIEngine,
    ThreePlusThreeEngine,
)
from ...clinical.pkpd import simulate_intrathecal_pkpd
from ...clinical.protocol import generate_clinical_protocol
from ..main import Fmt, header

AUTHOR_CITE = "Francesco Nobile (ORCID: 0009-0001-9477-3332)"
CLINICAL_DISCLAIMER = "RESEARCH PLANNING ONLY — NOT A CLINICAL DOSE RECOMMENDATION"


def _print_safety_banner() -> None:
    print(Fmt.yellow("!" * 78))
    print(Fmt.yellow("! RESEARCH PLANNING ONLY — NOT A CLINICAL DOSE RECOMMENDATION"))
    print(Fmt.dim("! Intrathecal dose calculations are mathematical translations based on"))
    print(Fmt.dim("! physiological CSF volume ratios. All human starting doses must be empirically"))
    print(Fmt.dim("! confirmed via GLP toxicology and authorized under an approved IND/CTA."))
    print(Fmt.yellow("!" * 78))
    print()


def _print_footer() -> None:
    print("\n" + "-" * 78)
    print(f"Cite: {AUTHOR_CITE}; run 'synaptik cite'")


# ---------------------------------------------------------------------------
# Subaction 1: Allometry
# ---------------------------------------------------------------------------


def cmd_clinical_allometry(args) -> int:
    noael_mg = getattr(args, "noael_mg", None)
    noael_mg_kg = getattr(args, "noael_mg_kg", None)
    nhp_weight = float(getattr(args, "nhp_weight", 3.0) or 3.0)
    nhp_vol = float(getattr(args, "nhp_volume", 15.0) or 15.0)
    adult_vol = float(getattr(args, "adult_volume", 140.0) or 140.0)
    pediatric_vol = float(getattr(args, "pediatric_volume", 65.0) or 65.0)
    sf = float(getattr(args, "safety_factor", 10.0) or 10.0)

    if noael_mg is None and noael_mg_kg is None:
        noael_mg = 30.0

    try:
        res_adult = calculate_allometric_dose(
            noael_mg=noael_mg,
            noael_mg_kg=noael_mg_kg,
            nhp_weight_kg=nhp_weight,
            nhp_volume_ml=nhp_vol,
            human_volume_ml=adult_vol,
            safety_factor=sf,
            cohort_name="adult",
            acknowledge_not_clinical=True,
        )
        res_pedia = calculate_allometric_dose(
            noael_mg=noael_mg,
            noael_mg_kg=noael_mg_kg,
            nhp_weight_kg=nhp_weight,
            nhp_volume_ml=nhp_vol,
            human_volume_ml=pediatric_vol,
            safety_factor=sf,
            cohort_name="pediatric",
            acknowledge_not_clinical=True,
        )
    except ValueError as e:
        print(Fmt.red(f"Allometry calculation error: {e}"))
        return 1

    if getattr(args, "json", False):
        out = {
            "author": AUTHOR_CITE,
            "noael_nhp_mg": res_adult.noael_nhp_mg,
            "nhp_csf_volume_ml": nhp_vol,
            "safety_factor": sf,
            "adult": res_adult.to_dict(),
            "pediatric": res_pedia.to_dict(),
            "methodology": "CSF volume allometric scaling (21 CFR Part 312 / FDA 2005)",
        }
        print(json.dumps(out, indent=2))
        return 0

    header("Intrathecal Allometric Dose Scaling: Nonclinical NHP to First-in-Human")
    print(f"Nonclinical Species         : Cynomolgus Macaque (NHP)")
    print(f"NHP NOAEL Dose              : {res_adult.noael_nhp_mg:.1f} mg")
    print(f"NHP CSF Volume              : {nhp_vol:.1f} mL")
    print(f"FDA Safety Factor           : {sf:.1f}x")
    print()
    print("=" * 78)
    print("Scalings Across Human Age Compartments")
    print("=" * 78)
    print(
        f"Adult Cohort (140 mL CSF)   : Ratio = {res_adult.scaling_ratio:.2f}x | "
        f"HED = {res_adult.hed_mg:.1f} mg | MRSD = {Fmt.green(f'{res_adult.mrsd_mg:.1f} mg')}"
    )
    print(
        f"Pediatric Infant (65 mL CSF): Ratio = {res_pedia.scaling_ratio:.2f}x | "
        f"HED = {res_pedia.hed_mg:.1f} mg | MRSD = {Fmt.green(f'{res_pedia.mrsd_mg:.1f} mg')}"
    )
    print()
    print(Fmt.dim("Regulatory Rationale: Intrathecal ASOs distribute strictly within CSF;"))
    print(Fmt.dim("systemic BSA scaling underestimates target neural tissue exposure."))
    _print_footer()
    return 0


# ---------------------------------------------------------------------------
# Subaction 2: PK/PD
# ---------------------------------------------------------------------------


def cmd_clinical_pkpd(args) -> int:
    dose_mg = float(getattr(args, "dose_mg", 50.0) or 50.0)
    route = str(getattr(args, "route", "bolus") or "bolus")
    csf_t12 = float(getattr(args, "csf_half_life", 24.0) or 24.0)
    tissue_t12 = float(getattr(args, "tissue_half_life", 135.0) or 135.0)
    days = int(getattr(args, "days", 360) or 360)
    infusion_h = float(getattr(args, "infusion_hours", 24.0) or 24.0)

    try:
        res = simulate_intrathecal_pkpd(
            dose_mg=dose_mg,
            route=route,
            csf_half_life_days=csf_t12,
            tissue_half_life_days=tissue_t12,
            duration_days=days,
            infusion_hours=infusion_h,
        )
    except ValueError as e:
        print(Fmt.red(f"PK/PD simulation error: {e}"))
        return 1

    if getattr(args, "json", False):
        print(json.dumps(res.to_dict(), indent=2))
        return 0

    header(f"Intrathecal PK/PD Simulation: {res.dose_mg:.1f} mg via {res.route.upper()}")
    print(f"CSF Half-Life (t1/2)        : {res.kinetics.csf_half_life_days:.1f} days (k_el = {res.kinetics.k_el_csf:.5f}/d)")
    print(f"Tissue Half-Life (t1/2)     : {res.kinetics.tissue_half_life_days:.1f} days (k_el = {res.kinetics.k_el_tissue:.5f}/d)")
    print(f"Peak CSF Concentration      : {res.exposure.cmax_csf_ug_ml:.1f} µg/mL (Tmax: {res.exposure.tmax_csf_days:.3f} d)")
    print(f"Peak CNS Tissue Exposure    : {res.exposure.cmax_tissue_ug_g:.1f} µg/g (Tmax: Day {res.exposure.tmax_tissue_days:.0f})")
    print(f"Total Tissue Exposure (AUC) : {res.exposure.auc_tissue_ug_day_g:.1f} µg·day/g")
    print()
    print("=" * 78)
    print("Day 90 Fluid Biomarker Response (Plasma NfL)")
    print("=" * 78)
    verdict_style = Fmt.green if res.day90_milestone_met else Fmt.red
    print(f"Day 90 Predicted NfL Level  : {res.day90_nfl_pg_ml:.1f} pg/mL")
    print(f"Day 90 Suppression          : {verdict_style(f'-{res.day90_nfl_reduction_pct:.1f}%')}")
    print(f"FDA Accelerated Milestone   : {verdict_style('MET (>30% reduction)' if res.day90_milestone_met else 'NOT MET (<=30% reduction)')}")

    if res.route == "catheter":
        print()
        print(Fmt.dim("Catheter Micro-infusion Advantage: CSF peak concentration flattened by >60%"))
        print(Fmt.dim("minimizing acute meningismus/radiculopathy risk while preserving tissue AUC."))

    _print_footer()
    return 0


# ---------------------------------------------------------------------------
# Subaction 3: Escalate
# ---------------------------------------------------------------------------


def cmd_clinical_escalate(args) -> int:
    design = str(getattr(args, "design", "boin") or "boin").lower()
    target_phi = float(getattr(args, "target_phi", 0.25) or 0.25)
    doses = getattr(args, "doses", None) or [10.0, 25.0, 50.0, 100.0, 150.0]
    history_str = getattr(args, "history", None)

    # Initialize engine
    if design == "boin":
        engine = BOINEscalationEngine(target_phi=target_phi, dose_levels=doses)
    elif design in ("3+3", "three_plus_three"):
        engine = ThreePlusThreeEngine(dose_levels=doses)
    elif design == "mtpi":
        engine = MTPIEngine(target_phi=target_phi, dose_levels=doses)
    else:
        print(Fmt.red(f"Unknown escalation design: {design}"))
        return 1

    # Parse cohort history if supplied
    eval_level = 1
    eval_n = 3
    eval_dlts = 0
    history_records: list[CohortOutcome] = []

    if history_str:
        try:
            parsed = json.loads(history_str)
            if isinstance(parsed, list) and len(parsed) > 0:
                for item in parsed:
                    rec = CohortOutcome(
                        dose_level=int(item.get("dose_level", 1)),
                        dose_mg=float(item.get("dose_mg", doses[int(item.get("dose_level", 1)) - 1])),
                        n_enrolled=int(item.get("n", item.get("n_patients", item.get("n_enrolled", 3)))),
                        n_dlts=int(item.get("dlts", item.get("n_dlts", 0))),
                    )
                    history_records.append(rec)
                last = history_records[-1]
                eval_level = last.dose_level
                eval_n = last.n_enrolled
                eval_dlts = last.n_dlts
        except Exception as e:
            print(Fmt.red(f"Invalid JSON format for --history: {e}"))
            return 1

    try:
        decision = engine.evaluate_cohort(
            dose_level=eval_level,
            n=eval_n,
            dlts=eval_dlts,
            history=history_records,
        )
    except ValueError as e:
        print(Fmt.red(f"Escalation evaluation error: {e}"))
        return 1

    if getattr(args, "json", False):
        out = decision.to_dict()
        out["history"] = [h.to_dict() for h in history_records]
        print(json.dumps(out, indent=2))
        return 0

    header(f"Adaptive Dose Escalation Engine: {decision.design.upper()} Design")
    print(f"Target Toxicity Rate (phi)  : {decision.target_phi:.2f}")
    if decision.design == "boin":
        print(f"Escalation Boundary (l_e)   : {decision.lambda_e:.4f}")
        print(f"De-escalation Boundary (l_d): {decision.lambda_d:.4f}")
    print(f"Current Evaluated Cohort    : Level {decision.current_dose_level} ({decision.current_dose_mg:.1f} mg)")
    print(f"Observed DLTs at Level      : {eval_dlts} / {eval_n} ({decision.observed_toxicity_rate:.1%})")
    print(f"Posterior Overdose Prob     : {decision.posterior_overdose_prob:.4f}")
    print()
    print("=" * 78)
    print("Bedside Recommendation")
    print("=" * 78)

    action_style = Fmt.green if decision.action.value == "escalate" else (Fmt.yellow if decision.action.value == "retain" else Fmt.red)
    print(f"Recommended Action          : {action_style(decision.action.value.upper())}")
    if decision.next_dose_level is not None:
        print(f"Next Cohort Dose            : Level {decision.next_dose_level} ({decision.next_dose_mg:.1f} mg)")
    if decision.mtd_declared_dose_level is not None:
        print(f"Declared MTD Level          : Level {decision.mtd_declared_dose_level} ({decision.mtd_declared_dose_mg:.1f} mg)")
    print(f"Clinical Rationale          : {decision.rationale}")

    if decision.is_dose_eliminated:
        print(Fmt.red(f"Safety Stopping Alert       : Dose levels {decision.eliminated_dose_levels} permanently eliminated."))

    _print_footer()
    return 0


# ---------------------------------------------------------------------------
# Subaction 4: Protocol
# ---------------------------------------------------------------------------


def cmd_clinical_protocol(args) -> int:
    gene = str(getattr(args, "gene", "STMN2") or "STMN2")
    disease = str(getattr(args, "disease", "als") or "als")
    noael_mg = float(getattr(args, "noael_mg", 30.0) or 30.0)
    design = str(getattr(args, "design", "boin") or "boin")
    target_phi = float(getattr(args, "target_phi", 0.25) or 0.25)
    doses = getattr(args, "doses", None) or [10.0, 25.0, 50.0, 100.0, 150.0]
    out_file = getattr(args, "output", None)

    try:
        proto = generate_clinical_protocol(
            target_gene=gene,
            disease=disease,
            noael_nhp_mg=noael_mg,
            design=design,
            target_phi=target_phi,
            planned_doses=doses,
        )
    except Exception as e:
        print(Fmt.red(f"Protocol generation error: {e}"))
        return 1

    if out_file:
        try:
            Path(out_file).write_text(proto.markdown, encoding="utf-8")
            print(Fmt.green(f"Protocol written to {out_file}"))
        except Exception as e:
            print(Fmt.red(f"Failed to write protocol to {out_file}: {e}"))
            return 1

    if getattr(args, "json", False):
        print(json.dumps(proto.to_dict(), indent=2))
        return 0

    header(f"FIH Clinical Protocol: {proto.target_gene} in {proto.indication.upper()}")
    print(f"Protocol Identifier         : {proto.protocol_id}")
    print(f"Investigational Product     : {proto.product_name}")
    print(f"Human Starting Dose (MRSD)  : {Fmt.green(f'{proto.mrsd_mg:.1f} mg')} (from {proto.noael_nhp_mg:.1f} mg NHP NOAEL)")
    print(f"Escalation Design           : {proto.escalation_design.upper()}")
    print(f"Biomarker Surrogate Gate    : {proto.primary_biomarker_milestone}")
    print()
    print("=" * 78)
    print("Bedside Dose-Limiting Toxicity (DLT) Criteria")
    print("=" * 78)
    for crit in proto.dlt_monitoring_criteria:
        print(f"  • {crit}")

    _print_footer()
    return 0


# ---------------------------------------------------------------------------
# Subparser Registration
# ---------------------------------------------------------------------------


def register(sub, common) -> None:
    """Register `synaptik clinical` and its four subactions."""
    desc = (
        "First-in-Human (FIH) Clinical Protocol & Dose Escalation Engine.\n\n"
        f"Author: {AUTHOR_CITE}\n"
        "Affiliation: Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School; University of Pavia\n"
        "GitHub: https://github.com/Fr4n0b/synaptik"
    )
    p = sub.add_parser(
        "clinical",
        parents=[common],
        help="First-in-Human trial protocol, CSF allometry, intrathecal PK/PD, and BOIN escalation",
        description=desc,
    )
    actions = p.add_subparsers(dest="clinical_action", required=True)

    # 1. Allometry
    s1 = actions.add_parser(
        "allometry",
        parents=[common],
        help="CSF volume allometric dose scaling from NHP NOAEL to human adult and pediatric",
    )
    s1.add_argument("--noael-mg", type=float, help="NHP absolute NOAEL dose in mg")
    s1.add_argument("--noael-mg-kg", type=float, help="NHP body weight NOAEL dose in mg/kg")
    s1.add_argument("--nhp-weight", type=float, default=3.0, help="NHP weight in kg (default: 3.0)")
    s1.add_argument("--nhp-volume", type=float, default=15.0, help="NHP CSF volume in mL (default: 15.0)")
    s1.add_argument("--adult-volume", type=float, default=140.0, help="Adult CSF volume in mL (default: 140.0)")
    s1.add_argument("--pediatric-volume", type=float, default=65.0, help="Pediatric CSF volume in mL (default: 65.0)")
    s1.add_argument("--safety-factor", type=float, default=10.0, help="FDA Safety factor (default: 10.0)")
    s1.add_argument("--json", action="store_true", help="Emit structured JSON")
    s1.set_defaults(func=cmd_clinical_allometry)

    # 2. PK/PD
    s2 = actions.add_parser(
        "pkpd",
        parents=[common],
        help="Coupled 2-compartment intrathecal PK/PD simulation (bolus vs catheter)",
    )
    s2.add_argument("--dose-mg", type=float, default=50.0, help="Administered dose in mg")
    s2.add_argument("--route", choices=["bolus", "catheter"], default="bolus", help="Administration route")
    s2.add_argument("--csf-half-life", type=float, default=24.0, help="CSF terminal half-life in days")
    s2.add_argument("--tissue-half-life", type=float, default=135.0, help="CNS tissue retention half-life in days")
    s2.add_argument("--days", type=int, default=360, help="Simulation duration in days")
    s2.add_argument("--infusion-hours", type=float, default=24.0, help="Catheter infusion duration in hours")
    s2.add_argument("--json", action="store_true", help="Emit structured JSON")
    s2.set_defaults(func=cmd_clinical_pkpd)

    # 3. Escalate
    s3 = actions.add_parser(
        "escalate",
        parents=[common],
        help="Adaptive dose escalation decisions under BOIN, 3+3, or mTPI",
    )
    s3.add_argument("--design", choices=["boin", "3+3", "mtpi"], default="boin", help="Escalation algorithm")
    s3.add_argument("--target-phi", type=float, default=0.25, help="Target toxicity rate phi")
    s3.add_argument("--doses", nargs="+", type=float, default=[10.0, 25.0, 50.0, 100.0, 150.0], help="Dose ladder in mg")
    s3.add_argument("--cohorts", type=int, default=6, help="Maximum number of cohorts")
    s3.add_argument("--cohort-size", type=int, default=3, help="Patients per cohort")
    s3.add_argument("--history", type=str, help="Cohort history as JSON string")
    s3.add_argument("--json", action="store_true", help="Emit structured JSON")
    s3.set_defaults(func=cmd_clinical_escalate)

    # 4. Protocol
    s4 = actions.add_parser(
        "protocol",
        parents=[common],
        help="Generate Phase 1 First-in-Human clinical trial protocol dossier",
    )
    s4.add_argument("--gene", default="STMN2", help="Target gene symbol")
    s4.add_argument("--disease", choices=["als", "ftd", "ad"], default="als", help="Disease indication")
    s4.add_argument("--noael-mg", type=float, default=30.0, help="NHP NOAEL dose in mg")
    s4.add_argument("--design", choices=["boin", "3+3", "mtpi"], default="boin", help="Escalation design")
    s4.add_argument("--target-phi", type=float, default=0.25, help="Target toxicity rate phi")
    s4.add_argument("--doses", nargs="+", type=float, default=[10.0, 25.0, 50.0, 100.0, 150.0], help="Planned doses in mg")
    s4.add_argument("--output", type=str, help="Export protocol to Markdown file")
    s4.add_argument("--json", action="store_true", help="Emit structured JSON")
    s4.set_defaults(func=cmd_clinical_protocol)
