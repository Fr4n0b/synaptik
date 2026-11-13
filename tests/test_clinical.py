"""
Comprehensive Test Suite: First-in-Human (FIH) Clinical Protocol & Dose Escalation Engine.

Covers:
1. CSF Volume Allometric Scaling (`synaptik.clinical.allometry`)
2. Intrathecal 2-Compartment Coupled PK/PD Simulator (`synaptik.clinical.pkpd`)
3. Adaptive Dose Escalation: 3+3, BOIN, mTPI (`synaptik.clinical.escalation`)
4. Bedside DLT Stopping Boundaries & Protocol Generator (`synaptik.clinical.protocol`)
5. CLI Command Interface (`synaptik.cli.commands.clinical`)
6. REST API Endpoints (`synaptik.ui.server`)

Principal Investigator & Author:
Francesco Nobile (ORCID: 0009-0001-9477-3332)
Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School;
University of Pavia.
"""

from __future__ import annotations

import argparse
import io
import json
import math
import sys
import tempfile
from pathlib import Path

from synaptik.clinical.allometry import (
    CSF_VOLUMES,
    AllometryResult,
    CSFCompartment,
    calculate_allometric_dose,
)
from synaptik.clinical.escalation import (
    BOINEscalationEngine,
    CohortOutcome,
    EscalationAction,
    EscalationDecision,
    EscalationDesign,
    MTPIEngine,
    ThreePlusThreeEngine,
    betai,
    calculate_beta_overdose_prob,
)
from synaptik.clinical.pkpd import (
    PKPDExposure,
    PKPDKinetics,
    PKPDPharmacodynamics,
    PKPDSimulationResult,
    simulate_intrathecal_pkpd,
)
from synaptik.clinical.protocol import (
    BiomarkerMilestoneResult,
    ClinicalProtocol,
    DLTCriteria,
    DLTEvaluationResult,
    DoseTier,
    evaluate_bedside_dlt,
    evaluate_nfl_milestone,
    generate_clinical_protocol,
)
from synaptik.cli.commands.clinical import (
    cmd_clinical_allometry,
    cmd_clinical_escalate,
    cmd_clinical_pkpd,
    cmd_clinical_protocol,
)
from synaptik.ui.server import create_app
from synaptik.store import Store


def check(desc: str, cond: bool) -> None:
    """Assertion check with descriptive output for test reporting."""
    if cond:
        print(f"  PASS  {desc}")
    else:
        print(f"  FAIL  {desc}")
        raise AssertionError(desc)


# ===========================================================================
# 1. Allometric Scaling Tests
# ===========================================================================


def test_allometry_requires_explicit_acknowledgment():
    """Verify that calling calculate_allometric_dose without acknowledgment raises ValueError."""
    try:
        calculate_allometric_dose(noael_mg=30.0)
        check("Missing acknowledge_not_clinical raises ValueError", False)
    except ValueError as e:
        check("Missing acknowledge_not_clinical raises ValueError", "RESEARCH PLANNING ONLY" in str(e))


def test_allometry_adult_scaling():
    """Verify standard NHP-to-adult CSF volume scaling and MRSD calculation."""
    res = calculate_allometric_dose(
        noael_mg=30.0,
        nhp_volume_ml=15.0,
        human_volume_ml=140.0,
        safety_factor=10.0,
        acknowledge_not_clinical=True,
    )
    check("Allometry returns AllometryResult instance", isinstance(res, AllometryResult))
    check("Adult scaling ratio is exactly 140/15 ~= 9.333", abs(res.scaling_ratio - 9.3333333) < 1e-4)
    check("Adult HED is 280.0 mg", abs(res.hed_mg - 280.0) < 1e-4)
    check("Adult MRSD is 28.0 mg (SF=10)", abs(res.mrsd_mg - 28.0) < 1e-4)
    check("Result dictionary contains author attribution", "Francesco Nobile" in res.author)
    check("Result carries research planning disclaimer", "RESEARCH PLANNING ONLY" in res.disclaimer)
    check("is_clinical_recommendation is False", res.is_clinical_recommendation is False)


def test_allometry_pediatric_scaling():
    """Verify pediatric infant (65 mL) and child (100 mL) CSF scaling."""
    res_infant = calculate_allometric_dose(
        noael_mg=30.0,
        nhp_volume_ml=15.0,
        human_volume_ml=65.0,
        safety_factor=10.0,
        cohort_name="pediatric_infant",
        acknowledge_not_clinical=True,
    )
    check("Infant scaling ratio is 65/15 ~= 4.333", abs(res_infant.scaling_ratio - 4.3333333) < 1e-4)
    check("Infant HED is 130.0 mg", abs(res_infant.hed_mg - 130.0) < 1e-4)
    check("Infant MRSD is 13.0 mg", abs(res_infant.mrsd_mg - 13.0) < 1e-4)

    res_child = calculate_allometric_dose(
        noael_mg=30.0,
        nhp_volume_ml=15.0,
        human_volume_ml=100.0,
        safety_factor=10.0,
        cohort_name="pediatric_child",
        acknowledge_not_clinical=True,
    )
    check("Child scaling ratio is 100/15 ~= 6.667", abs(res_child.scaling_ratio - 6.6666667) < 1e-4)
    check("Child HED is 200.0 mg", abs(res_child.hed_mg - 200.0) < 1e-4)
    check("Child MRSD is 20.0 mg", abs(res_child.mrsd_mg - 20.0) < 1e-4)

    check("Monotonicity across age groups: Infant MRSD < Child MRSD < Adult MRSD",
          res_infant.mrsd_mg < res_child.mrsd_mg < 28.0)


def test_allometry_weight_normalized_conversion():
    """Verify conversion from NHP body-weight dose (mg/kg) to absolute dose."""
    res = calculate_allometric_dose(
        noael_mg_kg=10.0,
        nhp_weight_kg=3.0,
        nhp_volume_ml=15.0,
        human_volume_ml=140.0,
        safety_factor=10.0,
        acknowledge_not_clinical=True,
    )
    check("Body-weight conversion: 10 mg/kg * 3 kg = 30 mg total", abs(res.noael_nhp_mg - 30.0) < 1e-4)
    check("Adult MRSD is 28.0 mg", abs(res.mrsd_mg - 28.0) < 1e-4)


def test_allometry_safety_factor_variation():
    """Verify configurable safety factor from 6.0 to 10.0."""
    res_sf6 = calculate_allometric_dose(
        noael_mg=30.0,
        human_volume_ml=140.0,
        safety_factor=6.0,
        acknowledge_not_clinical=True,
    )
    check("SF=6.0 yields MRSD = 46.667 mg", abs(res_sf6.mrsd_mg - 46.666667) < 1e-4)
    res_sf10 = calculate_allometric_dose(
        noael_mg=30.0,
        human_volume_ml=140.0,
        safety_factor=10.0,
        acknowledge_not_clinical=True,
    )
    check("Higher safety factor strictly reduces starting dose", res_sf10.mrsd_mg < res_sf6.mrsd_mg)


def test_allometry_validation_and_errors():
    """Verify input validation and error raising for degenerate parameters."""
    try:
        calculate_allometric_dose(noael_mg=0.0, acknowledge_not_clinical=True)
        check("Zero NOAEL raises ValueError", False)
    except ValueError:
        check("Zero NOAEL raises ValueError", True)

    try:
        calculate_allometric_dose(noael_mg=-15.0, acknowledge_not_clinical=True)
        check("Negative NOAEL raises ValueError", False)
    except ValueError:
        check("Negative NOAEL raises ValueError", True)

    try:
        calculate_allometric_dose(noael_mg=30.0, nhp_volume_ml=0.0, acknowledge_not_clinical=True)
        check("Zero NHP volume raises ValueError", False)
    except ValueError:
        check("Zero NHP volume raises ValueError", True)

    try:
        calculate_allometric_dose(noael_mg=30.0, human_volume_ml=-5.0, acknowledge_not_clinical=True)
        check("Negative human volume raises ValueError", False)
    except ValueError:
        check("Negative human volume raises ValueError", True)

    try:
        calculate_allometric_dose(noael_mg=30.0, safety_factor=0.0, acknowledge_not_clinical=True)
        check("Zero safety factor raises ValueError", False)
    except ValueError:
        check("Zero safety factor raises ValueError", True)

    try:
        calculate_allometric_dose(noael_mg=None, noael_mg_kg=None, acknowledge_not_clinical=True)
        check("Missing both NOAEL inputs raises ValueError", False)
    except ValueError:
        check("Missing both NOAEL inputs raises ValueError", True)


# ===========================================================================
# 2. Intrathecal PK/PD Simulator Tests
# ===========================================================================


def test_pkpd_bolus_vs_catheter_peak_concentration():
    """Verify catheter infusion flattens peak CSF concentration while preserving tissue exposure."""
    res_bolus = simulate_intrathecal_pkpd(dose_mg=50.0, route="bolus", duration_days=360)
    res_cath = simulate_intrathecal_pkpd(dose_mg=50.0, route="catheter", infusion_hours=24.0, duration_days=360)

    cmax_b = res_bolus.exposure.cmax_csf_ug_ml
    cmax_c = res_cath.exposure.cmax_csf_ug_ml

    check("Bolus peak CSF concentration is > 300 ug/mL (approx 357 ug/mL)", cmax_b > 300.0)
    check("Catheter peak CSF concentration is < 145 ug/mL", cmax_c < 145.0)
    check("Bolus Cmax is more than 2x catheter Cmax", cmax_b > cmax_c * 2.0)

    pct_red = ((cmax_b - cmax_c) / cmax_b) * 100.0
    check("Catheter peak reduction is >= 60%", pct_red >= 60.0)

    auc_b = res_bolus.exposure.auc_tissue_ug_day_g
    auc_c = res_cath.exposure.auc_tissue_ug_day_g
    auc_diff_pct = abs(auc_b - auc_c) / auc_b
    check("Tissue AUC between bolus and catheter is bioequivalent within 15%", auc_diff_pct < 0.15)


def test_pkpd_half_life_clearance_kinetics():
    """Verify terminal clearance kinetics follow nominal half-lives."""
    res = simulate_intrathecal_pkpd(dose_mg=50.0, csf_half_life_days=24.0, tissue_half_life_days=135.0)
    check("CSF elimination rate constant matches ln(2)/24", abs(res.kinetics.k_el_csf - (math.log(2.0)/24.0)) < 1e-5)
    check("Tissue elimination rate constant matches ln(2)/135", abs(res.kinetics.k_el_tissue - (math.log(2.0)/135.0)) < 1e-5)


def test_pkpd_mass_conservation_and_ode_stability():
    """Verify mass conservation and numerical stability across ODE solver."""
    res = simulate_intrathecal_pkpd(dose_mg=50.0, route="bolus", duration_days=360)
    all_non_negative = all(pt.a_csf_mg >= 0.0 and pt.a_tissue_mg >= 0.0 for pt in res.trajectory)
    check("Compartment drug masses are strictly non-negative", all_non_negative)

    all_conserved = all((pt.a_csf_mg + pt.a_tissue_mg) <= (50.0 + 1e-4) for pt in res.trajectory)
    check("Total active drug mass does not exceed administered 50 mg", all_conserved)

    p0 = res.trajectory[0]
    p_last = res.trajectory[-1]
    check("Terminal drug mass at Day 360 is substantially cleared (< 5 mg remaining)", (p_last.a_csf_mg + p_last.a_tissue_mg) < 5.0)


def test_pkpd_emax_hill_pharmacodynamic_coupling():
    """Verify target tissue engagement follows Hill pharmacodynamic model."""
    res = simulate_intrathecal_pkpd(dose_mg=50.0, duration_days=360)
    # Check that peak tissue engagement reaches near saturation (>95%)
    peak_eff = max(pt.target_engagement_pct for pt in res.trajectory)
    check("Therapeutic dose achieves near-maximal target engagement (>95%)", peak_eff > 95.0)

    # Check baseline t=0
    check("Baseline target engagement is 0% at t=0", res.trajectory[0].target_engagement_pct == 0.0)


def test_pkpd_day90_nfl_milestone_trajectory():
    """Verify evaluation of Day 90 plasma NfL surrogate milestone (>30% reduction)."""
    # Therapeutic 50 mg dose
    res_50 = simulate_intrathecal_pkpd(dose_mg=50.0, duration_days=360)
    check("50 mg dose meets Day 90 NfL milestone (>30% reduction)", res_50.day90_milestone_met is True)
    check("50 mg dose predicted Day 90 reduction exceeds 40%", res_50.day90_nfl_reduction_pct > 40.0)

    # Sub-therapeutic 0.5 mg dose
    res_sub = simulate_intrathecal_pkpd(dose_mg=0.5, duration_days=360)
    check("0.5 mg sub-therapeutic dose fails Day 90 NfL milestone", res_sub.day90_milestone_met is False)
    check("0.5 mg dose predicted Day 90 reduction is < 30%", res_sub.day90_nfl_reduction_pct < 30.0)


def test_pkpd_validation_and_errors():
    """Verify input argument validation on PK/PD simulator."""
    try:
        simulate_intrathecal_pkpd(dose_mg=0.0)
        check("Zero dose raises ValueError", False)
    except ValueError:
        check("Zero dose raises ValueError", True)

    try:
        simulate_intrathecal_pkpd(dose_mg=50.0, route="subcutaneous")
        check("Invalid route raises ValueError", False)
    except ValueError:
        check("Invalid route raises ValueError", True)

    try:
        simulate_intrathecal_pkpd(dose_mg=50.0, csf_half_life_days=0.0)
        check("Zero CSF half-life raises ValueError", False)
    except ValueError:
        check("Zero CSF half-life raises ValueError", True)

    try:
        simulate_intrathecal_pkpd(dose_mg=50.0, duration_days=-30)
        check("Negative duration raises ValueError", False)
    except ValueError:
        check("Negative duration raises ValueError", True)


# ===========================================================================
# 3. Adaptive Escalation Engine Tests
# ===========================================================================


def test_escalation_traditional_3plus3_decision_tree():
    """Verify all branching pathways of traditional 3+3 escalation."""
    engine = ThreePlusThreeEngine(dose_levels=[10.0, 25.0, 50.0, 100.0, 150.0])

    # 1. Level 1, 0/3 DLTs -> ESCALATE
    d1 = engine.evaluate_cohort(dose_level=1, n=3, dlts=0)
    check("3+3: 0/3 DLTs triggers ESCALATE to Level 2", d1.action == EscalationAction.ESCALATE and d1.next_dose_level == 2)

    # 2. Level 1, 1/3 DLTs -> RETAIN (expand to N=6)
    d2 = engine.evaluate_cohort(dose_level=1, n=3, dlts=1)
    check("3+3: 1/3 DLTs triggers RETAIN (expand cohort)", d2.action == EscalationAction.RETAIN and d2.next_dose_level == 1)

    # 3. Level 1 expanded, 1/6 DLTs -> ESCALATE
    d3 = engine.evaluate_cohort(dose_level=1, n=6, dlts=1)
    check("3+3: 1/6 DLTs triggers ESCALATE to Level 2", d3.action == EscalationAction.ESCALATE and d3.next_dose_level == 2)

    # 4. Level 1 expanded, 2/6 DLTs -> STOP_TRIAL
    d4 = engine.evaluate_cohort(dose_level=1, n=6, dlts=2)
    check("3+3: 2/6 DLTs at Level 1 halts trial (MTD not identified)", d4.action == EscalationAction.STOP_TRIAL and d4.mtd_declared_dose_level is None)

    # 5. Level 1, 2/3 DLTs -> STOP_TRIAL
    d5 = engine.evaluate_cohort(dose_level=1, n=3, dlts=2)
    check("3+3: 2/3 DLTs at Level 1 halts trial", d5.action == EscalationAction.STOP_TRIAL)

    # 6. Level 2, 2/3 DLTs -> STOP_TRIAL / declare Level 1 MTD
    d6 = engine.evaluate_cohort(dose_level=2, n=3, dlts=2)
    check("3+3: 2/3 DLTs at Level 2 declares Level 1 as MTD", d6.action == EscalationAction.STOP_TRIAL and d6.mtd_declared_dose_level == 1)


def test_escalation_boin_interval_boundaries():
    """Verify analytical calculation of BOIN boundaries lambda_e and lambda_d."""
    boin25 = BOINEscalationEngine(target_phi=0.25)
    check("BOIN target phi 0.25 lambda_e ~= 0.1968", abs(boin25.lambda_e - 0.1968) < 0.001)
    check("BOIN target phi 0.25 lambda_d ~= 0.2984", abs(boin25.lambda_d - 0.2984) < 0.001)
    check("BOIN ordering holds: lambda_e < phi < lambda_d", boin25.lambda_e < boin25.target_phi < boin25.lambda_d)

    boin30 = BOINEscalationEngine(target_phi=0.30)
    check("BOIN target phi 0.30 lambda_e ~= 0.2365", abs(boin30.lambda_e - 0.2365) < 0.001)
    check("BOIN target phi 0.30 lambda_d ~= 0.3585", abs(boin30.lambda_d - 0.3585) < 0.001)


def test_escalation_boin_decision_rules():
    """Verify BOIN decision logic based on observed DLT rates."""
    boin = BOINEscalationEngine(target_phi=0.25, dose_levels=[10.0, 25.0, 50.0, 100.0, 150.0])

    # 0/3 DLTs (p=0.0 <= 0.197) -> ESCALATE
    d1 = boin.evaluate_cohort(dose_level=1, n=3, dlts=0)
    check("BOIN: 0/3 DLTs recommends ESCALATE to Level 2", d1.action == EscalationAction.ESCALATE and d1.next_dose_level == 2)

    # 1/4 DLTs (p=0.25 in (0.197, 0.298]) -> RETAIN
    d2 = boin.evaluate_cohort(dose_level=2, n=4, dlts=1)
    check("BOIN: 1/4 DLTs recommends RETAIN Level 2", d2.action == EscalationAction.RETAIN and d2.next_dose_level == 2)

    # 1/3 DLTs (p=0.333 > 0.298) -> DE_ESCALATE
    d3 = boin.evaluate_cohort(dose_level=2, n=3, dlts=1)
    check("BOIN: 1/3 DLTs recommends DE_ESCALATE to Level 1", d3.action == EscalationAction.DE_ESCALATE and d3.next_dose_level == 1)

    # Maximum dose level (Level 5) with 0/3 DLTs -> RETAIN
    d_max = boin.evaluate_cohort(dose_level=5, n=3, dlts=0)
    check("BOIN at max dose with 0 DLTs retains Level 5", d_max.action == EscalationAction.RETAIN and d_max.next_dose_level == 5)


def test_escalation_boin_beta_posterior_safety_stopping():
    """Verify Bayesian conjugate posterior overdose safety stopping rule."""
    boin = BOINEscalationEngine(target_phi=0.25, dose_levels=[10.0, 25.0, 50.0, 100.0, 150.0])

    # 3/3 DLTs at Level 1: Beta(4, 1) posterior P(p > 0.25) = 1 - 0.25^4 = 0.9961 > 0.95
    d_toxic = boin.evaluate_cohort(dose_level=1, n=3, dlts=3)
    check("BOIN 3/3 DLTs at Level 1 triggers STOP_TRIAL", d_toxic.action == EscalationAction.STOP_TRIAL)
    check("BOIN 3/3 DLTs sets trial_stopped = True", d_toxic.trial_stopped is True)
    check("BOIN 3/3 DLTs posterior overdose prob > 0.99", d_toxic.posterior_overdose_prob > 0.99)
    check("BOIN 3/3 DLTs eliminates dose Level 1", d_toxic.is_dose_eliminated is True)

    # 0/3 DLTs at Level 1: Beta(1, 4) posterior P(p > 0.25) = 0.75^4 = 0.3164 < 0.95
    d_safe = boin.evaluate_cohort(dose_level=1, n=3, dlts=0)
    check("BOIN 0/3 DLTs does NOT eliminate dose", d_safe.is_dose_eliminated is False)
    check("BOIN 0/3 DLTs posterior overdose prob < 0.35", d_safe.posterior_overdose_prob < 0.35)


def test_escalation_mtpi_decision_intervals():
    """Verify mTPI Unit Probability Mass (UPM) calculations and transitions."""
    mtpi = MTPIEngine(target_phi=0.25, epsilon=0.05, dose_levels=[10.0, 25.0, 50.0, 100.0, 150.0])

    # 0/3 DLTs -> UPM_E highest -> ESCALATE
    d_esc = mtpi.evaluate_cohort(dose_level=1, n=3, dlts=0)
    check("mTPI: 0/3 DLTs recommends ESCALATE", d_esc.action == EscalationAction.ESCALATE and d_esc.next_dose_level == 2)

    # 1/3 DLTs -> UPM_T highest -> RETAIN
    d_ret = mtpi.evaluate_cohort(dose_level=2, n=3, dlts=1)
    check("mTPI: 1/3 DLTs recommends RETAIN", d_ret.action == EscalationAction.RETAIN and d_ret.next_dose_level == 2)

    # 3/3 DLTs -> Safety stopping / STOP_TRIAL
    d_stop = mtpi.evaluate_cohort(dose_level=1, n=3, dlts=3)
    check("mTPI: 3/3 DLTs triggers safety stopping", d_stop.action == EscalationAction.STOP_TRIAL)


def test_escalation_validation_and_errors():
    """Verify boundary input validation for escalation engines."""
    try:
        BOINEscalationEngine(target_phi=0.0)
        check("BOIN zero target phi raises ValueError", False)
    except ValueError:
        check("BOIN zero target phi raises ValueError", True)

    try:
        BOINEscalationEngine(target_phi=1.0)
        check("BOIN 1.0 target phi raises ValueError", False)
    except ValueError:
        check("BOIN 1.0 target phi raises ValueError", True)

    boin = BOINEscalationEngine(target_phi=0.25)
    try:
        boin.evaluate_cohort(dose_level=-1, n=3, dlts=0)
        check("Negative dose level raises ValueError", False)
    except ValueError:
        check("Negative dose level raises ValueError", True)

    try:
        boin.evaluate_cohort(dose_level=1, n=3, dlts=4)
        check("DLTs > n raises ValueError", False)
    except ValueError:
        check("DLTs > n raises ValueError", True)


# ===========================================================================
# 4. Bedside Protocol & DLT Stopping Tests
# ===========================================================================


def test_protocol_bedside_dlt_criteria_evaluation():
    """Verify deterministic bedside adverse event evaluation against 4 DLT criteria."""
    # 1. CSF Pleocytosis boundary: > 50 cells/uL
    res_50 = evaluate_bedside_dlt(csf_wbc_per_ul=50.0)
    check("WBC = 50.0 does NOT trigger DLT (requires >50)", res_50.is_dlt is False)
    res_51 = evaluate_bedside_dlt(csf_wbc_per_ul=51.0)
    check("WBC = 51.0 DOES trigger DLT", res_51.is_dlt is True and "CSF Pleocytosis" in res_51.triggered_criteria[0])

    # 2. Chemical Arachnoiditis
    res_arach = evaluate_bedside_dlt(chemical_arachnoiditis_present=True)
    check("Chemical arachnoiditis triggers DLT", res_arach.is_dlt is True)

    # 3. Sensory Radiculopathy: Grade >= 2 and duration > 7 days
    res_rad1 = evaluate_bedside_dlt(radiculopathy_grade=1, radiculopathy_duration_days=10.0)
    check("Grade 1 radiculopathy does NOT trigger DLT (requires Grade >= 2)", res_rad1.is_dlt is False)
    res_rad2_short = evaluate_bedside_dlt(radiculopathy_grade=2, radiculopathy_duration_days=5.0)
    check("Grade 2 radiculopathy for 5 days does NOT trigger DLT (requires >7 days)", res_rad2_short.is_dlt is False)
    res_rad2_long = evaluate_bedside_dlt(radiculopathy_grade=2, radiculopathy_duration_days=8.0)
    check("Grade 2 radiculopathy for 8 days DOES trigger DLT", res_rad2_long.is_dlt is True)

    # 4. Severe Thrombocytopenia: < 50,000 or > 50% drop
    res_plt_low = evaluate_bedside_dlt(platelet_count_per_ul=45000.0)
    check("Platelets < 50,000 triggers DLT", res_plt_low.is_dlt is True)
    res_plt_drop = evaluate_bedside_dlt(platelet_count_per_ul=60000.0, baseline_platelet_count=140000.0)
    check("Platelet drop > 50% (57.1% drop) triggers DLT", res_plt_drop.is_dlt is True)
    res_plt_ok = evaluate_bedside_dlt(platelet_count_per_ul=75000.0, baseline_platelet_count=100000.0)
    check("Platelets >= 50,000 and <= 50% drop does NOT trigger DLT", res_plt_ok.is_dlt is False)


def test_protocol_day90_nfl_biomarker_milestone():
    """Verify Day 90 fluid biomarker response milestone (>30.0% reduction)."""
    # 40% reduction -> Met
    m40 = evaluate_nfl_milestone(baseline_nfl_pg_ml=70.0, day90_nfl_pg_ml=42.0)
    check("40% NfL reduction meets Day 90 milestone", m40.milestone_met is True)
    check("40% reduction percentage is 40.0%", abs(m40.reduction_pct - 40.0) < 1e-4)

    # 20% reduction -> Not met
    m20 = evaluate_nfl_milestone(baseline_nfl_pg_ml=70.0, day90_nfl_pg_ml=56.0)
    check("20% NfL reduction fails Day 90 milestone", m20.milestone_met is False)

    # Exact 30.0% reduction -> Not met (strictly >30%)
    m30 = evaluate_nfl_milestone(baseline_nfl_pg_ml=70.0, day90_nfl_pg_ml=49.0)
    check("Exact 30.0% NfL reduction fails Day 90 milestone (requires >30%)", m30.milestone_met is False)


def test_protocol_document_generation():
    """Verify automated clinical trial protocol Markdown document generation."""
    proto = generate_clinical_protocol(
        target_gene="STMN2",
        disease="als",
        noael_nhp_mg=30.0,
        design="boin",
    )
    check("Protocol object is ClinicalProtocol instance", isinstance(proto, ClinicalProtocol))
    check("Protocol identifier is SYN-FIH-STMN2-001", proto.protocol_id == "SYN-FIH-STMN2-001")
    check("Protocol starting dose MRSD is 28.0 mg", abs(proto.mrsd_mg - 28.0) < 1e-4)
    check("Protocol Markdown contains title", "Phase 1 FIH Intrathecal Dose Escalation Protocol: STMN2" in proto.markdown)
    check("Protocol Markdown contains MRSD 28.0 mg", "28.0 mg" in proto.markdown)
    check("Protocol Markdown cites Francesco Nobile", "Francesco Nobile" in proto.markdown)
    check("Protocol Markdown contains verified ORCID", "0009-0001-9477-3332" in proto.markdown)


# ===========================================================================
# 5. CLI Command Tests
# ===========================================================================


def test_cli_clinical_allometry_action():
    """Verify CLI execution of `synaptik clinical allometry`."""
    args = argparse.Namespace(
        noael_mg=30.0,
        noael_mg_kg=None,
        nhp_weight=3.0,
        nhp_volume=15.0,
        adult_volume=140.0,
        pediatric_volume=65.0,
        safety_factor=10.0,
        json=True,
    )
    old_stdout = sys.stdout
    sys.stdout = io.StringIO()
    try:
        ret = cmd_clinical_allometry(args)
        out = sys.stdout.getvalue()
    finally:
        sys.stdout = old_stdout

    check("CLI allometry returns exit code 0", ret == 0)
    data = json.loads(out)
    check("CLI allometry JSON contains adult MRSD 28.0 mg", data["adult"]["mrsd_mg"] == 28.0)
    check("CLI allometry JSON contains author attribution", "Francesco Nobile" in data["author"])


def test_cli_clinical_pkpd_action():
    """Verify CLI execution of `synaptik clinical pkpd`."""
    args = argparse.Namespace(
        dose_mg=50.0,
        route="bolus",
        csf_half_life=24.0,
        tissue_half_life=135.0,
        days=360,
        infusion_hours=24.0,
        json=True,
    )
    old_stdout = sys.stdout
    sys.stdout = io.StringIO()
    try:
        ret = cmd_clinical_pkpd(args)
        out = sys.stdout.getvalue()
    finally:
        sys.stdout = old_stdout

    check("CLI pkpd returns exit code 0", ret == 0)
    data = json.loads(out)
    check("CLI pkpd JSON contains exposure metrics", "exposure" in data and "pharmacodynamics" in data)
    check("CLI pkpd Day 90 milestone is met", data["day90_milestone_met"] is True)


def test_cli_clinical_escalate_action():
    """Verify CLI execution of `synaptik clinical escalate`."""
    args = argparse.Namespace(
        design="boin",
        target_phi=0.25,
        doses=[10.0, 25.0, 50.0, 100.0, 150.0],
        history='[{"dose_level": 1, "n": 3, "dlts": 0}]',
        json=True,
    )
    old_stdout = sys.stdout
    sys.stdout = io.StringIO()
    try:
        ret = cmd_clinical_escalate(args)
        out = sys.stdout.getvalue()
    finally:
        sys.stdout = old_stdout

    check("CLI escalate returns exit code 0", ret == 0)
    data = json.loads(out)
    check("CLI escalate recommends escalate to Level 2", data["recommendation"] == "escalate" and data["recommended_dose_level"] == 2)


def test_cli_clinical_protocol_action():
    """Verify CLI execution of `synaptik clinical protocol`."""
    args = argparse.Namespace(
        gene="STMN2",
        disease="als",
        noael_mg=30.0,
        design="boin",
        target_phi=0.25,
        doses=[10.0, 25.0, 50.0, 100.0, 150.0],
        output=None,
        json=True,
    )
    old_stdout = sys.stdout
    sys.stdout = io.StringIO()
    try:
        ret = cmd_clinical_protocol(args)
        out = sys.stdout.getvalue()
    finally:
        sys.stdout = old_stdout

    check("CLI protocol returns exit code 0", ret == 0)
    data = json.loads(out)
    check("CLI protocol JSON contains protocol_id", "protocol_id" in data)
    check("CLI protocol JSON contains markdown", "markdown" in data and len(data["markdown"]) > 100)


# ===========================================================================
# 6. REST API Endpoint Tests
# ===========================================================================


def test_api_clinical_endpoints():
    """Verify all four `/api/clinical/*` REST endpoints on the Flask test client."""
    td = tempfile.TemporaryDirectory()
    db_path = str(Path(td.name) / "test_api_clinical.db")
    Store(db_path).close()

    app = create_app(db_path)
    client = app.test_client()

    # 1. /api/clinical/allometry
    res1 = client.post("/api/clinical/allometry", json={
        "noael_mg": 30.0,
        "nhp_volume_ml": 15.0,
        "adult_volume_ml": 140.0,
        "safety_factor": 10.0,
    })
    check("API POST /api/clinical/allometry returns 200 OK", res1.status_code == 200)
    d1 = res1.get_json()
    check("API allometry adult MRSD is 28.0 mg", d1["adult"]["mrsd_mg"] == 28.0)
    check("API allometry cites Francesco Nobile", "Francesco Nobile" in d1["author"])

    # Error path: negative dose
    res1_err = client.post("/api/clinical/allometry", json={"noael_mg": -5.0})
    check("API allometry negative dose returns 400 Bad Request", res1_err.status_code == 400)

    # 2. /api/clinical/pkpd
    res2 = client.post("/api/clinical/pkpd", json={
        "dose_mg": 50.0,
        "route": "bolus",
        "duration_days": 360,
    })
    check("API POST /api/clinical/pkpd returns 200 OK", res2.status_code == 200)
    d2 = res2.get_json()
    check("API pkpd peak CSF conc > 300 ug/mL", d2["exposure"]["cmax_csf_ug_ml"] > 300.0)
    check("API pkpd Day 90 surrogate milestone is met", d2["pharmacodynamics"]["clinical_surrogate_met"] is True)

    # Error path: invalid route
    res2_err = client.post("/api/clinical/pkpd", json={"dose_mg": 50.0, "route": "intramuscular"})
    check("API pkpd invalid route returns 400 Bad Request", res2_err.status_code == 400)

    # 3. /api/clinical/escalate
    res3 = client.post("/api/clinical/escalate", json={
        "design": "boin",
        "target_phi": 0.25,
        "current_level": 1,
        "history": [{"dose_level": 1, "n": 3, "dlts": 0}],
    })
    check("API POST /api/clinical/escalate returns 200 OK", res3.status_code == 200)
    d3 = res3.get_json()
    check("API escalate recommendation is escalate", d3["recommendation"] == "escalate")
    check("API escalate recommended dose level is 2", d3["recommended_dose_level"] == 2)

    # Error path: invalid design
    res3_err = client.post("/api/clinical/escalate", json={"design": "unknown_design"})
    check("API escalate invalid design returns 400 Bad Request", res3_err.status_code == 400)

    # 4. /api/clinical/protocol
    res4 = client.post("/api/clinical/protocol", json={
        "gene": "STMN2",
        "disease": "als",
        "noael_mg": 30.0,
        "design": "boin",
    })
    check("API POST /api/clinical/protocol returns 200 OK", res4.status_code == 200)
    d4 = res4.get_json()
    check("API protocol returns protocol_id SYN-FIH-STMN2-001", d4["protocol_id"] == "SYN-FIH-STMN2-001")
    check("API protocol returns markdown text", len(d4["markdown"]) > 100)


# ===========================================================================
# Top-Level Test Runner
# ===========================================================================


def run_tests():
    print("=" * 78)
    print("SYNAPTIK Clinical Protocol & Dose Escalation Engine Test Suite")
    print("Author: Francesco Nobile (ORCID: 0009-0001-9477-3332)")
    print("=" * 78)

    print("\n[1. Allometric Scaling Tests]")
    test_allometry_adult_scaling()
    test_allometry_pediatric_scaling()
    test_allometry_weight_normalized_conversion()
    test_allometry_safety_factor_variation()
    test_allometry_validation_and_errors()

    print("\n[2. Intrathecal PK/PD Simulator Tests]")
    test_pkpd_bolus_vs_catheter_peak_concentration()
    test_pkpd_half_life_clearance_kinetics()
    test_pkpd_mass_conservation_and_ode_stability()
    test_pkpd_emax_hill_pharmacodynamic_coupling()
    test_pkpd_day90_nfl_milestone_trajectory()
    test_pkpd_validation_and_errors()

    print("\n[3. Adaptive Dose Escalation Engine Tests]")
    test_escalation_traditional_3plus3_decision_tree()
    test_escalation_boin_interval_boundaries()
    test_escalation_boin_decision_rules()
    test_escalation_boin_beta_posterior_safety_stopping()
    test_escalation_mtpi_decision_intervals()
    test_escalation_validation_and_errors()

    print("\n[4. Bedside Protocol & DLT Stopping Tests]")
    test_protocol_bedside_dlt_criteria_evaluation()
    test_protocol_day90_nfl_biomarker_milestone()
    test_protocol_document_generation()

    print("\n[5. CLI Command Tests]")
    test_cli_clinical_allometry_action()
    test_cli_clinical_pkpd_action()
    test_cli_clinical_escalate_action()
    test_cli_clinical_protocol_action()

    print("\n[6. REST API Endpoint Tests]")
    test_api_clinical_endpoints()

    print("\n" + "=" * 78)
    print("ALL CLINICAL PROTOCOL & DOSE ESCALATION TESTS PASSED (100%)")
    print("=" * 78)


if __name__ == "__main__":
    run_tests()
