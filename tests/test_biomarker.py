"""Tests for the Translational Fluid Biomarker & Patient Stratification Panel (`synaptik biomarker`)."""

import unittest
from synaptik.biomarker.nfl import (
    BiomarkerAnalyzer,
    BiomarkerSubtype,
    PatientProgressionTier,
    compute_sample_size_two_arm,
)
from synaptik.workflows import evaluate_biomarker_trajectory


def check(desc, cond):
    if cond:
        print(f"  PASS  {desc}")
    else:
        print(f"  FAIL  {desc}")
        raise AssertionError(desc)


def run_tests():
    print("=" * 70)
    print("SYNAPTIK fluid biomarker & patient stratification test suite")
    print("=" * 70)

    analyzer = BiomarkerAnalyzer("als")

    # ------------------------------------------------------------------
    # 1. Longitudinal pharmacodynamics & clearance kinetics
    # ------------------------------------------------------------------
    print("\n[biomarker: pharmacodynamics & clearance kinetics]")
    res_stmn2 = analyzer.evaluate_biomarker_response("STMN2", in_vitro_rescue_pct=75.0)

    check("baseline plasma NfL is positive (got " + str(res_stmn2.baseline_plasma_nfl_pg_ml) + ")",
          res_stmn2.baseline_plasma_nfl_pg_ml > 30.0)
    check("maximal reduction matches scaling beta (~54%)",
          abs(res_stmn2.maximal_reduction_pct - (75.0 * 0.72)) < 0.2)

    traj = res_stmn2.longitudinal_trajectory
    check("week 0 has 0% reduction", traj[0].percent_reduction == 0.0)
    check("week 12 shows significant reduction (> 35%)", res_stmn2.predicted_week12_nfl < res_stmn2.baseline_plasma_nfl_pg_ml * 0.65)
    check("week 24 reaches near-plateau pharmacodynamic suppression",
          res_stmn2.predicted_week24_nfl < res_stmn2.predicted_week12_nfl)
    check("week 48 achieves maximal suppression",
          res_stmn2.predicted_week48_nfl <= res_stmn2.predicted_week24_nfl)

    # Monotonicity check
    is_monotonic = all(traj[i].predicted_nfl_pg_ml >= traj[i+1].predicted_nfl_pg_ml for i in range(len(traj)-1))
    check("longitudinal NfL trajectory decreases monotonically", is_monotonic)

    # ------------------------------------------------------------------
    # 2. Patient cohort stratification & trial sample size
    # ------------------------------------------------------------------
    print("\n[biomarker: patient cohort stratification & trial power]")
    strat = res_stmn2.stratification
    check("unstratified trial requires larger sample size than stratified cohort",
          strat.unstratified_n_per_arm > strat.stratified_n_per_arm)
    check("sample size reduction is >= 60% (got " + str(strat.sample_size_reduction_pct) + "%)",
          strat.sample_size_reduction_pct >= 60.0)
    check("trial feasibility verdict rendered", "FEASIBILITY" in strat.feasibility_verdict)

    # Fast progressor tier detection
    res_fast = analyzer.evaluate_biomarker_response("STMN2", in_vitro_rescue_pct=70.0, baseline_plasma_nfl=95.0)
    check("plasma NfL >= 75 pg/mL classified as fast progressor",
          res_fast.stratification.tier == PatientProgressionTier.FAST_PROGRESSOR)
    check("fast progressor expected decline rate > 1.0 pts/month",
          res_fast.stratification.expected_decline_rate_pts_month > 1.0)

    # Slow progressor tier detection
    res_slow = analyzer.evaluate_biomarker_response("STMN2", in_vitro_rescue_pct=70.0, baseline_plasma_nfl=25.0)
    check("plasma NfL < 40 pg/mL classified as slow progressor",
          res_slow.stratification.tier == PatientProgressionTier.SLOW_PROGRESSOR)

    # Sample size formula edge cases
    n_std = compute_sample_size_two_arm(0.30, 0.32)
    check("sample size calculation returns reasonable cohort size (got " + str(n_std) + ")", 15 < n_std < 50)

    # ------------------------------------------------------------------
    # 3. Regulatory surrogate endpoints & indications
    # ------------------------------------------------------------------
    print("\n[biomarker: regulatory surrogate milestone verdicts]")
    check("high rescue triggers substantial surrogate efficacy verdict",
          "SUBSTANTIAL SURROGATE EFFICACY" in res_stmn2.clinical_surrogate_verdict)
    check("verdict cites clinical precedent (VALOR / tofersen)",
          "tofersen" in res_stmn2.clinical_surrogate_verdict)

    res_low = analyzer.evaluate_biomarker_response("TEST", in_vitro_rescue_pct=15.0)
    check("low rescue triggers low surrogate efficacy warning",
          "LOW SURROGATE EFFICACY" in res_low.clinical_surrogate_verdict)

    # Indication support: FTD and AD
    analyzer_ftd = BiomarkerAnalyzer("ftd")
    res_ftd = analyzer_ftd.evaluate_biomarker_response("GRN", in_vitro_rescue_pct=80.0)
    check("FTD disease indication correctly assigned", res_ftd.disease == "FTD")
    check("FTD baseline NfL evaluated", res_ftd.baseline_plasma_nfl_pg_ml > 20.0)

    # ------------------------------------------------------------------
    # 4. Workflows API integration
    # ------------------------------------------------------------------
    print("\n[biomarker: workflows API serialization]")
    api_res = evaluate_biomarker_trajectory(disease="als", target_gene="STMN2", in_vitro_rescue_pct=70.0)
    check("workflow returns target_gene", api_res["target_gene"] == "STMN2")
    check("workflow returns maximal_reduction_pct", "maximal_reduction_pct" in api_res)
    check("workflow returns stratification dictionary", isinstance(api_res["stratification"], dict))
    check("workflow returns longitudinal trajectory list", len(api_res["longitudinal_trajectory"]) == 8)
    check("workflow returns actionable recommendations", len(api_res["actionable_trial_recommendations"]) >= 2)

    print("\n" + "=" * 70)
    print("All checks passed.")


if __name__ == "__main__":
    run_tests()
