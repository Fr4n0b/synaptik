"""
Patient Natural History Digital Twin & Bayesian Counterfactual Control Test Suite.

Run: python -m tests.test_twin
"""

from __future__ import annotations

import sys

from synaptik.twin import (
    PatientBaseline,
    VisitObservation,
    compute_conditional_progression,
    evaluate_treatment_efficacy,
    generate_counterfactual_cohort,
    generate_twin_fda_dossier,
)

FAILURES: list[str] = []


def check(cond: bool, msg: str) -> None:
    print(f"  {'PASS' if cond else 'FAIL'}  {msg}")
    if not cond:
        FAILURES.append(msg)


def test_baseline_and_distributions() -> None:
    print("\n[twin: baseline ingestion & PRO-ACT/Answer ALS covariance]")
    b = PatientBaseline("ALS-001", 55.0, "BULBAR", 10.0, 40.0)
    check(b.onset_site == "bulbar", "onset site normalized to lowercase")
    check(abs(b.baseline_slope - (-0.80)) < 1e-4, "pre-slope auto-computed as (40-48)/10 = -0.80")

    b_limb = PatientBaseline("ALS-002", 55.0, "limb", 10.0, 40.0, baseline_slope=-0.92, baseline_plasma_nfl=58.0)
    b_bulbar = PatientBaseline("ALS-003", 55.0, "bulbar", 10.0, 40.0, baseline_slope=-0.92, baseline_plasma_nfl=58.0)
    c_limb = compute_conditional_progression(b_limb)
    c_bulbar = compute_conditional_progression(b_bulbar)
    check(c_bulbar.mean_slope < c_limb.mean_slope, "bulbar onset progresses faster than limb")
    check(abs(c_limb.std_slope - 0.280) < 1e-4, "conditional residual SD is 0.280 pts/mo")


def test_cohort_simulation_and_bayesian() -> None:
    print("\n[twin: 1,000-twin simulation & Bayesian deceleration posterior]")
    b = PatientBaseline("ALS-004", 56.4, "bulbar", 14.5, 38.0, baseline_slope=-1.20, baseline_plasma_nfl=72.5)
    cohort = generate_counterfactual_cohort(b, n_twins=1000, duration_months=12, random_seed=42)
    check(cohort.n_twins == 1000 and len(cohort.trajectories) == 1000, "1,000 twin trajectories simulated")
    check(len(cohort.median_trajectory) == 13, "12-month trajectory has 13 monthly points (0..12)")
    check(all(0.0 <= x <= 48.0 for traj in cohort.trajectories for x in traj), "all scores bounded in [0, 48]")

    visits_stable = [
        VisitObservation(0.0, 38.0),
        VisitObservation(3.0, 38.0),
        VisitObservation(6.0, 38.0),
        VisitObservation(12.0, 38.0),
    ]
    eff = evaluate_treatment_efficacy(cohort, visits_stable)
    check(eff.p_deceleration >= 0.975, f"stabilized patient yields P(Delta>0)={eff.p_deceleration:.4f} >= 0.975")
    check("Tier A" in eff.efficacy_tier, f"classified as {eff.efficacy_tier}")

    dossier = generate_twin_fda_dossier(b, cohort, efficacy=eff)
    check("COMPUTATIONAL DRAFT" in dossier, "dossier carries computational draft banner")
    check("Francesco Nobile" in dossier and "0009-0001-9477-3332" in dossier, "dossier cites Francesco Nobile & ORCID")

    from synaptik.twin.twin_generator import evaluate_twin_calibration
    cal = evaluate_twin_calibration(b, n_held_out=200, n_twins=400, duration_months=12, random_seed=77)
    check(0.72 <= cal["mean_empirical_coverage"] <= 0.88,
          f"P10-P90 interval achieves ~80% empirical coverage (got {cal['mean_empirical_coverage']:.3f})")


def main() -> int:
    print("=" * 70)
    print("SYNAPTIK Patient Digital Twin Test Suite")
    print("=" * 70)
    test_baseline_and_distributions()
    test_cohort_simulation_and_bayesian()
    print("\n" + "=" * 70)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S):")
        for f in FAILURES:
            print("  -", f)
        return 1
    print("All digital twin checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
