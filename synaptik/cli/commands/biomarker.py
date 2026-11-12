"""`synaptik biomarker`: fluid biomarker trajectory modeling and patient trial stratification."""

from __future__ import annotations

import json
from dataclasses import asdict

from ...biomarker.nfl import BiomarkerAnalyzer, BiomarkerSubtype
from ..main import Fmt, header


def cmd_biomarker(args) -> int:
    disease = getattr(args, "disease", "als") or "als"
    target = getattr(args, "target", "STMN2") or "STMN2"
    rescue_pct = float(getattr(args, "rescue_pct", 70.0) or 70.0)
    baseline_nfl = getattr(args, "nfl_baseline", None)
    if baseline_nfl is not None:
        baseline_nfl = float(baseline_nfl)

    analyzer = BiomarkerAnalyzer(disease=disease)
    res = analyzer.evaluate_biomarker_response(
        target_gene=target,
        in_vitro_rescue_pct=rescue_pct,
        baseline_plasma_nfl=baseline_nfl,
    )

    if getattr(args, "json", False):
        out = {
            "disease": res.disease,
            "target_gene": res.target_gene,
            "in_vitro_rescue_pct": res.in_vitro_rescue_pct,
            "baseline_plasma_nfl_pg_ml": res.baseline_plasma_nfl_pg_ml,
            "biomarker": res.biomarker.value,
            "predicted_week12_nfl": res.predicted_week12_nfl,
            "predicted_week24_nfl": res.predicted_week24_nfl,
            "predicted_week48_nfl": res.predicted_week48_nfl,
            "maximal_reduction_pct": res.maximal_reduction_pct,
            "clinical_surrogate_verdict": res.clinical_surrogate_verdict,
            "stratification": asdict(res.stratification),
            "longitudinal_trajectory": [asdict(pt) for pt in res.longitudinal_trajectory],
            "actionable_trial_recommendations": res.actionable_trial_recommendations,
        }
        print(json.dumps(out, indent=2))
        return 0

    header(f"Translational Fluid Biomarker Panel: {res.target_gene} in {res.disease}")
    print(f"Target Biomarker Endpoint    : {res.biomarker.value.upper()} (Simoa ultrasensitive assay)")
    print(f"Assumed In Vitro Rescue      : {res.in_vitro_rescue_pct:.1f}%")
    print(f"Baseline Patient NfL Level   : {res.baseline_plasma_nfl_pg_ml:.1f} pg/mL")
    print(f"Maximal In Vivo Suppression  : {res.maximal_reduction_pct:.1f}%")
    print(f"Predicted Week 12 NfL        : {res.predicted_week12_nfl:.1f} pg/mL")
    print(f"Predicted Week 24 NfL        : {res.predicted_week24_nfl:.1f} pg/mL")
    print(f"Predicted Week 48 NfL        : {res.predicted_week48_nfl:.1f} pg/mL")

    print("\n" + "=" * 78)
    print("Regulatory Surrogate Milestone Verdict")
    print("=" * 78)
    if res.maximal_reduction_pct >= 40.0:
        print(Fmt.green(res.clinical_surrogate_verdict))
    elif res.maximal_reduction_pct >= 25.0:
        print(Fmt.yellow(res.clinical_surrogate_verdict))
    else:
        print(Fmt.red(res.clinical_surrogate_verdict))

    print("\nLongitudinal Pharmacodynamic Trajectory")
    print("=" * 78)
    print(f"{'Week':<8s} {'Predicted NfL (pg/mL)':<26s} {'Percent Reduction':<20s}")
    print("-" * 78)
    for pt in res.longitudinal_trajectory:
        print(f"Week {pt.week:<3d} {pt.predicted_nfl_pg_ml:>10.1f} pg/mL             {pt.percent_reduction:>8.1f}%")

    print("\nClinical Trial Patient Stratification (Power & Cohort Optimization)")
    print("=" * 78)
    strat = res.stratification
    print(f"Patient Progression Tier     : {strat.tier.value.upper()} ({strat.baseline_nfl_range})")
    print(f"Expected Disease Decline Rate: {strat.expected_decline_rate_pts_month:.2f} ALSFRS-R pts/month")
    print(f"Required Sample Size (All-Comers) : {strat.unstratified_n_per_arm} patients / arm")
    print(f"Required Sample Size (Stratified) : {strat.stratified_n_per_arm} patients / arm")
    print(f"Sample Size Reduction Efficiency  : {Fmt.green(str(strat.sample_size_reduction_pct) + '% saved')}")
    print(f"Trial Feasibility Rating          : {strat.feasibility_verdict}")

    if res.actionable_trial_recommendations:
        print("\nActionable Clinical Trial Recommendations")
        print("=" * 78)
        for rec in res.actionable_trial_recommendations:
            print(f"  -> {rec}")

    print("\n" + "-" * 78)
    print("Cite: Francesco Nobile (ORCID: 0009-0001-9477-3332); run 'synaptik cite'")
    return 0


def register(sub, common) -> None:
    s = sub.add_parser(
        "biomarker", parents=[common],
        help="model fluid biomarker trajectories (plasma NfL) and optimize clinical trial stratification",
    )
    s.add_argument("--disease", choices=["als", "ftd", "ad"], default="als",
                   help="disease indication (default: als)")
    s.add_argument("--target", metavar="GENE", default="STMN2",
                   help="therapeutic target gene symbol (default: STMN2)")
    s.add_argument("--rescue-pct", type=float, default=70.0,
                   help="in vitro phenotypic/axonal rescue percentage (default: 70.0)")
    s.add_argument("--nfl-baseline", type=float, metavar="PG_ML",
                   help="baseline patient plasma NfL concentration in pg/mL (default: indication median)")
    s.add_argument("--json", action="store_true", help="output results as JSON")
    s.set_defaults(func=cmd_biomarker)
