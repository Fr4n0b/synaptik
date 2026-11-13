"""
`synaptik twin`: Patient Natural History Digital Twin & Bayesian Counterfactual Control CLI.

Zero-PHI ephemeral execution: patient data is never persisted to disk or telemetry.

Principal Investigator & Author:
Francesco Nobile (ORCID: 0009-0001-9477-3332)
Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School;
University of Pavia.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import List, Optional

from ...twin.bayesian import evaluate_treatment_efficacy
from ...twin.dossier import generate_twin_fda_dossier
from ...twin.models import PatientBaseline, VisitObservation
from ...twin.twin_generator import generate_counterfactual_cohort
from ..main import Fmt, header

AUTHOR_CITE = "Francesco Nobile (ORCID: 0009-0001-9477-3332)"


def _parse_observed_scores(raw: Optional[str]) -> List[VisitObservation]:
    if not raw:
        return []
    visits: List[VisitObservation] = []
    for item in raw.split(","):
        item = item.strip()
        if not item:
            continue
        if ":" not in item:
            raise ValueError(f"Invalid visit format '{item}'. Expected 'month:score' (e.g. 0:38.0,3:37.5).")
        m_str, s_str = item.split(":", 1)
        visits.append(VisitObservation(month=float(m_str), alsfrs_score=float(s_str)))
    return visits


def _build_baseline(args) -> PatientBaseline:
    patient_val = getattr(args, "patient", None)
    patient = patient_val if isinstance(patient_val, str) and patient_val else "EPHEMERAL-PT"
    age = float(getattr(args, "age", 55.0) or 55.0)
    site = str(getattr(args, "site", "limb") or "limb")
    duration = float(getattr(args, "duration", 12.0) or 12.0)
    alsfrs = float(getattr(args, "alsfrs", 38.0) or 38.0)
    nfl_val = getattr(args, "nfl", None)
    nfl = float(nfl_val) if isinstance(nfl_val, (int, float, str)) else None
    slope_val = getattr(args, "slope", None)
    slope = float(slope_val) if isinstance(slope_val, (int, float, str)) else None
    return PatientBaseline(
        patient_id=patient,
        age_years=age,
        onset_site=site,
        disease_duration_months=duration,
        baseline_alsfrs_r=alsfrs,
        baseline_slope=slope,
        baseline_plasma_nfl=nfl,
    )


def cmd_twin(args) -> int:
    subcmd = getattr(args, "subcommand", None) or "simulate"
    as_json = bool(getattr(args, "json", False))
    months = int(getattr(args, "months", 18) or 18)
    twins = int(getattr(args, "twins", 1000) or 1000)
    seed_raw = getattr(args, "seed", 42)
    seed = int(seed_raw) if seed_raw is not None else None

    try:
        baseline = _build_baseline(args)
        cohort = generate_counterfactual_cohort(
            baseline,
            n_twins=twins,
            duration_months=months,
            random_seed=seed,
        )
    except ValueError as err:
        print(Fmt.red(f"Digital twin simulation error: {err}"))
        return 1

    if subcmd == "simulate":
        if as_json:
            print(json.dumps(cohort.to_dict(), indent=2))
        else:
            print(header("SYNAPTIK Patient Digital Twin — Counterfactual Simulation"))
            print(f"  Patient (Ephemeral): {baseline.patient_id} ({baseline.onset_site}, age {baseline.age_years:.1f})")
            print(f"  Twins Simulated:     {cohort.n_twins} over {cohort.duration_months} months")
            print(f"  Expected Slope:      {cohort.conditional_mean_slope:.3f} +/- {cohort.conditional_std_slope:.3f} pts/mo")
            print(f"  Month 0 -> {cohort.duration_months}:       {cohort.median_trajectory[0]:.1f} -> {cohort.median_trajectory[-1]:.1f} (P10-P90: {cohort.p10_trajectory[-1]:.1f}-{cohort.p90_trajectory[-1]:.1f})")
            print(f"\n  Author: {AUTHOR_CITE}")
        return 0

    if subcmd == "evaluate":
        raw_scores = getattr(args, "observed_scores", None)
        try:
            visits = _parse_observed_scores(raw_scores)
            eff = evaluate_treatment_efficacy(cohort, visits)
        except ValueError as err:
            print(Fmt.red(f"Evaluation error: {err}"))
            return 2

        if as_json:
            out = eff.to_dict()
            out["cohort"] = cohort.to_dict()
            print(json.dumps(out, indent=2))
        else:
            print(header("SYNAPTIK Digital Twin — Bayesian Efficacy Evaluation"))
            print(f"  Observed Slope:      {eff.observed_slope:.3f} +/- {eff.observed_slope_se:.3f} pts/mo")
            print(f"  Twin Median Slope:   {eff.median_counterfactual_slope:.3f} pts/mo")
            print(f"  Delta (Slowing):     {eff.treatment_effect_delta:+.3f} pts/mo")
            print(f"  P(Deceleration > 0): {eff.p_deceleration:.4f} (BF10 = {eff.bayes_factor:.1f})")
            print(f"  Efficacy Tier:       {eff.efficacy_tier}")
            print(f"\n  Author: {AUTHOR_CITE}")
        return 0

    if subcmd == "dossier":
        raw_scores = getattr(args, "observed_scores", None)
        eff = None
        if raw_scores:
            try:
                visits = _parse_observed_scores(raw_scores)
                if len(visits) >= 2:
                    eff = evaluate_treatment_efficacy(cohort, visits)
            except ValueError as err:
                print(Fmt.red(f"Dossier evaluation error: {err}"))
                return 1
        md = generate_twin_fda_dossier(baseline, cohort, efficacy=eff)
        output_path = getattr(args, "output", None)
        if output_path:
            Path(output_path).write_text(md, encoding="utf-8")
        if as_json:
            print(json.dumps({"author": AUTHOR_CITE, "markdown": md, "output": output_path}, indent=2))
        elif not output_path:
            print(md)
        else:
            print(f"Dossier written to {output_path} (Author: {AUTHOR_CITE})")
        return 0

    print(Fmt.red(f"Unknown twin subcommand: {subcmd}"))
    return 1


def register(sub, common) -> None:
    p = sub.add_parser(
        "twin",
        parents=[common],
        help="Patient Natural History Digital Twin & Bayesian Counterfactual Control",
    )
    sp = p.add_subparsers(dest="subcommand")

    for name, help_txt in (
        ("simulate", "Simulate 1,000 matched untreated digital twin trajectories"),
        ("evaluate", "Evaluate observed visits vs counterfactual cohort (Bayesian posterior)"),
        ("dossier", "Generate computational draft natural history dossier"),
    ):
        subp = sp.add_parser(name, parents=[common], help=help_txt)
        subp.add_argument("--patient", default="EPHEMERAL-PT")
        subp.add_argument("--age", type=float, default=55.0)
        subp.add_argument("--site", default="limb", choices=["bulbar", "limb", "spinal", "respiratory"])
        subp.add_argument("--duration", type=float, default=12.0)
        subp.add_argument("--alsfrs", type=float, default=38.0)
        subp.add_argument("--slope", type=float, default=None)
        subp.add_argument("--nfl", type=float, default=None)
        subp.add_argument("--months", type=int, default=18)
        subp.add_argument("--twins", type=int, default=1000)
        subp.add_argument("--seed", type=int, default=42)
        if name in ("evaluate", "dossier"):
            subp.add_argument("--observed-scores", default=None, help="Comma-separated month:score pairs")
        if name == "dossier":
            subp.add_argument("--output", default=None, help="Output Markdown path")

    p.set_defaults(func=cmd_twin)
