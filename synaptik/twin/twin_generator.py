"""
Monte Carlo Digital Twin Cohort Simulator (PRO-ACT / Answer ALS Reference).

Simulates N matched untreated counterfactual progression trajectories over a
target clinical follow-up horizon with realistic ALSFRS-R measurement noise,
ceiling clamping at 48.0, and an absorbing functional floor at 0.0.

Principal Investigator & Author:
Francesco Nobile (ORCID: 0009-0001-9477-3332)
Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School;
University of Pavia.
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from typing import Any, List, Optional

from .distributions import compute_conditional_progression
from .models import PatientBaseline


@dataclass
class CounterfactualCohort:
    """Simulated cohort of N untreated counterfactual digital twin trajectories."""
    baseline: PatientBaseline
    n_twins: int
    duration_months: int
    trajectories: List[List[float]]
    twin_slopes: List[float]
    median_trajectory: List[float]
    p10_trajectory: List[float]
    p90_trajectory: List[float]
    conditional_mean_slope: float
    conditional_std_slope: float
    author: str = "Francesco Nobile (ORCID: 0009-0001-9477-3332)"

    def to_dict(self) -> dict[str, Any]:
        return {
            "author": self.author,
            "patient_id": self.baseline.patient_id,
            "n_twins": self.n_twins,
            "duration_months": self.duration_months,
            "conditional_mean_slope": round(self.conditional_mean_slope, 4),
            "conditional_std_slope": round(self.conditional_std_slope, 4),
            "median_trajectory": [round(x, 2) for x in self.median_trajectory],
            "p10_trajectory": [round(x, 2) for x in self.p10_trajectory],
            "p90_trajectory": [round(x, 2) for x in self.p90_trajectory],
        }


def _quantile(sorted_vals: List[float], q: float) -> float:
    if not sorted_vals:
        return 0.0
    idx = q * (len(sorted_vals) - 1)
    lo = int(idx)
    hi = min(lo + 1, len(sorted_vals) - 1)
    frac = idx - lo
    return sorted_vals[lo] * (1.0 - frac) + sorted_vals[hi] * frac


def generate_counterfactual_cohort(
    baseline: PatientBaseline,
    n_twins: int = 1000,
    duration_months: int = 18,
    random_seed: Optional[int] = None,
    measurement_noise_sd: float = 0.85,
) -> CounterfactualCohort:
    """
    Generate N matched untreated counterfactual progression trajectories.

    Parameters
    ----------
    baseline : PatientBaseline
        Patient clinical baseline parameters.
    n_twins : int, default 1000
        Number of digital twin trajectories to simulate.
    duration_months : int, default 18
        Follow-up horizon in months (produces duration_months + 1 timepoints: 0..T).
    random_seed : int, optional
        Seed for bitwise-reproducible Monte Carlo simulation.
    measurement_noise_sd : float, default 0.85
        Per-visit clinical ALSFRS-R assessment noise standard deviation (pts).
    """
    if n_twins <= 0:
        raise ValueError("n_twins must be positive")
    if duration_months <= 0:
        raise ValueError("duration_months must be positive")

    cond = compute_conditional_progression(baseline)
    rng = random.Random(random_seed)

    trajectories: List[List[float]] = []
    twin_slopes: List[float] = []

    for _ in range(n_twins):
        slope = rng.gauss(cond.mean_slope, cond.std_slope)
        twin_slopes.append(slope)
        traj: List[float] = [round(baseline.baseline_alsfrs_r, 3)]
        absorbed = baseline.baseline_alsfrs_r <= 0.0

        for t in range(1, duration_months + 1):
            if absorbed:
                traj.append(0.0)
                continue
            noise = rng.gauss(0.0, measurement_noise_sd) if measurement_noise_sd > 0 else 0.0
            raw = baseline.baseline_alsfrs_r + slope * t + noise
            score = max(0.0, min(48.0, raw))
            if score <= 0.0:
                score = 0.0
                absorbed = True
            traj.append(round(score, 3))
        trajectories.append(traj)

    median_traj: List[float] = []
    p10_traj: List[float] = []
    p90_traj: List[float] = []
    for t in range(duration_months + 1):
        col = sorted(traj[t] for traj in trajectories)
        p10_traj.append(round(_quantile(col, 0.10), 3))
        median_traj.append(round(_quantile(col, 0.50), 3))
        p90_traj.append(round(_quantile(col, 0.90), 3))

    return CounterfactualCohort(
        baseline=baseline,
        n_twins=n_twins,
        duration_months=duration_months,
        trajectories=trajectories,
        twin_slopes=twin_slopes,
        median_trajectory=median_traj,
        p10_trajectory=p10_traj,
        p90_trajectory=p90_traj,
        conditional_mean_slope=cond.mean_slope,
        conditional_std_slope=cond.std_slope,
    )
