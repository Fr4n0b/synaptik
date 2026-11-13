"""
Bayesian Posterior Inference Engine for N-of-1 Treatment Efficacy.

Computes the posterior probability of disease deceleration P(Delta > 0 | Patient Data),
Bayes Factor (BF10), observed linear regression slope, and clinical efficacy tier
comparing treated longitudinal ALSFRS-R visits against the matched untreated digital twin cohort.

Principal Investigator & Author:
Francesco Nobile (ORCID: 0009-0001-9477-3332)
Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School;
University of Pavia.
"""

from __future__ import annotations

import math
import statistics
from dataclasses import dataclass
from typing import Any, List

from .models import VisitObservation
from .twin_generator import CounterfactualCohort


def _normal_cdf(z: float) -> float:
    """Standard normal cumulative distribution function Phi(z)."""
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


@dataclass
class BayesianEfficacyResult:
    """Bayesian counterfactual comparison outcome for an N-of-1 treated patient."""
    p_deceleration: float
    bayes_factor: float
    observed_slope: float
    observed_slope_se: float
    median_counterfactual_slope: float
    treatment_effect_delta: float
    efficacy_tier: str
    fda_regulatory_verdict: str
    author: str = "Francesco Nobile (ORCID: 0009-0001-9477-3332)"

    def to_dict(self) -> dict[str, Any]:
        return {
            "author": self.author,
            "p_deceleration": round(self.p_deceleration, 4),
            "bayes_factor": round(self.bayes_factor, 2),
            "observed_slope": round(self.observed_slope, 4),
            "observed_slope_se": round(self.observed_slope_se, 4),
            "median_counterfactual_slope": round(self.median_counterfactual_slope, 4),
            "treatment_effect_delta": round(self.treatment_effect_delta, 4),
            "efficacy_tier": self.efficacy_tier,
            "fda_regulatory_verdict": self.fda_regulatory_verdict,
        }


EfficacyEvaluationResult = BayesianEfficacyResult


def evaluate_treatment_efficacy(
    cohort: CounterfactualCohort,
    visits: List[VisitObservation],
) -> BayesianEfficacyResult:
    """
    Evaluate N-of-1 treatment efficacy by comparing observed longitudinal
    ALSFRS-R slope against the matched untreated counterfactual digital twin cohort.

    Requires at least 2 longitudinal visit observations with distinct timepoints.
    """
    if len(visits) < 2:
        raise ValueError("At least 2 longitudinal visits are required to estimate observed slope")

    ordered = sorted(visits, key=lambda v: v.month)
    xs = [v.month for v in ordered]
    ys = [v.alsfrs_score for v in ordered]

    n = len(xs)
    x_mean = sum(xs) / n
    y_mean = sum(ys) / n
    sxx = sum((x - x_mean) ** 2 for x in xs)
    if sxx <= 0.0:
        raise ValueError("Visit timepoints must span a non-zero duration")

    sxy = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys))
    obs_slope = sxy / sxx
    intercept = y_mean - obs_slope * x_mean

    if n > 2:
        residuals = [y - (intercept + obs_slope * x) for x, y in zip(xs, ys)]
        rss = sum(r * r for r in residuals)
        se_ols = math.sqrt((rss / (n - 2)) / sxx)
    else:
        # 2 visits: assign realistic clinical measurement uncertainty (SD ~ 0.85 pts)
        se_ols = math.sqrt(2.0 * (0.85 ** 2) / sxx)

    # Minimum floor for observation SE to reflect biological/clinical measurement variability
    obs_se = max(se_ols, 0.15)

    twin_median_slope = (
        statistics.median(cohort.twin_slopes)
        if cohort.twin_slopes
        else cohort.conditional_mean_slope
    )
    twin_std = cohort.conditional_std_slope if cohort.conditional_std_slope > 0 else 0.280

    # Positive delta means observed treated decline is slower (less negative) than counterfactual
    delta = obs_slope - twin_median_slope
    combined_sd = math.sqrt(obs_se ** 2 + twin_std ** 2)

    z = delta / combined_sd
    p_decel = _normal_cdf(z)
    p_decel_clamped = min(max(p_decel, 1e-6), 1.0 - 1e-6)
    bayes_factor = p_decel_clamped / (1.0 - p_decel_clamped)

    if p_decel >= 0.975:
        tier = "Tier A: Decisive Efficacy"
        verdict = "Strong Bayesian evidence of disease deceleration vs matched natural history control."
    elif p_decel >= 0.90:
        tier = "Tier B: Probable Slowing"
        verdict = "Moderate-to-strong Bayesian evidence of functional slowing; continue longitudinal follow-up."
    elif p_decel >= 0.25:
        tier = "Tier C: Inconclusive Signal"
        verdict = "Observed trajectory overlaps substantially with untreated counterfactual distribution."
    else:
        tier = "Tier D: Potential Acceleration / No Benefit"
        verdict = "Observed decline is faster than or indistinguishable from untreated natural history."

    return BayesianEfficacyResult(
        p_deceleration=p_decel,
        bayes_factor=bayes_factor,
        observed_slope=obs_slope,
        observed_slope_se=obs_se,
        median_counterfactual_slope=twin_median_slope,
        treatment_effect_delta=delta,
        efficacy_tier=tier,
        fda_regulatory_verdict=verdict,
    )
