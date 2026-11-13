"""
PRO-ACT / Answer ALS Reference Covariance & Conditional Progression Engine.

Models untreated ALSFRS-R progression slope conditioned on onset site,
baseline plasma neurofilament light chain (NfL), age at baseline, and
empirical Bayes shrinkage toward individual pre-treatment progression slope.

Principal Investigator & Author:
Francesco Nobile (ORCID: 0009-0001-9477-3332)
Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School;
University of Pavia.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .models import PatientBaseline

# Reference population medians derived from PRO-ACT / Answer ALS ALS cohorts
REFERENCE_POP_SLOPE = -0.92          # pts/month overall reference
REFERENCE_AGE_YEARS = 55.0
REFERENCE_NFL_PG_ML = 58.0
RESIDUAL_SLOPE_STD = 0.280           # pts/month conditional residual SD

SITE_ADJUSTMENTS = {
    "limb": 0.08,         # -0.92 + 0.08 = -0.84 pts/month
    "spinal": 0.08,       # equivalent to limb-onset spinal phenotype
    "bulbar": -0.28,      # -0.92 - 0.28 = -1.20 pts/month (~0.36 faster than limb)
    "respiratory": -0.35, # rapid respiratory phenotype
}


@dataclass
class ConditionalProgression:
    """Conditional untreated progression slope distribution for a matched patient."""
    mean_slope: float
    std_slope: float
    population_conditional_slope: float
    site_adjustment: float
    nfl_adjustment: float
    age_adjustment: float


def compute_conditional_progression(baseline: PatientBaseline) -> ConditionalProgression:
    """
    Compute conditional untreated ALSFRS-R slope distribution (pts/month)
    using PRO-ACT / Answer ALS covariate adjustments and empirical Bayes shrinkage.
    """
    site_adj = SITE_ADJUSTMENTS.get(baseline.onset_site, 0.08)

    nfl = baseline.baseline_plasma_nfl if baseline.baseline_plasma_nfl is not None else REFERENCE_NFL_PG_ML
    nfl_adj = -0.35 * math.log(nfl / REFERENCE_NFL_PG_ML)

    age_adj = -0.012 * (baseline.age_years - REFERENCE_AGE_YEARS)

    pop_cond = REFERENCE_POP_SLOPE + site_adj + nfl_adj + age_adj
    pre_slope = baseline.baseline_slope if baseline.baseline_slope is not None else pop_cond

    # Empirical Bayes shrinkage: 65% individual pre-slope, 35% conditional population prior
    mean_slope = 0.65 * pre_slope + 0.35 * pop_cond

    return ConditionalProgression(
        mean_slope=mean_slope,
        std_slope=RESIDUAL_SLOPE_STD,
        population_conditional_slope=pop_cond,
        site_adjustment=site_adj,
        nfl_adjustment=nfl_adj,
        age_adjustment=age_adj,
    )
