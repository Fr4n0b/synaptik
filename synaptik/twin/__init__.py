"""
Patient Natural History Digital Twin & Bayesian Counterfactual Control Engine.

Zero-PHI ephemeral simulation of 1,000 matched untreated progression trajectories
and Bayesian posterior estimation of treatment-induced disease deceleration.

Principal Investigator & Author:
Francesco Nobile (ORCID: 0009-0001-9477-3332)
"""

from .bayesian import (
    BayesianEfficacyResult,
    EfficacyEvaluationResult,
    evaluate_treatment_efficacy,
)
from .distributions import (
    ConditionalProgression,
    compute_conditional_progression,
)
from .dossier import generate_twin_fda_dossier
from .models import (
    PatientBaseline,
    TwinTrajectory,
    VisitObservation,
)
from .twin_generator import (
    CounterfactualCohort,
    generate_counterfactual_cohort,
)

__all__ = [
    "BayesianEfficacyResult",
    "ConditionalProgression",
    "CounterfactualCohort",
    "EfficacyEvaluationResult",
    "PatientBaseline",
    "TwinTrajectory",
    "VisitObservation",
    "compute_conditional_progression",
    "evaluate_treatment_efficacy",
    "generate_counterfactual_cohort",
    "generate_twin_fda_dossier",
]
