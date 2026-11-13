"""
Data Models for Patient Natural History Digital Twin & Bayesian Counterfactual Control.

All patient data structures are strictly ephemeral and in-memory (zero PHI persistence).

Principal Investigator & Author:
Francesco Nobile (ORCID: 0009-0001-9477-3332)
Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School;
University of Pavia.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

VALID_ONSET_SITES = {"bulbar", "limb", "spinal", "respiratory"}


@dataclass
class PatientBaseline:
    """Baseline clinical parameters for an N-of-1 ALS/neurodegenerative patient."""
    patient_id: str
    age_years: float
    onset_site: str
    disease_duration_months: float
    baseline_alsfrs_r: float
    baseline_slope: Optional[float] = None
    baseline_plasma_nfl: Optional[float] = None

    def __post_init__(self) -> None:
        if self.age_years <= 0.0 or self.age_years > 120.0:
            raise ValueError(f"Invalid age_years: {self.age_years}")
        if self.disease_duration_months <= 0.0:
            raise ValueError(f"disease_duration_months must be positive, got {self.disease_duration_months}")
        if not (0.0 <= self.baseline_alsfrs_r <= 48.0):
            raise ValueError(f"baseline_alsfrs_r must be in [0, 48], got {self.baseline_alsfrs_r}")
        if self.baseline_plasma_nfl is not None and self.baseline_plasma_nfl <= 0.0:
            raise ValueError(f"baseline_plasma_nfl must be positive, got {self.baseline_plasma_nfl}")

        self.onset_site = self.onset_site.strip().lower()
        if self.onset_site not in VALID_ONSET_SITES:
            raise ValueError(
                f"Invalid onset_site '{self.onset_site}'. Must be one of {sorted(VALID_ONSET_SITES)}."
            )

        if self.baseline_slope is None:
            self.baseline_slope = (self.baseline_alsfrs_r - 48.0) / self.disease_duration_months


@dataclass
class VisitObservation:
    """Single longitudinal treated clinical visit observation."""
    month: float
    alsfrs_score: float

    def __post_init__(self) -> None:
        if self.month < 0.0:
            raise ValueError(f"Visit month cannot be negative: {self.month}")
        if not (0.0 <= self.alsfrs_score <= 48.0):
            raise ValueError(f"alsfrs_score must be in [0, 48], got {self.alsfrs_score}")


@dataclass
class TwinTrajectory:
    """Single simulated untreated counterfactual patient trajectory."""
    twin_id: int
    latent_slope: float
    monthly_scores: List[float] = field(default_factory=list)
