"""Translational Fluid Biomarker & Patient Stratification Panel (`synaptik biomarker`).

Connects in vitro axonal rescue and molecular targets to human clinical fluid
biomarkers (Neurofilament Light Chain, NfL; plasma and CSF), modeling longitudinal
biomarker trajectories and clinical trial patient stratification.

Grounded in clinical trial and observational benchmarks:
- Miller et al. NEJM 2022 (VALOR trial: tofersen SOD1 ASO achieves ~60% plasma NfL suppression)
- Benatar et al. Ann Neurol 2018 (Pre-symptomatic NfL elevation prior to ALS onset)
- Lu et al. Lancet Neurol 2015 (NfL as a disease progression and survival predictor in ALS)
- FDA Accelerated Approval Guidance for Surrogate Endpoints (Plasma NfL in Motor Neuron Disease)
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Optional, Tuple


class BiomarkerSubtype(str, Enum):
    PLASMA_NFL = "plasma_nfl"       # Plasma Neurofilament Light (pg/mL, Simoa assay)
    CSF_NFL = "csf_nfl"             # CSF Neurofilament Light (pg/mL, ELISA/Simoa)
    PHOSPHO_TAU = "p_tau181"        # Phospho-tau 181 (pg/mL, Alzheimer's / FTD)


class PatientProgressionTier(str, Enum):
    FAST_PROGRESSOR = "fast"        # Loss > 1.2 ALSFRS-R pts/mo; baseline plasma NfL > 75 pg/mL
    INTERMEDIATE = "intermediate"   # Loss 0.6 - 1.2 pts/mo; baseline plasma NfL 40 - 75 pg/mL
    SLOW_PROGRESSOR = "slow"        # Loss < 0.6 pts/mo; baseline plasma NfL < 40 pg/mL


@dataclass
class BiomarkerTrajectoryPoint:
    week: int
    predicted_nfl_pg_ml: float
    percent_reduction: float


@dataclass
class TrialStratificationSummary:
    tier: PatientProgressionTier
    baseline_nfl_range: str
    expected_decline_rate_pts_month: float
    unstratified_n_per_arm: int
    stratified_n_per_arm: int
    sample_size_reduction_pct: float
    feasibility_verdict: str


@dataclass
class BiomarkerResponseResult:
    disease: str
    target_gene: str
    in_vitro_rescue_pct: float
    baseline_plasma_nfl_pg_ml: float
    biomarker: BiomarkerSubtype
    # Dynamics
    predicted_week12_nfl: float
    predicted_week24_nfl: float
    predicted_week48_nfl: float
    maximal_reduction_pct: float
    longitudinal_trajectory: List[BiomarkerTrajectoryPoint]
    # Clinical trial power
    stratification: TrialStratificationSummary
    clinical_surrogate_verdict: str
    actionable_trial_recommendations: List[str]


# Clinical baseline distributions for neurodegenerative indications (Plasma NfL pg/mL)
CLINICAL_NFL_BENCHMARKS: Dict[str, Dict[str, float]] = {
    "als": {
        "healthy_control_median": 12.0,
        "fast_progressor_median": 95.0,
        "intermediate_median": 58.0,
        "slow_progressor_median": 32.0,
        "clinical_cv_all_comers": 0.65,
        "clinical_cv_stratified": 0.32,
    },
    "ftd": {
        "healthy_control_median": 14.0,
        "fast_progressor_median": 85.0,
        "intermediate_median": 52.0,
        "slow_progressor_median": 28.0,
        "clinical_cv_all_comers": 0.60,
        "clinical_cv_stratified": 0.30,
    },
    "ad": {
        "healthy_control_median": 15.0,
        "fast_progressor_median": 45.0,
        "intermediate_median": 32.0,
        "slow_progressor_median": 22.0,
        "clinical_cv_all_comers": 0.50,
        "clinical_cv_stratified": 0.28,
    },
}


def compute_sample_size_two_arm(
    effect_reduction: float, # e.g. 0.30 for 30% slowing
    cv: float,               # coefficient of variation
    alpha: float = 0.05,
    power: float = 0.80,
) -> int:
    """Calculates sample size per arm using Cohen's d formula for two-sample t-test."""
    if effect_reduction <= 0 or cv <= 0:
        return 150
    # Cohen's d approx: (effect * mean) / (cv * mean) = effect / cv
    cohens_d = effect_reduction / cv
    # Normal approximation: n = 2 * (z_alpha/2 + z_beta)^2 / d^2
    z_alpha = 1.96  # for alpha = 0.05
    z_beta = 0.84   # for 80% power
    n = 2.0 * ((z_alpha + z_beta) ** 2) / (cohens_d ** 2)
    return max(12, int(math.ceil(n)))


class BiomarkerAnalyzer:
    """Translational fluid biomarker and patient cohort stratification engine."""

    def __init__(self, disease: str = "als"):
        self.disease = disease.lower().strip()
        self.benchmarks = CLINICAL_NFL_BENCHMARKS.get(self.disease, CLINICAL_NFL_BENCHMARKS["als"])

    def evaluate_biomarker_response(
        self,
        target_gene: str,
        in_vitro_rescue_pct: float = 65.0,
        baseline_plasma_nfl: Optional[float] = None,
        target_biomarker: BiomarkerSubtype = BiomarkerSubtype.PLASMA_NFL,
    ) -> BiomarkerResponseResult:
        """Projects in vitro axonal rescue into in vivo human plasma/CSF NfL trajectories."""
        gene_symbol = target_gene.upper().strip()
        rescue_pct = max(0.0, min(100.0, float(in_vitro_rescue_pct)))

        # Default baseline if not supplied: use intermediate clinical median
        if baseline_plasma_nfl is None or baseline_plasma_nfl <= 0:
            baseline_nfl = self.benchmarks["intermediate_median"]
        else:
            baseline_nfl = float(baseline_plasma_nfl)

        # Dish-to-human biomarker scaling factor (beta)
        # Accounting for biological clearance and systemic dilution
        scaling_beta = 0.72

        # Half-life tau of NfL clearance (weeks)
        tau_clearance_weeks = 8.0

        # Maximal achievable reduction percentage in vivo
        max_reduction_pct = rescue_pct * scaling_beta

        # Generate longitudinal trajectory (weeks 0, 4, 8, 12, 16, 24, 36, 48)
        weeks = [0, 4, 8, 12, 16, 24, 36, 48]
        trajectory: List[BiomarkerTrajectoryPoint] = []

        for w in weeks:
            if w == 0:
                current_nfl = baseline_nfl
                current_red = 0.0
            else:
                # Pharmacodynamic onset kinetics: 1 - exp(-t / tau)
                kinetic_fraction = 1.0 - math.exp(-w / tau_clearance_weeks)
                current_red = max_reduction_pct * kinetic_fraction
                current_nfl = baseline_nfl * (1.0 - (current_red / 100.0))

            trajectory.append(BiomarkerTrajectoryPoint(
                week=w,
                predicted_nfl_pg_ml=round(current_nfl, 1),
                percent_reduction=round(current_red, 1),
            ))

        w12 = next(p.predicted_nfl_pg_ml for p in trajectory if p.week == 12)
        w24 = next(p.predicted_nfl_pg_ml for p in trajectory if p.week == 24)
        w48 = next(p.predicted_nfl_pg_ml for p in trajectory if p.week == 48)

        # Stratification & Trial Power calculations
        # Determine patient tier based on baseline NfL
        if baseline_nfl >= 75.0:
            tier = PatientProgressionTier.FAST_PROGRESSOR
            baseline_range = "> 75 pg/mL"
            decline_rate = 1.35
        elif baseline_nfl >= 40.0:
            tier = PatientProgressionTier.INTERMEDIATE
            baseline_range = "40 - 75 pg/mL"
            decline_rate = 0.85
        else:
            tier = PatientProgressionTier.SLOW_PROGRESSOR
            baseline_range = "< 40 pg/mL"
            decline_rate = 0.40

        # Sample size calculations for 30% therapeutic slowing of progression
        cv_all = self.benchmarks["clinical_cv_all_comers"]
        cv_strat = self.benchmarks["clinical_cv_stratified"]

        n_unstrat = compute_sample_size_two_arm(effect_reduction=0.30, cv=cv_all)
        n_strat = compute_sample_size_two_arm(effect_reduction=0.30, cv=cv_strat)
        sample_reduction = round(((n_unstrat - n_strat) / n_unstrat) * 100.0, 1)

        strat_summary = TrialStratificationSummary(
            tier=tier,
            baseline_nfl_range=baseline_range,
            expected_decline_rate_pts_month=decline_rate,
            unstratified_n_per_arm=n_unstrat,
            stratified_n_per_arm=n_strat,
            sample_size_reduction_pct=sample_reduction,
            feasibility_verdict="HIGH CLINICAL TRIAL FEASIBILITY" if n_strat <= 60 else "MODERATE FEASIBILITY",
        )

        # Clinical regulatory surrogate verdict
        if max_reduction_pct >= 40.0:
            surrogate_verdict = (
                f"SUBSTANTIAL SURROGATE EFFICACY: Projected {max_reduction_pct:.1f}% NfL reduction exceeds "
                "the clinical surrogate milestone (e.g. tofersen VALOR trial ~50-60%), supporting accelerated approval eligibility."
            )
        elif max_reduction_pct >= 25.0:
            surrogate_verdict = (
                f"MODERATE SURROGATE EFFICACY: Projected {max_reduction_pct:.1f}% NfL reduction demonstrates "
                "biomarker engagement; evaluate combination therapy or increased exposure."
            )
        else:
            surrogate_verdict = (
                f"LOW SURROGATE EFFICACY: Projected {max_reduction_pct:.1f}% NfL reduction may be insufficient "
                "to overcome clinical trial noise without extreme cohort enrichment."
            )

        # Actionable recommendations
        recs: List[str] = [
            f"Enrich clinical trial with baseline {tier.value} progressors ({baseline_range}) to reduce required sample size from {n_unstrat} down to {n_strat} patients/arm (-{sample_reduction}%).",
            f"Set primary fluid pharmacodynamic endpoint at Week 24 (predicted NfL: {w24:.1f} pg/mL vs baseline {baseline_nfl:.1f} pg/mL).",
            "Monitor patient safety: acute transient NfL spikes during early dosing denote procedural or microglial stress rather than treatment failure.",
        ]

        return BiomarkerResponseResult(
            disease=self.disease.upper(),
            target_gene=gene_symbol,
            in_vitro_rescue_pct=rescue_pct,
            baseline_plasma_nfl_pg_ml=round(baseline_nfl, 1),
            biomarker=target_biomarker,
            predicted_week12_nfl=w12,
            predicted_week24_nfl=w24,
            predicted_week48_nfl=w48,
            maximal_reduction_pct=round(max_reduction_pct, 1),
            longitudinal_trajectory=trajectory,
            stratification=strat_summary,
            clinical_surrogate_verdict=surrogate_verdict,
            actionable_trial_recommendations=recs,
        )
