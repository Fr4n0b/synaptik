"""
Digital Twin Counterfactual Natural History Dossier Generator.

Formats matched untreated counterfactual simulations and Bayesian deceleration
posteriors into a structured Markdown report aligned with FDA Rare Disease Natural
History Guidance and 21 CFR Part 312 conventions (Computational Research Draft).

Principal Investigator & Author:
Francesco Nobile (ORCID: 0009-0001-9477-3332)
Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School;
University of Pavia.
"""

from __future__ import annotations

import datetime
from typing import Optional

from .bayesian import BayesianEfficacyResult
from .models import PatientBaseline
from .twin_generator import CounterfactualCohort


def generate_twin_fda_dossier(
    baseline: PatientBaseline,
    cohort: CounterfactualCohort,
    efficacy: Optional[BayesianEfficacyResult] = None,
) -> str:
    """
    Generate an audit-ready Markdown counterfactual natural history dossier
    (Computational Research Draft — Not a Regulatory Submission).
    """
    ts = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ")
    nfl_str = f"{baseline.baseline_plasma_nfl:.1f} pg/mL" if baseline.baseline_plasma_nfl is not None else "Not provided (population median 58.0 pg/mL prior used)"
    pre_slope_str = f"{baseline.baseline_slope:.3f} pts/mo" if baseline.baseline_slope is not None else "N/A"

    efficacy_section = ""
    if efficacy is not None:
        efficacy_section = f"""
## 3. Bayesian Treatment Deceleration Posterior
- **Observed Treated Slope**: {efficacy.observed_slope:.3f} ± {efficacy.observed_slope_se:.3f} ALSFRS-R pts/month
- **Counterfactual Untreated Median Slope**: {efficacy.median_counterfactual_slope:.3f} ALSFRS-R pts/month
- **Treatment Effect (Delta)**: {efficacy.treatment_effect_delta:+.3f} pts/month
- **Posterior Probability of Deceleration P(Delta > 0 | Data)**: {efficacy.p_deceleration:.3f} ({efficacy.p_deceleration * 100:.1f}%)
- **Bayes Factor (BF10)**: {efficacy.bayes_factor:.1f}
- **Efficacy Classification**: {efficacy.efficacy_tier}
- **Assessment**: {efficacy.fda_regulatory_verdict}
"""
    else:
        efficacy_section = """
## 3. Bayesian Treatment Deceleration Posterior
*No post-treatment longitudinal visits provided. Run `evaluate_treatment_efficacy()` with >= 2 follow-up visits to compute posterior probability of disease deceleration.*
"""

    t_mid = cohort.duration_months // 2
    t_end = cohort.duration_months

    md = f"""# N-of-1 Digital Twin Natural History & Counterfactual Control Report
**Patient Identifier (Ephemeral)**: {baseline.patient_id}  
**Generated On**: {ts}  
**Regulatory Format Alignment**: Structured in alignment with FDA Rare Disease Natural History Guidance and 21 CFR Part 312 external control conventions (Computational Research Draft)  

> [!CAUTION]
> **COMPUTATIONAL DRAFT — NOT A REGULATORY SUBMISSION**
> This document is a computational simulation report generated for preclinical and translational research planning. It is NOT an approved FDA regulatory filing and has not been independently audited by a biostatistician or regulatory authority.

---

## 1. Executive Summary
- **Baseline Age**: {baseline.age_years:.1f} years
- **Onset Site**: {baseline.onset_site}
- **Disease Duration**: {baseline.disease_duration_months:.1f} months
- **Baseline ALSFRS-R Score**: {baseline.baseline_alsfrs_r:.1f} / 48.0
- **Pre-Treatment Progression Slope**: {pre_slope_str}
- **Baseline Plasma NfL**: {nfl_str}

## 2. Counterfactual Cohort Simulation (PRO-ACT / Answer ALS Reference)
- **Simulated Untreated Digital Twins**: {cohort.n_twins} matched trajectories
- **Follow-up Horizon**: {cohort.duration_months} months
- **Conditional Untreated Slope**: {cohort.conditional_mean_slope:.3f} ± {cohort.conditional_std_slope:.3f} ALSFRS-R pts/month
- **Projected Untreated ALSFRS-R Trajectory**:
  - Month 0: {cohort.median_trajectory[0]:.1f} (P10–P90: {cohort.p10_trajectory[0]:.1f} – {cohort.p90_trajectory[0]:.1f})
  - Month {t_mid}: {cohort.median_trajectory[t_mid]:.1f} (P10–P90: {cohort.p10_trajectory[t_mid]:.1f} – {cohort.p90_trajectory[t_mid]:.1f})
  - Month {t_end}: {cohort.median_trajectory[t_end]:.1f} (P10–P90: {cohort.p10_trajectory[t_end]:.1f} – {cohort.p90_trajectory[t_end]:.1f})
{efficacy_section}
---

## 4. Attribution & Methodology
**Principal Investigator & Author**: Francesco Nobile (ORCID: [0009-0001-9477-3332](https://orcid.org/0009-0001-9477-3332))  
**Affiliation**: Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School; University of Pavia.  
**Privacy Guarantee**: Zero-PHI ephemeral in-memory computation. No patient identifiers or clinical measurements are persisted to disk.
"""
    return md
