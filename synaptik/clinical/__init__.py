"""
SYNAPTIK Clinical Protocol & Dose Escalation Engine.

First-in-Human (FIH) clinical protocol generator, cerebrospinal fluid (CSF) allometric
scaling, intrathecal 2-compartment PK/PD simulator, and adaptive Bayesian dose escalation.

Principal Investigator & Author:
Francesco Nobile (ORCID: 0009-0001-9477-3332)
Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School;
University of Pavia.
"""

from __future__ import annotations

from .allometry import (
    CSF_VOLUMES,
    AllometryResult,
    CSFCompartment,
    calculate_allometric_dose,
)
from .escalation import (
    BOINEscalationEngine,
    CohortOutcome,
    EscalationAction,
    EscalationDecision,
    EscalationDesign,
    MTPIEngine,
    ThreePlusThreeEngine,
    betai,
    calculate_beta_overdose_prob,
)
from .pkpd import (
    PKPDExposure,
    PKPDKinetics,
    PKPDPharmacodynamics,
    PKPDSimulationResult,
    PKPDTrajectoryPoint,
    simulate_intrathecal_pkpd,
)
from .protocol import (
    BiomarkerMilestoneResult,
    ClinicalProtocol,
    DLTCriteria,
    DLTEvaluationResult,
    DoseTier,
    evaluate_bedside_dlt,
    evaluate_nfl_milestone,
    generate_clinical_protocol,
)

AUTHOR = "Francesco Nobile"
ORCID = "0009-0001-9477-3332"
AFFILIATION = "Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School; University of Pavia"

__all__ = [
    # Allometry
    "CSFCompartment",
    "CSF_VOLUMES",
    "AllometryResult",
    "calculate_allometric_dose",
    # PK/PD
    "PKPDTrajectoryPoint",
    "PKPDKinetics",
    "PKPDExposure",
    "PKPDPharmacodynamics",
    "PKPDSimulationResult",
    "simulate_intrathecal_pkpd",
    # Escalation
    "EscalationAction",
    "EscalationDesign",
    "CohortOutcome",
    "EscalationDecision",
    "BOINEscalationEngine",
    "ThreePlusThreeEngine",
    "MTPIEngine",
    "betai",
    "calculate_beta_overdose_prob",
    # Protocol
    "DoseTier",
    "DLTCriteria",
    "DLTEvaluationResult",
    "evaluate_bedside_dlt",
    "BiomarkerMilestoneResult",
    "evaluate_nfl_milestone",
    "ClinicalProtocol",
    "generate_clinical_protocol",
    # Metadata
    "AUTHOR",
    "ORCID",
    "AFFILIATION",
]
