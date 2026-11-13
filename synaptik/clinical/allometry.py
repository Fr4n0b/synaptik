"""
CSF Volume Allometric Scaling for Intrathecal Therapeutics.

==============================================================================
CRITICAL CLINICAL & SAFETY NOTICE:
RESEARCH PLANNING ONLY — NOT A CLINICAL DOSE RECOMMENDATION.

Intrathecal drug administration in humans, particularly in pediatric and
infant populations, carries severe potential risks including chemical
arachnoiditis, pleocytosis, acute radiculopathy, and irreversible neurological
injury. The calculations in this module represent mathematical research
translations derived from physiological CSF volume literature (e.g., Cynomolgus
macaque ~15 mL vs human adult ~140 mL; adult vs pediatric CSF volume ratios).

This software does NOT provide medical advice, diagnosis, or treatment. It is
NOT an authorized prescription, clinical trial protocol clearance, or dose
recommendation. All clinical starting doses, escalation schemes, and safety
margins MUST be empirically determined via GLP toxicology studies and formally
reviewed and authorized under an Investigational New Drug (IND) or Clinical
Trial Authorisation (CTA) by an institutional IRB/ethics committee and national
regulatory health authorities (e.g., US FDA CDER/CBER, EMA).
==============================================================================

Principal Investigator & Author:
Francesco Nobile (ORCID: 0009-0001-9477-3332)
Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School;
University of Pavia.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any, Optional


class CSFCompartment(str, Enum):
    """Reference anatomical CSF compartments and populations."""
    NHP = "nhp"
    ADULT = "adult"
    PEDIATRIC_INFANT = "pediatric_infant"
    PEDIATRIC_CHILD = "pediatric_child"
    ADOLESCENT = "adolescent"


# Reference CSF subarachnoid volumes in milliliters (mL)
CSF_VOLUMES: dict[CSFCompartment, float] = {
    CSFCompartment.NHP: 15.0,
    CSFCompartment.ADULT: 140.0,
    CSFCompartment.PEDIATRIC_INFANT: 65.0,
    CSFCompartment.PEDIATRIC_CHILD: 100.0,
    CSFCompartment.ADOLESCENT: 130.0,
}


@dataclass
class AllometryResult:
    """Result of allometric dose translation from NHP NOAEL to human clinical dose.

    NOTE: Research planning calculation only. Not an authorized clinical dose.
    """
    noael_nhp_mg: float
    nhp_volume_ml: float
    human_volume_ml: float
    cohort_name: str
    scaling_ratio: float
    hed_mg: float
    safety_factor: float
    mrsd_mg: float
    mrsd_mg_kg: Optional[float] = None
    human_weight_kg: Optional[float] = None
    author: str = "Francesco Nobile (ORCID: 0009-0001-9477-3332)"
    disclaimer: str = "RESEARCH PLANNING ONLY — NOT A CLINICAL DOSE RECOMMENDATION"
    is_clinical_recommendation: bool = False

    def to_dict(self) -> dict[str, Any]:
        """Convert result to serializable dictionary."""
        d: dict[str, Any] = {
            "disclaimer": self.disclaimer,
            "is_clinical_recommendation": self.is_clinical_recommendation,
            "noael_nhp_mg": self.noael_nhp_mg,
            "nhp_volume_ml": self.nhp_volume_ml,
            "human_volume_ml": self.human_volume_ml,
            "cohort_name": self.cohort_name,
            "scaling_ratio": round(self.scaling_ratio, 4),
            "hed_mg": round(self.hed_mg, 4),
            "safety_factor": self.safety_factor,
            "mrsd_mg": round(self.mrsd_mg, 4),
            "author": self.author,
        }
        if self.mrsd_mg_kg is not None:
            d["mrsd_mg_kg"] = round(self.mrsd_mg_kg, 4)
        if self.human_weight_kg is not None:
            d["human_weight_kg"] = self.human_weight_kg
        return d


def calculate_allometric_dose(
    noael_mg: Optional[float] = None,
    noael_mg_kg: Optional[float] = None,
    nhp_weight_kg: float = 3.0,
    nhp_volume_ml: float = 15.0,
    human_volume_ml: float = 140.0,
    safety_factor: float = 10.0,
    human_weight_kg: Optional[float] = None,
    cohort_name: str = "adult",
    acknowledge_not_clinical: bool = False,
) -> AllometryResult:
    """
    Calculate Human Equivalent Dose (HED) and Maximum Recommended Starting Dose (MRSD)
    based on cerebrospinal fluid (CSF) volumetric allometric scaling.

    CRITICAL SAFETY REQUIREMENT:
    Preclinical research planning calculation only — NOT a clinical dose recommendation.
    Callers MUST pass acknowledge_not_clinical=True to acknowledge this.

    Parameters
    ----------
    noael_mg : float, optional
        NHP No Observed Adverse Effect Level (NOAEL) absolute dose in mg.
    noael_mg_kg : float, optional
        NHP NOAEL dose in mg/kg body weight.
    nhp_weight_kg : float, default 3.0
        Non-human primate body weight in kg.
    nhp_volume_ml : float, default 15.0
        NHP CSF subarachnoid volume in mL (reference: 12.0 - 15.0 mL).
    human_volume_ml : float, default 140.0
        Target human cohort CSF volume in mL (adult: 140-150 mL, infant: 65 mL, child: 100 mL).
    safety_factor : float, default 10.0
        FDA CDER/CBER safety factor applied to HED (standard: 10.0; range: 6.0 - 10.0).
    human_weight_kg : float, optional
        Patient body weight in kg for mg/kg normalized dose reporting.
    cohort_name : str, default "adult"
        Target cohort label ('adult', 'pediatric_infant', 'pediatric_child', etc.).
    acknowledge_not_clinical : bool, default False
        Required acknowledgment that this calculation is for preclinical research planning only.

    Returns
    -------
    AllometryResult
        Scaled dose metrics, HED, and MRSD.
    """
    if not acknowledge_not_clinical:
        raise ValueError(
            "RESEARCH PLANNING ONLY — NOT A CLINICAL DOSE RECOMMENDATION: "
            "calculate_allometric_dose() calculates mathematical allometric scaling for preclinical "
            "research planning only. It is not an authorized human dose recommendation. "
            "You must explicitly pass acknowledge_not_clinical=True to acknowledge this."
        )

    if noael_mg is None and noael_mg_kg is None:
        raise ValueError("Either noael_mg or noael_mg_kg must be provided")

    if noael_mg is not None:
        if noael_mg <= 0:
            raise ValueError("NOAEL dose must be positive")
        total_noael = float(noael_mg)
    else:
        assert noael_mg_kg is not None
        if noael_mg_kg <= 0:
            raise ValueError("NOAEL dose must be positive")
        if nhp_weight_kg <= 0:
            raise ValueError("NHP weight must be positive")
        total_noael = float(noael_mg_kg) * float(nhp_weight_kg)

    if nhp_volume_ml <= 0:
        raise ValueError("NHP CSF volume must be positive")
    if human_volume_ml <= 0:
        raise ValueError("Human CSF volume must be positive")
    if safety_factor <= 0:
        raise ValueError("Safety factor must be positive")

    scaling_ratio = human_volume_ml / nhp_volume_ml
    hed_mg = total_noael * scaling_ratio
    mrsd_mg = hed_mg / safety_factor

    mrsd_mg_kg: Optional[float] = None
    if human_weight_kg is not None:
        if human_weight_kg <= 0:
            raise ValueError("Human body weight must be positive")
        mrsd_mg_kg = mrsd_mg / human_weight_kg

    return AllometryResult(
        noael_nhp_mg=total_noael,
        nhp_volume_ml=nhp_volume_ml,
        human_volume_ml=human_volume_ml,
        cohort_name=cohort_name,
        scaling_ratio=scaling_ratio,
        hed_mg=hed_mg,
        safety_factor=safety_factor,
        mrsd_mg=mrsd_mg,
        mrsd_mg_kg=mrsd_mg_kg,
        human_weight_kg=human_weight_kg,
    )
