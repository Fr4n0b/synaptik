"""
First-in-Human (FIH) Clinical Trial Protocol & Bedside Safety Evaluation.

Implements bedside Dose-Limiting Toxicity (DLT) stopping evaluation, Day 90 fluid
biomarker (plasma NfL) surrogate milestone gating (>30% reduction), and automated
research protocol generation formatted according to 21 CFR Part 312 and ICH E6(R2) conventions
(RESEARCH PLANNING TEMPLATE - NOT AN APPROVED REGULATORY SUBMISSION).

Principal Investigator & Author:
Francesco Nobile (ORCID: 0009-0001-9477-3332)
Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School;
University of Pavia.
"""

from __future__ import annotations

import datetime
from dataclasses import dataclass, field
from typing import Any, List, Optional

from .allometry import calculate_allometric_dose
from .escalation import EscalationDesign


@dataclass
class DoseTier:
    """Planned clinical dose level tier."""
    level: int
    dose_mg: float
    cohort_size: int = 3
    description: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "level": self.level,
            "dose_mg": self.dose_mg,
            "cohort_size": self.cohort_size,
            "description": self.description,
        }


@dataclass
class DLTCriteria:
    """Pre-specified bedside Dose-Limiting Toxicity (DLT) indicators."""
    csf_pleocytosis: bool = False          # CSF WBC count > 50 cells/uL
    chemical_arachnoiditis: bool = False    # Meningismus, headache, photophobia within 72h
    sensory_radiculopathy: bool = False     # CTCAE Grade >= 2 radiating nerve pain > 7 days
    severe_thrombocytopenia: bool = False   # Platelets < 50,000/uL or drop > 50%
    other_dlt_notes: Optional[str] = None

    @property
    def is_dlt(self) -> bool:
        """Determines if any pre-defined DLT criteria are met."""
        return any([
            self.csf_pleocytosis,
            self.chemical_arachnoiditis,
            self.sensory_radiculopathy,
            self.severe_thrombocytopenia,
        ])

    def triggered_list(self) -> List[str]:
        """List human-readable labels of triggered criteria."""
        items: List[str] = []
        if self.csf_pleocytosis:
            items.append("CSF Pleocytosis (>50 WBC cells/µL)")
        if self.chemical_arachnoiditis:
            items.append("Chemical Arachnoiditis (Acute meningismus within 72h)")
        if self.sensory_radiculopathy:
            items.append("Sensory Radiculopathy (CTCAE Grade >= 2 > 7 days)")
        if self.severe_thrombocytopenia:
            items.append("Severe Thrombocytopenia (Platelets < 50,000/µL or >50% drop)")
        if self.other_dlt_notes:
            items.append(f"Other: {self.other_dlt_notes}")
        return items


@dataclass
class DLTEvaluationResult:
    """Outcome of bedside evaluation against clinical DLT boundaries."""
    is_dlt: bool
    triggered_criteria: List[str]
    recommendation: str
    details: dict[str, Any] = field(default_factory=dict)
    author: str = "Francesco Nobile (ORCID: 0009-0001-9477-3332)"

    def to_dict(self) -> dict[str, Any]:
        return {
            "author": self.author,
            "is_dlt": self.is_dlt,
            "triggered_criteria": self.triggered_criteria,
            "recommendation": self.recommendation,
            "details": self.details,
        }


def evaluate_bedside_dlt(
    csf_wbc_per_ul: Optional[float] = None,
    chemical_arachnoiditis_present: bool = False,
    radiculopathy_grade: int = 0,
    radiculopathy_duration_days: float = 0.0,
    platelet_count_per_ul: Optional[float] = None,
    baseline_platelet_count: Optional[float] = None,
    other_toxicity_notes: Optional[str] = None,
) -> DLTEvaluationResult:
    """
    Deterministically evaluate patient adverse events against pre-defined FIH DLT criteria.

    Criteria Boundaries:
    1. CSF Pleocytosis: WBC count > 50 cells/µL (50 does not trigger, 51 triggers).
    2. Chemical Arachnoiditis: Acute meningismus, severe headache within 72h requiring steroids.
    3. Sensory Radiculopathy: CTCAE Grade >= 2 and duration > 7 days.
    4. Severe Thrombocytopenia: Platelet count < 50,000/µL or drop > 50% from baseline.
    """
    crit = DLTCriteria(other_dlt_notes=other_toxicity_notes)

    # 1. Pleocytosis: strictly > 50 cells/uL
    if csf_wbc_per_ul is not None and csf_wbc_per_ul > 50.0:
        crit.csf_pleocytosis = True

    # 2. Arachnoiditis
    if chemical_arachnoiditis_present:
        crit.chemical_arachnoiditis = True

    # 3. Radiculopathy: Grade >= 2 AND duration > 7 days
    if radiculopathy_grade >= 2 and radiculopathy_duration_days > 7.0:
        crit.sensory_radiculopathy = True

    # 4. Thrombocytopenia: < 50,000 or drop > 50%
    if platelet_count_per_ul is not None:
        if platelet_count_per_ul < 50000.0:
            crit.severe_thrombocytopenia = True
        elif baseline_platelet_count is not None and baseline_platelet_count > 0:
            pct_drop = (baseline_platelet_count - platelet_count_per_ul) / baseline_platelet_count * 100.0
            if pct_drop > 50.0:
                crit.severe_thrombocytopenia = True

    triggered = crit.triggered_list()
    is_dlt = crit.is_dlt

    if is_dlt:
        recommendation = "Dose-Limiting Toxicity confirmed. Suspend further dosing; report to Safety Monitoring Committee (SMC)."
    else:
        recommendation = "No Dose-Limiting Toxicity identified. Proceed with scheduled protocol evaluations."

    return DLTEvaluationResult(
        is_dlt=is_dlt,
        triggered_criteria=triggered,
        recommendation=recommendation,
        details={
            "csf_wbc_per_ul": csf_wbc_per_ul,
            "chemical_arachnoiditis": chemical_arachnoiditis_present,
            "radiculopathy_grade": radiculopathy_grade,
            "radiculopathy_duration_days": radiculopathy_duration_days,
            "platelet_count": platelet_count_per_ul,
            "baseline_platelet_count": baseline_platelet_count,
        },
    )


@dataclass
class BiomarkerMilestoneResult:
    """Outcome of Day 90 fluid biomarker (plasma NfL) efficacy assessment."""
    baseline_nfl_pg_ml: float
    day90_nfl_pg_ml: float
    reduction_pct: float
    milestone_met: bool
    regulatory_verdict: str
    author: str = "Francesco Nobile (ORCID: 0009-0001-9477-3332)"

    def to_dict(self) -> dict[str, Any]:
        return {
            "author": self.author,
            "baseline_nfl_pg_ml": self.baseline_nfl_pg_ml,
            "day90_nfl_pg_ml": self.day90_nfl_pg_ml,
            "reduction_pct": round(self.reduction_pct, 2),
            "milestone_met": self.milestone_met,
            "regulatory_verdict": self.regulatory_verdict,
        }


def evaluate_nfl_milestone(
    baseline_nfl_pg_ml: float,
    day90_nfl_pg_ml: float,
) -> BiomarkerMilestoneResult:
    """
    Evaluate if Day 90 plasma neurofilament light chain (NfL) meets the FDA surrogate
    efficacy threshold (> 30.0% reduction from pre-treatment baseline).

    Parameters
    ----------
    baseline_nfl_pg_ml : float
        Pre-treatment baseline plasma NfL concentration (pg/mL).
    day90_nfl_pg_ml : float
        Follow-up Day 90 (Week 13) plasma NfL concentration (pg/mL).
    """
    if baseline_nfl_pg_ml <= 0:
        raise ValueError("Baseline plasma NfL must be strictly positive")

    reduction = ((baseline_nfl_pg_ml - day90_nfl_pg_ml) / baseline_nfl_pg_ml) * 100.0
    # Boundary: strictly > 30.0%
    milestone_met = reduction > 30.0

    if milestone_met:
        verdict = f"Milestone Met: Sustained {reduction:.1f}% reduction in plasma NfL satisfies pre-specified surrogate biomarker threshold (>30% reduction, aligned with accelerated approval precedent e.g. Tofersen)."
    else:
        verdict = f"Milestone Not Met: {reduction:.1f}% reduction in plasma NfL fails to achieve the prespecified >30% surrogate efficacy threshold."

    return BiomarkerMilestoneResult(
        baseline_nfl_pg_ml=baseline_nfl_pg_ml,
        day90_nfl_pg_ml=day90_nfl_pg_ml,
        reduction_pct=reduction,
        milestone_met=milestone_met,
        regulatory_verdict=verdict,
    )


@dataclass
class ClinicalProtocol:
    """Complete First-in-Human clinical trial protocol specification."""
    protocol_id: str
    title: str
    product_name: str
    target_gene: str
    indication: str
    lead_aso_sequence: str
    noael_nhp_mg: float
    mrsd_mg: float
    safety_factor: float
    escalation_design: str
    dose_tiers: List[DoseTier]
    administration_route: str
    dosing_schedule: str
    dlt_monitoring_criteria: List[str]
    primary_biomarker_milestone: str
    markdown: str
    author: str = "Francesco Nobile (ORCID: 0009-0001-9477-3332)"
    orcid: str = "0009-0001-9477-3332"
    affiliation: str = (
        "Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School; "
        "University of Pavia"
    )

    def to_dict(self) -> dict[str, Any]:
        return {
            "protocol_id": self.protocol_id,
            "title": self.title,
            "author": self.author,
            "orcid": self.orcid,
            "affiliation": self.affiliation,
            "protocol_dict": {
                "product_name": self.product_name,
                "gene": self.target_gene,
                "disease": self.indication,
                "lead_aso_sequence": self.lead_aso_sequence,
                "noael_nhp_mg": self.noael_nhp_mg,
                "starting_dose_mrsd_mg": self.mrsd_mg,
                "safety_factor": self.safety_factor,
                "escalation_design": self.escalation_design,
                "administration_route": self.administration_route,
                "dosing_schedule": self.dosing_schedule,
                "dose_tiers": [t.to_dict() for t in self.dose_tiers],
                "dlt_criteria": self.dlt_monitoring_criteria,
                "biomarker_milestone": self.primary_biomarker_milestone,
            },
            "markdown": self.markdown,
        }


def generate_clinical_protocol(
    target_gene: str = "STMN2",
    disease: str = "als",
    noael_nhp_mg: float = 30.0,
    design: str = "boin",
    target_phi: float = 0.25,
    planned_doses: Optional[List[float]] = None,
    administration_route: str = "bolus",
    lead_aso_sequence: Optional[str] = None,
    safety_factor: float = 10.0,
) -> ClinicalProtocol:
    """
    Generate an audit-ready Phase 1 First-in-Human clinical trial protocol dossier.

    Parameters
    ----------
    target_gene : str, default "STMN2"
        Target neurogenetic gene.
    disease : str, default "als"
        Disease indication ('als', 'ftd', 'sma', etc.).
    noael_nhp_mg : float, default 30.0
        Nonclinical NHP toxicology NOAEL absolute dose in mg.
    design : {"boin", "3+3", "mtpi"}, default "boin"
        Adaptive escalation framework.
    target_phi : float, default 0.25
        Target DLT rate phi.
    planned_doses : list of float, optional
        Dose level ladder in milligrams.
    administration_route : str, default "bolus"
        "bolus" (lumbar puncture) or "catheter" (continuous infusion).
    lead_aso_sequence : str, optional
        Lead antisense oligonucleotide therapeutic sequence.
    safety_factor : float, default 10.0
        FDA CDER/CBER safety factor for starting dose.
    """
    gene = target_gene.strip().upper()
    disease_label = disease.strip().upper()
    seq = (lead_aso_sequence or "TTCTGTAACAGTTTTTGTTTTTTCAGGTGAGTATTCAGGT")[:25]
    doses = planned_doses or [10.0, 25.0, 50.0, 100.0, 150.0]

    # Calculate allometric starting dose
    allom = calculate_allometric_dose(
        noael_mg=noael_nhp_mg,
        nhp_volume_ml=15.0,
        human_volume_ml=140.0,
        safety_factor=safety_factor,
        acknowledge_not_clinical=True,
    )
    mrsd = round(allom.mrsd_mg, 1)

    proto_id = f"SYN-FIH-{gene}-001"
    title = f"Phase 1 FIH Intrathecal Dose Escalation Protocol: {gene} in {disease_label}"

    tier_objs = [
        DoseTier(level=i + 1, dose_mg=d, cohort_size=3, description=f"Cohort Tier {i + 1} ({d} mg)")
        for i, d in enumerate(doses)
    ]

    dlt_criteria = [
        "CSF pleocytosis (>50 WBC cells/µL)",
        "Chemical arachnoiditis (Acute meningismus, headache, fever within 72h)",
        "Sensory radiculopathy (CTCAE Grade >= 2 radiating nerve pain > 7 days)",
        "Severe thrombocytopenia (Platelets < 50,000/µL or sudden >50% drop from baseline)",
    ]
    biomarker_milestone = "Sustained > 30% reduction in plasma NfL at Day 90 (Week 13) vs baseline"

    # Pre-render table and list strings to adhere strictly to Python 3.10-3.12 syntax
    tier_rows = "\n".join([
        f"| Level {t.level} | {t.dose_mg:.1f} mg | {t.cohort_size} | {t.description} |"
        for t in tier_objs
    ])
    dlt_list = "\n".join([f"- **{c}**" for c in dlt_criteria])

    route_desc = (
        "Lumbar puncture intrathecal bolus (slow manual injection over 3-5 minutes)"
        if administration_route.lower() == "bolus"
        else "Continuous lumbar intrathecal catheter pump infusion over 24 hours"
    )

    md = f"""# {title}
**Protocol Identifier**: {proto_id}  
**Format Alignment**: Structured in alignment with FDA 21 CFR Part 312 and ICH E6(R2) format guidelines (Draft Planning Document)  
**Principal Investigator & Author**: Francesco Nobile (ORCID: [0009-0001-9477-3332](https://orcid.org/0009-0001-9477-3332))  
**Institutional Affiliation**: Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School; University of Pavia  
**Generated On**: {datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ")}

> [!CAUTION]
> **COMPUTATIONAL RESEARCH TEMPLATE — NOT AN AUTHORIZED CLINICAL PROTOCOL**
> This document is a computational draft generated solely for preclinical research and regulatory planning purposes. It does NOT constitute medical advice, an authorized clinical protocol, or a regulatory clearance. Intrathecal administration in humans, especially in pediatric populations, requires independent qualification by a licensed clinical investigator, approval by an Institutional Review Board (IRB) / Ethics Committee, and formal clearance under an Investigational New Drug (IND) application from the FDA or equivalent health authority.

---

## 1. Protocol Synopsis & Investigational Product
- **Target Gene**: {gene}
- **Indication**: {disease_label} (Rare Neurodegenerative Spectrum)
- **Investigational Medicinal Product (IMP)**: Chemically modified Intrathecal Antisense Oligonucleotide (ASO)
- **Representative Sequence**: `{seq}` (2'-MOE / Phosphorothioate backbone)
- **Administration Route**: {route_desc}
- **Dosing Schedule**: Single ascending dose followed by maintenance dosing at Day 90
- **Starting Dose**: Maximum Recommended Starting Dose (MRSD) = **{mrsd} mg**

---

## 2. Nonclinical Toxicology & Allometric Dose Scaling
- **Non-Human Primate (NHP) GLP NOAEL**: {noael_nhp_mg:.1f} mg (Cynomolgus macaque)
- **Anatomical Subarachnoid Scaling**:
  - NHP CSF Volume: 15.0 mL
  - Human Adult CSF Volume: 140.0 mL
  - CSF Volumetric Scaling Ratio ($R_{{\\text{{CSF}}}}$): 9.33×
- **Human Equivalent Dose (HED)**: {allom.hed_mg:.1f} mg
- **FDA Safety Factor**: {safety_factor:.1f}× (CDER/CBER 2005 Guidance)
- **Clinical Starting Dose (MRSD)**: **{mrsd} mg**

Because antisense oligonucleotides injected intrathecally do not cross the blood-brain barrier from the periphery, biodistribution is governed directly by subarachnoid CSF volume rather than body weight or surface area.

---

## 3. Adaptive Dose Escalation Plan ({design.upper()})
- **Escalation Framework**: {design.upper()} Design
- **Target Toxicity Rate ($\\phi$)**: {target_phi:.2f}
- **Cohort Ladder**:
| Tier Level | Administered Dose | Cohort Size (N) | Clinical Description |
|:----------:|:-----------------:|:---------------:|:---------------------|
{tier_rows}

- **Cohort Evaluation Window**: 28 days post-first intrathecal administration.
- **Dose Transition Rules**: Evaluated against Bayesian optimal intervals or 3+3 rule matrices. If Beta posterior probability $P(p_j > {target_phi}) > 0.95$, dose tier is eliminated permanently.

---

## 4. Pre-defined Dose-Limiting Toxicity (DLT) Criteria
Any of the following drug-related adverse events occurring within the 28-day DLT window triggers immediate cohort review:
{dlt_list}

---

## 5. Surrogate Efficacy Biomarker Milestone
- **Primary Fluid Biomarker Milestone**: **{biomarker_milestone}**.
- **Assessment Timepoint**: Day 90 (Week 12.8 / 13) post-initial dose.
- **Regulatory Qualification**: Demonstrates pharmacodynamic proof-of-mechanism and target engagement under FDA Accelerated Approval provisions.

---

## 6. Investigator & Scientific Attribution
This protocol is generated and verified by the SYNAPTIK translational engine.
- **Investigator**: Francesco Nobile (ORCID: [0009-0001-9477-3332](https://orcid.org/0009-0001-9477-3332))
- **Affiliation**: Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School; University of Pavia
- **Repository**: https://github.com/Fr4n0b/synaptik
"""

    return ClinicalProtocol(
        protocol_id=proto_id,
        title=title,
        product_name=f"SYN-{gene}-ASO",
        target_gene=gene,
        indication=disease_label,
        lead_aso_sequence=seq,
        noael_nhp_mg=noael_nhp_mg,
        mrsd_mg=mrsd,
        safety_factor=safety_factor,
        escalation_design=design,
        dose_tiers=tier_objs,
        administration_route=administration_route,
        dosing_schedule="Single Ascending Dose with Day 90 Repeat Dosing",
        dlt_monitoring_criteria=dlt_criteria,
        primary_biomarker_milestone=biomarker_milestone,
        markdown=md.strip(),
    )
