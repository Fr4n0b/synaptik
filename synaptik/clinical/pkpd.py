"""
Coupled 2-Compartment Intrathecal Pharmacokinetic / Pharmacodynamic (PK/PD) Simulator.

Models lumbar bolus injection versus continuous catheter pump infusion, longitudinal
cerebrospinal fluid clearance (CSF t1/2 ~20-30 days), CNS neural tissue retention
(tissue t1/2 ~3-6 months), Emax Hill target engagement, and Day 90 fluid biomarker (NfL) response.

Principal Investigator & Author:
Francesco Nobile (ORCID: 0009-0001-9477-3332)
Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School;
University of Pavia.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, List, Optional


@dataclass
class PKPDTrajectoryPoint:
    """Individual longitudinal time point along the PK/PD trajectory."""
    day: float
    a_csf_mg: float
    c_csf_ug_ml: float
    a_tissue_mg: float
    c_tissue_ug_g: float
    target_engagement_pct: float
    plasma_nfl_pg_ml: float
    nfl_reduction_pct: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "day": round(self.day, 2),
            "a_csf_mg": round(self.a_csf_mg, 4),
            "c_csf_ug_ml": round(self.c_csf_ug_ml, 4),
            "a_tissue_mg": round(self.a_tissue_mg, 4),
            "c_tissue_ug_g": round(self.c_tissue_ug_g, 4),
            "target_engagement_pct": round(self.target_engagement_pct, 2),
            "plasma_nfl_pg_ml": round(self.plasma_nfl_pg_ml, 2),
            "nfl_reduction_pct": round(self.nfl_reduction_pct, 2),
        }


@dataclass
class PKPDKinetics:
    """Kinetic rate constants and half-lives."""
    csf_half_life_days: float
    tissue_half_life_days: float
    k_el_csf: float
    k_el_tissue: float
    k_uptake: float
    k_efflux: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "csf_half_life_days": self.csf_half_life_days,
            "tissue_half_life_days": self.tissue_half_life_days,
            "k_el_csf": round(self.k_el_csf, 5),
            "k_el_tissue": round(self.k_el_tissue, 5),
            "k_uptake": round(self.k_uptake, 4),
            "k_efflux": round(self.k_efflux, 4),
        }


@dataclass
class PKPDExposure:
    """Pharmacokinetic exposure and peak metrics."""
    cmax_csf_ug_ml: float
    tmax_csf_days: float
    cmax_tissue_ug_g: float
    tmax_tissue_days: float
    auc_csf_mg_day_l: float
    auc_tissue_ug_day_g: float

    @property
    def auc_csf_ug_day_ml(self) -> float:
        """Alias for auc_csf_mg_day_l (1 mg*day/L == 1 ug*day/mL)."""
        return self.auc_csf_mg_day_l

    @property
    def auc_tissue(self) -> float:
        """Convenience property for tissue AUC."""
        return self.auc_tissue_ug_day_g

    def to_dict(self) -> dict[str, Any]:
        return {
            "cmax_csf_ug_ml": round(self.cmax_csf_ug_ml, 2),
            "tmax_csf_days": round(self.tmax_csf_days, 4),
            "cmax_tissue_ug_g": round(self.cmax_tissue_ug_g, 2),
            "tmax_tissue_days": round(self.tmax_tissue_days, 1),
            "auc_csf_ug_day_ml": round(self.auc_csf_mg_day_l, 2),
            "auc_tissue_ug_day_g": round(self.auc_tissue_ug_day_g, 2),
        }


@dataclass
class PKPDPharmacodynamics:
    """Pharmacodynamic efficacy and Day 90 milestone results."""
    day90_tissue_ug_g: float
    day90_target_knockdown_pct: float
    day90_predicted_nfl_suppression_pct: float
    day90_nfl_pg_ml: float
    clinical_surrogate_met: bool

    def to_dict(self) -> dict[str, Any]:
        return {
            "day90_tissue_ug_g": round(self.day90_tissue_ug_g, 2),
            "day90_target_knockdown_pct": round(self.day90_target_knockdown_pct, 2),
            "day90_predicted_nfl_suppression_pct": round(self.day90_predicted_nfl_suppression_pct, 2),
            "day90_nfl_pg_ml": round(self.day90_nfl_pg_ml, 2),
            "clinical_surrogate_met": self.clinical_surrogate_met,
        }


@dataclass
class PKPDSimulationResult:
    """Consolidated outcome of the intrathecal PK/PD simulation."""
    dose_mg: float
    route: str
    infusion_duration_days: float
    csf_half_life_days: float
    tissue_half_life_days: float
    duration_days: int
    cmax_csf_ug_ml: float
    cmax_tissue_ug_g: float
    auc_csf_mg_day_l: float
    auc_tissue_ug_day_g: float
    day90_nfl_pg_ml: float
    day90_nfl_reduction_pct: float
    day90_milestone_met: bool
    kinetics: PKPDKinetics
    exposure: PKPDExposure
    pharmacodynamics: PKPDPharmacodynamics
    trajectory: List[PKPDTrajectoryPoint] = field(default_factory=list)
    author: str = "Francesco Nobile (ORCID: 0009-0001-9477-3332)"

    def to_dict(self) -> dict[str, Any]:
        return {
            "author": self.author,
            "dose_mg": self.dose_mg,
            "route": self.route,
            "duration_days": self.duration_days,
            "kinetics": self.kinetics.to_dict(),
            "exposure": self.exposure.to_dict(),
            "pharmacodynamics": self.pharmacodynamics.to_dict(),
            "day90_milestone_met": self.day90_milestone_met,
            "timepoints": [p.to_dict() for p in self.trajectory if p.day in {0, 1, 7, 14, 28, 60, 90, 180, 270, 360}],
        }


def simulate_intrathecal_pkpd(
    dose_mg: float = 50.0,
    route: str = "bolus",
    csf_half_life_days: float = 24.0,
    tissue_half_life_days: float = 135.0,
    duration_days: int = 360,
    infusion_hours: float = 24.0,
    baseline_nfl_pg_ml: float = 65.0,
    emax: float = 100.0,
    ec50_ug_g: float = 0.15,
    hill_gamma: float = 1.5,
    nfl_tau_days: float = 56.0,
    beta_scale: float = 0.72,
) -> PKPDSimulationResult:
    """
    Simulate intrathecal pharmacokinetics and pharmacodynamics using a coupled
    2-compartment differential equation model.

    Compartments:
    1. CSF Subarachnoid Space (V_csf = 140 mL = 0.140 L)
    2. CNS Neural Tissue - Spinal Cord & Cortex (Mass = 1400 g = 1.4 kg)

    Parameters
    ----------
    dose_mg : float, default 50.0
        Total administered intrathecal dose in milligrams.
    route : {"bolus", "catheter"}, default "bolus"
        Administration route: rapid lumbar bolus vs continuous catheter infusion.
    csf_half_life_days : float, default 24.0
        Terminal CSF elimination half-life in days (reference: 20-30 days).
    tissue_half_life_days : float, default 135.0
        Terminal neural tissue metabolic retention half-life in days (reference: 90-180 days).
    duration_days : int, default 360
        Simulation observation window in days.
    infusion_hours : float, default 24.0
        Catheter pump infusion duration in hours (used if route="catheter").
    baseline_nfl_pg_ml : float, default 65.0
        Pre-treatment baseline plasma NfL concentration.
    emax : float, default 100.0
        Maximal target engagement percentage.
    ec50_ug_g : float, default 0.15
        Tissue concentration producing 50% maximal effect.
    hill_gamma : float, default 1.5
        Hill cooperativity exponent.
    nfl_tau_days : float, default 56.0
        Systemic fluid turnover time constant for neurofilament light chain (8 weeks).
    beta_scale : float, default 0.72
        In vivo dish-to-human biomarker coupling scaling factor.

    Returns
    -------
    PKPDSimulationResult
        Complete trajectory, exposure metrics, and biomarker milestone evaluation.
    """
    if dose_mg <= 0:
        raise ValueError("Administered dose must be positive")
    route_clean = route.strip().lower()
    if route_clean not in ("bolus", "catheter"):
        raise ValueError("Unsupported administration route. Must be 'bolus' or 'catheter'")
    if csf_half_life_days <= 0:
        raise ValueError("CSF half-life must be positive")
    if tissue_half_life_days <= 0:
        raise ValueError("Tissue half-life must be positive")
    if duration_days <= 0:
        raise ValueError("Simulation duration must be positive")
    if route_clean == "catheter" and infusion_hours <= 0:
        raise ValueError("Catheter infusion duration must be positive")

    # Physiological compartment constants
    v_csf_l = 0.140        # 140 mL CSF subarachnoid volume
    v_tissue_g = 1400.0    # 1.4 kg neural tissue mass

    # Kinetic constants (1/day)
    k_el_csf = math.log(2.0) / csf_half_life_days
    k_el_tissue = math.log(2.0) / tissue_half_life_days
    k_uptake = 0.12        # tissue uptake rate constant
    k_efflux = 0.015       # tissue-to-CSF efflux rate constant

    kinetics = PKPDKinetics(
        csf_half_life_days=csf_half_life_days,
        tissue_half_life_days=tissue_half_life_days,
        k_el_csf=k_el_csf,
        k_el_tissue=k_el_tissue,
        k_uptake=k_uptake,
        k_efflux=k_efflux,
    )

    t_inf_days = (infusion_hours / 24.0) if route_clean == "catheter" else 0.0035

    # Catheter dampening factor on peak local CSF exposure
    # Continuous micro-infusion prevents acute bolus concentration spikes
    catheter_flattening = 0.28 if route_clean == "catheter" else 1.0

    # Differential equations
    def derivatives(t: float, a_c: float, a_t: float) -> tuple[float, float]:
        if route_clean == "bolus":
            r_in = 0.0
        else:
            r_in = (dose_mg / t_inf_days) if t <= t_inf_days else 0.0
        d_csf = r_in - (k_uptake + k_el_csf) * a_c + k_efflux * a_t
        d_tissue = k_uptake * a_c - (k_efflux + k_el_tissue) * a_t
        return d_csf, d_tissue

    # Initial conditions
    if route_clean == "bolus":
        a_csf = float(dose_mg)
        a_tissue = 0.0
    else:
        a_csf = 0.0
        a_tissue = 0.0

    # Time-stepping integration using Runge-Kutta 4th Order (RK4)
    dt = 0.1  # days
    steps = int(round(duration_days / dt))

    t = 0.0
    cmax_csf = (a_csf / v_csf_l) * catheter_flattening if route_clean == "bolus" else 0.0
    tmax_csf = 0.0035 if route_clean == "bolus" else 0.0
    cmax_tissue = 0.0
    tmax_tissue = 0.0
    auc_csf = 0.0
    auc_tissue = 0.0

    trajectory: List[PKPDTrajectoryPoint] = []

    # Helper for PD target engagement Hill model
    def compute_pd(c_tis: float) -> float:
        if c_tis <= 0.0:
            return 0.0
        nom = emax * (c_tis ** hill_gamma)
        denom = (ec50_ug_g ** hill_gamma) + (c_tis ** hill_gamma)
        return nom / denom

    # Record initial state t=0
    c_csf_0 = (a_csf / v_csf_l) * catheter_flattening
    c_tis_0 = (a_tissue * 1000.0) / v_tissue_g
    e_0 = compute_pd(c_tis_0)
    nfl_0 = baseline_nfl_pg_ml
    red_0 = 0.0

    trajectory.append(PKPDTrajectoryPoint(
        day=0.0,
        a_csf_mg=a_csf,
        c_csf_ug_ml=c_csf_0,
        a_tissue_mg=a_tissue,
        c_tissue_ug_g=c_tis_0,
        target_engagement_pct=e_0,
        plasma_nfl_pg_ml=nfl_0,
        nfl_reduction_pct=red_0,
    ))

    day90_point: Optional[PKPDTrajectoryPoint] = None

    for step in range(1, steps + 1):
        # RK4 intermediate evaluations
        k1_c, k1_t = derivatives(t, a_csf, a_tissue)
        k2_c, k2_t = derivatives(t + 0.5 * dt, a_csf + 0.5 * dt * k1_c, a_tissue + 0.5 * dt * k1_t)
        k3_c, k3_t = derivatives(t + 0.5 * dt, a_csf + 0.5 * dt * k2_c, a_tissue + 0.5 * dt * k2_t)
        k4_c, k4_t = derivatives(t + dt, a_csf + dt * k3_c, a_tissue + dt * k3_t)

        a_csf += (dt / 6.0) * (k1_c + 2.0 * k2_c + 2.0 * k3_c + k4_c)
        a_tissue += (dt / 6.0) * (k1_t + 2.0 * k2_t + 2.0 * k3_t + k4_t)

        # Ensure mass non-negativity and mass conservation
        a_csf = max(0.0, a_csf)
        a_tissue = max(0.0, a_tissue)
        if (a_csf + a_tissue) > (dose_mg + 1e-4):
            scale_fac = dose_mg / (a_csf + a_tissue)
            a_csf *= scale_fac
            a_tissue *= scale_fac

        t = round(step * dt, 4)

        c_csf = (a_csf / v_csf_l) * catheter_flattening
        c_tis = (a_tissue * 1000.0) / v_tissue_g

        if c_csf > cmax_csf:
            cmax_csf = c_csf
            tmax_csf = t

        if c_tis > cmax_tissue:
            cmax_tissue = c_tis
            tmax_tissue = t

        auc_csf += c_csf * dt
        auc_tissue += c_tis * dt

        pd_eff = compute_pd(c_tis)
        kinetic_delay = 1.0 - math.exp(-t / nfl_tau_days)
        nfl_red = pd_eff * beta_scale * kinetic_delay
        cur_nfl = baseline_nfl_pg_ml * (1.0 - (nfl_red / 100.0))

        pt = PKPDTrajectoryPoint(
            day=t,
            a_csf_mg=a_csf,
            c_csf_ug_ml=c_csf,
            a_tissue_mg=a_tissue,
            c_tissue_ug_g=c_tis,
            target_engagement_pct=pd_eff,
            plasma_nfl_pg_ml=cur_nfl,
            nfl_reduction_pct=nfl_red,
        )

        # Store integer or major milestone days for reporting
        if abs(t - round(t)) < 1e-4:
            trajectory.append(pt)

        if abs(t - 90.0) < (dt / 2.0):
            day90_point = pt

    if day90_point is None:
        # Fallback to closest point
        day90_point = min(trajectory, key=lambda p: abs(p.day - 90.0))

    milestone_met = day90_point.nfl_reduction_pct > 30.0

    exposure = PKPDExposure(
        cmax_csf_ug_ml=cmax_csf,
        tmax_csf_days=tmax_csf,
        cmax_tissue_ug_g=cmax_tissue,
        tmax_tissue_days=tmax_tissue,
        auc_csf_mg_day_l=auc_csf,
        auc_tissue_ug_day_g=auc_tissue,
    )

    pharmacodynamics = PKPDPharmacodynamics(
        day90_tissue_ug_g=day90_point.c_tissue_ug_g,
        day90_target_knockdown_pct=day90_point.target_engagement_pct,
        day90_predicted_nfl_suppression_pct=day90_point.nfl_reduction_pct,
        day90_nfl_pg_ml=day90_point.plasma_nfl_pg_ml,
        clinical_surrogate_met=milestone_met,
    )

    return PKPDSimulationResult(
        dose_mg=dose_mg,
        route=route_clean,
        infusion_duration_days=t_inf_days,
        csf_half_life_days=csf_half_life_days,
        tissue_half_life_days=tissue_half_life_days,
        duration_days=duration_days,
        cmax_csf_ug_ml=cmax_csf,
        cmax_tissue_ug_g=cmax_tissue,
        auc_csf_mg_day_l=auc_csf,
        auc_tissue_ug_day_g=auc_tissue,
        day90_nfl_pg_ml=day90_point.plasma_nfl_pg_ml,
        day90_nfl_reduction_pct=day90_point.nfl_reduction_pct,
        day90_milestone_met=milestone_met,
        kinetics=kinetics,
        exposure=exposure,
        pharmacodynamics=pharmacodynamics,
        trajectory=trajectory,
    )
