"""
Target Accessibility & Effective ASO Hybridization Affinity Engine.

Computes net thermodynamic binding free energy accounting for target RNA
secondary structure opening cost:
  Delta G_effective = Delta G_duplex + Delta G_open
and derives the effective dissociation constant Kd at physiological temperature:
  Kd = exp(Delta G_effective / (R * T))

Principal Investigator & Author:
Francesco Nobile (ORCID: 0009-0001-9477-3332)
Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School;
University of Pavia.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

R_KCAL_MOL_K = 0.001987204


@dataclass
class EffectiveBindingResult:
    """Effective thermodynamic hybridization outcome accounting for target folding."""
    delta_g_duplex: float
    delta_g_open: float
    delta_g_effective: float
    kd_effective_molar: float
    kd_effective_nm: float
    is_favorable: bool
    temperature_c: float = 37.0
    author: str = "Francesco Nobile (ORCID: 0009-0001-9477-3332)"

    def to_dict(self) -> dict[str, Any]:
        return {
            "author": self.author,
            "delta_g_duplex_kcal_mol": round(self.delta_g_duplex, 3),
            "delta_g_open_kcal_mol": round(self.delta_g_open, 3),
            "delta_g_effective_kcal_mol": round(self.delta_g_effective, 3),
            "kd_effective_molar": self.kd_effective_molar,
            "kd_effective_nm": self.kd_effective_nm,
            "is_favorable": self.is_favorable,
            "temperature_c": self.temperature_c,
        }


def compute_effective_binding(
    duplex_delta_g: float,
    delta_g_open: float,
    temperature_c: float = 37.0,
) -> EffectiveBindingResult:
    """
    Compute effective ASO binding free energy and Kd accounting for target unwinding penalty.

    Parameters
    ----------
    duplex_delta_g : float
        Hybridization free energy of ASO-RNA duplex (kcal/mol, negative).
    delta_g_open : float
        Thermodynamic opening penalty to unfold target RNA window (kcal/mol >= 0).
    temperature_c : float, default 37.0
        Physiological temperature in Celsius.
    """
    if delta_g_open < 0.0:
        raise ValueError("delta_g_open opening penalty must be non-negative (>= 0.0)")
    temp_kelvin = temperature_c + 273.15
    if temp_kelvin <= 0.0:
        raise ValueError("Absolute temperature must be positive (> -273.15 C)")

    dg_eff = duplex_delta_g + delta_g_open
    exponent = max(min(dg_eff / (R_KCAL_MOL_K * temp_kelvin), 500.0), -500.0)
    kd_molar = math.exp(exponent)
    kd_nm = kd_molar * 1e9

    return EffectiveBindingResult(
        delta_g_duplex=duplex_delta_g,
        delta_g_open=delta_g_open,
        delta_g_effective=round(dg_eff, 4),
        kd_effective_molar=kd_molar,
        kd_effective_nm=kd_nm,
        is_favorable=(dg_eff < 0.0),
        temperature_c=temperature_c,
    )
