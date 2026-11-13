"""
3D RNA Secondary Structure Biophysics & Steric Accessibility Module.

Provides pure-Python Zuker MFE folding, Turner 2004 thermodynamic parameters,
target opening penalty (Delta G_open), effective binding affinity, and
competitive spliceosomal displacement modeling.

Principal Investigator & Author:
Francesco Nobile (ORCID: 0009-0001-9477-3332)
"""

from .accessibility import EffectiveBindingResult, compute_effective_binding
from .competition import (
    COMPETITOR_KD_NM,
    CompetitionResult,
    DisplacementResult,
    SpliceosomeTarget,
    compute_steric_displacement,
    model_steric_displacement,
)
from .folding import (
    AccessibilityResult,
    AccessibilityTier,
    FoldingResult,
    compute_delta_g_open,
    compute_effective_affinity,
    compute_structure_features,
    compute_target_accessibility,
    fold_mfe_pure_python,
    fold_rna_mfe,
)
from .parameters import (
    get_hairpin_energy,
    get_stacking_energy,
    get_terminal_penalty,
    get_tetraloop_bonus,
)

__all__ = [
    "AccessibilityResult",
    "AccessibilityTier",
    "COMPETITOR_KD_NM",
    "CompetitionResult",
    "DisplacementResult",
    "EffectiveBindingResult",
    "FoldingResult",
    "SpliceosomeTarget",
    "compute_delta_g_open",
    "compute_effective_affinity",
    "compute_effective_binding",
    "compute_steric_displacement",
    "compute_structure_features",
    "compute_target_accessibility",
    "fold_mfe_pure_python",
    "fold_rna_mfe",
    "get_hairpin_energy",
    "get_stacking_energy",
    "get_terminal_penalty",
    "get_tetraloop_bonus",
    "model_steric_displacement",
]
