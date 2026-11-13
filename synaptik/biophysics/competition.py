"""
Spliceosomal Competitive Occlusion & Steric Displacement Engine.

Models the biophysical competition between an antisense oligonucleotide (ASO)
steric blocker and endogenous spliceosome machinery or RNA-binding proteins:
  1. U1 small nuclear RNP (U1 snRNP) at the 5' splice donor site.
  2. U2 small nuclear RNA auxiliary factor (U2AF65 / U2AF35) at polypyrimidine / 3' acceptor.
  3. Transactivation response DNA-binding protein 43 (TDP-43) at (UG)n consensus repeats.

Calculates equilibrium dissociation constant (Kd) ratios and predicts whether
physiological CSF concentration ([ASO] ~ 100-500 nM) achieves target displacement.

Principal Investigator & Author:
Francesco Nobile (ORCID: 0009-0001-9477-3332)
Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School;
University of Pavia.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum
from typing import Any, Optional

from .folding import compute_effective_affinity, compute_target_accessibility


class SpliceosomeTarget(str, Enum):
    U1_SN_RNP_DONOR = "u1_sn_rnp_donor"         # 5' splice donor (consensus: MAG|GURAGU)
    U2AF65_ACCEPTOR = "u2af65_acceptor"         # 3' polypyrimidine tract / acceptor
    TDP43_HOMOOLIGOMER = "tdp43_homooligomer"   # (UG)n cryptic repressor motif
    HNRNP_A1_SILENCER = "hnrnp_a1_silencer"     # Exonic / intronic splicing silencer (UAGGG)


# Literature equilibrium dissociation constants (Kd in nanomolar, nM)
# for endogenous human spliceosomal complexes binding pre-mRNA
ENDOGENOUS_KD_NM: dict[SpliceosomeTarget, float] = {
    SpliceosomeTarget.U1_SN_RNP_DONOR: 2.5,        # 1.0 - 5.0 nM (Roca & Krainer 2009)
    SpliceosomeTarget.U2AF65_ACCEPTOR: 12.0,       # 8.0 - 20.0 nM (Wu et al. 1999)
    SpliceosomeTarget.TDP43_HOMOOLIGOMER: 4.0,     # 2.0 - 8.0 nM (Polymenidou et al. 2011)
    SpliceosomeTarget.HNRNP_A1_SILENCER: 25.0,     # 15.0 - 40.0 nM (Burd & Dreyfuss 1994)
}

# Physiological nuclear concentrations in human motor neurons (estimated nM)
ENDOGENOUS_CONC_NM: dict[SpliceosomeTarget, float] = {
    SpliceosomeTarget.U1_SN_RNP_DONOR: 50.0,
    SpliceosomeTarget.U2AF65_ACCEPTOR: 40.0,
    SpliceosomeTarget.TDP43_HOMOOLIGOMER: 30.0,
    SpliceosomeTarget.HNRNP_A1_SILENCER: 60.0,
}


@dataclass
class CompetitionResult:
    """Quantitative outcome of ASO vs spliceosome competitive binding."""
    target_protein: SpliceosomeTarget
    target_protein_label: str
    endogenous_kd_nm: float
    endogenous_conc_nm: float
    aso_effective_delta_g: float
    aso_kd_nm: float
    aso_csf_conc_nm: float
    affinity_ratio: float               # Endogenous Kd / ASO Kd (>1 means ASO binds tighter)
    fractional_occupancy_aso: float     # Percentage [0-100%] of target sites occupied by ASO
    is_steric_displacement_favorable: bool
    verdict: str
    author: str = "Francesco Nobile (ORCID: 0009-0001-9477-3332)"

    def to_dict(self) -> dict[str, Any]:
        return {
            "author": self.author,
            "target_protein": self.target_protein.value,
            "target_protein_label": self.target_protein_label,
            "endogenous_kd_nm": self.endogenous_kd_nm,
            "endogenous_conc_nm": self.endogenous_conc_nm,
            "aso_effective_delta_g": self.aso_effective_delta_g,
            "aso_kd_nm": round(self.aso_kd_nm, 4),
            "aso_csf_conc_nm": self.aso_csf_conc_nm,
            "affinity_ratio": round(self.affinity_ratio, 2),
            "fractional_occupancy_aso_pct": round(self.fractional_occupancy_aso * 100.0, 1),
            "is_steric_displacement_favorable": self.is_steric_displacement_favorable,
            "verdict": self.verdict,
        }


def model_steric_displacement(
    delta_g_duplex: float,
    accessibility_penalty: float = 0.0,
    target_protein: SpliceosomeTarget = SpliceosomeTarget.U1_SN_RNP_DONOR,
    aso_csf_conc_nm: float = 150.0,
) -> CompetitionResult:
    """
    Model the thermodynamic equilibrium competition between an ASO and endogenous spliceosome.

    Computes fractional occupancy using competitive binding equilibrium equations:
      theta_ASO = ([ASO]/Kd_ASO) / (1 + [ASO]/Kd_ASO + [Protein]/Kd_Protein)

    Parameters
    ----------
    delta_g_duplex : float
        In silico duplex hybridization free energy in kcal/mol (negative value).
    accessibility_penalty : float, default 0.0
        Intramolecular RNA opening penalty in kcal/mol (positive value).
    target_protein : SpliceosomeTarget
        The endogenous complex being sterically blocked (e.g. U1 snRNP donor).
    aso_csf_conc_nm : float, default 150.0
        Clinical steady-state intrathecal concentration of ASO in nM (typical: 100 - 300 nM).

    Returns
    -------
    CompetitionResult
        Steric displacement metrics, occupancy percentage, and feasibility verdict.
    """
    if aso_csf_conc_nm <= 0:
        raise ValueError("ASO concentration must be strictly positive")

    dg_eff, aso_kd_nm = compute_effective_affinity(delta_g_duplex, accessibility_penalty)
    endo_kd = ENDOGENOUS_KD_NM.get(target_protein, 5.0)
    endo_conc = ENDOGENOUS_CONC_NM.get(target_protein, 40.0)

    # Affinity ratio: how many times tighter the ASO binds compared to the endogenous factor
    # affinity_ratio > 1.0 means ASO has sub-nanomolar or superior affinity
    affinity_ratio = endo_kd / max(aso_kd_nm, 1e-6)

    # Competitive binding fractional occupancy (theta)
    # theta_ASO = ( [ASO]/Kd_ASO ) / ( 1 + [ASO]/Kd_ASO + [Endo]/Kd_Endo )
    term_aso = aso_csf_conc_nm / max(aso_kd_nm, 1e-6)
    term_endo = endo_conc / max(endo_kd, 1e-6)
    theta_aso = term_aso / (1.0 + term_aso + term_endo)

    # Displacement is considered favorable if ASO occupies >= 75% of target pre-mRNAs
    is_favorable = theta_aso >= 0.75

    labels = {
        SpliceosomeTarget.U1_SN_RNP_DONOR: "U1 snRNP (5' Splice Donor)",
        SpliceosomeTarget.U2AF65_ACCEPTOR: "U2AF65 / U2AF35 (3' Splice Acceptor)",
        SpliceosomeTarget.TDP43_HOMOOLIGOMER: "TDP-43 Homooligomer (Cryptic Exon Repressor)",
        SpliceosomeTarget.HNRNP_A1_SILENCER: "hnRNP A1 (Exonic Splicing Silencer)",
    }
    label = labels.get(target_protein, target_protein.value)

    if is_favorable:
        verdict = (
            f"POTENT DISPLACEMENT: ASO achieves {theta_aso * 100:.1f}% target occupancy, "
            f"outcompeting endogenous {label} ({affinity_ratio:.1f}x higher affinity)."
        )
    else:
        verdict = (
            f"INSUFFICIENT OCCUPANCY: ASO achieves only {theta_aso * 100:.1f}% target occupancy. "
            f"Endogenous {label} retains significant binding; consider increasing sequence affinity or dose."
        )

    return CompetitionResult(
        target_protein=target_protein,
        target_protein_label=label,
        endogenous_kd_nm=endo_kd,
        endogenous_conc_nm=endo_conc,
        aso_effective_delta_g=dg_eff,
        aso_kd_nm=aso_kd_nm,
        aso_csf_conc_nm=aso_csf_conc_nm,
        affinity_ratio=affinity_ratio,
        fractional_occupancy_aso=theta_aso,
        is_steric_displacement_favorable=is_favorable,
        verdict=verdict,
    )


# Named spliceosomal / RBP competitor equilibrium Kd values (nM)
COMPETITOR_KD_NM: dict[str, float] = {
    "u1_snrnp": 25.0,
    "u2af65": 20.0,
    "u2af35": 40.0,
    "tardbp": 3.0,
    "tdp43": 3.0,
    "hnrnp_a1": 25.0,
}


@dataclass
class DisplacementResult:
    """Quantitative competitive displacement result for a named splice factor."""
    competitor_name: str
    competitor_kd_nm: float
    oligo_kd_nm: float
    displacement_ratio: float
    fractional_occupancy: float
    dominant_binder: str
    target_gene: Optional[str] = None
    author: str = "Francesco Nobile (ORCID: 0009-0001-9477-3332)"

    def to_dict(self) -> dict[str, Any]:
        return {
            "author": self.author,
            "competitor_name": self.competitor_name,
            "competitor_kd_nm": self.competitor_kd_nm,
            "oligo_kd_nm": self.oligo_kd_nm,
            "displacement_ratio": round(self.displacement_ratio, 4),
            "fractional_occupancy": round(self.fractional_occupancy, 4),
            "dominant_binder": self.dominant_binder,
            "target_gene": self.target_gene,
        }


def compute_steric_displacement(
    oligo_kd_nm: float,
    competitor_name: str,
    target_gene: Optional[str] = None,
    oligo_conc_nm: float = 150.0,
    competitor_conc_nm: float = 50.0,
) -> DisplacementResult:
    """
    Compute equilibrium competitive steric displacement of an endogenous RBP/splice factor.
    """
    if oligo_kd_nm <= 0.0:
        raise ValueError("oligo_kd_nm must be strictly positive")
    if oligo_conc_nm <= 0.0 or competitor_conc_nm < 0.0:
        raise ValueError("Concentrations must be positive")

    key = competitor_name.strip().lower().replace("-", "")
    if key not in COMPETITOR_KD_NM:
        raise ValueError(f"Unknown competitor '{competitor_name}'. Supported: {sorted(COMPETITOR_KD_NM)}")

    comp_kd = COMPETITOR_KD_NM[key]
    displacement_ratio = comp_kd / oligo_kd_nm

    term_oligo = oligo_conc_nm / oligo_kd_nm
    term_comp = competitor_conc_nm / comp_kd
    fractional_occupancy = term_oligo / (1.0 + term_oligo + term_comp)
    dominant = "oligo" if displacement_ratio > 1.0 else "competitor"

    return DisplacementResult(
        competitor_name=key,
        competitor_kd_nm=comp_kd,
        oligo_kd_nm=oligo_kd_nm,
        displacement_ratio=displacement_ratio,
        fractional_occupancy=fractional_occupancy,
        dominant_binder=dominant,
        target_gene=target_gene,
    )
