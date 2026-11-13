"""
Turner 2004 Nearest-Neighbor Biophysical Parameters for RNA Secondary Structure.

Defines Watson-Crick and wobble base stacking energies, terminal penalties,
hairpin loop initiation energies, and ultra-stable tetraloop thermodynamic bonuses
derived from Turner & Mathews (Nucleic Acids Res 2004).

Principal Investigator & Author:
Francesco Nobile (ORCID: 0009-0001-9477-3332)
Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School;
University of Pavia.
"""

from __future__ import annotations

import math
from typing import Dict, Tuple

# Turner 2004 base stacking free energies (kcal/mol at 37 C, 1 M NaCl)
# Key: (5'-pair1-3', 3'-pair2-5' antiparallel orientation)
TURNER_STACKING: Dict[Tuple[str, str], float] = {
    # Watson-Crick stacks
    ("GC", "CG"): -3.42,
    ("CG", "GC"): -2.36,
    ("GG", "CC"): -3.26,
    ("CC", "GG"): -3.26,
    ("GA", "CU"): -2.35,
    ("AG", "UC"): -2.08,
    ("CA", "GU"): -2.11,
    ("AC", "UG"): -2.24,
    ("AU", "UA"): -1.10,
    ("UA", "AU"): -1.33,
    ("AA", "UU"): -0.93,
    ("UU", "AA"): -0.93,
    # GU wobble stacks
    ("GG", "UU"): -1.40,
    ("UU", "GG"): -1.40,
    ("UG", "AC"): -2.11,
    ("GU", "CA"): -1.40,
    ("GA", "UU"): -1.30,
    ("AG", "UU"): -1.00,
    ("UG", "UG"): -0.50,
    ("GU", "GU"): -0.50,
}

# Standard Watson-Crick and wobble pairs
CANONICAL_RNA_PAIRS = {
    ("A", "U"), ("U", "A"),
    ("G", "C"), ("C", "G"),
    ("G", "U"), ("U", "G"),
}

# Empirical hairpin loop initiation free energies (kcal/mol at 37 C)
HAIRPIN_INITIATION: Dict[int, float] = {
    3: 5.70,
    4: 5.60,
    5: 5.50,
    6: 5.40,
    7: 5.50,
    8: 5.60,
    9: 5.70,
}

# Ultra-stable tetraloop bonuses (kcal/mol)
TETRALOOP_BONUSES: Dict[str, float] = {
    "GAAA": -3.00,
    "GCAA": -2.80,
    "UACG": -2.00,
    "GAGA": -2.50,
    "GGAA": -2.00,
    "CUUG": -2.80,
    "UNCG": -2.00,
}


def get_stacking_energy(pair1: str, pair2: str) -> float:
    """
    Return Turner 2004 stacking free energy for adjacent base pairs (kcal/mol).
    Returns >= 0.0 for non-canonical or destabilizing configurations.
    """
    p1 = pair1.upper().replace("T", "U")
    p2 = pair2.upper().replace("T", "U")
    if (p1, p2) in TURNER_STACKING:
        return TURNER_STACKING[(p1, p2)]
    # Check reverse complement symmetry
    if (p2, p1) in TURNER_STACKING:
        return TURNER_STACKING[(p2, p1)]
    return 0.0


def get_terminal_penalty(b1: str, b2: str) -> float:
    """
    Return Turner terminal AU penalty (+0.45 kcal/mol for terminal A-U pairs, 0.0 for G-C).
    """
    pair = (b1.upper().replace("T", "U"), b2.upper().replace("T", "U"))
    if pair in (("A", "U"), ("U", "A")):
        return 0.45
    return 0.0


def get_hairpin_energy(length: int) -> float:
    """
    Return hairpin loop initiation free energy (kcal/mol).
    Loops < 3 are sterically forbidden (returns infinity).
    Loops > 9 scale via Jacobson-Stockmayer polymer entropy scaling.
    """
    if length < 3:
        return float("inf")
    if length in HAIRPIN_INITIATION:
        return HAIRPIN_INITIATION[length]
    if length > 9:
        # Jacobson-Stockmayer formula: Delta G(L) = Delta G(6) + 1.0786 * ln(L / 6.0)
        return round(5.40 + 1.0786 * math.log(length / 6.0), 3)
    return 5.40


def get_tetraloop_bonus(sequence: str) -> float:
    """
    Return stabilizing free energy bonus for ultra-stable RNA tetraloops.
    """
    seq = sequence.upper().replace("T", "U")
    return TETRALOOP_BONUSES.get(seq, 0.0)
