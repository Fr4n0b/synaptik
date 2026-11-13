"""
RNA Secondary Structure & Thermodynamic Accessibility Engine.

Computes intramolecular RNA secondary structure folding, minimum free energy (MFE),
and the opening penalty (Delta G_open) required to unwind local target structures
for antisense oligonucleotide (ASO) hybridization.

Methodological Foundation:
  - Pure-Python Zuker-style dynamic programming (V and W matrices) with no external binaries.
  - Turner 2004 nearest-neighbor parameters for Watson-Crick (GC, AU) and wobble (GU) pairs,
    hairpin loop initiation energies, terminal AU penalties, and tetraloop bonuses.
  - Thermodynamic target opening penalty: Delta G_open = Delta G_constrained - Delta G_MFE >= 0.

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

from .parameters import (
    CANONICAL_RNA_PAIRS,
    get_hairpin_energy,
    get_stacking_energy,
    get_terminal_penalty,
    get_tetraloop_bonus,
)

INF = float("inf")


class AccessibilityTier(str, Enum):
    HIGHLY_ACCESSIBLE = "highly_accessible"          # Delta G_open < 1.5 kcal/mol
    MODERATELY_ACCESSIBLE = "moderately_accessible"  # 1.5 <= Delta G_open <= 4.0 kcal/mol
    STRUCTURED_OCCLUDED = "structured_occluded"      # Delta G_open > 4.0 kcal/mol


@dataclass
class FoldingResult:
    """Predicted intramolecular RNA secondary structure and thermodynamics."""
    sequence: str
    mfe_kcal_mol: float
    dot_bracket: str
    paired_bases: int
    unpaired_bases: int
    loop_count: int
    stem_count: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "sequence": self.sequence,
            "length": len(self.sequence),
            "mfe_kcal_mol": round(self.mfe_kcal_mol, 2),
            "dot_bracket": self.dot_bracket,
            "paired_bases": self.paired_bases,
            "unpaired_bases": self.unpaired_bases,
            "paired_fraction_pct": round((self.paired_bases / max(len(self.sequence), 1)) * 100.0, 1),
            "loop_count": self.loop_count,
            "stem_count": self.stem_count,
        }


@dataclass
class AccessibilityResult:
    """Thermodynamic accessibility of an ASO target window within structured RNA."""
    target_start: int
    target_end: int
    target_sequence: str
    full_sequence_len: int
    delta_g_full_mfe: float
    delta_g_constrained_mfe: float
    delta_g_open: float
    accessibility_tier: AccessibilityTier
    accessibility_penalty_kcal_mol: float
    author: str = "Francesco Nobile (ORCID: 0009-0001-9477-3332)"

    def to_dict(self) -> dict[str, Any]:
        return {
            "author": self.author,
            "target_window": f"{self.target_start}:{self.target_end}",
            "target_length": len(self.target_sequence),
            "target_sequence": self.target_sequence,
            "full_sequence_length": self.full_sequence_len,
            "delta_g_full_mfe": round(self.delta_g_full_mfe, 2),
            "delta_g_constrained_mfe": round(self.delta_g_constrained_mfe, 2),
            "delta_g_open": round(self.delta_g_open, 2),
            "accessibility_penalty_kcal_mol": round(self.accessibility_penalty_kcal_mol, 2),
            "accessibility_tier": self.accessibility_tier.value,
        }


def _internal_loop_penalty(size: int) -> float:
    """Empirical Turner 2004 bulge / internal loop initiation penalty (kcal/mol)."""
    if size <= 0:
        return 0.0
    if size == 1:
        return 3.8
    if size == 2:
        return 2.8
    if size == 3:
        return 3.2
    if size == 4:
        return 3.6
    return 3.6 + 1.0786 * math.log(size / 4.0)


def fold_rna_mfe(
    rna_sequence: str,
    min_hairpin_loop: int = 3,
    constrained_unpaired_indices: Optional[set[int]] = None,
) -> tuple[float, str]:
    """
    Predict RNA secondary structure Minimum Free Energy (MFE) using pure-Python
    Zuker-style dynamic programming (V and W matrices) with Turner 2004 parameters.

    Parameters
    ----------
    rna_sequence : str
        Target RNA (or DNA, T->U converted) sequence.
    min_hairpin_loop : int, default 3
        Minimum unpaired bases in a hairpin loop (steric limit >= 3).
    constrained_unpaired_indices : set of int, optional
        0-based indices forced to remain unpaired (for Delta G_open calculation).

    Returns
    -------
    tuple[float, str]
        (mfe_kcal_mol, dot_bracket_notation)
    """
    seq = rna_sequence.upper().replace("T", "U").strip()
    n = len(seq)
    if n < min_hairpin_loop + 2:
        return 0.0, "." * n

    constrained = constrained_unpaired_indices or set()

    # V[i][j]: MFE of seq[i..j] given that bases i and j pair
    V = [[INF] * n for _ in range(n)]
    # v_trace[i][j] stores decomposition: ("H",), ("S",), ("I", p, q), or ("M", k)
    v_trace: list[list[Optional[tuple[Any, ...]]]] = [[None] * n for _ in range(n)]

    # W[i][j]: unrestricted MFE of seq[i..j]
    W = [[0.0] * n for _ in range(n)]
    # w_trace[i][j]: ("LEFT",), ("RIGHT",), ("PAIR",), or ("SPLIT", k)
    w_trace: list[list[tuple[Any, ...]]] = [[("LEFT",)] * n for _ in range(n)]

    max_internal_loop = 8
    multiloop_init = 3.4

    for length in range(min_hairpin_loop + 2, n + 1):
        for i in range(0, n - length + 1):
            j = i + length - 1

            # 1. Compute V[i][j] if (i, j) can form a canonical pair
            if (
                i not in constrained
                and j not in constrained
                and (seq[i], seq[j]) in CANONICAL_RNA_PAIRS
            ):
                best_v = INF
                best_v_op: Optional[tuple[Any, ...]] = None

                # Option A: Hairpin loop closing (i, j)
                loop_len = j - i - 1
                if loop_len >= min_hairpin_loop:
                    e_hp = get_hairpin_energy(loop_len) + get_terminal_penalty(seq[i], seq[j])
                    if loop_len == 4:
                        e_hp += get_tetraloop_bonus(seq[i + 1 : j])
                    if e_hp < best_v:
                        best_v = e_hp
                        best_v_op = ("H",)

                # Option B: Stacked pair (i, j) closing (i+1, j-1)
                if i + 1 < j - 1 and V[i + 1][j - 1] < INF:
                    e_stack = (
                        get_stacking_energy(seq[i] + seq[i + 1], seq[j] + seq[j - 1])
                        + V[i + 1][j - 1]
                    )
                    if e_stack < best_v:
                        best_v = e_stack
                        best_v_op = ("S",)

                # Option C: Bulge / Internal loop (i, j) closing (p, q)
                for p in range(i + 1, min(i + max_internal_loop + 2, j - min_hairpin_loop)):
                    for q in range(
                        max(p + min_hairpin_loop + 1, j - max_internal_loop - 1 + (p - i - 1)),
                        j,
                    ):
                        if p == i + 1 and q == j - 1:
                            continue
                        loop_sz = (p - i - 1) + (j - q - 1)
                        if 1 <= loop_sz <= max_internal_loop and V[p][q] < INF:
                            e_int = (
                                _internal_loop_penalty(loop_sz)
                                + get_terminal_penalty(seq[i], seq[j])
                                + get_terminal_penalty(seq[p], seq[q])
                                + V[p][q]
                            )
                            if e_int < best_v:
                                best_v = e_int
                                best_v_op = ("I", p, q)

                # Option D: Multiloop bifurcation closing (i, j)
                if j - i > 2 * (min_hairpin_loop + 2):
                    for k in range(i + 2, j - 2):
                        e_multi = multiloop_init + W[i + 1][k] + W[k + 1][j - 1]
                        if e_multi < best_v:
                            best_v = e_multi
                            best_v_op = ("M", k)

                V[i][j] = best_v
                v_trace[i][j] = best_v_op

            # 2. Compute W[i][j]
            best_w = 0.0
            best_w_op: tuple[Any, ...] = ("LEFT",)

            if W[i][j - 1] < best_w:
                best_w = W[i][j - 1]
                best_w_op = ("LEFT",)

            if W[i + 1][j] < best_w:
                best_w = W[i + 1][j]
                best_w_op = ("RIGHT",)

            if V[i][j] < INF:
                e_ext = V[i][j] + get_terminal_penalty(seq[i], seq[j])
                if e_ext < best_w:
                    best_w = e_ext
                    best_w_op = ("PAIR",)

            for k in range(i, j):
                e_split = W[i][k] + W[k + 1][j]
                if e_split < best_w:
                    best_w = e_split
                    best_w_op = ("SPLIT", k)

            W[i][j] = best_w
            w_trace[i][j] = best_w_op

    structure = ["."] * n

    def backtrack_v(i: int, j: int) -> None:
        op = v_trace[i][j]
        if op is None:
            return
        structure[i] = "("
        structure[j] = ")"
        kind = op[0]
        if kind == "H":
            return
        elif kind == "S":
            backtrack_v(i + 1, j - 1)
        elif kind == "I":
            backtrack_v(op[1], op[2])
        elif kind == "M":
            backtrack_w(i + 1, op[1])
            backtrack_w(op[1] + 1, j - 1)

    def backtrack_w(i: int, j: int) -> None:
        if i >= j or j - i < min_hairpin_loop + 1:
            return
        op = w_trace[i][j]
        kind = op[0]
        if kind == "LEFT":
            if W[i][j] < 0.0:
                backtrack_w(i, j - 1)
        elif kind == "RIGHT":
            if W[i][j] < 0.0:
                backtrack_w(i + 1, j)
        elif kind == "PAIR":
            backtrack_v(i, j)
        elif kind == "SPLIT":
            backtrack_w(i, op[1])
            backtrack_w(op[1] + 1, j)

    mfe = min(0.0, W[0][n - 1])
    if mfe < 0.0:
        backtrack_w(0, n - 1)
    return round(mfe, 2), "".join(structure)


def fold_mfe_pure_python(rna_sequence: str) -> tuple[float, str]:
    """Alias for pure-Python Zuker MFE RNA secondary structure folding."""
    return fold_rna_mfe(rna_sequence)


def compute_structure_features(sequence: str, dot_bracket: str, mfe: float) -> FoldingResult:
    """Extract structural metrics from dot-bracket representation."""
    paired = dot_bracket.count("(") + dot_bracket.count(")")
    unpaired = dot_bracket.count(".")
    stems = dot_bracket.count("(")
    loops = len([x for x in dot_bracket.split(".") if x])
    return FoldingResult(
        sequence=sequence,
        mfe_kcal_mol=mfe,
        dot_bracket=dot_bracket,
        paired_bases=paired,
        unpaired_bases=unpaired,
        loop_count=max(loops, 0),
        stem_count=stems,
    )


def compute_target_accessibility(
    rna_sequence: str,
    target_start: int,
    target_end: int,
) -> AccessibilityResult:
    """
    Compute the thermodynamic opening penalty (Delta G_open) required for an ASO
    to bind a target window within a full RNA pre-mRNA transcript.
    """
    seq = rna_sequence.upper().replace("T", "U").strip()
    n = len(seq)
    if not (0 <= target_start < target_end <= n):
        raise ValueError(f"Invalid target window [{target_start}:{target_end}] for sequence length {n}")

    mfe_full, _ = fold_rna_mfe(seq)
    constrained_indices = set(range(target_start, target_end))
    mfe_constrained, _ = fold_rna_mfe(seq, constrained_unpaired_indices=constrained_indices)

    delta_g_open = max(0.0, mfe_constrained - mfe_full)

    if delta_g_open < 1.5:
        tier = AccessibilityTier.HIGHLY_ACCESSIBLE
    elif delta_g_open <= 4.0:
        tier = AccessibilityTier.MODERATELY_ACCESSIBLE
    else:
        tier = AccessibilityTier.STRUCTURED_OCCLUDED

    target_seq = seq[target_start:target_end]
    return AccessibilityResult(
        target_start=target_start,
        target_end=target_end,
        target_sequence=target_seq,
        full_sequence_len=n,
        delta_g_full_mfe=mfe_full,
        delta_g_constrained_mfe=mfe_constrained,
        delta_g_open=round(delta_g_open, 2),
        accessibility_tier=tier,
        accessibility_penalty_kcal_mol=round(delta_g_open, 2),
    )


def compute_delta_g_open(
    rna_sequence: str,
    binding_start: int,
    binding_length: int,
) -> float:
    """
    Convenience wrapper returning the opening penalty Delta G_open (kcal/mol >= 0)
    for a binding window [binding_start : binding_start + binding_length].
    """
    seq_len = len(rna_sequence.strip())
    if binding_start < 0 or binding_length < 0 or binding_start + binding_length > seq_len:
        raise ValueError(f"Invalid binding window [{binding_start}:{binding_start + binding_length}] for sequence length {seq_len}")
    if binding_length == 0:
        return 0.0
    res = compute_target_accessibility(
        rna_sequence,
        target_start=binding_start,
        target_end=binding_start + binding_length,
    )
    return res.delta_g_open


def compute_effective_affinity(
    delta_g_duplex: float,
    accessibility_penalty: float,
    temperature_celsius: float = 37.0,
) -> tuple[float, float]:
    """
    Calculate effective thermodynamic free energy and dissociation constant (Kd).
    Delta G_effective = Delta G_duplex + accessibility_penalty.
    """
    dg_eff = delta_g_duplex + accessibility_penalty
    temp_kelvin = temperature_celsius + 273.15
    r_kcal = 0.0019872

    exponent = max(min(dg_eff / (r_kcal * temp_kelvin), 50.0), -50.0)
    kd_molar = math.exp(exponent)
    kd_nm = kd_molar * 1e9
    return round(dg_eff, 2), round(kd_nm, 4)
