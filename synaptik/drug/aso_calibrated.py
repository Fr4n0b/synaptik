"""Empirically calibrated ASO gapmer knockdown potency model (v2 MLP).

Fitted on 131,566 non-neurogenetics unique (gene, sequence) gapmer ASO
measurements across 327 human genes from ASO Atlas (gene-grouped validation
early stopping, val Spearman rho = 0.3775) and evaluated on a strictly
held-out neurogenetics test set of 25,130 gapmers across 12 genes (SOD1,
MAPT, HTT, SNCA, ATXN2, APP, LRRK2, SMN1, FUS, C9ORF72, BACE1, APOE):

  - Uncalibrated biophysical heuristic: Spearman rho = 0.1652, ROC-AUC = 0.5623
  - V1 linear biophysical model:        Spearman rho = 0.3268, ROC-AUC = 0.6670
  - V2 positional+thermo non-linear MLP: Spearman rho = 0.3763, ROC-AUC = 0.6929

External zero-shot validation (completely independent datasets):
  - Master_ASO_KD_clean (n=19): Spearman rho = 0.2698
  - wetlab_KD (n=8):            Spearman rho = 0.4286

Scope and caveats
-----------------
Models sequence-intrinsic, positional, and thermodynamic knockdown potency
in vitro for RNase H1 gapmers from ASO primary sequence alone. Does NOT
account for target mRNA secondary-structure accessibility, nuclear vs
cytoplasmic localization, delivery/conjugation, or transcriptome off-targets.
"""

from __future__ import annotations

import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np

_DINUCS = (
    "AA", "AC", "AG", "AT",
    "CA", "CC", "CG", "CT",
    "GA", "GC", "GG", "GT",
    "TA", "TC", "TG", "TT",
)
_MOTIFS = (
    "GGGG", "TTTT", "AAAA", "CCCC", "GGG", "CCC",
    "CCAC", "TCCC", "GACT", "TCCT", "GGGA", "ATGG",
    "GCTG", "CAGC", "TGGA", "CG",
)

CALIBRATION_METADATA = {
    "dataset": "ASO Atlas (barneyhill/aso_atlas)",
    "model": "v2_positional_thermo_mlp_96_64_32_1",
    "n_train": 131706,
    "n_train_gapmer_only": 131566,
    "n_train_genes": 327,
    "n_test_neuro_gapmer": 25130,
    "n_test_neuro_genes": 12,
    "val_spearman_rho": 0.3775,
    "heldout_neuro_spearman_rho": 0.3763,
    "heldout_neuro_auroc_70pct": 0.6929,
    "baseline_neuro_spearman_rho": 0.1652,
    "baseline_neuro_auroc_70pct": 0.5623,
}


@lru_cache(maxsize=1)
def _load_mlp_weights() -> tuple[np.ndarray, ...]:
    data = json.loads(Path(__file__).with_name("aso_mlp_v2.json").read_text(encoding="utf-8"))
    return (
        np.array(data["mu"], dtype=np.float32),
        np.array(data["sd"], dtype=np.float32),
        np.array(data["W1"], dtype=np.float32),
        np.array(data["b1"], dtype=np.float32),
        np.array(data["W2"], dtype=np.float32),
        np.array(data["b2"], dtype=np.float32),
        np.array(data["W3"], dtype=np.float32),
        np.array(data["b3"], dtype=np.float32),
    )


def _base_fracs(sub: str) -> list[float]:
    m = max(1, len(sub))
    return [sub.count("A") / m, sub.count("C") / m, sub.count("G") / m, sub.count("T") / m]


def _max_run(seq: str) -> float:
    best, cur = 1, 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur += 1
            if cur > best:
                best = cur
        else:
            cur = 1
    return float(best)


def predict_empirical_knockdown_pct(
    oligo_seq: str,
    dg37: float,
    tm: float,
    hairpin_dg: float,
    self_dimer: float,
    cpg_count: int,
) -> float:
    """Predict expected in-vitro knockdown (%) from ASO sequence & thermodynamics."""
    from .aso import ASOChemistry, compute_thermodynamics

    seq = oligo_seq.strip().upper().replace("U", "T")
    n = len(seq)
    if n < 8:
        return 0.0

    gc = (seq.count("G") + seq.count("C")) / n
    w5 = seq[:5]
    w3 = seq[-5:]
    gap = seq[5:-5] if n > 10 else seq[2:-2]
    _, _, dg_w5, _ = compute_thermodynamics(w5, chemistry=ASOChemistry.MOE_PS)
    _, _, dg_w3, _ = compute_thermodynamics(w3, chemistry=ASOChemistry.MOE_PS)
    _, _, dg_gap, _ = compute_thermodynamics(gap, chemistry=ASOChemistry.GAPMER_5_10_5)

    bin_feats: list[float] = []
    for b in range(5):
        lo = (b * n) // 5
        hi = ((b + 1) * n) // 5
        bin_feats.extend(_base_fracs(seq[lo:hi]))

    end_onehot: list[float] = []
    for pos in (0, 1, 2, -3, -2, -1):
        ch = seq[pos]
        end_onehot.extend([1.0 if ch == b else 0.0 for b in "ACGT"])

    denom = max(1, n - 1)
    di_counts = {d: 0 for d in _DINUCS}
    for i in range(n - 1):
        d = seq[i : i + 2]
        if d in di_counts:
            di_counts[d] += 1
    di_freqs = [di_counts[d] / denom for d in _DINUCS]
    motif_hits = [float(seq.count(m)) for m in _MOTIFS]

    ent = 0.0
    for f in _base_fracs(seq):
        if f > 0:
            ent -= f * math.log2(f)
    pur_frac = (seq.count("A") + seq.count("G")) / n
    max_run = _max_run(seq)

    feats = np.array(
        [
            float(n), gc, gc * gc, dg37, tm, (tm - 62.0) ** 2,
            hairpin_dg, self_dimer, float(cpg_count),
            dg_w5, dg_w3, dg_gap, abs(dg_w5 - dg_w3),
            ent, pur_frac, max_run,
        ]
        + _base_fracs(seq)
        + bin_feats
        + end_onehot
        + di_freqs
        + motif_hits,
        dtype=np.float32,
    )

    mu, sd, W1, b1, W2, b2, W3, b3 = _load_mlp_weights()
    xs = (feats - mu) / sd
    a1 = np.maximum(0.0, xs @ W1 + b1)
    a2 = np.maximum(0.0, a1 @ W2 + b2)
    pred = float((a2 @ W3 + b3).ravel()[0])
    return round(max(0.0, min(100.0, pred)), 1)
