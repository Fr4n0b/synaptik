"""Empirically calibrated ASO gapmer knockdown potency model (v3 MLP).

Reproducibly trained via ``scripts/train_aso_mlp.py`` on non-neurogenetics
unique (gene, sequence) gapmer ASO measurements from ASO Atlas (gene-grouped
validation early stopping) with 162 features:
  - 20-position ASO base one-hot (10 from 5', 10 from 3')
  - Regional wing/gap thermodynamics, dinucleotide frequencies, sequence motifs
  - Local target mRNA accessibility (dG_open, stem-loop pairing probability,
    and +-50 nt flanking AU/GC composition) mapped via Ensembl GRCh38 cDNA

Run ``python scripts/train_aso_mlp.py`` to reproduce training and held-out
neurogenetics evaluation metrics with 95% bootstrap confidence intervals.
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
_PAIR_E = {
    ("G", "C"): -2.4, ("C", "G"): -2.4,
    ("A", "T"): -1.1, ("T", "A"): -1.1,
    ("G", "T"): -0.5, ("T", "G"): -0.5,
}
_COMP_TRANS = str.maketrans("ACGTU", "TGCAA")

CALIBRATION_METADATA = {
    "dataset": "ASO Atlas (barneyhill/aso_atlas) + Ensembl GRCh38 cDNA flanks",
    "model": "v3_accessibility_pos20_mlp_162_128_64_1",
    "n_features": 162,
    "training_script": "scripts/train_aso_mlp.py",
}


@lru_cache(maxsize=1)
def _load_mlp_weights() -> tuple[np.ndarray, ...]:
    data = json.loads(Path(__file__).with_name("aso_mlp_v3.json").read_text(encoding="utf-8"))
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


def compute_target_accessibility(
    target_site: str,
    flank_5p: str = "",
    flank_3p: str = "",
) -> list[float]:
    """Compute local target mRNA accessibility (dG_open & stem-loop pairing) from +-50nt flanks."""
    if not flank_5p and not flank_3p:
        return [0.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.5]

    f5 = flank_5p[-35:]
    f3 = flank_3p[:35]
    win = f5 + target_site + f3
    t_start = len(f5)
    t_end = t_start + len(target_site)
    L = len(win)

    f5_gc = (flank_5p.count("G") + flank_5p.count("C")) / max(1, len(flank_5p)) if flank_5p else 0.5
    f3_gc = (flank_3p.count("G") + flank_3p.count("C")) / max(1, len(flank_3p)) if flank_3p else 0.5
    f5_prox_au = (flank_5p[-15:].count("A") + flank_5p[-15:].count("T")) / max(1, len(flank_5p[-15:])) if flank_5p else 0.5
    f3_prox_au = (flank_3p[:15].count("A") + flank_3p[:15].count("T")) / max(1, len(flank_3p[:15])) if flank_3p else 0.5

    min_dg = 0.0
    paired_mask = [False] * len(target_site)
    total_stem_e = 0.0

    for stem_len in (6, 5, 4):
        for i in range(max(0, t_start - 20), min(L - stem_len, t_end)):
            s1 = win[i : i + stem_len]
            j_min = i + stem_len + 3
            j_max = min(L - stem_len, i + 36)
            for j in range(j_min, j_max + 1):
                ov1 = max(0, min(i + stem_len, t_end) - max(i, t_start))
                ov2 = max(0, min(j + stem_len, t_end) - max(j, t_start))
                if ov1 == 0 and ov2 == 0:
                    continue
                e = 3.4
                mismatches = 0
                for k in range(stem_len):
                    pe = _PAIR_E.get((s1[k], win[j + stem_len - 1 - k]))
                    if pe is None:
                        mismatches += 1
                        if mismatches > 1:
                            break
                        e += 0.8
                    else:
                        e += pe
                if mismatches <= 1 and e < -1.5:
                    if e < min_dg:
                        min_dg = e
                    total_stem_e += -e
                    if ov1 > 0:
                        for p in range(max(i, t_start), min(i + stem_len, t_end)):
                            paired_mask[p - t_start] = True
                    if ov2 > 0:
                        for p in range(max(j, t_start), min(j + stem_len, t_end)):
                            paired_mask[p - t_start] = True

    dg_open = max(0.0, -min_dg)
    paired_frac = sum(paired_mask) / max(1, len(target_site))
    stem_density = min(10.0, total_stem_e / max(1, len(target_site)))
    acc_prob = 1.0 / (1.0 + math.exp((dg_open - 3.0) / 1.5))

    pinch_dg = 0.0
    if len(f5) >= 12 and len(f3) >= 12:
        for i in range(len(f5) - 5):
            s1 = f5[i : i + 5]
            for j in range(len(f3) - 5):
                e = 2.0
                mm = 0
                for k in range(5):
                    pe = _PAIR_E.get((s1[k], f3[j + 4 - k]))
                    if pe is None:
                        mm += 1
                        if mm > 0:
                            break
                    else:
                        e += pe
                if mm == 0 and e < pinch_dg:
                    pinch_dg = e

    return [
        1.0,
        f5_gc,
        f3_gc,
        f5_prox_au,
        f3_prox_au,
        dg_open,
        paired_frac,
        stem_density,
        -pinch_dg,
        acc_prob,
    ]


def predict_empirical_knockdown_pct(
    oligo_seq: str,
    dg37: float,
    tm: float,
    hairpin_dg: float,
    self_dimer: float,
    cpg_count: int,
    flank_5p: str = "",
    flank_3p: str = "",
) -> float:
    """Predict expected in-vitro knockdown (%) from ASO sequence, thermodynamics, and target mRNA context."""
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

    pos_onehot: list[float] = []
    for pos in list(range(10)) + list(range(-10, 0)):
        ch = seq[pos] if (pos >= 0 and pos < n) or (pos < 0 and abs(pos) <= n) else "N"
        pos_onehot.extend([1.0 if ch == b else 0.0 for b in "ACGT"])

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

    target_site = seq.translate(_COMP_TRANS)[::-1]
    acc_feats = compute_target_accessibility(
        target_site,
        flank_5p.strip().upper().replace("U", "T"),
        flank_3p.strip().upper().replace("U", "T"),
    )

    feats = np.array(
        [
            float(n), gc, gc * gc, dg37, tm, (tm - 62.0) ** 2,
            hairpin_dg, self_dimer, float(cpg_count),
            dg_w5, dg_w3, dg_gap, abs(dg_w5 - dg_w3),
            ent, pur_frac, max_run,
        ]
        + _base_fracs(seq)
        + bin_feats
        + pos_onehot
        + di_freqs
        + motif_hits
        + acc_feats,
        dtype=np.float32,
    )

    mu, sd, W1, b1, W2, b2, W3, b3 = _load_mlp_weights()
    xs = (feats - mu) / sd
    a1 = np.maximum(0.0, xs @ W1 + b1)
    a2 = np.maximum(0.0, a1 @ W2 + b2)
    pred = float((a2 @ W3 + b3).ravel()[0])
    return round(max(0.0, min(100.0, pred)), 1)
