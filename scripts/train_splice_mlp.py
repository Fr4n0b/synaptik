"""Reproducible training and 5-fold exon-grouped CV for SYNAPTIK splice MLP (v2).

Reads MFASS (`~/.synaptik/benchmarks/mfass_snv_clean.txt`), logs SHA256
checksum, hyperparameters, and exon-grouped 5-fold splits to disk, trains the
18-feature MLP (`18 -> 32 -> 1`) with fixed seed=42, and reports out-of-fold
AUROC and AUPRC with 95% bootstrap CIs.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.stats import mannwhitneyu

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from synaptik.sequence.splice import (
    _BRANCH_MOTIFS,
    _MFASS_ESE_HEXAMERS,
    _MFASS_ESS_HEXAMERS,
    _NOVA_MOTIFS,
    _PTBP_MOTIFS,
    _RBFOX_MOTIFS,
    _TDP43_MOTIFS,
    _count_hexamers,
    _count_tuple_motifs,
    scan_splice_sites,
)

SEED = 42
HIDDEN = 32
EPOCHS = 40
BATCH_SIZE = 512
LR = 5e-3
L2 = 1e-3

DEFAULT_MFASS = Path.home() / ".synaptik" / "benchmarks" / "mfass_snv_clean.txt"
SPLITS_OUT = Path.home() / ".synaptik" / "cache" / "mfass_cv_splits.json"
WEIGHTS_OUT = REPO / "synaptik" / "sequence" / "splice_mlp_v2.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def featurize_snv(ref_u: str, alt_u: str) -> np.ndarray:
    ref_donor = dict(scan_splice_sites(ref_u, "donor"))
    alt_donor = dict(scan_splice_sites(alt_u, "donor"))
    ref_acceptor = dict(scan_splice_sites(ref_u, "acceptor"))
    alt_acceptor = dict(scan_splice_sites(alt_u, "acceptor"))

    diffs = [i for i, (a, b) in enumerate(zip(ref_u, alt_u)) if a != b]
    vpos = diffs[0] if diffs else len(ref_u) // 2
    wlo, whi = max(0, vpos - 15), min(len(alt_u), vpos + 16)

    d_pos = ({p for p in ref_donor if wlo <= p < whi} | {p for p in alt_donor if wlo <= p < whi}) or (set(ref_donor) | set(alt_donor))
    a_pos = ({p for p in ref_acceptor if wlo <= p < whi} | {p for p in alt_acceptor if wlo <= p < whi}) or (set(ref_acceptor) | set(alt_acceptor))

    max_d_loss = max((ref_donor.get(p, 0.0) - alt_donor.get(p, 0.0) for p in d_pos), default=0.0)
    max_d_gain = max((alt_donor.get(p, 0.0) - ref_donor.get(p, 0.0) for p in d_pos), default=0.0)
    max_a_loss = max((ref_acceptor.get(p, 0.0) - alt_acceptor.get(p, 0.0) for p in a_pos), default=0.0)
    max_a_gain = max((alt_acceptor.get(p, 0.0) - ref_acceptor.get(p, 0.0) for p in a_pos), default=0.0)
    max_ref_d = max((ref_donor.get(p, 0.0) for p in d_pos), default=0.0)
    max_ref_a = max((ref_acceptor.get(p, 0.0) for p in a_pos), default=0.0)

    lo, hi = max(0, vpos - 6), min(len(ref_u), vpos + 7)
    rw, aw = ref_u[lo:hi], alt_u[lo:hi]
    ese_loss = float(_count_hexamers(rw, _MFASS_ESE_HEXAMERS) - _count_hexamers(aw, _MFASS_ESE_HEXAMERS))
    ess_gain = float(_count_hexamers(aw, _MFASS_ESS_HEXAMERS) - _count_hexamers(rw, _MFASS_ESS_HEXAMERS))
    tdp_delta = float(_count_tuple_motifs(aw, _TDP43_MOTIFS) - _count_tuple_motifs(rw, _TDP43_MOTIFS))
    ptbp_delta = float(_count_tuple_motifs(aw, _PTBP_MOTIFS) - _count_tuple_motifs(rw, _PTBP_MOTIFS))
    rbfox_delta = float(_count_tuple_motifs(aw, _RBFOX_MOTIFS) - _count_tuple_motifs(rw, _RBFOX_MOTIFS))
    nova_delta = float(_count_tuple_motifs(aw, _NOVA_MOTIFS) - _count_tuple_motifs(rw, _NOVA_MOTIFS))
    branch_loss = float(max(0, _count_tuple_motifs(rw, _BRANCH_MOTIFS) - _count_tuple_motifs(aw, _BRANCH_MOTIFS)))

    plo, phi = max(0, vpos - 10), min(len(ref_u), vpos + 11)
    ppt_ref = sum(1 for c in ref_u[plo:phi] if c in "CT") / max(1, phi - plo)
    ppt_alt = sum(1 for c in alt_u[plo:phi] if c in "CT") / max(1, phi - plo)
    ppt_delta = ppt_alt - ppt_ref

    pyr_loss = 1.0 if (vpos < len(ref_u) and ref_u[vpos] in "CT" and vpos < len(alt_u) and alt_u[vpos] in "AG") else 0.0
    gc_delta = float((alt_u.count("G") + alt_u.count("C")) - (ref_u.count("G") + ref_u.count("C")))
    rel_pos = abs(vpos - len(ref_u) / 2.0) / max(1.0, len(ref_u) / 2.0)
    is_gt_ag_core = 1.0 if (max_d_loss >= 6.0 or max_a_loss >= 6.0) else 0.0

    return np.array([
        max(0.0, max_d_loss), max(0.0, max_d_gain),
        max(0.0, max_a_loss), max(0.0, max_a_gain),
        max_ref_d, max_ref_a,
        ese_loss, ess_gain, tdp_delta, ptbp_delta, rbfox_delta, nova_delta, branch_loss,
        ppt_delta, pyr_loss, gc_delta, rel_pos, is_gt_ag_core,
    ], dtype=np.float32)


def auprc(y: np.ndarray, scores: np.ndarray) -> float:
    order = np.argsort(-scores)
    ys = y[order]
    tp = np.cumsum(ys)
    prec = tp / np.arange(1, len(ys) + 1)
    return float(prec[ys > 0.5].sum() / max(1, ys.sum()))


def auroc_ci(y: np.ndarray, scores: np.ndarray, n_boot: int = 1000, seed: int = 0) -> tuple[float, float, float]:
    pos = scores[y > 0.5]
    neg = scores[y < 0.5]
    u, _ = mannwhitneyu(pos, neg, alternative="two-sided")
    point = float(u) / (len(pos) * len(neg))
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        sp = pos[rng.integers(0, len(pos), len(pos))]
        sn = neg[rng.integers(0, len(neg), len(neg))]
        ub, _ = mannwhitneyu(sp, sn, alternative="two-sided")
        vals.append(float(ub) / (len(sp) * len(sn)))
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return point, float(lo), float(hi)


def fit_mlp(X_tr: np.ndarray, y_tr: np.ndarray, seed: int):
    rng = np.random.default_rng(seed)
    mu = X_tr.mean(axis=0)
    sd = np.where(X_tr.std(axis=0) < 1e-6, 1.0, X_tr.std(axis=0))
    Xs = (X_tr - mu) / sd
    d = Xs.shape[1]
    W1 = rng.normal(0, math.sqrt(2.0 / d), (d, HIDDEN)).astype(np.float32)
    b1 = np.zeros(HIDDEN, dtype=np.float32)
    W2 = rng.normal(0, math.sqrt(2.0 / HIDDEN), (HIDDEN, 1)).astype(np.float32)
    b2 = np.array([math.log(max(1e-4, y_tr.mean()) / max(1e-4, 1.0 - y_tr.mean()))], dtype=np.float32)
    params = [W1, b1, W2, b2]
    m_a = [np.zeros_like(p) for p in params]
    v_a = [np.zeros_like(p) for p in params]
    step = 0
    for _ in range(EPOCHS):
        perm = rng.permutation(len(Xs))
        for i in range(0, len(Xs), BATCH_SIZE):
            idx = perm[i : i + BATCH_SIZE]
            xb, yb = Xs[idx], y_tr[idx, None]
            z1 = xb @ W1 + b1
            a1 = np.maximum(0.0, z1)
            logit = np.clip(a1 @ W2 + b2, -20.0, 20.0)
            p = 1.0 / (1.0 + np.exp(-logit))
            err = (p - yb) / len(xb)
            gW2 = a1.T @ err + L2 * W2
            gb2 = err.sum(axis=0)
            dz1 = (err @ W2.T) * (z1 > 0)
            gW1 = xb.T @ dz1 + L2 * W1
            gb1 = dz1.sum(axis=0)
            step += 1
            for k, (pr, gr) in enumerate(zip(params, [gW1, gb1, gW2, gb2])):
                m_a[k] = 0.9 * m_a[k] + 0.1 * gr
                v_a[k] = 0.999 * v_a[k] + 0.001 * (gr * gr)
                pr -= LR * (m_a[k] / (1.0 - 0.9**step)) / (np.sqrt(v_a[k] / (1.0 - 0.999**step)) + 1e-8)
    return mu, sd, W1, b1, W2, b2


def main() -> None:
    if not DEFAULT_MFASS.exists():
        raise FileNotFoundError(f"Missing MFASS file: {DEFAULT_MFASS}")
    print(f"Dataset: {DEFAULT_MFASS}")
    print(f"SHA256 : {sha256_file(DEFAULT_MFASS)}")
    print(f"Config : seed={SEED}, hidden={HIDDEN}, epochs={EPOCHS}, batch={BATCH_SIZE}, lr={LR}, l2={L2}")

    X_rows, y_rows, exons = [], [], []
    with open(DEFAULT_MFASS, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["category"] == "natural":
                continue
            ref_u = r["nat_seq"].strip().upper()
            alt_u = r["original_seq"].strip().upper()
            if len(ref_u) != len(alt_u) or len(ref_u) < 20:
                continue
            X_rows.append(featurize_snv(ref_u, alt_u))
            y_rows.append(1.0 if r["strong_lof"].strip().upper() == "TRUE" else 0.0)
            exons.append(r["ensembl_id"])

    X = np.vstack(X_rows)
    y = np.array(y_rows, dtype=np.float32)
    exons_arr = np.array(exons)
    uniq_exons = np.unique(exons_arr)
    rng = np.random.default_rng(SEED)
    rng.shuffle(uniq_exons)
    exon_to_fold = {e: int(i % 5) for i, e in enumerate(uniq_exons)}
    folds = np.array([exon_to_fold[e] for e in exons_arr])

    SPLITS_OUT.parent.mkdir(parents=True, exist_ok=True)
    SPLITS_OUT.write_text(json.dumps({"seed": SEED, "exon_to_fold": exon_to_fold}), encoding="utf-8")

    oof = np.zeros(len(y), dtype=np.float32)
    for f in range(5):
        tr, te = folds != f, folds == f
        mu, sd, W1, b1, W2, b2 = fit_mlp(X[tr], y[tr], SEED + f)
        xs = (X[te] - mu) / sd
        a1 = np.maximum(0.0, xs @ W1 + b1)
        oof[te] = 1.0 / (1.0 + np.exp(-np.clip((a1 @ W2 + b2).ravel(), -20.0, 20.0)))

    auc_pt, auc_lo, auc_hi = auroc_ci(y, oof)
    ap = auprc(y, oof)
    print(f"5-Fold Exon-Grouped OOF (n={len(y)} SNVs, {int(y.sum())} pos / {len(y)-int(y.sum())} neg, {len(uniq_exons)} exons):")
    print(f"  AUROC = {auc_pt:.4f} (95% CI {auc_lo:.4f} - {auc_hi:.4f})")
    print(f"  AUPRC = {ap:.4f} (random baseline = {y.mean():.4f})")

    # Fit final production weights on full dataset with fixed seed
    mu, sd, W1, b1, W2, b2 = fit_mlp(X, y, SEED)
    WEIGHTS_OUT.write_text(json.dumps({
        "mu": [round(float(v), 6) for v in mu],
        "sd": [round(float(v), 6) for v in sd],
        "W1": [[round(float(v), 6) for v in row] for row in W1],
        "b1": [round(float(v), 6) for v in b1],
        "W2": [[round(float(v), 6) for v in row] for row in W2],
        "b2": [round(float(v), 6) for v in b2],
    }), encoding="utf-8")
    print(f"Saved clean weights (without baked metrics) to {WEIGHTS_OUT}")


if __name__ == "__main__":
    main()
