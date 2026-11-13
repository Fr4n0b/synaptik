"""Reproducible training and evaluation for SYNAPTIK ASO v3 MLP.

Reads ASO Atlas (`~/.synaptik/benchmarks/aso_atlas.pkl`) and Ensembl GRCh38 cDNA
(`~/.synaptik/data/Homo_sapiens.GRCh38.cdna.all.fa.gz`), logs SHA256 checksums,
hyperparameters, and the zero-leakage gene split to disk, trains the 162-feature
MLP (`162 -> 128 -> 64 -> 1`) with fixed seed=42, and reports held-out
neurogenetics Spearman rho and ROC-AUC (>=70% knockdown) with 95% bootstrap CIs.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import pickle
import re
import sys
from pathlib import Path

import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from synaptik.drug.aso import (
    ASOChemistry,
    compute_thermodynamics,
    count_cpg_motifs,
    predict_hairpin_delta_g,
    predict_self_dimer_score,
)
from synaptik.drug.aso_calibrated import (
    _DINUCS,
    _MOTIFS,
    _base_fracs,
    _max_run,
    compute_target_accessibility,
)

SEED = 42
HIDDEN1 = 128
HIDDEN2 = 64
EPOCHS = 30
BATCH_SIZE = 512
LR = 2e-3
L2 = 1e-3

NEURO_TEST_GENES = sorted({
    "SOD1", "MAPT", "HTT", "SNCA", "ATXN2", "APP",
    "LRRK2", "SMN1", "FUS", "C9ORF72", "BACE1", "APOE",
})

ASO_ATLAS_PATH = Path.home() / ".synaptik" / "benchmarks" / "aso_atlas.pkl"
FASTA_PATH = Path.home() / ".synaptik" / "data" / "Homo_sapiens.GRCh38.cdna.all.fa.gz"
SPLITS_OUT = Path.home() / ".synaptik" / "cache" / "aso_v3_splits.json"
WEIGHTS_OUT = REPO / "synaptik" / "drug" / "aso_mlp_v3.json"

_TRANS = str.maketrans("ACGTU", "TGCAA")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class _Stub:
    def __setstate__(self, s):
        if isinstance(s, dict):
            self.__dict__.update(s)


class _Unpickler(pickle.Unpickler):
    def find_class(self, mod, name):
        if mod.startswith("src"):
            return type(name, (_Stub,), {})
        return super().find_class(mod, name)


def _featurize(seq: str, target_site: str, flank_5p: str, flank_3p: str) -> np.ndarray:
    n = len(seq)
    gc = (seq.count("G") + seq.count("C")) / n
    _, _, dg37, tm = compute_thermodynamics(seq, chemistry=ASOChemistry.GAPMER_5_10_5)
    hp = predict_hairpin_delta_g(seq)
    sd_score = predict_self_dimer_score(seq)
    cpg, _ = count_cpg_motifs(seq)

    w5, w3 = seq[:5], seq[-5:]
    gap = seq[5:-5] if n > 10 else seq[2:-2]
    _, _, dg_w5, _ = compute_thermodynamics(w5, chemistry=ASOChemistry.MOE_PS)
    _, _, dg_w3, _ = compute_thermodynamics(w3, chemistry=ASOChemistry.MOE_PS)
    _, _, dg_gap, _ = compute_thermodynamics(gap, chemistry=ASOChemistry.GAPMER_5_10_5)

    bin_feats: list[float] = []
    for b in range(5):
        bin_feats.extend(_base_fracs(seq[(b * n) // 5 : ((b + 1) * n) // 5]))

    pos_onehot: list[float] = []
    for pos in list(range(10)) + list(range(-10, 0)):
        ch = seq[pos] if (0 <= pos < n) or (pos < 0 and abs(pos) <= n) else "N"
        pos_onehot.extend([1.0 if ch == b else 0.0 for b in "ACGT"])

    denom = max(1, n - 1)
    di_counts = {d: 0 for d in _DINUCS}
    for i in range(n - 1):
        d = seq[i : i + 2]
        if d in di_counts:
            di_counts[d] += 1
    di_freqs = [di_counts[d] / denom for d in _DINUCS]
    motif_hits = [float(seq.count(m)) for m in _MOTIFS]

    ent = sum(-f * math.log2(f) for f in _base_fracs(seq) if f > 0)
    pur_frac = (seq.count("A") + seq.count("G")) / n
    max_run = _max_run(seq)
    acc_feats = compute_target_accessibility(target_site, flank_5p, flank_3p)

    return np.array(
        [
            float(n), gc, gc * gc, dg37, tm, (tm - 62.0) ** 2,
            hp, sd_score, float(cpg),
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


def auroc_ci(y_true: np.ndarray, y_score: np.ndarray, n_boot: int = 1000, seed: int = 0) -> tuple[float, float, float]:
    pos = y_score[y_true >= 70.0]
    neg = y_score[y_true < 70.0]
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


def main() -> None:
    if not ASO_ATLAS_PATH.exists() or not FASTA_PATH.exists():
        raise FileNotFoundError("Missing ASO Atlas or Ensembl GRCh38 cDNA FASTA.")

    print(f"ASO Atlas : {ASO_ATLAS_PATH} (SHA256={sha256_file(ASO_ATLAS_PATH)})")
    print(f"Transcript: {FASTA_PATH} (SHA256={sha256_file(FASTA_PATH)})")
    print(f"Config    : seed={SEED}, arch=162->{HIDDEN1}->{HIDDEN2}->1, epochs={EPOCHS}, batch={BATCH_SIZE}, lr={LR}, l2={L2}")

    with open(ASO_ATLAS_PATH, "rb") as f:
        df = _Unpickler(f).load()
    df = df[df["target_gene"].notna() & df["aso_sequence_5_to_3"].notna() & df["inhibition_percent"].notna()].copy()
    df["gene"] = df["target_gene"].astype(str).str.strip().str.upper()
    df["seq"] = df["aso_sequence_5_to_3"].astype(str).str.strip().str.upper().str.replace("U", "T")
    df = df[~df["steric_blocking"].fillna(False).astype(bool)].copy()
    df = df[df["seq"].str.match(r"^[ACGT]{14,25}$")].copy()
    agg = df.groupby(["gene", "seq"], as_index=False)["inhibition_percent"].median()
    genes_needed = set(agg["gene"].unique())

    gene_txs: dict[str, list[str]] = {}
    gene_re = re.compile(r"gene_symbol:(\S+)")
    cur_gene, cur_seq = None, []
    with gzip.open(FASTA_PATH, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                if cur_gene and cur_gene in genes_needed and cur_seq:
                    gene_txs.setdefault(cur_gene, []).append("".join(cur_seq))
                m = gene_re.search(line)
                cur_gene = m.group(1).upper() if m else None
                cur_seq = []
            elif cur_gene and cur_gene in genes_needed:
                cur_seq.append(line.strip().upper())
        if cur_gene and cur_gene in genes_needed and cur_seq:
            gene_txs.setdefault(cur_gene, []).append("".join(cur_seq))
    for g in gene_txs:
        gene_txs[g].sort(key=len, reverse=True)

    X_list = []
    for g, s in zip(agg["gene"], agg["seq"]):
        ts = s.translate(_TRANS)[::-1]
        f5, f3 = "", ""
        for tx in gene_txs.get(g, ()):
            idx = tx.find(ts)
            if idx != -1:
                f5 = tx[max(0, idx - 50) : idx]
                f3 = tx[idx + len(ts) : idx + len(ts) + 50]
                break
        X_list.append(_featurize(s, ts, f5, f3))

    X = np.vstack(X_list)
    y = agg["inhibition_percent"].to_numpy(dtype=np.float32)
    genes = agg["gene"].to_numpy()

    is_test = np.isin(genes, NEURO_TEST_GENES)
    X_tr_all, y_tr_all, g_tr_all = X[~is_test], y[~is_test], genes[~is_test]
    X_te, y_te = X[is_test], y[is_test]

    rng = np.random.default_rng(SEED)
    uniq_tr = np.unique(g_tr_all)
    rng.shuffle(uniq_tr)
    val_genes = sorted(uniq_tr[: max(1, int(0.15 * len(uniq_tr)))].tolist())
    train_genes = sorted(set(uniq_tr.tolist()) - set(val_genes))
    is_val = np.isin(g_tr_all, val_genes)
    X_tr, y_tr = X_tr_all[~is_val], y_tr_all[~is_val]
    X_val, y_val = X_tr_all[is_val], y_tr_all[is_val]

    SPLITS_OUT.parent.mkdir(parents=True, exist_ok=True)
    SPLITS_OUT.write_text(json.dumps({
        "seed": SEED,
        "train_genes": train_genes,
        "val_genes": val_genes,
        "test_genes": NEURO_TEST_GENES,
    }, indent=2), encoding="utf-8")

    mu = X_tr.mean(axis=0)
    sd = np.where(X_tr.std(axis=0) < 1e-6, 1.0, X_tr.std(axis=0))
    X_tr_s, X_val_s, X_te_s = (X_tr - mu) / sd, (X_val - mu) / sd, (X_te - mu) / sd

    d_in = X_tr_s.shape[1]
    W1 = rng.normal(0, math.sqrt(2.0 / d_in), (d_in, HIDDEN1)).astype(np.float32)
    b1 = np.zeros(HIDDEN1, dtype=np.float32)
    W2 = rng.normal(0, math.sqrt(2.0 / HIDDEN1), (HIDDEN1, HIDDEN2)).astype(np.float32)
    b2 = np.zeros(HIDDEN2, dtype=np.float32)
    W3 = rng.normal(0, math.sqrt(2.0 / HIDDEN2), (HIDDEN2, 1)).astype(np.float32)
    b3 = np.array([float(y_tr.mean())], dtype=np.float32)
    params = [W1, b1, W2, b2, W3, b3]
    m_a = [np.zeros_like(p) for p in params]
    v_a = [np.zeros_like(p) for p in params]

    best_val_rho, best_w, step = -1.0, None, 0
    for epoch in range(1, EPOCHS + 1):
        perm = rng.permutation(len(X_tr_s))
        for i in range(0, len(X_tr_s), BATCH_SIZE):
            idx = perm[i : i + BATCH_SIZE]
            xb, yb = X_tr_s[idx], y_tr[idx, None]
            z1 = xb @ W1 + b1
            a1 = np.maximum(0.0, z1)
            z2 = a1 @ W2 + b2
            a2 = np.maximum(0.0, z2)
            err = (a2 @ W3 + b3 - yb) * (2.0 / len(xb))
            gW3, gb3 = a2.T @ err + L2 * W3, err.sum(axis=0)
            dz2 = (err @ W3.T) * (z2 > 0)
            gW2, gb2 = a1.T @ dz2 + L2 * W2, dz2.sum(axis=0)
            dz1 = (dz2 @ W2.T) * (z1 > 0)
            gW1, gb1 = xb.T @ dz1 + L2 * W1, dz1.sum(axis=0)
            step += 1
            for k, (pr, gr) in enumerate(zip(params, [gW1, gb1, gW2, gb2, gW3, gb3])):
                m_a[k] = 0.9 * m_a[k] + 0.1 * gr
                v_a[k] = 0.999 * v_a[k] + 0.001 * (gr * gr)
                pr -= LR * (m_a[k] / (1.0 - 0.9**step)) / (np.sqrt(v_a[k] / (1.0 - 0.999**step)) + 1e-8)

        va1 = np.maximum(0.0, X_val_s @ W1 + b1)
        vpred = (np.maximum(0.0, va1 @ W2 + b2) @ W3 + b3).ravel()
        vrho = float(spearmanr(y_val, vpred).statistic)
        if vrho > best_val_rho:
            best_val_rho = vrho
            best_w = [p.copy() for p in params]

    W1, b1, W2, b2, W3, b3 = best_w
    ta1 = np.maximum(0.0, X_te_s @ W1 + b1)
    tpred = (np.maximum(0.0, ta1 @ W2 + b2) @ W3 + b3).ravel()
    rho = float(spearmanr(y_te, tpred).statistic)
    auc_pt, auc_lo, auc_hi = auroc_ci(y_te, tpred)
    print(f"Held-out neurogenetics (n={len(y_te)}, >=70% pos={int((y_te>=70).sum())}):")
    print(f"  Spearman rho = {rho:.4f}")
    print(f"  ROC-AUC      = {auc_pt:.4f} (95% CI {auc_lo:.4f} - {auc_hi:.4f})")

    WEIGHTS_OUT.write_text(json.dumps({
        "model": "v3_accessibility_pos20_mlp",
        "n_features": int(d_in),
        "mu": [round(float(v), 6) for v in mu],
        "sd": [round(float(v), 6) for v in sd],
        "W1": [[round(float(v), 6) for v in row] for row in W1],
        "b1": [round(float(v), 6) for v in b1],
        "W2": [[round(float(v), 6) for v in row] for row in W2],
        "b2": [round(float(v), 6) for v in b2],
        "W3": [[round(float(v), 6) for v in row] for row in W3],
        "b3": [round(float(v), 6) for v in b3],
    }), encoding="utf-8")
    print(f"Saved clean weights (without baked metrics) to {WEIGHTS_OUT}")


if __name__ == "__main__":
    main()
