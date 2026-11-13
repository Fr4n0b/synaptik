"""Reproducible temporal cross-validation for SYNAPTIK EvidenceWeighter.

Reads cached Europe PMC full-text corpus (`~/.synaptik/cache/epmc_oa_1k_papers.jsonl`),
logs SHA256 checksum and temporal cutoff splits, excludes raw publication volume
(`n_records`) and does NOT penalize `null_ratio`, and reports out-of-sample
Cutoff 2022 AUROC with 95% bootstrap CIs vs. `record_count`, `enrichment`, and
`p_value`. Refuses to declare a winner when confidence intervals overlap.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from synaptik.extract.fulltext import extract_records_from_fulltext
from synaptik.infer.benchmark import (
    RetrospectiveBenchmark,
    auroc,
    bootstrap_auroc_ci,
    bootstrap_auroc_difference,
)
from synaptik.normalize.labs import LabResolver
from synaptik.normalize.matcher import PhenotypeMatcher
from synaptik.schema import SourceMetadata

SEED = 42
L2_REG = 0.05
TRAIN_CUTOFFS = (2018, 2019, 2020, 2021)
TEST_CUTOFF = 2022
FEAT_NAMES = ("n_labs", "n_donor_lines", "n_clones", "is_isogenic_pair", "assay_diversity")
LOG_COLS = {"n_labs", "n_donor_lines", "n_clones", "assay_diversity"}

CORPUS_CACHE = Path.home() / ".synaptik" / "cache" / "epmc_oa_1k_papers.jsonl"
SPLITS_OUT = Path.home() / ".synaptik" / "cache" / "evidence_temporal_splits.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def feat_vec(d: dict[str, float]) -> list[float]:
    return [
        math.log1p(max(0.0, float(d.get(k, 0.0)))) if k in LOG_COLS else float(d.get(k, 0.0))
        for k in FEAT_NAMES
    ]


def main() -> None:
    if not CORPUS_CACHE.exists():
        raise FileNotFoundError(f"Missing cached OA corpus: {CORPUS_CACHE}")

    print(f"Corpus : {CORPUS_CACHE}")
    print(f"SHA256 : {sha256_file(CORPUS_CACHE)}")
    print(f"Config : seed={SEED}, l2={L2_REG}, features={FEAT_NAMES}, train_cutoffs={TRAIN_CUTOFFS}, test_cutoff={TEST_CUTOFF}")

    papers = []
    with open(CORPUS_CACHE, encoding="utf-8") as f:
        for line in f:
            if line.strip():
                obj = json.loads(line)
                papers.append((SourceMetadata.model_validate(obj["meta"]), obj["sections"]))

    lr = LabResolver()
    lr.fit([m for m, _ in papers])
    for m, _ in papers:
        m.lab_id = lr.assign(m)

    pm = PhenotypeMatcher()
    records = [r for m, s in papers for r in extract_records_from_fulltext(m, s, matcher=pm)]
    print(f"Loaded {len(papers)} OA papers -> {len(records)} ExperimentRecords")

    bench = RetrospectiveBenchmark(min_records_per_phenotype=3)
    folds = {cy: bench.run_gene_phenotype(records, cy, min_records=2) for cy in (*TRAIN_CUTOFFS, TEST_CUTOFF)}

    SPLITS_OUT.write_text(json.dumps({
        "seed": SEED,
        "train_cutoffs": list(TRAIN_CUTOFFS),
        "test_cutoff": TEST_CUTOFF,
        "fold_counts": {str(cy): len(folds[cy].cases) for cy in folds},
    }, indent=2), encoding="utf-8")

    X_tr = np.array([feat_vec(c.features) for cy in TRAIN_CUTOFFS for c in folds[cy].cases], dtype=np.float64)
    y_tr = np.array([1.0 if c.replicated else 0.0 for cy in TRAIN_CUTOFFS for c in folds[cy].cases], dtype=np.float64)

    test_cases = folds[TEST_CUTOFF].cases
    X_te = np.array([feat_vec(c.features) for c in test_cases], dtype=np.float64)
    y_te = [bool(c.replicated) for c in test_cases]

    def loss(p: np.ndarray) -> float:
        w, b = p[:-1], p[-1]
        z = np.clip(X_tr @ w + b, -20.0, 20.0)
        pr = 1.0 / (1.0 + np.exp(-z))
        return float(-np.mean(y_tr * np.log(pr + 1e-12) + (1.0 - y_tr) * np.log(1.0 - pr + 1e-12)) + L2_REG * np.sum(w * w))

    opt = minimize(loss, np.zeros(len(FEAT_NAMES) + 1), method="BFGS")
    w_fit, b_fit = opt.x[:-1], float(opt.x[-1])
    print("\nLearned rigor-only weights (no n_records, no null_ratio penalty):")
    for k, w in zip(FEAT_NAMES, w_fit):
        print(f"  {k:20s}: {w:+.4f}")
    print(f"  {'intercept':20s}: {b_fit:+.4f}")

    p_te = 1.0 / (1.0 + np.exp(-(X_te @ w_fit + b_fit)))
    n_pos = sum(y_te)
    n_neg = len(y_te) - n_pos
    print(f"\nOut-of-sample Cutoff {TEST_CUTOFF} (n={len(test_cases)}, {n_pos} pos / {n_neg} neg):")

    baselines = {
        "rigor_only": [float(v) for v in p_te],
        "record_count": [float(c.n_records) for c in test_cases],
        "enrichment": [float(c.enrichment) for c in test_cases],
        "p_value": [-math.log10(max(c.p_value, 1e-300)) for c in test_cases],
    }
    for name, sc in baselines.items():
        sp = [sc[i] for i, y in enumerate(y_te) if y]
        sn = [sc[i] for i, y in enumerate(y_te) if not y]
        val = auroc(sp, sn)
        ci = bootstrap_auroc_ci(sp, sn, seed=SEED)
        print(f"  {name:14s}: AUROC={val:.4f} (95% CI {ci[0]:.4f}-{ci[1]:.4f})")

    diff_ci = bootstrap_auroc_difference(y_te, baselines["rigor_only"], baselines["record_count"], seed=SEED)
    print(f"  diff (rigor_only - record_count): 95% CI {diff_ci[0]:+.4f} to {diff_ci[1]:+.4f}")
    if diff_ci[0] <= 0.0 <= diff_ci[1]:
        print("  Verdict: Confidence intervals overlap zero — no statistically resolvable winner.")


if __name__ == "__main__":
    main()
