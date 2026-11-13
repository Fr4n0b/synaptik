"""Empirically calibrated ASO gapmer knockdown potency model.

Fitted on 131,706 non-neurogenetics unique (gene, sequence) ASO measurements
across 327 human genes from ASO Atlas (5-fold gene-grouped CV selected L2
alpha=200.0, CV Spearman rho=0.2709) and evaluated on a strictly held-out
neurogenetics test set of 26,844 ASOs across 15 genes (SOD1, MAPT, HTT, SNCA,
ATXN2, APP, LRRK2, SMN1, FUS, C9ORF72, STMN2, BACE1, APOE, GRN, SCN1A):

  - Uncalibrated biophysical heuristic: Spearman rho = 0.1652, ROC-AUC = 0.5623
  - Calibrated biophysical model:     Spearman rho = 0.3268, ROC-AUC = 0.6670

Scope and caveats
-----------------
This score models sequence-intrinsic and thermodynamic knockdown potency in
vitro from ASO primary sequence alone. It does NOT account for target mRNA
secondary-structure accessibility, nuclear vs cytoplasmic localization,
delivery/conjugation, or transcriptome-wide off-target hybridization.
"""

from __future__ import annotations

_DINUCS = (
    "AA", "AC", "AG", "AT",
    "CA", "CC", "CG", "CT",
    "GA", "GC", "GG", "GT",
    "TA", "TC", "TG", "TT",
)

CALIBRATION_METADATA = {
    "dataset": "ASO Atlas (barneyhill/aso_atlas)",
    "n_train": 131706,
    "n_train_genes": 327,
    "n_test_neuro": 26844,
    "n_test_neuro_genes": 15,
    "cv_spearman_rho": 0.2709,
    "heldout_neuro_spearman_rho": 0.3268,
    "heldout_neuro_auroc_70pct": 0.6670,
    "baseline_neuro_spearman_rho": 0.1652,
    "baseline_neuro_auroc_70pct": 0.5623,
}

_MU = (
    18.165239, 0.472215, 0.242234, 0.263861, 0.263924, 0.219026, 0.253189,
    -23.051393, 59.138629, 115.475951, 0.001526, 0.072001, 0.144154, 0.210233,
    -2.044263, 3.073892, 0.359596, 0.483184, 0.500569, 0.456457, 0.083717,
    0.079359, 0.048927, 0.070111, 0.057542, 0.071651, 0.065161, 0.082288,
    0.071989, 0.021442, 0.075937, 0.056398, 0.057355, 0.054707, 0.050605,
    0.056053, 0.066564, 0.068832, 0.073365,
)

_SD = (
    1.993697, 0.138732, 0.136192, 0.111569, 0.117002, 0.11016, 0.11652,
    5.255346, 10.358017, 135.066932, 0.039036, 0.25849, 0.351246, 0.407474,
    1.600795, 0.803528, 0.651321, 0.242035, 0.239874, 0.172451, 0.276963,
    0.270298, 0.215716, 0.075854, 0.055023, 0.059536, 0.059645, 0.061467,
    0.078614, 0.039041, 0.060285, 0.055221, 0.056277, 0.061469, 0.052167,
    0.059347, 0.058209, 0.058581, 0.080371,
)

_INTERCEPT = 44.535659
_WEIGHTS = (
    -2.171427, 5.116264, -17.129301, -1.015785, -5.097857, 4.861157, 1.495752,
    -20.596454, -14.986423, -0.562786, -0.385045, -0.198934, -0.81645, -0.452737,
    1.8845, -0.49632, 0.072449, 0.96806, 3.096347, 0.597065, -0.144788,
    -0.359465, 0.285506, -4.710588, -0.929686, -2.931187, 0.570816, 1.398693,
    -2.56053, -1.663482, 0.704723, 0.882507, 2.380185, -0.858665, 4.257558,
    0.209922, 2.313501, -0.9573, 3.031786,
)


def predict_empirical_knockdown_pct(
    oligo_seq: str,
    dg37: float,
    tm: float,
    hairpin_dg: float,
    self_dimer: float,
    cpg_count: int,
) -> float:
    """Predict expected in-vitro knockdown (%) from ASO sequence & thermodynamics."""
    seq = oligo_seq.strip().upper().replace("U", "T")
    n = len(seq)
    if n < 8:
        return 0.0

    gc = (seq.count("G") + seq.count("C")) / n
    a_frac = seq.count("A") / n
    t_frac = seq.count("T") / n
    g_frac = seq.count("G") / n
    c_frac = seq.count("C") / n
    has_g4 = 1.0 if "GGGG" in seq else 0.0
    has_t4 = 1.0 if "TTTT" in seq else 0.0
    has_ggg = 1.0 if "GGG" in seq else 0.0
    has_ccc = 1.0 if "CCC" in seq else 0.0
    end5_gc = (seq[:4].count("G") + seq[:4].count("C")) / 4.0
    end3_gc = (seq[-4:].count("G") + seq[-4:].count("C")) / 4.0
    mid = seq[4:-4]
    mid_gc = (mid.count("G") + mid.count("C")) / max(1, len(mid))
    ccac = 1.0 if "CCAC" in seq else 0.0
    tccc = 1.0 if "TCCC" in seq else 0.0
    gact = 1.0 if "GACT" in seq else 0.0

    denom = max(1, n - 1)
    di_counts = {d: 0 for d in _DINUCS}
    for i in range(n - 1):
        d = seq[i : i + 2]
        if d in di_counts:
            di_counts[d] += 1
    di_freqs = [di_counts[d] / denom for d in _DINUCS]

    feats = [
        float(n), gc, gc * gc, a_frac, t_frac, g_frac, c_frac,
        dg37, tm, (tm - 62.0) ** 2,
        has_g4, has_t4, has_ggg, has_ccc, hairpin_dg, self_dimer, float(cpg_count),
        end5_gc, end3_gc, mid_gc, ccac, tccc, gact,
    ] + di_freqs

    pred = _INTERCEPT
    for val, mu, sd, w in zip(feats, _MU, _SD, _WEIGHTS):
        pred += ((val - mu) / sd) * w
    return round(max(0.0, min(100.0, pred)), 1)
