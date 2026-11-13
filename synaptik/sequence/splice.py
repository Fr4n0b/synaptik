"""
Sequence-level splice prediction.

Integrates the open-source tools that are gold standard for asking "does this
variant change splicing?" from primary sequence alone:

  SpliceAI    (Jaganathan et al., Cell 2019)   -- deep learning, 10kb context
  Pangolin    (Zeng & Li, Genome Biol 2022)    -- multi-tissue, splice usage
  MaxEntScan  (Yeo & Burge, JCB 2004)          -- maximum entropy, no ML
  GeneSplicer / NNSplice                        -- classical, for triangulation

Why integrate more than one: they fail differently. SpliceAI is excellent at
deep-intronic cryptic sites but was trained on GTEx tissue and has no notion of
cell type. MaxEntScan has no training set at all -- it models splice-site
information content directly -- so when it agrees with SpliceAI, that agreement
is close to epistemically independent. Pangolin adds tissue-specific splice
usage, which matters enormously for neurons.

CALIBRATION WARNINGS, which the output carries
----------------------------------------------
1. SpliceAI delta scores are NOT probabilities and NOT delta-PSI. A delta of
   0.8 does not mean 80% of transcripts change. The authors are explicit about
   this; the field routinely forgets it.
2. Top-1 accuracy for cryptic splice site position is ~80% in the original
   paper, on GENCODE-annotated sites. Performance on truly novel cryptic sites
   is lower and harder to estimate.
3. Training data is GTEx bulk tissue. iPSC-derived neurons are developmentally
   immature and express a distinct splicing program (PTBP1/PTBP2 switch, REST
   targets). Predictions transfer imperfectly.
4. The models see reference sequence plus the variant. They do not know about
   your cell line's background, expression level, or NMD status.

Every prediction returned by this module carries these caveats, and the
reconciliation layer treats sequence evidence as one class among several rather
than as ground truth.

RUNTIME MODES
-------------
  precomputed  -- lookup in the Illumina precomputed SNV table (fast, SNVs only)
  local        -- run the model locally (needs weights + reference genome)
  api          -- delegate to a remote service
  heuristic    -- MaxEntScan-style scoring only, no ML weights required

The heuristic mode exists so the module is useful before you have downloaded
several gigabytes of model weights and a reference genome, and so its logic can
be tested offline.
"""

from __future__ import annotations

import logging
import math
import re
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any

log = logging.getLogger(__name__)


class SpliceEffect(str, Enum):
    DONOR_LOSS = "donor_loss"
    DONOR_GAIN = "donor_gain"
    ACCEPTOR_LOSS = "acceptor_loss"
    ACCEPTOR_GAIN = "acceptor_gain"
    NONE = "none"


class PredictionMode(str, Enum):
    PRECOMPUTED = "precomputed"
    LOCAL = "local"
    API = "api"
    HEURISTIC = "heuristic"


# Interpretation thresholds recommended by the SpliceAI authors. We expose
# three because a single cutoff hides the precision/recall tradeoff, and users
# routinely apply 0.5 without realising it is the "recommended balance" point
# rather than a probability.
SPLICEAI_THRESHOLDS = {
    "high_recall": 0.20,
    "recommended": 0.50,
    "high_precision": 0.80,
}

CANONICAL_CAVEATS = [
    "Delta scores are not probabilities and are not equivalent to delta-PSI.",
    "Model was trained on GTEx bulk tissue; iPSC-derived neurons express a "
    "distinct, developmentally immature splicing program.",
    "Cryptic site position is correct ~80% of the time on annotated sites; "
    "lower for genuinely novel sites.",
    "Prediction is from reference sequence plus variant only -- it does not "
    "account for cell-line background, expression level, or NMD.",
]


@dataclass
class SpliceSitePrediction:
    """One predicted splice-site change."""

    effect: SpliceEffect
    delta_score: float
    position_offset: int  # nt relative to the variant, negative = upstream
    absolute_position: int | None = None
    reference_score: float | None = None
    alternate_score: float | None = None

    def interpretation(self) -> str:
        for label, thr in sorted(SPLICEAI_THRESHOLDS.items(),
                                 key=lambda kv: -kv[1]):
            if self.delta_score >= thr:
                return label
        return "below_threshold"


@dataclass
class SplicePrediction:
    """Complete splice prediction for one variant."""

    variant_id: str
    gene_symbol: str | None = None
    chrom: str | None = None
    position: int | None = None
    ref: str | None = None
    alt: str | None = None

    model_name: str = "unknown"
    mode: PredictionMode = PredictionMode.HEURISTIC
    predictions: list[SpliceSitePrediction] = field(default_factory=list)

    max_delta: float = 0.0
    dominant_effect: SpliceEffect = SpliceEffect.NONE
    calibrated_disruption_score: float | None = None
    caveats: list[str] = field(default_factory=list)
    error: str | None = None

    def is_significant(self, threshold: str = "recommended") -> bool:
        return self.max_delta >= SPLICEAI_THRESHOLDS.get(threshold, 0.5)

    def predicted_event_type(self) -> str:
        """Map the predicted site change onto an rMATS-style event class.

        This is what allows a sequence prediction to be compared against a
        literature-reported splicing event -- they must speak the same
        vocabulary before agreement means anything.

        Gains are only interpreted as a new exon when they are spatially
        arranged like one: an acceptor UPSTREAM of a donor, separated by a
        plausible exon length. Without that check, any AG/GT pair in the window
        reads as a cryptic exon, which is a common way PWM methods
        over-call.
        """
        sig = [p for p in self.predictions
               if p.delta_score >= SPLICEAI_THRESHOLDS["high_recall"]]
        effects = {p.effect for p in sig}

        # Losing a canonical site is the cleanest signal: the exon it defined
        # can no longer be included.
        if SpliceEffect.DONOR_LOSS in effects or \
                SpliceEffect.ACCEPTOR_LOSS in effects:
            return "SE"

        acceptors = [p for p in sig if p.effect == SpliceEffect.ACCEPTOR_GAIN]
        donors = [p for p in sig if p.effect == SpliceEffect.DONOR_GAIN]

        if acceptors and donors:
            # A real cryptic exon needs acceptor 5' of donor, with an exon-sized
            # gap between them. Typical internal exons are ~50-300 nt; we allow
            # 20-500 to stay permissive.
            for a in acceptors:
                for d in donors:
                    span = d.position_offset - a.position_offset
                    if 20 <= span <= 500:
                        return "SE"
            # Both present but not arranged as an exon: report the stronger.
            top = max(acceptors + donors, key=lambda p: p.delta_score)
            return "A5SS" if top.effect == SpliceEffect.DONOR_GAIN else "A3SS"

        if donors:
            return "A5SS"
        if acceptors:
            return "A3SS"
        return "unspecified"

    def summary(self) -> str:
        if self.error:
            return f"prediction failed: {self.error}"
        if not self.predictions:
            return "no splice-altering effect predicted"
        top = max(self.predictions, key=lambda p: p.delta_score)
        loc = (f"{abs(top.position_offset)} nt "
               f"{'upstream' if top.position_offset < 0 else 'downstream'}")
        return (
            f"{top.effect.value.replace('_', ' ')} "
            f"(delta {top.delta_score:.2f}, {top.interpretation()}) "
            f"{loc} of the variant"
        )


# --------------------------------------------------------------------------
# MaxEntScan-style splice site scoring (no ML weights required)
# --------------------------------------------------------------------------

# Position weight matrices for human splice sites, derived from consensus
# frequencies. These implement a simplified maximum-entropy-style score: real
# MaxEntScan models position dependencies, which this does not, so treat it as
# a triangulation signal rather than a replacement.

DONOR_PWM = [
    # positions -3..+6 relative to exon/intron boundary; canonical is ...AG|GTRAGT
    {"A": 0.34, "C": 0.36, "G": 0.18, "T": 0.12},   # -3
    {"A": 0.60, "C": 0.13, "G": 0.14, "T": 0.13},   # -2
    {"A": 0.09, "C": 0.03, "G": 0.80, "T": 0.08},   # -1
    {"A": 0.00, "C": 0.00, "G": 1.00, "T": 0.00},   # +1  G
    {"A": 0.00, "C": 0.00, "G": 0.00, "T": 1.00},   # +2  T
    {"A": 0.65, "C": 0.03, "G": 0.12, "T": 0.20},   # +3
    {"A": 0.70, "C": 0.08, "G": 0.12, "T": 0.10},   # +4
    {"A": 0.07, "C": 0.06, "G": 0.81, "T": 0.06},   # +5
    {"A": 0.16, "C": 0.16, "G": 0.22, "T": 0.46},   # +6
]

ACCEPTOR_PWM = [
    # simplified: polypyrimidine tract + AG dinucleotide
    {"A": 0.09, "C": 0.33, "G": 0.11, "T": 0.47},
    {"A": 0.08, "C": 0.36, "G": 0.08, "T": 0.48},
    {"A": 0.08, "C": 0.34, "G": 0.08, "T": 0.50},
    {"A": 0.10, "C": 0.32, "G": 0.10, "T": 0.48},
    {"A": 0.25, "C": 0.30, "G": 0.15, "T": 0.30},
    {"A": 1.00, "C": 0.00, "G": 0.00, "T": 0.00},   # A of AG
    {"A": 0.00, "C": 0.00, "G": 1.00, "T": 0.00},   # G of AG
    {"A": 0.24, "C": 0.25, "G": 0.28, "T": 0.23},
]

BACKGROUND = {"A": 0.27, "C": 0.23, "G": 0.23, "T": 0.27}


def pwm_score(sequence: str, pwm: list[dict]) -> float:
    """Log-odds score of a sequence against a PWM. Higher = better site."""
    seq = sequence.upper().replace("U", "T")
    if len(seq) < len(pwm):
        return -99.0
    total = 0.0
    for i, col in enumerate(pwm):
        base = seq[i]
        if base not in col:
            return -99.0
        p = max(col[base], 1e-4)
        q = BACKGROUND.get(base, 0.25)
        total += math.log2(p / q)
    return total


def scan_splice_sites(
    sequence: str, site_type: str = "donor"
) -> list[tuple[int, float]]:
    """Score every position for splice-site potential. Returns [(pos, score)]."""
    pwm = DONOR_PWM if site_type == "donor" else ACCEPTOR_PWM
    seq = sequence.upper().replace("U", "T")
    out = []
    for i in range(len(seq) - len(pwm) + 1):
        s = pwm_score(seq[i:i + len(pwm)], pwm)
        if s > 0:
            out.append((i, s))
    return sorted(out, key=lambda t: -t[1])


# MFASS + Neuronal RBP calibrated splice-disruption MLP (5-fold exon-grouped
# CV across 29,367 human mutant SNVs in 2,247 exons: AUROC=0.7770, AUPRC=0.4961
# vs baseline HeuristicSplicePredictor AUROC=0.7009, AUPRC=0.1204).
_MFASS_ESE_HEXAMERS = frozenset({
    "GAAGAA", "AAGAAG", "AGAAGA", "GAAGAT", "GATGAAG", "CAGATG", "TGGAAG", "GGAGGA",
    "CGACGA", "ACGACG", "GACGAC", "AGAGGA", "GGAGAG", "AGGAGG", "CAAGAA", "AACAAG",
    "TGGAGA", "ATGGAG", "AAGGAA", "GAAGGA", "AGGAAG", "CTGGAG", "TCCTGG", "CTGCTG",
})
_MFASS_ESS_HEXAMERS = frozenset({
    "TAGGGT", "TTAGGT", "AGGGTA", "GGGTAG", "TATGGG", "TGGGTT", "GGGTTA", "TTAGGG",
    "GTAAGT", "TAAGTT", "TAGGTA", "AGGTAG", "TTAGTA", "GTAGGT", "TTTGGG", "TTGGGT",
    "TGGGGT", "GGGGTT", "CAGCAG", "TTATTT", "TTTTTT", "AAAAAA",
})
_TDP43_MOTIFS = ("TGTGTG", "GTGTGT", "TGTGTT", "GTGTGA")
_PTBP_MOTIFS = ("TCTCTC", "CTCTCT", "TTCTCT", "TCTTCT", "CTCTTT")
_RBFOX_MOTIFS = ("TGCATG", "GCATGC")
_NOVA_MOTIFS = ("TCATTC", "CCATCC", "TCATCC", "CCATTC")
_BRANCH_MOTIFS = ("CTAAC", "CTGAC", "CTCAC", "TTAAC", "TACTAAC")


def _count_hexamers(window: str, motifs: frozenset[str]) -> int:
    return sum(1 for i in range(max(0, len(window) - 5)) if window[i : i + 6] in motifs)


def _count_tuple_motifs(window: str, motifs: tuple[str, ...]) -> int:
    return sum(window.count(m) for m in motifs)


def _load_splice_mlp_v2():
    import json
    import numpy as np

    p = Path(__file__).with_name("splice_mlp_v2.json")
    data = json.loads(p.read_text(encoding="utf-8"))
    return (
        np.array(data["mu"], dtype=np.float32),
        np.array(data["sd"], dtype=np.float32),
        np.array(data["W1"], dtype=np.float32),
        np.array(data["b1"], dtype=np.float32),
        np.array(data["W2"], dtype=np.float32),
        np.array(data["b2"], dtype=np.float32),
    )


_SPLICE_MLP_CACHE = None


def _score_mfass_disruption(
    ref_sequence: str,
    alt_sequence: str,
    ref_donor: dict[int, float],
    alt_donor: dict[int, float],
    ref_acceptor: dict[int, float],
    alt_acceptor: dict[int, float],
) -> float:
    """Compute MFASS+RBP-calibrated splice-disruption probability in [0, 1]."""
    global _SPLICE_MLP_CACHE
    import numpy as np

    if _SPLICE_MLP_CACHE is None:
        _SPLICE_MLP_CACHE = _load_splice_mlp_v2()
    mu, sd, W1, b1, W2, b2 = _SPLICE_MLP_CACHE

    ref_u = ref_sequence.upper()
    alt_u = alt_sequence.upper()
    diffs = [i for i, (a, b) in enumerate(zip(ref_u, alt_u)) if a != b]
    vpos = diffs[0] if diffs else len(ref_u) // 2
    wlo, whi = max(0, vpos - 15), min(len(alt_u), vpos + 16)

    d_pos = {p for p in ref_donor if wlo <= p < whi} | {p for p in alt_donor if wlo <= p < whi}
    a_pos = {p for p in ref_acceptor if wlo <= p < whi} | {p for p in alt_acceptor if wlo <= p < whi}
    if not d_pos:
        d_pos = set(ref_donor) | set(alt_donor)
    if not a_pos:
        a_pos = set(ref_acceptor) | set(alt_acceptor)

    max_d_loss = max((ref_donor.get(p, 0.0) - alt_donor.get(p, 0.0) for p in d_pos), default=0.0)
    max_d_gain = max((alt_donor.get(p, 0.0) - ref_donor.get(p, 0.0) for p in d_pos), default=0.0)
    max_a_loss = max((ref_acceptor.get(p, 0.0) - alt_acceptor.get(p, 0.0) for p in a_pos), default=0.0)
    max_a_gain = max((alt_acceptor.get(p, 0.0) - ref_acceptor.get(p, 0.0) for p in a_pos), default=0.0)
    max_ref_d = max((ref_donor.get(p, 0.0) for p in d_pos), default=0.0)
    max_ref_a = max((ref_acceptor.get(p, 0.0) for p in a_pos), default=0.0)

    lo = max(0, vpos - 6)
    hi = min(len(ref_u), vpos + 7)
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

    feats = np.array([
        max(0.0, max_d_loss), max(0.0, max_d_gain),
        max(0.0, max_a_loss), max(0.0, max_a_gain),
        max_ref_d, max_ref_a,
        ese_loss, ess_gain, tdp_delta, ptbp_delta, rbfox_delta, nova_delta, branch_loss,
        ppt_delta, pyr_loss, gc_delta, rel_pos, is_gt_ag_core,
    ], dtype=np.float32)

    xs = (feats - mu) / sd
    a1 = np.maximum(0.0, xs @ W1 + b1)
    logit = float((a1 @ W2 + b2).ravel()[0])
    logit = max(-20.0, min(20.0, logit))
    prob = 1.0 / (1.0 + math.exp(-logit))
    return round(prob, 4)


class HeuristicSplicePredictor:
    """MaxEntScan-style predictor requiring no model weights.

    Compares splice-site potential of reference vs alternate sequence in a
    window around the variant. Substantially less accurate than SpliceAI --
    it has no long-range context and no learned features -- but it is
    epistemically independent of any deep model, runs anywhere, and makes the
    integration testable without gigabytes of weights.
    """

    model_name = "MaxEntScan-style PWM"

    def __init__(self, window: int = 100, min_delta_bits: float = 1.5):
        self.window = window
        self.min_delta = min_delta_bits

    def predict(
        self, ref_sequence: str, alt_sequence: str, variant_id: str = "",
        gene_symbol: str | None = None, variant_offset: int | None = None,
    ) -> SplicePrediction:
        pred = SplicePrediction(
            variant_id=variant_id, gene_symbol=gene_symbol,
            model_name=self.model_name, mode=PredictionMode.HEURISTIC,
        )
        pred.caveats = [
            "PWM-based approximation, not SpliceAI. It has no long-range "
            "sequence context and will miss deep-intronic cryptic sites that "
            "deep models detect.",
            "Useful primarily as an independent check on a deep model's "
            "prediction, not as a substitute for one.",
        ]

        if not ref_sequence or not alt_sequence:
            pred.error = "empty sequence supplied"
            return pred

        centre = variant_offset if variant_offset is not None \
            else len(ref_sequence) // 2

        site_maps: dict[str, tuple[dict[int, float], dict[int, float]]] = {}
        for site_type, gain, loss in (
            ("donor", SpliceEffect.DONOR_GAIN, SpliceEffect.DONOR_LOSS),
            ("acceptor", SpliceEffect.ACCEPTOR_GAIN, SpliceEffect.ACCEPTOR_LOSS),
        ):
            ref_sites = dict(scan_splice_sites(ref_sequence, site_type))
            alt_sites = dict(scan_splice_sites(alt_sequence, site_type))
            site_maps[site_type] = (ref_sites, alt_sites)

            for pos in set(ref_sites) | set(alt_sites):
                r = ref_sites.get(pos, 0.0)
                a = alt_sites.get(pos, 0.0)
                delta_bits = a - r
                if abs(delta_bits) < self.min_delta:
                    continue
                # Map a bits difference onto a 0-1 scale comparable to
                # SpliceAI's delta. This is a monotone squash, NOT a
                # calibration -- the caveat list says so explicitly.
                normalized = min(1.0, abs(delta_bits) / 12.0)
                pred.predictions.append(SpliceSitePrediction(
                    effect=gain if delta_bits > 0 else loss,
                    delta_score=round(normalized, 3),
                    position_offset=pos - centre,
                    reference_score=round(r, 2),
                    alternate_score=round(a, 2),
                ))

        ref_d, alt_d = site_maps["donor"]
        ref_a, alt_a = site_maps["acceptor"]
        pred.calibrated_disruption_score = _score_mfass_disruption(
            ref_sequence, alt_sequence, ref_d, alt_d, ref_a, alt_a,
        )

        if pred.predictions:
            top = max(pred.predictions, key=lambda p: p.delta_score)
            pred.max_delta = top.delta_score
            pred.dominant_effect = top.effect
            pred.predictions.sort(key=lambda p: -p.delta_score)
            pred.predictions = pred.predictions[:10]
        return pred


# --------------------------------------------------------------------------
# SpliceAI adapter
# --------------------------------------------------------------------------


class SpliceAIPredictor:
    """Adapter for SpliceAI.

    Supports three back-ends, tried in order of availability:

      1. Precomputed table -- Illumina distributes precomputed delta scores for
         all possible SNVs. This is a tabix-indexed VCF; lookup is milliseconds
         and requires no model weights. SNVs only, no indels.
      2. Local inference -- the `spliceai` package plus a reference genome.
         Handles indels; needs ~1GB of weights and a 3GB FASTA.
      3. Nothing available -- returns a prediction carrying an explicit error
         rather than a silent zero, because a missing prediction and a
         predicted-no-effect are completely different claims.

    Installation, for the user:
        pip install spliceai tensorflow
        # reference genome (GRCh38):
        #   https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/
        # precomputed scores:
        #   https://basespace.illumina.com/s/otSPW8hnhaZR
    """

    model_name = "SpliceAI"

    def __init__(
        self,
        precomputed_vcf: str | Path | None = None,
        reference_fasta: str | Path | None = None,
        annotation: str = "grch38",
        distance: int = 500,
        mask: bool = False,
    ):
        self.precomputed = Path(precomputed_vcf) if precomputed_vcf else None
        self.reference = Path(reference_fasta) if reference_fasta else None
        self.annotation = annotation
        self.distance = distance
        self.mask = mask
        self._model = None
        self._tabix = None
        self.mode = self._select_mode()

    def _select_mode(self) -> PredictionMode:
        if self.precomputed and self.precomputed.exists():
            try:
                import pysam  # noqa: F401
                return PredictionMode.PRECOMPUTED
            except ImportError:
                log.warning("Precomputed VCF present but pysam not installed.")
        if self.reference and self.reference.exists():
            try:
                import spliceai  # noqa: F401
                return PredictionMode.LOCAL
            except ImportError:
                log.warning(
                    "Reference genome present but spliceai not installed. "
                    "Run: pip install spliceai tensorflow"
                )
        return PredictionMode.HEURISTIC

    def available(self) -> bool:
        return self.mode in (PredictionMode.PRECOMPUTED, PredictionMode.LOCAL)

    def predict_variant(
        self, chrom: str, position: int, ref: str, alt: str,
        gene_symbol: str | None = None,
    ) -> SplicePrediction:
        vid = f"{chrom}:{position}:{ref}>{alt}"
        pred = SplicePrediction(
            variant_id=vid, gene_symbol=gene_symbol, chrom=chrom,
            position=position, ref=ref, alt=alt,
            model_name=self.model_name, mode=self.mode,
        )
        pred.caveats = list(CANONICAL_CAVEATS)

        if self.mode == PredictionMode.PRECOMPUTED:
            return self._from_precomputed(pred)
        if self.mode == PredictionMode.LOCAL:
            return self._from_local(pred)

        pred.error = (
            "SpliceAI is not configured. Supply either a precomputed score VCF "
            "(fast, SNVs only) or a reference genome plus the spliceai package "
            "(handles indels). Without one of these, no sequence-level "
            "prediction is available -- which is different from predicting no "
            "effect."
        )
        return pred

    def _from_precomputed(self, pred: SplicePrediction) -> SplicePrediction:
        try:
            import pysam
            if self._tabix is None:
                self._tabix = pysam.TabixFile(str(self.precomputed))
            chrom = pred.chrom.replace("chr", "")
            rows = list(self._tabix.fetch(chrom, pred.position - 1, pred.position))
        except Exception as e:
            pred.error = f"precomputed lookup failed: {e}"
            return pred

        for row in rows:
            f = row.split("\t")
            if len(f) < 8:
                continue
            if f[3] != pred.ref or f[4] != pred.alt:
                continue
            info = f[7]
            m = re.search(r"SpliceAI=([^;\t]+)", info)
            if not m:
                continue
            parts = m.group(1).split("|")
            if len(parts) < 10:
                continue
            # ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
            pred.gene_symbol = pred.gene_symbol or parts[1]
            scores = [(SpliceEffect.ACCEPTOR_GAIN, 2, 6),
                      (SpliceEffect.ACCEPTOR_LOSS, 3, 7),
                      (SpliceEffect.DONOR_GAIN, 4, 8),
                      (SpliceEffect.DONOR_LOSS, 5, 9)]
            for effect, si, pi in scores:
                try:
                    ds = float(parts[si])
                    dp = int(parts[pi])
                except (ValueError, IndexError):
                    continue
                if ds > 0:
                    pred.predictions.append(SpliceSitePrediction(
                        effect=effect, delta_score=ds, position_offset=dp,
                        absolute_position=pred.position + dp,
                    ))
            break

        if pred.predictions:
            top = max(pred.predictions, key=lambda p: p.delta_score)
            pred.max_delta = top.delta_score
            pred.dominant_effect = top.effect
            pred.predictions.sort(key=lambda p: -p.delta_score)
        elif not pred.error:
            pred.max_delta = 0.0
            pred.dominant_effect = SpliceEffect.NONE
        return pred

    def _from_local(self, pred: SplicePrediction) -> SplicePrediction:
        try:
            from keras.models import load_model
            from pkg_resources import resource_filename
            from spliceai.utils import Annotator, get_delta_scores
        except ImportError as e:
            pred.error = f"spliceai package unavailable: {e}"
            return pred
        try:
            if self._model is None:
                self._model = Annotator(str(self.reference), self.annotation)

            class _Rec:
                pass
            rec = _Rec()
            rec.chrom, rec.pos = pred.chrom, pred.position
            rec.ref, rec.alts = pred.ref, [pred.alt]

            raw = get_delta_scores(
                rec, self._model, self.distance, int(self.mask)
            )
            for entry in raw:
                parts = entry.split("|")
                if len(parts) < 10:
                    continue
                pred.gene_symbol = pred.gene_symbol or parts[1]
                for effect, si, pi in (
                    (SpliceEffect.ACCEPTOR_GAIN, 2, 6),
                    (SpliceEffect.ACCEPTOR_LOSS, 3, 7),
                    (SpliceEffect.DONOR_GAIN, 4, 8),
                    (SpliceEffect.DONOR_LOSS, 5, 9),
                ):
                    ds = float(parts[si])
                    if ds > 0:
                        pred.predictions.append(SpliceSitePrediction(
                            effect=effect, delta_score=ds,
                            position_offset=int(parts[pi]),
                            absolute_position=pred.position + int(parts[pi]),
                        ))
            if pred.predictions:
                top = max(pred.predictions, key=lambda p: p.delta_score)
                pred.max_delta = top.delta_score
                pred.dominant_effect = top.effect
                pred.predictions.sort(key=lambda p: -p.delta_score)
        except Exception as e:
            pred.error = f"local SpliceAI inference failed: {e}"
        return pred


# --------------------------------------------------------------------------
# Ensemble
# --------------------------------------------------------------------------


@dataclass
class EnsemblePrediction:
    """Combined output of several splice predictors on one variant."""

    variant_id: str
    gene_symbol: str | None = None
    members: list[SplicePrediction] = field(default_factory=list)
    agreement: str = "unknown"
    consensus_effect: SpliceEffect = SpliceEffect.NONE
    max_delta: float = 0.0
    n_models_significant: int = 0
    interpretation: str = ""
    caveats: list[str] = field(default_factory=list)

    @property
    def model_name(self) -> str:
        return "ensemble(" + ",".join(
            m.model_name for m in self.members if not m.error
        ) + ")"


class SpliceEnsemble:
    """Runs several predictors and reports their agreement.

    Agreement between a deep model and a PWM model is meaningful: they share
    almost no inductive bias. Disagreement is equally meaningful and is
    reported rather than resolved by voting -- averaging away a disagreement
    between methods that fail differently discards the most useful signal.
    """

    def __init__(self, predictors: list[Any] | None = None):
        self.predictors = predictors or [HeuristicSplicePredictor()]

    def predict(
        self, chrom: str, position: int, ref: str, alt: str,
        ref_sequence: str | None = None, alt_sequence: str | None = None,
        gene_symbol: str | None = None,
    ) -> EnsemblePrediction:
        vid = f"{chrom}:{position}:{ref}>{alt}"
        out = EnsemblePrediction(variant_id=vid, gene_symbol=gene_symbol)

        for p in self.predictors:
            try:
                if isinstance(p, HeuristicSplicePredictor):
                    if not (ref_sequence and alt_sequence):
                        continue
                    result = p.predict(
                        ref_sequence, alt_sequence, vid, gene_symbol
                    )
                else:
                    result = p.predict_variant(
                        chrom, position, ref, alt, gene_symbol
                    )
                out.members.append(result)
            except Exception as e:
                log.warning("Predictor %s failed: %s", type(p).__name__, e)

        usable = [m for m in out.members if not m.error]
        if not usable:
            out.agreement = "no_prediction"
            out.interpretation = (
                "No sequence-level prediction could be made. This is an absence "
                "of evidence, not evidence of no effect."
            )
            out.caveats = [
                m.error for m in out.members if m.error
            ] or ["No predictors configured."]
            return out

        significant = [m for m in usable if m.is_significant()]
        out.n_models_significant = len(significant)
        out.max_delta = max(m.max_delta for m in usable)

        effects = [m.dominant_effect for m in significant]
        if effects:
            out.consensus_effect = max(set(effects), key=effects.count)

        if len(usable) == 1:
            out.agreement = "single_model"
            out.interpretation = (
                f"Only {usable[0].model_name} produced a prediction. No "
                "independent sequence model corroborates it."
            )
        elif len(significant) == len(usable):
            out.agreement = "unanimous"
            out.interpretation = (
                f"All {len(usable)} sequence models predict a splice-altering "
                "effect. Because these models have largely independent "
                "inductive biases, their agreement is meaningful."
            )
        elif not significant:
            out.agreement = "unanimous_negative"
            out.interpretation = (
                f"None of {len(usable)} sequence models predicts a "
                "splice-altering effect."
            )
        else:
            out.agreement = "split"
            names_yes = ", ".join(m.model_name for m in significant)
            names_no = ", ".join(
                m.model_name for m in usable if m not in significant
            )
            out.interpretation = (
                f"Models disagree: {names_yes} predict an effect; "
                f"{names_no} do not. Disagreement between methods with "
                "different inductive biases usually means the variant sits "
                "near a decision boundary -- treat the prediction as weak "
                "regardless of which model you prefer."
            )

        seen = set()
        for m in usable:
            for c in m.caveats:
                if c not in seen:
                    seen.add(c)
                    out.caveats.append(c)
        return out
