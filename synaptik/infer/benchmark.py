"""
Retrospective benchmark: do SYNAPTIK's confidence calls predict the future?

A tool that ranks evidence makes an implicit, testable claim: associations it
rates highly should hold up, and associations it flags should not. This module
tests that claim on the corpus itself, with a temporal hold-out.

  1. Freeze the corpus at a cutoff year. Compute every phenotype-pathway
     association using only records published up to that year.
  2. Look at what happened afterwards. An association counts as
     INDEPENDENTLY REPLICATED if, after the cutoff, a research group that had
     not contributed to it before reports the same pathway for the same
     phenotype in the same direction.
  3. Ask whether the pre-cutoff ranking predicted that outcome, and compare
     against the rankings a conventional tool would have produced
     (enrichment, p-value, raw paper count).

The headline number is the AUROC of each ranking for predicting later
independent replication. If evidence-structure ranking does not beat
enrichment or p-value on a real corpus, the tool's central premise is wrong
for that corpus, and this module is how you would find out.

What this benchmark does NOT show
---------------------------------
* Replication in the literature is not truth. A popular wrong idea gets
  "replicated"; an unfashionable right one does not get tested. The outcome
  measures whether a claim was independently re-observed, which is the best
  proxy available from literature alone, and no more than that.
* Attention is a confounder: heavily studied pathways are more likely to be
  reported again regardless of evidence quality. The raw-count baseline is
  included precisely to expose this -- a ranking only adds value to the
  extent it beats "more papers predicts more papers".
* Associations whose phenotype has no post-cutoff records at all are
  unobservable and are excluded, not counted as failures.
* On the synthetic demo corpus the numbers demonstrate the procedure only.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, field

from scipy import stats

from ..schema import Direction, ExperimentRecord
from .association import AssociationEngine, PathwayAssociation
from .weighting import EvidenceWeighter

MIN_CLASS_SIZE = 10  # below this per class, an AUROC is anecdote


@dataclass
class BenchmarkCase:
    phenotype_key: str
    pathway_label: str
    tier: str
    rank_score: float
    enrichment: float
    p_value: float
    n_records: int
    n_labs_before: int
    replicated: bool
    contradicted: bool
    n_new_labs: int


@dataclass
class RankerScore:
    name: str
    auroc: float | None
    description: str = ""


@dataclass
class BenchmarkResult:
    cutoff_year: int
    n_records_before: int = 0
    n_records_after: int = 0
    n_phenotypes: int = 0
    cases: list[BenchmarkCase] = field(default_factory=list)
    n_replicated: int = 0
    rankers: list[RankerScore] = field(default_factory=list)
    tier_table: list[dict] = field(default_factory=list)
    underpowered: bool = True
    best_baseline: str | None = None
    difference_vs_best_baseline: float | None = None
    difference_ci95: tuple[float, float] | None = None
    verdict: str = ""
    caveats: list[str] = field(default_factory=list)

    def auroc(self, name: str) -> float | None:
        for r in self.rankers:
            if r.name == name:
                return r.auroc
        return None


def auroc(scores_pos: list[float], scores_neg: list[float]) -> float | None:
    """Area under the ROC curve via the Mann-Whitney U statistic.

    Ties count half. Returns None when either class is empty, because an
    AUROC is undefined there and 0.5 would be a fabricated answer.
    """
    if not scores_pos or not scores_neg:
        return None
    u, _ = stats.mannwhitneyu(scores_pos, scores_neg, alternative="two-sided")
    return float(u) / (len(scores_pos) * len(scores_neg))


def bootstrap_auroc_difference(
    labels: list[bool], a: list[float], b: list[float],
    n_boot: int = 2000, seed: int = 0,
) -> tuple[float, float] | None:
    """95% bootstrap interval for AUROC(a) - AUROC(b), resampling cases.

    Paired: both rankings are scored on the same resampled associations, so
    the interval reflects their difference rather than their separate noise.
    Deterministic for a given seed, so a benchmark run is reproducible.
    """
    import numpy as np

    y = np.asarray(labels, dtype=bool)
    xa, xb = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    n = len(y)
    if n == 0 or y.all() or not y.any():
        return None

    def auc(scores, yy):
        n_pos = int(yy.sum())
        n_neg = len(yy) - n_pos
        if n_pos == 0 or n_neg == 0:
            return None
        ranks = stats.rankdata(scores)
        return (ranks[yy].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)

    rng = np.random.default_rng(seed)
    diffs = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        ua, ub = auc(xa[idx], y[idx]), auc(xb[idx], y[idx])
        if ua is not None and ub is not None:
            diffs.append(ua - ub)
    if len(diffs) < n_boot // 2:
        return None
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    return (float(lo), float(hi))


class RetrospectiveBenchmark:
    def __init__(
        self,
        weighter: EvidenceWeighter | None = None,
        min_records_per_phenotype: int = 5,
    ):
        self.min_records = min_records_per_phenotype
        self._weighter = weighter

    def run(
        self, records: list[ExperimentRecord], cutoff_year: int
    ) -> BenchmarkResult:
        res = BenchmarkResult(cutoff_year=cutoff_year)
        dated = [r for r in records if r.source.year]
        before = [r for r in dated if r.source.year <= cutoff_year]
        after = [r for r in dated if r.source.year > cutoff_year]
        res.n_records_before, res.n_records_after = len(before), len(after)

        undated = len(records) - len(dated)
        if undated:
            res.caveats.append(
                f"{undated} record(s) have no publication year and were "
                "excluded.")
        if not before or not after:
            res.verdict = (
                "Cannot run: the cutoff leaves no records on one side. Choose "
                "a cutoff year inside the range covered by the corpus.")
            return res

        # Weight as a reader AT the cutoff would have: recency relative to
        # the cutoff, not to today. Retraction status is the one unavoidable
        # leak -- a paper retracted after the cutoff is already flagged.
        weighter = self._weighter or EvidenceWeighter(current_year=cutoff_year)
        engine = AssociationEngine(weighter)

        before_by_ph: dict[str, list[ExperimentRecord]] = defaultdict(list)
        after_by_ph: dict[str, list[ExperimentRecord]] = defaultdict(list)
        for r in before:
            before_by_ph[r.outcome.phenotype.key()].append(r)
        for r in after:
            after_by_ph[r.outcome.phenotype.key()].append(r)

        for ph, fg in before_by_ph.items():
            if len(fg) < self.min_records:
                continue
            future = after_by_ph.get(ph)
            if not future:
                continue  # unobservable, not a failure
            res.n_phenotypes += 1
            for assoc in engine.associate(fg, before):
                res.cases.append(self._evaluate(ph, assoc, fg, future))

        self._score(res)
        return res

    # ---------------- outcome ----------------

    @staticmethod
    def _evaluate(
        ph: str, assoc: PathwayAssociation,
        fg: list[ExperimentRecord], future: list[ExperimentRecord],
    ) -> BenchmarkCase:
        prior_labs = {
            r.source.lab_id for r in fg
            if r.source.lab_id
            and any(p.key() == assoc.pathway_key for p in r.pathways)
        }
        dominant = assoc.dominant_direction
        directional = dominant in (Direction.INCREASE.value,
                                   Direction.DECREASE.value)

        confirming, contradicting = set(), set()
        for r in future:
            lab = r.source.lab_id
            if not lab or lab in prior_labs:
                continue  # the same group again is not independent
            if r.source.is_retracted:
                continue
            if not any(p.key() == assoc.pathway_key for p in r.pathways):
                continue
            d = r.outcome.direction.value
            if r.is_negative_result():
                contradicting.add(lab)
            elif not directional or d == dominant:
                confirming.add(lab)
            elif d in (Direction.INCREASE.value, Direction.DECREASE.value):
                contradicting.add(lab)

        return BenchmarkCase(
            phenotype_key=ph, pathway_label=assoc.pathway_label,
            tier=assoc.confidence_tier(), rank_score=assoc.rank_score(),
            enrichment=assoc.enrichment, p_value=assoc.p_fisher,
            n_records=assoc.n_records, n_labs_before=assoc.n_labs,
            replicated=bool(confirming),
            contradicted=bool(contradicting) and not confirming,
            n_new_labs=len(confirming),
        )

    # ---------------- scoring ----------------

    @staticmethod
    def _score(res: BenchmarkResult) -> None:
        cases = res.cases
        if not cases:
            res.verdict = (
                "No evaluable associations: no phenotype has enough records "
                "before the cutoff together with any records after it.")
            return
        pos = [c for c in cases if c.replicated]
        neg = [c for c in cases if not c.replicated]
        res.n_replicated = len(pos)
        res.underpowered = min(len(pos), len(neg)) < MIN_CLASS_SIZE

        def ranker(name, fn, description):
            return RankerScore(
                name, auroc([fn(c) for c in pos], [fn(c) for c in neg]),
                description)

        res.rankers = [
            ranker("synaptik", lambda c: c.rank_score,
                   "evidence-structure rank score"),
            ranker("enrichment", lambda c: c.enrichment,
                   "odds ratio alone"),
            ranker("p_value", lambda c: -math.log10(max(c.p_value, 1e-300)),
                   "Fisher p-value alone"),
            ranker("record_count", lambda c: float(c.n_records),
                   "number of prior records (attention baseline)"),
        ]

        order = ["strong", "moderate", "weak", "single-lab", "insufficient"]
        by_tier: dict[str, list[BenchmarkCase]] = defaultdict(list)
        for c in cases:
            by_tier[c.tier].append(c)
        for t in order:
            group = by_tier.get(t, [])
            if not group:
                continue
            res.tier_table.append({
                "tier": t, "n": len(group),
                "replicated": sum(1 for c in group if c.replicated),
                "replication_rate": sum(1 for c in group if c.replicated)
                / len(group),
                "contradicted": sum(1 for c in group if c.contradicted),
            })

        res.caveats += [
            "Independent re-observation in the literature is a proxy for "
            "truth, not truth: popular ideas are re-reported whether or not "
            "they are right.",
            "Compare against the record-count baseline. A ranking is only "
            "informative to the extent it beats 'more papers predicts more "
            "papers'.",
        ]
        if res.underpowered:
            res.caveats.insert(0, (
                f"Underpowered: {len(pos)} replicated vs {len(neg)} "
                f"non-replicated associations. With fewer than "
                f"{MIN_CLASS_SIZE} in either class the AUROCs below are "
                "illustrative only."))

        s = res.auroc("synaptik")
        base = {n: res.auroc(n) for n in ("enrichment", "p_value",
                                          "record_count")}
        usable = {k: v for k, v in base.items() if v is not None}
        if s is None or not usable:
            res.verdict = (
                "AUROC undefined: every evaluable association had the same "
                "outcome.")
            return
        best_name, best = max(usable.items(), key=lambda kv: kv[1])
        fns = {
            "enrichment": lambda c: c.enrichment,
            "p_value": lambda c: -math.log10(max(c.p_value, 1e-300)),
            "record_count": lambda c: float(c.n_records),
        }
        diff = s - best
        ci = bootstrap_auroc_difference(
            [c.replicated for c in cases],
            [c.rank_score for c in cases],
            [fns[best_name](c) for c in cases],
        )
        res.best_baseline = best_name
        res.difference_vs_best_baseline = round(diff, 3)
        res.difference_ci95 = ci
        ci_txt = (f"difference {diff:+.2f}, 95% CI {ci[0]:+.2f} to "
                  f"{ci[1]:+.2f}" if ci else f"difference {diff:+.2f}")

        # A verdict of "better" or "worse" is only issued when the bootstrap
        # interval excludes zero. Anything else is reported as no detectable
        # difference, however the point estimates happen to fall.
        if ci and ci[0] > 0:
            res.verdict = (
                f"Evidence-structure ranking predicted later independent "
                f"replication better than the best conventional baseline "
                f"({best_name}): AUROC {s:.2f} vs {best:.2f} ({ci_txt}).")
        elif ci and ci[1] < 0:
            res.verdict = (
                f"Evidence-structure ranking predicted later replication "
                f"WORSE than the {best_name} baseline: AUROC {s:.2f} vs "
                f"{best:.2f} ({ci_txt}). The tiers should not be trusted on "
                "this corpus until that is understood.")
        else:
            res.verdict = (
                f"No detectable difference between evidence-structure "
                f"ranking and the best baseline ({best_name}): AUROC "
                f"{s:.2f} vs {best:.2f} ({ci_txt}). "
                + ("The corpus is too small to tell them apart."
                   if res.underpowered else
                   "On this corpus the ranking adds nothing measurable over "
                   "that baseline."))
