"""
Experiment design briefs.

Experiment lookup answers "has this been done?". A design brief answers the
next question: "given what has been done, how should I run mine?"

It condenses the closest prior experiments into:

  * the consensus protocol (maturity, co-culture, markers, replicate numbers)
  * the effect sizes that were reported, and what they imply for sample size
  * the weaknesses of the existing evidence that a new experiment could fix
  * a checklist of pitfalls specific to the readout

Two statistical positions are built in, because they are the two most common
ways an iPSC experiment ends up uninterpretable:

1. THE UNIT OF REPLICATION IS THE DIFFERENTIATION, NOT THE CELL.
   Three hundred neurons from one differentiation are one observation of the
   biology plus 299 observations of that dish. Sample sizes here are always in
   independent differentiations (or independent lines), and the brief says so.

2. PUBLISHED EFFECT SIZES ARE INFLATED.
   Small studies that reach significance over-estimate the true effect
   (the winner's curse). Powering a study on the published median therefore
   under-powers it. The brief reports the sample size for the published effect
   AND for a halved effect, and recommends planning on the latter.

The power calculation needs a variability estimate that papers almost never
report (between-differentiation CV). Rather than invent one, the brief shows
the answer across a range of plausible CVs, or uses the lab's own value when
supplied.
"""

from __future__ import annotations

import math
import statistics
from collections import Counter
from dataclasses import dataclass, field

from scipy import stats

from ..schema import AssayClass, Direction, ExperimentQuery, ORTHOGONALITY_GROUPS
from .lookup import LookupResult, MatchTier

DEFAULT_CVS = (0.15, 0.25, 0.35)

# Readout-specific pitfalls. Each is something that routinely invalidates the
# assay in iPSC-derived neurons, not generic advice.
ASSAY_PITFALLS: dict[AssayClass, list[str]] = {
    AssayClass.IMAGING_SYNAPTIC_PUNCTA: [
        "Normalise puncta to dendrite length (MAP2), not to field of view: "
        "neurite density differs between lines and will masquerade as a "
        "synaptic phenotype.",
        "Count colocalised pre- and postsynaptic puncta, not one marker "
        "alone; single-marker puncta include non-synaptic protein clusters.",
        "Fix acquisition settings and the thresholding method before "
        "unblinding, and image blind to genotype.",
        "Match plating density and astrocyte batch across genotypes; both "
        "change synapse number independently of the mutation.",
    ],
    AssayClass.ELECTROPHYS_PATCH: [
        "Report series resistance and exclusion criteria; cell selection is "
        "the largest source of bias in patch-clamp from mixed cultures.",
        "mEPSC frequency reflects both release probability and synapse "
        "number; pair it with an imaging readout before interpreting.",
        "Record across several differentiations, interleaving genotypes on "
        "the same day.",
    ],
    AssayClass.ELECTROPHYS_MEA: [
        "Network activity depends strongly on plating density and age; "
        "normalise to active electrodes and compare the same DIV.",
        "Wells are not independent of their plate: include plate as a factor "
        "and replicate across differentiations.",
    ],
    AssayClass.CALCIUM_IMAGING: [
        "Indicator expression level alters apparent transient amplitude; "
        "check it does not differ between genotypes.",
    ],
    AssayClass.SPLICING_ANALYSIS: [
        "Confirm mutant and control are matched for maturity (e.g. the "
        "PTBP1-to-PTBP2 switch): a maturation offset produces hundreds of "
        "splicing differences unrelated to the mutation.",
        "Validate key events by RT-PCR; short-read dPSI below ~0.1 is "
        "usually within noise.",
        "If the event is predicted to trigger NMD, measure protein as well "
        "as transcript, and consider an NMD-inhibited sample.",
    ],
    AssayClass.RNASEQ_BULK: [
        "Differences in neuronal purity between lines dominate bulk "
        "expression; quantify composition (marker genes) per sample.",
    ],
    AssayClass.WESTERN_BLOT: [
        "Normalise to total protein stain rather than a single housekeeping "
        "band, and stay inside the linear range.",
    ],
    AssayClass.VIABILITY: [
        "Distinguish cell death from failure to differentiate or attach; "
        "count neurons at plating as a baseline.",
    ],
    AssayClass.IMAGING_MORPHOLOGY: [
        "Trace blind to genotype, and analyse sparse-labelled cells so "
        "neurites can be assigned to a soma.",
    ],
}

GENERAL_PITFALLS = [
    "Use isogenic controls where the question is about a specific variant; "
    "between-donor variability otherwise exceeds most disease effects.",
    "Treat each independent differentiation as one replicate. Cells, fields "
    "and wells from the same differentiation are technical replicates.",
]


@dataclass
class PowerEstimate:
    effect_percent: float
    cv: float
    cohens_d: float
    n_per_group: int | None
    n_paired: int | None = None
    paired_rho: float = 0.7
    label: str = ""

    def describe(self) -> str:
        n = (f"{self.n_per_group} differentiations per group"
             if self.n_per_group else "more than 200 per group (impractical)")
        paired_str = f" [paired rho={self.paired_rho:.2g}: {self.n_paired} pairs]" if self.n_paired else ""
        return (f"{abs(self.effect_percent):.0f}% effect, CV "
                f"{self.cv:.0%}: {n}{paired_str}")


@dataclass
class DesignBrief:
    query: ExperimentQuery
    n_exact: int = 0
    n_close: int = 0
    n_labs: int = 0
    basis: str = ""

    # consensus protocol
    div_median: int | None = None
    div_range: tuple[int, int] | None = None
    differentiation_protocols: list[tuple[str, int]] = field(default_factory=list)
    cocultures: list[tuple[str, int]] = field(default_factory=list)
    markers: list[tuple[str, int]] = field(default_factory=list)
    statistical_tests: list[tuple[str, int]] = field(default_factory=list)
    n_replicates_median: int | None = None
    n_replicates_range: tuple[int, int] | None = None
    isogenic_fraction: float | None = None

    # effects
    consensus_direction: str = "unknown"
    direction_counts: dict[str, int] = field(default_factory=dict)
    effect_sizes_percent: list[float] = field(default_factory=list)
    median_effect_percent: float | None = None
    n_null_results: int = 0
    power: list[PowerEstimate] = field(default_factory=list)
    recommended_n: int | None = None
    recommended_n_paired: int | None = None
    recommended_n_basis: str = ""
    paired_rho: float = 0.7
    # (n per group, CV, smallest detectable % effect) for feasible designs
    detectable: list[tuple[int, float, float]] = field(default_factory=list)
    impractical: bool = False

    # judgement
    weaknesses: list[str] = field(default_factory=list)
    recommendations: list[str] = field(default_factory=list)
    pitfalls: list[str] = field(default_factory=list)
    orthogonal_assays: list[str] = field(default_factory=list)
    key_papers: list[str] = field(default_factory=list)
    value_statement: str = ""


def sample_size_two_group(
    cohens_d: float, alpha: float = 0.05, power: float = 0.8, max_n: int = 200
) -> int | None:
    """Per-group n for a two-sided two-sample t-test, by exact noncentral t.

    Returns None when more than `max_n` per group would be needed -- for a
    differentiation-level replicate that is not an experiment anyone will run,
    and saying so is more useful than printing 1,400.
    """
    d = abs(cohens_d)
    if d <= 0 or not math.isfinite(d):
        return None
    for n in range(2, max_n + 1):
        df = 2 * n - 2
        nc = d * math.sqrt(n / 2.0)
        t_crit = stats.t.ppf(1 - alpha / 2.0, df)
        achieved = (1 - stats.nct.cdf(t_crit, df, nc)
                    + stats.nct.cdf(-t_crit, df, nc))
        if achieved >= power:
            return n
    return None


def sample_size_paired(
    cohens_d: float, rho: float = 0.7, alpha: float = 0.05, power: float = 0.8, max_n: int = 100
) -> int | None:
    """Number of pairs for a two-sided paired t-test by exact noncentral t.

    In human iPSC modeling, differentiating mutant and isogenic control side by side
    in each batch (parallel differentiation rounds) removes shared batch-to-batch
    variability. rho is the inter-batch correlation (typically 0.6 to 0.85).
    """
    d = abs(cohens_d)
    if d <= 0 or not math.isfinite(d) or rho >= 1.0 or rho < 0:
        return None
    denom = math.sqrt(2.0 * (1.0 - rho))
    if denom <= 0:
        return None
    dz = d / denom
    for n in range(2, max_n + 1):
        df = n - 1
        nc = dz * math.sqrt(n)
        t_crit = stats.t.ppf(1 - alpha / 2.0, df)
        achieved = (1 - stats.nct.cdf(t_crit, df, nc)
                    + stats.nct.cdf(-t_crit, df, nc))
        if achieved >= power:
            return n
    return None


def detectable_effect_d(
    n_per_group: int, alpha: float = 0.05, power: float = 0.8
) -> float | None:
    """Smallest standardized effect detectable with n per group (bisection)."""
    if n_per_group < 2:
        return None
    df = 2 * n_per_group - 2
    t_crit = stats.t.ppf(1 - alpha / 2.0, df)

    def achieved(d: float) -> float:
        nc = d * math.sqrt(n_per_group / 2.0)
        p = 1 - stats.nct.cdf(t_crit, df, nc) + stats.nct.cdf(-t_crit, df, nc)
        # The noncentral t CDF loses precision at large noncentrality and can
        # return NaN; power there is 1 for every practical purpose.
        return 1.0 if math.isnan(p) else float(p)

    # d = 8 is far beyond any real effect; keeping the bracket tight also
    # keeps the noncentral t in its numerically stable range.
    lo, hi = 0.0, 8.0
    if achieved(hi) < power:
        return None
    for _ in range(50):
        mid = (lo + hi) / 2.0
        if achieved(mid) >= power:
            hi = mid
        else:
            lo = mid
    return hi


FEASIBLE_N = (3, 4, 6, 8)
PRACTICAL_LIMIT = 12  # differentiations per group most labs can actually run


def _top(counter: Counter, n: int = 5) -> list[tuple[str, int]]:
    return [(k, v) for k, v in counter.most_common(n) if k]


class DesignAdvisor:
    """Builds a design brief from the result of an experiment lookup."""

    def __init__(self, alpha: float = 0.05, power: float = 0.8):
        self.alpha = alpha
        self.power = power

    def brief(
        self, lookup: LookupResult, cv: float | None = None, rho: float = 0.7,
    ) -> DesignBrief:
        q = lookup.query
        b = DesignBrief(query=q)

        exact = lookup.by_tier(MatchTier.EXACT)
        close = lookup.by_tier(MatchTier.CLOSE)
        b.n_exact, b.n_close = len(exact), len(close)
        basis = exact + close
        if not basis:
            related = lookup.by_tier(MatchTier.RELATED)
            basis = related
            b.basis = (
                f"No exact or close prior experiment. Protocol details below "
                f"are drawn from {len(related)} RELATED experiment(s) and are "
                "weaker guidance: they differ from your design in system, "
                "maturity or readout."
                if related else
                "No comparable prior experiment in the corpus. There is no "
                "published protocol or effect size to anchor this design."
            )
        else:
            b.basis = (
                f"Based on {b.n_exact} exact and {b.n_close} close prior "
                "experiment(s).")

        recs = [m.record for m in basis]
        b.n_labs = len({r.source.lab_id for r in recs if r.source.lab_id})
        self._protocol(b, recs)
        self._effects(b, recs, cv, anchored=bool(exact or close), rho=rho)
        self._weaknesses(b, recs, lookup)
        self._pitfalls(b, recs)
        b.key_papers = [m.summary_line() for m in basis[:6]]
        b.value_statement = self._value(b, lookup.conflicting)
        return b

    # ---------------- protocol ----------------

    @staticmethod
    def _protocol(b: DesignBrief, recs: list) -> None:
        divs = [r.cell_model.days_in_vitro for r in recs
                if r.cell_model.days_in_vitro is not None]
        if divs:
            b.div_median = int(statistics.median(divs))
            b.div_range = (min(divs), max(divs))
        b.differentiation_protocols = _top(Counter(
            r.cell_model.differentiation_protocol for r in recs
            if r.cell_model.differentiation_protocol))
        b.cocultures = _top(Counter(
            r.cell_model.coculture for r in recs if r.cell_model.coculture))
        b.markers = _top(Counter(
            m.upper() for r in recs for m in r.outcome.markers_used), 8)
        b.statistical_tests = _top(Counter(
            r.outcome.statistical_test for r in recs
            if r.outcome.statistical_test), 4)
        ns = [r.outcome.n_biological_replicates for r in recs
              if r.outcome.n_biological_replicates is not None]
        if ns:
            b.n_replicates_median = int(statistics.median(ns))
            b.n_replicates_range = (min(ns), max(ns))
        iso = [r.cell_model.is_isogenic_pair for r in recs
               if r.cell_model.is_isogenic_pair is not None]
        if iso:
            b.isogenic_fraction = sum(1 for x in iso if x) / len(iso)

    # ---------------- effects and power ----------------

    def _effects(
        self, b: DesignBrief, recs: list, cv: float | None, anchored: bool, rho: float = 0.7,
    ) -> None:
        dirs = Counter(r.outcome.direction.value for r in recs
                       if r.outcome.direction != Direction.UNSPECIFIED)
        b.direction_counts = dict(dirs)
        b.n_null_results = dirs.get(Direction.NO_CHANGE.value, 0)
        directional = {k: v for k, v in dirs.items()
                       if k != Direction.NO_CHANGE.value}
        if directional:
            top, n = max(directional.items(), key=lambda kv: kv[1])
            share = n / sum(directional.values())
            b.consensus_direction = top if share >= 0.7 else "disputed"

        # Only percent-change effects are commensurable enough to pool.
        b.effect_sizes_percent = [
            r.outcome.effect_size for r in recs
            if r.outcome.effect_size is not None
            and r.outcome.effect_size_unit == "percent_change"
            and r.outcome.direction != Direction.NO_CHANGE
        ]
        if not b.effect_sizes_percent or not anchored:
            return
        b.median_effect_percent = statistics.median(
            abs(x) for x in b.effect_sizes_percent)

        published = b.median_effect_percent
        shrunken = published / 2.0
        cvs = (cv,) if cv else DEFAULT_CVS
        b.paired_rho = rho
        for label, eff in (("published median", published),
                           ("halved (winner's-curse allowance)", shrunken)):
            for c in cvs:
                d = (eff / 100.0) / c
                n_two = sample_size_two_group(d, self.alpha, self.power)
                n_pr = sample_size_paired(d, rho=rho, alpha=self.alpha, power=self.power)
                b.power.append(PowerEstimate(
                    effect_percent=eff, cv=c, cohens_d=round(d, 2),
                    n_per_group=n_two,
                    n_paired=n_pr,
                    paired_rho=rho,
                    label=label,
                ))

        ref_cv = cv if cv else 0.25
        d_shrunk = (shrunken / 100.0) / ref_cv
        b.recommended_n = sample_size_two_group(d_shrunk, self.alpha, self.power)
        b.recommended_n_paired = sample_size_paired(d_shrunk, rho=rho, alpha=self.alpha, power=self.power)
        b.recommended_n_basis = (
            f"{self.power:.0%} power, two-sided alpha {self.alpha}, to detect "
            f"half the published median effect ({shrunken:.0f}%) with a "
            f"between-differentiation CV of {ref_cv:.0%}"
            + ("" if cv else
               " (an assumption: supply your own CV for a real estimate)")
            + ".")
        b.impractical = (b.recommended_n is None
                         or b.recommended_n > PRACTICAL_LIMIT)
        for n in FEASIBLE_N:
            d = detectable_effect_d(n, self.alpha, self.power)
            if d is not None:
                b.detectable.append((n, ref_cv, round(d * ref_cv * 100.0, 1)))

    # ---------------- weaknesses ----------------

    @staticmethod
    def _weaknesses(b: DesignBrief, recs: list, lookup: LookupResult) -> None:
        if not recs:
            b.recommendations.append(
                "Nothing comparable is published. Run a small pilot (3 "
                "differentiations) to estimate variability before committing "
                "to a full design.")
            return

        if b.n_labs <= 1:
            b.weaknesses.append(
                "All prior evidence comes from one research group; it has not "
                "been independently replicated.")
            b.recommendations.append(
                "An independent replication is publishable and useful in "
                "itself. Follow the original protocol closely for the "
                "primary readout before adding variations.")
        if b.consensus_direction == "disputed" or lookup.conflicting:
            b.weaknesses.append(
                "Prior reports disagree on the direction of the effect ("
                + ", ".join(f"{v} {k.replace('_', ' ')}"
                            for k, v in sorted(b.direction_counts.items(),
                                               key=lambda kv: -kv[1]))
                + ").")
            b.recommendations.append(
                "Identify what differs between the conflicting reports "
                "(maturity, co-culture, controls) and vary exactly that "
                "factor within one experiment.")
        if b.n_null_results == 0 and len(recs) >= 4:
            b.weaknesses.append(
                f"No null result among {len(recs)} prior reports; the "
                "published effect sizes are probably inflated.")
        if b.isogenic_fraction is not None and b.isogenic_fraction < 0.5:
            b.weaknesses.append(
                f"Only {b.isogenic_fraction:.0%} of prior experiments used "
                "isogenic controls.")
            b.recommendations.append(
                "Using an isogenic pair would make your result stronger than "
                "most of the existing evidence.")
        if b.n_replicates_median is not None and b.n_replicates_median <= 3:
            b.weaknesses.append(
                f"Prior studies used a median of {b.n_replicates_median} "
                "biological replicates, too few to estimate the effect size "
                "reliably.")
        quantified = sum(1 for r in recs if r.outcome.is_quantified)
        if quantified < len(recs) / 2:
            b.weaknesses.append(
                f"Only {quantified} of {len(recs)} prior reports quantified "
                "the effect.")
        if b.impractical and b.median_effect_percent is not None:
            # Telling someone to run 37 differentiations is not advice. Say
            # what a feasible design can and cannot detect, and how to make
            # the experiment more sensitive instead of merely bigger.
            feasible = next((d for d in b.detectable if d[0] == 6), None)
            b.weaknesses.append(
                "An unpaired design powered for a realistic (halved) effect "
                "would need more differentiations than is practical"
                + (f" ({b.recommended_n} per group)." if b.recommended_n
                   else " (more than 200 per group)."))
            if b.recommended_n_paired and b.recommended_n_paired <= 15:
                b.recommendations.append(
                    f"RESCUE VIA PAIRING: Differentiating mutant and isogenic control side by side "
                    f"in each batch (paired design, rho={b.paired_rho:.2g}) cuts the required "
                    f"sample size down to {b.recommended_n_paired} pairs!")
            else:
                b.recommendations.append(
                    "Reduce variance rather than add replicates: differentiate "
                    "mutant and isogenic control side by side in each batch and "
                    "analyse them as PAIRS. Pairing removes batch-to-batch "
                    "variability, which is usually the largest component.")
            if feasible:
                b.recommendations.append(
                    f"With {feasible[0]} unpaired differentiations per group "
                    f"you can detect effects of about {feasible[2]:.0f}% or "
                    "larger; state that limit in advance so a null result is "
                    "interpretable.")
            b.recommendations.append(
                "Measure your own between-differentiation CV in a pilot and "
                "re-run this brief with --cv; the estimate above rests on an "
                "assumed value.")
        elif b.recommended_n and b.n_replicates_median and \
                b.recommended_n > b.n_replicates_median:
            b.recommendations.append(
                f"Plan for about {b.recommended_n} independent "
                f"differentiations per group rather than the "
                f"{b.n_replicates_median} typical of prior work.")

        used = {ORTHOGONALITY_GROUPS.get(r.assay_class, "other") for r in recs}
        candidates = {
            "morphological": "synaptic imaging",
            "functional_electrical": "electrophysiology (patch-clamp or MEA)",
            "proteomic": "protein-level measurement (western blot / proteomics)",
            "transcriptomic": "transcript-level measurement",
        }
        b.orthogonal_assays = [
            label for grp, label in candidates.items() if grp not in used]
        if len(used - {"other"}) <= 1 and b.orthogonal_assays:
            b.weaknesses.append(
                "Prior evidence uses a single assay family, so the finding "
                "could be specific to that measurement.")
            b.recommendations.append(
                "Add an orthogonal readout: "
                + "; ".join(b.orthogonal_assays[:2]) + ".")

    @staticmethod
    def _pitfalls(b: DesignBrief, recs: list) -> None:
        assay = b.query.assay_class
        if assay is None and recs:
            assay = Counter(r.assay_class for r in recs).most_common(1)[0][0]
        b.pitfalls = list(ASSAY_PITFALLS.get(assay, [])) + GENERAL_PITFALLS

    @staticmethod
    def _value(b: DesignBrief, conflicting: bool = False) -> str:
        if b.n_exact + b.n_close == 0:
            return ("Unexplored: no closely comparable experiment is in the "
                    "corpus. Verify against PubMed, then treat as novel.")
        if b.n_labs <= 1:
            return ("First independent replication: the existing evidence "
                    "comes from a single group.")
        if b.consensus_direction == "disputed" or conflicting:
            return ("Tie-breaker: the literature is split, and a well-"
                    "controlled experiment could resolve it.")
        return (
            f"Confirmatory: {b.n_labs} groups have reported this. A new "
            "experiment adds most by extending it (your variant, an "
            "orthogonal readout, a rescue) rather than repeating it.")
