"""
Phenotype -> pathway association inference.

This is the analytical core. The central claim it makes is narrow and honest:

    "Among published experiments that measured phenotype P, pathway W is
     mentioned/implicated more often than chance, supported by N independent
     labs across K orthogonal assay families, with the following robustness
     properties."

It does NOT claim that W causes P. Literature co-occurrence cannot establish
that, and any tool that reports a p-value implying it is lying to the user.

What we compute, and why each piece is necessary:

  enrichment      -- weighted odds ratio of pathway W among P-records vs the
                     background corpus. The headline number.
  p_fisher        -- Fisher exact test on the unweighted contingency table,
                     with Benjamini-Hochberg FDR across all tested pathways.
  n_labs          -- distinct research groups. The honest replication count.
  n_ortho_groups  -- distinct assay families. Guards against a result that is
                     an artifact of one measurement technique.
  leave_one_lab_out -- does the association survive removing the single most
                     influential lab? If not, it is one group's finding.
  direction_consistency -- do the reports agree on the sign of the effect?
                     Conflicting directions are informative, not noise.
  negative_evidence -- count of explicit null results. Papers that looked and
                     found nothing are the strongest available correction for
                     publication bias.
  publication_bias_flag -- heuristics for a literature that has only ever
                     reported positive findings.

Every returned association carries these fields. The UI surfaces them all,
because a single score would let a user mistake "well studied" for "well
established" -- which is the exact error this tool exists to prevent.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from dataclasses import dataclass, field

from scipy import stats

from ..schema import Direction, ExperimentRecord
from .weighting import EvidenceWeighter


@dataclass
class SupportingRecord:
    """One record backing an association, with its weight and provenance."""

    record_id: str
    source_id: str
    title: str
    journal: str | None
    year: int | None
    lab_id: str | None
    direction: str
    effect_size: float | None
    p_value: float | None
    assay_class: str
    ortho_group: str
    weight: float
    quote: str | None
    is_negative: bool


@dataclass
class PathwayAssociation:
    pathway_key: str
    pathway_label: str

    # headline statistics
    enrichment: float = 1.0
    weighted_support: float = 0.0
    p_fisher: float = 1.0
    q_value: float = 1.0

    # contingency counts
    n_records: int = 0
    n_papers: int = 0
    n_labs: int = 0
    n_ortho_groups: int = 0
    n_background_records: int = 0

    # robustness
    leave_one_lab_out_enrichment: float | None = None
    lab_concentration: float = 0.0  # Herfindahl index over labs, 0..1
    survives_lolo: bool = False

    # direction
    direction_counts: dict[str, int] = field(default_factory=dict)
    direction_consistency: float = 0.0
    dominant_direction: str = "unspecified"

    # negative evidence / bias
    n_negative_records: int = 0
    negative_fraction: float = 0.0
    publication_bias_flag: str = ""

    # effect sizes where reported
    median_effect_size: float | None = None
    effect_size_iqr: tuple[float, float] | None = None
    n_with_effect_size: int = 0

    # temporal
    first_year: int | None = None
    last_year: int | None = None

    supporting: list[SupportingRecord] = field(default_factory=list)

    def rank_score(self) -> float:
        """Ranking score. Enrichment and support must BOTH be present.

        Sorting by weighted support alone promotes ubiquitous pathways
        (apoptosis, cell signaling) that appear in every paper regardless of
        phenotype -- they accumulate support without being specific. Sorting by
        enrichment alone promotes noise from pathways with two records.

        We therefore take log-enrichment (specificity) times a saturating
        function of weighted support (evidence mass), and hard-gate on
        enrichment > 1: a depleted pathway is never "top-associated" no matter
        how much has been written about it.
        """
        if self.enrichment <= 1.0:
            return 0.0
        specificity = math.log2(self.enrichment)
        mass = math.log1p(max(0.0, self.weighted_support))
        # Independent replication is the scarcest signal; reward it explicitly.
        replication = 1.0 + 0.5 * math.log1p(max(0, self.n_labs - 1))
        orthogonality = 1.0 + 0.25 * max(0, self.n_ortho_groups - 1)
        return specificity * mass * replication * orthogonality

    def confidence_tier(self) -> str:
        """A blunt summary label. Deliberately conservative.

        The thresholds encode a view: nothing from a single lab is 'strong',
        no matter how many papers or how good the journal.
        """
        # Weighted support gates the upper tiers. Without this, evidence that
        # has been zeroed out (retracted) or heavily discounted (review
        # assertions, unquantified) can still reach 'moderate' on lab count
        # alone -- the counts say two groups, the weights say almost nothing.
        if self.weighted_support < 1.0:
            return "insufficient"
        if self.n_labs >= 3 and self.n_ortho_groups >= 2 and self.survives_lolo \
                and self.q_value < 0.05 and self.direction_consistency >= 0.7 \
                and self.weighted_support >= 3.0:
            return "strong"
        if self.n_labs >= 2 and self.q_value < 0.1 \
                and self.direction_consistency >= 0.7 \
                and self.weighted_support >= 2.0:
            return "moderate"
        if self.n_labs >= 2:
            return "weak"
        return "single-lab"

    def caveats(self) -> list[str]:
        """Human-readable warnings. Shown next to every result."""
        out = []
        if self.n_labs < 2:
            out.append(
                "All evidence comes from one research group -- not independently replicated."
            )
        if self.n_ortho_groups < 2:
            out.append(
                "All evidence uses one assay family -- may reflect a methodological artifact."
            )
        if not self.survives_lolo and self.n_labs >= 2:
            out.append(
                "Association disappears when the most prolific lab is excluded."
            )
        if self.direction_consistency < 0.7 and self.n_records >= 3:
            out.append(
                f"Reported directions are split ({self.direction_counts}) -- "
                "the literature does not agree on the sign of this effect."
            )
        if self.publication_bias_flag:
            out.append(self.publication_bias_flag)
        if self.n_with_effect_size == 0:
            out.append("No quantitative effect sizes reported in any supporting study.")
        if self.lab_concentration > 0.6:
            out.append(
                f"Evidence is concentrated in few labs (concentration index "
                f"{self.lab_concentration:.2f})."
            )
        return out


class AssociationEngine:
    """Computes weighted phenotype->pathway associations over a record set."""

    def __init__(
        self,
        weighter: EvidenceWeighter | None = None,
        min_records: int = 2,
        include_negative_in_support: bool = False,
    ):
        self.weighter = weighter or EvidenceWeighter()
        self.min_records = min_records
        self.include_negative = include_negative_in_support

    # ---------------- main entry ----------------

    def associate(
        self,
        foreground: list[ExperimentRecord],
        background: list[ExperimentRecord],
        direction_filter: Direction | None = None,
    ) -> list[PathwayAssociation]:
        """
        foreground: records measuring the phenotype of interest
        background: the rest of the corpus, defining the null expectation
        direction_filter: restrict to records reporting this direction

        Returning the background-conditioned enrichment rather than a raw count
        is what stops ubiquitously-mentioned pathways (mTOR, apoptosis) from
        topping every query.
        """
        fg = foreground
        if direction_filter is not None:
            fg = [
                r for r in foreground
                if r.outcome.direction == direction_filter or r.is_negative_result()
            ]
        if not fg:
            return []

        fg_ids = {r.record_id for r in fg}
        bg = [r for r in background if r.record_id not in fg_ids]

        fg_by_pathway = self._index_by_pathway(fg)
        bg_by_pathway = self._index_by_pathway(bg)

        n_fg, n_bg = len(fg), len(bg)
        results: list[PathwayAssociation] = []

        for pkey, recs in fg_by_pathway.items():
            if len(recs) < self.min_records:
                continue
            assoc = self._build_association(
                pkey, recs, len(bg_by_pathway.get(pkey, [])), n_fg, n_bg
            )
            if assoc:
                results.append(assoc)

        self._apply_fdr(results)
        results.sort(key=lambda a: (-a.rank_score(), a.q_value))
        return results

    # ---------------- construction ----------------

    @staticmethod
    def _index_by_pathway(
        records: list[ExperimentRecord],
    ) -> dict[str, list[ExperimentRecord]]:
        idx: dict[str, list[ExperimentRecord]] = defaultdict(list)
        for r in records:
            seen = set()
            for p in r.pathways:
                k = p.key()
                if k not in seen:
                    seen.add(k)
                    idx[k].append(r)
        return idx

    def _build_association(
        self,
        pkey: str,
        recs: list[ExperimentRecord],
        bg_hits: int,
        n_fg: int,
        n_bg: int,
    ) -> PathwayAssociation | None:
        label = next(
            (p.ontology_label or p.raw for r in recs for p in r.pathways
             if p.key() == pkey),
            pkey,
        )
        a = PathwayAssociation(pathway_key=pkey, pathway_label=label)

        positives = [r for r in recs if not r.is_negative_result()]
        negatives = [r for r in recs if r.is_negative_result()]
        support_set = recs if self.include_negative else positives
        if not support_set:
            support_set = recs

        # --- weights and supporting records ---
        total_w = 0.0
        for r in support_set:
            wb = self.weighter.weight(r)
            total_w += wb.composite
            a.supporting.append(
                SupportingRecord(
                    record_id=r.record_id,
                    source_id=r.source.source_id,
                    title=r.source.title,
                    journal=r.source.journal,
                    year=r.source.year,
                    lab_id=r.source.lab_id,
                    direction=r.outcome.direction.value,
                    effect_size=r.outcome.effect_size,
                    p_value=r.outcome.p_value,
                    assay_class=r.assay_class.value,
                    ortho_group=r.orthogonality_group(),
                    weight=round(wb.composite, 4),
                    quote=r.provenance.quote,
                    is_negative=r.is_negative_result(),
                )
            )
        a.weighted_support = total_w

        # --- counts ---
        a.n_records = len(recs)
        a.n_papers = len({r.source.source_id for r in recs})
        labs = [r.source.lab_id for r in recs if r.source.lab_id]
        a.n_labs = len(set(labs))
        a.n_ortho_groups = len({r.orthogonality_group() for r in recs})
        a.n_background_records = bg_hits
        a.n_negative_records = len(negatives)
        a.negative_fraction = len(negatives) / len(recs) if recs else 0.0

        # --- enrichment + Fisher ---
        a_ct = len(recs)
        b_ct = max(0, n_fg - a_ct)
        c_ct = bg_hits
        d_ct = max(0, n_bg - bg_hits)
        a.enrichment = self._odds_ratio(a_ct, b_ct, c_ct, d_ct)
        try:
            _, a.p_fisher = stats.fisher_exact(
                [[a_ct, b_ct], [c_ct, d_ct]], alternative="greater"
            )
        except ValueError:
            a.p_fisher = 1.0

        # --- lab concentration + leave-one-lab-out ---
        a.lab_concentration = self._herfindahl(labs)
        a.leave_one_lab_out_enrichment, a.survives_lolo = self._lolo(
            recs, labs, n_fg, bg_hits, n_bg
        )

        # --- direction consistency ---
        dc = Counter(r.outcome.direction.value for r in recs
                     if r.outcome.direction != Direction.UNSPECIFIED)
        a.direction_counts = dict(dc)
        directional = {k: v for k, v in dc.items() if k != Direction.NO_CHANGE.value}
        if directional:
            total_dir = sum(directional.values())
            top_dir, top_n = max(directional.items(), key=lambda kv: kv[1])
            a.dominant_direction = top_dir
            a.direction_consistency = top_n / total_dir
        else:
            a.dominant_direction = (
                Direction.NO_CHANGE.value if dc else Direction.UNSPECIFIED.value
            )
            a.direction_consistency = 1.0 if dc else 0.0

        # --- effect sizes ---
        es = [r.outcome.effect_size for r in recs if r.outcome.effect_size is not None]
        a.n_with_effect_size = len(es)
        if es:
            es_sorted = sorted(es)
            a.median_effect_size = float(
                stats.scoreatpercentile(es_sorted, 50)
            )
            if len(es) >= 4:
                a.effect_size_iqr = (
                    float(stats.scoreatpercentile(es_sorted, 25)),
                    float(stats.scoreatpercentile(es_sorted, 75)),
                )

        # --- temporal ---
        years = [r.source.year for r in recs if r.source.year]
        if years:
            a.first_year, a.last_year = min(years), max(years)

        a.publication_bias_flag = self._bias_flag(a, recs)
        return a

    # ---------------- statistics helpers ----------------

    @staticmethod
    def _odds_ratio(a: int, b: int, c: int, d: int) -> float:
        """Haldane-Anscombe corrected odds ratio (adds 0.5 to each cell).

        The correction prevents division-by-zero and infinite ORs when a
        pathway appears in the foreground but never in the background, which
        is common for specific pathways in a focused corpus.
        """
        return ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))

    @staticmethod
    def _herfindahl(labs: list[str]) -> float:
        """Concentration of evidence across labs. 1.0 = all from one lab."""
        if not labs:
            return 1.0
        counts = Counter(labs)
        total = sum(counts.values())
        return sum((c / total) ** 2 for c in counts.values())

    def _lolo(
        self,
        recs: list[ExperimentRecord],
        labs: list[str],
        n_fg: int,
        bg_hits: int,
        n_bg: int,
    ) -> tuple[float | None, bool]:
        """Leave-one-lab-out: drop the most prolific lab, recompute enrichment.

        This is the single most useful robustness check for this data. A real
        biological association should survive it; a one-group narrative
        typically does not.
        """
        if not labs or len(set(labs)) < 2:
            return (None, False)
        top_lab = Counter(labs).most_common(1)[0][0]
        kept = [r for r in recs if r.source.lab_id != top_lab]
        if not kept:
            return (None, False)
        a_ct = len(kept)
        b_ct = max(0, n_fg - a_ct)
        or_lolo = self._odds_ratio(a_ct, b_ct, bg_hits, max(0, n_bg - bg_hits))
        try:
            _, p = stats.fisher_exact(
                [[a_ct, b_ct], [bg_hits, max(0, n_bg - bg_hits)]],
                alternative="greater",
            )
        except ValueError:
            p = 1.0
        return (or_lolo, bool(or_lolo > 1.2 and p < 0.05))

    @staticmethod
    def _bias_flag(a: PathwayAssociation, recs: list[ExperimentRecord]) -> str:
        """Heuristic publication-bias warnings."""
        if a.n_records >= 5 and a.n_negative_records == 0:
            return (
                f"No null results among {a.n_records} reports -- literature may be "
                "filtered for positive findings."
            )
        small_n = [
            r for r in recs
            if r.outcome.n_biological_replicates is not None
            and r.outcome.n_biological_replicates <= 3
        ]
        sig = [r for r in recs if r.outcome.p_value is not None and r.outcome.p_value < 0.05]
        if len(small_n) >= 3 and len(sig) >= 3 and len(sig) == len(
            [r for r in recs if r.outcome.p_value is not None]
        ):
            return (
                "All reported p-values are significant despite small samples -- "
                "possible selective reporting."
            )
        return ""

    @staticmethod
    def _apply_fdr(results: list[PathwayAssociation]) -> None:
        """Benjamini-Hochberg across all tested pathways.

        Necessary because we test hundreds of pathways per query; uncorrected
        p-values would produce a wall of spurious hits.
        """
        if not results:
            return
        m = len(results)
        indexed = sorted(enumerate(results), key=lambda t: t[1].p_fisher)
        prev_q = 1.0
        for rank in range(m - 1, -1, -1):
            i, assoc = indexed[rank]
            q = min(prev_q, assoc.p_fisher * m / (rank + 1))
            assoc.q_value = min(1.0, q)
            prev_q = q


# --------------------------------------------------------------------------
# Cross-phenotype coherence
# --------------------------------------------------------------------------


def phenotype_coherence(
    associations_a: list[PathwayAssociation],
    associations_b: list[PathwayAssociation],
    top_n: int = 30,
) -> dict:
    """Compare pathway profiles of two phenotypes.

    Useful concretely: if your imaging phenotype and your electrophysiology
    phenotype converge on the same pathways, that is genuine orthogonal
    support. If they do not, one of them may be an artifact.
    """
    a_top = {x.pathway_key: x for x in associations_a[:top_n]}
    b_top = {x.pathway_key: x for x in associations_b[:top_n]}
    shared = set(a_top) & set(b_top)
    union = set(a_top) | set(b_top)
    return {
        "jaccard": len(shared) / len(union) if union else 0.0,
        "shared_pathways": [
            {
                "pathway": a_top[k].pathway_label,
                "enrichment_a": round(a_top[k].enrichment, 2),
                "enrichment_b": round(b_top[k].enrichment, 2),
                "direction_a": a_top[k].dominant_direction,
                "direction_b": b_top[k].dominant_direction,
                "consistent_direction": a_top[k].dominant_direction
                == b_top[k].dominant_direction,
            }
            for k in sorted(shared, key=lambda k: -a_top[k].weighted_support)
        ],
        "only_a": [a_top[k].pathway_label for k in set(a_top) - shared][:10],
        "only_b": [b_top[k].pathway_label for k in set(b_top) - shared][:10],
    }
