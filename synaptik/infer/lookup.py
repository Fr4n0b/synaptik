"""
Component 1: "Has this experiment already been done?"

The naive version of this is keyword search, which answers a different and
less useful question ("what papers mention these words"). What the user needs
is a graded answer along the dimensions that determine whether a prior result
actually applies to their planned experiment:

  EXACT       -- same system, same perturbation, same assay, same maturity.
                 The experiment has been done. Read these papers first.
  CLOSE       -- differs in one non-critical dimension (e.g. DIV 21 vs 25,
                 or NGN2 vs small-molecule induction). Likely informative.
  RELATED     -- same perturbation, different readout; or same readout,
                 different system. Useful context, not a direct answer.
  ADJACENT    -- same gene or pathway, otherwise different. Background.

Certain mismatches are disqualifying regardless of other similarity, because
the results are not comparable:
  - different species (human iPSC vs mouse primary)
  - different maturity bin (DIV14 vs DIV49 synapse counts)
  - opposite perturbation direction (knockout vs overexpression)

These are encoded as hard blocks rather than soft penalties, because a
similarity score that lets a mouse primary-neuron study rank as "exact match"
for a human iPSC experiment is worse than no answer at all.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from enum import Enum

from ..schema import (
    AssayClass, CellModelClass, ExperimentQuery, ExperimentRecord,
    PerturbationType,
)
from .weighting import EvidenceWeighter

# Perturbation types that push in opposite directions. Comparing a knockout
# result to an overexpression result as if they were the same experiment is a
# category error, so we block the match.
OPPOSING_PERTURBATIONS = {
    frozenset({PerturbationType.KNOCKOUT, PerturbationType.OVEREXPRESSION}),
    frozenset({PerturbationType.KNOCKDOWN, PerturbationType.OVEREXPRESSION}),
}

# Perturbations grouped by the direction they push gene dosage/function.
# Enumerating pairs is not enough: a patient loss-of-function mutation versus
# an overexpression construct is just as opposed as knockout versus
# overexpression, and pairwise listing misses it. Classifying each type and
# blocking on class mismatch catches the whole family.
PERTURBATION_VALENCE: dict[PerturbationType, str] = {
    PerturbationType.KNOCKOUT: "loss",
    PerturbationType.KNOCKDOWN: "loss",
    PerturbationType.MUTATION_PATIENT: "variant",   # usually loss, not always
    PerturbationType.MUTATION_ENGINEERED: "variant",
    PerturbationType.OVEREXPRESSION: "gain",
    PerturbationType.ASO: "loss",
    PerturbationType.ISOFORM_SWITCH: "shift",
    PerturbationType.SMALL_MOLECULE: "modulator",
    PerturbationType.NONE: "none",
    PerturbationType.OTHER: "unknown",
}

# Valence pairs that must never be treated as comparable experiments.
BLOCKED_VALENCE_PAIRS = {
    frozenset({"loss", "gain"}),
    frozenset({"variant", "gain"}),
}


def opposing(a: PerturbationType, b: PerturbationType) -> bool:
    """True if two perturbations push in incompatible directions."""
    if frozenset({a, b}) in OPPOSING_PERTURBATIONS:
        return True
    va = PERTURBATION_VALENCE.get(a, "unknown")
    vb = PERTURBATION_VALENCE.get(b, "unknown")
    return frozenset({va, vb}) in BLOCKED_VALENCE_PAIRS

# Cell models close enough that a result transfers with a caveat.
MODEL_NEIGHBORHOODS = {
    CellModelClass.IPSC_NEURON_NGN2: {
        CellModelClass.IPSC_NEURON_SMALL_MOLECULE: 0.75,
        CellModelClass.IPSC_ORGANOID: 0.45,
        CellModelClass.IPSC_NPC: 0.30,
        CellModelClass.PRIMARY_RODENT_NEURON: 0.25,
    },
    CellModelClass.IPSC_NEURON_SMALL_MOLECULE: {
        CellModelClass.IPSC_NEURON_NGN2: 0.75,
        CellModelClass.IPSC_ORGANOID: 0.45,
    },
}

MATURITY_ADJACENCY = {
    ("early", "mid"): 0.5, ("mid", "late"): 0.6, ("late", "extended"): 0.7,
}


class MatchTier(str, Enum):
    EXACT = "exact"
    CLOSE = "close"
    RELATED = "related"
    ADJACENT = "adjacent"


@dataclass
class ProtocolSummary:
    """What a prior study actually did -- the practical payload of a match."""

    differentiation_protocol: str | None = None
    days_in_vitro: int | None = None
    coculture: str | None = None
    markers: list[str] = field(default_factory=list)
    n_biological_replicates: int | None = None
    statistical_test: str | None = None
    is_isogenic_pair: bool | None = None
    control_description: str | None = None
    dose: str | None = None
    duration: str | None = None


@dataclass
class ExperimentMatch:
    record: ExperimentRecord
    tier: MatchTier
    similarity: float
    weight: float
    matched_on: list[str] = field(default_factory=list)
    differs_on: list[str] = field(default_factory=list)

    def protocol(self) -> ProtocolSummary:
        r = self.record
        return ProtocolSummary(
            differentiation_protocol=r.cell_model.differentiation_protocol,
            days_in_vitro=r.cell_model.days_in_vitro,
            coculture=r.cell_model.coculture,
            markers=r.outcome.markers_used,
            n_biological_replicates=r.outcome.n_biological_replicates,
            statistical_test=r.outcome.statistical_test,
            is_isogenic_pair=r.cell_model.is_isogenic_pair,
            control_description=r.perturbation.control_description,
            dose=r.perturbation.dose,
            duration=r.perturbation.duration,
        )

    def summary_line(self) -> str:
        r = self.record
        auth = r.source.first_author_surname() or "?"
        yr = r.source.year or "n.d."
        es = (
            f"{r.outcome.effect_size:+.0f}{'%' if r.outcome.effect_size_unit == 'percent_change' else ''}"
            if r.outcome.effect_size is not None else r.outcome.direction.value
        )
        return (
            f"{auth} et al. {yr} ({r.source.journal or 'unknown'}): "
            f"{r.outcome.phenotype.raw} {es}"
        )


@dataclass
class LookupResult:
    query: ExperimentQuery
    matches: list[ExperimentMatch]
    n_exact: int = 0
    n_close: int = 0
    verdict: str = ""
    consensus: str = ""
    conflicting: bool = False

    def by_tier(self, tier: MatchTier) -> list[ExperimentMatch]:
        return [m for m in self.matches if m.tier == tier]


class ExperimentLookup:
    def __init__(self, weighter: EvidenceWeighter | None = None):
        self.weighter = weighter or EvidenceWeighter()

    def lookup(
        self,
        query: ExperimentQuery,
        corpus: list[ExperimentRecord],
        max_results: int = 60,
    ) -> LookupResult:
        matches: list[ExperimentMatch] = []
        for rec in corpus:
            m = self._score(query, rec)
            if m:
                matches.append(m)

        matches.sort(key=lambda m: (-m.similarity, -m.weight))
        matches = matches[:max_results]

        res = LookupResult(query=query, matches=matches)
        res.n_exact = sum(1 for m in matches if m.tier == MatchTier.EXACT)
        res.n_close = sum(1 for m in matches if m.tier == MatchTier.CLOSE)
        self._verdict(res)
        return res

    # ---------------- scoring ----------------

    def _score(
        self, q: ExperimentQuery, r: ExperimentRecord
    ) -> ExperimentMatch | None:
        matched: list[str] = []
        differs: list[str] = []
        score = 0.0
        max_score = 0.0
        blocked = False

        # --- cell model (weight 3) ---
        if q.cell_model_class:
            max_score += 3.0
            if r.cell_model.model_class == q.cell_model_class:
                score += 3.0
                matched.append(f"cell model ({q.cell_model_class.value})")
            else:
                neigh = MODEL_NEIGHBORHOODS.get(q.cell_model_class, {})
                sim = neigh.get(r.cell_model.model_class, 0.0)
                if sim == 0.0:
                    blocked = True
                    differs.append(
                        f"different system ({r.cell_model.model_class.value})"
                    )
                else:
                    score += 3.0 * sim
                    differs.append(
                        f"system ({r.cell_model.model_class.value} vs "
                        f"{q.cell_model_class.value})"
                    )

        # --- maturity (weight 2) -- hard-blocking when far apart ---
        if q.days_in_vitro is not None:
            max_score += 2.0
            qb, rb = q.maturity_bin(), r.cell_model.maturity_bin()
            if rb == "unknown":
                differs.append("DIV not reported")
            elif qb == rb:
                score += 2.0
                matched.append(f"maturity ({rb}, DIV~{r.cell_model.days_in_vitro})")
            else:
                adj = MATURITY_ADJACENCY.get(
                    (qb, rb), MATURITY_ADJACENCY.get((rb, qb), 0.0)
                )
                if adj == 0.0:
                    blocked = True
                    differs.append(f"incomparable maturity ({rb} vs {qb})")
                else:
                    score += 2.0 * adj
                    differs.append(f"maturity ({rb} vs {qb})")

        # --- perturbation type (weight 2) ---
        if q.perturbation_type:
            max_score += 2.0
            if r.perturbation.ptype == q.perturbation_type:
                score += 2.0
                matched.append(f"perturbation ({q.perturbation_type.value})")
            elif opposing(q.perturbation_type, r.perturbation.ptype):
                blocked = True
                differs.append(
                    f"OPPOSING perturbation ({r.perturbation.ptype.value} vs "
                    f"{q.perturbation_type.value}) -- effect sign not comparable"
                )
            else:
                differs.append(f"perturbation ({r.perturbation.ptype.value})")

        # --- perturbation target (weight 3) ---
        if q.perturbation_target:
            max_score += 3.0
            tgt = (r.perturbation.target.raw if r.perturbation.target else "")
            gene_hits = {g.raw.upper() for g in r.genes}
            qt = q.perturbation_target.upper()
            if tgt and tgt.upper() == qt:
                score += 3.0
                matched.append(f"target gene ({tgt})")
            elif qt in gene_hits:
                score += 1.5
                matched.append(f"gene mentioned ({q.perturbation_target})")
            else:
                differs.append(f"target ({tgt or 'none'})")

        # --- exact variant (weight 4) -- the strongest possible signal ---
        if q.variant_hgvs:
            max_score += 4.0
            if r.perturbation.variant_hgvs == q.variant_hgvs:
                score += 4.0
                matched.append(f"exact variant ({q.variant_hgvs})")

        # --- assay (weight 3) ---
        if q.assay_class:
            max_score += 3.0
            if r.assay_class == q.assay_class:
                score += 3.0
                matched.append(f"assay ({q.assay_class.value})")
            elif r.orthogonality_group() == _ortho(q.assay_class):
                score += 1.5
                matched.append(f"related assay ({r.assay_class.value})")
            else:
                differs.append(f"assay ({r.assay_class.value})")

        # --- markers (weight 2) ---
        if q.markers:
            max_score += 2.0
            qm = {m.upper() for m in q.markers}
            rm = {m.upper() for m in r.outcome.markers_used}
            if qm and rm:
                overlap = len(qm & rm) / len(qm)
                score += 2.0 * overlap
                if overlap >= 0.99:
                    matched.append(f"markers ({', '.join(sorted(qm & rm))})")
                elif overlap > 0:
                    matched.append(f"some markers ({', '.join(sorted(qm & rm))})")
                else:
                    differs.append(f"markers ({', '.join(sorted(rm))})")

        # --- phenotype (weight 3) ---
        if q.phenotype:
            max_score += 3.0
            qp = q.phenotype.lower()
            rp = r.outcome.phenotype.raw.lower()
            rid = r.outcome.phenotype.ontology_id
            if qp == rp or (rid and qp == rid.lower()):
                score += 3.0
                matched.append(f"phenotype ({r.outcome.phenotype.raw})")
            elif qp in rp or rp in qp:
                score += 2.0
                matched.append(f"similar phenotype ({r.outcome.phenotype.raw})")
            else:
                differs.append(f"phenotype ({r.outcome.phenotype.raw})")

        if max_score == 0.0:
            return None
        similarity = score / max_score

        # A match cannot be "exact" unless the READOUT agrees. Without this,
        # a query that specifies system + DIV + target but no phenotype scores
        # 1.0 against every record from a matching system -- so a viability
        # assay is reported as an exact match for a splicing experiment. The
        # readout is what makes two experiments the same experiment.
        readout_specified = bool(q.phenotype or q.assay_class)
        readout_agrees = False
        if q.phenotype:
            qp = q.phenotype.lower()
            rp = r.outcome.phenotype.raw.lower()
            rid = (r.outcome.phenotype.ontology_id or "").lower()
            readout_agrees = qp == rp or qp == rid or qp in rp or rp in qp
        elif q.assay_class:
            # Assay equality alone is a weak proxy for "same readout": the same
            # instrument measures many different things. We additionally require
            # that the record's phenotype is one this assay actually produces,
            # which filters records whose assay and phenotype were mis-paired
            # during extraction.
            readout_agrees = (
                r.assay_class == q.assay_class
                and _phenotype_fits_assay(r.outcome.phenotype.ontology_id,
                                          q.assay_class)
            )

        if blocked:
            # Blocked matches are demoted, never promoted above 'related'.
            similarity = min(similarity, 0.45)

        if similarity < 0.20:
            return None

        if blocked:
            tier = MatchTier.RELATED if similarity >= 0.35 else MatchTier.ADJACENT
        elif similarity >= 0.85 and (readout_agrees or not readout_specified):
            tier = MatchTier.EXACT
        elif similarity >= 0.65:
            tier = MatchTier.CLOSE
        elif similarity >= 0.40:
            tier = MatchTier.RELATED
        else:
            tier = MatchTier.ADJACENT

        # Cap at CLOSE when the user named a readout and this record used a
        # different one, however well everything else lines up.
        if readout_specified and not readout_agrees and tier == MatchTier.EXACT:
            tier = MatchTier.CLOSE
            differs.append("different readout")

        return ExperimentMatch(
            record=r,
            tier=tier,
            similarity=round(similarity, 3),
            weight=round(self.weighter.weight(r).composite, 3),
            matched_on=matched,
            differs_on=differs,
        )

    # ---------------- verdict ----------------

    @staticmethod
    def _verdict(res: LookupResult) -> None:
        strong = [m for m in res.matches
                  if m.tier in (MatchTier.EXACT, MatchTier.CLOSE)]
        if not strong:
            related = res.by_tier(MatchTier.RELATED)
            res.verdict = (
                "No prior study matches this experiment closely. "
                + (f"{len(related)} related studies provide partial context."
                   if related else "The literature appears not to cover this.")
            )
            return

        labs = {m.record.source.lab_id for m in strong if m.record.source.lab_id}
        dirs = defaultdict(float)
        for m in strong:
            dirs[m.record.outcome.direction.value] += m.weight

        if dirs:
            top_dir, top_w = max(dirs.items(), key=lambda kv: kv[1])
            total_w = sum(dirs.values()) or 1.0
            frac = top_w / total_w
            res.conflicting = frac < 0.7 and len(dirs) > 1
            res.consensus = (
                f"{top_dir} ({frac:.0%} of weighted evidence)"
                if not res.conflicting
                else f"CONFLICTING: {dict((k, round(v, 2)) for k, v in dirs.items())}"
            )

        n_lab = len(labs)
        res.verdict = (
            f"{res.n_exact} exact and {res.n_close} close prior experiments "
            f"from {n_lab} independent lab{'s' if n_lab != 1 else ''}. "
            + (
                "Direction of effect is disputed -- read the conflicting reports "
                "before designing your version."
                if res.conflicting
                else f"Consensus: {res.consensus}."
            )
        )
        if n_lab == 1 and (res.n_exact + res.n_close) > 1:
            res.verdict += (
                " Note: all matches come from a single group, so this is "
                "not independent replication."
            )


def _ortho(assay: AssayClass) -> str:
    from ..schema import ORTHOGONALITY_GROUPS
    return ORTHOGONALITY_GROUPS.get(assay, "other")


def _phenotype_fits_assay(
    phenotype_id: str | None, assay: AssayClass
) -> bool:
    """Is this phenotype something the named assay can actually measure?

    Catches extraction errors where an assay and an outcome from different
    experiments in the same paper were paired into one record. A viability
    number attributed to a splicing assay is a mis-parse, and treating it as
    an exact match for a splicing experiment would be actively misleading.

    Unknown phenotypes pass, because we would rather admit an uncertain match
    than hide a real one.
    """
    if not phenotype_id:
        return True
    from ..normalize.phenotype_ontology import BY_ID

    term = BY_ID.get(phenotype_id)
    if not term or not term.typical_assays:
        return True
    return assay.value in term.typical_assays
