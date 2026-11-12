"""
Sequence <-> literature reconciliation for splicing.

This is the module that does something no existing tool does for iPSC-neuron
work: it takes what a sequence model predicts about a variant and asks whether
the published literature actually observed that.

Four outcomes, each meaning something different and each requiring a different
response from the researcher:

  CORROBORATED   Model predicts an effect; literature reports the same effect
                 in the same gene, same direction, same event class. The
                 mechanism is sequence-explicable AND observed. This is the
                 strongest evidential state available short of doing the
                 experiment yourself.

  PREDICTED_ONLY Model predicts an effect; nobody has looked. This is the
                 highest-value experimental opportunity in the system --
                 a testable, sequence-grounded hypothesis with no prior art.
                 These are surfaced prominently and deliberately.

  OBSERVED_ONLY  Literature reports an effect; sequence does not predict it.
                 Either the mechanism is indirect (the variant acts through a
                 regulator, not a splice site), or the literature report is
                 wrong, or the model missed it. All three are informative.

  CONTRADICTED   Model predicts one effect; literature reports a different or
                 opposite one. Something is wrong. Finding out what is usually
                 more valuable than either claim taken alone.

Why this matters
----------------
Splice prediction tools tell you what SHOULD happen. The literature tells you
what someone SAW happen. These are different epistemic objects and the gap
between them is where the biology hides. A variant that SpliceAI calls at 0.9
but that nobody has ever reported as splice-altering is not a boring negative
-- it is either an untested prediction or a systematic blind spot in how the
field assays that gene.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from enum import Enum

from ..evidence.reconcile import (
    EvidenceClass, EvidenceItem, EvidenceReconciler, ReconciliationResult,
)
from ..infer.splicing import NormalizedSplicingEvent
from ..schema import Direction, SplicingEventType
from .splice import (
    SPLICEAI_THRESHOLDS, EnsemblePrediction, SpliceEffect, SplicePrediction,
)

log = logging.getLogger(__name__)


class ReconciliationState(str, Enum):
    CORROBORATED = "corroborated"
    PREDICTED_ONLY = "predicted_only"
    OBSERVED_ONLY = "observed_only"
    CONTRADICTED = "contradicted"
    BOTH_SILENT = "both_silent"


# Which predicted splice effect should produce which observed direction of
# inclusion. A lost donor or acceptor removes a canonical exon -> inclusion
# goes DOWN. A gained site creates a new exon or shifts a boundary ->
# inclusion of the new form goes UP.
EFFECT_TO_EXPECTED_DIRECTION = {
    SpliceEffect.DONOR_LOSS: Direction.DECREASE,
    SpliceEffect.ACCEPTOR_LOSS: Direction.DECREASE,
    SpliceEffect.DONOR_GAIN: Direction.INCREASE,
    SpliceEffect.ACCEPTOR_GAIN: Direction.INCREASE,
}


@dataclass
class SpliceReconciliation:
    """One variant, reconciled between prediction and literature."""

    variant_id: str
    gene_symbol: str | None = None
    state: ReconciliationState = ReconciliationState.BOTH_SILENT

    predicted_effect: SpliceEffect = SpliceEffect.NONE
    predicted_delta: float = 0.0
    predicted_event_type: str = "unspecified"
    predicted_direction: Direction = Direction.UNSPECIFIED
    prediction_available: bool = False

    observed_events: list[NormalizedSplicingEvent] = field(default_factory=list)
    n_observing_papers: int = 0
    n_observing_labs: int = 0
    observed_direction: Direction = Direction.UNSPECIFIED
    observed_event_types: list[str] = field(default_factory=list)

    event_type_match: bool = False
    direction_match: bool = False

    reconciliation: ReconciliationResult | None = None
    narrative: str = ""
    opportunity: str = ""
    caveats: list[str] = field(default_factory=list)
    priority: float = 0.0

    def is_novel_hypothesis(self) -> bool:
        """A confident prediction with no published observation."""
        return (
            self.state == ReconciliationState.PREDICTED_ONLY
            and self.predicted_delta >= SPLICEAI_THRESHOLDS["recommended"]
        )


class SpliceLiteratureReconciler:
    """Reconciles splice predictions against observed splicing events."""

    def __init__(
        self,
        significance_threshold: float = SPLICEAI_THRESHOLDS["recommended"],
        weak_threshold: float = SPLICEAI_THRESHOLDS["high_recall"],
    ):
        self.threshold = significance_threshold
        self.weak = weak_threshold
        self.reconciler = EvidenceReconciler()

    def reconcile(
        self,
        prediction: SplicePrediction | EnsemblePrediction | None,
        observed: list[NormalizedSplicingEvent],
        gene_symbol: str | None = None,
        variant_id: str = "",
        lab_ids: dict[str, str] | None = None,
    ) -> SpliceReconciliation:
        rec = SpliceReconciliation(
            variant_id=variant_id or getattr(prediction, "variant_id", ""),
            gene_symbol=gene_symbol or getattr(prediction, "gene_symbol", None),
        )

        self._load_prediction(rec, prediction)
        self._load_observations(rec, observed, lab_ids or {})
        rec.state = self._classify(rec)
        rec.reconciliation = self._build_evidence(rec)
        rec.narrative = self._narrate(rec)
        rec.opportunity = self._opportunity(rec)
        rec.priority = self._priority(rec)
        rec.caveats = self._caveats(rec, prediction)
        return rec

    # ---------------- loading ----------------

    def _load_prediction(self, rec: SpliceReconciliation, prediction) -> None:
        if prediction is None:
            return
        err = getattr(prediction, "error", None)
        if err:
            return

        if isinstance(prediction, EnsemblePrediction):
            usable = [m for m in prediction.members if not m.error]
            if not usable:
                return
            rec.prediction_available = True
            rec.predicted_delta = prediction.max_delta
            rec.predicted_effect = prediction.consensus_effect
            best = max(usable, key=lambda m: m.max_delta)
            rec.predicted_event_type = best.predicted_event_type()
        else:
            rec.prediction_available = True
            rec.predicted_delta = prediction.max_delta
            rec.predicted_effect = prediction.dominant_effect
            rec.predicted_event_type = prediction.predicted_event_type()

        rec.predicted_direction = EFFECT_TO_EXPECTED_DIRECTION.get(
            rec.predicted_effect, Direction.UNSPECIFIED
        )

    def _load_observations(
        self, rec: SpliceReconciliation,
        observed: list[NormalizedSplicingEvent],
        lab_ids: dict[str, str],
    ) -> None:
        if rec.gene_symbol:
            target = rec.gene_symbol.upper()
            observed = [e for e in observed if e.gene.upper() == target]
        rec.observed_events = observed
        if not observed:
            return

        rec.n_observing_papers = len({e.source_id for e in observed})
        rec.n_observing_labs = len(
            {lab_ids.get(e.source_id) for e in observed
             if lab_ids.get(e.source_id)}
        ) or rec.n_observing_papers

        rec.observed_event_types = sorted(
            {e.event_type.value for e in observed
             if e.event_type != SplicingEventType.UNSPECIFIED}
        )

        dirs = [e.direction for e in observed
                if e.direction in (Direction.INCREASE, Direction.DECREASE)]
        if dirs:
            rec.observed_direction = max(set(dirs), key=dirs.count)

        rec.event_type_match = (
            rec.predicted_event_type in rec.observed_event_types
            if rec.observed_event_types and
            rec.predicted_event_type != "unspecified" else False
        )
        rec.direction_match = (
            rec.predicted_direction == rec.observed_direction
            and rec.predicted_direction != Direction.UNSPECIFIED
        )

    # ---------------- classification ----------------

    def _classify(self, rec: SpliceReconciliation) -> ReconciliationState:
        predicts = rec.prediction_available and rec.predicted_delta >= self.threshold
        predicts_weakly = (
            rec.prediction_available and rec.predicted_delta >= self.weak
        )
        observes = bool(rec.observed_events)

        if not predicts_weakly and not observes:
            return ReconciliationState.BOTH_SILENT
        if predicts and not observes:
            return ReconciliationState.PREDICTED_ONLY
        if observes and not predicts_weakly:
            return ReconciliationState.OBSERVED_ONLY
        if predicts_weakly and observes:
            # Both say something. Do they say the same thing?
            if rec.direction_match or rec.event_type_match:
                return ReconciliationState.CORROBORATED
            if (
                rec.predicted_direction != Direction.UNSPECIFIED
                and rec.observed_direction != Direction.UNSPECIFIED
                and rec.predicted_direction != rec.observed_direction
            ):
                return ReconciliationState.CONTRADICTED
            # Same gene, both report something, but the specifics do not line
            # up -- treat as corroborated-at-gene-level with a caveat rather
            # than claiming contradiction we cannot substantiate.
            return ReconciliationState.CORROBORATED
        return ReconciliationState.BOTH_SILENT

    # ---------------- evidence assembly ----------------

    def _build_evidence(self, rec: SpliceReconciliation) -> ReconciliationResult:
        items: list[EvidenceItem] = []

        # When the two sources disagree about the SPECIFICS, they are not both
        # supporting one claim -- they are supporting incompatible ones. The
        # claim we reconcile against is the model's, so a literature report
        # that contradicts the model's direction counts as refuting it.
        # Without this, a contradiction is scored as convergence, which is
        # exactly backwards.
        contradicted = rec.state == ReconciliationState.CONTRADICTED

        if rec.prediction_available:
            supports = rec.predicted_delta >= self.threshold
            items.append(EvidenceItem(
                source_class=EvidenceClass.SEQUENCE,
                source_name="splice prediction",
                claim_direction="supports" if supports else "refutes",
                magnitude=rec.predicted_delta,
                magnitude_unit="delta_score",
                confidence=min(0.85, 0.35 + 0.5 * rec.predicted_delta),
                provenance=rec.variant_id,
            ))

        if rec.observed_events:
            quantified = [e for e in rec.observed_events
                          if e.delta_psi is not None]
            median_dpsi = None
            if quantified:
                vals = sorted(abs(e.delta_psi) for e in quantified)
                median_dpsi = vals[len(vals) // 2]
            items.append(EvidenceItem(
                source_class=EvidenceClass.LITERATURE,
                source_name="published splicing analyses",
                claim_direction="refutes" if contradicted else "supports",
                magnitude=median_dpsi,
                magnitude_unit="abs_delta_psi",
                confidence=0.5 + 0.1 * min(4, rec.n_observing_labs),
                n_independent=rec.n_observing_labs,
                provenance=f"{rec.n_observing_papers} papers",
                caveats=(
                    [f"Reports {rec.observed_direction.value} where the model "
                     f"predicts {rec.predicted_direction.value}."]
                    if contradicted else []
                ),
            ))

        claim = (
            f"Variant {rec.variant_id} alters splicing of "
            f"{rec.gene_symbol or 'the target gene'} as the sequence model "
            f"predicts"
        )
        return self.reconciler.reconcile(claim, items)

    # ---------------- narrative ----------------

    def _narrate(self, rec: SpliceReconciliation) -> str:
        g = rec.gene_symbol or "this gene"
        if rec.state == ReconciliationState.CORROBORATED:
            detail = []
            if rec.event_type_match:
                detail.append(f"event class ({rec.predicted_event_type})")
            if rec.direction_match:
                detail.append(f"direction ({rec.observed_direction.value})")
            agree = " and ".join(detail) if detail else "the gene affected"
            return (
                f"Sequence prediction and published observation agree on "
                f"{agree}. The model predicts "
                f"{rec.predicted_effect.value.replace('_', ' ')} "
                f"(delta {rec.predicted_delta:.2f}) and "
                f"{rec.n_observing_papers} paper(s) from "
                f"{rec.n_observing_labs} group(s) report splicing changes in "
                f"{g}. These are independent lines of evidence: one from "
                f"primary sequence, one from observation."
            )
        if rec.state == ReconciliationState.PREDICTED_ONLY:
            return (
                f"The sequence model predicts "
                f"{rec.predicted_effect.value.replace('_', ' ')} in {g} "
                f"(delta {rec.predicted_delta:.2f}), but no paper in the corpus "
                f"reports a splicing change for this gene. Nobody appears to "
                f"have looked."
            )
        if rec.state == ReconciliationState.OBSERVED_ONLY:
            return (
                f"{rec.n_observing_papers} paper(s) report splicing changes in "
                f"{g}, but the sequence model predicts no splice-site effect "
                f"for this variant (delta {rec.predicted_delta:.2f}). The most "
                f"likely explanations are that the mechanism is indirect -- "
                f"acting through a splicing regulator rather than a splice "
                f"site -- or that the observed change is driven by a different "
                f"variant."
            )
        if rec.state == ReconciliationState.CONTRADICTED:
            return (
                f"Direct conflict. The model predicts inclusion should "
                f"{rec.predicted_direction.value} while the literature reports "
                f"it {rec.observed_direction.value} in {g}. One of these is "
                f"wrong, and determining which would be genuinely informative."
            )
        return (
            f"Neither the sequence model nor the literature indicates a "
            f"splicing effect for this variant in {g}."
        )

    def _opportunity(self, rec: SpliceReconciliation) -> str:
        if rec.state == ReconciliationState.PREDICTED_ONLY:
            return (
                "HIGH-VALUE UNTESTED PREDICTION. A minigene assay or targeted "
                "RT-PCR in your line would test a specific, sequence-grounded "
                "hypothesis that nobody has published. If it holds, it is novel; "
                "if it does not, that is a calibration datapoint on the model."
            )
        if rec.state == ReconciliationState.CONTRADICTED:
            return (
                "Resolving this contradiction is worth more than either claim. "
                "Check first whether the published work used a different "
                "transcript isoform or genome build -- coordinate mismatches "
                "explain many apparent contradictions."
            )
        if rec.state == ReconciliationState.OBSERVED_ONLY:
            return (
                "Look for an indirect mechanism: is a splicing regulator "
                "(TDP-43, PTBP1, SRSF1) implicated in the same papers? If so, "
                "the variant may act on the regulator rather than on this "
                "gene's splice sites directly."
            )
        if rec.state == ReconciliationState.CORROBORATED:
            return (
                "Well-supported from two independent directions. Reasonable to "
                "build on, though it has still not been tested in your specific "
                "line and maturity state."
            )
        return ""

    @staticmethod
    def _priority(rec: SpliceReconciliation) -> float:
        """Experimental priority. Untested confident predictions rank highest.

        This ordering is a deliberate stance: the tool is most useful when it
        points at things nobody has done, not when it confirms what everyone
        already knows.
        """
        if rec.state == ReconciliationState.PREDICTED_ONLY:
            return 3.0 * rec.predicted_delta
        if rec.state == ReconciliationState.CONTRADICTED:
            return 2.5 * max(rec.predicted_delta, 0.5)
        if rec.state == ReconciliationState.OBSERVED_ONLY:
            return 1.5 * min(1.0, rec.n_observing_labs / 3.0)
        if rec.state == ReconciliationState.CORROBORATED:
            return 1.0 * rec.predicted_delta
        return 0.0

    @staticmethod
    def _caveats(rec: SpliceReconciliation, prediction) -> list[str]:
        out: list[str] = []
        if prediction is not None:
            out.extend(getattr(prediction, "caveats", [])[:4])
        if rec.state == ReconciliationState.PREDICTED_ONLY:
            out.append(
                "Absence from the corpus is not absence from the literature: "
                "the paper may exist but be paywalled, or may not have been "
                "harvested. Check manually before claiming novelty."
            )
        if rec.observed_events:
            uncertain = [e for e in rec.observed_events
                         if e.sign_confidence != "high"]
            if uncertain:
                out.append(
                    f"{len(uncertain)} observed event(s) have uncertain dPSI "
                    "sign, so the direction comparison may be unreliable."
                )
        if rec.gene_symbol is None:
            out.append(
                "No gene symbol resolved, so literature matching was done "
                "without gene filtering and may include unrelated events."
            )
        return out


# --------------------------------------------------------------------------
# Batch analysis
# --------------------------------------------------------------------------


@dataclass
class VariantPanel:
    """Reconciliation across a set of variants, ranked by opportunity."""

    reconciliations: list[SpliceReconciliation] = field(default_factory=list)

    def novel_hypotheses(self) -> list[SpliceReconciliation]:
        return sorted(
            [r for r in self.reconciliations if r.is_novel_hypothesis()],
            key=lambda r: -r.priority,
        )

    def contradictions(self) -> list[SpliceReconciliation]:
        return [r for r in self.reconciliations
                if r.state == ReconciliationState.CONTRADICTED]

    def corroborated(self) -> list[SpliceReconciliation]:
        return [r for r in self.reconciliations
                if r.state == ReconciliationState.CORROBORATED]

    def by_priority(self) -> list[SpliceReconciliation]:
        return sorted(self.reconciliations, key=lambda r: -r.priority)

    def summary(self) -> dict:
        from collections import Counter
        states = Counter(r.state.value for r in self.reconciliations)
        return {
            "n_variants": len(self.reconciliations),
            "states": dict(states),
            "n_novel_hypotheses": len(self.novel_hypotheses()),
            "n_contradictions": len(self.contradictions()),
            "n_corroborated": len(self.corroborated()),
        }
