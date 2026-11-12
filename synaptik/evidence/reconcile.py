"""
Multi-modal evidence integration.

This is the conceptual heart of SYNAPTIK v2. Everything else feeds it.

The problem it solves: a claim about biology can be supported by evidence of
fundamentally different kinds -- what people have published, what a sequence
model predicts, what a population database shows, what a functional screen
measured. These are not interchangeable. They fail in different ways, they are
wrong about different things, and averaging them into one score destroys the
information that makes them jointly useful.

The insight this module is built on: evidence from independent epistemic
sources is worth far more when it AGREES than the sum of its parts, and worth
far more when it DISAGREES than either part alone. A literature report of exon
skipping plus a SpliceAI prediction of a disrupted donor site is not "two
pieces of evidence." It is a mechanistic hypothesis with a sequence-level
explanation. Conversely, a strong literature consensus that the sequence model
flatly contradicts is either a model failure or a field-wide error -- and
either way it is the most interesting thing on the page.

So we do not fuse. We RECONCILE, and we report the reconciliation.

Epistemic independence classes
------------------------------
Two pieces of evidence are independent only if their errors are uncorrelated.
Twenty papers from a field with a shared assumption are NOT independent. A
sequence model trained on the same gene annotations that guided those papers is
only partly independent. We track this explicitly:

  LITERATURE   -- what has been observed and published
  SEQUENCE     -- what a model predicts from primary sequence
  POPULATION   -- what human genetic variation reveals (constraint, absence)
  FUNCTIONAL   -- what direct perturbation screens measured
  STRUCTURAL   -- what protein/RNA structure implies
  CLINICAL     -- what patient phenotypes show

Sources within a class share failure modes. Across classes they largely do not.
The convergence score rewards cross-class agreement and is deliberately
unimpressed by twenty concordant literature reports.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Any


class EvidenceClass(str, Enum):
    """Epistemic source classes. Errors correlate within, not across."""

    LITERATURE = "literature"
    SEQUENCE = "sequence"
    POPULATION = "population"
    FUNCTIONAL = "functional"
    STRUCTURAL = "structural"
    CLINICAL = "clinical"
    EXPRESSION = "expression"


# How much a class can be trusted on its own, before convergence. These are
# not quality judgements -- they encode how directly each class bears on a
# causal claim about a molecular mechanism.
CLASS_PRIOR: dict[EvidenceClass, float] = {
    EvidenceClass.FUNCTIONAL: 1.00,   # direct perturbation, most causal
    EvidenceClass.CLINICAL: 0.85,     # real phenotypes, confounded
    EvidenceClass.LITERATURE: 0.75,   # observed, but publication-filtered
    EvidenceClass.POPULATION: 0.70,   # strong for constraint, indirect for mechanism
    EvidenceClass.SEQUENCE: 0.65,     # predictive, not observational
    EvidenceClass.EXPRESSION: 0.60,   # correlational by nature
    EvidenceClass.STRUCTURAL: 0.55,   # inferential, often modelled
}


class Concordance(str, Enum):
    CONVERGENT = "convergent"        # independent classes agree
    CONSISTENT = "consistent"        # agree, but not independent
    SILENT = "silent"                # one class has nothing to say
    TENSION = "tension"              # partial disagreement
    CONTRADICTORY = "contradictory"  # direct conflict across classes


@dataclass
class EvidenceItem:
    """One piece of evidence bearing on one claim."""

    source_class: EvidenceClass
    source_name: str            # "SpliceAI", "Europe PMC", "gnomAD"
    claim_direction: str        # "supports", "refutes", "neutral"
    magnitude: float | None = None      # effect size in the source's own units
    magnitude_unit: str | None = None
    confidence: float = 0.5             # 0-1, source-reported or derived
    n_independent: int = 1              # independent observations within source
    provenance: str = ""                # identifier for auditing
    caveats: list[str] = field(default_factory=list)
    raw: dict[str, Any] = field(default_factory=dict)

    def weight(self) -> float:
        """Prior-scaled confidence, with sublinear credit for replication.

        Sublinear because the tenth paper from the same field adds much less
        than the second: they share assumptions, reagents, and reviewers.
        """
        prior = CLASS_PRIOR.get(self.source_class, 0.5)
        replication = 1.0 + 0.4 * math.log1p(max(0, self.n_independent - 1))
        return prior * self.confidence * replication


@dataclass
class ReconciliationResult:
    """The output: not a fused score, but a structured account of agreement."""

    claim: str
    items: list[EvidenceItem] = field(default_factory=list)

    concordance: Concordance = Concordance.SILENT
    n_classes_supporting: int = 0
    n_classes_refuting: int = 0
    classes_present: list[str] = field(default_factory=list)
    classes_absent: list[str] = field(default_factory=list)

    convergence_score: float = 0.0
    support_weight: float = 0.0
    refute_weight: float = 0.0

    narrative: str = ""
    tensions: list[str] = field(default_factory=list)
    what_would_settle_it: list[str] = field(default_factory=list)

    def is_actionable(self) -> bool:
        """Enough cross-class support to justify bench time."""
        return (
            self.concordance == Concordance.CONVERGENT
            and self.n_classes_supporting >= 2
            and self.convergence_score >= 1.5
        )

    def is_interesting_disagreement(self) -> bool:
        """Disagreement between well-supported independent classes.

        These are flagged prominently because they are where the field is
        actually wrong about something, and therefore where a well-designed
        experiment has the most value.
        """
        return (
            self.concordance in (Concordance.CONTRADICTORY, Concordance.TENSION)
            and self.n_classes_supporting >= 1
            and self.n_classes_refuting >= 1
        )


class EvidenceReconciler:
    """Reconciles evidence across epistemic classes for a single claim."""

    def __init__(self, min_classes_for_convergence: int = 2):
        self.min_classes = min_classes_for_convergence

    def reconcile(self, claim: str, items: list[EvidenceItem]) -> ReconciliationResult:
        res = ReconciliationResult(claim=claim, items=list(items))
        if not items:
            res.narrative = "No evidence of any kind bears on this claim."
            res.classes_absent = [c.value for c in EvidenceClass]
            return res

        # Group by class; within-class evidence is pooled, not counted twice.
        by_class: dict[EvidenceClass, list[EvidenceItem]] = {}
        for it in items:
            by_class.setdefault(it.source_class, []).append(it)

        supporting: set[EvidenceClass] = set()
        refuting: set[EvidenceClass] = set()

        for cls, group in by_class.items():
            sup = sum(i.weight() for i in group if i.claim_direction == "supports")
            ref = sum(i.weight() for i in group if i.claim_direction == "refutes")
            res.support_weight += sup
            res.refute_weight += ref
            # A class counts as taking a position only if it leans clearly.
            if sup > ref * 1.5 and sup > 0.3:
                supporting.add(cls)
            elif ref > sup * 1.5 and ref > 0.3:
                refuting.add(cls)

        res.n_classes_supporting = len(supporting)
        res.n_classes_refuting = len(refuting)
        res.classes_present = sorted(c.value for c in by_class)
        res.classes_absent = sorted(
            c.value for c in EvidenceClass if c not in by_class
        )

        res.concordance = self._classify(supporting, refuting, by_class)
        res.convergence_score = self._score(supporting, refuting, by_class)
        res.tensions = self._describe_tensions(by_class, supporting, refuting)
        res.narrative = self._narrate(res, supporting, refuting)
        res.what_would_settle_it = self._next_steps(res, by_class)
        return res

    # ---------------- internals ----------------

    @staticmethod
    def _classify(
        supporting: set, refuting: set, by_class: dict
    ) -> Concordance:
        if supporting and refuting:
            # Both directions represented across independent classes.
            return (
                Concordance.CONTRADICTORY
                if len(supporting) >= 1 and len(refuting) >= 1
                and min(len(supporting), len(refuting)) >= 1
                and abs(len(supporting) - len(refuting)) <= 1
                else Concordance.TENSION
            )
        if len(supporting) >= 2:
            return Concordance.CONVERGENT
        if len(supporting) == 1:
            return Concordance.CONSISTENT
        if len(refuting) >= 2:
            return Concordance.CONVERGENT  # convergent refutation is convergence
        if len(refuting) == 1:
            return Concordance.CONSISTENT
        return Concordance.SILENT

    def _score(self, supporting: set, refuting: set, by_class: dict) -> float:
        """Convergence score.

        Deliberately superlinear in the NUMBER OF CLASSES and only logarithmic
        in the weight within a class. Two independent classes agreeing beats
        one class with twice the evidence -- which is the entire point.
        """
        if not supporting and not refuting:
            return 0.0
        dominant = supporting if len(supporting) >= len(refuting) else refuting
        if not dominant:
            return 0.0

        mass = 0.0
        for cls in dominant:
            group = by_class.get(cls, [])
            w = sum(i.weight() for i in group)
            mass += math.log1p(w)

        breadth = len(dominant) ** 1.5
        # Contradiction discounts the score without erasing it.
        opposition = len(refuting if dominant is supporting else supporting)
        penalty = 1.0 / (1.0 + 0.8 * opposition)
        return round(mass * breadth * penalty, 3)

    @staticmethod
    def _describe_tensions(
        by_class: dict, supporting: set, refuting: set
    ) -> list[str]:
        out = []
        for cls in refuting:
            for other in supporting:
                out.append(
                    f"{cls.value} evidence refutes what {other.value} evidence "
                    f"supports."
                )
        # Within-class disagreement is also worth surfacing.
        for cls, group in by_class.items():
            sup = [i for i in group if i.claim_direction == "supports"]
            ref = [i for i in group if i.claim_direction == "refutes"]
            if sup and ref:
                out.append(
                    f"{cls.value} evidence is internally split "
                    f"({len(sup)} supporting, {len(ref)} refuting)."
                )
        return out

    @staticmethod
    def _narrate(
        res: ReconciliationResult, supporting: set, refuting: set
    ) -> str:
        sup_names = ", ".join(sorted(c.value for c in supporting))
        ref_names = ", ".join(sorted(c.value for c in refuting))

        if res.concordance == Concordance.CONVERGENT and supporting:
            return (
                f"Independent lines of evidence converge: {sup_names} all "
                f"support this. Because these classes fail in different ways, "
                f"agreement between them is substantially stronger than the "
                f"same weight of evidence from any one of them."
            )
        if res.concordance == Concordance.CONVERGENT and refuting:
            return (
                f"Independent lines of evidence converge AGAINST this claim "
                f"({ref_names}). Treat it as likely false rather than merely "
                f"unsupported."
            )
        if res.concordance == Concordance.CONTRADICTORY:
            return (
                f"Direct contradiction: {sup_names} support this while "
                f"{ref_names} refute it. One of these is wrong, and finding "
                f"out which is more informative than either taken alone. This "
                f"is a high-value experiment."
            )
        if res.concordance == Concordance.TENSION:
            return (
                f"Partial disagreement between {sup_names} and {ref_names}. "
                f"Often means the effect is conditional on something neither "
                f"line of evidence has controlled for."
            )
        if res.concordance == Concordance.CONSISTENT:
            only = sup_names or ref_names
            return (
                f"Supported by {only} alone. No independent line of evidence "
                f"either corroborates or contradicts it, so the claim rests "
                f"entirely on that source's assumptions."
            )
        return "No evidence class takes a clear position on this claim."

    @staticmethod
    def _next_steps(res: ReconciliationResult, by_class: dict) -> list[str]:
        """Concrete suggestions for resolving the current epistemic state.

        This is the output a researcher can act on: not "here is a score" but
        "here is the specific gap, and here is what would close it."
        """
        steps: list[str] = []
        absent = set(res.classes_absent)

        if EvidenceClass.SEQUENCE.value in absent:
            steps.append(
                "No sequence-level prediction available. Run the variant "
                "through SpliceAI/Pangolin to see whether the primary sequence "
                "predicts the reported effect."
            )
        if EvidenceClass.FUNCTIONAL.value in absent:
            steps.append(
                "No direct functional perturbation data. A CRISPR or ASO "
                "experiment in your own line would supply the missing causal "
                "class -- this is usually the highest-value next experiment."
            )
        if EvidenceClass.POPULATION.value in absent:
            steps.append(
                "No population-genetics constraint data. Check gnomAD "
                "constraint scores and whether the variant class is depleted "
                "in healthy populations."
            )
        if res.concordance == Concordance.CONTRADICTORY:
            steps.append(
                "Resolve the contradiction directly: identify what differs "
                "between the conflicting sources (cell model, maturity, assay) "
                "and design an experiment that varies exactly that."
            )
        if res.concordance == Concordance.CONSISTENT:
            steps.append(
                "Seek corroboration from a second, independent evidence class "
                "before building on this."
            )
        return steps


# --------------------------------------------------------------------------
# Convenience constructors for each source type
# --------------------------------------------------------------------------


def from_literature_association(assoc, direction: str = "supports") -> EvidenceItem:
    """Build an EvidenceItem from a PathwayAssociation."""
    caveats = list(getattr(assoc, "caveats", lambda: [])())
    return EvidenceItem(
        source_class=EvidenceClass.LITERATURE,
        source_name="published literature",
        claim_direction=direction,
        magnitude=getattr(assoc, "enrichment", None),
        magnitude_unit="odds_ratio",
        confidence=_tier_confidence(getattr(assoc, "confidence_tier", lambda: "")()),
        n_independent=getattr(assoc, "n_labs", 1),
        provenance=f"{getattr(assoc, 'n_papers', 0)} papers",
        caveats=caveats,
    )


def _tier_confidence(tier: str) -> float:
    return {
        "strong": 0.9, "moderate": 0.7, "weak": 0.45,
        "single-lab": 0.3, "insufficient": 0.15,
    }.get(tier, 0.4)


def from_splice_prediction(pred, threshold: float = 0.5) -> EvidenceItem:
    """Build an EvidenceItem from a splice-model prediction."""
    score = getattr(pred, "max_delta", 0.0) or 0.0
    return EvidenceItem(
        source_class=EvidenceClass.SEQUENCE,
        source_name=getattr(pred, "model_name", "splice model"),
        claim_direction="supports" if score >= threshold else "refutes",
        magnitude=score,
        magnitude_unit="delta_score",
        # Sequence models are well calibrated in the middle of their range and
        # overconfident at the extremes; we compress accordingly.
        confidence=min(0.85, 0.35 + 0.5 * score),
        n_independent=1,
        provenance=getattr(pred, "variant_id", ""),
        caveats=list(getattr(pred, "caveats", [])),
    )


def from_constraint(gene_symbol: str, pli: float | None,
                    loeuf: float | None) -> EvidenceItem:
    """Build an EvidenceItem from gnomAD constraint metrics."""
    caveats = []
    if pli is None and loeuf is None:
        caveats.append("No constraint data available for this gene.")
        conf = 0.1
        direction = "neutral"
    else:
        # LOEUF < 0.35 indicates strong selection against LoF.
        constrained = (loeuf is not None and loeuf < 0.35) or (
            pli is not None and pli > 0.9
        )
        direction = "supports" if constrained else "neutral"
        conf = 0.75 if constrained else 0.4
        if loeuf is not None and loeuf > 1.0:
            direction = "refutes"
            caveats.append(
                "Gene tolerates loss of function in healthy populations, "
                "which argues against a strong haploinsufficiency mechanism."
            )
    return EvidenceItem(
        source_class=EvidenceClass.POPULATION,
        source_name="gnomAD constraint",
        claim_direction=direction,
        magnitude=loeuf,
        magnitude_unit="LOEUF",
        confidence=conf,
        provenance=gene_symbol,
        caveats=caveats,
        raw={"pLI": pli, "LOEUF": loeuf},
    )
