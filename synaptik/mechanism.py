"""
Integrated mechanism analysis.

The top-level question this answers: given a gene, a variant, or a phenotype,
what does ALL available evidence -- published, predicted, population,
functional, clinical -- jointly say, and where do the sources disagree?

This is the piece that makes SYNAPTIK v2 different from a literature-mining
tool. Instead of asking "what has been published about X", it asks "what is
true about X, according to every independent line of evidence, and what would
resolve the remaining uncertainty".

Output design principle: the most valuable thing this can produce is not a
confident answer. It is a well-characterised disagreement, or a confident
prediction nobody has tested. Those are surfaced above corroborated findings,
because corroborated findings are what you would have found anyway.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

from .evidence.reconcile import (
    Concordance, EvidenceClass, EvidenceItem, EvidenceReconciler,
    ReconciliationResult, from_literature_association,
)
from .infer.splicing import SplicingAnalyzer
from .omics.evidence_sources import (
    ClinicalEvidence, ClinVarClient, DepMapClient, ExpressionClient,
    FunctionalEvidence, GeneConstraint, GnomADClient,
)
from .schema import Direction
from .sequence.reconcile_splice import (
    SpliceLiteratureReconciler, SpliceReconciliation, VariantPanel,
)
from .sequence.splice import SpliceEnsemble

log = logging.getLogger(__name__)


@dataclass
class MechanismHypothesis:
    """One candidate mechanism, with all supporting and refuting evidence."""

    mechanism: str
    gene_symbol: str | None = None
    pathway: str | None = None

    reconciliation: ReconciliationResult | None = None
    evidence_items: list[EvidenceItem] = field(default_factory=list)

    constraint: GeneConstraint | None = None
    clinical: ClinicalEvidence | None = None
    functional: FunctionalEvidence | None = None
    splice: SpliceReconciliation | None = None

    priority: float = 0.0
    verdict: str = ""
    experiments: list[str] = field(default_factory=list)

    def concordance(self) -> str:
        return (self.reconciliation.concordance.value
                if self.reconciliation else "unknown")

    def n_independent_classes(self) -> int:
        return (self.reconciliation.n_classes_supporting
                if self.reconciliation else 0)


@dataclass
class GeneDossier:
    """Everything the system knows about one gene, across all evidence classes."""

    gene_symbol: str

    constraint: GeneConstraint | None = None
    clinical: ClinicalEvidence | None = None
    functional: FunctionalEvidence | None = None
    expression: Any = None

    n_literature_records: int = 0
    n_papers: int = 0
    n_labs: int = 0
    literature_phenotypes: list[str] = field(default_factory=list)
    splicing_events: int = 0
    implicated_regulators: list[str] = field(default_factory=list)

    hypotheses: list[MechanismHypothesis] = field(default_factory=list)
    contradictions: list[str] = field(default_factory=list)
    gaps: list[str] = field(default_factory=list)
    summary: str = ""

    def evidence_classes_present(self) -> list[str]:
        out = []
        if self.n_literature_records:
            out.append("literature")
        if self.constraint and not self.constraint.error:
            out.append("population")
        if self.clinical and not self.clinical.error:
            out.append("clinical")
        if self.functional and not self.functional.error:
            out.append("functional")
        if self.expression and not getattr(self.expression, "error", None):
            out.append("expression")
        if self.splicing_events:
            out.append("sequence")
        return out


class MechanismAnalyst:
    """Integrates every evidence class for a gene, variant, or phenotype."""

    def __init__(
        self,
        store,
        session=None,
        gnomad: GnomADClient | None = None,
        clinvar: ClinVarClient | None = None,
        depmap: DepMapClient | None = None,
        expression: ExpressionClient | None = None,
        splice_ensemble: SpliceEnsemble | None = None,
        offline: bool = False,
    ):
        self.store = store
        self.session = session
        self.offline = offline
        self.gnomad = gnomad
        self.clinvar = clinvar
        self.depmap = depmap
        self.expression = expression
        self.splice = splice_ensemble or SpliceEnsemble()
        self.reconciler = EvidenceReconciler()
        self.splice_reconciler = SpliceLiteratureReconciler()
        self.splicing_analyzer = SplicingAnalyzer()

    # ---------------- gene dossier ----------------

    def gene_dossier(self, gene_symbol: str) -> GeneDossier:
        """Assemble every evidence class for one gene."""
        d = GeneDossier(gene_symbol=gene_symbol.upper())
        self._load_literature(d)
        self._load_external(d)
        d.hypotheses = self._build_hypotheses(d)
        d.contradictions = self._find_contradictions(d)
        d.gaps = self._find_gaps(d)
        d.summary = self._summarize(d)
        return d

    def _load_literature(self, d: GeneDossier) -> None:
        g = d.gene_symbol
        rows = self.store.conn.execute(
            """SELECT DISTINCT r.record_id FROM records r
               LEFT JOIN record_genes rg ON rg.record_id = r.record_id
               WHERE UPPER(rg.gene_label) = ?
                  OR UPPER(r.perturbation_target_key) LIKE ?
                  OR r.record_id IN (
                       SELECT record_id FROM splicing_events
                       WHERE UPPER(gene_key) LIKE ?)""",
            (g, f"%{g.lower()}%", f"%{g.lower()}%"),
        ).fetchall()
        ids = [r["record_id"] for r in rows]
        if not ids:
            return

        records = [self.store.get_record(i) for i in ids]
        records = [r for r in records if r]
        d.n_literature_records = len(records)
        d.n_papers = len({r.source.source_id for r in records})
        d.n_labs = len({r.source.lab_id for r in records if r.source.lab_id})
        d.literature_phenotypes = sorted({
            r.outcome.phenotype.ontology_label or r.outcome.phenotype.raw
            for r in records
        })[:12]

        events = self.splicing_analyzer.collect_events(records)
        d.splicing_events = len(events)
        if events:
            labs = {r.source.source_id: r.source.lab_id for r in records}
            regs = self.splicing_analyzer.infer_regulators(events, records)
            d.implicated_regulators = [r.rbp for r in regs]

    def _load_external(self, d: GeneDossier) -> None:
        if self.offline:
            return
        g = d.gene_symbol
        if self.gnomad:
            try:
                d.constraint = self.gnomad.constraint(g)
            except Exception as e:
                log.warning("gnomAD lookup failed for %s: %s", g, e)
        if self.clinvar:
            try:
                d.clinical = self.clinvar.gene_summary(g)
            except Exception as e:
                log.warning("ClinVar lookup failed for %s: %s", g, e)
        if self.depmap:
            try:
                d.functional = self.depmap.gene_effect(g)
            except Exception as e:
                log.warning("DepMap lookup failed for %s: %s", g, e)
        if self.expression:
            try:
                d.expression = self.expression.profile(g)
            except Exception as e:
                log.warning("Expression lookup failed for %s: %s", g, e)

    def _build_hypotheses(self, d: GeneDossier) -> list[MechanismHypothesis]:
        """Reconcile the loss-of-function mechanism across all classes."""
        items: list[EvidenceItem] = []

        if d.n_literature_records:
            items.append(EvidenceItem(
                source_class=EvidenceClass.LITERATURE,
                source_name="published records",
                claim_direction="supports",
                confidence=min(0.85, 0.4 + 0.1 * min(4, d.n_labs)),
                n_independent=max(1, d.n_labs),
                provenance=f"{d.n_papers} papers",
                caveats=(
                    ["All literature evidence comes from a single group."]
                    if d.n_labs <= 1 else []
                ),
            ))
        if d.constraint and not d.constraint.error:
            items.append(d.constraint.to_evidence(claim_is_lof_mechanism=True))
        if d.clinical and not d.clinical.error:
            items.append(d.clinical.to_evidence())
        if d.functional and not d.functional.error:
            items.append(d.functional.to_evidence())
        if d.expression and not getattr(d.expression, "error", None):
            items.append(d.expression.to_evidence())

        hyp = MechanismHypothesis(
            mechanism=f"Loss of {d.gene_symbol} function causes a neuronal phenotype",
            gene_symbol=d.gene_symbol,
            evidence_items=items,
            constraint=d.constraint,
            clinical=d.clinical,
            functional=d.functional,
        )
        hyp.reconciliation = self.reconciler.reconcile(hyp.mechanism, items)
        hyp.priority = hyp.reconciliation.convergence_score
        hyp.verdict = hyp.reconciliation.narrative
        hyp.experiments = list(hyp.reconciliation.what_would_settle_it)
        return [hyp]

    @staticmethod
    def _find_contradictions(d: GeneDossier) -> list[str]:
        """Cross-class conflicts that a literature-only view would never see."""
        out = []
        lit_claims_lof = d.n_literature_records > 0

        if lit_claims_lof and d.constraint and d.constraint.constraint_class() == "tolerant":
            out.append(
                f"The literature reports phenotypes from {d.gene_symbol} "
                f"perturbation, but the gene is tolerant of loss of function "
                f"in healthy populations (LOEUF {d.constraint.loeuf:.2f})"
                + (f", including {d.constraint.n_homozygous_lof} homozygous "
                   f"LoF individuals" if d.constraint.n_homozygous_lof else "")
                + ". Either the published mechanism is not simple "
                "haploinsufficiency, the effect is subtle enough to escape "
                "selection, or it is context-dependent."
            )
        if lit_claims_lof and d.expression and \
                getattr(d.expression, "neuron_expressed", None) is False:
            out.append(
                f"The literature reports neuronal phenotypes for "
                f"{d.gene_symbol}, but the gene does not appear to be "
                f"appreciably expressed in neurons. A cell-autonomous "
                f"mechanism is questionable."
            )
        if d.functional and d.functional.is_common_essential and lit_claims_lof:
            out.append(
                f"{d.gene_symbol} is commonly essential across cell lines. "
                f"Published neuronal phenotypes may be secondary to impaired "
                f"viability rather than a neuron-specific mechanism -- check "
                f"whether the studies controlled for cell death."
            )
        if d.clinical and d.clinical.n_vus > 20 and \
                not d.clinical.has_clinical_support() and lit_claims_lof:
            out.append(
                f"Research literature implicates {d.gene_symbol}, but ClinVar "
                f"contains {d.clinical.n_vus} variants of uncertain "
                f"significance and no pathogenic classifications. Clinical "
                f"genetics has not corroborated the research claim."
            )
        return out

    @staticmethod
    def _find_gaps(d: GeneDossier) -> list[str]:
        present = set(d.evidence_classes_present())
        gaps = []
        if "population" not in present:
            gaps.append(
                "No population-constraint data. gnomAD would indicate whether "
                "loss of this gene is tolerated in healthy humans."
            )
        if "functional" not in present:
            gaps.append(
                "No functional perturbation data. This is the most causal "
                "evidence class and its absence is the main reason claims here "
                "remain correlational."
            )
        if "clinical" not in present:
            gaps.append(
                "No clinical genetics data. ClinVar/OMIM would show whether "
                "patient variants support the mechanism."
            )
        if "sequence" not in present:
            gaps.append(
                "No sequence-level prediction. If you have specific variants, "
                "running them through the splice ensemble adds an independent "
                "evidence class."
            )
        if d.n_labs <= 1 and d.n_literature_records:
            gaps.append(
                "All literature evidence comes from one research group; no "
                "independent replication exists in the corpus."
            )
        return gaps

    @staticmethod
    def _summarize(d: GeneDossier) -> str:
        classes = d.evidence_classes_present()
        if not classes:
            return (
                f"No evidence of any kind for {d.gene_symbol} in the current "
                f"knowledge base."
            )
        parts = [
            f"{d.gene_symbol}: evidence from {len(classes)} class(es) "
            f"({', '.join(classes)})."
        ]
        if d.n_literature_records:
            parts.append(
                f"{d.n_literature_records} records from {d.n_papers} papers "
                f"across {d.n_labs} group(s)."
            )
        if d.constraint and not d.constraint.error:
            parts.append(d.constraint.interpretation())
        if d.contradictions:
            parts.append(
                f"{len(d.contradictions)} cross-class contradiction(s) found -- "
                f"see below; these are the most informative findings here."
            )
        elif len(classes) >= 3:
            parts.append(
                "Multiple independent evidence classes are available and none "
                "contradict each other."
            )
        return " ".join(parts)

    # ---------------- variant analysis ----------------

    def analyze_variants(
        self,
        variants: list[tuple[str, int, str, str]],
        gene_symbol: str | None = None,
    ) -> VariantPanel:
        """Reconcile a set of variants against the literature.

        variants: list of (chrom, position, ref, alt)
        """
        panel = VariantPanel()

        observed = []
        lab_ids: dict[str, str] = {}
        if gene_symbol:
            rows = self.store.conn.execute(
                "SELECT record_id FROM splicing_events WHERE UPPER(gene_key) LIKE ?",
                (f"%{gene_symbol.lower()}%",),
            ).fetchall()
            recs = [self.store.get_record(r["record_id"]) for r in rows]
            recs = [r for r in recs if r]
            observed = self.splicing_analyzer.collect_events(recs)
            lab_ids = {r.source.source_id: r.source.lab_id for r in recs
                       if r.source.lab_id}

        for chrom, pos, ref, alt in variants:
            try:
                pred = self.splice.predict(
                    chrom, pos, ref, alt, gene_symbol=gene_symbol
                )
            except Exception as e:
                log.warning("Splice prediction failed for %s:%s: %s",
                            chrom, pos, e)
                pred = None
            rec = self.splice_reconciler.reconcile(
                pred, observed, gene_symbol=gene_symbol,
                variant_id=f"{chrom}:{pos}:{ref}>{alt}", lab_ids=lab_ids,
            )
            panel.reconciliations.append(rec)
        return panel

    # ---------------- discovery ----------------

    def find_untested_predictions(
        self, gene_symbols: list[str], variants_by_gene: dict[str, list],
    ) -> list[SpliceReconciliation]:
        """Surface confident predictions with no published observation.

        This is the discovery mode: it looks for places where sequence models
        make a confident claim that nobody has tested. Those are, by
        construction, the experiments with the best ratio of novelty to risk.
        """
        out: list[SpliceReconciliation] = []
        for gene in gene_symbols:
            variants = variants_by_gene.get(gene, [])
            if not variants:
                continue
            panel = self.analyze_variants(variants, gene_symbol=gene)
            out.extend(panel.novel_hypotheses())
        return sorted(out, key=lambda r: -r.priority)
