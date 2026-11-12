"""
High-level analysis API.

One object that the CLI, the web UI, and any notebook all call. Keeping the
orchestration here rather than in the interfaces means the reasoning is
identical regardless of how you reach it -- there is one definition of "what
does this phenotype mean", not three that drift apart.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

from .drug.repurposing import (
    DrugCandidate, MoADirection, RepurposingEngine, curated_targets,
    required_moa_direction,
)
from .infer.association import (
    AssociationEngine, PathwayAssociation, phenotype_coherence,
)
from .infer.lookup import ExperimentLookup, LookupResult, MatchTier
from .infer.splicing import (
    RegulatorHypothesis, SplicingAnalyzer, aso_strategy, nmd_likelihood,
)
from .infer.weighting import EvidenceWeighter
from .normalize.matcher import PhenotypeMatcher
from .normalize.phenotype_ontology import BY_ID, ancestors, commensurable
from .schema import Direction, ExperimentQuery, ExperimentRecord
from .store import Store

log = logging.getLogger(__name__)


@dataclass
class PhenotypeAnalysis:
    """Complete answer to 'what does this phenotype mean?'"""

    phenotype_query: str
    resolved_phenotype: str | None = None
    resolved_label: str | None = None
    direction: Direction = Direction.UNSPECIFIED

    n_matching_records: int = 0
    n_papers: int = 0
    n_labs: int = 0
    corpus_size: int = 0

    associations: list[PathwayAssociation] = field(default_factory=list)
    regulators: list[RegulatorHypothesis] = field(default_factory=list)
    coherence: dict = field(default_factory=dict)
    drug_candidates: list[DrugCandidate] = field(default_factory=list)

    warnings: list[str] = field(default_factory=list)
    interpretation: str = ""

    def top(self, n: int = 10) -> list[PathwayAssociation]:
        return self.associations[:n]

    def strong_only(self) -> list[PathwayAssociation]:
        return [a for a in self.associations
                if a.confidence_tier() in ("strong", "moderate")]


class Analyst:
    """The main analysis interface."""

    def __init__(
        self,
        store: Store,
        weighter: EvidenceWeighter | None = None,
        repurposing: RepurposingEngine | None = None,
    ):
        self.store = store
        self.weighter = weighter or EvidenceWeighter()
        self.matcher = PhenotypeMatcher()
        self.assoc_engine = AssociationEngine(self.weighter)
        self.lookup_engine = ExperimentLookup(self.weighter)
        self.splicing = SplicingAnalyzer()
        self.repurposing = repurposing

    # ---------------- component 1 ----------------

    def lookup_experiment(
        self, query: ExperimentQuery, max_results: int = 60
    ) -> LookupResult:
        """Has this experiment been done? Restricts the corpus intelligently.

        Loading every record for every lookup does not scale past a few hundred
        thousand records, so we pre-filter in SQL on the highest-selectivity
        fields before scoring in Python.
        """
        clauses, params = [], []
        if query.perturbation_target:
            clauses.append(
                "(perturbation_target_key LIKE ? OR record_id IN "
                "(SELECT record_id FROM record_genes WHERE gene_key LIKE ?))"
            )
            like = f"%{query.perturbation_target.lower()}%"
            params += [like, like]
        if query.variant_hgvs:
            clauses.append("variant_hgvs = ?")
            params.append(query.variant_hgvs)
        if query.assay_class and not clauses:
            clauses.append("assay_class = ?")
            params.append(query.assay_class.value)

        where = " OR ".join(clauses) if clauses else ""
        corpus = list(self.store.iter_records(where, tuple(params)))
        if not corpus:
            corpus = list(self.store.iter_records())
        return self.lookup_engine.lookup(query, corpus, max_results)

    # ---------------- component 2 ----------------

    def analyze_phenotype(
        self,
        phenotype: str,
        direction: Direction = Direction.UNSPECIFIED,
        include_ancestors: bool = True,
        include_drugs: bool = False,
        max_drug_pathways: int = 3,
    ) -> PhenotypeAnalysis:
        """What does this phenotype mean? The main entry point."""
        res = PhenotypeAnalysis(phenotype_query=phenotype, direction=direction)

        term = self.matcher.match(phenotype)
        res.resolved_phenotype = term.ontology_id
        res.resolved_label = term.ontology_label

        keys = self._phenotype_keys(term, phenotype, include_ancestors)
        foreground = self.store.records_for_phenotype(keys)

        if not foreground:
            res.warnings.append(
                f"No records match '{phenotype}'. Either the corpus does not "
                "cover it, or the phenotype vocabulary needs a synonym for how "
                "you phrased it."
            )
            return res

        if direction != Direction.UNSPECIFIED:
            foreground = [
                r for r in foreground
                if r.outcome.direction == direction or r.is_negative_result()
            ]

        all_records = list(self.store.iter_records())
        res.corpus_size = len(all_records)
        res.n_matching_records = len(foreground)
        res.n_papers = len({r.source.source_id for r in foreground})
        res.n_labs = len({r.source.lab_id for r in foreground if r.source.lab_id})

        res.associations = self.assoc_engine.associate(
            foreground, all_records, direction_filter=None
        )

        events = self.splicing.collect_events(foreground)
        if events:
            res.regulators = self.splicing.infer_regulators(events, foreground)

        self._add_warnings(res, foreground)
        res.interpretation = self._interpret(res)

        if include_drugs and self.repurposing:
            res.drug_candidates = self._drug_candidates(
                res, max_pathways=max_drug_pathways
            )
        return res

    def _phenotype_keys(
        self, term, raw: str, include_ancestors: bool
    ) -> list[str]:
        keys: list[str] = []
        if term.ontology_id:
            keys.append(term.ontology_id)
            if include_ancestors:
                keys.extend(ancestors(term.ontology_id))
            # Terms curated into the same pool_group are commensurable and
            # would otherwise fragment the evidence across near-synonyms.
            for tid in BY_ID:
                if tid not in keys and commensurable(term.ontology_id, tid):
                    keys.append(tid)
        keys.append(f"raw:{raw.strip().lower()}")
        return keys

    # ---------------- interpretation ----------------

    @staticmethod
    def _add_warnings(res: PhenotypeAnalysis, foreground: list) -> None:
        if res.n_labs <= 1:
            res.warnings.append(
                "All matching evidence comes from one research group. Nothing "
                "below should be treated as independently established."
            )
        if res.n_matching_records < 5:
            res.warnings.append(
                f"Only {res.n_matching_records} records match. Associations "
                "computed on this little evidence are unstable."
            )
        n_neg = sum(1 for r in foreground if r.is_negative_result())
        if n_neg == 0 and res.n_matching_records >= 10:
            res.warnings.append(
                "No null results anywhere in the matching literature. Either "
                "this phenotype is unusually robust, or negative findings were "
                "not published."
            )
        if not res.resolved_phenotype:
            res.warnings.append(
                "Phenotype did not map to the controlled vocabulary; matching "
                "fell back to exact text. Consider adding a synonym."
            )

    @staticmethod
    def _interpret(res: PhenotypeAnalysis) -> str:
        if not res.associations:
            return (
                f"{res.n_matching_records} records measure this phenotype, but "
                "no pathway is enriched above background. This usually means "
                "the papers describe the phenotype without attributing it to "
                "a mechanism."
            )
        strong = res.strong_only()
        top = res.associations[0]
        lines = [
            f"{res.n_matching_records} records from {res.n_papers} papers "
            f"({res.n_labs} independent groups) measure this phenotype.",
        ]
        if strong:
            names = ", ".join(a.pathway_label for a in strong[:3])
            lines.append(
                f"{len(strong)} pathway association(s) reach moderate confidence "
                f"or better: {names}."
            )
        else:
            lines.append(
                f"No association reaches moderate confidence. The strongest, "
                f"'{top.pathway_label}', is tiered '{top.confidence_tier()}' -- "
                "treat it as a hypothesis to test, not a finding."
            )
        if res.regulators:
            r = res.regulators[0]
            lines.append(
                f"Splicing evidence points to {r.rbp} "
                f"({r.n_supporting_events} events, {r.n_labs} labs)."
            )
        return " ".join(lines)

    # ---------------- drugs ----------------

    def _drug_candidates(
        self, res: PhenotypeAnalysis, max_pathways: int = 3
    ) -> list[DrugCandidate]:
        out: list[DrugCandidate] = []
        pool = res.strong_only() or res.associations[:max_pathways]

        for assoc in pool[:max_pathways]:
            pdir = (
                Direction.DECREASE if assoc.dominant_direction == "decrease"
                else Direction.INCREASE if assoc.dominant_direction == "increase"
                else Direction.UNSPECIFIED
            )
            required = required_moa_direction(res.direction, pdir)
            if required == MoADirection.UNKNOWN:
                continue

            for symbol in curated_targets(assoc.pathway_label)[:5]:
                try:
                    cands = self.repurposing.candidates_for_target(
                        symbol, required,
                        pathway_key=assoc.pathway_key,
                        pathway_label=assoc.pathway_label,
                        evidence_tier=assoc.confidence_tier(),
                    )
                except Exception as e:
                    log.warning("Drug lookup failed for %s: %s", symbol, e)
                    continue
                out.extend(cands[:5])

        if out:
            corpus = list(self.store.iter_records())
            self.repurposing.cross_reference_corpus(out, corpus)
        out.sort(key=lambda c: -c.rank_score())
        return out[:40]

    # ---------------- splicing ----------------

    def analyze_splicing(
        self, gene: str | None = None, regulator: str | None = None
    ) -> dict[str, Any]:
        """Splicing-focused view: events, regulators, and ASO feasibility."""
        clauses, params = [], []
        if gene:
            clauses.append("gene_key LIKE ?")
            params.append(f"%{gene.lower()}%")
        if regulator:
            clauses.append("regulator_key LIKE ?")
            params.append(f"%{regulator.lower()}%")
        where = " AND ".join(clauses) if clauses else "1=1"

        rows = self.store.conn.execute(
            f"SELECT record_id FROM splicing_events WHERE {where}", tuple(params)
        ).fetchall()
        ids = [r["record_id"] for r in rows]
        if not ids:
            return {"events": [], "regulators": [], "message":
                    "No splicing events match that filter in the corpus."}

        records = [self.store.get_record(i) for i in ids]
        records = [r for r in records if r]
        events = self.splicing.collect_events(records)
        regs = self.splicing.infer_regulators(events, records)

        strategies = []
        for e in events[:20]:
            s = aso_strategy(e.event_type, e.direction)
            level, expl = nmd_likelihood(e.event_type, e.direction)
            strategies.append({
                "gene": e.gene,
                "event_type": e.event_type.value,
                "delta_psi": e.delta_psi,
                "magnitude": e.magnitude_bin(),
                "sign_confidence": e.sign_confidence,
                "aso_feasible": s.feasible,
                "aso_approach": s.approach,
                "aso_precedent": s.precedent,
                "nmd_risk": level,
                "nmd_note": expl,
                "source_id": e.source_id,
            })

        return {
            "n_events": len(events),
            "events": strategies,
            "regulators": regs,
            "genes": sorted({e.gene for e in events}),
        }

    # ---------------- cross-checks ----------------

    def compare_phenotypes(self, pheno_a: str, pheno_b: str) -> dict:
        """Do two phenotypes implicate the same biology?

        The practical use: you observe a puncta phenotype and an ephys
        phenotype. If they converge on the same pathways, that is orthogonal
        support. If not, one may be an artifact.
        """
        a = self.analyze_phenotype(pheno_a)
        b = self.analyze_phenotype(pheno_b)
        coh = phenotype_coherence(a.associations, b.associations)
        coh["phenotype_a"] = {
            "query": pheno_a, "n_records": a.n_matching_records,
            "n_labs": a.n_labs,
        }
        coh["phenotype_b"] = {
            "query": pheno_b, "n_records": b.n_matching_records,
            "n_labs": b.n_labs,
        }
        j = coh["jaccard"]
        coh["verdict"] = (
            "Strong convergence -- these phenotypes implicate overlapping "
            "biology, which is genuine orthogonal support."
            if j >= 0.4 else
            "Partial overlap -- some shared pathways, but the profiles differ."
            if j >= 0.15 else
            "Little overlap. Either these are mechanistically distinct, or one "
            "is poorly covered in the corpus."
        )
        return coh

    def corpus_report(self) -> dict:
        """Health check on the knowledge base itself."""
        s = self.store.stats()
        rows = self.store.conn.execute(
            """SELECT phenotype_key, COUNT(*) c FROM records
               GROUP BY phenotype_key ORDER BY c DESC LIMIT 20"""
        ).fetchall()
        top_phen = [(r["phenotype_key"], r["c"]) for r in rows]

        rows = self.store.conn.execute(
            """SELECT pathway_label, COUNT(*) c FROM record_pathways
               GROUP BY pathway_key ORDER BY c DESC LIMIT 20"""
        ).fetchall()
        top_path = [(r["pathway_label"], r["c"]) for r in rows]

        rows = self.store.conn.execute(
            """SELECT lab_id, COUNT(*) c FROM sources WHERE lab_id IS NOT NULL
               GROUP BY lab_id ORDER BY c DESC LIMIT 10"""
        ).fetchall()
        top_labs = [(r["lab_id"], r["c"]) for r in rows]

        unmapped = self.store.conn.execute(
            "SELECT COUNT(*) c FROM records WHERE phenotype_key LIKE 'raw:%'"
        ).fetchone()["c"]

        return {
            **s,
            "unmapped_phenotype_records": unmapped,
            "unmapped_fraction": unmapped / max(1, s["records"]),
            "top_phenotypes": top_phen,
            "top_pathways": top_path,
            "most_prolific_labs": top_labs,
            "negative_result_fraction": s["negative_records"] / max(1, s["records"]),
            "quantified_fraction": s["quantified_records"] / max(1, s["records"]),
        }
