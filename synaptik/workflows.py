"""
Shared entry points for the lab-facing workflows.

The CLI and the web interface both call these, so a deposit, a design brief or
a benchmark means exactly the same thing however it is reached. Each workflow
has a `*_to_dict` serializer used for `--json` output and for the HTTP API.
"""

from __future__ import annotations

from .analyst import Analyst
from .infer.benchmark import BenchmarkResult, RetrospectiveBenchmark
from .infer.design import DesignAdvisor, DesignBrief
from .infer.splicing import NormalizedSplicingEvent, SplicingAnalyzer
from .labdata import (
    LabDataReconciler, ReconciledPhenotype, SplicingDepositReport,
    SplicingTable, parse_phenotype_table, reconcile_phenotypes,
)
from .schema import ExperimentQuery
from .store import Store


def corpus_splicing_events(
    store: Store,
) -> tuple[list[NormalizedSplicingEvent], dict[str, str]]:
    """Every published splicing event in the knowledge base, plus lab ids."""
    rows = store.conn.execute("SELECT record_id FROM splicing_events").fetchall()
    records = [store.get_record(r["record_id"]) for r in rows]
    records = [r for r in records if r and not r.source.is_retracted]
    events = SplicingAnalyzer().collect_events(records)
    labs = {r.source.source_id: r.source.lab_id for r in records
            if r.source.lab_id}
    return events, labs


def deposit_splicing(
    store: Store, table: SplicingTable,
    max_fdr: float = 0.05, min_abs_dpsi: float = 0.1,
) -> SplicingDepositReport:
    events, labs = corpus_splicing_events(store)
    return LabDataReconciler(max_fdr, min_abs_dpsi).reconcile_splicing(
        table, events, labs)


def deposit_phenotypes(
    store: Store, source: str, text: str | None = None
) -> tuple[list[ReconciledPhenotype], list[str]]:
    observations, warnings = parse_phenotype_table(source, text=text)
    return reconcile_phenotypes(observations, store), warnings


def design_brief(
    store: Store, query: ExperimentQuery, cv: float | None = None,
    rho: float = 0.7, power: float = 0.8, alpha: float = 0.05,
) -> DesignBrief:
    lookup = Analyst(store).lookup_experiment(query, max_results=200)
    return DesignAdvisor(alpha=alpha, power=power).brief(lookup, cv=cv, rho=rho)


from .sequence.benchmark_spliceai import (
    SpliceBenchmarkResult, run_spliceai_benchmark,
)
from .drug.aso import ASOCandidate, ASODesigner, ASODesignResult, ASOMode


def run_benchmark(store: Store, cutoff_year: int) -> BenchmarkResult:
    return RetrospectiveBenchmark().run(list(store.iter_records()), cutoff_year)


def run_splice_benchmark(threshold: float = 0.20) -> SpliceBenchmarkResult:
    return run_spliceai_benchmark(threshold=threshold)


def design_aso(
    target_sequence: str,
    mode: str = "splice_switch",
    target_gene: str = "",
    target_region: str = "splice site",
    top_k: int = 5,
) -> ASODesignResult:
    """In silico design of synthesis-ready therapeutic antisense oligonucleotides."""
    aso_mode = ASOMode.GAPMER if mode.lower() == "gapmer" else ASOMode.SPLICE_SWITCH
    return ASODesigner().design(
        target_sequence=target_sequence,
        mode=aso_mode,
        target_gene=target_gene,
        target_region=target_region,
        top_k=top_k,
    )


# --------------------------------------------------------------------------
# serializers
# --------------------------------------------------------------------------


def splicing_report_to_dict(rep: SplicingDepositReport, limit: int = 200) -> dict:
    return {
        "summary": rep.summary,
        "format": rep.table.source_format,
        "n_events": rep.n_events,
        "n_genes_tested": rep.n_genes_tested,
        "n_significant": rep.n_significant,
        "n_rows_skipped": rep.table.n_rows_skipped,
        "skip_reasons": rep.table.skip_reasons,
        "thresholds": rep.thresholds,
        "corpus_events": rep.corpus_events,
        "corpus_genes": rep.corpus_genes,
        "counts": rep.counts(),
        "warnings": rep.warnings,
        "events": [{
            "gene": r.event.gene,
            "event_type": r.event.event_type.value,
            "delta_psi": r.event.delta_psi,
            "fdr": r.event.fdr,
            "coordinates": r.event.coordinates,
            "exon_length": r.event.exon_length,
            "call": r.call.value,
            "priority": round(r.priority, 3),
            "n_papers": r.n_papers,
            "n_labs": r.n_labs,
            "literature_directions": r.literature_directions,
            "literature_event_types": r.literature_event_types,
            "regulators_reported": r.regulators_reported,
            "source_ids": r.source_ids,
            "nmd_risk": r.nmd_risk,
            "nmd_note": r.nmd_note,
            "aso_feasible": r.aso_feasible,
            "aso_approach": r.aso_approach,
            "note": r.note,
        } for r in rep.reconciled[:limit]],
        "not_reproduced": [{
            "gene": n.gene, "n_papers": n.n_papers, "n_labs": n.n_labs,
            "literature_direction": n.literature_direction,
            "regulators": n.regulators,
            "your_largest_dpsi": n.user_best_dpsi,
            "your_fdr": n.user_best_fdr,
            "source_ids": n.source_ids,
        } for n in rep.not_reproduced[:limit]],
        "regulators": [{
            "rbp": g.rbp,
            "targets_in_corpus": g.n_targets_in_corpus,
            "targets_tested": g.n_targets_tested,
            "targets_changed": g.n_targets_changed,
            "overlap_genes": g.overlap_genes,
            "odds_ratio": g.odds_ratio,
            "p_value": g.p_value,
            "q_value": g.q_value,
            "direction_agreement": g.direction_agreement,
            "supporting_labs": g.n_supporting_labs,
            "canonical_effect": g.canonical_effect,
            "caveats": g.caveats(),
        } for g in rep.regulators],
    }


def phenotype_report_to_dict(
    rows: list[ReconciledPhenotype], warnings: list[str]
) -> dict:
    return {
        "warnings": warnings,
        "observations": [{
            "phenotype": r.observation.phenotype,
            "resolved_label": r.resolved_label,
            "gene": r.observation.gene,
            "your_direction": r.observation.direction.value,
            "call": r.call,
            "n_records": r.n_records,
            "n_labs": r.n_labs,
            "literature_directions": r.literature_directions,
            "note": r.note,
            "source_ids": r.source_ids,
        } for r in rows],
    }


def design_brief_to_dict(b: DesignBrief) -> dict:
    return {
        "basis": b.basis,
        "value": b.value_statement,
        "n_exact": b.n_exact, "n_close": b.n_close, "n_labs": b.n_labs,
        "protocol": {
            "div_median": b.div_median, "div_range": b.div_range,
            "differentiation": b.differentiation_protocols,
            "coculture": b.cocultures, "markers": b.markers,
            "statistical_tests": b.statistical_tests,
            "n_replicates_median": b.n_replicates_median,
            "n_replicates_range": b.n_replicates_range,
            "isogenic_fraction": b.isogenic_fraction,
        },
        "effects": {
            "consensus_direction": b.consensus_direction,
            "direction_counts": b.direction_counts,
            "n_quantified": len(b.effect_sizes_percent),
            "median_effect_percent": b.median_effect_percent,
            "n_null_results": b.n_null_results,
        },
        "power": [{
            "scenario": p.label, "effect_percent": p.effect_percent,
            "cv": p.cv, "cohens_d": p.cohens_d, "n_per_group": p.n_per_group,
            "n_paired": p.n_paired, "paired_rho": p.paired_rho,
        } for p in b.power],
        "recommended_n_per_group": b.recommended_n,
        "recommended_n_paired": b.recommended_n_paired,
        "paired_rho": b.paired_rho,
        "recommended_n_basis": b.recommended_n_basis,
        "impractical_unpaired": b.impractical,
        "detectable_effects": [
            {"n_per_group": n, "cv": cv, "min_effect_percent": eff}
            for n, cv, eff in b.detectable],
        "weaknesses": b.weaknesses,
        "recommendations": b.recommendations,
        "pitfalls": b.pitfalls,
        "orthogonal_assays": b.orthogonal_assays,
        "key_papers": b.key_papers,
    }


def benchmark_to_dict(res: BenchmarkResult) -> dict:
    return {
        "cutoff_year": res.cutoff_year,
        "n_records_before": res.n_records_before,
        "n_records_after": res.n_records_after,
        "n_phenotypes": res.n_phenotypes,
        "n_associations": len(res.cases),
        "n_replicated": res.n_replicated,
        "underpowered": res.underpowered,
        "verdict": res.verdict,
        "best_baseline": res.best_baseline,
        "difference_vs_best_baseline": res.difference_vs_best_baseline,
        "difference_ci95": res.difference_ci95,
        "rankers": [{"name": r.name, "auroc": r.auroc,
                     "description": r.description} for r in res.rankers],
        "tiers": res.tier_table,
        "caveats": res.caveats,
    }


def spliceai_benchmark_to_dict(res: SpliceBenchmarkResult) -> dict:
    return {
        "n_cases": res.n_cases,
        "n_positives": res.n_positives,
        "n_negatives": res.n_negatives,
        "threshold": res.threshold,
        "sensitivity": res.sensitivity,
        "specificity": res.specificity,
        "precision": res.precision,
        "accuracy": res.accuracy,
        "auroc": res.auroc,
        "details": res.details,
    }


def aso_result_to_dict(res: ASODesignResult) -> dict:
    from dataclasses import asdict
    return {
        "target_gene": res.target_gene,
        "target_region": res.target_region,
        "target_sequence": res.target_sequence,
        "mode": res.mode.value,
        "chemistry": res.chemistry.value,
        "n_windows_evaluated": res.n_windows_evaluated,
        "n_passed_filters": res.n_passed_filters,
        "precedent": res.precedent,
        "caveats": res.caveats,
        "candidates": [asdict(c) for c in res.candidates],
    }
