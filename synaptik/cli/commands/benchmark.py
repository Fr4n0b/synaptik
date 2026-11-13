"""`synaptik benchmark`: do the confidence tiers predict later replication?"""

from __future__ import annotations

import json

from ...store import Store
from ...workflows import (
    benchmark_to_dict,
    run_benchmark,
    run_splice_benchmark,
    spliceai_benchmark_to_dict,
)
from ..main import Fmt, header


def cmd_benchmark(args) -> int:
    if getattr(args, "all", False):
        from ...drug.benchmark_aso import run_aso_tiling_benchmark
        from ...twin.models import PatientBaseline
        from ...twin.twin_generator import evaluate_twin_calibration

        threshold = getattr(args, "threshold", 0.20) or 0.20
        splice_res = run_splice_benchmark(threshold=threshold)
        aso_res = run_aso_tiling_benchmark()
        twin_cal = evaluate_twin_calibration(
            PatientBaseline("BENCH-CAL", 56.0, "bulbar", 12.0, 38.0, baseline_slope=-1.1, baseline_plasma_nfl=68.0),
            n_held_out=200,
            n_twins=400,
            duration_months=12,
            random_seed=42,
        )

        temporal_payload = None
        store = Store(args.db)
        years = [r["year"] for r in store.conn.execute(
            "SELECT year FROM records WHERE year IS NOT NULL ORDER BY year")]
        if years:
            cutoff = args.cutoff if args.cutoff is not None else years[int(len(years) * 0.7)]
            temporal_payload = benchmark_to_dict(run_benchmark(store, cutoff))

        scorecard = {
            "splice_benchmark": spliceai_benchmark_to_dict(splice_res),
            "aso_tiling_benchmark": aso_res.to_dict(),
            "digital_twin_calibration": twin_cal,
            "temporal_holdout": temporal_payload,
        }
        if args.json:
            print(json.dumps(scorecard, indent=2, default=str))
            return 0

        header("SYNAPTIK Empirical Validation Scorecard")
        print(f"  1. Splice Benchmark (n={splice_res.n_cases}): "
              f"Sens={splice_res.sensitivity:.1%}, Spec={splice_res.specificity:.1%}, ROC-AUC={splice_res.auroc:.3f}")
        print(f"  2. ASO Tiling Benchmark (n={aso_res.n_loci} loci): "
              f"Hit@5(bio)={aso_res.hit_rate_at_5_biophysics:.1%}, MRR(bio)={aso_res.mrr_biophysics:.3f} "
              f"vs seq-only Hit@5={aso_res.hit_rate_at_5_sequence_only:.1%}, MRR={aso_res.mrr_sequence_only:.3f}")
        print(f"  3. Digital Twin Calibration (n=200 held-out): "
              f"P10-P90 empirical coverage={twin_cal['mean_empirical_coverage']:.1%} (nominal 80.0%)")
        if temporal_payload:
            syn_ranker = next((r for r in temporal_payload.get("rankers", []) if r["name"] == "synaptik"), None)
            auroc_str = f"{syn_ranker['auroc']:.3f}" if syn_ranker and syn_ranker.get("auroc") is not None else "undefined"
            print(f"  4. Literature Temporal Hold-Out (cutoff={temporal_payload['cutoff_year']}): "
                  f"synaptik AUROC={auroc_str}")
        else:
            print("  4. Literature Temporal Hold-Out: skipped (no dated records in DB)")
        return 0

    if getattr(args, "aso", False):
        from ...drug.benchmark_aso import run_aso_tiling_benchmark
        aso_res = run_aso_tiling_benchmark()
        if args.json:
            print(json.dumps(aso_res.to_dict(), indent=2))
            return 0
        header("Published ASO Tiling Benchmark (4 Neurogenetics Loci)")
        print(f"  Biophysics-coupled : Hit@5={aso_res.hit_rate_at_5_biophysics:.1%}, MRR={aso_res.mrr_biophysics:.3f}")
        print(f"  Sequence-only      : Hit@5={aso_res.hit_rate_at_5_sequence_only:.1%}, MRR={aso_res.mrr_sequence_only:.3f}")
        print()
        for o in aso_res.outcomes:
            print(f"  {o.gene:<8s} rank(bio)={str(o.best_rank_with_biophysics):<5s} "
                  f"rank(seq)={str(o.best_rank_sequence_only):<5s} "
                  f"top_score={o.top_candidate_score:<5.1f} dG_open={o.top_delta_g_open}")
        return 0

    if getattr(args, "spliceai", False):
        threshold = getattr(args, "threshold", 0.20) or 0.20
        res = run_splice_benchmark(threshold=threshold)
        if args.json:
            print(json.dumps(spliceai_benchmark_to_dict(res), indent=2, default=str))
            return 0

        header("Neuronal Splice Prediction Benchmark (Gold Standard)")
        print(Fmt.dim(
            "Evaluation against curated ALS/FTD and neurodegenerative splice variants\n"
            "STMN2, UNC13A, SMN2, ELP1, COL6A1, and matched neutral controls."
        ))
        print(f"  total variants evaluated : {res.n_cases} ({res.n_positives} positive, {res.n_negatives} negative)")
        print(f"  decision threshold (score): {res.threshold:.2f}")
        print()
        header("Validation Metrics")
        print(f"  Sensitivity (Recall)     : {res.sensitivity:6.1%}")
        print(f"  Specificity              : {res.specificity:6.1%}")
        print(f"  Precision (PPV)          : {res.precision:6.1%}")
        print(f"  Balanced Accuracy        : {res.accuracy:6.1%}")
        print(f"  ROC-AUC                  : {res.auroc:6.3f}")
        print()
        header("Variant-Level Predictions")
        print(f"  {'Gene':<8s} {'Variant':<22s} {'Mechanism':<35s} {'Score':>6s} {'Obs':>5s} {'Pred':>5s}")
        for d in res.details:
            obs = "POS" if d["is_splice_altering"] else "NEG"
            pred = "POS" if d["predicted_positive"] else "NEG"
            match = Fmt.green("✓") if obs == pred else Fmt.red("✗")
            mech = (d["mechanism"][:32] + "...") if len(d["mechanism"]) > 35 else d["mechanism"]
            print(f"  {d['gene']:<8s} {d['variant_id']:<22s} {mech:<35s} {d['score']:6.2f} {obs:>5s} {pred:>5s} {match}")
        return 0

    store = Store(args.db)
    cutoff = args.cutoff
    if cutoff is None:
        years = [r["year"] for r in store.conn.execute(
            "SELECT year FROM records WHERE year IS NOT NULL ORDER BY year")]
        if not years:
            print(Fmt.red("No dated records in the knowledge base."))
            return 1
        # Default: the year below which ~70% of records fall.
        cutoff = years[int(len(years) * 0.7)]

    res = run_benchmark(store, cutoff)
    if args.json:
        print(json.dumps(benchmark_to_dict(res), indent=2, default=str))
        return 0

    header(f"Retrospective benchmark, cutoff {res.cutoff_year}")
    print(Fmt.dim(
        "Associations computed from records up to the cutoff, scored against "
        "whether a NEW research group later reported the same thing."))
    print(f"  records before / after   {res.n_records_before} / "
          f"{res.n_records_after}")
    print(f"  phenotypes evaluated     {res.n_phenotypes}")
    print(f"  associations evaluated   {len(res.cases)} "
          f"({res.n_replicated} later replicated)")
    print()
    style = Fmt.yellow if res.underpowered else Fmt.bold
    print(style(res.verdict))

    if res.rankers:
        header("Predicting later independent replication (AUROC)")
        for r in res.rankers:
            val = "undefined" if r.auroc is None else f"{r.auroc:.3f}"
            mark = Fmt.bold if r.name == "synaptik" else str
            print(mark(f"  {r.name:14s} {val:>9s}") +
                  Fmt.dim(f"   {r.description}"))
        print(Fmt.dim("  0.5 = chance. The comparison that matters is "
                      "against record_count."))
    if res.tier_table:
        header("Replication rate by confidence tier")
        print(f"  {'tier':14s} {'n':>4s} {'replicated':>11s} "
              f"{'rate':>6s} {'contradicted':>13s}")
        for t in res.tier_table:
            print(f"  {t['tier']:14s} {t['n']:4d} {t['replicated']:11d} "
                  f"{t['replication_rate']:6.0%} {t['contradicted']:13d}")
    if res.caveats:
        print()
        for c in res.caveats:
            print(Fmt.yellow(f"  ! {c}"))
    return 0


def register(sub, common) -> None:
    s = sub.add_parser(
        "benchmark", parents=[common],
        help="test whether confidence tiers predict later independent "
             "replication (temporal hold-out), or validate sequence predictors")
    s.add_argument("--cutoff", type=int,
                   help="freeze the corpus at this year (default: the year "
                        "below which ~70%% of records fall)")
    s.add_argument("--spliceai", action="store_true",
                   help="run neuronal splice prediction benchmark against curated gold-standard variants")
    s.add_argument("--aso", action="store_true",
                   help="run ASO tiling benchmark against published neurogenetics loci")
    s.add_argument("--all", action="store_true",
                   help="run full empirical validation scorecard (splice, ASO, twin calibration, temporal hold-out)")
    s.add_argument("--threshold", type=float, default=0.20,
                   help="splice prediction score cutoff (default: 0.20)")
    s.add_argument("--json", action="store_true")
    s.set_defaults(func=cmd_benchmark)
