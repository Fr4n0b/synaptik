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
    s.add_argument("--threshold", type=float, default=0.20,
                   help="splice prediction score cutoff (default: 0.20)")
    s.add_argument("--json", action="store_true")
    s.set_defaults(func=cmd_benchmark)
