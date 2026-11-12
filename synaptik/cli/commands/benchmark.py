"""`synaptik benchmark`: do the confidence tiers predict later replication?"""

from __future__ import annotations

import json

from ...store import Store
from ...workflows import benchmark_to_dict, run_benchmark
from ..main import Fmt, header


def cmd_benchmark(args) -> int:
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
             "replication (temporal hold-out)")
    s.add_argument("--cutoff", type=int,
                   help="freeze the corpus at this year (default: the year "
                        "below which ~70%% of records fall)")
    s.add_argument("--json", action="store_true")
    s.set_defaults(func=cmd_benchmark)
