"""`synaptik deposit`: reconcile your own results against the literature."""

from __future__ import annotations

import json

from ...labdata import LabCall, LabDataError, parse_splicing_table
from ...store import Store
from ...workflows import (
    deposit_phenotypes, deposit_splicing, phenotype_report_to_dict,
    splicing_report_to_dict,
)
from ..main import Fmt, header

CALL_STYLE = {
    "novel": Fmt.green, "discordant": Fmt.red, "consistent": Fmt.cyan,
    "same_gene_other_event": Fmt.yellow, "contested": Fmt.yellow,
}


def cmd_deposit_splicing(args) -> int:
    try:
        table = parse_splicing_table(
            args.file, perturbed=args.perturbed, event_type=args.event_type)
    except LabDataError as e:
        print(Fmt.red(f"Cannot use this table: {e}"))
        return 1

    rep = deposit_splicing(Store(args.db), table, args.fdr, args.min_dpsi)
    if args.json:
        print(json.dumps(splicing_report_to_dict(rep), indent=2, default=str))
        return 0

    header("Your splicing results vs the literature")
    print(Fmt.dim(
        f"{table.source_format} table, {rep.n_events} events in "
        f"{rep.n_genes_tested} genes; thresholds FDR <= {args.fdr}, "
        f"|dPSI| >= {args.min_dpsi}; corpus has {rep.corpus_events} published "
        f"events in {rep.corpus_genes} genes"))
    print(Fmt.dim("Your data is analysed in memory and is not stored."))
    print()
    print(rep.summary)
    for w in rep.warnings:
        print(Fmt.yellow(f"  ! {w}"))

    sections = [
        (LabCall.DISCORDANT, "CONTRADICT the literature - check these first"),
        (LabCall.NOVEL, "NOVEL - no published event in this gene"),
        (LabCall.SAME_GENE, "Known gene, different event"),
        (LabCall.CONSISTENT, "CONSISTENT with published reports"),
    ]
    for call, title in sections:
        group = rep.by_call(call)
        if not group:
            continue
        header(f"{title}  ({len(group)})")
        for r in group[: args.limit]:
            e = r.event
            style = CALL_STYLE.get(call.value, str)
            fdr = f"FDR {e.fdr:.2g}" if e.fdr is not None else "FDR n/a"
            print(style(Fmt.bold(e.gene)) +
                  f"  {e.event_type.value}  dPSI {e.delta_psi:+.2f}  {fdr}")
            print(f"  {r.note}")
            if r.n_papers:
                print(Fmt.dim(f"  literature: {r.n_papers} paper(s), "
                              f"{r.n_labs} group(s); "
                              f"{', '.join(r.source_ids[:4])}"))
            if r.regulators_reported:
                print(Fmt.dim("  published regulator(s): "
                              + ", ".join(r.regulators_reported)))
            if r.nmd_risk == "high":
                print(Fmt.red(f"  frame: {r.nmd_note}"))
            if r.aso_feasible and call in (LabCall.NOVEL, LabCall.DISCORDANT):
                print(Fmt.blue(f"  ASO: {r.aso_approach[:110]}"))
        if len(group) > args.limit:
            print(Fmt.dim(f"  ... and {len(group) - args.limit} more "
                          "(use --limit or --json)"))

    if rep.regulators:
        header("Regulators whose published targets you hit")
        print(Fmt.dim("Fisher test against the genes tested in YOUR table."))
        for g in rep.regulators[:8]:
            flag = Fmt.green if g.q_value < 0.1 else Fmt.dim
            agreement = ("n/a" if g.direction_agreement is None
                         else f"{g.direction_agreement:.0%}")
            print(flag(Fmt.bold(g.rbp)) +
                  f"  {g.n_targets_changed}/{g.n_targets_tested} tested "
                  f"targets changed  OR {g.odds_ratio}  q {g.q_value:.2g}  "
                  f"direction agreement {agreement}")
            print(Fmt.dim(f"  {', '.join(g.overlap_genes)}"))
            for c in g.caveats():
                print(Fmt.yellow(f"  ! {c}"))

    if rep.not_reproduced:
        header(f"Published events you tested but did not see  "
               f"({len(rep.not_reproduced)})")
        print(Fmt.dim("Genes present in your table with no significant "
                      "change. Genes absent from your table are not listed: "
                      "absence is not a negative result."))
        for n in rep.not_reproduced[: args.limit]:
            dpsi = ("n/a" if n.user_best_dpsi is None
                    else f"{n.user_best_dpsi:+.2f}")
            print(Fmt.bold(n.gene) +
                  f"  published: {n.literature_direction} "
                  f"({n.n_papers} paper(s), {n.n_labs} group(s))"
                  f"  your largest dPSI: {dpsi}")
    return 0


def cmd_deposit_phenotypes(args) -> int:
    try:
        rows, warnings = deposit_phenotypes(Store(args.db), args.file)
    except LabDataError as e:
        print(Fmt.red(f"Cannot use this table: {e}"))
        return 1
    if args.json:
        print(json.dumps(phenotype_report_to_dict(rows, warnings),
                         indent=2, default=str))
        return 0

    header("Your phenotypes vs the literature")
    print(Fmt.dim("Your data is analysed in memory and is not stored."))
    for w in warnings:
        print(Fmt.yellow(f"  ! {w}"))
    print()
    for r in rows:
        ob = r.observation
        style = CALL_STYLE.get(r.call, str)
        print(style(Fmt.bold(f"[{r.call}]")) +
              f"  {ob.phenotype}  ({ob.direction.value}"
              + (f", {ob.gene}" if ob.gene else "") + ")")
        if r.resolved_label:
            print(Fmt.dim(f"  resolved to: {r.resolved_label}"))
        print(f"  {r.note}")
        if r.source_ids:
            print(Fmt.dim(f"  {', '.join(r.source_ids[:5])}"))
        print()
    return 0


def register(sub, common) -> None:
    p = sub.add_parser(
        "deposit", parents=[common],
        help="reconcile your own (unpublished) results against the literature")
    kinds = p.add_subparsers(dest="kind", required=True)

    s = kinds.add_parser(
        "splicing", parents=[common],
        help="an rMATS table, or a generic gene/delta_psi/fdr table")
    s.add_argument("file")
    s.add_argument(
        "--perturbed", required=True,
        choices=["sample1", "sample2", "positive", "negative"],
        help="which side is the perturbed condition. rMATS: sample1 or "
             "sample2 (IncLevelDifference = sample1 - sample2). Generic "
             "tables: 'positive' if positive delta_psi means more inclusion "
             "in the perturbed condition")
    s.add_argument("--event-type", choices=["SE", "RI", "A5SS", "A3SS", "MXE"],
                   help="inferred from an rMATS file name when omitted")
    s.add_argument("--fdr", type=float, default=0.05)
    s.add_argument("--min-dpsi", type=float, default=0.1)
    s.add_argument("--limit", type=int, default=15)
    s.add_argument("--json", action="store_true")
    s.set_defaults(func=cmd_deposit_splicing)

    s = kinds.add_parser(
        "phenotypes", parents=[common],
        help="a CSV with phenotype and direction columns")
    s.add_argument("file")
    s.add_argument("--json", action="store_true")
    s.set_defaults(func=cmd_deposit_phenotypes)
