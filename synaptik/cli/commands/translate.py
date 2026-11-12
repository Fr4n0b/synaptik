"""`synaptik translate`: bridge in vitro iPSC data to human post-mortem single-cell atlases."""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

from ...labdata import parse_splicing_table
from ...schema import Direction
from ...translate.atlas import evaluate_translational_concordance
from ..main import Fmt, header


def cmd_translate(args) -> int:
    disease = getattr(args, "disease", "als") or "als"
    entries: list[tuple[str, Direction | None, float | None]] = []

    if getattr(args, "splicing", None):
        source = args.splicing
        perturbed = getattr(args, "perturbed", "sample1") or "sample1"
        try:
            table = parse_splicing_table(source, perturbed=perturbed)
            for ev in table.events:
                direct = Direction.INCREASE if ev.delta_psi > 0 else Direction.DECREASE
                entries.append((ev.gene, direct, ev.delta_psi))
        except Exception as e:
            print(Fmt.red(f"Error reading splicing table: {e}"))
            return 1
    elif getattr(args, "genes", None):
        for g in args.genes:
            entries.append((g, None, None))
    else:
        print(Fmt.red("Error: Supply either --splicing <file> or --genes <gene1 gene2 ...>."))
        return 1

    if not entries:
        print(Fmt.yellow("No genes found to evaluate."))
        return 0

    res = evaluate_translational_concordance(entries, disease=disease)

    if getattr(args, "json", False):
        data = {
            "disease": res.disease.value,
            "n_genes_evaluated": res.n_genes_evaluated,
            "n_concordant": res.n_concordant,
            "n_discordant": res.n_discordant,
            "n_developmental_gap": res.n_developmental_gap,
            "n_cell_type_mismatch": res.n_cell_type_mismatch,
            "n_uncharacterized": res.n_uncharacterized,
            "concordance_score_pct": res.concordance_score_pct,
            "translational_verdict": res.translational_verdict,
            "enriched_cell_types": [
                {"cell_type": ct, "count": cnt, "percent": pct}
                for ct, cnt, pct in res.enriched_cell_types
            ],
            "recommendations": res.recommendations,
            "gene_calls": [asdict(c) for c in res.gene_calls],
        }
        print(json.dumps(data, indent=2, default=str))
        return 0

    header(f"Dish-to-Brain Translational Concordance: {res.disease.value}")
    print(Fmt.dim(
        "Cross-referencing in vitro findings against human post-mortem single-nucleus "
        "RNA-seq (NYGC ALS, SEA-AD, Allen Brain Cell Atlas)."
    ))
    print(f"  Genes / Events Evaluated    : {res.n_genes_evaluated}")
    print(f"  Post-Mortem Concordant      : {Fmt.green(str(res.n_concordant))}")
    print(f"  Post-Mortem Discordant      : {Fmt.red(str(res.n_discordant)) if res.n_discordant > 0 else '0'}")
    print(f"  Developmental Gaps          : {Fmt.yellow(str(res.n_developmental_gap))}")
    print(f"  Uncharacterized In Vivo     : {res.n_uncharacterized}")
    print(f"  Translational Concordance   : {Fmt.bold(f'{res.concordance_score_pct:.1f}%')}")
    print()

    header("Translational Verdict")
    style = Fmt.green if res.concordance_score_pct >= 75 else (Fmt.yellow if res.concordance_score_pct >= 50 else Fmt.red)
    print(style(res.translational_verdict))
    print()

    if res.enriched_cell_types:
        header("Human In Vivo Cell-Type Specificity")
        for ct, cnt, pct in res.enriched_cell_types[:4]:
            print(f"  {ct:<42s} : {cnt} marker(s) ({pct:.1f}%)")
        print()

    header("Gene-Level Translational Mapping")
    print(f"  {'Gene':<10s} {'Status':<16s} {'Cell Type':<32s} {'In Vitro':>9s} {'Post-Mortem':>12s} {'Evidence'}")
    for c in res.gene_calls:
        if c.status.value == "concordant":
            tag = Fmt.green("CONCORDANT")
        elif c.status.value == "discordant":
            tag = Fmt.red("DISCORDANT")
        elif c.status.value == "developmental_gap":
            tag = Fmt.yellow("DEV_GAP")
        else:
            tag = Fmt.dim("UNMAPPED")

        ev_short = c.evidence_source.split(";")[0][:30]
        print(f"  {c.gene:<10s} {tag:<25s} {c.primary_cell_type[:30]:<32s} {c.in_vitro_direction:>9s} {c.post_mortem_direction:>12s}  {ev_short}")

    if res.recommendations:
        print()
        header("Translational Recommendations")
        for r in res.recommendations:
            print(f"  -> {r}")

    print("\n" + "-" * 78)
    print("Cite: Francesco Nobile (ORCID: 0009-0001-9477-3332); run 'synaptik cite'")
    return 0


def register(sub, common) -> None:
    s = sub.add_parser(
        "translate", parents=[common],
        help="evaluate in vitro iPSC findings against human post-mortem single-cell brain atlases",
    )
    s.add_argument("--splicing", metavar="PATH",
                   help="path to an rMATS, MAJIQ, or LeafCutter differential splicing table")
    s.add_argument("--genes", nargs="+", metavar="GENE",
                   help="list of gene symbols to evaluate")
    s.add_argument("--disease", choices=["als", "ftd", "ad"], default="als",
                   help="target human post-mortem disease atlas (default: als)")
    s.add_argument("--perturbed", default="sample1",
                   help="which sample/group is perturbed (for splicing tables: sample1/sample2/positive/negative)")
    s.add_argument("--json", action="store_true", help="output results as JSON")
    s.set_defaults(func=cmd_translate)
