"""`synaptik design`: a design brief for a planned experiment."""

from __future__ import annotations

import json

from ...schema import (
    AssayClass, CellModelClass, ExperimentQuery, PerturbationType,
)
from ...store import Store
from ...workflows import design_brief, design_brief_to_dict
from ..main import Fmt, header


def _fmt_counts(pairs) -> str:
    return ", ".join(f"{k} ({v})" for k, v in pairs)


def cmd_design(args) -> int:
    q = ExperimentQuery(
        cell_model_class=CellModelClass(args.model) if args.model else None,
        days_in_vitro=args.div,
        perturbation_type=PerturbationType(args.perturbation)
        if args.perturbation else None,
        perturbation_target=args.target,
        variant_hgvs=args.variant,
        assay_class=AssayClass(args.assay) if args.assay else None,
        markers=args.markers or [],
        phenotype=args.phenotype,
    )
    if not any([q.cell_model_class, q.perturbation_target, q.assay_class,
                q.phenotype, q.variant_hgvs]):
        print(Fmt.red("Describe the experiment: give at least --target, "
                      "--assay or --phenotype."))
        return 1
    if args.cv is not None and not (0 < args.cv < 2):
        print(Fmt.red("--cv is a fraction, e.g. 0.25 for 25%."))
        return 1
    rho = args.rho if args.rho is not None else 0.7
    if not (0 <= rho < 1):
        print(Fmt.red("--rho must be between 0 and 0.99 (correlation between paired batches)."))
        return 1

    b = design_brief(Store(args.db), q, cv=args.cv, rho=rho, power=args.power)
    if args.json:
        print(json.dumps(design_brief_to_dict(b), indent=2, default=str))
        return 0

    header("Experiment design brief")
    print(b.basis)
    print(Fmt.bold("What your experiment would add: ") + b.value_statement)

    if b.n_exact + b.n_close or b.div_median is not None:
        header("How it has been done")
        if b.div_median is not None:
            print(f"  maturity          DIV {b.div_median} "
                  f"(range {b.div_range[0]}-{b.div_range[1]})")
        if b.differentiation_protocols:
            print(f"  differentiation   {_fmt_counts(b.differentiation_protocols)}")
        if b.cocultures:
            print(f"  co-culture        {_fmt_counts(b.cocultures)}")
        if b.markers:
            print(f"  markers           {_fmt_counts(b.markers)}")
        if b.n_replicates_median is not None:
            print(f"  replicates        median {b.n_replicates_median} "
                  f"(range {b.n_replicates_range[0]}-"
                  f"{b.n_replicates_range[1]})")
        if b.isogenic_fraction is not None:
            print(f"  isogenic controls {b.isogenic_fraction:.0%} of prior "
                  "experiments")
        if b.statistical_tests:
            print(f"  statistics        {_fmt_counts(b.statistical_tests)}")

    if b.direction_counts:
        header("What was found")
        counts = ", ".join(
            f"{v} {k.replace('_', ' ')}" for k, v in
            sorted(b.direction_counts.items(), key=lambda kv: -kv[1]))
        print(f"  direction         {b.consensus_direction}  ({counts})")
        if b.median_effect_percent is not None:
            print(f"  effect size       median {b.median_effect_percent:.0f}% "
                  f"({len(b.effect_sizes_percent)} quantified report(s))")
        print(f"  null results      {b.n_null_results}")

    if b.power:
        header("Sample size (independent differentiations per group / isogenic pairs)")
        print(Fmt.dim(f"  two-sample & paired t-tests, {args.power:.0%} power, "
                      "two-sided alpha 0.05"))
        last = None
        for p in b.power:
            if p.label != last:
                print(Fmt.bold(f"  {p.label}"))
                last = p.label
            n_unp = str(p.n_per_group) if p.n_per_group else ">200"
            n_pr = str(p.n_paired) if p.n_paired else ">100"
            print(f"    {abs(p.effect_percent):5.0f}% effect   CV "
                  f"{p.cv:.0%}   unpaired n = {n_unp:<5s} | paired (rho={p.paired_rho:.2g}) = {n_pr} pairs")
        print()
        if b.impractical:
            n = b.recommended_n if b.recommended_n else ">200"
            print(Fmt.yellow(Fmt.bold(
                f"  A realistic effect needs {n} per group unpaired: not practical.")))
            if b.recommended_n_paired:
                print(Fmt.green(Fmt.bold(
                    f"  PAIRING RESCUE: Differentiating isogenic pairs in parallel (rho={b.paired_rho:.2g}) needs only {b.recommended_n_paired} pairs!")))
        elif b.recommended_n:
            print(Fmt.green(Fmt.bold(
                f"  Plan for {b.recommended_n} per group unpaired" +
                (f" (or {b.recommended_n_paired} isogenic pairs)." if b.recommended_n_paired else "."))))
        print(Fmt.dim(f"  {b.recommended_n_basis}"))
        if b.detectable:
            print()
            print(Fmt.bold("  What a feasible design can detect"))
            for n, cv, eff in b.detectable:
                print(f"    n = {n} per group   CV {cv:.0%}   "
                      f"effects of {eff:.0f}% or larger")
        print(Fmt.dim("  Cells, fields and wells from one differentiation "
                      "are technical replicates and do not count toward n."))

    if b.weaknesses:
        header("Weaknesses in the existing evidence")
        for w in b.weaknesses:
            print(Fmt.yellow(f"  - {w}"))
    if b.recommendations:
        header("Recommendations")
        for r in b.recommendations:
            print(Fmt.blue(f"  -> {r}"))
    if b.pitfalls:
        header("Pitfalls for this readout")
        for p in b.pitfalls:
            print(f"  - {p}")
    if b.key_papers:
        header("Read first")
        for k in b.key_papers:
            print(Fmt.dim(f"  {k}"))
    return 0


def register(sub, common) -> None:
    s = sub.add_parser(
        "design", parents=[common],
        help="design brief for a planned experiment: protocol, sample size, "
             "pitfalls")
    s.add_argument("--model", choices=[c.value for c in CellModelClass])
    s.add_argument("--div", type=int, help="days in vitro")
    s.add_argument("--perturbation",
                   choices=[c.value for c in PerturbationType])
    s.add_argument("--target", help="gene or protein perturbed")
    s.add_argument("--variant", help="HGVS variant string")
    s.add_argument("--assay", choices=[c.value for c in AssayClass])
    s.add_argument("--markers", nargs="+")
    s.add_argument("--phenotype")
    s.add_argument("--cv", type=float,
                   help="your between-differentiation coefficient of "
                        "variation as a fraction (e.g. 0.25)")
    s.add_argument("--rho", type=float, default=0.7,
                   help="expected inter-batch correlation in a paired isogenic "
                        "design (default 0.7, typical 0.6-0.85)")
    s.add_argument("--power", type=float, default=0.8)
    s.add_argument("--json", action="store_true")
    s.set_defaults(func=cmd_design)
