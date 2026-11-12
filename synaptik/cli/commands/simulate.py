"""`synaptik simulate`: causal genomic counterfactual simulation for therapeutic target discovery."""

from __future__ import annotations

import json
from dataclasses import asdict

from ...simulate.counterfactual import (
    CounterfactualSimulator,
    InterventionType,
    BENCHMARK_DISEASE_STATES,
)
from ..main import Fmt, header


def cmd_simulate(args) -> int:
    disease = getattr(args, "disease", "als") or "als"
    disease = disease.lower().strip()
    alpha = float(getattr(args, "alpha", 0.8) or 0.8)
    as_json = getattr(args, "json", False)

    simulator = CounterfactualSimulator()

    # Mode 1: High-throughput target screen
    if getattr(args, "screen", False):
        try:
            results = simulator.screen(disease=disease, efficacy_alpha=alpha)
        except Exception as e:
            if as_json:
                print(json.dumps({"error": str(e)}))
            else:
                print(Fmt.red(f"Error during screening: {e}"))
            return 1

        if as_json:
            out = {
                "disease": disease.upper(),
                "efficacy_alpha": alpha,
                "n_targets_screened": len(results),
                "ranked_targets": [
                    {
                        "target_gene": r.target_gene,
                        "intervention": r.intervention.value,
                        "composite_rescue_score": r.composite_rescue_score,
                        "rescue_efficiency_pct": r.rescue_efficiency_pct,
                        "cosine_rescue_alignment": r.cosine_rescue_alignment,
                        "collateral_toxicity_index": r.collateral_toxicity_index,
                        "is_clinically_viable": r.is_clinically_viable,
                        "verdict": r.verdict,
                    }
                    for r in results
                ],
            }
            print(json.dumps(out, indent=2))
            return 0

        header(f"Causal Target Screen: {disease.upper()} Transcriptomic Rescue (alpha={alpha:.2f})")
        print(f"Evaluated {len(results)} perturbation targets against {disease.upper()} baseline vector.\n")

        print(f"{'Target':<10s} {'Intervention':<16s} {'Score':<8s} {'Rescue %':<10s} {'Cosine':<8s} {'Toxicity':<10s} {'Viability':<12s}")
        print("-" * 76)

        for r in results:
            if r.composite_rescue_score >= 70.0:
                score_str = Fmt.green(f"{r.composite_rescue_score:>5.1f}")
            elif r.composite_rescue_score >= 45.0:
                score_str = Fmt.yellow(f"{r.composite_rescue_score:>5.1f}")
            else:
                score_str = Fmt.red(f"{r.composite_rescue_score:>5.1f}")

            viab_str = Fmt.green("VIABLE") if r.is_clinically_viable else Fmt.red("NON-VIABLE")
            int_str = r.intervention.value
            print(f"{r.target_gene:<10s} {int_str:<16s} {score_str:<17s} {r.rescue_efficiency_pct:>6.1f}%   {r.cosine_rescue_alignment:>6.2f}   {r.collateral_toxicity_index:>6.2f}     {viab_str}")

        print("\nTop Lead Recommendations:")
        viable_leads = [r for r in results if r.is_clinically_viable]
        if viable_leads:
            top = viable_leads[0]
            print(f"  * Highest-ranking therapeutic intervention: {top.target_gene} ({top.intervention.value})")
            print(f"    Rescue Score: {top.composite_rescue_score}/100 | Cosine: {top.cosine_rescue_alignment:.2f} | Collateral Risk: {top.collateral_toxicity_index:.2f}")
            for rec in top.recommendations:
                print(f"    -> {rec}")
        else:
            print("  ! No viable candidates passed toxicity and efficiency filters.")
        return 0

    # Mode 2: Targeted single perturbation simulation
    target = getattr(args, "target", None)
    if not target:
        if as_json:
            print(json.dumps({"error": "Target gene required. Specify --target <GENE> or --screen."}))
        else:
            print(Fmt.red("Error: Target gene required. Specify --target <GENE> or use --screen to rank all candidates."))
        return 1

    intervention_str = getattr(args, "intervention", "knockdown") or "knockdown"
    try:
        intervention = InterventionType(intervention_str.lower().strip())
    except ValueError:
        intervention = InterventionType.KNOCKDOWN

    try:
        res = simulator.simulate(
            disease=disease,
            target_gene=target,
            intervention=intervention,
            efficacy_alpha=alpha,
        )
    except Exception as e:
        if as_json:
            print(json.dumps({"error": str(e)}))
        else:
            print(Fmt.red(f"Error during simulation: {e}"))
        return 1

    if as_json:
        out = {
            "disease": res.disease,
            "target_gene": res.target_gene,
            "intervention": res.intervention.value,
            "efficacy_alpha": res.efficacy_alpha,
            "disease_magnitude": res.disease_magnitude,
            "post_intervention_magnitude": res.post_intervention_magnitude,
            "rescue_efficiency_pct": res.rescue_efficiency_pct,
            "cosine_rescue_alignment": res.cosine_rescue_alignment,
            "collateral_toxicity_index": res.collateral_toxicity_index,
            "composite_rescue_score": res.composite_rescue_score,
            "is_clinically_viable": res.is_clinically_viable,
            "verdict": res.verdict,
            "key_gene_trajectories": res.key_gene_trajectories,
            "pathway_impacts": res.pathway_impacts,
            "recommendations": res.recommendations,
        }
        print(json.dumps(out, indent=2))
        return 0

    header(f"Causal Counterfactual Simulation: {res.target_gene} in {res.disease}")
    print(f"Perturbation Type            : {res.intervention.value}")
    print(f"Perturbation Efficacy (alpha): {res.efficacy_alpha:.2f}")
    print(f"Baseline Disease Magnitude   : {res.disease_magnitude:.3f}")
    print(f"Post-Intervention Magnitude  : {res.post_intervention_magnitude:.3f}")
    print(f"Vector Rescue Efficiency     : {res.rescue_efficiency_pct:.1f}%")
    print(f"Cosine Directional Alignment : {res.cosine_rescue_alignment:.3f} (-1 to +1)")
    print(f"Collateral Toxicity Index    : {res.collateral_toxicity_index:.2f} (0=clean, 1=lethal)")
    print(f"Composite Rescue Score       : {res.composite_rescue_score:.1f} / 100")
    print(f"Clinical Viability Status    : {'VIABLE' if res.is_clinically_viable else 'NON-VIABLE'}")

    print("\n" + "=" * 78)
    print("Translational Verdict")
    print("=" * 78)
    if res.is_clinically_viable and res.composite_rescue_score >= 70:
        print(Fmt.green(res.verdict))
    elif res.is_clinically_viable:
        print(Fmt.yellow(res.verdict))
    else:
        print(Fmt.red(res.verdict))

    print("\nTranscriptomic Trajectory (Key Dysregulated Genes)")
    print("=" * 78)
    print(f"{'Gene':<12s} {'Baseline log2FC':<18s} {'Shift (Delta)':<18s} {'Post-Intervention':<18s}")
    print("-" * 78)
    for traj in res.key_gene_trajectories[:12]:
        g = traj["gene"]
        base = traj["baseline_log2fc"]
        delta = traj["delta_log2fc"]
        post = traj["post_log2fc"]
        print(f"{g:<12s} {base:>+10.2f}         {delta:>+10.2f}         {post:>+10.2f}")

    if res.recommendations:
        print("\nTranslational Recommendations")
        print("=" * 78)
        for rec in res.recommendations:
            print(f"  -> {rec}")

    return 0


def register(sub, common) -> None:
    s = sub.add_parser(
        "simulate", parents=[common],
        help="causal genomic counterfactual simulation for target discovery and phenotype rescue",
    )
    s.add_argument("--disease", choices=["als", "ftd", "ad"], default="als",
                   help="target disease state vector (default: als)")
    s.add_argument("--target", metavar="GENE",
                   help="perturbation target gene symbol (e.g. ATXN2, STMN2, UNC13A, PIKFYVE)")
    s.add_argument("--intervention", choices=["knockdown", "knockout", "overexpression", "aso_block"],
                   default="knockdown", help="perturbation modality (default: knockdown)")
    s.add_argument("--alpha", type=float, default=0.8,
                   help="perturbation efficacy scaling factor between 0.0 and 1.0 (default: 0.8)")
    s.add_argument("--screen", action="store_true",
                   help="screen and rank all cataloged perturbation targets against disease state")
    s.add_argument("--json", action="store_true", help="output results as JSON")
    s.set_defaults(func=cmd_simulate)
