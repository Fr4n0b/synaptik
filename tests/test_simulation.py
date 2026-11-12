"""Tests for the Causal Genomic Counterfactual Simulator (`synaptik simulate`)."""

import unittest
from synaptik.simulate.counterfactual import (
    CounterfactualSimulator,
    InterventionType,
    BENCHMARK_DISEASE_STATES,
)
from synaptik.workflows import simulate_intervention, screen_counterfactual_targets


def check(desc, cond):
    if cond:
        print(f"  PASS  {desc}")
    else:
        print(f"  FAIL  {desc}")
        raise AssertionError(desc)


def run_tests():
    print("=" * 70)
    print("SYNAPTIK genomic counterfactual simulation test suite")
    print("=" * 70)

    sim = CounterfactualSimulator()

    # ------------------------------------------------------------------
    # 1. Vector mathematics and scaling
    # ------------------------------------------------------------------
    print("\n[simulate: vector mathematics & scaling]")
    res_full = sim.simulate("als", "ATXN2", InterventionType.KNOCKDOWN, efficacy_alpha=1.0)
    res_half = sim.simulate("als", "ATXN2", InterventionType.KNOCKDOWN, efficacy_alpha=0.5)

    check("disease baseline magnitude is positive", res_full.disease_magnitude > 4.0)
    check("cosine rescue alignment is strongly positive for ATXN2 knockdown", res_full.cosine_rescue_alignment > 0.40)
    check("half efficacy produces smaller delta than full efficacy",
          res_half.post_intervention_magnitude > res_full.post_intervention_magnitude)
    check("rescue efficiency pct is higher with full efficacy",
          res_full.rescue_efficiency_pct > res_half.rescue_efficiency_pct)

    # ------------------------------------------------------------------
    # 2. Curated ALS therapeutic rescues
    # ------------------------------------------------------------------
    print("\n[simulate: ALS therapeutic rescues]")
    # ATXN2 Knockdown
    res_atxn2 = sim.simulate("als", "ATXN2", InterventionType.KNOCKDOWN, efficacy_alpha=0.8)
    check("ATXN2 knockdown predicted clinically viable", res_atxn2.is_clinically_viable)
    check("ATXN2 composite score >= 60 (got " + str(res_atxn2.composite_rescue_score) + ")",
          res_atxn2.composite_rescue_score >= 60.0)
    check("ATXN2 collateral toxicity index is low (< 0.20)", res_atxn2.collateral_toxicity_index < 0.20)
    check("verdict reflects strong rescue", "STRONG RESCUE" in res_atxn2.verdict)

    # STMN2 ASO rescue
    res_stmn2 = sim.simulate("als", "STMN2", InterventionType.ASO_STERIC_BLOCK, efficacy_alpha=0.9)
    check("STMN2 ASO block predicted clinically viable", res_stmn2.is_clinically_viable)
    stmn2_traj = next(t for t in res_stmn2.key_gene_trajectories if t["gene"] == "STMN2")
    check("STMN2 baseline was severely negative in ALS (-2.2)", stmn2_traj["baseline_log2fc"] == -2.2)
    check("STMN2 delta shifts strongly positive (+1.89)", stmn2_traj["delta_log2fc"] > 1.8)
    check("STMN2 post-intervention state is rescued close to 0", abs(stmn2_traj["post_log2fc"]) < 0.5)

    # UNC13A ASO rescue
    res_unc13a = sim.simulate("als", "UNC13A", InterventionType.ASO_STERIC_BLOCK, efficacy_alpha=0.85)
    check("UNC13A ASO block predicted clinically viable", res_unc13a.is_clinically_viable)
    check("UNC13A composite score >= 55 (got " + str(res_unc13a.composite_rescue_score) + ")",
          res_unc13a.composite_rescue_score >= 55.0)

    # ------------------------------------------------------------------
    # 3. Lethal & toxic counterfactuals
    # ------------------------------------------------------------------
    print("\n[simulate: lethal & toxic negative controls]")
    # POLR2A Knockdown (Pan-essential lethal control)
    res_polr2a = sim.simulate("als", "POLR2A", InterventionType.KNOCKDOWN, efficacy_alpha=0.8)
    check("POLR2A knockdown flagged as non-viable", not res_polr2a.is_clinically_viable)
    check("POLR2A collateral toxicity index is critical (> 0.90)", res_polr2a.collateral_toxicity_index > 0.90)
    check("POLR2A verdict reports viability failure", "VIABILITY FAILURE" in res_polr2a.verdict)

    # TARDBP Overexpression (Known toxic aggregation phenotype)
    res_tardbp = sim.simulate("als", "TARDBP", InterventionType.OVEREXPRESSION, efficacy_alpha=0.8)
    check("TARDBP overexpression flagged as non-viable", not res_tardbp.is_clinically_viable)
    check("TARDBP overexpression collateral risk is elevated (> 0.50)", res_tardbp.collateral_toxicity_index >= 0.50)

    # ------------------------------------------------------------------
    # 4. High-throughput screening & ranking
    # ------------------------------------------------------------------
    print("\n[simulate: high-throughput counterfactual target screening]")
    screen_results = sim.screen("als", efficacy_alpha=0.8)
    check("screen evaluates multiple cataloged targets", len(screen_results) >= 6)
    check("screen results sorted descending by composite rescue score",
          screen_results[0].composite_rescue_score >= screen_results[1].composite_rescue_score)

    top_target = screen_results[0]
    check("top target is clinically viable", top_target.is_clinically_viable)
    check("top target is an established ALS rescue candidate (STMN2 or ATXN2)",
          top_target.target_gene in ("STMN2", "ATXN2", "UNC13A"))

    bottom_target = screen_results[-1]
    check("bottom target is non-viable lethal/toxic control", not bottom_target.is_clinically_viable)
    check("bottom target is POLR2A or TARDBP", bottom_target.target_gene in ("POLR2A", "TARDBP"))

    # ------------------------------------------------------------------
    # 5. Workflows and API integration
    # ------------------------------------------------------------------
    print("\n[simulate: workflows API serialization]")
    api_single = simulate_intervention("als", "ATXN2", intervention="knockdown", efficacy_alpha=0.8)
    check("workflow simulate_intervention returns target_gene", api_single["target_gene"] == "ATXN2")
    check("workflow simulate_intervention returns is_clinically_viable", api_single["is_clinically_viable"] is True)
    check("workflow simulate_intervention contains key_gene_trajectories", len(api_single["key_gene_trajectories"]) > 0)

    api_screen = screen_counterfactual_targets("als", efficacy_alpha=0.8)
    check("workflow screen_counterfactual_targets returns list", isinstance(api_screen, list))
    check("workflow screen returns sorted entries",
          api_screen[0]["composite_rescue_score"] >= api_screen[-1]["composite_rescue_score"])

    print("\n" + "=" * 70)
    print("All checks passed.")


if __name__ == "__main__":
    run_tests()
