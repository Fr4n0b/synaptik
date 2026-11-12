"""Tests for Preclinical Study Design Registry & Multi-Site Meta-Analysis (`synaptik protocol`)."""

import unittest
from synaptik.protocol.registry import (
    create_preclinical_protocol,
    run_random_effects_meta_analysis,
    SiteReplicationResult,
    generate_96well_randomized_layout,
)
from synaptik.workflows import generate_study_protocol, synthesize_multisite_replication


def check(desc, cond):
    if cond:
        print(f"  PASS  {desc}")
    else:
        print(f"  FAIL  {desc}")
        raise AssertionError(desc)


def run_tests():
    print("=" * 70)
    print("SYNAPTIK preclinical protocol registry & meta-analysis test suite")
    print("=" * 70)

    # ------------------------------------------------------------------
    # 1. Study design protocol generation & cryptographic hashing
    # ------------------------------------------------------------------
    print("\n[protocol: pre-registration design & cryptographic hashing]")
    proto1 = create_preclinical_protocol("TARDBP", disease="als", paired_n=8, rho=0.75)
    proto2 = create_preclinical_protocol("TARDBP", disease="als", paired_n=8, rho=0.75)
    proto3 = create_preclinical_protocol("STMN2", disease="als", paired_n=8, rho=0.75)

    check("SHA-256 fingerprint is 64 characters hex", len(proto1.protocol_fingerprint_sha256) == 64)
    check("protocol fingerprint is strictly deterministic",
          proto1.protocol_fingerprint_sha256 == proto2.protocol_fingerprint_sha256)
    check("different target gene produces different fingerprint",
          proto1.protocol_fingerprint_sha256 != proto3.protocol_fingerprint_sha256)
    check("blinding verification token generated", proto1.blinding_verification_token.startswith("SYN-"))
    check("all 7 methodological rigor checklist items verified",
          len(proto1.rigor_checklist) == 7 and all(proto1.rigor_checklist.values()))

    # Plate layout checks
    layout = proto1.sample_plate_layout
    check("sample plate layout maps 16 wells for 8 pairs", len(layout) == 16)
    # Edge effect check: no A or H rows, no col 1 or 12
    wells = list(layout.keys())
    has_edge = any(w[0] in ("A", "H") or int(w[1:]) in (1, 12) for w in wells)
    check("perimeter edge wells reserved as evaporation barrier", not has_edge)
    # Pairing balance: 8 mutants and 8 controls
    n_mut = sum(1 for v in layout.values() if "Mutant" in v)
    n_ctrl = sum(1 for v in layout.values() if "IsogenicControl" in v)
    check("plate layout has exact balanced allocation (8 mutants, 8 controls)", n_mut == 8 and n_ctrl == 8)

    # ------------------------------------------------------------------
    # 2. Multi-site random-effects meta-analysis
    # ------------------------------------------------------------------
    print("\n[protocol: multi-site meta-analysis synthesis]")
    # Concordant replication dataset (low heterogeneity)
    concordant_sites = [
        SiteReplicationResult(site_id="Site_MGH", n_pairs=10, mean_difference=-28.5, standard_error=3.2),
        SiteReplicationResult(site_id="Site_Stanford", n_pairs=8, mean_difference=-31.0, standard_error=4.1),
        SiteReplicationResult(site_id="Site_UCL", n_pairs=12, mean_difference=-27.0, standard_error=2.8),
        SiteReplicationResult(site_id="Site_Karolinska", n_pairs=8, mean_difference=-29.4, standard_error=3.8),
    ]

    meta_concordant = run_random_effects_meta_analysis(concordant_sites, target_gene="TARDBP")
    check("evaluated 4 replication sites", meta_concordant.k_sites == 4)
    check("total paired replicates is 38", meta_concordant.total_n_pairs == 38)
    check("pooled effect is significantly negative (~ -28.6)", -32.0 < meta_concordant.pooled_effect_random < -26.0)
    check("p-value demonstrates high statistical significance (p < 0.001)", meta_concordant.p_value < 0.001)
    check("concordant sites yield low I^2 heterogeneity (< 25%)", meta_concordant.i_squared_pct < 25.0)
    check("replicability confirmed as statistically robust", meta_concordant.is_statistically_robust)
    check("leave-one-out sensitivity covers all 4 sites", len(meta_concordant.leave_one_out_sensitivity) == 4)

    # Divergent replication dataset (high heterogeneity)
    divergent_sites = [
        SiteReplicationResult(site_id="Site_A", n_pairs=10, mean_difference=-45.0, standard_error=2.5),
        SiteReplicationResult(site_id="Site_B", n_pairs=10, mean_difference=-5.0, standard_error=2.5),
        SiteReplicationResult(site_id="Site_C", n_pairs=10, mean_difference=-22.0, standard_error=2.5),
    ]
    meta_divergent = run_random_effects_meta_analysis(divergent_sites, target_gene="UNC13A")
    check("divergent sites trigger high I^2 heterogeneity (> 75%)", meta_divergent.i_squared_pct > 75.0)
    check("between-study variance tau^2 is substantial (> 100)", meta_divergent.between_study_variance_tau2 > 100.0)
    check("verdict reports substantial heterogeneity", "SUBSTANTIAL HETEROGENEITY" in meta_divergent.heterogeneity_verdict)

    # ------------------------------------------------------------------
    # 3. Workflows API integration
    # ------------------------------------------------------------------
    print("\n[protocol: workflows API serialization]")
    api_proto = generate_study_protocol("TARDBP", disease="als", paired_n=8)
    check("workflow generate_study_protocol returns target_gene", api_proto["target_gene"] == "TARDBP")
    check("workflow returns protocol_fingerprint_sha256", "protocol_fingerprint_sha256" in api_proto)
    check("workflow returns sample_plate_layout", len(api_proto["sample_plate_layout"]) == 16)

    sites_json = [
        {"site_id": "Lab1", "n_pairs": 8, "mean_difference": -15.0, "standard_error": 2.0},
        {"site_id": "Lab2", "n_pairs": 8, "mean_difference": -17.5, "standard_error": 2.2},
    ]
    api_meta = synthesize_multisite_replication(sites_json, target_gene="STMN2")
    check("workflow synthesize_multisite_replication returns k_sites", api_meta["k_sites"] == 2)
    check("workflow returns pooled_effect_random", api_meta["pooled_effect_random"] < -10.0)

    print("\n" + "=" * 70)
    print("All checks passed.")


if __name__ == "__main__":
    run_tests()
