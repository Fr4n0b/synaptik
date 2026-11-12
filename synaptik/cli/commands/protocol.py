"""`synaptik protocol`: preclinical study design pre-registration and multi-site meta-analysis."""

from __future__ import annotations

import json
from dataclasses import asdict

from ...protocol.registry import (
    PreclinicalProtocol,
    SiteReplicationResult,
    create_preclinical_protocol,
    run_random_effects_meta_analysis,
)
from ..main import Fmt, header


def cmd_protocol(args) -> int:
    target = getattr(args, "target", "TARDBP") or "TARDBP"
    disease = getattr(args, "disease", "als") or "als"
    cell_model = getattr(args, "model", "ipsc_motor_neuron") or "ipsc_motor_neuron"
    div = int(getattr(args, "div", 28) or 28)
    assay = getattr(args, "assay", "imaging_synaptic_puncta") or "imaging_synaptic_puncta"
    endpoint = getattr(args, "endpoint", "SYN1/PSD95 colocalized puncta density") or "SYN1/PSD95 colocalized puncta density"
    pairs = int(getattr(args, "pairs", 8) or 8)
    rho = float(getattr(args, "rho", 0.75) or 0.75)

    # Mode 1: Multi-site meta-analysis synthesis
    if getattr(args, "meta_analysis", None):
        meta_file = args.meta_analysis
        sites: list[SiteReplicationResult] = []
        try:
            with open(meta_file, "r", encoding="utf-8") as f:
                data = json.load(f)
                for item in data.get("sites", []):
                    sites.append(SiteReplicationResult(
                        site_id=item["site_id"],
                        n_pairs=int(item["n_pairs"]),
                        mean_difference=float(item["mean_difference"]),
                        standard_error=float(item["standard_error"]),
                    ))
        except Exception as e:
            if getattr(args, "json", False):
                print(json.dumps({"error": f"Failed reading meta-analysis file: {e}"}))
            else:
                print(Fmt.red(f"Error reading meta-analysis file: {e}"))
            return 1

        if len(sites) < 2:
            print(Fmt.red("Error: Meta-analysis requires at least 2 independent sites in JSON file."))
            return 1

        res = run_random_effects_meta_analysis(sites, target_gene=target, endpoint=endpoint)

        if getattr(args, "json", False):
            print(json.dumps(asdict(res), indent=2))
            return 0

        header(f"Multi-Site Meta-Analysis Synthesis: {res.target_gene} ({res.endpoint})")
        print(f"Replication Sites Evaluated  : {res.k_sites} independent laboratories")
        print(f"Total Paired Replicates      : {res.total_n_pairs} matched isogenic pairs")
        print(f"Random-Effects Pooled Effect : {res.pooled_effect_random:+.2f} (SE: {res.pooled_se_random:.3f})")
        print(f"95% Confidence Interval      : [{res.ci95_lower:+.2f}, {res.ci95_upper:+.2f}]")
        print(f"Z-Score / P-Value            : z = {res.z_score:.2f}, p = {res.p_value:.4e}")
        print(f"Heterogeneity I^2 Index      : {res.i_squared_pct:.1f}% (Cochran Q = {res.cochrans_q:.2f}, df = {res.df})")
        print(f"Between-Study Variance tau^2 : {res.between_study_variance_tau2:.3f}")

        print("\n" + "=" * 78)
        print("Meta-Analytic Replicability Verdict")
        print("=" * 78)
        if res.is_statistically_robust:
            print(Fmt.green("REPLICABILITY CONFIRMED: Robust pooled effect across multi-site laboratories with acceptable heterogeneity."))
        else:
            print(Fmt.yellow("REPLICABILITY CAVEAT: Moderate/high heterogeneity or sensitivity to single-site exclusion detected."))
        print(f"  -> {res.heterogeneity_verdict}")

        print("\nLeave-One-Site-Out Sensitivity Analysis")
        print("=" * 78)
        print(f"{'Omitted Site':<24s} {'Pooled Effect (Excl)':<24s} {'Shift (Delta)':<18s}")
        print("-" * 78)
        for item in res.leave_one_out_sensitivity:
            site = item["omitted_site"]
            p_excl = item["pooled_effect_excluding_site"]
            d_shift = item["delta_from_overall"]
            print(f"{site:<24s} {p_excl:>+12.2f}            {d_shift:>+10.2f}")

        print("\n" + "-" * 78)
        print("Cite: Francesco Nobile (ORCID: 0009-0001-9477-3332); run 'synaptik cite'")
        return 0

    # Mode 2: Generate pre-registration study protocol
    protocol = create_preclinical_protocol(
        target_gene=target,
        disease=disease,
        cell_model=cell_model,
        maturity_div=div,
        assay_family=assay,
        primary_endpoint=endpoint,
        paired_n=pairs,
        rho=rho,
    )

    if getattr(args, "json", False):
        print(json.dumps(asdict(protocol), indent=2))
        return 0

    header(f"Preclinical Study Design Protocol: {protocol.target_gene} in {protocol.disease}")
    print(f"Cryptographic SHA-256 Hash   : {protocol.protocol_fingerprint_sha256}")
    print(f"Blinding Verification Token  : {protocol.blinding_verification_token}")
    print(f"Target Gene / Readout Assay  : {protocol.target_gene} | {protocol.assay_family}")
    print(f"Cell Model & Maturity        : {protocol.cell_model} (DIV {protocol.maturity_div})")
    print(f"Pre-Specified Primary Readout: {protocol.primary_endpoint} (direction: {protocol.primary_direction})")
    print(f"Matched Sample Size (Pairs)  : {protocol.recommended_n_pairs} isogenic pairs (assumed rho = {protocol.paired_correlation_rho:.2f})")
    print(f"Minimal Detectable Effect    : {protocol.minimal_detectable_effect_pct:.1f}%")
    print(f"Blinding Architecture        : {protocol.blinding_scheme.value.upper()}")
    print(f"Randomization Allocation     : {protocol.randomization_strategy}")
    print(f"Plate Evaporation Barrier    : Perimeter wells reserved for PBS (edge effect mitigation)")

    print("\n" + "=" * 78)
    print("Pre-Registration Methodological Rigor Checklist")
    print("=" * 78)
    for check_name, passed in protocol.rigor_checklist.items():
        status = Fmt.green("VERIFIED") if passed else Fmt.red("MISSING")
        fmt_name = check_name.replace("_", " ").title()
        print(f"  [{status}] {fmt_name}")

    print("\nRandomized Block 96-Well Plate Allocation (Excerpt)")
    print("=" * 78)
    items = list(protocol.sample_plate_layout.items())[:12]
    for well, assignment in items:
        print(f"  Well {well:<4s} -> {assignment}")

    print(f"\n  [OK] Pre-registration fingerprint generated. Export with --json to archive in repository.")
    print("  [Cite: Francesco Nobile (ORCID: 0009-0001-9477-3332); run 'synaptik cite']")
    return 0


def register(sub, common) -> None:
    s = sub.add_parser(
        "protocol", parents=[common],
        help="generate pre-registered preclinical study protocols and synthesize multi-site meta-analyses",
    )
    s.add_argument("--target", metavar="GENE", default="TARDBP",
                   help="target gene symbol (default: TARDBP)")
    s.add_argument("--disease", choices=["als", "ftd", "ad"], default="als",
                   help="disease model indication (default: als)")
    s.add_argument("--model", default="ipsc_motor_neuron",
                   help="cell model system (default: ipsc_motor_neuron)")
    s.add_argument("--div", type=int, default=28,
                   help="maturity days in vitro (default: 28)")
    s.add_argument("--assay", default="imaging_synaptic_puncta",
                   help="experimental assay family (default: imaging_synaptic_puncta)")
    s.add_argument("--endpoint", default="SYN1/PSD95 colocalized puncta density",
                   help="primary quantitative endpoint description")
    s.add_argument("--pairs", type=int, default=8,
                   help="number of matched isogenic pairs (default: 8)")
    s.add_argument("--rho", type=float, default=0.75,
                   help="paired isogenic correlation coefficient (default: 0.75)")
    s.add_argument("--meta-analysis", metavar="PATH",
                   help="path to a JSON file containing multi-site replication study results")
    s.add_argument("--json", action="store_true", help="output results as JSON")
    s.set_defaults(func=cmd_protocol)
