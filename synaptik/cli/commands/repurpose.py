"""`synaptik repurpose`: Small-molecule drug repurposing and phenotypic rescue engine."""

from __future__ import annotations

import json

from ...drug.repurposing import MoADirection
from ...drug.repurposing_catalog import (
    repurposed_compound_to_dict,
    screen_repurposing_candidates,
)
from ..main import Fmt, header
from .cite import AUTHOR, ORCID_URL


def cmd_repurpose(args) -> int:
    target_gene = args.target.strip().upper()
    min_phase = getattr(args, "min_phase", 0)
    bbb_only = getattr(args, "bbb_only", False)
    direction_raw = getattr(args, "direction", "activator")

    candidates = screen_repurposing_candidates(
        target_gene=target_gene,
        required_direction=direction_raw,
        min_clinical_phase=min_phase,
        bbb_only=bbb_only,
    )

    if getattr(args, "json", False):
        payload = {
            "target_gene": target_gene,
            "min_clinical_phase": min_phase,
            "bbb_only": bbb_only,
            "total_candidates_found": len(candidates),
            "candidates": [repurposed_compound_to_dict(c) for c in candidates],
            "citation": {
                "author": AUTHOR,
                "orcid": "0009-0001-9477-3332",
                "orcid_url": ORCID_URL,
            },
        }
        print(json.dumps(payload, indent=2))
        return 0

    print(header(f"Small-Molecule Drug Repurposing: {target_gene} Target Pathway"))
    print(f"  Target: {target_gene} | Desired MoA: {direction_raw.upper()} | Min Phase: {min_phase} | CNS BBB Filter: {'Active' if bbb_only else 'All'}\n")

    if not candidates:
        print(f"  No candidates matched the specified criteria for target '{target_gene}'.")
        return 0

    print(Fmt.bold(f"{'DRUG CANDIDATE':<24} {'PHASE':<16} {'CNS MPO':<10} {'BBB PENETRANCE':<16} {'PREDICTED RESCUE'}"))
    print("-" * 88)

    for c in candidates:
        bbb_verdict = Fmt.green("Plausible") if c.bbb_plausible else Fmt.yellow("Dish Only")
        rescue_str = f"{c.predicted_rescue_efficacy_pct:.1f}%"
        print(
            f"{c.name[:23]:<24} "
            f"{c.phase_label[:15]:<16} "
            f"{c.cns_mpo_score:>4.1f} / 6.0   "
            f"{bbb_verdict:<25} "
            f"{rescue_str}"
        )

    lead = candidates[0]
    print(f"\n{Fmt.bold('Lead Repurposing Candidate:')} {lead.name} ({lead.chembl_id})")
    print(f"  * Clinical Status: {lead.phase_label}")
    print(f"  * Primary Target: {lead.target_name} ({lead.target_symbol})")
    print(f"  * Mechanism of Action: {lead.mechanism_summary}")
    print(f"  * CNS MPO Druggability Score: {lead.cns_mpo_score:.2f} / 6.00 (MW={lead.mw_da} Da, clogP={lead.clogp}, TPSA={lead.tpsa_a2} Å²)")
    print(f"  * Literature Precedent: {lead.literature_precedent}")

    print(Fmt.dim(f"\nCitation: {AUTHOR} (ORCID: {ORCID_URL}). Use 'synaptik cite' to export BibTeX/APA."))
    return 0


def register(sub, common) -> None:
    s = sub.add_parser(
        "repurpose", parents=[common],
        help="screen and rank brain-penetrant small-molecule drugs for phenotype rescue",
    )
    s.add_argument("--target", required=True, metavar="GENE", help="target gene or pathway (e.g. STMN2, SMN2, MTOR, SOD1)")
    s.add_argument("--direction", choices=["activator", "inhibitor", "modulator"], default="activator", help="desired MoA direction (default: activator)")
    s.add_argument("--min-phase", type=int, default=0, choices=[0, 1, 2, 3, 4], help="minimum clinical phase (4=Approved, default: 0)")
    s.add_argument("--bbb-only", action="store_true", help="restrict strictly to blood-brain barrier penetrant compounds")
    s.add_argument("--json", action="store_true", help="output candidate molecules as JSON")
    s.set_defaults(func=cmd_repurpose)

