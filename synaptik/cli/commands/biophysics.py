"""
`synaptik biophysics`: 3D RNA Secondary Structure Folding & Steric Accessibility CLI.

Principal Investigator & Author:
Francesco Nobile (ORCID: 0009-0001-9477-3332)
Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School;
University of Pavia.
"""

from __future__ import annotations

import json

from ...biophysics.competition import compute_steric_displacement
from ...biophysics.folding import compute_target_accessibility, fold_rna_mfe
from ..main import Fmt, header

AUTHOR_CITE = "Francesco Nobile (ORCID: 0009-0001-9477-3332)"


def cmd_biophysics(args) -> int:
    subcmd = getattr(args, "subcommand", None) or "open"
    as_json = bool(getattr(args, "json", False))

    if subcmd == "open":
        seq = getattr(args, "target_sequence", None) or getattr(args, "sequence", None) or ""
        start = int(getattr(args, "oligo_start", 0) or 0)
        length = int(getattr(args, "oligo_len", 18) or 18)
        if not seq:
            print(Fmt.red("Error: target_sequence / --sequence is required."))
            return 1
        try:
            acc = compute_target_accessibility(seq, start, start + length)
            mfe, db = fold_rna_mfe(seq)
        except ValueError as err:
            print(Fmt.red(f"Biophysics error: {err}"))
            return 1

        if as_json:
            out = acc.to_dict()
            out["dot_bracket"] = db
            print(json.dumps(out, indent=2))
        else:
            print(header("SYNAPTIK 3D RNA Biophysics — Target Accessibility"))
            print(f"  Target window:      [{start}:{start + length}] ({acc.target_sequence})")
            print(f"  Full MFE:           {acc.delta_g_full_mfe:.2f} kcal/mol  {db}")
            print(f"  Constrained MFE:    {acc.delta_g_constrained_mfe:.2f} kcal/mol")
            print(f"  Opening Penalty:    {acc.delta_g_open:.2f} kcal/mol ({acc.accessibility_tier.value})")
            print(f"\n  Author: {AUTHOR_CITE}")
        return 0

    if subcmd == "displace":
        oligo_kd = float(getattr(args, "oligo_kd", 0.1) or 0.1)
        competitor = str(getattr(args, "competitor", "u1_snrnp") or "u1_snrnp")
        gene = getattr(args, "gene", None)
        oligo_conc = float(getattr(args, "oligo_conc", 150.0) or 150.0)
        comp_conc = float(getattr(args, "competitor_conc", 50.0) or 50.0)
        try:
            res = compute_steric_displacement(
                oligo_kd_nm=oligo_kd,
                competitor_name=competitor,
                target_gene=gene,
                oligo_conc_nm=oligo_conc,
                competitor_conc_nm=comp_conc,
            )
        except ValueError as err:
            print(Fmt.red(f"Displacement error: {err}"))
            return 1

        if as_json:
            print(json.dumps(res.to_dict(), indent=2))
        else:
            print(header("SYNAPTIK 3D RNA Biophysics — Spliceosome Steric Displacement"))
            print(f"  Competitor:         {res.competitor_name} (Kd = {res.competitor_kd_nm} nM)")
            print(f"  ASO Kd:             {res.oligo_kd_nm} nM")
            print(f"  Displacement Ratio: {res.displacement_ratio:.2f}x")
            print(f"  Fractional Occ.:    {res.fractional_occupancy * 100:.1f}%")
            print(f"  Dominant Binder:    {res.dominant_binder}")
            print(f"\n  Author: {AUTHOR_CITE}")
        return 0

    if subcmd == "scan":
        seq = getattr(args, "sequence", None) or getattr(args, "target_sequence", None) or ""
        oligo_len = int(getattr(args, "oligo_len", 18) or 18)
        gene = getattr(args, "gene", "TARGET") or "TARGET"
        if len(seq) < oligo_len:
            print(Fmt.red(f"Error: sequence length ({len(seq)}) < oligo_len ({oligo_len})."))
            return 1
        windows = []
        for i in range(0, len(seq) - oligo_len + 1):
            acc = compute_target_accessibility(seq, i, i + oligo_len)
            windows.append({
                "start": i,
                "end": i + oligo_len,
                "target": acc.target_sequence,
                "delta_g_open": acc.delta_g_open,
                "tier": acc.accessibility_tier.value,
            })
        windows.sort(key=lambda w: w["delta_g_open"])
        if as_json:
            print(json.dumps({"author": AUTHOR_CITE, "gene": gene, "windows": windows}, indent=2))
        else:
            print(header(f"SYNAPTIK Biophysics Accessibility Scan — {gene}"))
            for w in windows[:5]:
                print(f"  [{w['start']:2d}:{w['end']:2d}] {w['target']}  dG_open={w['delta_g_open']:5.2f} kcal/mol  ({w['tier']})")
            print(f"\n  Author: {AUTHOR_CITE}")
        return 0

    print(Fmt.red(f"Unknown biophysics subcommand: {subcmd}"))
    return 1


def register(sub, common) -> None:
    p = sub.add_parser(
        "biophysics",
        parents=[common],
        help="3D RNA MFE folding, opening penalty (Delta G_open), and spliceosomal steric displacement",
    )
    sp = p.add_subparsers(dest="subcommand")

    p_open = sp.add_parser("open", parents=[common], help="Compute target opening penalty Delta G_open")
    p_open.add_argument("--target-sequence", "--sequence", dest="target_sequence", required=True)
    p_open.add_argument("--oligo-start", type=int, default=0)
    p_open.add_argument("--oligo-len", type=int, default=18)

    p_disp = sp.add_parser("displace", parents=[common], help="Model ASO vs spliceosome steric competition")
    p_disp.add_argument("--gene", default=None)
    p_disp.add_argument("--oligo-kd", type=float, default=0.1)
    p_disp.add_argument("--competitor", default="u1_snrnp")
    p_disp.add_argument("--oligo-conc", type=float, default=150.0)
    p_disp.add_argument("--competitor-conc", type=float, default=50.0)

    p_scan = sp.add_parser("scan", parents=[common], help="Scan target RNA for accessible ASO windows")
    p_scan.add_argument("--gene", default="TARGET")
    p_scan.add_argument("--sequence", required=True)
    p_scan.add_argument("--oligo-len", type=int, default=18)

    p.set_defaults(func=cmd_biophysics)
