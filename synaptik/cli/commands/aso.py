"""`synaptik aso`: in silico therapeutic antisense oligonucleotide design."""

from __future__ import annotations

import json
from dataclasses import asdict

from ...drug.aso import ASODesigner, ASOMode
from ..main import Fmt, header


def cmd_aso(args) -> int:
    seq = (args.sequence or "").strip().upper()
    if not seq:
        print(Fmt.red("Error: A target sequence is required. Supply via --sequence <ACGTU...>."))
        return 1

    mode = ASOMode.GAPMER if args.mode == "gapmer" else ASOMode.SPLICE_SWITCH
    gene = (args.gene or "").strip().upper()
    region = args.region or ("splice junction" if mode == ASOMode.SPLICE_SWITCH else "transcript target")
    top_k = args.top or 3

    designer = ASODesigner()
    res = designer.design(
        target_sequence=seq,
        mode=mode,
        target_gene=gene,
        target_region=region,
        top_k=top_k,
        transcriptome=getattr(args, "transcriptome", None),
        max_mismatches=getattr(args, "max_mismatches", None),
    )

    if args.json:
        data = {
            "target_gene": res.target_gene,
            "target_region": res.target_region,
            "target_sequence": res.target_sequence,
            "mode": res.mode.value,
            "chemistry": res.chemistry.value,
            "n_windows_evaluated": res.n_windows_evaluated,
            "n_passed_filters": res.n_passed_filters,
            "precedent": res.precedent,
            "caveats": res.caveats,
            "candidates": [asdict(c) for c in res.candidates],
        }
        print(json.dumps(data, indent=2, default=str))
        return 0

    header(f"In Silico ASO Engineering: {gene or 'Target'} ({res.mode.value})")
    print(Fmt.dim(
        f"Architecture: {res.chemistry.value} | Clinical Precedent: {res.precedent}\n"
        f"Evaluated {res.n_windows_evaluated} tiling windows across {len(seq)} nt target sequence."
    ))
    print(f"  Passed biophysical & specificity filters: {res.n_passed_filters} / {res.n_windows_evaluated}")
    print()

    if not res.candidates:
        print(Fmt.yellow("No candidates passed all biophysical and off-target screening filters."))
        return 0

    header(f"Top {len(res.candidates)} Candidate ASOs (for experimental validation)")
    print(f"  {'Rank':<5s} {'ASO Sequence (5->3)':<24s} {'Tm (C)':>7s} {'dG_37':>8s} {'GC%':>6s} {'Score':>7s}  IDT Order Format")
    for rank, c in enumerate(res.candidates, start=1):
        score_fmt = Fmt.green(f"{c.composite_score:5.1f}") if c.composite_score >= 80 else f"{c.composite_score:5.1f}"
        print(f"  #{rank:<4d} {c.oligo_sequence_5_to_3:<24s} {c.melting_temp_c:7.1f} {c.gibbs_free_energy_37c:8.2f} {c.gc_content_pct:5.1f}% {score_fmt:>7s} {c.vendor_idt_format}")

    print()
    header("Lead Candidate Specifications")
    lead = res.candidates[0]
    print(f"  Lead Sequence (5'->3') : {Fmt.bold(lead.oligo_sequence_5_to_3)}")
    print(f"  Target RNA (5'->3')     : {lead.target_sequence_5_to_3} (pos {lead.target_start}-{lead.target_end})")
    print(f"  Duplex Thermodynamics   : Tm = {lead.melting_temp_c} C, dH = {lead.delta_h_kcal_mol} kcal/mol, dS = {lead.delta_s_cal_mol_k} cal/(mol*K), dG_37 = {lead.gibbs_free_energy_37c} kcal/mol")
    print(f"  Biophysical Liabilities : G-quadruplex={lead.has_quadruplex_risk}, Poly-T/U={lead.has_poly_u_risk}, Hairpin dG={lead.hairpin_delta_g} kcal/mol, 3'-Dimer={lead.self_dimer_score}")
    if lead.delta_g_open_kcal_mol is not None:
        print(f"  3D Target Accessibility : dG_open={lead.delta_g_open_kcal_mol:.2f} kcal/mol, dG_eff={lead.effective_dg_kcal_mol:.2f} kcal/mol ({lead.accessibility_tier})")
    sp = lead.specificity
    if lead.specificity_screened:
        counts = sp.counts_by_mismatch()
        summary = (", ".join(f"{n}mm: {c}" for n, c in counts.items())
                   if counts else "none found")
        style = Fmt.red if any(n <= 1 for n in counts) else (
            Fmt.yellow if counts else Fmt.green)
        print(f"  Off-Target Screen       : {style(summary)}  "
              f"[{sp.n_transcripts_searched:,} transcripts, "
              f"{sp.n_bases_searched:,} bp searched, "
              f"<={sp.max_mismatches_searched} mismatches]")
        for h in sp.hits[:5]:
            print(Fmt.dim(f"                            - {h.describe()}"))
    else:
        print(f"  Off-Target Screen       : {Fmt.red('NOT SCREENED')} "
              "(no --transcriptome given; no specificity evidence)")
    print(f"  Molecular Weight        : {lead.molecular_weight_da:,.1f} Da (sodium salt)")
    print(f"  CSF Half-Life (class)   : ~{lead.chemistry_class_csf_half_life_days} days "
          f"for {lead.chemistry.value} chemistry "
          + Fmt.dim("(literature class value, not a per-sequence prediction)"))
    print(f"  Vendor Order String     : {lead.vendor_idt_format}")

    print()
    print(Fmt.bold("Before ordering this oligo"))
    for item in (
        "Off-target screen against a current full transcriptome "
        "(BLAST/Bowtie), not only this tool.",
        "Confirm the intended splicing or knockdown effect in cells "
        "(RT-PCR / RNA-seq) before drawing any conclusion.",
        "Transcriptome-wide RNA-seq of treated cells to detect off-target "
        "effects this sequence search cannot predict.",
        "Any use beyond cultured cells requires toxicology, regulatory and "
        "ethics review. This tool does not assess any of those.",
    ):
        print(f"  [ ] {item}")

    if lead.penalties_and_warnings:
        print()
        print(Fmt.yellow("Warnings:"))
        for w in lead.penalties_and_warnings:
            print(Fmt.yellow(f"  ! {w}"))

    if args.export:
        try:
            with open(args.export, "w", encoding="utf-8") as f:
                f.write("rank,oligo_sequence_5_to_3,target_sequence_5_to_3,length,tm_c,delta_g_37,gc_pct,score,molecular_weight_da,vendor_idt_format\n")
                for r, c in enumerate(res.candidates, start=1):
                    f.write(f"{r},{c.oligo_sequence_5_to_3},{c.target_sequence_5_to_3},{c.length},{c.melting_temp_c},{c.gibbs_free_energy_37c},{c.gc_content_pct},{c.composite_score},{c.molecular_weight_da},{c.vendor_idt_format}\n")
            print()
            print(Fmt.green(f"Exported {len(res.candidates)} candidates to {args.export}"))
        except OSError as e:
            print(Fmt.red(f"Failed to export: {e}"))

    if getattr(args, "svg", None) and res.candidates:
        from ...sequence.visualize import SplicingDiagramConfig, export_splicing_architecture_svg
        cfg = SplicingDiagramConfig(
            gene_symbol=gene or "Target",
            lead_aso_sequence=lead.oligo_sequence_5_to_3,
            melting_temp_c=lead.melting_temp_c,
        )
        try:
            export_splicing_architecture_svg(cfg, args.svg)
            print()
            print(Fmt.green(f"Exported publication-grade SVG splicing diagram to {args.svg}"))
        except OSError as e:
            print(Fmt.red(f"Failed to export SVG diagram: {e}"))

    print("\n" + "-" * 78)
    print("Cite: Francesco Nobile (ORCID: 0009-0001-9477-3332); run 'synaptik cite'")
    return 0


def register(sub, common) -> None:
    s = sub.add_parser(
        "aso", parents=[common],
        help="in silico engineering of synthesis-ready antisense oligonucleotides (ASOs)",
    )
    s.add_argument("--sequence", required=True,
                   help="target RNA or DNA sequence (e.g. 40-100 nt around splice site or variant)")
    s.add_argument("--mode", choices=["splice_switch", "gapmer"], default="splice_switch",
                   help="therapeutic mechanism (splice_switch for steric blocking, gapmer for RNase H degradation)")
    s.add_argument("--gene", help="target gene symbol (e.g. STMN2, SOD1, UNC13A)")
    s.add_argument("--region", help="description of target region (e.g. 'cryptic donor', 'exon 7 ISS')")
    s.add_argument("--top", type=int, default=3, help="number of top candidate ASOs to display (default: 3)")
    s.add_argument("--export", metavar="PATH", help="export candidates to CSV file")
    s.add_argument("--svg", metavar="PATH", help="export publication-grade SVG splicing and ASO architecture diagram")
    s.add_argument("--transcriptome", metavar="FASTA",
                   help="reference transcripts (FASTA, .gz ok) to screen "
                        "candidates against. Without this no specificity "
                        "claim is made.")
    s.add_argument("--max-mismatches", type=int,
                   help="off-target mismatch tolerance (default: 3 for "
                        "gapmers, 2 for splice-switching)")
    s.add_argument("--json", action="store_true", help="output candidate specifications as JSON")
    s.set_defaults(func=cmd_aso)
