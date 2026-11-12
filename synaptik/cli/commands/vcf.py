"""`synaptik vcf`: Ingest patient genomic VCF files, prioritize cryptic splice variants, and trigger N-of-1 therapeutic translation."""

from __future__ import annotations

import json
from pathlib import Path

from ...sequence.vcf import PatientVCFParser, variant_to_dict
from ..main import Fmt, header
from .cite import AUTHOR, ORCID_URL


def cmd_vcf(args) -> int:
    vcf_path = Path(args.vcf_file)
    if not vcf_path.exists():
        print(f"Error: VCF file not found at {vcf_path}")
        return 1

    target_genes = [args.target.upper()] if getattr(args, "target", None) else None
    parser = PatientVCFParser(target_genes=target_genes)
    all_variants = parser.parse_file(vcf_path)

    min_score = getattr(args, "min_splice_score", 0.20)
    prioritized = parser.prioritize_actionable_variants(all_variants, min_splice_score=min_score)

    if getattr(args, "json", False):
        payload = {
            "vcf_file": str(vcf_path),
            "total_variants_parsed": len(all_variants),
            "actionable_candidates_count": len(prioritized),
            "candidates": [variant_to_dict(v) for v in prioritized],
            "citation": {
                "author": AUTHOR,
                "orcid": "0009-0001-9477-3332",
                "orcid_url": ORCID_URL,
            },
        }
        print(json.dumps(payload, indent=2))
        return 0

    print(header(f"Patient Genomic VCF Splice Prioritization: {vcf_path.name}"))
    print(f"  Parsed {len(all_variants)} total variants from patient VCF.")
    print(f"  Identified {len(prioritized)} high-priority actionable splice candidates (min delta >= {min_score:.2f}).\n")

    if not prioritized:
        print("  No high-priority cryptic splice variants detected in neuro target catalog.")
        return 0

    print(Fmt.bold(f"{'VARIANT ID':<18} {'GENE':<10} {'CLINICAL SIG':<25} {'SPLICE EFFECT':<22} {'DELTA':<8} {'PRIORITY'}"))
    print("-" * 95)

    for v in prioritized:
        delta_str = f"{v.spliceai_max_delta:.2f}"
        print(
            f"{v.variant_id:<18} "
            f"{v.gene_symbol or 'Unknown':<10} "
            f"{v.clinical_significance[:24]:<25} "
            f"{v.splice_effect:<22} "
            f"{delta_str:<8} "
            f"{v.priority_score:>6.1f} / 100"
        )

    top = prioritized[0]
    print(f"\n{Fmt.bold('Top Actionable Therapeutic Target:')} {top.gene_symbol} ({top.variant_id})")
    print(f"  * Clinical Status: {top.clinical_significance}")
    print(f"  * Splice Alteration: {top.splice_effect} (SpliceAI delta = {top.spliceai_max_delta:.2f})")
    if top.gnomad_af is not None:
        print(f"  * Population Frequency (gnomAD): {top.gnomad_af:.6f}")
    if top.cadd_phred is not None:
        print(f"  * Pathogenicity Score (CADD Phred): {top.cadd_phred:.1f}")

    print(f"\n{Fmt.green('Ready to translate into personalized ASO & preclinical trial protocol:')}")
    print(f"  Run: {Fmt.bold(f'synaptik nof1 --gene {top.gene_symbol} --patient PATIENT-VCF --dossier IND_Dossier_{top.gene_symbol}.md')}")

    print(Fmt.dim(f"\nCitation: {AUTHOR} (ORCID: {ORCID_URL}). Use 'synaptik cite' to export BibTeX/APA."))
    return 0


def register(sub, common) -> None:
    s = sub.add_parser(
        "vcf", parents=[common],
        help="ingest patient genomic VCF, prioritize cryptic splice variants, and trigger translational therapy",
    )
    s.add_argument("vcf_file", metavar="VCF_FILE", help="path to patient VCF file (.vcf or .vcf.gz)")
    s.add_argument("--disease", choices=["als", "ftd", "ad"], default="als", help="disease context (default: als)")
    s.add_argument("--target", metavar="GENE", help="restrict to specific target gene (e.g. STMN2, UNC13A)")
    s.add_argument("--min-splice-score", type=float, default=0.20, help="minimum SpliceAI delta score (default: 0.20)")
    s.add_argument("--json", action="store_true", help="output prioritized variants as JSON")
    s.set_defaults(func=cmd_vcf)

