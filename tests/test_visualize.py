"""Tests for Splicing & ASO Architecture Vector Visualization (`synaptik/sequence/visualize.py`)."""

from pathlib import Path
from synaptik.sequence.visualize import (
    SplicingDiagramConfig,
    export_splicing_architecture_svg,
    generate_splicing_architecture_svg,
)


def check(desc, cond):
    if cond:
        print(f"  PASS  {desc}")
    else:
        print(f"  FAIL  {desc}")
        raise AssertionError(desc)


def run_tests():
    print("=" * 70)
    print("SYNAPTIK Splicing & ASO Architecture SVG Visualizer test suite")
    print("=" * 70)

    # ------------------------------------------------------------------
    # 1. Default SVG generation
    # ------------------------------------------------------------------
    print("\n[default svg generation]")
    cfg = SplicingDiagramConfig()
    svg = generate_splicing_architecture_svg(cfg)

    check("svg starts with <svg", svg.strip().startswith("<svg"))
    check("svg ends with </svg>", svg.strip().endswith("</svg>"))
    check("svg contains xmlns attribute", 'xmlns="http://www.w3.org/2000/svg"' in svg)
    check("svg contains viewBox", f'viewBox="0 0 {cfg.width} {cfg.height}"' in svg)

    # ------------------------------------------------------------------
    # 2. Structural biological elements present
    # ------------------------------------------------------------------
    print("\n[biological structural elements]")
    check("svg contains gene symbol STMN2", "STMN2" in svg)
    check("svg contains upstream Exon 1", "Exon 1" in svg)
    check("svg contains downstream Exon 2", "Exon 2" in svg)
    check("svg contains Cryptic Exon 2a", "Cryptic Exon 2a" in svg)
    check("svg contains splice acceptor AG", ">AG<" in svg)
    check("svg contains splice donor GT", ">GT<" in svg)
    check("svg contains canonical splicing arc", "Canonical Splicing" in svg)
    check("svg contains aberrant cryptic splice arc", "Aberrant Cryptic Splice" in svg)
    check("svg contains ASO steric-blocker label", "ASO Steric Blocker" in svg)
    check("svg contains lead oligo sequence", cfg.lead_aso_sequence in svg)
    check("svg contains author citation Nobile", "Nobile" in svg)
    check("svg contains ORCID identifier", "0009-0001-9477-3332" in svg)

    # ------------------------------------------------------------------
    # 3. Custom configuration rendering
    # ------------------------------------------------------------------
    print("\n[custom configuration rendering]")
    custom_cfg = SplicingDiagramConfig(
        gene_symbol="UNC13A",
        upstream_exon_name="Exon 20",
        downstream_exon_name="Exon 21",
        cryptic_exon_name="Cryptic Exon 20b",
        lead_aso_sequence="GCTACAGTTATCCCTCCA",
        melting_temp_c=68.5,
    )
    custom_svg = generate_splicing_architecture_svg(custom_cfg)
    check("custom svg contains UNC13A", "UNC13A" in custom_svg)
    check("custom svg contains Exon 20", "Exon 20" in custom_svg)
    check("custom svg contains Exon 21", "Exon 21" in custom_svg)
    check("custom svg contains Cryptic Exon 20b", "Cryptic Exon 20b" in custom_svg)
    check("custom svg contains custom ASO sequence", "GCTACAGTTATCCCTCCA" in custom_svg)

    # ------------------------------------------------------------------
    # 4. Export to file on disk
    # ------------------------------------------------------------------
    print("\n[file export]")
    out_file = Path(__file__).parent / "test_out_splicing.svg"
    try:
        p = export_splicing_architecture_svg(custom_cfg, out_file)
        check("exported SVG file exists", p.exists())
        check("exported SVG file size > 1000 bytes", p.stat().st_size > 1000)
    finally:
        if out_file.exists():
            out_file.unlink()

    print("\n" + "=" * 70)
    print("ALL SVG VISUALIZATION TESTS PASSED")
    print("=" * 70)


if __name__ == "__main__":
    run_tests()
