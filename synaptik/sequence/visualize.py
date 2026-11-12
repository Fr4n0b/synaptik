"""Publication-Grade Splicing & Therapeutic ASO Architecture Vector Visualizer (SVG).

Generates high-resolution, scalable vector graphics (SVG) illustrating
pre-mRNA exon-intron architecture, cryptic exon inclusion events,
sashimi-style splice junction arcs, and steric-blocking antisense
oligonucleotide (ASO) binding footprints.
"""

from __future__ import annotations

import html
from dataclasses import dataclass
from pathlib import Path


@dataclass
class SplicingDiagramConfig:
    width: int = 880
    height: int = 340
    gene_symbol: str = "STMN2"
    upstream_exon_name: str = "Exon 1"
    downstream_exon_name: str = "Exon 2"
    cryptic_exon_name: str = "Cryptic Exon 2a"
    delta_psi_pct: float = 85.0
    lead_aso_sequence: str = "CTGAATACTCACCTGAAA"
    lead_aso_chemistry: str = "2'-MOE / PS"
    melting_temp_c: float = 66.0
    author_citation: str = "Nobile, F. (2026). SYNAPTIK. ORCID: 0009-0001-9477-3332."


def generate_splicing_architecture_svg(config: SplicingDiagramConfig) -> str:
    """Construct an audit-ready, publication-grade SVG diagram."""
    w = config.width
    h = config.height

    # Gene / Exon geometry
    intron_y = 140
    intron_start_x = 70
    intron_end_x = w - 70

    exon_w = 90
    exon_h = 44
    exon_y = intron_y - exon_h // 2

    up_exon_x = intron_start_x
    down_exon_x = intron_end_x - exon_w

    cryptic_w = 110
    cryptic_h = 44
    cryptic_x = (intron_start_x + intron_end_x) // 2 - cryptic_w // 2

    # ASO footprint geometry
    aso_w = 160
    aso_h = 28
    aso_x = cryptic_x + (cryptic_w - aso_w) // 2
    aso_y = intron_y + 45

    # Arc controls
    # Canonical arc: up_exon_x + exon_w -> down_exon_x (top arc)
    c_start_x = up_exon_x + exon_w
    c_end_x = down_exon_x
    c_mid_x = (c_start_x + c_end_x) // 2
    c_ctrl_y = 40

    # Cryptic inclusion arc: up_exon_x + exon_w -> cryptic_x
    cr_start_x = up_exon_x + exon_w
    cr_end_x = cryptic_x
    cr_mid_x = (cr_start_x + cr_end_x) // 2
    cr_ctrl_y = 85

    svg_lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="100%" height="100%" style="font-family: -apple-system, BlinkMacSystemFont, \'Segoe UI\', Roboto, Helvetica, Arial, sans-serif; background: #0f172a; border-radius: 8px;">',
        "  <defs>",
        '    <linearGradient id="canonGrad" x1="0%" y1="0%" x2="100%" y2="100%">',
        '      <stop offset="0%" stop-color="#3b82f6" />',
        '      <stop offset="100%" stop-color="#1d4ed8" />',
        "    </linearGradient>",
        '    <linearGradient id="crypticGrad" x1="0%" y1="0%" x2="100%" y2="100%">',
        '      <stop offset="0%" stop-color="#ef4444" />',
        '      <stop offset="100%" stop-color="#b91c1c" />',
        "    </linearGradient>",
        '    <linearGradient id="asoGrad" x1="0%" y1="0%" x2="100%" y2="100%">',
        '      <stop offset="0%" stop-color="#f59e0b" />',
        '      <stop offset="100%" stop-color="#d97706" />',
        "    </linearGradient>",
        '    <filter id="glow" x="-20%" y="-20%" width="140%" height="140%">',
        '      <feDropShadow dx="0" dy="2" stdDeviation="4" flood-color="#000000" flood-opacity="0.4" />',
        "    </filter>",
        "  </defs>",
        "",
        "  <!-- Header & Target Metadata -->",
        f'  <text x="30" y="32" fill="#f8fafc" font-size="16" font-weight="700">Pre-mRNA Splicing &amp; Therapeutic ASO Steric-Block Footprint</text>',
        f'  <text x="30" y="52" fill="#94a3b8" font-size="12">Target Locus: <tspan fill="#38bdf8" font-weight="600">{html.escape(config.gene_symbol)}</tspan> | Cryptic Exclusion ΔΨ: <tspan fill="#f87171" font-weight="600">+{config.delta_psi_pct:.1f}%</tspan> | Lead Oligo Tm: <tspan fill="#fbbf24" font-weight="600">{config.melting_temp_c:.1f}°C</tspan></text>',
        "",
        "  <!-- Intron Line -->",
        f'  <line x1="{intron_start_x}" y1="{intron_y}" x2="{intron_end_x}" y2="{intron_y}" stroke="#475569" stroke-width="3" stroke-dasharray="4 4" />',
        "",
        "  <!-- Canonical Splicing Arc (Top, Blue) -->",
        f'  <path d="M {c_start_x} {intron_y} Q {c_mid_x} {c_ctrl_y} {c_end_x} {intron_y}" fill="none" stroke="#60a5fa" stroke-width="2.5" stroke-opacity="0.85" />',
        f'  <text x="{c_mid_x}" y="{c_ctrl_y - 8}" fill="#93c5fd" font-size="11" font-weight="600" text-anchor="middle">Canonical Splicing (Normal Full-Length Transcript)</text>',
        "",
        "  <!-- Cryptic Aberrant Splicing Arc (Red) -->",
        f'  <path d="M {cr_start_x} {intron_y} Q {cr_mid_x} {cr_ctrl_y} {cr_end_x} {intron_y}" fill="none" stroke="#f87171" stroke-width="2" stroke-dasharray="3 3" />',
        f'  <text x="{cr_mid_x}" y="{cr_ctrl_y + 14}" fill="#fca5a5" font-size="10" text-anchor="middle">Aberrant Cryptic Splice</text>',
        "",
        "  <!-- Upstream Canonical Exon -->",
        f'  <rect x="{up_exon_x}" y="{exon_y}" width="{exon_w}" height="{exon_h}" rx="6" fill="url(#canonGrad)" stroke="#60a5fa" stroke-width="1.5" filter="url(#glow)" />',
        f'  <text x="{up_exon_x + exon_w // 2}" y="{intron_y + 5}" fill="#ffffff" font-size="13" font-weight="600" text-anchor="middle">{html.escape(config.upstream_exon_name)}</text>',
        "",
        "  <!-- Downstream Canonical Exon -->",
        f'  <rect x="{down_exon_x}" y="{exon_y}" width="{exon_w}" height="{exon_h}" rx="6" fill="url(#canonGrad)" stroke="#60a5fa" stroke-width="1.5" filter="url(#glow)" />',
        f'  <text x="{down_exon_x + exon_w // 2}" y="{intron_y + 5}" fill="#ffffff" font-size="13" font-weight="600" text-anchor="middle">{html.escape(config.downstream_exon_name)}</text>',
        "",
        "  <!-- Cryptic Exon (Pathology) -->",
        f'  <rect x="{cryptic_x}" y="{exon_y}" width="{cryptic_w}" height="{cryptic_h}" rx="6" fill="url(#crypticGrad)" stroke="#f87171" stroke-width="1.5" filter="url(#glow)" />',
        f'  <text x="{cryptic_x + cryptic_w // 2}" y="{intron_y - 2}" fill="#ffffff" font-size="12" font-weight="700" text-anchor="middle">{html.escape(config.cryptic_exon_name)}</text>',
        f'  <text x="{cryptic_x + cryptic_w // 2}" y="{intron_y + 14}" fill="#fecaca" font-size="9" text-anchor="middle">Premature Poly-A / Stop</text>',
        "",
        "  <!-- Splice Donor/Acceptor Tags -->",
        f'  <text x="{cryptic_x - 12}" y="{intron_y + 4}" fill="#fca5a5" font-size="10" font-weight="700" text-anchor="end">AG</text>',
        f'  <text x="{cryptic_x + cryptic_w + 12}" y="{intron_y + 4}" fill="#fca5a5" font-size="10" font-weight="700" text-anchor="start">GT</text>',
        "",
        "  <!-- Steric-Blocking ASO Footprint -->",
        f'  <rect x="{aso_x}" y="{aso_y}" width="{aso_w}" height="{aso_h}" rx="5" fill="url(#asoGrad)" stroke="#fbbf24" stroke-width="1.5" filter="url(#glow)" />',
        f'  <text x="{aso_x + aso_w // 2}" y="{aso_y + 18}" fill="#ffffff" font-size="11" font-weight="700" text-anchor="middle">ASO Steric Blocker (18-mer)</text>',
        "",
        "  <!-- ASO Nucleotide Sequence Details -->",
        f'  <text x="{aso_x + aso_w // 2}" y="{aso_y + 42}" fill="#fde68a" font-family="monospace" font-size="11" font-weight="600" text-anchor="middle">5\'-{html.escape(config.lead_aso_sequence)}-\'3</text>',
        f'  <text x="{aso_x + aso_w // 2}" y="{aso_y + 57}" fill="#cbd5e1" font-size="10" text-anchor="middle">Chemistry: {html.escape(config.lead_aso_chemistry)} (Gap-free splice-masking)</text>',
        "",
        "  <!-- Steric Block Inhibit Arrow -->",
        f'  <line x1="{aso_x + aso_w // 2}" y1="{aso_y}" x2="{aso_x + aso_w // 2}" y2="{intron_y + exon_h // 2 + 6}" stroke="#fbbf24" stroke-width="2" />',
        f'  <line x1="{aso_x + aso_w // 2 - 8}" y1="{intron_y + exon_h // 2 + 6}" x2="{aso_x + aso_w // 2 + 8}" y2="{intron_y + exon_h // 2 + 6}" stroke="#fbbf24" stroke-width="2.5" />',
        "",
        "  <!-- Scale & Footer Citation -->",
        f'  <text x="30" y="{h - 18}" fill="#64748b" font-size="11">{html.escape(config.author_citation)}</text>',
        f'  <text x="{w - 30}" y="{h - 18}" fill="#64748b" font-size="11" text-anchor="end">SYNAPTIK Splicing CAD Engine</text>',
        "</svg>",
    ]

    return "\n".join(svg_lines)


def export_splicing_architecture_svg(config: SplicingDiagramConfig, out_path: str | Path) -> Path:
    """Render and save the SVG diagram to a file."""
    svg_content = generate_splicing_architecture_svg(config)
    p = Path(out_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(svg_content, encoding="utf-8")
    return p
