"""Turnkey PI-ready Lab Audit dossier: `synaptik audit --gene <GENE>`."""

from __future__ import annotations

import json
from pathlib import Path

from ...drug.aso import ASODesigner, ASOMode
from ...infer.association import AssociationEngine
from ...normalize.matcher import PhenotypeMatcher
from ...schema import Direction
from ...store import Store


def _build_audit_dict(
    gene: str,
    phenotype: str | None,
    db_path: str,
    target_seq: str | None,
    mode: str,
) -> dict:
    gene_u = gene.strip().upper()
    records = []
    db_exists = Path(db_path).exists()
    if db_exists:
        store = Store(db_path)
        try:
            all_recs = list(store.iter_records())
        finally:
            store.close()
        for r in all_recs:
            t = (r.perturbation.target.raw if r.perturbation.target else "").strip().upper()
            if t != gene_u:
                continue
            if phenotype:
                pm = PhenotypeMatcher()
                q_term = pm.match(phenotype)
                q_key = q_term.key() if q_term else phenotype.strip().lower()
                r_key = r.outcome.phenotype.key().lower()
                if q_key.lower() not in r_key and phenotype.strip().lower() not in r.outcome.phenotype.raw.lower():
                    continue
            records.append(r)

    labs = {r.source.lab_id for r in records if r.source.lab_id}
    papers = {r.source.source_id for r in records}
    donors = [r.cell_model.n_donor_lines for r in records if r.cell_model.n_donor_lines is not None]
    clones = [r.cell_model.n_clones for r in records if r.cell_model.n_clones is not None]
    n_isogenic = sum(1 for r in records if r.cell_model.is_isogenic_pair is True)
    n_null = sum(1 for r in records if r.outcome.direction == Direction.NO_CHANGE)
    assays = sorted({r.assay.assay_class.value for r in records})

    top_assocs = []
    if records:
        engine = AssociationEngine()
        pheno_label = phenotype or (records[0].outcome.phenotype.ontology_label or records[0].outcome.phenotype.raw)
        assocs = engine.rank_for_phenotype(pheno_label, records)
        for a in assocs[:5]:
            top_assocs.append({
                "pathway": a.pathway_label,
                "tier": a.tier.value,
                "n_labs": a.n_labs,
                "n_records": a.n_records,
                "survives_lolo": a.survives_leave_one_lab_out,
                "direction_consistency": round(a.direction_consistency, 3),
                "caveats": a.caveats,
            })

    aso_candidates = []
    if target_seq:
        designer = ASODesigner()
        aso_mode = ASOMode.SPLICE_SWITCH if mode == "splice" else ASOMode.GAPMER
        res = designer.design(target_seq, mode=aso_mode, target_gene=gene_u, top_k=5)
        for c in res.candidates[:5]:
            aso_candidates.append({
                "oligo_5_to_3": c.oligo_sequence_5_to_3,
                "length": c.length,
                "tm_c": c.melting_temp_c,
                "dg37_kcal_mol": c.gibbs_free_energy_37c,
                "composite_score": c.composite_score,
                "empirical_knockdown_pct": c.empirical_knockdown_pct,
                "specificity_screened": c.specificity_screened,
            })

    caveats = []
    if not db_exists:
        caveats.append(f"Local corpus DB not found at {db_path}; literature section is empty. Run `synaptik harvest` first.")
    elif not records:
        caveats.append(f"No records matched gene={gene_u}" + (f" phenotype={phenotype}" if phenotype else "") + " in local store.")
    if records and len(labs) < 2:
        caveats.append("Single-lab evidence only (<2 independent research groups).")
    if records and not donors:
        caveats.append("No donor line counts (n_donor_lines) extracted; verify multi-donor replication manually.")
    if records and n_isogenic == 0:
        caveats.append("No CRISPR-corrected isogenic controls flagged in extracted records.")

    return {
        "gene": gene_u,
        "phenotype_filter": phenotype,
        "literature_audit": {
            "n_records": len(records),
            "n_papers": len(papers),
            "n_independent_labs": len(labs),
            "max_donor_lines_reported": max(donors) if donors else None,
            "max_clones_reported": max(clones) if clones else None,
            "n_isogenic_controlled_records": n_isogenic,
            "n_null_results": n_null,
            "assay_classes": assays,
            "top_pathway_associations": top_assocs,
        },
        "aso_shortlist": aso_candidates,
        "caveats": caveats,
    }


def _format_markdown(d: dict) -> str:
    lit = d["literature_audit"]
    lines = [
        f"# SYNAPTIK Lab Audit Dossier: {d['gene']}",
        "",
        f"**Phenotype filter**: {d['phenotype_filter'] or 'All phenotypes'}",
        "",
        "## 1. Literature Evidence-Structure Audit",
        f"- **Records / Papers**: {lit['n_records']} records across {lit['n_papers']} papers",
        f"- **Independent Research Groups**: {lit['n_independent_labs']}",
        f"- **Max Donor Lines Reported**: {lit['max_donor_lines_reported'] if lit['max_donor_lines_reported'] is not None else 'Not reported'}",
        f"- **Max Clones Reported**: {lit['max_clones_reported'] if lit['max_clones_reported'] is not None else 'Not reported'}",
        f"- **Isogenic-Controlled Records**: {lit['n_isogenic_controlled_records']}",
        f"- **Null / Unchanged Results Retained**: {lit['n_null_results']}",
        f"- **Assay Diversity**: {', '.join(lit['assay_classes']) if lit['assay_classes'] else 'None'}",
    ]
    if lit["top_pathway_associations"]:
        lines.append("")
        lines.append("### Top Pathway Associations")
        for a in lit["top_pathway_associations"]:
            lines.append(
                f"- **{a['pathway']}** — tier=`{a['tier']}`, labs={a['n_labs']}, "
                f"records={a['n_records']}, LOLO_survives={a['survives_lolo']}, "
                f"consistency={a['direction_consistency']}"
            )
    if d["aso_shortlist"]:
        lines.append("")
        lines.append("## 2. Top ASO Candidates")
        for idx, c in enumerate(d["aso_shortlist"], 1):
            kd_str = f", emp_KD={c['empirical_knockdown_pct']}%" if c["empirical_knockdown_pct"] is not None else ""
            lines.append(
                f"{idx}. `5'-{c['oligo_5_to_3']}-3'` ({c['length']}mer, Tm={c['tm_c']}C, "
                f"dG37={c['dg37_kcal_mol']} kcal/mol, score={c['composite_score']}{kd_str}, "
                f"transcriptome_screened={c['specificity_screened']})"
            )
    if d["caveats"]:
        lines.append("")
        lines.append("## Caveats")
        for cv in d["caveats"]:
            lines.append(f"- {cv}")
    return "\n".join(lines) + "\n"


def cmd_audit(args) -> int:
    report = _build_audit_dict(
        gene=args.gene,
        phenotype=args.phenotype,
        db_path=args.db,
        target_seq=args.target_seq,
        mode=args.mode,
    )
    if args.json:
        text = json.dumps(report, indent=2)
    else:
        text = _format_markdown(report)

    if args.out:
        Path(args.out).write_text(text, encoding="utf-8")
        print(f"Wrote audit dossier to {args.out}")
    else:
        print(text)
    return 0


def register(sub, common) -> None:
    p = sub.add_parser(
        "audit",
        parents=[common],
        help="turnkey PI-ready gene/phenotype evidence & ASO audit dossier",
    )
    p.add_argument("--gene", required=True, help="target gene symbol (e.g. STMN2, TARDBP, SOD1)")
    p.add_argument("--phenotype", default=None, help="optional phenotype filter")
    p.add_argument("--target-seq", default=None, help="optional target RNA/DNA sequence for ASO shortlist")
    p.add_argument("--mode", choices=["gapmer", "splice"], default="gapmer", help="ASO mode when --target-seq is given")
    p.add_argument("--out", default=None, help="write dossier to file instead of stdout")
    p.add_argument("--json", action="store_true", help="output JSON instead of Markdown")
    p.set_defaults(func=cmd_audit)
