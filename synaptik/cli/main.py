"""
SYNAPTIK command-line interface.

Design principle: every command prints the evidence structure, not just the
answer. A user who runs `synaptik phenotype "SYN1 puncta"` and sees a ranked
list without lab counts and caveats would draw exactly the wrong conclusions,
so the caveats are not optional output.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from pathlib import Path

from ..analyst import Analyst
from ..drug.repurposing import (
    ChEMBLClient, OpenTargetsClient, RepurposingEngine,
)
from ..extract.llm import LLMExtractor
from ..infer.lookup import MatchTier
from ..infer.weighting import EvidenceWeighter
from ..ingest.http import PoliteSession
from ..pipeline import CORPUS_QUERIES, Pipeline
from ..schema import (
    AssayClass, CellModelClass, Direction, ExperimentQuery, PerturbationType,
)
from ..store import Store

DEFAULT_DB = os.environ.get("SYNAPTIK_DB", "data/synaptik.db")


# --------------------------------------------------------------------------
# output helpers
# --------------------------------------------------------------------------

class Fmt:
    """ANSI formatting, disabled when not a TTY or when NO_COLOR is set."""

    enabled = sys.stdout.isatty() and not os.environ.get("NO_COLOR")

    @classmethod
    def _w(cls, code: str, s: str) -> str:
        return f"\033[{code}m{s}\033[0m" if cls.enabled else s

    @classmethod
    def bold(cls, s): return cls._w("1", s)
    @classmethod
    def dim(cls, s): return cls._w("2", s)
    @classmethod
    def red(cls, s): return cls._w("31", s)
    @classmethod
    def green(cls, s): return cls._w("32", s)
    @classmethod
    def yellow(cls, s): return cls._w("33", s)
    @classmethod
    def blue(cls, s): return cls._w("34", s)
    @classmethod
    def cyan(cls, s): return cls._w("36", s)


TIER_STYLE = {
    "strong": Fmt.green, "moderate": Fmt.cyan, "weak": Fmt.yellow,
    "single-lab": Fmt.red, "insufficient": Fmt.dim,
}


def hr(char: str = "-", width: int = 78) -> None:
    print(Fmt.dim(char * width))


def header(text: str) -> None:
    print()
    print(Fmt.bold(text))
    hr("=")


# --------------------------------------------------------------------------
# commands
# --------------------------------------------------------------------------


def cmd_init(args) -> int:
    store = Store(args.db)
    print(f"Initialized knowledge base at {Fmt.bold(args.db)}")
    for k, v in store.stats().items():
        print(f"  {k}: {v}")
    return 0


def cmd_preflight(args) -> int:
    """Validate live API responses against the parser contracts."""
    store = Store(args.db)
    session = PoliteSession(store=store, contact_email=args.email)
    pipe = Pipeline(store, session=session, contact_email=args.email,
                    ncbi_api_key=args.ncbi_key, strict_contracts=False)
    header("Preflight: API contract check")
    print(Fmt.dim("Fetching a few records from each API and checking that the "
                  "fields\nour parsers depend on are still present."))
    print()
    ok, report = pipe.preflight(progress=print)
    print()
    print(report)
    print()
    if ok:
        print(Fmt.green("PASS - safe to harvest."))
        return 0
    print(Fmt.red("FAIL - do not harvest until this is fixed."))
    return 1


def cmd_contracts(args) -> int:
    """Show or export the API contracts without touching the network."""
    from ..ingest.contracts import ALL_CONTRACTS, dump_contracts, self_test

    if args.export:
        dump_contracts(args.export)
        print(f"Contracts written to {args.export}")
        return 0
    if args.self_test:
        header("Contract self-test (offline, against fixtures)")
        ok = self_test()
        print()
        print(Fmt.green("PASS") if ok else Fmt.red("FAIL"))
        return 0 if ok else 1

    header("API response contracts")
    for key, c in ALL_CONTRACTS.items():
        print(Fmt.bold(f"{c.endpoint}") + Fmt.dim(f"   [{key}]"))
        if c.root_path:
            print(Fmt.dim(f"  records at: {c.root_path}"))
        for f_ in c.fields:
            tag = Fmt.red("required") if f_.required else Fmt.dim("optional")
            print(f"  {f_.path:42s} {tag}" +
                  (Fmt.dim(f"  {f_.note}") if f_.note else ""))
        print()
    return 0


def cmd_harvest(args) -> int:
    store = Store(args.db)
    session = PoliteSession(store=store, contact_email=args.email)
    pipe = Pipeline(
        store, session=session, contact_email=args.email,
        ncbi_api_key=args.ncbi_key, scimago_csv=args.scimago,
        strict_contracts=not args.allow_drift,
    )

    if not args.skip_preflight:
        ok, report = pipe.preflight(progress=lambda m: print(Fmt.dim(m)))
        if not ok:
            print()
            print(report)
            print()
            print(Fmt.red("Preflight failed. Harvest aborted."))
            print(Fmt.dim("Re-run with --skip-preflight to override, but "
                          "expect records with missing fields."))
            return 1
        print(Fmt.green("preflight OK"))
    queries = CORPUS_QUERIES
    if args.query:
        queries = {"custom": args.query}
    elif args.only:
        missing = [k for k in args.only if k not in CORPUS_QUERIES]
        if missing:
            print(Fmt.red(f"Unknown query set(s): {', '.join(missing)}"))
            print(f"Available: {', '.join(CORPUS_QUERIES)}")
            return 1
        queries = {k: CORPUS_QUERIES[k] for k in args.only}

    header("Harvesting corpus metadata")
    st = pipe.harvest(queries, max_per_query=args.max_per_query, progress=print)
    print()
    if st.failed():
        print(Fmt.red(st.summary()))
        print(Fmt.dim(f"HTTP: {session.stats()}"))
        return 1
    print(st.summary())
    print(Fmt.dim(f"HTTP: {session.stats()}"))

    if args.fulltext:
        header("Fetching full text")
        st2 = pipe.fetch_full_texts(limit=args.fulltext_limit, progress=print)
        print(st2.summary())
    return 0


def cmd_extract(args) -> int:
    store = Store(args.db)
    if not os.environ.get("ANTHROPIC_API_KEY") and not args.dry_run:
        print(Fmt.red("ANTHROPIC_API_KEY is not set."))
        print("Set it, or use --dry-run to see how many papers are queued.")
        return 1

    if args.dry_run:
        n = store.conn.execute(
            """SELECT COUNT(*) c FROM sources s JOIN documents d
               ON d.source_id = s.source_id
               WHERE s.source_id NOT IN (SELECT DISTINCT source_id FROM records)
                 AND s.is_review = 0"""
        ).fetchone()["c"]
        print(f"{n} papers queued for extraction.")
        print(Fmt.dim(f"Estimated cost at ~$0.02/paper: ${n * 0.02:.2f}"))
        return 0

    extractor = LLMExtractor(model=args.model)
    pipe = Pipeline(store, extractor=extractor, scimago_csv=args.scimago)
    header("Extracting experiment records")
    st = pipe.extract_records(limit=args.limit, progress=print)
    print()
    print(st.summary())
    if st.unmapped_phenotypes:
        print()
        print(Fmt.yellow("Most frequent unmapped phenotypes:"))
        print(Fmt.dim("  (add these as synonyms in normalize/phenotype_ontology.py)"))
        for term, count in st.unmapped_phenotypes[:15]:
            print(f"  {count:4d}  {term}")
    return 0


def cmd_renormalize(args) -> int:
    """Re-apply the phenotype vocabulary to existing records.

    Free and fast: it does NOT re-run the LLM. Use after adding synonyms to
    normalize/phenotype_ontology.py so previously-unmapped records become
    visible to analysis.
    """
    from ..pipeline import Pipeline

    store = Store(args.db)
    pipe = Pipeline(store, scimago_csv=args.scimago, strict_contracts=False)

    before = store.conn.execute(
        "SELECT COUNT(*) c FROM records WHERE phenotype_key LIKE 'raw:%'"
    ).fetchone()["c"]
    total = store.conn.execute("SELECT COUNT(*) c FROM records").fetchone()["c"]

    if not total:
        print(Fmt.yellow("No records in the knowledge base yet."))
        return 0

    header("Re-normalizing phenotype vocabulary")
    print(Fmt.dim(f"{total} records; {before} currently unmapped "
                  f"({before / max(1, total):.0%})"))
    print(Fmt.dim("This does not call the LLM and costs nothing."))
    print()

    n = 0
    for rec in list(store.iter_records()):
        old_id = rec.record_id
        old_key = rec.outcome.phenotype.key()
        updated = pipe.normalize_record(rec)
        if updated.outcome.phenotype.key() != old_key:
            n += 1
        # Record IDs are content hashes that include the phenotype key, so a
        # successful re-mapping produces a NEW id. Inserting without deleting
        # the old row leaves both behind and silently inflates every count in
        # the system, which would corrupt the replication statistics.
        if updated.record_id != old_id:
            store.conn.execute(
                "DELETE FROM records WHERE record_id=?", (old_id,)
            )
        store.insert_record(updated)
    store.conn.commit()

    after = store.conn.execute(
        "SELECT COUNT(*) c FROM records WHERE phenotype_key LIKE 'raw:%'"
    ).fetchone()["c"]

    print(f"  {n} record(s) changed mapping")
    print(f"  unmapped: {before} -> {Fmt.bold(str(after))} "
          f"({after / max(1, total):.0%})")
    if after < before:
        print(Fmt.green(f"  recovered {before - after} record(s) for analysis"))
    elif before and after == before:
        print(Fmt.yellow(
            "  No change. Check that your synonyms match the exact phrasings "
            "reported by `synaptik report` under unmapped phenotypes."))

    if after:
        print()
        print(Fmt.bold("Most frequent still-unmapped phenotypes:"))
        rows = store.conn.execute(
            """SELECT phenotype_raw, COUNT(*) c FROM records
               WHERE phenotype_key LIKE 'raw:%'
               GROUP BY phenotype_key ORDER BY c DESC LIMIT 15"""
        ).fetchall()
        for r in rows:
            print(f"  {r['c']:5d}  {r['phenotype_raw']}")
        print()
        print(Fmt.dim("Add these to normalize/phenotype_ontology.py, then "
                      "re-run this command."))
    return 0


def cmd_lookup(args) -> int:
    store = Store(args.db)
    an = Analyst(store, EvidenceWeighter())

    q = ExperimentQuery(
        cell_model_class=CellModelClass(args.model) if args.model else None,
        days_in_vitro=args.div,
        perturbation_type=PerturbationType(args.perturbation)
        if args.perturbation else None,
        perturbation_target=args.target,
        variant_hgvs=args.variant,
        assay_class=AssayClass(args.assay) if args.assay else None,
        markers=args.markers or [],
        phenotype=args.phenotype,
    )

    header("Has this experiment been done?")
    parts = []
    if args.model: parts.append(f"system={args.model}")
    if args.div: parts.append(f"DIV={args.div}")
    if args.target: parts.append(f"target={args.target}")
    if args.assay: parts.append(f"assay={args.assay}")
    if args.phenotype: parts.append(f"phenotype='{args.phenotype}'")
    print(Fmt.dim("query: " + ", ".join(parts)))
    print()

    res = an.lookup_experiment(q, max_results=args.limit)

    if args.json:
        print(json.dumps({
            "verdict": res.verdict,
            "n_exact": res.n_exact, "n_close": res.n_close,
            "conflicting": res.conflicting,
            "matches": [{
                "source_id": m.record.source.source_id,
                "tier": m.tier.value, "similarity": m.similarity,
                "weight": m.weight, "summary": m.summary_line(),
                "matched_on": m.matched_on, "differs_on": m.differs_on,
                "protocol": m.protocol().__dict__,
            } for m in res.matches],
        }, indent=2, default=str))
        return 0

    style = Fmt.red if res.conflicting else Fmt.bold
    print(style(res.verdict))
    print()

    for tier in (MatchTier.EXACT, MatchTier.CLOSE, MatchTier.RELATED,
                 MatchTier.ADJACENT):
        group = res.by_tier(tier)
        if not group:
            continue
        label = {
            MatchTier.EXACT: "EXACT MATCHES - this experiment has been done",
            MatchTier.CLOSE: "CLOSE - differs in one dimension",
            MatchTier.RELATED: "RELATED - useful context, not a direct answer",
            MatchTier.ADJACENT: "ADJACENT - background",
        }[tier]
        print(Fmt.bold(label))
        for m in group[: args.per_tier]:
            print(f"  {m.summary_line()}")
            print(Fmt.dim(f"    {m.record.source.source_id}  "
                          f"sim={m.similarity:.2f} weight={m.weight:.2f}"))
            p = m.protocol()
            proto = []
            if p.differentiation_protocol:
                proto.append(p.differentiation_protocol)
            if p.days_in_vitro:
                proto.append(f"DIV{p.days_in_vitro}")
            if p.coculture:
                proto.append(p.coculture)
            if p.markers:
                proto.append("markers: " + "/".join(p.markers[:4]))
            if p.n_biological_replicates:
                proto.append(f"n={p.n_biological_replicates}")
            if p.is_isogenic_pair:
                proto.append("isogenic")
            if proto:
                print(Fmt.dim(f"    protocol: {'; '.join(proto)}"))
            if m.differs_on:
                print(Fmt.yellow(f"    differs: {'; '.join(m.differs_on)}"))
        if len(group) > args.per_tier:
            print(Fmt.dim(f"  ... and {len(group) - args.per_tier} more"))
        print()
    return 0


def cmd_phenotype(args) -> int:
    store = Store(args.db)
    repurposer = None
    if args.drugs:
        session = PoliteSession(store=store, contact_email=args.email)
        repurposer = RepurposingEngine(
            chembl=ChEMBLClient(session),
            opentargets=OpenTargetsClient(session),
        )
    an = Analyst(store, EvidenceWeighter(), repurposing=repurposer)

    direction = Direction(args.direction) if args.direction else Direction.UNSPECIFIED
    res = an.analyze_phenotype(
        args.phenotype, direction=direction, include_drugs=args.drugs
    )

    if args.json:
        print(json.dumps({
            "phenotype": res.phenotype_query,
            "resolved": res.resolved_phenotype,
            "label": res.resolved_label,
            "n_records": res.n_matching_records,
            "n_papers": res.n_papers, "n_labs": res.n_labs,
            "interpretation": res.interpretation,
            "warnings": res.warnings,
            "associations": [{
                "pathway": a.pathway_label, "tier": a.confidence_tier(),
                "enrichment": round(a.enrichment, 2), "q_value": a.q_value,
                "n_records": a.n_records, "n_papers": a.n_papers,
                "n_labs": a.n_labs, "n_assay_families": a.n_ortho_groups,
                "direction": a.dominant_direction,
                "direction_consistency": round(a.direction_consistency, 2),
                "survives_leave_one_lab_out": a.survives_lolo,
                "n_negative_results": a.n_negative_records,
                "median_effect_size": a.median_effect_size,
                "caveats": a.caveats(),
                "supporting": [s.__dict__ for s in a.supporting[:5]],
            } for a in res.top(args.limit)],
            "regulators": [{
                "rbp": r.rbp, "events": r.n_supporting_events,
                "labs": r.n_labs, "signature": r.signature_match,
                "genes": r.example_genes, "notes": r.therapeutic_notes,
            } for r in res.regulators],
            "drugs": [{
                "name": d.name, "chembl_id": d.chembl_id,
                "target": d.target_symbol, "moa": d.moa_direction.value,
                "stage": d.stage_label(), "potency_nm": d.potency_nm,
                "cns": d.cns.summary(), "rationale": d.rationale,
                "warnings": d.warnings,
                "prior_evidence": d.prior_neuronal_evidence,
            } for d in res.drug_candidates],
        }, indent=2, default=str))
        return 0

    header(f"Phenotype analysis: {res.phenotype_query}")
    if res.resolved_label:
        print(Fmt.dim(f"resolved to: {res.resolved_label} ({res.resolved_phenotype})"))
    print(f"{res.n_matching_records} records | {res.n_papers} papers | "
          f"{Fmt.bold(str(res.n_labs))} independent groups | "
          f"corpus {res.corpus_size}")
    print()
    print(res.interpretation)

    if res.warnings:
        print()
        for w in res.warnings:
            print(Fmt.yellow(f"  ! {w}"))

    if not res.associations:
        return 0

    header("Associated pathways")
    print(f"{'pathway':30s} {'tier':12s} {'enrich':>7s} {'q':>9s} "
          f"{'labs':>5s} {'assays':>7s} {'dir':>9s}")
    hr()
    for a in res.top(args.limit):
        tier = a.confidence_tier()
        style = TIER_STYLE.get(tier, str)
        print(f"{a.pathway_label[:29]:30s} {style(f'{tier:12s}')} "
              f"{a.enrichment:7.1f} {a.q_value:9.2e} {a.n_labs:5d} "
              f"{a.n_ortho_groups:7d} {a.dominant_direction:>9s}")

    if args.verbose:
        for a in res.top(args.limit):
            print()
            print(Fmt.bold(f"{a.pathway_label}  [{a.confidence_tier()}]"))
            print(f"  {a.n_records} records, {a.n_papers} papers, {a.n_labs} labs, "
                  f"{a.n_ortho_groups} assay families")
            print(f"  enrichment {a.enrichment:.1f}x (q={a.q_value:.2e}), "
                  f"direction consistency {a.direction_consistency:.0%}")
            if a.leave_one_lab_out_enrichment is not None:
                print(f"  leave-one-lab-out: {a.leave_one_lab_out_enrichment:.1f}x "
                      f"({'survives' if a.survives_lolo else Fmt.red('FAILS')})")
            if a.median_effect_size is not None:
                iqr = (f", IQR {a.effect_size_iqr[0]:.0f} to "
                       f"{a.effect_size_iqr[1]:.0f}" if a.effect_size_iqr else "")
                print(f"  median effect size {a.median_effect_size:.0f}{iqr} "
                      f"({a.n_with_effect_size} quantified)")
            if a.n_negative_records:
                print(f"  {a.n_negative_records} null result(s) reported")
            for c in a.caveats():
                print(Fmt.yellow(f"  ! {c}"))
            print(Fmt.dim("  supporting:"))
            for s in a.supporting[:3]:
                print(Fmt.dim(f"    {s.source_id} {s.journal or '?'} {s.year or ''} "
                              f"[{s.direction}] w={s.weight:.2f}"))

    if res.regulators:
        header("Splicing regulators implicated")
        for r in res.regulators:
            print(Fmt.bold(f"{r.rbp}") +
                  f"  {r.n_supporting_events} events, {r.n_papers} papers, "
                  f"{r.n_labs} labs")
            print(f"  {r.signature_match}")
            if r.example_genes:
                print(Fmt.dim(f"  genes: {', '.join(r.example_genes)}"))
            for n in r.therapeutic_notes[:2]:
                print(Fmt.blue(f"  -> {n}"))
            for c in r.caveats():
                print(Fmt.yellow(f"  ! {c}"))
            print()

    if res.drug_candidates:
        header("Candidate compounds")
        print(Fmt.dim("Direction-matched only: compounds whose mechanism opposes "
                      "the observed change."))
        print()
        for d in res.drug_candidates[: args.drug_limit]:
            print(Fmt.bold(f"{d.name}") + Fmt.dim(f"  ({d.chembl_id})"))
            print(f"  {d.target_symbol} {d.moa_direction.value} | "
                  f"{d.stage_label()} | " +
                  (f"{d.potency_nm:.0f} nM {d.potency_type}"
                   if d.potency_nm else "potency unknown"))
            print(f"  pathway: {d.pathway_label} [{d.pathway_evidence_tier}]")
            print(Fmt.dim(f"  {d.cns.summary()}"))
            print(Fmt.dim(f"  {d.rationale}"))
            for n in d.cns.notes:
                print(Fmt.yellow(f"  ! {n}"))
            for w in d.warnings:
                print(Fmt.red(f"  ! {w}"))
            for e in d.prior_neuronal_evidence:
                print(Fmt.green(f"  * already tested: {e}"))
            print()
    return 0


def cmd_splicing(args) -> int:
    store = Store(args.db)
    an = Analyst(store)
    res = an.analyze_splicing(gene=args.gene, regulator=args.regulator)

    if args.json:
        print(json.dumps({
            "n_events": res.get("n_events", 0),
            "genes": res.get("genes", []),
            "events": res.get("events", []),
            "regulators": [{
                "rbp": r.rbp, "events": r.n_supporting_events,
                "labs": r.n_labs, "consistency": round(r.consistency, 2),
                "signature": r.signature_match, "genes": r.example_genes,
                "therapeutic_notes": r.therapeutic_notes,
                "caveats": r.caveats(),
            } for r in res.get("regulators", [])],
        }, indent=2, default=str))
        return 0

    header("Splicing analysis")
    if not res.get("events"):
        print(res.get("message", "No events found."))
        return 0

    print(f"{res['n_events']} events across {len(res['genes'])} genes")
    print(Fmt.dim(f"genes: {', '.join(res['genes'][:20])}"))

    if res["regulators"]:
        header("Implicated regulators")
        for r in res["regulators"]:
            print(Fmt.bold(r.rbp) +
                  f"  {r.n_supporting_events} events, {r.n_labs} labs, "
                  f"consistency {r.consistency:.0%}")
            print(f"  {r.signature_match}")
            if r.profile:
                print(Fmt.dim(f"  motifs: {', '.join(r.profile.motifs)}"))
                print(Fmt.dim(f"  {r.profile.canonical_effect}"))
            for n in r.therapeutic_notes:
                print(Fmt.blue(f"  -> {n}"))
            for c in r.caveats():
                print(Fmt.yellow(f"  ! {c}"))
            print()

    header("Events and therapeutic entry points")
    for e in res["events"][: args.limit]:
        dpsi = f"{e['delta_psi']:+.2f}" if e["delta_psi"] is not None else "n/a"
        print(Fmt.bold(f"{e['gene']}") +
              f"  {e['event_type']}  dPSI={dpsi}  ({e['magnitude']})")
        if e["sign_confidence"] != "high":
            print(Fmt.yellow("  ! dPSI sign is uncertain -- verify against the figure"))
        nmd_style = Fmt.red if e["nmd_risk"] == "high" else Fmt.dim
        print(nmd_style(f"  NMD risk: {e['nmd_risk']} -- {e['nmd_note'][:90]}"))
        if e["aso_feasible"]:
            print(Fmt.green(f"  ASO: {e['aso_approach'][:100]}"))
            if e["aso_precedent"]:
                print(Fmt.dim(f"       precedent: {e['aso_precedent']}"))
        print(Fmt.dim(f"  {e['source_id']}"))
        print()
    return 0


def cmd_compare(args) -> int:
    store = Store(args.db)
    an = Analyst(store)
    res = an.compare_phenotypes(args.phenotype_a, args.phenotype_b)

    if args.json:
        print(json.dumps(res, indent=2, default=str))
        return 0

    header(f"Coherence: '{args.phenotype_a}' vs '{args.phenotype_b}'")
    a, b = res["phenotype_a"], res["phenotype_b"]
    print(f"  A: {a['n_records']} records, {a['n_labs']} labs")
    print(f"  B: {b['n_records']} records, {b['n_labs']} labs")
    jacc = res['jaccard']
    print(f"  Jaccard overlap: {Fmt.bold(f'{jacc:.2f}')}")
    print()
    print(res["verdict"])

    if res["shared_pathways"]:
        print()
        print(Fmt.bold("Shared pathways"))
        for s in res["shared_pathways"][:12]:
            mark = (Fmt.green("consistent") if s["consistent_direction"]
                    else Fmt.red("OPPOSITE directions"))
            print(f"  {s['pathway'][:34]:36s} "
                  f"{s['enrichment_a']:6.1f}x / {s['enrichment_b']:6.1f}x  {mark}")
    if res["only_a"]:
        print()
        print(Fmt.dim(f"Only in A: {', '.join(res['only_a'][:6])}"))
    if res["only_b"]:
        print(Fmt.dim(f"Only in B: {', '.join(res['only_b'][:6])}"))
    return 0


def cmd_gene(args) -> int:
    """Full multi-class evidence dossier for one gene."""
    from ..mechanism import MechanismAnalyst
    from ..omics.evidence_sources import (
        ClinVarClient, ExpressionClient, GnomADClient,
    )

    store = Store(args.db)
    gnomad = clinvar = expr = None
    if not args.offline:
        sess = PoliteSession(store=store, contact_email=args.email)
        gnomad, clinvar, expr = (
            GnomADClient(sess), ClinVarClient(sess), ExpressionClient(sess)
        )
    ma = MechanismAnalyst(store, gnomad=gnomad, clinvar=clinvar,
                          expression=expr, offline=args.offline)
    d = ma.gene_dossier(args.gene)

    if args.json:
        print(json.dumps({
            "gene": d.gene_symbol, "summary": d.summary,
            "evidence_classes": d.evidence_classes_present(),
            "n_records": d.n_literature_records, "n_papers": d.n_papers,
            "n_labs": d.n_labs, "splicing_events": d.splicing_events,
            "regulators": d.implicated_regulators,
            "phenotypes": d.literature_phenotypes,
            "contradictions": d.contradictions, "gaps": d.gaps,
            "hypotheses": [{
                "mechanism": h.mechanism, "concordance": h.concordance(),
                "convergence_score": h.reconciliation.convergence_score
                if h.reconciliation else 0,
                "verdict": h.verdict, "experiments": h.experiments,
            } for h in d.hypotheses],
        }, indent=2, default=str))
        return 0

    header(f"Evidence dossier: {d.gene_symbol}")
    print(d.summary)
    print()
    classes = d.evidence_classes_present()
    print(Fmt.bold("Evidence classes present: ") +
          (", ".join(classes) if classes else Fmt.red("none")))
    if d.n_literature_records:
        print(f"  literature   {d.n_literature_records} records, "
              f"{d.n_papers} papers, {Fmt.bold(str(d.n_labs))} groups")
    if d.constraint and not d.constraint.error:
        print(f"  population   {d.constraint.interpretation()[:120]}")
    if d.clinical and not d.clinical.error:
        print(f"  clinical     {d.clinical.interpretation()[:120]}")
    if d.functional and not d.functional.error:
        print(f"  functional   {d.functional.interpretation()[:120]}")
    if d.splicing_events:
        print(f"  splicing     {d.splicing_events} events; regulators: "
              f"{', '.join(d.implicated_regulators) or 'none inferred'}")

    if d.contradictions:
        header("Cross-class contradictions")
        print(Fmt.dim("These are the most informative findings here: one "
                      "evidence class disagrees with another."))
        print()
        for c in d.contradictions:
            print(Fmt.red("  ! ") + c)
            print()

    for h in d.hypotheses:
        header(f"Mechanism: {h.mechanism}")
        style = {"convergent": Fmt.green, "contradictory": Fmt.red,
                 "tension": Fmt.yellow}.get(h.concordance(), Fmt.dim)
        print(style(f"  {h.concordance().upper()}") +
              f"   convergence score "
              f"{h.reconciliation.convergence_score if h.reconciliation else 0}")
        print(f"  {h.verdict}")
        if h.experiments:
            print()
            print(Fmt.bold("  What would settle it:"))
            for e in h.experiments:
                print(Fmt.blue(f"    -> {e}"))

    if d.gaps:
        header("Evidence gaps")
        for g in d.gaps:
            print(Fmt.yellow(f"  - {g}"))
    return 0


def cmd_variant(args) -> int:
    """Reconcile variants against the literature."""
    from ..mechanism import MechanismAnalyst
    from ..sequence.splice import SpliceAIPredictor, SpliceEnsemble

    store = Store(args.db)
    preds = []
    if args.spliceai_vcf or args.reference:
        preds.append(SpliceAIPredictor(precomputed_vcf=args.spliceai_vcf,
                                       reference_fasta=args.reference))
    ma = MechanismAnalyst(store, splice_ensemble=SpliceEnsemble(preds or None),
                          offline=True)

    variants = []
    for v in args.variants:
        try:
            chrom, pos, change = v.split(":")
            ref, alt = change.split(">")
            variants.append((chrom, int(pos), ref, alt))
        except ValueError:
            print(Fmt.red(f"Cannot parse '{v}'. Use chr:pos:REF>ALT"))
            return 1

    panel = ma.analyze_variants(variants, gene_symbol=args.gene)

    if args.json:
        print(json.dumps({
            "summary": panel.summary(),
            "variants": [{
                "variant": r.variant_id, "gene": r.gene_symbol,
                "state": r.state.value, "priority": round(r.priority, 2),
                "predicted_delta": r.predicted_delta,
                "predicted_effect": r.predicted_effect.value,
                "n_observing_papers": r.n_observing_papers,
                "n_observing_labs": r.n_observing_labs,
                "narrative": r.narrative, "opportunity": r.opportunity,
                "caveats": r.caveats,
            } for r in panel.by_priority()],
        }, indent=2, default=str))
        return 0

    header("Sequence vs literature reconciliation")
    s_ = panel.summary()
    print(f"{s_['n_variants']} variant(s): {s_['states']}")
    print()
    novel = panel.novel_hypotheses()
    if novel:
        print(Fmt.green(Fmt.bold(
            f"{len(novel)} UNTESTED PREDICTION(S) -- highest experimental value")))
        print()
    for r in panel.by_priority():
        style = {"predicted_only": Fmt.green, "contradicted": Fmt.red,
                 "corroborated": Fmt.cyan}.get(r.state.value, Fmt.dim)
        print(style(Fmt.bold(f"{r.variant_id}")) +
              f"  [{r.state.value}]  priority {r.priority:.2f}")
        print(f"  {r.narrative}")
        if r.opportunity:
            print(Fmt.blue(f"  -> {r.opportunity}"))
        for c in r.caveats[:3]:
            print(Fmt.dim(f"  ! {c}"))
        print()
    return 0


def cmd_report(args) -> int:
    store = Store(args.db)
    an = Analyst(store)
    rep = an.corpus_report()

    if args.json:
        print(json.dumps(rep, indent=2, default=str))
        return 0

    header("Knowledge base report")
    print(f"  sources          {rep['sources']}")
    print(f"  documents        {rep['documents']}")
    print(f"  records          {rep['records']}")
    print(f"  distinct labs    {rep['distinct_labs']}")
    print(f"  splicing events  {rep['splicing_events']}")
    print(f"  retracted        {rep['retracted']}")
    print()
    q = rep["quantified_fraction"]
    n = rep["negative_result_fraction"]
    u = rep["unmapped_fraction"]
    print(f"  quantified       {q:.0%}" +
          ("" if q > 0.4 else Fmt.yellow("   (low: many records lack effect sizes)")))
    print(f"  null results     {n:.0%}" +
          ("" if n > 0.05 else Fmt.yellow("   (low: publication bias likely)")))
    print(f"  unmapped pheno   {u:.0%}" +
          ("" if u < 0.25 else Fmt.yellow("   (high: extend the vocabulary)")))

    print()
    print(Fmt.bold("Most frequent phenotypes"))
    for k, c in rep["top_phenotypes"][:10]:
        print(f"  {c:5d}  {k}")
    print()
    print(Fmt.bold("Most frequent pathways"))
    for k, c in rep["top_pathways"][:10]:
        print(f"  {c:5d}  {k}")
    print()
    print(Fmt.bold("Most prolific groups") +
          Fmt.dim("  (check these for over-representation)"))
    for k, c in rep["most_prolific_labs"][:6]:
        print(f"  {c:5d}  {k}")
    return 0


def cmd_ui(args) -> int:
    """Launch the web interface."""
    from ..ui.server import create_app

    if not Path(args.db).exists():
        print(Fmt.red(f"No knowledge base at {args.db}"))
        print(Fmt.dim("Run `synaptik demo` or `synaptik harvest` first."))
        return 1
    app = create_app(args.db, enable_drugs=args.drugs, contact_email=args.email)
    print(f"SYNAPTIK at {Fmt.bold(f'http://{args.host}:{args.port}')}")
    print(Fmt.dim(f"knowledge base: {args.db}"))
    if args.host != "127.0.0.1":
        print(Fmt.red("WARNING: bound to a non-local address with no "
                      "authentication."))
    app.run(host=args.host, port=args.port, debug=False)
    return 0


def cmd_demo(args) -> int:
    """Build the demo knowledge base."""
    # tests/ is not part of the installed package, so locate it relative to
    # the source tree and fail with a useful message if it is absent.
    repo_root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(repo_root))
    try:
        from tests.demo_data import build_demo_store
    except ImportError:
        print(Fmt.red("Demo data builder not found."))
        print(Fmt.dim(f"Looked in {repo_root / 'tests'}. The demo requires the "
                      "source tree;\nrun it from a clone rather than an "
                      "installed wheel."))
        return 1

    print(f"Building demo knowledge base at {args.db} ...")
    store = build_demo_store(args.db)
    for k, v in store.stats().items():
        print(f"  {k}: {v}")
    print()
    print("Try:")
    print(f"  synaptik phenotype 'SYN1 puncta density' --db {args.db} -v")
    print(f"  synaptik splicing --db {args.db}")
    print(f"  synaptik lookup --target TARDBP --assay imaging_synaptic_puncta "
          f"--div 21 --db {args.db}")
    return 0


# --------------------------------------------------------------------------
# argument parsing
# --------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="synaptik",
        description="Literature-grounded experiment lookup and "
                    "phenotype-to-pathway inference for iPSC-neuron research.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  synaptik demo                                    build a demo knowledge base
  synaptik harvest --email you@mgh.harvard.edu     harvest the corpus
  synaptik extract --dry-run                       estimate extraction cost
  synaptik phenotype "reduced SYN1 puncta" -v      analyze a phenotype
  synaptik lookup --target SCN2A --div 21          check if it has been done
  synaptik splicing --gene STMN2                   splicing-focused view
  synaptik compare "SYN1 puncta" "mEPSC frequency" cross-check two readouts
""",
    )
    p.add_argument("--db", default=DEFAULT_DB,
                   help=f"database path (default {DEFAULT_DB})")
    p.add_argument("--verbose-log", action="store_true",
                   help="enable debug logging")

    # Every subcommand also accepts --db, so `synaptik demo --db x` works as
    # naturally as `synaptik --db x demo`. Requiring the global position is a
    # usability trap that produces a confusing argparse error.
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--db", default=None, help=argparse.SUPPRESS)
    common.add_argument("--verbose-log", action="store_true",
                        help=argparse.SUPPRESS)

    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("init", parents=[common], help="create an empty knowledge base")
    s.set_defaults(func=cmd_init)

    s = sub.add_parser("demo", parents=[common], help="build a demo knowledge base with fake data")
    s.set_defaults(func=cmd_demo)

    s = sub.add_parser("preflight",
                       help="check that the APIs still return what we parse")
    s.add_argument("--email")
    s.add_argument("--ncbi-key")
    s.set_defaults(func=cmd_preflight)

    s = sub.add_parser("contracts", parents=[common], help="show or test the API field contracts")
    s.add_argument("--export", metavar="PATH", help="write contracts to JSON")
    s.add_argument("--self-test", action="store_true",
                   help="validate contracts against offline fixtures")
    s.set_defaults(func=cmd_contracts)

    s = sub.add_parser("harvest", parents=[common], help="fetch paper metadata and full text")
    s.add_argument("--skip-preflight", action="store_true",
                   help="skip the API contract check (not recommended)")
    s.add_argument("--allow-drift", action="store_true",
                   help="warn instead of aborting when a contract is violated")
    s.add_argument("--email", help="contact email (required by NCBI/EuropePMC etiquette)")
    s.add_argument("--ncbi-key", help="NCBI API key (raises the rate limit 3x)")
    s.add_argument("--scimago", help="path to a SCImago journal rank CSV")
    s.add_argument("--query", help="custom Europe PMC query (overrides the defaults)")
    s.add_argument("--only", nargs="+", help=f"run only these query sets: "
                                             f"{', '.join(CORPUS_QUERIES)}")
    s.add_argument("--max-per-query", type=int, default=2000)
    s.add_argument("--fulltext", action="store_true", help="also fetch full text")
    s.add_argument("--fulltext-limit", type=int)
    s.set_defaults(func=cmd_harvest)

    s = sub.add_parser("extract", parents=[common], help="run LLM extraction over fetched papers")
    s.add_argument("--model", default="claude-sonnet-4-6")
    s.add_argument("--limit", type=int)
    s.add_argument("--scimago")
    s.add_argument("--dry-run", action="store_true", help="show the queue and cost estimate")
    s.set_defaults(func=cmd_extract)

    s = sub.add_parser("renormalize", parents=[common],
                       help="re-apply the phenotype vocabulary (free, no LLM)")
    s.add_argument("--scimago")
    s.set_defaults(func=cmd_renormalize)

    s = sub.add_parser("lookup", parents=[common], help="has this experiment been done?")
    s.add_argument("--model", choices=[c.value for c in CellModelClass])
    s.add_argument("--div", type=int, help="days in vitro")
    s.add_argument("--perturbation", choices=[c.value for c in PerturbationType])
    s.add_argument("--target", help="gene or protein perturbed")
    s.add_argument("--variant", help="HGVS variant string")
    s.add_argument("--assay", choices=[c.value for c in AssayClass])
    s.add_argument("--markers", nargs="+")
    s.add_argument("--phenotype")
    s.add_argument("--limit", type=int, default=60)
    s.add_argument("--per-tier", type=int, default=6)
    s.add_argument("--json", action="store_true")
    s.set_defaults(func=cmd_lookup)

    s = sub.add_parser("phenotype", parents=[common], help="what does this phenotype mean?")
    s.add_argument("phenotype")
    s.add_argument("--direction", choices=[d.value for d in Direction])
    s.add_argument("--limit", type=int, default=15)
    s.add_argument("-v", "--verbose", action="store_true",
                   help="show evidence structure and caveats per pathway")
    s.add_argument("--drugs", action="store_true", help="query ChEMBL for compounds")
    s.add_argument("--drug-limit", type=int, default=10)
    s.add_argument("--email")
    s.add_argument("--json", action="store_true")
    s.set_defaults(func=cmd_phenotype)

    s = sub.add_parser("splicing", parents=[common], help="splicing events, regulators, ASO strategies")
    s.add_argument("--gene")
    s.add_argument("--regulator")
    s.add_argument("--limit", type=int, default=15)
    s.add_argument("--json", action="store_true")
    s.set_defaults(func=cmd_splicing)

    s = sub.add_parser("compare", parents=[common], help="do two phenotypes implicate the same biology?")
    s.add_argument("phenotype_a")
    s.add_argument("phenotype_b")
    s.add_argument("--json", action="store_true")
    s.set_defaults(func=cmd_compare)

    s = sub.add_parser("gene", parents=[common],
                       help="multi-class evidence dossier for a gene")
    s.add_argument("gene")
    s.add_argument("--offline", action="store_true",
                   help="skip external APIs; literature evidence only")
    s.add_argument("--email")
    s.add_argument("--json", action="store_true")
    s.set_defaults(func=cmd_gene)

    s = sub.add_parser("variant", parents=[common],
                       help="reconcile variants against the literature")
    s.add_argument("variants", nargs="+", metavar="chr:pos:REF>ALT")
    s.add_argument("--gene", help="gene symbol for literature matching")
    s.add_argument("--spliceai-vcf", help="precomputed SpliceAI score VCF")
    s.add_argument("--reference", help="reference genome FASTA")
    s.add_argument("--json", action="store_true")
    s.set_defaults(func=cmd_variant)

    s = sub.add_parser("ui", parents=[common], help="launch the web interface")
    s.add_argument("--host", default="127.0.0.1",
                   help="bind address; the default keeps this local-only")
    s.add_argument("--port", type=int, default=8420)
    s.add_argument("--drugs", action="store_true",
                   help="enable ChEMBL lookups (requires network)")
    s.add_argument("--email")
    s.set_defaults(func=cmd_ui)

    # Newer commands live in cli/commands/ (one module each).
    from .commands import register_all
    register_all(sub, common)

    s = sub.add_parser("report", parents=[common], help="knowledge base health check")
    s.add_argument("--json", action="store_true")
    s.set_defaults(func=cmd_report)

    return p


def main(argv=None) -> int:
    parser = build_parser()
    argv_list = list(argv) if argv is not None else sys.argv[1:]
    args = parser.parse_args(argv_list)

    # argparse lets the subparser's default (None) clobber a --db given before
    # the subcommand. Recover the user's intent from the raw argv: if --db
    # appears anywhere, that value wins; otherwise fall back to the default.
    if getattr(args, "db", None) is None:
        explicit = None
        for i, tok in enumerate(argv_list):
            if tok == "--db" and i + 1 < len(argv_list):
                explicit = argv_list[i + 1]
            elif tok.startswith("--db="):
                explicit = tok.split("=", 1)[1]
        args.db = explicit or DEFAULT_DB

    logging.basicConfig(
        level=logging.DEBUG if args.verbose_log else logging.WARNING,
        format="%(levelname)s %(name)s: %(message)s",
    )
    try:
        return args.func(args)
    except KeyboardInterrupt:
        print("\nInterrupted.")
        return 130
    except Exception as e:
        if args.verbose_log:
            raise
        print(Fmt.red(f"Error: {e}"))
        print(Fmt.dim("Run with --verbose-log for a traceback."))
        return 1


if __name__ == "__main__":
    sys.exit(main())
