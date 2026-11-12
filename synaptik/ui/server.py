"""
Web server.

Read-only over the knowledge base. Bind to localhost by default -- this holds
an institution's literature analysis and has no authentication layer, so
exposing it to a network is a deliberate act requiring --host.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

from flask import Flask, jsonify, request, send_from_directory

from ..analyst import Analyst
from ..drug.repurposing import ChEMBLClient, OpenTargetsClient, RepurposingEngine
from ..infer.lookup import MatchTier
from ..infer.weighting import EvidenceWeighter
from ..ingest.http import PoliteSession
from ..schema import (
    AssayClass, CellModelClass, Direction, ExperimentQuery, PerturbationType,
)
from ..store import Store

log = logging.getLogger(__name__)
STATIC_DIR = Path(__file__).parent / "static"


def create_app(db_path: str, enable_drugs: bool = False,
               contact_email: str | None = None) -> Flask:
    app = Flask(__name__, static_folder=str(STATIC_DIR))
    app.config["DB_PATH"] = db_path
    app.config["ENABLE_DRUGS"] = enable_drugs
    app.config["EMAIL"] = contact_email

    def get_analyst() -> Analyst:
        """One Store per request: SQLite connections are not thread-safe."""
        store = Store(app.config["DB_PATH"])
        repurposer = None
        if app.config["ENABLE_DRUGS"]:
            session = PoliteSession(
                store=store, contact_email=app.config["EMAIL"]
            )
            repurposer = RepurposingEngine(
                chembl=ChEMBLClient(session),
                opentargets=OpenTargetsClient(session),
            )
        return Analyst(store, EvidenceWeighter(), repurposing=repurposer)

    # ---------------- static ----------------

    @app.route("/")
    def index():
        return send_from_directory(STATIC_DIR, "index.html")

    @app.route("/<path:filename>")
    def static_files(filename):
        return send_from_directory(STATIC_DIR, filename)

    # ---------------- api ----------------

    @app.route("/api/stats")
    def api_stats():
        an = get_analyst()
        try:
            return jsonify(an.corpus_report())
        finally:
            an.store.close()

    @app.route("/api/vocab")
    def api_vocab():
        """Vocabulary the UI needs to build its form controls."""
        from ..normalize.phenotype_ontology import TERMS

        return jsonify({
            "phenotypes": [
                {"id": t.term_id, "label": t.label,
                 "synonyms": list(t.synonyms)[:5], "axis": t.measurement_axis}
                for t in TERMS if t.measurement_axis != "composite"
            ],
            "assays": [c.value for c in AssayClass],
            "models": [c.value for c in CellModelClass],
            "perturbations": [c.value for c in PerturbationType],
            "directions": [d.value for d in Direction],
        })

    @app.route("/api/phenotype", methods=["POST"])
    def api_phenotype():
        body = request.get_json(force=True, silent=True) or {}
        phenotype = (body.get("phenotype") or "").strip()
        if not phenotype:
            return jsonify({"error": "phenotype is required"}), 400

        direction = Direction.UNSPECIFIED
        if body.get("direction"):
            try:
                direction = Direction(body["direction"])
            except ValueError:
                pass

        an = get_analyst()
        try:
            res = an.analyze_phenotype(
                phenotype, direction=direction,
                include_drugs=bool(body.get("include_drugs"))
                and app.config["ENABLE_DRUGS"],
            )
            return jsonify(_serialize_phenotype(res))
        finally:
            an.store.close()

    @app.route("/api/lookup", methods=["POST"])
    def api_lookup():
        body = request.get_json(force=True, silent=True) or {}

        def enum_or_none(cls, key):
            v = body.get(key)
            if not v:
                return None
            try:
                return cls(v)
            except ValueError:
                return None

        q = ExperimentQuery(
            cell_model_class=enum_or_none(CellModelClass, "model"),
            days_in_vitro=_int_or_none(body.get("div")),
            perturbation_type=enum_or_none(PerturbationType, "perturbation"),
            perturbation_target=(body.get("target") or "").strip() or None,
            variant_hgvs=(body.get("variant") or "").strip() or None,
            assay_class=enum_or_none(AssayClass, "assay"),
            markers=[m.strip() for m in (body.get("markers") or []) if m.strip()],
            phenotype=(body.get("phenotype") or "").strip() or None,
        )
        if not any([q.cell_model_class, q.perturbation_target, q.assay_class,
                    q.phenotype, q.variant_hgvs]):
            return jsonify({"error": "specify at least one search dimension"}), 400

        an = get_analyst()
        try:
            res = an.lookup_experiment(q, max_results=80)
            return jsonify({
                "verdict": res.verdict,
                "consensus": res.consensus,
                "conflicting": res.conflicting,
                "n_exact": res.n_exact,
                "n_close": res.n_close,
                "tiers": {
                    tier.value: [_serialize_match(m) for m in res.by_tier(tier)]
                    for tier in MatchTier
                },
            })
        finally:
            an.store.close()

    @app.route("/api/splicing", methods=["POST"])
    def api_splicing():
        body = request.get_json(force=True, silent=True) or {}
        an = get_analyst()
        try:
            res = an.analyze_splicing(
                gene=(body.get("gene") or "").strip() or None,
                regulator=(body.get("regulator") or "").strip() or None,
            )
            return jsonify({
                "n_events": res.get("n_events", 0),
                "genes": res.get("genes", []),
                "events": res.get("events", []),
                "message": res.get("message"),
                "regulators": [{
                    "rbp": r.rbp,
                    "n_events": r.n_supporting_events,
                    "n_papers": r.n_papers,
                    "n_labs": r.n_labs,
                    "consistency": round(r.consistency, 2),
                    "signature": r.signature_match,
                    "genes": r.example_genes,
                    "notes": r.therapeutic_notes,
                    "caveats": r.caveats(),
                    "motifs": list(r.profile.motifs) if r.profile else [],
                    "canonical_effect": r.profile.canonical_effect
                    if r.profile else None,
                    "disease_links": list(r.profile.disease_links)
                    if r.profile else [],
                } for r in res.get("regulators", [])],
            })
        finally:
            an.store.close()

    @app.route("/api/compare", methods=["POST"])
    def api_compare():
        body = request.get_json(force=True, silent=True) or {}
        a = (body.get("phenotype_a") or "").strip()
        b = (body.get("phenotype_b") or "").strip()
        if not a or not b:
            return jsonify({"error": "two phenotypes are required"}), 400
        an = get_analyst()
        try:
            return jsonify(an.compare_phenotypes(a, b))
        finally:
            an.store.close()

    @app.errorhandler(500)
    def handle_500(e):
        log.exception("Server error")
        return jsonify({"error": "Something went wrong on the server. "
                                 "Check the server log for details."}), 500

    return app


# --------------------------------------------------------------------------
# serialization
# --------------------------------------------------------------------------


def _int_or_none(v):
    try:
        return int(v) if v not in (None, "") else None
    except (TypeError, ValueError):
        return None


def _serialize_phenotype(res) -> dict:
    return {
        "phenotype": res.phenotype_query,
        "resolved_id": res.resolved_phenotype,
        "resolved_label": res.resolved_label,
        "n_records": res.n_matching_records,
        "n_papers": res.n_papers,
        "n_labs": res.n_labs,
        "corpus_size": res.corpus_size,
        "interpretation": res.interpretation,
        "warnings": res.warnings,
        "associations": [{
            "pathway": a.pathway_label,
            "pathway_key": a.pathway_key,
            "tier": a.confidence_tier(),
            "rank_score": round(a.rank_score(), 2),
            "enrichment": round(a.enrichment, 2),
            "q_value": a.q_value,
            "p_fisher": a.p_fisher,
            "n_records": a.n_records,
            "n_papers": a.n_papers,
            "n_labs": a.n_labs,
            "n_assay_families": a.n_ortho_groups,
            "direction": a.dominant_direction,
            "direction_counts": a.direction_counts,
            "direction_consistency": round(a.direction_consistency, 2),
            "survives_lolo": a.survives_lolo,
            "lolo_enrichment": (round(a.leave_one_lab_out_enrichment, 2)
                                if a.leave_one_lab_out_enrichment else None),
            "lab_concentration": round(a.lab_concentration, 2),
            "n_negative": a.n_negative_records,
            "median_effect_size": a.median_effect_size,
            "effect_size_iqr": a.effect_size_iqr,
            "n_with_effect_size": a.n_with_effect_size,
            "first_year": a.first_year,
            "last_year": a.last_year,
            "caveats": a.caveats(),
            "supporting": [{
                "source_id": s.source_id,
                "title": s.title,
                "journal": s.journal,
                "year": s.year,
                "direction": s.direction,
                "effect_size": s.effect_size,
                "p_value": s.p_value,
                "assay": s.assay_class,
                "weight": s.weight,
                "quote": s.quote,
                "is_negative": s.is_negative,
            } for s in a.supporting[:8]],
        } for a in res.associations],
        "regulators": [{
            "rbp": r.rbp,
            "n_events": r.n_supporting_events,
            "n_labs": r.n_labs,
            "signature": r.signature_match,
            "genes": r.example_genes,
            "notes": r.therapeutic_notes,
            "caveats": r.caveats(),
        } for r in res.regulators],
        "drugs": [{
            "name": d.name,
            "chembl_id": d.chembl_id,
            "target": d.target_symbol,
            "moa": d.moa_direction.value,
            "stage": d.stage_label(),
            "potency_nm": d.potency_nm,
            "potency_type": d.potency_type,
            "cns_summary": d.cns.summary(),
            "cns_notes": d.cns.notes,
            "pathway": d.pathway_label,
            "tier": d.pathway_evidence_tier,
            "rationale": d.rationale,
            "warnings": d.warnings,
            "prior_evidence": d.prior_neuronal_evidence,
        } for d in res.drug_candidates],
    }


def _serialize_match(m) -> dict:
    r = m.record
    p = m.protocol()
    return {
        "source_id": r.source.source_id,
        "pmid": r.source.pmid,
        "doi": r.source.doi,
        "title": r.source.title,
        "journal": r.source.journal,
        "year": r.source.year,
        "quartile": r.source.sjr_quartile,
        "is_retracted": r.source.is_retracted,
        "is_preprint": r.source.is_preprint,
        "first_author": r.source.first_author_surname(),
        "similarity": m.similarity,
        "weight": m.weight,
        "tier": m.tier.value,
        "matched_on": m.matched_on,
        "differs_on": m.differs_on,
        "phenotype": r.outcome.phenotype.raw,
        "direction": r.outcome.direction.value,
        "effect_size": r.outcome.effect_size,
        "effect_unit": r.outcome.effect_size_unit,
        "p_value": r.outcome.p_value,
        "assay": r.assay_class.value,
        "protocol": {
            "differentiation": p.differentiation_protocol,
            "div": p.days_in_vitro,
            "coculture": p.coculture,
            "markers": p.markers,
            "n_replicates": p.n_biological_replicates,
            "statistical_test": p.statistical_test,
            "isogenic": p.is_isogenic_pair,
            "control": p.control_description,
            "dose": p.dose,
            "duration": p.duration,
        },
        "quote": r.provenance.quote,
    }


def main() -> int:
    import argparse

    ap = argparse.ArgumentParser(description="SYNAPTIK web interface")
    ap.add_argument("--db", default=os.environ.get("SYNAPTIK_DB",
                                                   "data/synaptik.db"))
    ap.add_argument("--host", default="127.0.0.1",
                    help="bind address; the default keeps this local-only")
    ap.add_argument("--port", type=int, default=8420)
    ap.add_argument("--drugs", action="store_true",
                    help="enable ChEMBL lookups (requires network)")
    ap.add_argument("--email", help="contact email for API etiquette")
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO)
    app = create_app(args.db, enable_drugs=args.drugs, contact_email=args.email)
    print(f"SYNAPTIK running at http://{args.host}:{args.port}")
    print(f"knowledge base: {args.db}")
    if args.host != "127.0.0.1":
        print("WARNING: bound to a non-local address with no authentication.")
    app.run(host=args.host, port=args.port, debug=args.debug)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
