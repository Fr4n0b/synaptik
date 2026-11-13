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

from flask import Flask, Response, jsonify, request, send_from_directory

from ..analyst import Analyst
from ..drug.repurposing import ChEMBLClient, OpenTargetsClient, RepurposingEngine
from ..infer.lookup import MatchTier
from ..infer.weighting import EvidenceWeighter
from ..ingest.http import PoliteSession
from ..schema import (
    AssayClass, CellModelClass, Direction, ExperimentQuery, PerturbationType,
)
from ..labdata import LabDataError, parse_splicing_table
from ..store import Store
from ..workflows import (
    benchmark_to_dict, deposit_phenotypes, deposit_splicing, design_brief,
    design_brief_to_dict, phenotype_report_to_dict, run_benchmark,
    splicing_report_to_dict,
)

log = logging.getLogger(__name__)
STATIC_DIR = Path(__file__).parent / "static"


def create_app(db_path: str, enable_drugs: bool = False,
               contact_email: str | None = None) -> Flask:
    app = Flask(__name__, static_folder=str(STATIC_DIR))
    app.config["DB_PATH"] = db_path
    # Uploaded tables are parsed in memory; cap them so a stray BAM file
    # cannot exhaust the server.
    app.config["MAX_CONTENT_LENGTH"] = 64 * 1024 * 1024
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

    def query_from(body: dict) -> ExperimentQuery:
        def enum_or_none(cls, key):
            v = body.get(key)
            if not v:
                return None
            try:
                return cls(v)
            except ValueError:
                return None

        return ExperimentQuery(
            cell_model_class=enum_or_none(CellModelClass, "model"),
            days_in_vitro=_int_or_none(body.get("div")),
            perturbation_type=enum_or_none(PerturbationType, "perturbation"),
            perturbation_target=(body.get("target") or "").strip() or None,
            variant_hgvs=(body.get("variant") or "").strip() or None,
            assay_class=enum_or_none(AssayClass, "assay"),
            markers=[m.strip() for m in (body.get("markers") or [])
                     if m.strip()],
            phenotype=(body.get("phenotype") or "").strip() or None,
        )

    @app.route("/api/lookup", methods=["POST"])
    def api_lookup():
        body = request.get_json(force=True, silent=True) or {}

        q = query_from(body)
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

    @app.route("/api/design", methods=["POST"])
    def api_design():
        body = request.get_json(force=True, silent=True) or {}
        q = query_from(body)
        if not any([q.cell_model_class, q.perturbation_target, q.assay_class,
                    q.phenotype, q.variant_hgvs]):
            return jsonify({"error": "describe the experiment first: give at "
                                     "least a target, assay or readout"}), 400
        cv = None
        if body.get("cv") not in (None, ""):
            try:
                cv = float(body["cv"])
            except (TypeError, ValueError):
                return jsonify({"error": "CV must be a number"}), 400
            if cv > 2:      # entered as a percentage
                cv = cv / 100.0
            if not 0 < cv < 2:
                return jsonify({"error": "CV must be between 0 and 200%"}), 400
        store = Store(app.config["DB_PATH"])
        try:
            return jsonify(design_brief_to_dict(design_brief(store, q, cv=cv)))
        finally:
            store.close()

    def _uploaded_text():
        """Decode an uploaded table in memory. Nothing is written to disk."""
        f = request.files.get("file")
        if f is None or not f.filename:
            return None, None
        raw = f.read()
        return f.filename, raw.decode("utf-8-sig", errors="replace")

    @app.route("/api/deposit/splicing", methods=["POST"])
    def api_deposit_splicing():
        name, text = _uploaded_text()
        if text is None:
            return jsonify({"error": "choose a file to upload"}), 400
        form = request.form
        try:
            max_fdr = float(form.get("fdr") or 0.05)
            min_dpsi = float(form.get("min_dpsi") or 0.1)
        except ValueError:
            return jsonify({"error": "thresholds must be numbers"}), 400
        try:
            table = parse_splicing_table(
                name, perturbed=form.get("perturbed") or "",
                event_type=form.get("event_type") or None, text=text)
        except LabDataError as e:
            return jsonify({"error": str(e)}), 400
        store = Store(app.config["DB_PATH"])
        try:
            rep = deposit_splicing(store, table, max_fdr, min_dpsi)
            return jsonify(splicing_report_to_dict(rep))
        finally:
            store.close()

    @app.route("/api/deposit/phenotypes", methods=["POST"])
    def api_deposit_phenotypes():
        name, text = _uploaded_text()
        if text is None:
            return jsonify({"error": "choose a file to upload"}), 400
        store = Store(app.config["DB_PATH"])
        try:
            rows, warnings = deposit_phenotypes(store, name, text=text)
            return jsonify(phenotype_report_to_dict(rows, warnings))
        except LabDataError as e:
            return jsonify({"error": str(e)}), 400
        finally:
            store.close()

    @app.route("/api/benchmark", methods=["POST"])
    def api_benchmark():
        body = request.get_json(force=True, silent=True) or {}
        store = Store(app.config["DB_PATH"])
        try:
            cutoff = _int_or_none(body.get("cutoff"))
            if cutoff is None:
                years = [r["year"] for r in store.conn.execute(
                    "SELECT year FROM records WHERE year IS NOT NULL "
                    "ORDER BY year")]
                if not years:
                    return jsonify({"error": "no dated records"}), 400
                cutoff = years[int(len(years) * 0.7)]
            return jsonify(benchmark_to_dict(run_benchmark(store, cutoff)))
        finally:
            store.close()

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

    @app.route("/api/nof1", methods=["POST"])
    def api_nof1():
        body = request.get_json(force=True, silent=True) or {}
        gene = (body.get("gene") or "").strip().upper()
        if not gene:
            return jsonify({"error": "target gene symbol is required"}), 400
        patient = (body.get("patient") or "PATIENT-N1").strip()
        disease = (body.get("disease") or "als").strip().lower()
        variant = (body.get("variant") or "cryptic_splice_variant").strip()
        seq = body.get("sequence") or None
        nfl = float(body.get("nfl") or 65.0)
        div = int(body.get("div") or 28)

        from ..workflows import run_nof1_workflow
        dossier_dict, dossier_md = run_nof1_workflow(
            patient_id=patient,
            target_gene=gene,
            disease=disease,
            variant_id=variant,
            target_sequence=seq,
            baseline_nfl=nfl,
            div=div,
        )
        return jsonify({"dossier": dossier_dict, "markdown": dossier_md})

    @app.route("/api/aso", methods=["POST"])
    def api_aso():
        body = request.get_json(force=True, silent=True) or {}
        seq = (body.get("sequence") or "").strip().upper()
        if not seq:
            return jsonify({"error": "RNA target sequence is required"}), 400
        gene = (body.get("gene") or "TARGET").strip().upper()
        mode = (body.get("mode") or "splice_switch").strip()
        top_k = int(body.get("top_k") or body.get("length") or 5)
        from ..workflows import aso_result_to_dict, design_aso
        res = design_aso(target_sequence=seq, target_gene=gene, mode=mode, top_k=top_k)
        return jsonify(aso_result_to_dict(res))

    @app.route("/api/translate", methods=["POST"])
    def api_translate():
        body = request.get_json(force=True, silent=True) or {}
        genes = body.get("genes") or []
        disease = (body.get("disease") or "als").strip().lower()
        if not genes:
            return jsonify({"error": "list of gene symbols is required"}), 400
        entries = [(str(g).strip().upper(), None, None) for g in genes]
        from ..workflows import translate_findings, translate_result_to_dict
        res = translate_findings(gene_entries=entries, disease=disease)
        return jsonify(translate_result_to_dict(res))

    @app.route("/api/simulate", methods=["POST"])
    def api_simulate():
        body = request.get_json(force=True, silent=True) or {}
        disease = (body.get("disease") or "als").strip().lower()
        if body.get("screen"):
            from ..workflows import screen_counterfactual_targets, simulation_result_to_dict
            results = screen_counterfactual_targets(disease=disease)
            return jsonify({"results": [simulation_result_to_dict(r) for r in results]})
        target = (body.get("target") or "STMN2").strip().upper()
        pert = (body.get("intervention") or "aso_steric_block").strip()
        alpha = float(body.get("alpha") or 0.85)
        from ..workflows import simulate_intervention
        res = simulate_intervention(disease=disease, target_gene=target, intervention=pert, efficacy_alpha=alpha)
        return jsonify(res)

    @app.route("/api/biomarker", methods=["POST"])
    def api_biomarker():
        body = request.get_json(force=True, silent=True) or {}
        disease = (body.get("disease") or "als").strip().lower()
        target = (body.get("target") or "STMN2").strip().upper()
        rescue_pct = float(body.get("rescue_pct") or 75.0)
        baseline_nfl = float(body.get("baseline_nfl") or 65.0)
        from ..workflows import evaluate_biomarker_trajectory
        res = evaluate_biomarker_trajectory(
            target_gene=target,
            disease=disease,
            in_vitro_rescue_pct=rescue_pct,
            baseline_plasma_nfl=baseline_nfl,
        )
        return jsonify(res)

    @app.route("/api/protocol", methods=["POST"])
    def api_protocol():
        body = request.get_json(force=True, silent=True) or {}
        if body.get("sites"):
            from ..workflows import synthesize_multisite_replication
            target = (body.get("target") or "STMN2").strip().upper()
            endpoint = (body.get("endpoint") or "Primary phenotypic rescue").strip()
            return jsonify(synthesize_multisite_replication(body["sites"], target_gene=target, endpoint=endpoint))
        target = (body.get("target") or "TARDBP").strip().upper()
        disease = (body.get("disease") or "als").strip().lower()
        cell_model = (body.get("cell_model") or "ipsc_motor_neuron").strip()
        div = int(body.get("div") or 28)
        pairs = int(body.get("pairs") or 8)
        rho = float(body.get("rho") or 0.75)
        from ..workflows import generate_study_protocol
        return jsonify(generate_study_protocol(target_gene=target, disease=disease, cell_model=cell_model, maturity_div=div, paired_n=pairs, rho=rho))

    @app.route("/api/vcf", methods=["POST"])
    def api_vcf():
        # Supports JSON payload { "vcf_text": "...", "target": "..." } or multipart form upload
        body = request.get_json(force=True, silent=True) or {}
        vcf_text = body.get("vcf_text")
        target_gene = body.get("target") or request.form.get("target") or None

        if not vcf_text:
            name, text = _uploaded_text()
            vcf_text = text

        if not vcf_text:
            return jsonify({"error": "Provide VCF text or upload a VCF file"}), 400

        from ..workflows import run_vcf_prioritization_workflow
        res = run_vcf_prioritization_workflow(vcf_text, target_gene=target_gene)
        return jsonify(res)

    @app.route("/api/repurpose", methods=["POST"])
    def api_repurpose():
        body = request.get_json(force=True, silent=True) or {}
        target = (body.get("target") or "STMN2").strip().upper()
        direction = (body.get("direction") or "activator").strip().lower()
        min_phase = int(body.get("min_phase") or 0)
        bbb_only = bool(body.get("bbb_only") or False)

        from ..workflows import run_drug_repurposing_workflow
        cands = run_drug_repurposing_workflow(
            target_gene=target,
            required_direction=direction,
            min_clinical_phase=min_phase,
            bbb_only=bbb_only,
        )
        return jsonify({"target": target, "direction": direction, "candidates": cands})

    @app.route("/api/visualize/splicing", methods=["GET", "POST"])
    def api_visualize_splicing():
        body = request.get_json(force=True, silent=True) or {}
        gene = (request.args.get("gene") or body.get("gene") or "STMN2").strip().upper()
        lead_aso = request.args.get("aso") or body.get("lead_aso_sequence") or None
        delta_psi = float(request.args.get("delta_psi") or body.get("delta_psi") or 85.0)

        from ..workflows import generate_splicing_svg_workflow
        svg_xml = generate_splicing_svg_workflow(gene_symbol=gene, lead_aso_sequence=lead_aso, delta_psi_pct=delta_psi)
        return Response(svg_xml, mimetype="image/svg+xml")

    @app.route("/api/nof1/preind", methods=["POST"])
    def api_nof1_preind():
        body = request.get_json(force=True, silent=True) or {}
        gene = (body.get("gene") or "STMN2").strip().upper()
        patient = (body.get("patient") or "PATIENT-N1").strip()
        disease = (body.get("disease") or "als").strip().lower()
        variant = (body.get("variant") or "cryptic_splice_variant").strip()
        seq = body.get("sequence") or None
        nfl = float(body.get("nfl") or 65.0)
        div = int(body.get("div") or 28)

        from ..workflows import run_preind_workflow
        preind_md = run_preind_workflow(
            patient_id=patient,
            target_gene=gene,
            disease=disease,
            variant_id=variant,
            target_sequence=seq,
            baseline_nfl=nfl,
            div=div,
        )
        return jsonify({"gene": gene, "patient": patient, "markdown": preind_md})

    @app.route("/api/cite", methods=["GET", "POST"])
    def api_cite():
        from ..cli.commands.cite import AFFILIATION, AUTHOR, ORCID, ORCID_URL, REPO_URL, get_apa, get_bibtex, get_nature
        from .. import __version__ as version
        return jsonify({
            "author": AUTHOR,
            "orcid": ORCID,
            "orcid_url": ORCID_URL,
            "affiliation": AFFILIATION,
            "version": version,
            "repository": REPO_URL,
            "citations": {
                "bibtex": get_bibtex(version),
                "apa": get_apa(version),
                "nature": get_nature(version),
            },
        })

    # ---------------- clinical endpoints ----------------

    @app.route("/api/clinical/allometry", methods=["POST"])
    def api_clinical_allometry():
        body = request.get_json(force=True, silent=True)
        if body is None:
            return jsonify({"error": "Missing JSON request body"}), 400

        noael_mg = body.get("noael_mg")
        noael_mg_kg = body.get("noael_mg_kg")
        if noael_mg is None and noael_mg_kg is None:
            noael_mg = 30.0
        try:
            noael_val = float(noael_mg) if noael_mg is not None else None
            noael_kg_val = float(noael_mg_kg) if noael_mg_kg is not None else None
            nhp_wt = float(body.get("nhp_weight_kg", 3.0) or 3.0)
            nhp_vol = float(body.get("nhp_volume_ml", 15.0) or 15.0)
            adult_vol = float(body.get("adult_volume_ml", 140.0) or 140.0)
            pedia_vol = float(body.get("pediatric_volume_ml", 65.0) or 65.0)
            sf = float(body.get("safety_factor", 10.0) or 10.0)
        except (TypeError, ValueError) as err:
            return jsonify({"error": f"Invalid numerical parameter: {err}"}), 400

        from ..clinical.allometry import calculate_allometric_dose
        try:
            adult_res = calculate_allometric_dose(
                noael_mg=noael_val,
                noael_mg_kg=noael_kg_val,
                nhp_weight_kg=nhp_wt,
                nhp_volume_ml=nhp_vol,
                human_volume_ml=adult_vol,
                safety_factor=sf,
                cohort_name="adult",
                acknowledge_not_clinical=True,
            )
            pedia_res = calculate_allometric_dose(
                noael_mg=noael_val,
                noael_mg_kg=noael_kg_val,
                nhp_weight_kg=nhp_wt,
                nhp_volume_ml=nhp_vol,
                human_volume_ml=pedia_vol,
                safety_factor=sf,
                cohort_name="pediatric",
                acknowledge_not_clinical=True,
            )
        except ValueError as err:
            return jsonify({"error": str(err)}), 400

        return jsonify({
            "disclaimer": "RESEARCH PLANNING ONLY — NOT A CLINICAL DOSE RECOMMENDATION",
            "is_clinical_recommendation": False,
            "author": "Francesco Nobile (ORCID: 0009-0001-9477-3332)",
            "noael_mg": adult_res.noael_nhp_mg,
            "nhp_csf_volume_ml": nhp_vol,
            "safety_factor": sf,
            "adult": {
                "csf_volume_ml": adult_res.human_volume_ml,
                "scaling_ratio": round(adult_res.scaling_ratio, 4),
                "hed_mg": round(adult_res.hed_mg, 4),
                "mrsd_mg": round(adult_res.mrsd_mg, 4),
            },
            "pediatric": {
                "csf_volume_ml": pedia_res.human_volume_ml,
                "scaling_ratio": round(pedia_res.scaling_ratio, 4),
                "hed_mg": round(pedia_res.hed_mg, 4),
                "mrsd_mg": round(pedia_res.mrsd_mg, 4),
            },
            "methodology": "CSF volume allometric scaling (Preclinical Research Planning / FDA 2005 Guidance)",
        })

    @app.route("/api/clinical/pkpd", methods=["POST"])
    def api_clinical_pkpd():
        body = request.get_json(force=True, silent=True)
        if body is None:
            return jsonify({"error": "Missing JSON request body"}), 400

        try:
            dose_mg = float(body.get("dose_mg", 50.0) or 50.0)
            route = str(body.get("route", "bolus") or "bolus").strip().lower()
            csf_t12 = float(body.get("csf_half_life_days", 24.0) or 24.0)
            tis_t12 = float(body.get("tissue_half_life_days", 135.0) or 135.0)
            days = int(body.get("duration_days", 360) or 360)
            infusion_h = float(body.get("infusion_hours", 24.0) or 24.0)
        except (TypeError, ValueError) as err:
            return jsonify({"error": f"Invalid numerical parameter: {err}"}), 400

        if route not in ("bolus", "catheter"):
            return jsonify({"error": "Unsupported administration route. Must be 'bolus' or 'catheter'"}), 400

        from ..clinical.pkpd import simulate_intrathecal_pkpd
        try:
            sim = simulate_intrathecal_pkpd(
                dose_mg=dose_mg,
                route=route,
                csf_half_life_days=csf_t12,
                tissue_half_life_days=tis_t12,
                duration_days=days,
                infusion_hours=infusion_h,
            )
        except ValueError as err:
            return jsonify({"error": str(err)}), 400

        return jsonify({
            "author": "Francesco Nobile (ORCID: 0009-0001-9477-3332)",
            "dose_mg": sim.dose_mg,
            "route": sim.route,
            "kinetics": sim.kinetics.to_dict(),
            "exposure": sim.exposure.to_dict(),
            "pharmacodynamics": sim.pharmacodynamics.to_dict(),
            "timepoints": [
                {
                    "day": p.day,
                    "csf_conc_ug_ml": round(p.c_csf_ug_ml, 2),
                    "tissue_conc_ug_g": round(p.c_tissue_ug_g, 2),
                    "pd_effect_pct": round(p.target_engagement_pct, 1),
                    "nfl_reduction_pct": round(p.nfl_reduction_pct, 1),
                }
                for p in sim.trajectory
                if p.day in {0, 1, 7, 14, 28, 60, 90, 180, 270, 360}
            ],
        })

    @app.route("/api/clinical/escalate", methods=["POST"])
    def api_clinical_escalate():
        body = request.get_json(force=True, silent=True)
        if body is None:
            return jsonify({"error": "Missing JSON request body"}), 400

        design = str(body.get("design", "boin") or "boin").strip().lower()
        if design not in ("boin", "3+3", "mtpi"):
            return jsonify({"error": f"Invalid escalation design: {design}. Must be 'boin', '3+3', or 'mtpi'"}), 400

        try:
            target_phi = float(body.get("target_phi", 0.25) or 0.25)
            dose_levels = [float(x) for x in (body.get("dose_levels") or [10.0, 25.0, 50.0, 100.0, 150.0])]
            current_level = int(body.get("current_level", 1) or 1)
        except (TypeError, ValueError) as err:
            return jsonify({"error": f"Invalid parameter format: {err}"}), 400

        if target_phi <= 0.0 or target_phi >= 1.0:
            return jsonify({"error": "Target toxicity rate phi must be strictly between 0 and 1"}), 400

        if current_level < 1 or current_level > len(dose_levels):
            return jsonify({"error": f"current_level {current_level} is out of bounds for dose_levels (1..{len(dose_levels)})"}), 400

        from ..clinical.escalation import BOINEscalationEngine, CohortOutcome, MTPIEngine, ThreePlusThreeEngine

        if design == "boin":
            engine = BOINEscalationEngine(target_phi=target_phi, dose_levels=dose_levels)
        elif design == "3+3":
            engine = ThreePlusThreeEngine(dose_levels=dose_levels)
        else:
            engine = MTPIEngine(target_phi=target_phi, dose_levels=dose_levels)

        history_raw = body.get("history") or []
        history_records: list[CohortOutcome] = []
        eval_n = int(body.get("n_patients", 3) or 3)
        eval_dlts = int(body.get("dlts", 0) or 0)

        if isinstance(history_raw, list) and len(history_raw) > 0:
            for item in history_raw:
                try:
                    lvl = int(item.get("dose_level", 1))
                    mg = float(item.get("dose_mg", dose_levels[lvl - 1] if 1 <= lvl <= len(dose_levels) else 10.0))
                    n_pts = int(item.get("n_patients", item.get("n", item.get("n_enrolled", 3))))
                    d_count = int(item.get("dlts", item.get("n_dlts", 0)))
                    history_records.append(CohortOutcome(dose_level=lvl, dose_mg=mg, n_enrolled=n_pts, n_dlts=d_count))
                except Exception:
                    pass
            if history_records:
                eval_level = history_records[-1].dose_level
                eval_n = history_records[-1].n_enrolled
                eval_dlts = history_records[-1].n_dlts
                current_level = eval_level

        try:
            decision = engine.evaluate_cohort(
                dose_level=current_level,
                n=eval_n,
                dlts=eval_dlts,
                history=history_records,
            )
        except ValueError as err:
            return jsonify({"error": str(err)}), 400

        rec_str = decision.action.value
        if rec_str == "retain":
            rec_str = "stay"

        return jsonify({
            "author": "Francesco Nobile (ORCID: 0009-0001-9477-3332)",
            "design": decision.design,
            "target_phi": decision.target_phi,
            "lambda_e": round(decision.lambda_e, 3),
            "lambda_d": round(decision.lambda_d, 3),
            "recommendation": rec_str,
            "recommended_dose_level": decision.next_dose_level,
            "recommended_dose_mg": decision.next_dose_mg,
            "rationale": decision.rationale,
            "is_trial_stopped": decision.trial_stopped,
            "safety_stopping_triggered": decision.trial_stopped or decision.is_dose_eliminated,
            "dose_evaluations": [
                {
                    "dose_level": h.dose_level,
                    "dose_mg": h.dose_mg,
                    "n_patients": h.n_enrolled,
                    "dlts": h.n_dlts,
                    "observed_dlt_rate": round(h.n_dlts / h.n_enrolled, 3) if h.n_enrolled > 0 else 0.0,
                    "posterior_excess_toxicity_prob": round(decision.posterior_overdose_prob, 3),
                    "eliminated": h.dose_level in decision.eliminated_dose_levels,
                }
                for h in history_records
            ],
        })

    @app.route("/api/clinical/protocol", methods=["POST"])
    def api_clinical_protocol():
        body = request.get_json(force=True, silent=True)
        if body is None:
            return jsonify({"error": "Missing JSON request body"}), 400

        gene = (body.get("gene") or "STMN2").strip().upper()
        disease = (body.get("disease") or "als").strip().lower()
        try:
            noael_mg = float(body.get("noael_mg", 30.0) or 30.0)
            design = str(body.get("design", "boin") or "boin").strip().lower()
            target_phi = float(body.get("target_phi", 0.25) or 0.25)
            doses = [float(x) for x in (body.get("doses") or [10.0, 25.0, 50.0, 100.0, 150.0])]
            route = str(body.get("route", "bolus") or "bolus").strip().lower()
            sf = float(body.get("safety_factor", 10.0) or 10.0)
        except (TypeError, ValueError) as err:
            return jsonify({"error": f"Invalid parameter format: {err}"}), 400

        from ..clinical.protocol import generate_clinical_protocol
        try:
            proto = generate_clinical_protocol(
                target_gene=gene,
                disease=disease,
                noael_nhp_mg=noael_mg,
                design=design,
                target_phi=target_phi,
                planned_doses=doses,
                administration_route=route,
                safety_factor=sf,
            )
        except Exception as err:
            return jsonify({"error": str(err)}), 400

        return jsonify(proto.to_dict())

    @app.route("/api/biophysics", methods=["POST"])
    def api_biophysics():
        body = request.get_json(force=True, silent=True)
        if body is None:
            return jsonify({"error": "Missing JSON request body"}), 400
        seq = (body.get("sequence") or body.get("target_sequence") or "").strip()
        if not seq:
            return jsonify({"error": "sequence is required"}), 400
        try:
            start = int(body.get("target_start", 0) or 0)
            end = int(body.get("target_end", len(seq)) or len(seq))
            duplex_dg = float(body.get("duplex_delta_g", -22.0) or -22.0)
            competitor = str(body.get("competitor", "u1_snrnp") or "u1_snrnp")
        except (TypeError, ValueError) as err:
            return jsonify({"error": f"Invalid parameter: {err}"}), 400

        from ..biophysics.accessibility import compute_effective_binding
        from ..biophysics.competition import compute_steric_displacement
        from ..biophysics.folding import compute_target_accessibility, fold_rna_mfe

        try:
            mfe, db = fold_rna_mfe(seq)
            acc = compute_target_accessibility(seq, start, end)
            eff = compute_effective_binding(duplex_dg, acc.delta_g_open)
            disp = compute_steric_displacement(max(eff.kd_effective_nm, 1e-6), competitor_name=competitor)
        except ValueError as err:
            return jsonify({"error": str(err)}), 400

        return jsonify({
            "author": "Francesco Nobile (ORCID: 0009-0001-9477-3332)",
            "mfe_kcal_mol": mfe,
            "dot_bracket": db,
            "accessibility": acc.to_dict(),
            "effective_binding": eff.to_dict(),
            "steric_displacement": disp.to_dict(),
        })

    @app.route("/api/twin", methods=["POST"])
    def api_twin():
        body = request.get_json(force=True, silent=True)
        if body is None:
            return jsonify({"error": "Missing JSON request body"}), 400
        try:
            patient_id = str(body.get("patient_id", "EPHEMERAL-PT") or "EPHEMERAL-PT")
            age = float(body.get("age_years", 55.0) or 55.0)
            site = str(body.get("onset_site", "limb") or "limb")
            duration = float(body.get("disease_duration_months", 12.0) or 12.0)
            alsfrs = float(body.get("baseline_alsfrs_r", 38.0) or 38.0)
            nfl = float(body["baseline_plasma_nfl"]) if body.get("baseline_plasma_nfl") is not None else None
            slope = float(body["baseline_slope"]) if body.get("baseline_slope") is not None else None
            months = int(body.get("duration_months", 18) or 18)
            twins = int(body.get("n_twins", 1000) or 1000)
            seed = int(body.get("random_seed", 42)) if body.get("random_seed") is not None else 42
        except (TypeError, ValueError) as err:
            return jsonify({"error": f"Invalid parameter: {err}"}), 400

        from ..twin.bayesian import evaluate_treatment_efficacy
        from ..twin.models import PatientBaseline, VisitObservation
        from ..twin.twin_generator import generate_counterfactual_cohort

        try:
            baseline = PatientBaseline(
                patient_id=patient_id,
                age_years=age,
                onset_site=site,
                disease_duration_months=duration,
                baseline_alsfrs_r=alsfrs,
                baseline_slope=slope,
                baseline_plasma_nfl=nfl,
            )
            cohort = generate_counterfactual_cohort(baseline, n_twins=twins, duration_months=months, random_seed=seed)
            out = {"cohort": cohort.to_dict(), "author": "Francesco Nobile (ORCID: 0009-0001-9477-3332)"}
            raw_visits = body.get("visits")
            if isinstance(raw_visits, list) and len(raw_visits) >= 2:
                visits = [VisitObservation(month=float(v["month"]), alsfrs_score=float(v["alsfrs_score"])) for v in raw_visits]
                eff = evaluate_treatment_efficacy(cohort, visits)
                out["efficacy"] = eff.to_dict()
        except (ValueError, KeyError, TypeError) as err:
            return jsonify({"error": str(err)}), 400

        return jsonify(out)

    @app.errorhandler(413)
    def handle_413(e):
        return jsonify({"error": "That file is larger than 64 MB. Upload the "
                                 "results table, not raw data."}), 413

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
