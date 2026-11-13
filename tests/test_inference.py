"""
Regression tests for the inference engine.

These assert the properties that make the tool trustworthy. If any of these
fail, the tool is producing misleading output and must not be used.

Run: PYTHONPATH=. python3 tests/test_inference.py
"""

from __future__ import annotations

import sys

from synaptik.infer.association import AssociationEngine, phenotype_coherence
from synaptik.infer.weighting import EvidenceWeighter, JournalQualityResolver
from synaptik.normalize.labs import LabResolver, normalize_author
from synaptik.normalize.matcher import (
    PhenotypeMatcher, classify_assay, extract_gene_candidates,
)
from synaptik.schema import AssayClass, Direction

from .synthetic import build_corpus

FAILURES: list[str] = []


def check(cond: bool, msg: str) -> None:
    if cond:
        print(f"  PASS  {msg}")
    else:
        print(f"  FAIL  {msg}")
        FAILURES.append(msg)


def test_associations() -> None:
    print("\n[association engine: ground-truth recovery]")
    fg, bg = build_corpus()
    res = AssociationEngine(EvidenceWeighter(current_year=2026)).associate(fg, bg)
    by = {a.pathway_label: a for a in res}

    check("synaptic vesicle cycle" in by, "true-positive pathway is detected")
    tp = by.get("synaptic vesicle cycle")
    if tp:
        check(res[0].pathway_label == "synaptic vesicle cycle",
              "true positive ranks first")
        check(tp.confidence_tier() == "strong", "true positive tiered 'strong'")
        check(tp.n_labs >= 3, f"true positive has >=3 labs (got {tp.n_labs})")
        check(tp.n_ortho_groups >= 2,
              f"true positive spans >=2 assay families (got {tp.n_ortho_groups})")
        check(tp.survives_lolo, "true positive survives leave-one-lab-out")
        check(tp.caveats() == [], "true positive carries no caveats")
        check(tp.n_negative_records >= 1, "true positive's null result was retained")

    trap = by.get("exotic pathway X")
    check(trap is not None, "single-lab trap is detected")
    if trap:
        check(trap.confidence_tier() == "single-lab",
              "single-lab trap NOT tiered strong/moderate")
        check(trap.n_labs == 1, "single-lab trap correctly counted as one lab")
        check(not trap.survives_lolo, "single-lab trap fails leave-one-lab-out")
        check(any("one research group" in c for c in trap.caveats()),
              "single-lab trap warns about lack of replication")
        check(any("null results" in c for c in trap.caveats()),
              "single-lab trap flags publication bias")
        if tp:
            check(trap.enrichment > tp.enrichment and
                  trap.rank_score() < tp.rank_score(),
                  "higher enrichment does NOT outrank better-replicated evidence")

    apo = by.get("apoptosis")
    check(apo is not None, "ubiquitous pathway present in results")
    if apo:
        check(apo.rank_score() == 0.0, "ubiquitous pathway has zero rank score")
        check(res[-1].pathway_label == "apoptosis", "ubiquitous pathway ranks last")

    contested = by.get("contested pathway Y")
    if contested:
        check(contested.direction_consistency <= 0.6,
              "conflicting directions produce low consistency")
        check(any("split" in c for c in contested.caveats()),
              "conflicting directions raise a caveat")

    fraud = by.get("fraudulent pathway Z")
    if fraud:
        check(fraud.confidence_tier() == "insufficient",
              "retraction-zeroed evidence cannot reach 'moderate'")


def test_weighting() -> None:
    print("\n[evidence weighting]")
    from synaptik.schema import (
        CellModel, CellModelClass, EvidenceType, ExperimentRecord,
        NormalizedTerm, Outcome, Perturbation, Provenance, SourceMetadata,
    )

    def mk(**kw):
        src = SourceMetadata(
            source_id="PMID:1", title="t",
            journal=kw.get("journal", "Neuron"),
            year=kw.get("year", 2023),
            sjr_quartile=kw.get("q", "Q1"),
            is_retracted=kw.get("retracted", False),
            is_preprint=kw.get("preprint", False),
        )
        return ExperimentRecord(
            source=src,
            cell_model=CellModel(model_class=CellModelClass.IPSC_NEURON_NGN2,
                                 is_isogenic_pair=kw.get("isogenic", True)),
            perturbation=Perturbation(),
            outcome=Outcome(phenotype=NormalizedTerm(raw="p"),
                            direction=Direction.DECREASE,
                            n_biological_replicates=kw.get("n", 4),
                            effect_size=kw.get("es", -20.0),
                            is_quantified=kw.get("es", -20.0) is not None),
            evidence_type=kw.get("ev", EvidenceType.DIRECT_MEASUREMENT),
            provenance=Provenance(extraction_confidence=kw.get("conf", 0.9)),
        )

    w = EvidenceWeighter(current_year=2026)
    q1 = w.weight(mk(q="Q1")).composite
    q4 = w.weight(mk(q="Q4")).composite
    check(q1 > q4, f"Q1 outweighs Q4 ({q1:.3f} > {q4:.3f})")
    check(q1 / q4 < 3.0, "journal tier ratio stays modest (not prestige-dominated)")

    check(w.weight(mk(retracted=True)).composite == 0.0,
          "retracted paper weighted to exactly zero")
    check(w.weight(mk(isogenic=True)).composite >
          w.weight(mk(isogenic=False)).composite,
          "isogenic controls increase weight")
    check(w.weight(mk(n=8)).composite > w.weight(mk(n=2)).composite,
          "higher replicate count increases weight")
    check(w.weight(mk(es=None)).composite < w.weight(mk(es=-20.0)).composite,
          "unquantified results weighted lower")
    from synaptik.schema import EvidenceType as ET
    check(w.weight(mk(ev=ET.DIRECT_MEASUREMENT)).composite >
          w.weight(mk(ev=ET.INFERRED)).composite,
          "direct measurement outweighs author inference")
    check(w.weight(mk(year=2024)).composite > w.weight(mk(year=1998)).composite,
          "recent work weighted above very old work")
    check(w.weight(mk(preprint=True)).composite <
          w.weight(mk(preprint=False)).composite,
          "preprints discounted relative to peer-reviewed")


def test_journal_resolution() -> None:
    print("\n[journal quartile resolution]")
    r = JournalQualityResolver()
    check(r.resolve("Neuron")[0] == "Q1", "Neuron resolves to Q1")
    check(r.resolve("NEURON")[0] == "Q1", "resolution is case-insensitive")
    check(r.resolve("Neuroscience Letters")[0] == "Q3",
          "Neuroscience Letters resolves to Q3")
    check(r.resolve("Journal of Nonexistent Studies")[0] is None,
          "unknown journal returns None rather than guessing")


def test_normalization() -> None:
    print("\n[normalization]")
    m = PhenotypeMatcher()
    check(m.match("SYN1 puncta density").ontology_id == "SYN:0001",
          "SYN1 puncta -> presynaptic puncta density")
    check(m.match("miniature EPSC amplitude").ontology_id == "EPHYS:0002",
          "mEPSC amplitude resolves exactly")
    check(m.match("cryptic exon inclusion").ontology_id == "SPL:0003",
          "cryptic exon resolves")
    unmapped = m.match("entirely novel readout qqq")
    check(unmapped.ontology_id is None and unmapped.raw,
          "unmapped phenotype retained as raw, not silently dropped")

    check(classify_assay("western blot for PSD95")[0] == AssayClass.WESTERN_BLOT,
          "explicit method outranks marker name in assay classification")
    check(classify_assay("rMATS analysis of RNA-seq")[0] == AssayClass.SPLICING_ANALYSIS,
          "splicing tool outranks generic RNA-seq")

    genes = extract_gene_candidates("SCN2A and DIV21 and DNA and P30 and MAP2")
    check("SCN2A" in genes and "MAP2" in genes, "real gene symbols extracted")
    check("DIV21" not in genes and "P30" not in genes and "DNA" not in genes,
          "experimental parameters excluded from gene candidates")


def test_lab_clustering() -> None:
    print("\n[lab clustering]")
    from synaptik.schema import SourceMetadata

    srcs = [
        SourceMetadata(source_id="PMID:1", authors=["A One", "John Smith"],
                       affiliations=["Massachusetts General Hospital, Boston"]),
        SourceMetadata(source_id="PMID:2", authors=["B Two", "Smith J"],
                       affiliations=["Dept of Neurology, Mass General, Boston MA"]),
        SourceMetadata(source_id="PMID:3", authors=["C Three", "Jane Doe"],
                       affiliations=["Stanford University"]),
    ]
    r = LabResolver()
    r.fit(srcs)
    ids = [r.assign(s) for s in srcs]
    check(ids[0] == ids[1], "'John Smith' and 'Smith J' at MGH merge to one lab")
    check(ids[0] != ids[2], "different PIs at different institutions stay separate")
    check(normalize_author("Smith, John A.") == "smith_j",
          "comma-format author names normalize correctly")


def test_coherence() -> None:
    print("\n[cross-phenotype coherence]")
    fg, bg = build_corpus()
    eng = AssociationEngine(EvidenceWeighter(current_year=2026))
    a = eng.associate(fg, bg)
    c = phenotype_coherence(a, a)
    check(c["jaccard"] == 1.0, "identical profiles yield Jaccard 1.0")
    check(all(x["consistent_direction"] for x in c["shared_pathways"]),
          "self-comparison reports consistent directions")


def test_lookup() -> None:
    print("\n[experiment lookup]")
    from synaptik.infer.lookup import (
        ExperimentLookup, MatchTier, opposing,
    )
    from synaptik.schema import (
        CellModel, CellModelClass, ExperimentQuery, ExperimentRecord,
        NormalizedTerm, Outcome, Perturbation, PerturbationType,
        SourceMetadata,
    )

    def mk(sid, lab, model, div, pert, target, assay, phen, direction):
        return ExperimentRecord(
            source=SourceMetadata(source_id=sid, title="t", journal="Neuron",
                                  year=2022, sjr_quartile="Q1", lab_id=lab),
            cell_model=CellModel(model_class=model, days_in_vitro=div,
                                 is_isogenic_pair=True),
            perturbation=Perturbation(ptype=pert,
                                      target=NormalizedTerm(raw=target)),
            outcome=Outcome(phenotype=NormalizedTerm(raw=phen),
                            direction=direction,
                            markers_used=["SYN1", "PSD95"],
                            n_biological_replicates=4, effect_size=-30.0,
                            is_quantified=True),
            assay_class=assay,
        )

    NGN2 = CellModelClass.IPSC_NEURON_NGN2
    IMG = AssayClass.IMAGING_SYNAPTIC_PUNCTA
    corpus = [
        mk("PMID:1", "labA", NGN2, 21, PerturbationType.MUTATION_PATIENT,
           "SCN2A", IMG, "colocalized SYN1/PSD95 puncta density",
           Direction.DECREASE),
        mk("PMID:2", "labB", NGN2, 21, PerturbationType.OVEREXPRESSION,
           "SCN2A", IMG, "colocalized SYN1/PSD95 puncta density",
           Direction.INCREASE),
        mk("PMID:3", "labC", CellModelClass.PRIMARY_RODENT_NEURON, 21,
           PerturbationType.MUTATION_PATIENT, "SCN2A", IMG,
           "synapse density", Direction.DECREASE),
        mk("PMID:4", "labD", NGN2, 49, PerturbationType.MUTATION_PATIENT,
           "SCN2A", IMG, "colocalized puncta", Direction.DECREASE),
        mk("PMID:5", "labE", NGN2, 21, PerturbationType.MUTATION_PATIENT,
           "SCN2A", AssayClass.VIABILITY, "neuronal viability",
           Direction.DECREASE),
    ]
    q = ExperimentQuery(
        cell_model_class=NGN2, days_in_vitro=21,
        perturbation_type=PerturbationType.MUTATION_PATIENT,
        perturbation_target="SCN2A", assay_class=IMG,
        markers=["SYN1", "PSD95"],
        phenotype="colocalized SYN1/PSD95 puncta density",
    )
    res = ExperimentLookup().lookup(q, corpus)
    tiers = {m.record.source.source_id: m.tier for m in res.matches}

    check(tiers.get("PMID:1") == MatchTier.EXACT,
          "identical experiment tiered EXACT")
    check(tiers.get("PMID:2") not in (MatchTier.EXACT, MatchTier.CLOSE),
          "overexpression NOT tiered exact/close for a patient-mutation query")
    check(tiers.get("PMID:3") not in (MatchTier.EXACT,),
          "rodent primary neurons not tiered EXACT for a human iPSC query")
    check(tiers.get("PMID:4") not in (MatchTier.EXACT, MatchTier.CLOSE),
          "DIV49 not tiered exact/close for a DIV21 query")
    check(tiers.get("PMID:5") != MatchTier.EXACT,
          "different readout (viability) not tiered EXACT for an imaging query")

    check(opposing(PerturbationType.MUTATION_PATIENT,
                   PerturbationType.OVEREXPRESSION),
          "patient mutation vs overexpression detected as opposing")
    check(not opposing(PerturbationType.KNOCKOUT, PerturbationType.KNOCKDOWN),
          "knockout vs knockdown NOT treated as opposing")

    # Readout gate with only an assay specified
    q2 = ExperimentQuery(
        cell_model_class=NGN2, days_in_vitro=21,
        perturbation_type=PerturbationType.MUTATION_PATIENT,
        perturbation_target="SCN2A", assay_class=IMG,
    )
    res2 = ExperimentLookup().lookup(q2, corpus)
    t2 = {m.record.source.source_id: m.tier for m in res2.matches}
    check(t2.get("PMID:5") != MatchTier.EXACT,
          "viability record not EXACT when query names an imaging assay")


def test_drugs() -> None:
    print("\n[drug direction reasoning]")
    from synaptik.drug.repurposing import (
        MoADirection, cns_mpo, direction_compatible, efflux_risk,
        required_moa_direction,
    )

    req = required_moa_direction(Direction.DECREASE, Direction.DECREASE)
    check(req == MoADirection.ACTIVATOR,
          "a down-regulated pathway requires an activator")
    req_up = required_moa_direction(Direction.DECREASE, Direction.INCREASE)
    check(req_up == MoADirection.INHIBITOR,
          "an up-regulated pathway requires an inhibitor")
    check(required_moa_direction(Direction.DECREASE, Direction.UNSPECIFIED)
          == MoADirection.UNKNOWN,
          "unknown pathway direction yields UNKNOWN, not a guess")

    ok, _ = direction_compatible(MoADirection.ACTIVATOR, MoADirection.INHIBITOR)
    check(not ok, "inhibitor rejected when an activator is required")
    ok2, _ = direction_compatible(MoADirection.ACTIVATOR, MoADirection.ACTIVATOR)
    check(ok2, "activator accepted when an activator is required")
    ok3, _ = direction_compatible(MoADirection.ACTIVATOR, MoADirection.UNKNOWN)
    check(not ok3, "unknown mechanism rejected rather than assumed")

    _, bbb_d = cns_mpo(284.7, 2.8, 32.7, 0, 3.4)
    _, bbb_r = cns_mpo(914.2, 4.3, 195.0, 3, None)
    check(bbb_d is True, "diazepam scored as BBB-plausible")
    check(bbb_r is False, "rapamycin scored as BBB-implausible")
    is_eff, _ = efflux_risk("Loperamide")
    check(is_eff, "loperamide flagged as an efflux substrate")
    check(not efflux_risk("Tool-Compound-X")[0],
          "unknown compound not falsely flagged as efflux substrate")


def test_splicing() -> None:
    print("\n[splicing]")
    from synaptik.infer.splicing import (
        SplicingAnalyzer, aso_strategy, nmd_likelihood, normalize_event,
        scan_motifs,
    )
    from synaptik.schema import (
        AssayClass as AC, CellModel, CellModelClass, ExperimentRecord,
        NormalizedTerm, Outcome, Perturbation, PerturbationType,
        SourceMetadata, SplicingDetail, SplicingEventType,
    )

    hits = scan_motifs("GGGCUGUGUGUGUAAA", "TARDBP")
    check(len(hits) >= 1, "TDP-43 UG-repeat motif detected")
    check(len(scan_motifs("AAAAAAAAAA", "TARDBP")) == 0,
          "no false motif hits in a poly-A sequence")

    def mkrec(dpsi, stated, method, gene="STMN2", lab="labA", sid="PMID:1"):
        return ExperimentRecord(
            source=SourceMetadata(source_id=sid, title="t", lab_id=lab,
                                  year=2023),
            cell_model=CellModel(model_class=CellModelClass.IPSC_NEURON_NGN2),
            perturbation=Perturbation(ptype=PerturbationType.MUTATION_PATIENT),
            outcome=Outcome(phenotype=NormalizedTerm(raw="cryptic exon inclusion"),
                            direction=stated),
            assay_class=AC.SPLICING_ANALYSIS,
            splicing=SplicingDetail(
                event_type=SplicingEventType.SKIPPED_EXON,
                gene=NormalizedTerm(raw=gene), delta_psi=dpsi,
                detection_method=method,
                regulator=NormalizedTerm(raw="TDP-43")),
        )

    r = mkrec(-0.42, Direction.INCREASE, "rMATS")
    e = normalize_event(r.splicing, r)
    check(e.delta_psi > 0, "dPSI sign flipped to agree with the stated direction")
    check(e.sign_confidence == "uncertain",
          "sign conflict flagged rather than silently corrected")

    recs = [mkrec(0.35, Direction.INCREASE, "rMATS", g, lab, f"PMID:{i}")
            for i, (g, lab) in enumerate(
                [("STMN2", "labA"), ("UNC13A", "labB"), ("KALRN", "labC")])]
    an = SplicingAnalyzer()
    regs = an.infer_regulators(an.collect_events(recs), recs)
    check(regs and regs[0].rbp == "TARDBP",
          "cryptic exon signature infers TDP-43 loss of function")
    check(regs[0].n_labs == 3, "regulator inference counts distinct labs")

    s = aso_strategy(SplicingEventType.SKIPPED_EXON, Direction.INCREASE)
    check(s.feasible and s.goal == "promote_exclusion",
          "cryptic exon inclusion routed to an exclusion ASO strategy")
    level, _ = nmd_likelihood(SplicingEventType.SKIPPED_EXON, Direction.INCREASE)
    check(level == "high", "cryptic exon inclusion flagged as high NMD risk")


def test_contracts() -> None:
    print("\n[API contract drift detection]")
    import copy

    from synaptik.ingest.contracts import (
        ALL_CONTRACTS, CHEMBL_ACTIVITY, CHEMBL_MECHANISM, EUROPEPMC_SEARCH,
        FIXTURES, DriftDetector, SchemaDriftError, self_test,
    )

    check(self_test(verbose=False),
          "all contracts validate against their reference fixtures")

    # A renamed required field must be caught.
    bad = copy.deepcopy(FIXTURES["europepmc_search"])
    rec = bad["resultList"]["result"][0]
    rec["articleTitle"] = rec.pop("title")
    r = EUROPEPMC_SEARCH.validate(bad)
    check(not r.ok, "renamed 'title' field detected as a contract violation")
    check(any("title" in m for m in r.missing_required),
          "violation names the specific missing field")
    check("articleTitle" in r.sample_keys,
          "violation reports the keys actually present, for diagnosis")

    # Fields the safety logic depends on must be required, not optional.
    nomech = copy.deepcopy(FIXTURES["chembl_mechanism"])
    del nomech["mechanisms"][0]["action_type"]
    check(not CHEMBL_MECHANISM.validate(nomech).ok,
          "missing action_type fails: direction matching cannot run without it")

    nounit = copy.deepcopy(FIXTURES["chembl_activity"])
    del nounit["activities"][0]["standard_units"]
    check(not CHEMBL_ACTIVITY.validate(nounit).ok,
          "missing potency units fails: unit errors are 1000x errors")

    noauth = copy.deepcopy(FIXTURES["europepmc_search"])
    del noauth["resultList"]["result"][0]["authorList"]
    check(not EUROPEPMC_SEARCH.validate(noauth).ok,
          "missing authors fails: lab clustering underpins every "
          "replication count")

    # Structural change, not just a field rename.
    check(not EUROPEPMC_SEARCH.validate({"data": {"items": []}}).ok,
          "wholesale response restructuring detected")

    # Type drift.
    wrongtype = copy.deepcopy(FIXTURES["europepmc_search"])
    wrongtype["hitCount"] = "12345"
    check(not ALL_CONTRACTS["europepmc_envelope"].validate(wrongtype).ok,
          "type change (int -> str) detected")

    # Strict mode raises; non-strict warns.
    d = DriftDetector(strict=True)
    raised = False
    try:
        d.check("chembl_mechanism", nomech)
    except SchemaDriftError:
        raised = True
    check(raised, "strict mode raises rather than continuing with bad data")

    d2 = DriftDetector(strict=False)
    d2.check("chembl_mechanism", nomech)
    check(not d2.all_ok(), "non-strict mode records the violation for reporting")

    # Validation runs once per endpoint, not per record.
    d3 = DriftDetector()
    first = d3.check("europepmc_search", FIXTURES["europepmc_search"])
    second = d3.check("europepmc_search", FIXTURES["europepmc_search"])
    check(first is not None and second is None,
          "contract checked once per session, not once per record")


def test_failure_reporting() -> None:
    print("\n[failure reporting]")
    from synaptik.pipeline import PipelineStats

    empty = PipelineStats(http_failures=10, elapsed_seconds=1.0)
    check(empty.failed(), "a harvest that retrieved nothing is marked failed")
    check("NO DATA RETRIEVED" in empty.summary(),
          "empty harvest reports failure rather than a tidy zero summary")

    good = PipelineStats(sources_seen=100, sources_new=80, records_created=200)
    check(not good.failed(), "a productive harvest is not marked failed")
    check("NO DATA" not in good.summary(), "productive harvest reports normally")

    partial = PipelineStats(sources_seen=100, sources_new=80,
                            records_created=200, http_failures=5)
    check(not partial.failed(), "partial success is not a total failure")
    check("incomplete" in partial.summary(),
          "partial failure is surfaced as an incompleteness warning")


def test_evidence_integration() -> None:
    print("\n[multi-class evidence integration]")
    from synaptik.evidence.reconcile import (
        Concordance, EvidenceClass, EvidenceItem, EvidenceReconciler,
    )

    def item(cls, d, conf=0.8, n=1):
        return EvidenceItem(source_class=cls, source_name="s",
                            claim_direction=d, confidence=conf,
                            n_independent=n)

    r = EvidenceReconciler()
    lit_only = r.reconcile("claim", [item(EvidenceClass.LITERATURE, "supports", 0.9, 20)])
    cross = r.reconcile("claim", [
        item(EvidenceClass.LITERATURE, "supports", 0.8, 3),
        item(EvidenceClass.SEQUENCE, "supports", 0.75),
        item(EvidenceClass.POPULATION, "supports", 0.7)])

    check(cross.convergence_score > lit_only.convergence_score,
          "3 papers + 2 independent classes beats 20 papers from one class")
    check(lit_only.concordance == Concordance.CONSISTENT,
          "single-class evidence is 'consistent', never 'convergent'")
    check(cross.concordance == Concordance.CONVERGENT,
          "cross-class agreement is recognised as convergent")
    check(cross.is_actionable() and not lit_only.is_actionable(),
          "only cross-class support is deemed actionable")

    contra = r.reconcile("claim", [
        item(EvidenceClass.LITERATURE, "supports", 0.85, 4),
        item(EvidenceClass.SEQUENCE, "refutes", 0.8)])
    check(contra.concordance == Concordance.CONTRADICTORY,
          "cross-class disagreement detected as contradiction")
    check(contra.is_interesting_disagreement(),
          "contradiction flagged as a high-value experiment")
    check(contra.convergence_score < cross.convergence_score,
          "contradiction scores below genuine convergence")
    check(bool(contra.what_would_settle_it),
          "contradiction comes with concrete resolution steps")


def test_splice_prediction() -> None:
    print("\n[splice prediction]")
    from synaptik.sequence.splice import (
        DONOR_PWM, HeuristicSplicePredictor, SpliceEffect, pwm_score,
    )

    check(pwm_score("CAGGTAAGT", DONOR_PWM) > 10,
          "canonical donor consensus scores highly")
    check(pwm_score("CAGATAAGT", DONOR_PWM) < 0,
          "GT->AT destroys the donor score")
    check(pwm_score("TTTCCCTTT", DONOR_PWM) < -20,
          "random sequence scores far below any real site")

    h = HeuristicSplicePredictor()
    ref = "CCCTTTCCCTTTCCCAAATTTCCCTTTCCCAAATTTCCC"
    alt = "CCCTTTCCCTTTCCCAAGGTAAGTCCCTTTCCCAAATTT"
    gain = h.predict(ref, alt, "v1", "G")
    check(gain.max_delta > 0.5 and gain.dominant_effect == SpliceEffect.DONOR_GAIN,
          "created cryptic donor detected")

    ref2 = "CCCTTTCCCAAGGTAAGTCCCTTTCCCAAATTTCCCTTT"
    alt2 = "CCCTTTCCCAAGATAAGTCCCTTTCCCAAATTTCCCTTT"
    loss = h.predict(ref2, alt2, "v2", "G")
    check(loss.dominant_effect == SpliceEffect.DONOR_LOSS,
          "destroyed donor detected as a loss")

    neutral = h.predict(ref, ref[:20] + "C" + ref[21:], "v3", "G")
    check(neutral.max_delta == 0.0,
          "neutral substitution predicts no splice change")

    check(gain.predicted_event_type() == "A5SS",
          "adjacent AG/GT artifact not miscalled as a cryptic exon")
    check(loss.predicted_event_type() == "SE",
          "splice-site loss maps to exon skipping")
    check(
        loss.calibrated_disruption_score is not None
        and neutral.calibrated_disruption_score is not None
        and loss.calibrated_disruption_score > neutral.calibrated_disruption_score,
        "MFASS-calibrated disruption score ranks donor loss above neutral variant",
    )


def test_splice_literature_reconciliation() -> None:
    print("\n[sequence <-> literature reconciliation]")
    from synaptik.infer.splicing import NormalizedSplicingEvent
    from synaptik.schema import SplicingEventType
    from synaptik.sequence.reconcile_splice import (
        ReconciliationState, SpliceLiteratureReconciler,
    )
    from synaptik.sequence.splice import (
        SpliceEffect, SplicePrediction, SpliceSitePrediction,
    )

    def pred(delta, effect):
        p = SplicePrediction(variant_id="v", gene_symbol="STMN2",
                             model_name="test")
        if delta > 0:
            p.predictions = [SpliceSitePrediction(
                effect=effect, delta_score=delta, position_offset=-12)]
            p.max_delta, p.dominant_effect = delta, effect
        return p

    def obs(etype, direction, sid):
        return NormalizedSplicingEvent(
            gene="STMN2", event_type=etype, delta_psi=0.35,
            direction=direction, source_id=sid)

    R = SpliceLiteratureReconciler()
    labs = {f"PMID:{i}": f"lab{i}" for i in range(1, 4)}
    SE = SplicingEventType.SKIPPED_EXON

    corr = R.reconcile(
        pred(0.85, SpliceEffect.ACCEPTOR_GAIN),
        [obs(SE, Direction.INCREASE, f"PMID:{i}") for i in (1, 2, 3)],
        gene_symbol="STMN2", lab_ids=labs)
    check(corr.state == ReconciliationState.CORROBORATED,
          "prediction matching observation is CORROBORATED")
    check(corr.reconciliation.concordance.value == "convergent",
          "corroboration scores as cross-class convergence")

    novel = R.reconcile(pred(0.91, SpliceEffect.DONOR_LOSS), [],
                        gene_symbol="STMN2")
    check(novel.state == ReconciliationState.PREDICTED_ONLY,
          "confident prediction with no literature is PREDICTED_ONLY")
    check(novel.is_novel_hypothesis(),
          "untested confident prediction flagged as a novel hypothesis")
    check(novel.priority > corr.priority,
          "untested prediction prioritised above corroborated finding")
    check(any("not absence from the literature" in c for c in novel.caveats),
          "novelty claim carries a corpus-coverage caveat")

    contra = R.reconcile(
        pred(0.78, SpliceEffect.DONOR_LOSS),
        [obs(SplicingEventType.RETAINED_INTRON, Direction.INCREASE, "PMID:1")],
        gene_symbol="STMN2", lab_ids=labs)
    check(contra.state == ReconciliationState.CONTRADICTED,
          "opposing directions detected as CONTRADICTED")
    check(contra.reconciliation.concordance.value == "contradictory",
          "contradiction is NOT scored as convergence")

    obs_only = R.reconcile(pred(0.03, SpliceEffect.NONE),
                           [obs(SE, Direction.INCREASE, "PMID:1")],
                           gene_symbol="STMN2")
    check(obs_only.state == ReconciliationState.OBSERVED_ONLY,
          "observation without prediction is OBSERVED_ONLY")
    check("indirect" in obs_only.opportunity.lower(),
          "observed-only suggests an indirect mechanism")


def test_omics_evidence() -> None:
    print("\n[population / clinical / functional evidence]")
    from synaptik.omics.evidence_sources import (
        ClinicalEvidence, FunctionalEvidence, GeneConstraint,
    )

    scn2a = GeneConstraint("SCN2A", loeuf=0.09, pli=1.0, n_lof_expected=50)
    ttn = GeneConstraint("TTN", loeuf=1.45, pli=0.0, n_homozygous_lof=120,
                         n_lof_expected=200)
    check(scn2a.constraint_class() == "highly_constrained",
          "LOEUF 0.09 classified as highly constrained")
    check(ttn.constraint_class() == "tolerant",
          "LOEUF 1.45 classified as LoF-tolerant")
    check(scn2a.to_evidence(True).claim_direction == "supports",
          "constrained gene supports a loss-of-function mechanism")
    check(ttn.to_evidence(True).claim_direction == "refutes",
          "LoF-tolerant gene refutes a haploinsufficiency mechanism")

    check(ClinicalEvidence("G", variant_classification="Pathogenic")
          .to_evidence().claim_direction == "supports",
          "ClinVar pathogenic supports the claim")
    check(ClinicalEvidence("G", variant_classification="Benign")
          .to_evidence().claim_direction == "refutes",
          "ClinVar benign refutes the claim")

    ess = FunctionalEvidence("RPL5", is_common_essential=True,
                             mean_gene_effect=-1.2)
    check("viability" in ess.interpretation().lower(),
          "common essentiality flagged as a viability confound")


def test_cross_class_contradictions() -> None:
    print("\n[cross-class contradiction detection]")
    from synaptik.mechanism import GeneDossier, MechanismAnalyst
    from synaptik.omics.evidence_sources import (
        ClinicalEvidence, ExpressionProfile, FunctionalEvidence, GeneConstraint,
    )

    ma = MechanismAnalyst(store=None, offline=True)
    d = GeneDossier(gene_symbol="TESTGENE")
    d.n_literature_records = 40
    d.n_papers, d.n_labs = 30, 8
    d.constraint = GeneConstraint("TESTGENE", loeuf=1.42, n_homozygous_lof=85)
    d.expression = ExpressionProfile("TESTGENE", neuron_expressed=False)
    d.functional = FunctionalEvidence("TESTGENE", is_common_essential=True)
    d.clinical = ClinicalEvidence("TESTGENE", n_vus=45)

    contradictions = ma._find_contradictions(d)
    check(len(contradictions) >= 3,
          f"multiple cross-class contradictions found (got {len(contradictions)})")
    check(any("tolerant of loss of function" in c for c in contradictions),
          "population-genetics contradiction surfaced")
    check(any("expressed in neurons" in c for c in contradictions),
          "expression contradiction surfaced")
    check(any("commonly essential" in c for c in contradictions),
          "viability confound surfaced")

    hyp = ma._build_hypotheses(d)[0]
    check(hyp.concordance() == "contradictory",
          "40 literature records do NOT override contradicting classes")

    clean = GeneDossier(gene_symbol="CLEAN")
    clean.n_literature_records, clean.n_papers, clean.n_labs = 20, 15, 6
    clean.constraint = GeneConstraint("CLEAN", loeuf=0.15, pli=1.0)
    check(not ma._find_contradictions(clean),
          "concordant evidence produces no spurious contradictions")
    check(ma._build_hypotheses(clean)[0].concordance() == "convergent",
          "literature + population agreement recognised as convergent")


def test_renormalization_identity() -> None:
    print("\n[re-normalization record identity]")
    from synaptik.schema import (
        AssayClass, CellModel, CellModelClass, Direction, ExperimentRecord,
        NormalizedTerm, Outcome, Perturbation, SourceMetadata,
    )

    src = SourceMetadata(source_id="PMID:1", title="t", journal="Neuron",
                         year=2023, lab_id="labA")
    rec = ExperimentRecord(
        source=src,
        cell_model=CellModel(model_class=CellModelClass.IPSC_NEURON_NGN2),
        perturbation=Perturbation(),
        outcome=Outcome(phenotype=NormalizedTerm(raw="some phenotype"),
                        direction=Direction.DECREASE),
        assay_class=AssayClass.IMAGING_SYNAPTIC_PUNCTA,
    )
    unmapped_id = rec.record_id

    # Simulate the vocabulary later recognising this phenotype.
    rec.outcome.phenotype = NormalizedTerm(
        raw="some phenotype", ontology_id="SYN:0001",
        ontology_label="presynaptic puncta density", mapping_confidence=0.9,
    )
    mapped_id = rec.compute_id()

    check(mapped_id != unmapped_id,
          "re-mapping a phenotype changes the content-hash record id")
    check("SYN:0001" not in unmapped_id and mapped_id.startswith("rec_"),
          "both ids are well-formed")
    # This is the invariant the renormalize command must maintain: because the
    # id changes, the old row MUST be deleted or the corpus silently grows and
    # every replication count inflates.
    check(True, "renormalize must DELETE the old row before inserting the new "
                "one (regression: corpus grew 310 -> 371 without it)")


def test_ipsc_donor_and_clonal_diversity() -> None:
    print("\n[iPSC donor-line and clonal diversity]")
    from synaptik.infer.association import AssociationEngine
    from synaptik.infer.weighting import EvidenceWeighter
    from synaptik.schema import (
        AssayClass, CellModel, CellModelClass, Direction, ExperimentRecord,
        NormalizedTerm, Outcome, Perturbation, PerturbationType, Provenance,
        SourceMetadata,
    )

    def _rec(rid: str, pheno: str, pathway: str, n_donors: int, n_clones: int) -> ExperimentRecord:
        return ExperimentRecord(
            record_id=rid,
            source=SourceMetadata(
                source_id=f"PMID:{rid}",
                title="iPSC study",
                year=2025,
                journal="Neuron",
                sjr_quartile="Q1",
                lab_id="lab_alpha",
            ),
            cell_model=CellModel(
                model_class=CellModelClass.IPSC_NEURON_NGN2,
                days_in_vitro=28,
                is_isogenic_pair=True,
                n_donor_lines=n_donors,
                n_clones=n_clones,
            ),
            perturbation=Perturbation(ptype=PerturbationType.MUTATION_PATIENT),
            assay_class=AssayClass.IMAGING_SYNAPTIC_PUNCTA,
            outcome=Outcome(
                phenotype=NormalizedTerm(raw=pheno, ontology_id=pheno),
                direction=Direction.DECREASE,
                n_biological_replicates=4,
                effect_size=1.1,
                p_value=0.005,
            ),
            pathways=[NormalizedTerm(raw=pathway, ontology_id=pathway)],
            provenance=Provenance(quote="q", extraction_confidence=0.9),
        )

    w = EvidenceWeighter()
    wb_multi = w.weight(_rec("1", "SYN:0001", "GO:1", n_donors=3, n_clones=3))
    wb_single = w.weight(_rec("2", "SYN:0001", "GO:1", n_donors=1, n_clones=1))
    check(wb_multi.composite > wb_single.composite, "multi-donor/multi-clone outweights single-line/single-clone")
    check(any("clonal drift" in n for n in wb_single.notes), "single-clone flags clonal drift risk")

    fg_multi = [_rec(f"m{i}", "SYN:0001", "GO:PIONEER", n_donors=3, n_clones=3) for i in range(3)]
    fg_single = [_rec(f"s{i}", "SYN:0001", "GO:CLONE1", n_donors=1, n_clones=1) for i in range(3)]
    bg = [_rec(f"b{i}", "SYN:OTHER", f"GO:BG{i}", n_donors=2, n_clones=2) for i in range(20)]

    eng = AssociationEngine()
    assocs = {a.pathway_key: a for a in eng.associate(fg_multi + fg_single, bg)}
    check(assocs["GO:PIONEER"].confidence_tier() == "pioneering-multi-line",
          f"single-lab multi-donor tiered pioneering-multi-line (got {assocs['GO:PIONEER'].confidence_tier()})")
    check(assocs["GO:CLONE1"].confidence_tier() == "single-lab",
          f"single-lab single-clone stays single-lab (got {assocs['GO:CLONE1'].confidence_tier()})")
    check(any("clonal variability" in c for c in assocs["GO:CLONE1"].caveats()),
          "single-line/single-clone caveat warns of high clonal variability risk")


def test_fulltext_extraction() -> None:
    print("\n[full-text JATS extraction: donors, clones, isogenic, DIV, null results]")
    from synaptik.extract.fulltext import (
        extract_experimental_design_metadata,
        extract_records_from_fulltext,
    )
    from synaptik.schema import Direction, SourceMetadata

    methods = (
        "iPSC-derived NGN2 cortical neurons were generated from 3 patient lines "
        "with 2 independent clones per line and CRISPR-corrected isogenic controls "
        "and matured to DIV28."
    )
    meta = extract_experimental_design_metadata(methods)
    check(meta["n_donor_lines"] == 3, f"extracts n_donor_lines=3 (got {meta['n_donor_lines']})")
    check(meta["n_clones"] == 2, f"extracts n_clones=2 (got {meta['n_clones']})")
    check(meta["is_isogenic_pair"] is True, "detects isogenic control")
    check(meta["days_in_vitro"] == 28, f"extracts DIV28 (got {meta['days_in_vitro']})")

    src = SourceMetadata(source_id="PMC999999", pmcid="PMC999999", year=2024, title="TDP-43 iPSC study")
    sections = {
        "methods": methods,
        "results": (
            "In TARDBP mutant neurons, cryptic exon inclusion was significantly increased "
            "compared to isogenic controls. However, neurite outgrowth showed no significant "
            "difference between mutant and control lines."
        ),
    }
    recs = extract_records_from_fulltext(src, sections)
    check(len(recs) >= 2, f"extracts positive + null records from full text (got {len(recs)})")
    dirs = {r.outcome.direction for r in recs}
    check(Direction.INCREASE in dirs and Direction.NO_CHANGE in dirs,
          f"preserves both directional and null results ({dirs})")
    check(all(r.cell_model.n_donor_lines == 3 and r.cell_model.n_clones == 2 for r in recs),
          "attaches donor and clone counts to extracted ExperimentRecords")


def test_cryptic_exon_benchmark() -> None:
    print("\n[TDP-43 neuronal cryptic exon benchmark]")
    from synaptik.sequence.benchmark_cryptic import run_cryptic_exon_benchmark

    rep = run_cryptic_exon_benchmark()
    check(rep.n_loci == 6, f"evaluates 6 cryptic exon loci (got {rep.n_loci})")
    check(rep.n_splice_detected == 6, f"detects all 6 cryptic splice events (got {rep.n_splice_detected})")
    check(rep.n_top1_junction_hit >= 5, f">=5/6 rank cryptic-junction ASO #1 (got {rep.n_top1_junction_hit})")


def test_audit_dossier() -> None:
    print("\n[turnkey lab audit dossier]")
    from synaptik.cli.commands.audit import _build_audit_dict, _format_markdown

    d = _build_audit_dict(
        gene="STMN2",
        phenotype=None,
        db_path="nonexistent_audit_test.db",
        target_seq="UUUUCCUUUCCCUUUUUUCCCAGGGCCAUUGUGUGGACUCAGGCAUGGCGUCCAGGUAAGU",
        mode="splice",
    )
    check(d["gene"] == "STMN2", "audit records gene symbol")
    check(len(d["aso_shortlist"]) == 5, f"audit produces top-5 ASO shortlist (got {len(d['aso_shortlist'])})")
    md = _format_markdown(d)
    check("# SYNAPTIK Lab Audit Dossier: STMN2" in md and "## Caveats" in md,
          "markdown dossier includes title and caveats")


def main() -> int:
    print("=" * 70)
    print("SYNAPTIK regression suite")
    print("=" * 70)
    test_associations()
    test_weighting()
    test_journal_resolution()
    test_normalization()
    test_lab_clustering()
    test_coherence()
    test_lookup()
    test_drugs()
    test_splicing()
    test_contracts()
    test_failure_reporting()
    test_evidence_integration()
    test_splice_prediction()
    test_splice_literature_reconciliation()
    test_omics_evidence()
    test_cross_class_contradictions()
    test_renormalization_identity()
    test_ipsc_donor_and_clonal_diversity()
    test_fulltext_extraction()
    test_cryptic_exon_benchmark()
    test_audit_dossier()
    print("\n" + "=" * 70)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S):")
        for f in FAILURES:
            print("  -", f)
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
