"""
Regression tests for the lab-facing workflows: bring-your-own-data,
experiment design briefs, the retrospective benchmark, and phenotype pooling.

As elsewhere in this suite, the assertions are scientific properties. Each
one names a way the tool could mislead a researcher, and fails if it does.

Run: python3 -m tests.test_lab_workflows
"""

from __future__ import annotations

import sys

from synaptik.infer.benchmark import (
    RetrospectiveBenchmark, auroc, bootstrap_auroc_difference,
)
from synaptik.infer.design import (
    DesignAdvisor, detectable_effect_d, sample_size_paired, sample_size_two_group,
)
from synaptik.infer.lookup import ExperimentLookup
from synaptik.infer.splicing import NormalizedSplicingEvent
from synaptik.sequence.benchmark_spliceai import run_spliceai_benchmark
from synaptik.labdata import (
    LabCall, LabDataError, LabDataReconciler, UserSplicingEvent,
    frame_consequence, parse_phenotype_table, parse_splicing_table,
)
from synaptik.normalize.matcher import PhenotypeMatcher
from synaptik.normalize.phenotype_ontology import commensurable
from synaptik.schema import (
    AssayClass, CellModel, CellModelClass, Direction, ExperimentQuery,
    ExperimentRecord, NormalizedTerm, Outcome, Perturbation, PerturbationType,
    Provenance, SourceMetadata, SplicingEventType,
)

FAILURES: list[str] = []


def check(cond: bool, msg: str) -> None:
    print(f"  {'PASS' if cond else 'FAIL'}  {msg}")
    if not cond:
        FAILURES.append(msg)


RMATS_HEADER = (
    "ID\tGeneID\tgeneSymbol\tchr\tstrand\texonStart_0base\texonEnd\t"
    "upstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\tID\tIJC_SAMPLE_1\t"
    "SJC_SAMPLE_1\tIJC_SAMPLE_2\tSJC_SAMPLE_2\tIncFormLen\tSkipFormLen\t"
    "PValue\tFDR\tIncLevel1\tIncLevel2\tIncLevelDifference"
)


def rmats(rows: list[tuple]) -> str:
    """rows: (gene, dpsi, fdr, exon_length)"""
    lines = [RMATS_HEADER]
    for i, (gene, dpsi, fdr, length) in enumerate(rows):
        lines.append("\t".join(map(str, [
            i, f'"ENSG{i}"', f'"{gene}"', "chr1", "+", 1000, 1000 + length,
            500, 600, 2000, 2100, i, "10,10", "5,5", "5,5", "10,10", 200, 100,
            fdr / 10, fdr, "0.6,0.6", "0.4,0.4", dpsi])))
    return "\n".join(lines) + "\n"


def lit(gene, direction, sid, etype=SplicingEventType.SKIPPED_EXON,
        regulator="TDP-43", conf="high"):
    return NormalizedSplicingEvent(
        gene=gene, event_type=etype, delta_psi=0.3, direction=direction,
        source_id=sid, regulator=regulator, sign_confidence=conf)


# --------------------------------------------------------------------------


def test_splicing_table_parsing() -> None:
    print("\n[lab data: parsing]")
    text = rmats([("STMN2", 0.40, 0.001, 120), ("KALRN", -0.30, 0.01, 97)])

    raised = False
    try:
        parse_splicing_table("SE.MATS.JC.txt", perturbed="", text=text)
    except LabDataError as e:
        raised = "no safe default" in str(e)
    check(raised, "rMATS table without a stated perturbed group is REFUSED, "
                  "not guessed")

    t1 = parse_splicing_table("SE.MATS.JC.txt", perturbed="sample1", text=text)
    t2 = parse_splicing_table("SE.MATS.JC.txt", perturbed="sample2", text=text)
    d1 = {e.gene: e.delta_psi for e in t1.events}
    d2 = {e.gene: e.delta_psi for e in t2.events}
    check(d1["STMN2"] > 0 and d2["STMN2"] < 0,
          "declaring sample2 as perturbed flips every dPSI sign")
    check(t1.events[0].event_type == SplicingEventType.SKIPPED_EXON,
          "event type inferred from the rMATS file name")
    check(t1.events[0].gene == "STMN2",
          "quoted rMATS gene symbols are unquoted")
    check(t1.events[0].exon_length == 120,
          "exon length derived from the table's coordinates")

    bad = rmats([("A", 40.0, 0.01, 120), ("B", 0.2, 0.01, 120),
                 ("C", "NA", 0.01, 120)])
    t = parse_splicing_table("SE.MATS.JC.txt", perturbed="sample1", text=bad)
    check(len(t.events) == 1 and t.n_rows_skipped == 2,
          "percent-scale and missing dPSI rows are skipped, not rescaled or "
          "zero-filled")
    check(any("outside [-1, 1]" in w for w in t.warnings),
          "skipped out-of-range rows are reported to the user")

    generic = "gene,delta_psi\nSTMN2,0.4\nKALRN,-0.3\n"
    g = parse_splicing_table("mine.csv", perturbed="positive", text=generic)
    check(len(g.events) == 2 and g.events[0].fdr is None,
          "generic table without FDR parses with FDR left as None")
    check(any("No FDR column" in w for w in g.warnings),
          "missing FDR is flagged: significance cannot be assessed")
    check(any("Event type unknown" in w for w in g.warnings),
          "unknown event type is flagged, not assumed to be SE")

    raised = False
    try:
        parse_splicing_table("x.csv", perturbed="positive",
                             text="foo,bar\n1,2\n")
    except LabDataError:
        raised = True
    check(raised, "a table without gene/dPSI columns is rejected with an "
                  "explanation")

    # MAJIQ VOILA parsing test
    majiq_text = (
        "Gene Name\tGene ID\tLSV ID\tE(dPSI) per LSV junction\tP(|dPSI| > 0.20) per LSV junction\tLSV Type\tJunctions coords\n"
        "STMN2\tENSG1\tlsv_1\t0.420;0.050\t0.00001;0.85000\ts\tchr8:79611214-79611334;chr8:79611214-79615000\n"
        "UNC13A\tENSG2\tlsv_2\t0.310\t0.00002\tt\tchr19:17642414-17642511\n"
    )
    m1 = parse_splicing_table("voila.tsv", perturbed="sample1", text=majiq_text)
    check(len(m1.events) == 2, "MAJIQ VOILA table parses successfully")
    check(m1.events[0].gene == "STMN2" and abs(m1.events[0].delta_psi - 0.420) < 1e-4,
          "MAJIQ multi-junction event selects max |E(dPSI)|")
    check(m1.events[0].exon_length == 120, "MAJIQ junction coordinates parsed for exon length")

    # LeafCutter parsing test
    leaf_text = (
        "cluster\tgenes\tlogef\tp\tp.adjust\n"
        "chr8:79611214:79611334:clu_1_+\tSTMN2\t0.420\t1e-10\t1e-09\n"
        "chr19:17642414:17642511:clu_2_+\tUNC13A\t-0.310\t2e-07\t2e-06\n"
    )
    l1 = parse_splicing_table("leafcutter_ds.tsv", perturbed="positive", text=leaf_text)
    check(len(l1.events) == 2, "LeafCutter differential splicing table parses successfully")
    check(l1.events[0].gene == "STMN2" and l1.events[0].exon_length == 120,
          "LeafCutter cluster coordinates parsed for coordinates and length")
    check(l1.events[0].fdr == 1e-09, "LeafCutter p.adjust parsed as FDR")


def test_splicing_reconciliation() -> None:
    print("\n[lab data: reconciliation with the literature]")
    rows = [("STMN2", 0.42, 1e-6, 120), ("UNC13A", 0.31, 1e-4, 97),
            ("KALRN", -0.28, 1e-3, 120), ("NEWGENE", 0.35, 1e-4, 143),
            ("AGRN", 0.02, 0.7, 120), ("WEAK", 0.06, 0.001, 120),
            ("NOTSIG", 0.4, 0.4, 120)]
    rows += [(f"BKG{i}", 0.01, 0.9, 120) for i in range(60)]
    table = parse_splicing_table(
        "SE.MATS.JC.txt", perturbed="sample1", text=rmats(rows))
    corpus = [
        lit("STMN2", Direction.INCREASE, "PMID:1"),
        lit("STMN2", Direction.INCREASE, "PMID:2"),
        lit("UNC13A", Direction.INCREASE, "PMID:3"),
        lit("KALRN", Direction.INCREASE, "PMID:4"),
        lit("AGRN", Direction.INCREASE, "PMID:5"),
        lit("ABSENT", Direction.INCREASE, "PMID:6"),
    ]
    labs = {f"PMID:{i}": f"lab{i}" for i in range(1, 7)}
    rep = LabDataReconciler().reconcile_splicing(table, corpus, labs)
    calls = {r.event.gene: r.call for r in rep.reconciled}

    check(calls.get("STMN2") == LabCall.CONSISTENT,
          "same gene, class and direction as the literature -> consistent")
    check(calls.get("KALRN") == LabCall.DISCORDANT,
          "opposite direction to the literature -> discordant")
    check(calls.get("NEWGENE") == LabCall.NOVEL,
          "gene with no published event -> novel")
    check("WEAK" not in calls and "NOTSIG" not in calls,
          "events failing the effect-size or FDR threshold are not called")
    check(rep.reconciled[0].call in (LabCall.DISCORDANT, LabCall.NOVEL),
          "contradictions and novel events are ranked above confirmations")
    novel = next(r for r in rep.reconciled if r.event.gene == "NEWGENE")
    check("Check PubMed" in novel.note,
          "a 'novel' call carries the corpus-coverage caveat")

    nr = {n.gene for n in rep.not_reproduced}
    check("AGRN" in nr,
          "published event in a gene tested-but-unchanged is reported as not "
          "reproduced")
    check("ABSENT" not in nr,
          "a gene absent from the user's table is NOT reported as a negative "
          "result")

    tdp = next((g for g in rep.regulators if g.rbp == "TARDBP"), None)
    check(tdp is not None and tdp.n_targets_tested == 4
          and tdp.n_targets_changed == 3,
          "regulator enrichment counts only targets tested in the user's table")
    check(tdp is not None and tdp.q_value < 0.05,
          "published TDP-43 targets are enriched among the changed genes")
    check(tdp is not None and tdp.direction_agreement is not None
          and abs(tdp.direction_agreement - 2 / 3) < 1e-6,
          "direction agreement with the published targets is reported")

    empty = LabDataReconciler().reconcile_splicing(table, [], {})
    check(all(r.call == LabCall.NOVEL for r in empty.reconciled)
          and any("reflects the corpus" in w for w in empty.warnings),
          "with an empty corpus, 'novel' is explained as a corpus limitation")

    typed = lit("STMN2", Direction.INCREASE, "PMID:9",
                etype=SplicingEventType.RETAINED_INTRON)
    other = LabDataReconciler().reconcile_splicing(table, [typed], {})
    c2 = {r.event.gene: r.call for r in other.reconciled}
    check(c2.get("STMN2") == LabCall.SAME_GENE,
          "a different event class in the same gene is not called consistent")


def test_frame_consequence() -> None:
    print("\n[lab data: reading-frame consequence]")
    SE = SplicingEventType.SKIPPED_EXON
    shift = UserSplicingEvent("G", SE, 0.3, exon_length=143)
    inframe = UserSplicingEvent("G", SE, 0.3, exon_length=120)
    nocoord = UserSplicingEvent("G", SE, 0.3)
    check(frame_consequence(shift)[0] == "high",
          "exon length not divisible by three flagged as frame-shifting")
    check(frame_consequence(inframe)[0] == "low",
          "in-frame cassette exon is NOT flagged as an NMD risk")
    check(frame_consequence(nocoord)[0] == "unknown",
          "without coordinates the frame consequence is unknown, not assumed")


def test_phenotype_table() -> None:
    print("\n[lab data: phenotype table]")
    text = ("phenotype,direction,gene\nSYN1 puncta density,decrease,TARDBP\n"
            "mEPSC frequency,sideways,TARDBP\nviability,unchanged,\n")
    rows, warnings = parse_phenotype_table("p.csv", text=text)
    check(len(rows) == 2, "rows with an unreadable direction are skipped")
    check(any("not understood" in w for w in warnings),
          "an unreadable direction is reported, never read as 'no change'")
    check(rows[1].direction == Direction.NO_CHANGE and rows[1].gene is None,
          "'unchanged' parses as an explicit null result")


# --------------------------------------------------------------------------


def mk_record(sid, lab, direction, es=None, n=3, div=21, isogenic=True,
              assay=AssayClass.IMAGING_SYNAPTIC_PUNCTA, year=2022,
              pathways=("p",), phenotype="colocalized puncta density"):
    return ExperimentRecord(
        source=SourceMetadata(source_id=sid, title="t", journal="Neuron",
                              year=year, sjr_quartile="Q1", lab_id=lab,
                              authors=["A B", "C D"]),
        cell_model=CellModel(model_class=CellModelClass.IPSC_NEURON_NGN2,
                             days_in_vitro=div, is_isogenic_pair=isogenic,
                             coculture="mouse astrocytes"),
        perturbation=Perturbation(ptype=PerturbationType.MUTATION_PATIENT,
                                  target=NormalizedTerm(raw="SCN2A")),
        outcome=Outcome(
            phenotype=NormalizedTerm(raw=phenotype, ontology_id="SYN:0003"),
            direction=direction, effect_size=es,
            effect_size_unit="percent_change" if es is not None else None,
            n_biological_replicates=n, markers_used=["SYN1", "PSD95"],
            is_quantified=es is not None),
        assay_class=assay,
        pathways=[NormalizedTerm(raw=p, ontology_id=p) for p in pathways],
        provenance=Provenance(quote="q" * 20, extraction_confidence=0.9),
    )


QUERY = ExperimentQuery(
    cell_model_class=CellModelClass.IPSC_NEURON_NGN2, days_in_vitro=21,
    perturbation_type=PerturbationType.MUTATION_PATIENT,
    perturbation_target="SCN2A",
    assay_class=AssayClass.IMAGING_SYNAPTIC_PUNCTA,
    markers=["SYN1", "PSD95"], phenotype="colocalized puncta density")


def brief_for(corpus, cv=None):
    return DesignAdvisor().brief(ExperimentLookup().lookup(QUERY, corpus), cv=cv)


def test_power() -> None:
    print("\n[design: power calculation]")
    got = [sample_size_two_group(d) for d in (0.5, 0.8, 1.0, 2.0)]
    check(got == [64, 26, 17, 6],
          f"sample sizes match the standard tables for d=0.5/0.8/1/2 "
          f"(got {got})")
    check(sample_size_two_group(0.05) is None,
          "an effect needing >200 per group returns None, not a huge number")
    check(sample_size_two_group(0.0) is None,
          "a zero effect has no sample size")
    ds = [detectable_effect_d(n) for n in (3, 4, 6, 8, 16)]
    check(all(a > b for a, b in zip(ds, ds[1:])),
          "detectable effect shrinks monotonically as n grows")
    check(abs(detectable_effect_d(64) - 0.5) < 0.02,
          "detectable effect at n=64 inverts the sample-size calculation")

    # Paired isogenic power analysis
    check(sample_size_paired(0.8, rho=0.7) == 10,
          "paired sample size d=0.8, rho=0.7 yields n=10 pairs")
    check(sample_size_paired(0.8, rho=0.8) == 8,
          "higher inter-batch correlation rho=0.8 reduces required pairs to 8")
    check(sample_size_paired(0.8, rho=0.0) in (26, 27),
          "rho=0 in paired test matches independent two-group requirement (~26-27)")
    check(sample_size_paired(0.0) is None and sample_size_paired(0.8, rho=1.0) is None,
          "invalid effect size or correlation in paired power returns None")


def test_design_brief() -> None:
    print("\n[design: brief]")
    single_lab = [mk_record(f"PMID:{i}", "labA", Direction.DECREASE, es=-40.0)
                  for i in range(5)]
    b = brief_for(single_lab)
    check(b.n_exact == 5 and b.n_labs == 1, "prior experiments counted by lab")
    check(any("one research group" in w for w in b.weaknesses),
          "single-lab prior evidence is named as a weakness")
    check("replication" in b.value_statement.lower(),
          "a single-lab finding is framed as worth replicating")
    check(any("No null result" in w for w in b.weaknesses),
          "all-positive prior evidence is flagged as likely inflated")
    check(b.median_effect_percent == 40.0, "median published effect reported")

    published = [p for p in b.power if p.label == "published median"]
    halved = [p for p in b.power if p.label.startswith("halved")]
    check(len(published) == 3 and len(halved) == 3,
          "sample size is shown across a range of CVs when none is supplied")
    check(all(h.n_per_group is None or h.n_per_group > p.n_per_group
              for h, p in zip(halved, published)),
          "planning on the halved effect always needs more replicates")
    check(all(p.n_paired is not None and (p.n_per_group is None or p.n_paired < p.n_per_group)
              for p in published),
          "paired isogenic design requires fewer differentiations across all CVs")
    check(b.recommended_n_paired is not None and b.recommended_n_paired < b.recommended_n,
          "paired isogenic design recommendation requires fewer differentiations than two-group")
    check("assumption" in b.recommended_n_basis,
          "the recommended n states that its CV is an assumption")

    own = brief_for(single_lab, cv=0.10)
    check(len(own.power) == 2 and "assumption" not in own.recommended_n_basis,
          "a user-supplied CV replaces the assumed range")
    check(own.recommended_n < b.recommended_n,
          "a lower measured CV lowers the required sample size")

    small = [mk_record(f"PMID:{i}", f"lab{i}", Direction.DECREASE, es=-12.0)
             for i in range(4)]
    s = brief_for(small)
    check(s.impractical and any("PAIRS" in r for r in s.recommendations),
          "an impractical sample size leads to a paired-design "
          "recommendation, not a bare number")
    check(any("state that limit in advance" in r for r in s.recommendations),
          "the brief states what a feasible design can detect")

    split = [mk_record("PMID:1", "lab1", Direction.DECREASE, es=-30.0),
             mk_record("PMID:2", "lab2", Direction.INCREASE, es=25.0),
             mk_record("PMID:3", "lab3", Direction.DECREASE, es=-20.0),
             mk_record("PMID:4", "lab4", Direction.INCREASE, es=30.0)]
    d = brief_for(split)
    check(d.consensus_direction == "disputed"
          and "Tie-breaker" in d.value_statement,
          "a split literature is framed as a tie-breaking opportunity")

    none = brief_for([])
    check(none.n_exact == 0 and not none.power
          and "Unexplored" in none.value_statement,
          "with no prior work the brief invents no protocol or effect size")
    check(any("pilot" in r for r in none.recommendations),
          "with no prior work a pilot is recommended")
    check(any("dendrite length" in p for p in b.pitfalls),
          "readout-specific pitfalls are included for synaptic imaging")
    check(any("independent differentiation" in p for p in b.pitfalls),
          "the unit-of-replication warning is always present")

    no_es = [mk_record(f"PMID:{i}", f"lab{i}", Direction.DECREASE)
             for i in range(3)]
    check(not brief_for(no_es).power,
          "no power estimate is produced when no effect size was reported")


# --------------------------------------------------------------------------


def test_benchmark() -> None:
    print("\n[benchmark: temporal hold-out]")
    check(auroc([3, 4, 5], [0, 1, 2]) == 1.0, "perfect ranking has AUROC 1")
    check(auroc([0, 1, 2], [3, 4, 5]) == 0.0, "inverted ranking has AUROC 0")
    check(auroc([1, 1], [1, 1]) == 0.5, "ties score 0.5")
    check(auroc([], [1, 2]) is None,
          "AUROC with an empty class is undefined, not 0.5")

    labels = [True] * 30 + [False] * 30
    good = [1.0] * 30 + [0.0] * 30
    noise = [(i * 7919) % 13 for i in range(60)]
    ci = bootstrap_auroc_difference(labels, good, noise)
    check(ci is not None and ci[0] > 0,
          "bootstrap interval excludes zero for a clearly better ranking")
    ci0 = bootstrap_auroc_difference(labels, noise, noise)
    check(ci0 is not None and ci0[0] <= 0 <= ci0[1],
          "bootstrap interval includes zero for identical rankings")
    check(bootstrap_auroc_difference(labels, good, noise)
          == bootstrap_auroc_difference(labels, good, noise),
          "the bootstrap is deterministic for a given seed")

    # Before the cutoff: 'real' is reported by 3 labs, 'solo' by one lab many
    # times. After it: a NEW lab confirms 'real'; only the SAME lab repeats
    # 'solo'; a new lab reports the opposite of 'flip'.
    recs = []
    for i, lab in enumerate(["l1", "l2", "l3"]):
        recs.append(mk_record(f"A{i}", lab, Direction.DECREASE, es=-30.0,
                              year=2018, pathways=("real",)))
    for i in range(6):
        recs.append(mk_record(f"B{i}", "solo_lab", Direction.DECREASE,
                              es=-40.0, year=2019, pathways=("solo",)))
    for i, lab in enumerate(["l4", "l5"]):
        recs.append(mk_record(f"C{i}", lab, Direction.DECREASE, es=-25.0,
                              year=2019, pathways=("flip",)))
    for i in range(40):
        recs.append(mk_record(f"BG{i}", f"bg{i % 10}", Direction.DECREASE,
                              year=2017, pathways=("background",),
                              phenotype="other"))
        recs[-1].outcome.phenotype = NormalizedTerm(raw="other",
                                                    ontology_id="OTHER:1")
    recs.append(mk_record("F1", "new_lab", Direction.DECREASE, es=-28.0,
                          year=2023, pathways=("real",)))
    recs.append(mk_record("F2", "solo_lab", Direction.DECREASE, es=-41.0,
                          year=2023, pathways=("solo",)))
    recs.append(mk_record("F3", "another_new", Direction.INCREASE, es=20.0,
                          year=2024, pathways=("flip",)))

    res = RetrospectiveBenchmark(min_records_per_phenotype=3).run(recs, 2020)
    by = {c.pathway_label: c for c in res.cases}
    check(by["real"].replicated,
          "confirmation by a new lab counts as independent replication")
    check(not by["solo"].replicated,
          "the SAME lab publishing again does NOT count as replication")
    check(by["flip"].contradicted and not by["flip"].replicated,
          "a new lab reporting the opposite direction counts as contradiction")
    check(res.underpowered and "illustrative" in " ".join(res.caveats),
          "a small benchmark is labelled underpowered")
    check(res.n_records_before == 51 and res.n_records_after == 3,
          "records are split strictly at the cutoff year")

    leak = RetrospectiveBenchmark(min_records_per_phenotype=3).run(recs, 2030)
    check(not leak.cases and "Cannot run" in leak.verdict,
          "a cutoff with nothing after it is refused, not scored")


# --------------------------------------------------------------------------


def test_pooling_and_vocabulary() -> None:
    print("\n[phenotype pooling and vocabulary]")
    check(not commensurable("EPHYS:0001", "EPHYS:0003"),
          "mEPSC and mIPSC frequency are NOT pooled (excitatory vs inhibitory)")
    check(not commensurable("SYN:0001", "MORPH:0003"),
          "synaptic puncta and dendritic spines are NOT pooled")
    check(not commensurable("SYN:0001", "SYN:0005"),
          "excitatory and inhibitory puncta are NOT pooled")
    check(commensurable("SYN:0001", "SYN:0003"),
          "pre-synaptic and colocalised excitatory puncta density ARE pooled")

    m = PhenotypeMatcher()
    want = {
        "weighted mean firing rate": "EPHYS:0006",
        "mfr": "EPHYS:0006",
        "mean firing rate (mfr)": "EPHYS:0006",
        "network burst duration": "EPHYS:0007",
        "mea burst duration": "EPHYS:0007",
        "synchrony index": "EPHYS:0008",
        "spike synchrony": "EPHYS:0008",
        "pairwise synchrony": "EPHYS:0008",
        "correlation index": "EPHYS:0008",
        "resting membrane potential": "EPHYS:0009",
        "input resistance": "EPHYS:0010",
        "rheobase": "EPHYS:0011",
        "VGAT puncta density": "SYN:0005",
        "gephyrin puncta": "SYN:0005",
        "Mander's coefficient": "SYN:0006",
        "mushroom spine density": "MORPH:0005",
        "spine density": "MORPH:0003",
        "network burst frequency": "EPHYS:0004",
        "mitochondrial membrane potential": "CELL:0002",
    }
    for text, tid in want.items():
        got = m.match(text).ontology_id
        check(got == tid, f"'{text}' -> {tid} (got {got})")

    ids = {m.match(t).ontology_id
           for t in ("rheobase", "input resistance",
                     "resting membrane potential")}
    check(len(ids) == 3,
          "rheobase, input resistance and resting potential are separate "
          "terms: their directions mean opposite things for excitability")

    for text in ("NMDA receptor current", "capability index"):
        check(m.match(text).ontology_id is None,
              f"short synonym does not match inside another word: '{text}'")
    check(m.match("NMD target abundance").ontology_id == "SPL:0005",
          "short synonym still matches as a whole word")


def test_spliceai_benchmark() -> None:
    print("\n[sequence: neuronal splice benchmark]")
    from synaptik.sequence.benchmark_spliceai import NEURONAL_SPLICE_BENCHMARK
    check(all(c.ref_seq != c.alt_seq for c in NEURONAL_SPLICE_BENCHMARK),
          "all positive and negative benchmark cases carry genuine sequence substitutions")
    res = run_spliceai_benchmark(threshold=0.20)
    check(res.n_cases == 20, f"benchmark suite evaluates 20 gold-standard cases (got {res.n_cases})")
    check(res.n_positives == 10 and res.n_negatives == 10,
          f"benchmark suite contains 10 true positives and 10 matched controls (got {res.n_positives}/{res.n_negatives})")
    check(res.sensitivity >= 0.80, f"sensitivity >= 80% (got {res.sensitivity:.1%})")
    check(res.specificity >= 0.80, f"specificity >= 80% (got {res.specificity:.1%})")
    check(res.auroc >= 0.90, f"ROC-AUC >= 0.90 (got {res.auroc:.3f})")


def test_heuristic_abstract_extractor() -> None:
    print("\n[extract: deterministic abstract extractor]")
    from synaptik.extract.heuristic import extract_records_from_abstract
    from synaptik.schema import Direction, SourceMetadata

    src = SourceMetadata(
        source_id="PMID:99999999",
        title="TDP-43 loss in iPSC motor neurons",
        year=2023,
        journal="Neuron",
        authors=["Jane Doe", "John Smith"],
        corresponding_author="John Smith",
        affiliations=["Massachusetts General Hospital"],
        lab_id="lab_smith",
    )
    abstract = (
        "Human iPSC-derived NGN2 neurons lacking TARDBP showed markedly decreased "
        "presynaptic puncta density and impaired synaptic transmission. "
        "Conversely, cryptic exon inclusion increased in STMN2 splicing assays."
    )
    recs = extract_records_from_abstract(src, abstract)
    check(len(recs) >= 1, f"extracted sentence-grounded records from abstract (got {len(recs)})")
    check(any(r.perturbation.target and r.perturbation.target.ontology_label == "TARDBP" and r.outcome.direction == Direction.DECREASE for r in recs),
          "captured TARDBP decrease direction")
    check(all(r.provenance.quote in abstract for r in recs), "every record carries a verbatim abstract quote")


def test_next_stage_upgrades() -> None:
    print("\n[next-stage: pathway mapping, gene-phenotype benchmark, extraction F1, accessibility, audit]")
    from synaptik.cli.commands.audit import _build_audit_dict
    from synaptik.drug.aso_calibrated import compute_target_accessibility
    from synaptik.extract.benchmark_extraction import run_extraction_benchmark
    from synaptik.extract.fulltext import extract_records_from_fulltext
    from synaptik.normalize.phenotype_ontology import TERMS, canonical_pathways_for_phenotype
    from synaptik.schema import Direction, SourceMetadata

    # 1. Every PhenotypeTerm maps to >=1 canonical GO/Reactome pathway
    unmapped = [t.term_id for t in TERMS if not canonical_pathways_for_phenotype(t.term_id)]
    check(not unmapped, f"all PhenotypeTerms map to canonical GO/Reactome pathways (missing: {unmapped})")

    # 2. Fulltext extractor populates pathways, assay_class, and gene-phenotype benchmark works
    recs = []
    for i, (yr, lab) in enumerate([(2019, "lab_a"), (2020, "lab_b"), (2023, "lab_c")]):
        src = SourceMetadata(source_id=f"PMID:{1000+i}", year=yr, lab_id=lab)
        secs = {
            "methods": "NGN2 cortical neurons from 3 ALS donor lines and isogenic controls (2 clones) at DIV28.",
            "results": "TARDBP knockdown significantly decreased presynaptic puncta density.",
        }
        recs.extend(extract_records_from_fulltext(src, secs))
    check(len(recs) == 3 and all(r.pathways for r in recs), "fulltext extractor populates canonical pathways")
    gp_res = RetrospectiveBenchmark().run_gene_phenotype(recs, cutoff_year=2020, min_records=2)
    check(len(gp_res.cases) == 1 and gp_res.cases[0].replicated, "gene-phenotype temporal replication tracks new-lab confirmation")

    # 3. Gold-standard extraction benchmark
    ext_res = run_extraction_benchmark()
    check(all(m.f1 >= 0.85 for m in ext_res.metrics), "all extraction benchmark fields achieve F1 >= 0.85")

    # 4. Target mRNA accessibility features
    acc = compute_target_accessibility("GCGCGCGCGCGCGCGCGCGC", "ATATATATATATATATATAT", "GCGCGCGCGCGCGCGCGCGC")
    check(len(acc) == 10 and acc[0] == 1.0 and 0.0 <= acc[9] <= 1.0, "target accessibility returns bounded 10-element feature vector")

    # 5. Audit dict handles missing DB cleanly
    audit = _build_audit_dict("STMN2", None, "nonexistent_audit_test.db", None, "gapmer", live=False)
    check(audit["gene"] == "STMN2" and audit["caveats"], "audit handles missing DB gracefully with caveat")


def main() -> int:
    print("=" * 70)
    print("SYNAPTIK lab-workflow regression suite")
    print("=" * 70)
    test_splicing_table_parsing()
    test_splicing_reconciliation()
    test_frame_consequence()
    test_phenotype_table()
    test_power()
    test_design_brief()
    test_benchmark()
    test_pooling_and_vocabulary()
    test_spliceai_benchmark()
    test_heuristic_abstract_extractor()
    test_next_stage_upgrades()
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
