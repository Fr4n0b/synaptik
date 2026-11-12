"""
Builds a realistic demo knowledge base.

This is not a toy: it models the actual structure of the iPSC-neuron
literature, including the pathologies that make naive analysis fail --
one dominant lab, publication bias toward positive findings, a retracted
paper, conflicting reports, and phenotypes described in inconsistent language.

Used for the end-to-end integration test and as a sandbox the user can explore
before harvesting real data.
"""

from __future__ import annotations

import random

from synaptik.pipeline import Pipeline
from synaptik.schema import (
    AssayClass, CellModel, CellModelClass, Direction, EvidenceType,
    ExperimentRecord, NormalizedTerm, Outcome, Perturbation, PerturbationType,
    Provenance, SourceMetadata, SplicingDetail, SplicingEventType,
)
from synaptik.store import Store

JOURNALS = {
    "Q1": ["Neuron", "Nature Neuroscience", "Cell Stem Cell", "eLife",
           "Molecular Psychiatry", "Brain", "Nature Communications"],
    "Q2": ["Human Molecular Genetics", "Scientific Reports", "Molecular Brain",
           "Journal of Neurochemistry"],
    "Q3": ["Neuroscience Letters", "Brain Research"],
}

INSTITUTIONS = [
    "Massachusetts General Hospital, Boston, MA",
    "Broad Institute, Cambridge, MA",
    "Stanford University, Stanford, CA",
    "UCSF, San Francisco, CA",
    "Columbia University, New York, NY",
    "Karolinska Institutet, Stockholm, Sweden",
    "UCL Institute of Neurology, London, UK",
    "RIKEN Center for Brain Science, Wako, Japan",
    "Johns Hopkins University, Baltimore, MD",
    "Cedars-Sinai Medical Center, Los Angeles, CA",
]

SURNAMES = ["Chen", "Okafor", "Nakamura", "Petrov", "Silva", "Ahmed",
            "Lindqvist", "Rossi", "Dubois", "Kowalski", "Tanaka", "Mbeki"]

# Phrasings the same phenotype appears under, to exercise normalization.
PHENOTYPE_PHRASINGS = {
    "SYN:0003": [
        "colocalized SYN1/PSD95 puncta density", "synapse density",
        "synaptic puncta density", "number of synapses",
        "co-localized pre- and postsynaptic puncta",
    ],
    "SYN:0001": [
        "SYN1 puncta density", "synapsin puncta", "presynaptic puncta density",
        "VGLUT1 puncta density",
    ],
    "EPHYS:0001": ["mEPSC frequency", "miniature EPSC frequency",
                   "spontaneous EPSC frequency"],
    "EPHYS:0002": ["mEPSC amplitude", "miniature EPSC amplitude"],
    "MORPH:0002": ["dendritic branching", "Sholl analysis", "arborization"],
    "SPL:0003": ["cryptic exon inclusion", "cryptic splicing"],
    "SPL:0002": ["intron retention", "retained intron"],
    "CELL:0001": ["neuronal viability", "cell survival", "cell death"],
}


def build_demo_store(path: str = "data/demo.db", seed: int = 42) -> Store:
    rng = random.Random(seed)
    store = Store(path)
    pipe = Pipeline(store)

    sources: list[SourceMetadata] = []
    records: list[ExperimentRecord] = []
    pmid = 30000000

    def mk_source(lab_idx, year, quartile, retracted=False, preprint=False,
                  review=False):
        nonlocal pmid
        pmid += rng.randint(1, 400)
        inst = INSTITUTIONS[lab_idx % len(INSTITUTIONS)]
        pi = f"{rng.choice(SURNAMES)} {chr(65 + lab_idx % 26)}"
        first = f"{rng.choice(SURNAMES)} {chr(65 + rng.randint(0,25))}"
        return SourceMetadata(
            source_id=f"PMID:{pmid}", pmid=str(pmid),
            title=f"Study of neuronal phenotypes in patient iPSC neurons ({pmid})",
            journal=rng.choice(JOURNALS[quartile]), year=year,
            sjr_quartile=quartile, is_retracted=retracted, is_preprint=preprint,
            is_review=review,
            authors=[first, f"Middle {chr(65+rng.randint(0,25))}", pi],
            affiliations=[inst], corresponding_author=pi,
            open_access=True, full_text_available=True,
            citation_count=rng.randint(0, 300),
        )

    def mk_record(src, pheno_id, direction, pathways, assay, **kw):
        phrasing = rng.choice(PHENOTYPE_PHRASINGS.get(pheno_id, [pheno_id]))
        es = kw.get("effect_size")
        return ExperimentRecord(
            source=src,
            cell_model=CellModel(
                model_class=kw.get("model", CellModelClass.IPSC_NEURON_NGN2),
                days_in_vitro=kw.get("div", rng.choice([21, 21, 24, 28])),
                is_isogenic_pair=kw.get("isogenic", True),
                differentiation_protocol="doxycycline-inducible NGN2",
                coculture=kw.get("coculture", "mouse astrocytes"),
                species="human",
            ),
            perturbation=Perturbation(
                ptype=kw.get("ptype", PerturbationType.MUTATION_PATIENT),
                target=NormalizedTerm(raw=kw.get("target", "TARDBP")),
                agent=kw.get("agent"),
                variant_hgvs=kw.get("variant"),
                control_description="isogenic corrected control",
            ),
            outcome=Outcome(
                phenotype=NormalizedTerm(raw=phrasing),
                direction=direction,
                effect_size=es,
                effect_size_unit="percent_change" if es is not None else None,
                p_value=kw.get("p", round(rng.uniform(0.0005, 0.045), 4)),
                n_biological_replicates=kw.get("n", rng.choice([3, 3, 4, 4, 6])),
                n_cells=rng.randint(30, 400),
                statistical_test=rng.choice(["t-test", "one-way ANOVA",
                                             "Mann-Whitney U"]),
                markers_used=kw.get("markers", ["SYN1", "PSD95", "MAP2"]),
                is_quantified=es is not None,
            ),
            assay_class=assay,
            assay_raw=kw.get("assay_raw", assay.value.replace("_", " ")),
            evidence_type=kw.get("ev", EvidenceType.DIRECT_MEASUREMENT),
            pathways=[NormalizedTerm(raw=p) for p in pathways],
            genes=[NormalizedTerm(raw=g) for g in kw.get("genes", [])],
            splicing=kw.get("splicing"),
            provenance=Provenance(
                section="results",
                quote="synthetic supporting sentence for the demo corpus",
                extraction_confidence=round(rng.uniform(0.7, 0.95), 2),
            ),
        )

    # ---- SIGNAL 1: synaptic vesicle cycle <-> presynaptic puncta loss ----
    # Well replicated: 7 labs, 3 assay families, one honest null.
    for lab in range(7):
        for rep in range(rng.randint(1, 3)):
            q = "Q1" if lab < 4 else rng.choice(["Q1", "Q2"])
            src = mk_source(lab, 2015 + lab + rep, q)
            sources.append(src)
            assay = [AssayClass.IMAGING_SYNAPTIC_PUNCTA,
                     AssayClass.ELECTROPHYS_PATCH,
                     AssayClass.PROTEOMICS][(lab + rep) % 3]
            pheno = ("SYN:0001" if assay == AssayClass.IMAGING_SYNAPTIC_PUNCTA
                     else "EPHYS:0001" if assay == AssayClass.ELECTROPHYS_PATCH
                     else "SYN:0001")
            records.append(mk_record(
                src, pheno, Direction.DECREASE,
                ["synaptic vesicle cycle", "neurotransmitter release",
                 "apoptosis"],
                assay, effect_size=-rng.uniform(20, 45),
                genes=["SYN1", "SNAP25", "VAMP2"],
            ))
    src = mk_source(11, 2023, "Q1")
    sources.append(src)
    records.append(mk_record(
        src, "SYN:0001", Direction.NO_CHANGE, ["synaptic vesicle cycle"],
        AssayClass.IMAGING_SYNAPTIC_PUNCTA, effect_size=None, p=0.72,
    ))

    # ---- TRAP: one prolific lab, one assay, always positive ----
    for j in range(9):
        src = mk_source(3, 2017 + j % 6, "Q1")  # same institution index = same lab
        src.authors[-1] = "Prolific P"
        src.corresponding_author = "Prolific P"
        sources.append(src)
        records.append(mk_record(
            src, "SYN:0001", Direction.DECREASE,
            ["exotic kinase cascade", "apoptosis"],
            AssayClass.IMAGING_SYNAPTIC_PUNCTA,
            effect_size=-rng.uniform(35, 50), n=3, p=0.002,
        ))

    # ---- CONFLICT: mTOR reported both ways ----
    for j, d in enumerate([Direction.DECREASE, Direction.INCREASE,
                           Direction.INCREASE, Direction.DECREASE,
                           Direction.NO_CHANGE]):
        src = mk_source(4 + j, 2018 + j, rng.choice(["Q1", "Q2"]))
        sources.append(src)
        records.append(mk_record(
            src, "SYN:0001", d, ["mTOR signaling", "autophagy"],
            AssayClass.IMAGING_SYNAPTIC_PUNCTA,
            effect_size=(None if d == Direction.NO_CHANGE
                         else rng.uniform(15, 30) * (1 if d == Direction.INCREASE else -1)),
        ))

    # ---- RETRACTED ----
    src = mk_source(9, 2014, "Q1", retracted=True)
    sources.append(src)
    records.append(mk_record(
        src, "SYN:0001", Direction.DECREASE, ["fabricated pathway"],
        AssayClass.IMAGING_SYNAPTIC_PUNCTA, effect_size=-78.0, p=0.00001,
    ))

    # ---- SPLICING: TDP-43 cryptic exon signature across labs ----
    cryptic_genes = ["STMN2", "UNC13A", "KALRN", "AGRN", "POLDIP3"]
    for j, gene in enumerate(cryptic_genes):
        for lab_off in range(2 if j < 3 else 1):
            src = mk_source(j + lab_off, 2019 + j, "Q1")
            sources.append(src)
            records.append(mk_record(
                src, "SPL:0003", Direction.INCREASE,
                ["RNA splicing", "nonsense-mediated decay"],
                AssayClass.SPLICING_ANALYSIS,
                assay_raw="rMATS analysis of RNA-seq",
                target="TARDBP", genes=[gene, "TARDBP"],
                effect_size=None,
                splicing=SplicingDetail(
                    event_type=SplicingEventType.SKIPPED_EXON,
                    gene=NormalizedTerm(raw=gene),
                    delta_psi=round(rng.uniform(0.25, 0.55), 2),
                    detection_method="rMATS",
                    regulator=NormalizedTerm(raw="TDP-43"),
                    genome_build="GRCh38",
                ),
            ))

    # ---- SPLICING: intron retention, different regulator ----
    for j in range(3):
        src = mk_source(6 + j, 2020 + j, "Q2")
        sources.append(src)
        records.append(mk_record(
            src, "SPL:0002", Direction.INCREASE, ["RNA splicing"],
            AssayClass.SPLICING_ANALYSIS, assay_raw="MAJIQ on RNA-seq",
            target="SRSF1", genes=["SRSF1"], effect_size=None,
            splicing=SplicingDetail(
                event_type=SplicingEventType.RETAINED_INTRON,
                gene=NormalizedTerm(raw=f"GENE{j}"),
                delta_psi=round(rng.uniform(0.15, 0.35), 2),
                detection_method="MAJIQ",
                regulator=NormalizedTerm(raw="SRSF1"),
            ),
        ))

    # ---- DRUG RESCUE experiments ----
    for agent, path, eff in [("rapamycin", "mTOR signaling", 25.0),
                             ("TG003", "RNA splicing", 18.0),
                             ("BDNF", "BDNF signaling", 30.0)]:
        for lab_off in range(2):
            src = mk_source(2 + lab_off, 2021 + lab_off, "Q1")
            sources.append(src)
            records.append(mk_record(
                src, "SYN:0001", Direction.INCREASE, [path],
                AssayClass.IMAGING_SYNAPTIC_PUNCTA,
                ptype=PerturbationType.SMALL_MOLECULE, agent=agent,
                target=path.split()[0], effect_size=eff,
                ev=EvidenceType.RESCUE,
            ))

    # ---- BACKGROUND: broad corpus so enrichment has a real denominator ----
    other_phenos = ["MORPH:0002", "CELL:0001", "EPHYS:0002", "MORPH:0001"]
    common_paths = ["apoptosis", "cell cycle", "oxidative stress",
                    "protein synthesis", "calcium signaling",
                    "cytoskeletal organization", "mitochondrial function"]
    for j in range(260):
        lab = j % 40 + 12
        src = mk_source(lab, 2014 + j % 12,
                        rng.choice(["Q1", "Q2", "Q2", "Q3"]),
                        preprint=(j % 37 == 0))
        sources.append(src)
        paths = rng.sample(common_paths, rng.randint(1, 3))
        if j % 22 == 0:
            paths.append("synaptic vesicle cycle")
        if j % 19 == 0:
            paths.append("mTOR signaling")
        records.append(mk_record(
            src, rng.choice(other_phenos),
            rng.choice([Direction.DECREASE, Direction.INCREASE,
                        Direction.NO_CHANGE]),
            paths, rng.choice(list(AssayClass)),
            effect_size=(None if rng.random() < 0.3
                         else round(rng.uniform(-40, 40), 1)),
            model=rng.choice([
                CellModelClass.IPSC_NEURON_NGN2,
                CellModelClass.IPSC_NEURON_NGN2,
                CellModelClass.IPSC_NEURON_SMALL_MOLECULE,
                CellModelClass.PRIMARY_RODENT_NEURON,
            ]),
        ))

    # ---- persist through the real normalization path ----
    pipe.labs.fit(sources)
    with store.tx():
        for s in sources:
            s.lab_id = pipe.labs.assign(s)
            pipe.journals.annotate(s)
            store.upsert_source(s)
            store.index_source_text(s.source_id, s.title, "")

    normalized = []
    for r in records:
        r.source.lab_id = pipe.labs.assign(r.source)
        normalized.append(pipe.normalize_record(r))
    store.insert_records(normalized)

    return store


if __name__ == "__main__":
    import sys

    path = sys.argv[1] if len(sys.argv) > 1 else "data/demo.db"
    s = build_demo_store(path)
    print(f"Demo store built at {path}")
    for k, v in s.stats().items():
        print(f"  {k}: {v}")
