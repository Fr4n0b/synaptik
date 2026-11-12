"""
Synthetic corpus with known ground truth.

Because we cannot validate the statistics against the real literature without
network access, we generate a corpus where the right answer is known by
construction, and assert the engine recovers it. Three planted signals:

  TRUE POSITIVE   -- 'synaptic vesicle cycle' is genuinely associated with
                     presynaptic puncta loss: many labs, multiple assay
                     families, consistent direction. Must rank highly and be
                     tiered 'strong'.

  SINGLE-LAB TRAP -- 'exotic pathway X' appears in many papers, all from one
                     lab, one assay. Must NOT be tiered 'strong', must fail
                     leave-one-lab-out, and must carry the single-group caveat.

  UBIQUITOUS      -- 'apoptosis' appears everywhere including the background.
                     Must NOT be enriched despite a high raw count.

If the engine cannot separate these three, it is not fit for purpose.
"""

from __future__ import annotations

import random

from synaptik.schema import (
    AssayClass, CellModel, CellModelClass, Direction, EvidenceType,
    ExperimentRecord, NormalizedTerm, Outcome, Perturbation, PerturbationType,
    Provenance, SourceMetadata,
)

JOURNALS_Q1 = ["Neuron", "Nature Neuroscience", "Cell Stem Cell", "eLife"]
JOURNALS_Q2 = ["Human Molecular Genetics", "Scientific Reports", "Molecular Brain"]
JOURNALS_Q3 = ["Neuroscience Letters", "Brain Research"]


def _src(i, lab, journal, year, quartile, retracted=False, preprint=False):
    return SourceMetadata(
        source_id=f"PMID:{i}",
        pmid=str(i),
        title=f"Study {i} of neuronal phenotypes",
        journal=journal,
        year=year,
        sjr_quartile=quartile,
        lab_id=lab,
        is_retracted=retracted,
        is_preprint=preprint,
        authors=[f"Author{i}", f"PI_{lab}"],
    )


def _rec(src, phenotype, direction, pathways, assay, n_rep=4, es=None, p=0.01,
         isogenic=True, ev=EvidenceType.DIRECT_MEASUREMENT, conf=0.85):
    return ExperimentRecord(
        source=src,
        cell_model=CellModel(
            model_class=CellModelClass.IPSC_NEURON_NGN2,
            days_in_vitro=21,
            is_isogenic_pair=isogenic,
            differentiation_protocol="doxycycline-inducible NGN2",
        ),
        perturbation=Perturbation(
            ptype=PerturbationType.MUTATION_PATIENT,
            target=NormalizedTerm(raw="TESTGENE"),
        ),
        outcome=Outcome(
            phenotype=NormalizedTerm(raw=phenotype, ontology_id=phenotype),
            direction=direction,
            effect_size=es,
            effect_size_unit="percent_change" if es is not None else None,
            p_value=p,
            n_biological_replicates=n_rep,
            is_quantified=es is not None,
        ),
        assay_class=assay,
        assay_raw=assay.value,
        evidence_type=ev,
        pathways=[NormalizedTerm(raw=x, ontology_id=x) for x in pathways],
        provenance=Provenance(quote="synthetic supporting span for validation",
                              extraction_confidence=conf),
    )


def build_corpus(seed: int = 7) -> tuple[list[ExperimentRecord], list[ExperimentRecord]]:
    """Returns (foreground, background). Foreground = presynaptic puncta records."""
    rng = random.Random(seed)
    fg: list[ExperimentRecord] = []
    bg: list[ExperimentRecord] = []
    i = 1000

    # --- TRUE POSITIVE: synaptic vesicle cycle, 6 labs, 3 assay families ---
    true_labs = [f"lab_true_{k}" for k in range(6)]
    assays = [
        AssayClass.IMAGING_SYNAPTIC_PUNCTA,
        AssayClass.ELECTROPHYS_PATCH,
        AssayClass.PROTEOMICS,
    ]
    for k, lab in enumerate(true_labs):
        for j in range(2):
            i += 1
            q = "Q1" if k < 4 else "Q2"
            journal = rng.choice(JOURNALS_Q1 if q == "Q1" else JOURNALS_Q2)
            fg.append(_rec(
                _src(i, lab, journal, 2016 + k, q),
                "SYN:0001", Direction.DECREASE,
                ["synaptic vesicle cycle", "apoptosis"],
                assays[(k + j) % 3],
                es=-30.0 - k, p=0.004,
            ))
    # one honest null result for the true pathway
    i += 1
    fg.append(_rec(
        _src(i, "lab_true_9", "eLife", 2022, "Q1"),
        "SYN:0001", Direction.NO_CHANGE,
        ["synaptic vesicle cycle"], AssayClass.IMAGING_SYNAPTIC_PUNCTA,
        es=None, p=0.61,
    ))

    # --- SINGLE-LAB TRAP: one lab, 8 papers, one assay, always positive ---
    for j in range(8):
        i += 1
        fg.append(_rec(
            _src(i, "lab_prolific", rng.choice(JOURNALS_Q1), 2018 + j % 5, "Q1"),
            "SYN:0001", Direction.DECREASE,
            ["exotic pathway X", "apoptosis"],
            AssayClass.IMAGING_SYNAPTIC_PUNCTA,
            es=-40.0, p=0.001, n_rep=3,
        ))

    # --- CONFLICTING: reported in both directions by different labs ---
    for j, (lab, d) in enumerate(
        [("lab_c1", Direction.DECREASE), ("lab_c2", Direction.INCREASE),
         ("lab_c3", Direction.INCREASE), ("lab_c4", Direction.DECREASE)]
    ):
        i += 1
        fg.append(_rec(
            _src(i, lab, rng.choice(JOURNALS_Q2), 2019 + j, "Q2"),
            "SYN:0001", d, ["contested pathway Y"],
            AssayClass.IMAGING_SYNAPTIC_PUNCTA, es=20.0, p=0.03,
        ))

    # --- RETRACTED paper supporting a bogus pathway ---
    i += 1
    fg.append(_rec(
        _src(i, "lab_bad", "Nature Neuroscience", 2015, "Q1", retracted=True),
        "SYN:0001", Direction.DECREASE, ["fraudulent pathway Z"],
        AssayClass.IMAGING_SYNAPTIC_PUNCTA, es=-70.0, p=0.0001,
    ))
    i += 1
    fg.append(_rec(
        _src(i, "lab_bad2", "Brain Research", 2016, "Q3"),
        "SYN:0001", Direction.DECREASE, ["fraudulent pathway Z"],
        AssayClass.IMAGING_SYNAPTIC_PUNCTA, es=-65.0, p=0.001,
    ))

    # --- BACKGROUND: apoptosis everywhere, vesicle cycle rare ---
    for j in range(300):
        i += 1
        lab = f"lab_bg_{j % 40}"
        pathways = ["apoptosis"]
        if j % 25 == 0:
            pathways.append("synaptic vesicle cycle")
        if j % 30 == 0:
            pathways.append("contested pathway Y")
        bg.append(_rec(
            _src(i, lab, rng.choice(JOURNALS_Q1 + JOURNALS_Q2 + JOURNALS_Q3),
                 2015 + j % 10, rng.choice(["Q1", "Q2", "Q3"])),
            "OTHER:0001", Direction.DECREASE, pathways,
            rng.choice(list(AssayClass)), es=-10.0, p=0.04,
        ))

    return fg, bg
