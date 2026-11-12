"""
Controlled vocabulary for iPSC-neuron phenotypes.

Why a custom vocabulary: GO, HPO, and CL each cover part of this space but
none covers "SYN1 puncta density measured by confocal in NGN2 neurons" as a
single addressable concept. GO has "synapse organization" (a process, not a
measurement); HPO is clinical. We therefore define a measurement-level
vocabulary and cross-reference it to GO/HPO where a mapping exists.

Each term carries:
  - synonyms: surface forms seen in the literature, used for matching
  - go_xref / hpo_xref: cross-references for downstream pathway joins
  - measurement_axis: what is physically being counted, so that two terms on
    the same axis can be recognized as commensurable
  - polarity_hint: whether "increase" is usually reported as a gain or loss of
    function, used only for display, never for inference

Extending this file is the main way to adapt the tool to a new phenotype
domain. Add a term, re-run normalization, done.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class PhenotypeTerm:
    term_id: str
    label: str
    synonyms: tuple[str, ...] = ()
    measurement_axis: str = ""
    # Terms sharing a pool_group are pooled when analysing a phenotype. This
    # is an explicit, curated statement that they measure the same biological
    # quantity. It is deliberately NOT derived from measurement_axis: mEPSC
    # and mIPSC frequency are both event rates, and are not the same thing.
    pool_group: str = ""
    go_xref: tuple[str, ...] = ()
    hpo_xref: tuple[str, ...] = ()
    typical_assays: tuple[str, ...] = ()
    parent: str | None = None
    notes: str = ""


# --------------------------------------------------------------------------
# Synaptic phenotypes -- the core of the user's use case
# --------------------------------------------------------------------------

TERMS: list[PhenotypeTerm] = [
    PhenotypeTerm(
        "SYN:0001", "presynaptic puncta density",
        (
            "synapsin puncta", "syn1 puncta", "synapsin-1 puncta density",
            "presynaptic puncta", "vglut1 puncta", "vglut1 puncta density",
            "synaptophysin puncta", "syp puncta", "bassoon puncta",
            "presynaptic terminal density", "presynaptic marker density",
            "synaptic vesicle protein puncta",
        ),
        measurement_axis="puncta_count_per_length",
        pool_group="excitatory_synapse_number",
        go_xref=("GO:0048786",),  # presynaptic active zone
        typical_assays=("imaging_synaptic_puncta",),
        parent="SYN:0000",
    ),
    PhenotypeTerm(
        "SYN:0002", "postsynaptic puncta density",
        (
            "psd95 puncta", "psd-95 puncta density", "homer1 puncta",
            "shank3 puncta", "postsynaptic density puncta",
            "postsynaptic marker density", "dlg4 puncta", "homer puncta density",
        ),
        measurement_axis="puncta_count_per_length",
        pool_group="excitatory_synapse_number",
        go_xref=("GO:0014069",),  # postsynaptic density
        typical_assays=("imaging_synaptic_puncta",),
        parent="SYN:0000",
    ),
    PhenotypeTerm(
        "SYN:0003", "colocalized synaptic puncta density",
        (
            "colocalized puncta", "co-localized synaptic puncta",
            "syn1/psd95 colocalization", "pre-post colocalization",
            "synaptic puncta density", "synapse number", "synapse density",
            "synapse count", "number of synapses", "synaptic density",
            "colocalized pre- and postsynaptic puncta",
        ),
        measurement_axis="puncta_count_per_length",
        pool_group="excitatory_synapse_number",
        go_xref=("GO:0045202",),  # synapse
        typical_assays=("imaging_synaptic_puncta",),
        parent="SYN:0000",
        notes="Most rigorous imaging readout of synapse number; requires both markers.",
    ),
    PhenotypeTerm(
        "SYN:0004", "synaptic puncta size",
        ("puncta size", "puncta area", "synaptic puncta intensity",
         "puncta volume", "punctum area", "synapse size"),
        measurement_axis="puncta_morphometry",
        typical_assays=("imaging_synaptic_puncta",),
        parent="SYN:0000",
    ),
    PhenotypeTerm(
        "SYN:0005", "inhibitory synaptic puncta density",
        ("vgat puncta", "vgat puncta density", "gephyrin puncta",
         "gephyrin puncta density", "gad65 puncta", "gad67 puncta",
         "inhibitory synapse density", "inhibitory synapse number",
         "gabaergic synapse density", "inhibitory puncta"),
        measurement_axis="puncta_count_per_length",
        go_xref=("GO:0060077",),  # inhibitory synapse
        typical_assays=("imaging_synaptic_puncta",),
        parent="SYN:0000",
        notes="Kept separate from excitatory puncta: the two can move in "
              "opposite directions and pooling them hides E/I imbalance.",
    ),
    PhenotypeTerm(
        "SYN:0006", "synaptic marker colocalization coefficient",
        ("manders coefficient", "mander's coefficient", "manders overlap",
         "manders colocalization", "pearson colocalization",
         "pearson's correlation coefficient", "colocalization coefficient",
         "colocalisation coefficient", "overlap coefficient"),
        measurement_axis="colocalization_coefficient",
        typical_assays=("imaging_synaptic_puncta",),
        parent="SYN:0000",
        notes="An intensity-correlation measure, not a synapse count; "
              "sensitive to background and thresholding.",
    ),
    PhenotypeTerm(
        "SYN:0000", "synaptic phenotype", ("synaptic abnormality",),
        measurement_axis="composite",
        go_xref=("GO:0045202",),
    ),
    # ---- electrophysiology ----
    PhenotypeTerm(
        "EPHYS:0001", "mEPSC frequency",
        ("mepsc frequency", "miniature epsc frequency", "spontaneous epsc frequency",
         "sepsc frequency", "mini frequency", "excitatory event frequency"),
        measurement_axis="event_rate",
        go_xref=("GO:0060079",),
        typical_assays=("electrophys_patch",),
        parent="EPHYS:0000",
        notes="Classically interpreted as presynaptic release probability or synapse number.",
    ),
    PhenotypeTerm(
        "EPHYS:0002", "mEPSC amplitude",
        ("mepsc amplitude", "miniature epsc amplitude", "sepsc amplitude",
         "mini amplitude", "quantal amplitude"),
        measurement_axis="event_amplitude",
        typical_assays=("electrophys_patch",),
        parent="EPHYS:0000",
        notes="Classically interpreted as postsynaptic receptor content.",
    ),
    PhenotypeTerm(
        "EPHYS:0003", "mIPSC frequency",
        ("mipsc frequency", "miniature ipsc frequency", "sipsc frequency",
         "inhibitory event frequency"),
        measurement_axis="event_rate",
        typical_assays=("electrophys_patch",),
        parent="EPHYS:0000",
    ),
    # ---- multielectrode array ----
    PhenotypeTerm(
        "EPHYS:0004", "network burst frequency",
        ("network burst frequency", "network burst rate", "burst frequency",
         "bursting frequency", "mea burst frequency", "burst rate",
         "number of network bursts", "network bursts per minute"),
        measurement_axis="event_rate",
        typical_assays=("electrophys_mea",),
        parent="EPHYS:0000",
    ),
    PhenotypeTerm(
        "EPHYS:0006", "mean firing rate",
        ("mean firing rate", "weighted mean firing rate", "firing rate",
         "spike rate", "spike frequency", "spontaneous firing rate",
         "spontaneous spiking"),
        measurement_axis="event_rate",
        typical_assays=("electrophys_mea",),
        parent="EPHYS:0000",
        notes="Depends strongly on plating density and age; compare the "
              "same DIV and normalise to active electrodes.",
    ),
    PhenotypeTerm(
        "EPHYS:0007", "burst duration",
        ("burst duration", "network burst duration", "mean burst duration",
         "burst length"),
        measurement_axis="duration",
        typical_assays=("electrophys_mea",),
        parent="EPHYS:0000",
    ),
    PhenotypeTerm(
        "EPHYS:0008", "network synchrony",
        ("synchrony index", "network synchrony", "synchronicity",
         "synchrony", "area under cross-correlation",
         "area under normalized cross-correlation", "network synchronization"),
        measurement_axis="synchrony_index",
        typical_assays=("electrophys_mea", "calcium_imaging"),
        parent="EPHYS:0000",
    ),
    # ---- intrinsic properties (patch-clamp) ----
    # These are separate terms on purpose. They were once pooled as
    # "intrinsic excitability", but their directions mean opposite things: a
    # HIGHER rheobase is LOWER excitability, a higher input resistance is
    # higher excitability. Pooling them makes direction consistency
    # meaningless.
    PhenotypeTerm(
        "EPHYS:0005", "intrinsic excitability",
        ("intrinsic excitability", "excitability", "ap firing",
         "action potential firing", "evoked firing",
         "evoked action potentials", "number of action potentials",
         "f-i curve", "firing frequency upon current injection"),
        measurement_axis="evoked_firing",
        typical_assays=("electrophys_patch",),
        parent="EPHYS:0000",
    ),
    PhenotypeTerm(
        "EPHYS:0009", "resting membrane potential",
        ("resting membrane potential", "resting potential",
         "membrane potential"),
        measurement_axis="voltage",
        typical_assays=("electrophys_patch",),
        parent="EPHYS:0000",
        notes="'Increase' is ambiguous for a negative quantity; check "
              "whether a paper means depolarised or hyperpolarised.",
    ),
    PhenotypeTerm(
        "EPHYS:0010", "input resistance",
        ("input resistance", "membrane resistance"),
        measurement_axis="resistance",
        typical_assays=("electrophys_patch",),
        parent="EPHYS:0000",
        notes="Falls as neurons mature; a difference between lines can "
              "reflect a maturation offset.",
    ),
    PhenotypeTerm(
        "EPHYS:0011", "rheobase",
        ("rheobase", "rheobase current", "current threshold",
         "minimum current to elicit an action potential"),
        measurement_axis="current",
        typical_assays=("electrophys_patch",),
        parent="EPHYS:0000",
        notes="Higher rheobase means LOWER excitability.",
    ),
    PhenotypeTerm(
        "EPHYS:0012", "action potential threshold",
        ("action potential threshold", "ap threshold", "spike threshold",
         "voltage threshold"),
        measurement_axis="voltage",
        typical_assays=("electrophys_patch",),
        parent="EPHYS:0000",
    ),
    PhenotypeTerm(
        "EPHYS:0013", "mIPSC amplitude",
        ("mipsc amplitude", "miniature ipsc amplitude", "sipsc amplitude"),
        measurement_axis="event_amplitude",
        typical_assays=("electrophys_patch",),
        parent="EPHYS:0000",
    ),
    PhenotypeTerm(
        "EPHYS:0014", "membrane capacitance",
        ("membrane capacitance", "cell capacitance", "whole-cell capacitance"),
        measurement_axis="capacitance",
        typical_assays=("electrophys_patch",),
        parent="EPHYS:0000",
        notes="A proxy for cell size; changes confound current-density "
              "comparisons.",
    ),
    PhenotypeTerm(
        "EPHYS:0000", "electrophysiological phenotype", ("ephys phenotype",),
        measurement_axis="composite",
    ),
    # ---- morphology ----
    PhenotypeTerm(
        "MORPH:0001", "neurite length",
        ("total neurite length", "neurite outgrowth", "axon length",
         "process length", "neurite extension", "total process length"),
        measurement_axis="length",
        go_xref=("GO:0031175",),
        typical_assays=("imaging_morphology",),
        parent="MORPH:0000",
    ),
    PhenotypeTerm(
        "MORPH:0002", "dendritic branching",
        ("branching complexity", "sholl analysis", "dendritic complexity",
         "branch points", "arborization", "neurite branching",
         "number of branch points", "dendritic arborization"),
        measurement_axis="branch_count",
        go_xref=("GO:0016358",),
        typical_assays=("imaging_morphology",),
        parent="MORPH:0000",
    ),
    PhenotypeTerm(
        "MORPH:0003", "dendritic spine density",
        ("spine density", "spine number", "dendritic spines",
         "dendritic spine number", "spines per 10"),
        measurement_axis="puncta_count_per_length",
        go_xref=("GO:0060996",),
        typical_assays=("imaging_morphology",),
        parent="MORPH:0000",
    ),
    PhenotypeTerm(
        "MORPH:0005", "dendritic spine morphology",
        ("spine morphology", "spine morphometry", "mushroom spines",
         "mushroom spine density", "mushroom spine fraction",
         "stubby spines", "thin spines", "filopodia", "spine head width",
         "spine head diameter", "spine length", "spine maturity",
         "spine classification"),
        measurement_axis="spine_morphometry",
        go_xref=("GO:0060997",),  # dendritic spine morphogenesis
        typical_assays=("imaging_morphology",),
        parent="MORPH:0000",
        notes="Subtype distribution, not number. Classic spine classes are "
              "rare in iPSC-derived neurons before late maturity.",
    ),
    PhenotypeTerm(
        "MORPH:0004", "soma size",
        ("cell body size", "soma area", "somatic size", "cell size"),
        measurement_axis="area",
        typical_assays=("imaging_morphology",),
        parent="MORPH:0000",
    ),
    PhenotypeTerm(
        "MORPH:0000", "morphological phenotype", ("morphology defect",),
        measurement_axis="composite",
    ),
    # ---- survival / stress ----
    PhenotypeTerm(
        "CELL:0001", "neuronal viability",
        ("cell viability", "cell survival", "neuronal survival", "cell death",
         "apoptosis", "caspase activation", "neurodegeneration", "cytotoxicity",
         "ldh release", "propidium iodide"),
        measurement_axis="fraction_viable",
        go_xref=("GO:0043524",),
        typical_assays=("viability",),
    ),
    PhenotypeTerm(
        "CELL:0002", "mitochondrial function",
        ("mitochondrial membrane potential", "oxygen consumption rate", "ocr",
         "atp production", "mitochondrial respiration", "seahorse",
         "mitochondrial dysfunction", "mitochondrial morphology"),
        measurement_axis="metabolic_rate",
        go_xref=("GO:0007005",),
    ),
    PhenotypeTerm(
        "CELL:0003", "oxidative stress",
        ("ros levels", "reactive oxygen species", "oxidative damage",
         "lipid peroxidation", "glutathione levels"),
        measurement_axis="concentration",
        go_xref=("GO:0006979",),
    ),
    PhenotypeTerm(
        "CELL:0004", "protein aggregation",
        ("aggregation", "protein aggregates", "inclusion bodies",
         "insoluble fraction", "proteostasis", "aggregate formation"),
        measurement_axis="aggregate_count",
        go_xref=("GO:0043933",),
    ),
    PhenotypeTerm(
        "CELL:0005", "autophagy flux",
        ("autophagic flux", "lc3 puncta", "lc3-ii levels", "p62 levels",
         "autophagosome number", "lysosomal function"),
        measurement_axis="flux",
        go_xref=("GO:0006914",),
    ),
    PhenotypeTerm(
        "CELL:0006", "calcium signaling",
        ("calcium transients", "ca2+ transients", "calcium imaging",
         "calcium influx", "intracellular calcium", "spontaneous calcium activity",
         "calcium oscillations", "calcium handling"),
        measurement_axis="event_rate",
        go_xref=("GO:0019722",),
        typical_assays=("calcium_imaging",),
    ),
    PhenotypeTerm(
        "CELL:0007", "ER stress",
        ("endoplasmic reticulum stress", "unfolded protein response", "upr",
         "chop induction", "bip levels", "xbp1 splicing"),
        measurement_axis="marker_level",
        go_xref=("GO:0034976",),
    ),
    # ---- splicing ----
    PhenotypeTerm(
        "SPL:0001", "exon skipping",
        ("exon skipping", "skipped exon", "exon exclusion", "cassette exon",
         "exon inclusion", "alternative exon usage", "percent spliced in",
         "psi change", "delta psi"),
        measurement_axis="psi",
        go_xref=("GO:0000380",),
        typical_assays=("splicing_analysis",),
        parent="SPL:0000",
    ),
    PhenotypeTerm(
        "SPL:0002", "intron retention",
        ("intron retention", "retained intron", "unspliced transcript",
         "intronic reads"),
        measurement_axis="psi",
        go_xref=("GO:0000380",),
        typical_assays=("splicing_analysis",),
        parent="SPL:0000",
    ),
    PhenotypeTerm(
        "SPL:0003", "cryptic exon inclusion",
        ("cryptic exon", "cryptic splicing", "de novo exon",
         "nonconserved cryptic exon", "cryptic exon inclusion"),
        measurement_axis="psi",
        typical_assays=("splicing_analysis",),
        parent="SPL:0000",
        notes="Hallmark of TDP-43 loss of nuclear function.",
    ),
    PhenotypeTerm(
        "SPL:0004", "alternative polyadenylation",
        ("alternative polyadenylation", "apa", "3'utr shortening",
         "polyadenylation site usage", "utr lengthening"),
        measurement_axis="pau",
        typical_assays=("splicing_analysis",),
        parent="SPL:0000",
    ),
    PhenotypeTerm(
        "SPL:0005", "nonsense-mediated decay",
        ("nmd", "nonsense-mediated decay", "nmd target", "ptc-containing transcript",
         "transcript degradation", "nonsense mediated mrna decay"),
        measurement_axis="transcript_level",
        go_xref=("GO:0000184",),
        parent="SPL:0000",
    ),
    PhenotypeTerm(
        "SPL:0006", "isoform switching",
        ("isoform switch", "isoform ratio", "transcript isoform usage",
         "differential isoform expression", "isoform abundance"),
        measurement_axis="isoform_fraction",
        typical_assays=("splicing_analysis",),
        parent="SPL:0000",
    ),
    PhenotypeTerm(
        "SPL:0007", "splicing factor mislocalization",
        ("nuclear depletion", "cytoplasmic mislocalization", "rbp mislocalization",
         "tdp-43 mislocalization", "nuclear clearance", "nucleocytoplasmic transport"),
        measurement_axis="localization_ratio",
        go_xref=("GO:0006406",),
        parent="SPL:0000",
    ),
    PhenotypeTerm(
        "SPL:0000", "splicing phenotype",
        ("splicing defect", "aberrant splicing", "mis-splicing", "missplicing",
         "dysregulated splicing", "splicing alteration", "splicing dysregulation"),
        measurement_axis="composite",
        go_xref=("GO:0008380",),
    ),
    # ---- transcriptomic / proteomic ----
    PhenotypeTerm(
        "EXPR:0001", "gene expression change",
        ("differential expression", "mrna levels", "transcript abundance",
         "gene expression", "upregulation", "downregulation", "expression level"),
        measurement_axis="transcript_level",
        typical_assays=("rnaseq_bulk", "qpcr"),
    ),
    PhenotypeTerm(
        "EXPR:0002", "protein level change",
        ("protein levels", "protein abundance", "protein expression",
         "western blot quantification", "immunoreactivity"),
        measurement_axis="protein_level",
        typical_assays=("western_blot", "proteomics"),
    ),
    # ---- maturation / identity ----
    PhenotypeTerm(
        "DIFF:0001", "neuronal maturation",
        ("maturation", "maturity", "developmental stage", "maturation markers",
         "immature phenotype", "neuronal maturity"),
        measurement_axis="composite",
        go_xref=("GO:0048666",),
    ),
    PhenotypeTerm(
        "DIFF:0002", "differentiation efficiency",
        ("differentiation efficiency", "conversion efficiency", "purity",
         "map2-positive cells", "tuj1-positive", "percent neurons",
         "neuronal yield", "induction efficiency"),
        measurement_axis="fraction_positive",
        typical_assays=("imaging_morphology",),
    ),
    PhenotypeTerm(
        "DIFF:0003", "cell identity marker expression",
        ("marker expression", "subtype identity", "cortical layer identity",
         "excitatory identity", "vglut expression", "ctip2 expression"),
        measurement_axis="marker_level",
    ),
]

BY_ID: dict[str, PhenotypeTerm] = {t.term_id: t for t in TERMS}

# Synonym -> term_id. Longest-first matching happens in the matcher.
SYNONYM_INDEX: dict[str, str] = {}
for _t in TERMS:
    SYNONYM_INDEX[_t.label.lower()] = _t.term_id
    for _s in _t.synonyms:
        SYNONYM_INDEX[_s.lower()] = _t.term_id


def ancestors(term_id: str) -> list[str]:
    """Walk the parent chain. Used to roll up sparse evidence."""
    out, cur, seen = [], BY_ID.get(term_id), set()
    while cur and cur.parent and cur.parent not in seen:
        seen.add(cur.parent)
        out.append(cur.parent)
        cur = BY_ID.get(cur.parent)
    return out


def commensurable(a: str, b: str) -> bool:
    """True if two phenotype terms are curated as the same quantity.

    Decided by an explicit `pool_group`, never by measurement axis. Sharing an
    axis only means two things are counted the same way: mEPSC and mIPSC
    frequency are both event rates, synaptic puncta and dendritic spines are
    both counts per length, and pooling either pair would merge distinct
    biology and corrupt the direction-consistency statistics.
    """
    ta, tb = BY_ID.get(a), BY_ID.get(b)
    if not ta or not tb:
        return False
    return bool(ta.pool_group) and ta.pool_group == tb.pool_group
