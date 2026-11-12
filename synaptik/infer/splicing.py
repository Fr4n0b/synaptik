"""
Splicing-specific analysis layer.

Splicing gets its own module because the generic phenotype->pathway machinery
loses the information that matters most in this domain. A splicing phenotype
is not just "aberrant splicing" -- it is a specific event, in a specific gene,
at specific coordinates, with a direction (dPSI) and, critically, a REGULATOR:
the RNA-binding protein whose loss or gain produced it.

That regulator is the actionable node. If your mutation causes mis-splicing,
the therapeutic question is rarely "which drug hits this exon" -- it is:

  (a) which RBP's activity is perturbed, and can it be restored?
  (b) can an ASO block or promote this specific event directly?
  (c) is this event a known NMD trigger, meaning the real defect is loss of
      the protein rather than production of a bad isoform?

This module answers those. It also implements RBP motif proximity scoring,
which is what lets you go from "exon X is skipped" to "and the nearby motifs
are TDP-43 binding sites" -- a hypothesis you can test.

Coordinate handling note: dPSI sign conventions differ between tools (rMATS
reports IncLevelDifference as control-minus-treatment in some versions,
treatment-minus-control in others). We normalize to a single convention --
positive dPSI means MORE inclusion in the perturbed condition -- and record
which tool produced each value so the normalization is auditable.
"""

from __future__ import annotations

import logging
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field

from ..schema import (
    Direction, ExperimentRecord, SplicingDetail, SplicingEventType,
)

log = logging.getLogger(__name__)


# --------------------------------------------------------------------------
# RNA-binding protein reference
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class RBPProfile:
    """A splicing regulator, its binding preference, and its disease context."""

    symbol: str
    motifs: tuple[str, ...]  # IUPAC consensus, 5'->3' on the sense strand
    binding_position: str  # where binding typically acts
    canonical_effect: str  # what binding usually does to inclusion
    disease_links: tuple[str, ...] = ()
    notes: str = ""


# Curated from CLIP-seq and functional studies. Motif consensus sequences use
# IUPAC codes. This set covers the RBPs most relevant to neurodegeneration and
# neurodevelopmental disease; it is extensible by adding entries here.
RBP_REFERENCE: dict[str, RBPProfile] = {
    "TARDBP": RBPProfile(
        "TARDBP", ("UGUGUGU", "GUGUGUG", "UGUGU"),
        "deep intronic, downstream of the regulated exon",
        "represses cryptic exon inclusion; loss causes cryptic exon inclusion",
        ("ALS", "FTD", "LATE"),
        "TDP-43. Nuclear depletion is the central event in ~97% of ALS. "
        "Loss-of-function signature is inclusion of nonconserved cryptic exons "
        "(STMN2, UNC13A are the best-characterized targets).",
    ),
    "FUS": RBPProfile(
        "FUS", ("GUGGU", "GGUG", "AUUAAA"),
        "broad; promoter-proximal and intronic",
        "context-dependent; both represses and promotes inclusion",
        ("ALS", "FTD"),
        "Mutations cluster in the NLS and cause cytoplasmic mislocalization.",
    ),
    "HNRNPA1": RBPProfile(
        "HNRNPA1", ("UAGGGA", "UAGGGU", "UAG"),
        "exonic and intronic splicing silencers",
        "represses inclusion; antagonizes SR proteins",
        ("ALS", "multisystem proteinopathy"),
    ),
    "PTBP1": RBPProfile(
        "PTBP1", ("CUCUCU", "UCUUC", "CUCU"),
        "polypyrimidine tract, upstream of the 3' splice site",
        "represses inclusion by competing with U2AF",
        ("neuronal differentiation",),
        "PTBP1->PTBP2 switch drives neuronal maturation; a major confound in "
        "iPSC-neuron splicing comparisons across DIV.",
    ),
    "PTBP2": RBPProfile(
        "PTBP2", ("CUCUCU", "UCUU"),
        "polypyrimidine tract",
        "represses inclusion; neuronal paralog of PTBP1",
        ("neuronal maturation",),
    ),
    "NOVA1": RBPProfile(
        "NOVA1", ("YCAY", "UCAUUCAU"),
        "intronic, position determines effect",
        "promotes inclusion when downstream, represses when upstream",
        ("POMA", "epilepsy"),
    ),
    "NOVA2": RBPProfile(
        "NOVA2", ("YCAY", "UCAU"),
        "intronic clusters",
        "position-dependent; regulates synaptic gene splicing",
        ("neurodevelopmental disorder",),
    ),
    "RBFOX1": RBPProfile(
        "RBFOX1", ("UGCAUG", "GCAUG"),
        "intronic; downstream promotes, upstream represses",
        "promotes inclusion when bound downstream of the exon",
        ("autism", "epilepsy", "schizophrenia"),
        "A2BP1. Strong autism association; regulates large neuronal splicing "
        "programs including synaptic genes.",
    ),
    "RBFOX2": RBPProfile(
        "RBFOX2", ("UGCAUG", "GCAUG"),
        "intronic",
        "position-dependent, same logic as RBFOX1",
        ("neurodevelopment",),
    ),
    "SRSF1": RBPProfile(
        "SRSF1", ("GGAGGA", "AGGAC", "GAAGAA"),
        "exonic splicing enhancers",
        "promotes inclusion; recruits U1/U2AF",
        ("cancer", "SMA modifier"),
    ),
    "SRSF2": RBPProfile(
        "SRSF2", ("GGNG", "SSNG"),
        "exonic splicing enhancers",
        "promotes inclusion",
        ("MDS", "leukemia"),
    ),
    "HNRNPH1": RBPProfile(
        "HNRNPH1", ("GGGGG", "GGGA"),
        "G-runs, intronic and exonic",
        "context-dependent",
        ("neurodevelopmental disorder",),
    ),
    "MBNL1": RBPProfile(
        "MBNL1", ("YGCY", "CGCU", "UGCU"),
        "intronic",
        "promotes adult splicing patterns; sequestered in DM1",
        ("myotonic dystrophy",),
        "Sequestration by CUG repeat RNA causes the DM1 splicing signature.",
    ),
    "CELF1": RBPProfile(
        "CELF1", ("UGUU", "GUGUGU"),
        "intronic",
        "antagonizes MBNL; promotes fetal splicing patterns",
        ("myotonic dystrophy",),
    ),
    "QKI": RBPProfile(
        "QKI", ("ACUAAY", "UAAY"),
        "intronic",
        "regulates myelination-associated splicing",
        ("schizophrenia", "glioma"),
    ),
    "ELAVL3": RBPProfile(
        "ELAVL3", ("UUUUU", "AUUUA"),
        "3'UTR and intronic",
        "neuronal; stabilizes transcripts and regulates splicing",
        ("neurodegeneration",),
    ),
    "SFPQ": RBPProfile(
        "SFPQ", ("UGUGAA", "GUAGUA"),
        "long intron, co-transcriptional",
        "required for long-intron transcription in neurons",
        ("ALS", "FTD"),
    ),
    "MATR3": RBPProfile(
        "MATR3", ("AUCUU", "UCUU"),
        "intronic",
        "represses inclusion; interacts with TDP-43",
        ("ALS", "distal myopathy"),
    ),
    "HNRNPA2B1": RBPProfile(
        "HNRNPA2B1", ("UAGG", "AGGUAA"),
        "splicing silencers",
        "represses inclusion",
        ("multisystem proteinopathy", "ALS"),
    ),
}

# IUPAC ambiguity expansion for motif search.
IUPAC = {
    "A": "A", "C": "C", "G": "G", "U": "U", "T": "T",
    "R": "[AG]", "Y": "[CUT]", "S": "[GC]", "W": "[AUT]",
    "K": "[GUT]", "M": "[AC]", "B": "[CGUT]", "D": "[AGUT]",
    "H": "[ACUT]", "V": "[ACG]", "N": "[ACGUT]",
}


def motif_to_regex(motif: str) -> re.Pattern:
    """Compile an IUPAC consensus into a regex matching both RNA and DNA."""
    parts = [IUPAC.get(c.upper(), re.escape(c)) for c in motif]
    return re.compile("".join(parts), re.IGNORECASE)


def scan_motifs(sequence: str, rbp: str) -> list[tuple[int, str]]:
    """Find all occurrences of an RBP's motifs. Returns [(position, motif)].

    Sequence may be DNA or RNA; T/U are treated interchangeably.
    """
    profile = RBP_REFERENCE.get(rbp.upper())
    if not profile or not sequence:
        return []
    seq = sequence.upper().replace("T", "U")
    hits: list[tuple[int, str]] = []
    for motif in profile.motifs:
        pat = motif_to_regex(motif)
        for m in pat.finditer(seq):
            hits.append((m.start(), motif))
    return sorted(set(hits))


def motif_density(sequence: str, rbp: str, window: int = 100) -> float:
    """Motif hits per 100 nt. A crude but useful enrichment signal."""
    if not sequence:
        return 0.0
    hits = scan_motifs(sequence, rbp)
    return len(hits) * window / max(len(sequence), 1)


# --------------------------------------------------------------------------
# Event normalization
# --------------------------------------------------------------------------

# rMATS versions have differed in the sign convention of IncLevelDifference.
# We standardize on: positive dPSI = MORE inclusion in the perturbed/mutant
# condition. Anything we cannot confidently orient is flagged rather than
# silently flipped, because a sign error here inverts the biology.
TOOL_SIGN_CONVENTION = {
    "rmats": "sample1_minus_sample2",
    "majiq": "condition_minus_control",
    "leafcutter": "unspecified",
    "suppa": "condition_minus_control",
    "vast-tools": "condition_minus_control",
    "rt-pcr": "unspecified",
}


@dataclass
class NormalizedSplicingEvent:
    gene: str
    event_type: SplicingEventType
    delta_psi: float | None
    direction: Direction
    coordinates: str | None = None
    genome_build: str | None = None
    detection_method: str | None = None
    regulator: str | None = None
    source_id: str = ""
    record_id: str = ""
    sign_confidence: str = "high"  # high | uncertain
    notes: list[str] = field(default_factory=list)

    def is_cryptic(self) -> bool:
        return self.event_type == SplicingEventType.SKIPPED_EXON and \
            self.direction == Direction.INCREASE

    def magnitude_bin(self) -> str:
        """|dPSI| bins. Below 0.1 is usually noise in short-read data."""
        if self.delta_psi is None:
            return "unknown"
        a = abs(self.delta_psi)
        if a < 0.10:
            return "negligible"
        if a < 0.20:
            return "modest"
        if a < 0.40:
            return "substantial"
        return "large"


def normalize_event(
    detail: SplicingDetail, record: ExperimentRecord
) -> NormalizedSplicingEvent:
    """Convert a raw SplicingDetail into a sign-normalized event."""
    method = (detail.detection_method or "").strip().lower()
    conv = TOOL_SIGN_CONVENTION.get(method, "unspecified")

    dpsi = detail.delta_psi
    sign_conf = "high"
    notes: list[str] = []

    if dpsi is not None and conv == "unspecified":
        sign_conf = "uncertain"
        notes.append(
            f"Sign convention for '{detail.detection_method or 'unknown method'}' "
            "is not standardized; verify the direction against the source figure."
        )

    # Prefer the stated outcome direction when it conflicts with the dPSI sign,
    # because prose is written by the author and the number may be transcribed
    # with either convention.
    stated = record.outcome.direction
    if dpsi is not None and stated in (Direction.INCREASE, Direction.DECREASE):
        implied = Direction.INCREASE if dpsi > 0 else Direction.DECREASE
        if implied != stated:
            sign_conf = "uncertain"
            notes.append(
                f"Reported dPSI ({dpsi:+.2f}) implies {implied.value} but the text "
                f"states {stated.value}. Using the stated direction."
            )
            dpsi = -dpsi

    direction = stated
    if direction == Direction.UNSPECIFIED and dpsi is not None:
        direction = Direction.INCREASE if dpsi > 0 else Direction.DECREASE

    return NormalizedSplicingEvent(
        gene=detail.gene.raw if detail.gene else "unknown",
        event_type=detail.event_type,
        delta_psi=dpsi,
        direction=direction,
        coordinates=detail.coordinates,
        genome_build=detail.genome_build,
        detection_method=detail.detection_method,
        regulator=detail.regulator.raw if detail.regulator else None,
        source_id=record.source.source_id,
        record_id=record.record_id,
        sign_confidence=sign_conf,
        notes=notes,
    )


# --------------------------------------------------------------------------
# Regulator inference
# --------------------------------------------------------------------------


@dataclass
class RegulatorHypothesis:
    """A candidate RBP explaining an observed splicing signature."""

    rbp: str
    n_supporting_events: int
    n_papers: int
    n_labs: int
    consistency: float  # fraction of events in the expected direction
    signature_match: str  # description of why this RBP fits
    profile: RBPProfile | None = None
    example_genes: list[str] = field(default_factory=list)
    therapeutic_notes: list[str] = field(default_factory=list)
    score: float = 0.0

    def caveats(self) -> list[str]:
        out = []
        if self.n_labs < 2:
            out.append("Supported by a single research group.")
        if self.n_supporting_events < 3:
            out.append(
                f"Only {self.n_supporting_events} events support this regulator."
            )
        if self.consistency < 0.7:
            out.append(
                "Events are not directionally consistent with a simple "
                "loss- or gain-of-function model."
            )
        return out


class SplicingAnalyzer:
    """Infers regulators and therapeutic entry points from splicing evidence."""

    def __init__(self, min_events: int = 2):
        self.min_events = min_events

    def collect_events(
        self, records: list[ExperimentRecord]
    ) -> list[NormalizedSplicingEvent]:
        out = []
        for r in records:
            if r.splicing:
                out.append(normalize_event(r.splicing, r))
        return out

    def infer_regulators(
        self, events: list[NormalizedSplicingEvent],
        records: list[ExperimentRecord] | None = None,
    ) -> list[RegulatorHypothesis]:
        """Rank RBPs that could explain the observed event set.

        Two sources of signal:
          - explicit: the paper named a regulator
          - implicit: the event signature matches a known RBP's loss/gain
            phenotype (e.g. cryptic exon inclusion -> TDP-43 loss)
        """
        lab_by_source: dict[str, str] = {}
        if records:
            for r in records:
                if r.source.lab_id:
                    lab_by_source[r.source.source_id] = r.source.lab_id

        buckets: dict[str, list[NormalizedSplicingEvent]] = defaultdict(list)
        for e in events:
            if e.regulator:
                key = e.regulator.upper().replace("TDP-43", "TARDBP")
                key = re.sub(r"[^A-Z0-9]", "", key)
                if key in RBP_REFERENCE:
                    buckets[key].append(e)

        # Implicit signature: cryptic exon inclusion is the TDP-43 loss hallmark
        cryptic = [e for e in events if e.is_cryptic()]
        if len(cryptic) >= 2 and "TARDBP" not in buckets:
            buckets["TARDBP"].extend(cryptic)

        hypotheses: list[RegulatorHypothesis] = []
        for rbp, evs in buckets.items():
            if len(evs) < self.min_events:
                continue
            profile = RBP_REFERENCE.get(rbp)
            papers = {e.source_id for e in evs}
            labs = {lab_by_source.get(s) for s in papers if lab_by_source.get(s)}

            dirs = Counter(e.direction.value for e in evs
                           if e.direction != Direction.UNSPECIFIED)
            consistency = (
                max(dirs.values()) / sum(dirs.values()) if dirs else 0.0
            )

            sig = self._signature_description(rbp, evs, profile)
            hyp = RegulatorHypothesis(
                rbp=rbp,
                n_supporting_events=len(evs),
                n_papers=len(papers),
                n_labs=len(labs),
                consistency=consistency,
                signature_match=sig,
                profile=profile,
                example_genes=sorted({e.gene for e in evs})[:8],
                therapeutic_notes=self._therapeutic_notes(rbp, evs),
            )
            hyp.score = (
                len(evs) * (1.0 + 0.5 * max(0, len(labs) - 1)) * (0.5 + consistency)
            )
            hypotheses.append(hyp)

        hypotheses.sort(key=lambda h: -h.score)
        return hypotheses

    @staticmethod
    def _signature_description(
        rbp: str, evs: list[NormalizedSplicingEvent], profile: RBPProfile | None
    ) -> str:
        cryptic_n = sum(1 for e in evs if e.is_cryptic())
        ri_n = sum(1 for e in evs
                   if e.event_type == SplicingEventType.RETAINED_INTRON)
        inc = sum(1 for e in evs if e.direction == Direction.INCREASE)
        dec = sum(1 for e in evs if e.direction == Direction.DECREASE)

        if rbp == "TARDBP" and cryptic_n >= 2:
            return (
                f"{cryptic_n} cryptic exon inclusion events -- the canonical "
                "signature of TDP-43 nuclear loss of function."
            )
        if ri_n >= 2:
            return f"{ri_n} intron retention events, consistent with reduced " \
                   f"{rbp} splicing activity."
        if inc > dec * 2 and profile:
            return (
                f"Predominantly increased inclusion ({inc} up vs {dec} down); "
                f"{rbp} canonically {profile.canonical_effect}, so this fits "
                "loss of repression."
            )
        if dec > inc * 2 and profile:
            return (
                f"Predominantly decreased inclusion ({dec} down vs {inc} up); "
                f"consistent with loss of {rbp}-dependent enhancement."
            )
        return f"{len(evs)} events associated with {rbp} in the literature."

    @staticmethod
    def _therapeutic_notes(
        rbp: str, evs: list[NormalizedSplicingEvent]
    ) -> list[str]:
        notes: list[str] = []
        if rbp == "TARDBP":
            notes.append(
                "Restoring nuclear TDP-43 function is not directly druggable, "
                "but individual cryptic exons ARE targetable with ASOs -- the "
                "STMN2 cryptic exon is the proof of concept."
            )
            genes = {e.gene.upper() for e in evs}
            if "STMN2" in genes:
                notes.append(
                    "STMN2 cryptic exon present: ASO-mediated exclusion restores "
                    "full-length STMN2 and axonal regeneration in iPSC neurons."
                )
            if "UNC13A" in genes:
                notes.append(
                    "UNC13A cryptic exon present: a major ALS GWAS locus; the "
                    "risk haplotype increases cryptic inclusion."
                )
        if rbp in ("PTBP1", "PTBP2"):
            notes.append(
                "PTBP1/2 levels change dramatically during neuronal maturation. "
                "Before interpreting this as a disease phenotype, confirm your "
                "mutant and control neurons are matched for DIV and maturation "
                "markers -- this is the most common confound in the field."
            )
        if rbp in ("SRSF1", "SRSF2"):
            notes.append(
                "SR protein activity is regulated by SRPK1/CLK kinases, which "
                "are druggable (e.g. CLK inhibitors such as TG003)."
            )
        if rbp == "MBNL1":
            notes.append(
                "MBNL sequestration is reversible by small molecules or ASOs "
                "targeting the repeat RNA."
            )
        notes.append(
            "Any splicing event with a defined exon boundary is in principle "
            "ASO-addressable; feasibility depends on the event type and whether "
            "exclusion or inclusion is desired."
        )
        return notes


# --------------------------------------------------------------------------
# ASO design feasibility
# --------------------------------------------------------------------------


@dataclass
class ASOStrategy:
    event_type: SplicingEventType
    goal: str  # "promote_inclusion" | "promote_exclusion"
    feasible: bool
    approach: str
    target_region: str
    precedent: str = ""
    caveats: list[str] = field(default_factory=list)


def aso_strategy(
    event_type: SplicingEventType, current_direction: Direction
) -> ASOStrategy:
    """What ASO approach would correct this event?

    Splice-switching ASOs work by steric blocking. Blocking a splice site or a
    splicing enhancer promotes exclusion; blocking a silencer promotes
    inclusion. Which one you want depends on the direction of the defect.
    """
    goal = (
        "promote_exclusion" if current_direction == Direction.INCREASE
        else "promote_inclusion"
    )

    if event_type == SplicingEventType.SKIPPED_EXON:
        if goal == "promote_exclusion":
            return ASOStrategy(
                event_type, goal, True,
                "Steric-blocking ASO targeting the 3' or 5' splice site of the "
                "aberrantly included exon, or an exonic splicing enhancer within it.",
                "Splice acceptor/donor of the cryptic exon, or ESE within it.",
                precedent="STMN2 cryptic exon exclusion in ALS iPSC neurons; "
                          "clinically validated chemistry (2'-MOE, PMO).",
                caveats=[
                    "Requires the exon boundary to be mapped precisely.",
                    "Screen a tiling series -- efficacy is highly position-dependent.",
                ],
            )
        return ASOStrategy(
            event_type, goal, True,
            "Steric-blocking ASO targeting an intronic splicing silencer (ISS) "
            "to de-repress inclusion.",
            "ISS elements in the flanking introns, typically within 100 nt of "
            "the splice site.",
            precedent="Nusinersen (SMN2 exon 7 inclusion) -- the archetypal "
                      "approved splice-switching ASO.",
            caveats=[
                "Silencer position must be identified, usually by tiling.",
                "Inclusion is generally harder to achieve than exclusion.",
            ],
        )

    if event_type == SplicingEventType.RETAINED_INTRON:
        return ASOStrategy(
            event_type, goal, goal == "promote_exclusion",
            "ASOs against retained introns are less reliable; intron retention "
            "often reflects a co-transcriptional or export defect rather than a "
            "local splice-site problem.",
            "5' splice site of the retained intron, if a weak site is identified.",
            caveats=[
                "Consider whether the retention is a consequence rather than a "
                "cause -- correcting it may not restore function.",
                "Retained-intron transcripts are frequently NMD substrates; "
                "check whether the real defect is transcript loss.",
            ],
        )

    if event_type == SplicingEventType.ALT_POLYA:
        return ASOStrategy(
            event_type, goal, True,
            "ASOs blocking the proximal polyadenylation signal can shift usage "
            "to the distal site.",
            "AAUAAA hexamer and downstream GU-rich element of the proximal site.",
            caveats=["Effects on transcript stability must be checked separately."],
        )

    if event_type in (SplicingEventType.ALT_5_SS, SplicingEventType.ALT_3_SS):
        return ASOStrategy(
            event_type, goal, True,
            "Block the undesired splice site to redirect the spliceosome to the "
            "alternative site.",
            "The aberrantly used 5' or 3' splice site.",
            caveats=["Competing sites must be close enough that redirection is "
                     "the favoured outcome."],
        )

    return ASOStrategy(
        event_type, goal, False,
        "No standard ASO approach for this event class.",
        "n/a",
        caveats=["Consider small-molecule splicing modulators or regulator-level "
                 "intervention instead."],
    )


# --------------------------------------------------------------------------
# NMD assessment
# --------------------------------------------------------------------------


def nmd_likelihood(
    event_type: SplicingEventType, direction: Direction,
    frame_disrupting: bool | None = None,
) -> tuple[str, str]:
    """Is this event likely to trigger nonsense-mediated decay?

    This matters because it changes the therapeutic goal entirely. If the
    aberrant transcript is degraded, the defect is loss of protein and the aim
    is to restore normal splicing. If it is translated, you may also be dealing
    with a toxic or dominant-negative product.
    """
    if frame_disrupting is False:
        return ("low", "Event preserves reading frame; transcript is likely "
                       "translated rather than degraded.")

    high_risk = {
        SplicingEventType.RETAINED_INTRON,
        SplicingEventType.SKIPPED_EXON,
    }
    if event_type in high_risk and direction == Direction.INCREASE:
        return (
            "high",
            "Inclusion of a cryptic exon or retained intron frequently "
            "introduces a premature termination codon, making the transcript an "
            "NMD substrate. The functional consequence is usually LOSS of the "
            "encoded protein -- measure protein level, not just transcript.",
        )
    if event_type == SplicingEventType.SKIPPED_EXON and \
            direction == Direction.DECREASE:
        return (
            "moderate",
            "Exon skipping disrupts frame in roughly two thirds of cases. "
            "Check whether the skipped exon length is a multiple of three.",
        )
    return ("unknown", "Frame impact not determinable from the event type alone; "
                       "check the exon length and reading frame.")
