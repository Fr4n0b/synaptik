"""
Canonical schema for SYNAPTIK.

The central object is the ExperimentRecord: one (paper, biological system,
perturbation, assay, measured outcome) tuple. A single paper typically yields
5-50 of these. Making this the unit of record -- rather than the paper -- is
what allows "has this experiment been done?" to be answered precisely.

Design notes:
  - Every field that can be ontology-normalized carries both the raw surface
    string (what the paper said) and the normalized ID. We never discard the
    raw text; normalization is lossy and reviewers need to audit it.
  - Effect direction and magnitude are separated. Many papers report direction
    without magnitude; we must not silently coerce that into a zero.
  - Provenance is mandatory. Every record knows which sentence(s) it came from.
"""

from __future__ import annotations

import hashlib
from datetime import date
from enum import Enum
from typing import Any

from pydantic import BaseModel, Field, field_validator


# --------------------------------------------------------------------------
# Controlled vocabularies
# --------------------------------------------------------------------------


class Direction(str, Enum):
    """Direction of a measured change relative to the stated control."""

    INCREASE = "increase"
    DECREASE = "decrease"
    NO_CHANGE = "no_change"  # explicitly tested, explicitly null -- precious
    BIPHASIC = "biphasic"
    UNSPECIFIED = "unspecified"


class EvidenceType(str, Enum):
    """How direct is the link between the perturbation and the outcome."""

    DIRECT_MEASUREMENT = "direct_measurement"
    RESCUE = "rescue"  # perturbation reverses a phenotype
    CORRELATION = "correlation"
    INFERRED = "inferred"  # author's interpretation, not measured
    REVIEW_ASSERTION = "review_assertion"  # secondary source


class CellModelClass(str, Enum):
    IPSC_NEURON_NGN2 = "ipsc_neuron_ngn2"
    IPSC_NEURON_SMALL_MOLECULE = "ipsc_neuron_small_molecule"
    IPSC_NPC = "ipsc_npc"
    IPSC_ORGANOID = "ipsc_organoid"
    PRIMARY_RODENT_NEURON = "primary_rodent_neuron"
    IMMORTALIZED_LINE = "immortalized_line"
    IN_VIVO_RODENT = "in_vivo_rodent"
    POST_MORTEM_HUMAN = "post_mortem_human"
    OTHER = "other"


class AssayClass(str, Enum):
    """Assay families. Orthogonality between families is what we score."""

    IMAGING_SYNAPTIC_PUNCTA = "imaging_synaptic_puncta"
    IMAGING_MORPHOLOGY = "imaging_morphology"
    ELECTROPHYS_PATCH = "electrophys_patch"
    ELECTROPHYS_MEA = "electrophys_mea"
    CALCIUM_IMAGING = "calcium_imaging"
    RNASEQ_BULK = "rnaseq_bulk"
    RNASEQ_SINGLE_CELL = "rnaseq_single_cell"
    SPLICING_ANALYSIS = "splicing_analysis"
    PROTEOMICS = "proteomics"
    WESTERN_BLOT = "western_blot"
    QPCR = "qpcr"
    VIABILITY = "viability"
    BIOCHEMICAL_BINDING = "biochemical_binding"
    OTHER = "other"


# Assay families that constitute genuinely independent lines of evidence.
# Two records from the same group are NOT orthogonal even if the family differs
# in name -- see infer/orthogonality.py.
ORTHOGONALITY_GROUPS: dict[AssayClass, str] = {
    AssayClass.IMAGING_SYNAPTIC_PUNCTA: "morphological",
    AssayClass.IMAGING_MORPHOLOGY: "morphological",
    AssayClass.ELECTROPHYS_PATCH: "functional_electrical",
    AssayClass.ELECTROPHYS_MEA: "functional_electrical",
    AssayClass.CALCIUM_IMAGING: "functional_electrical",
    AssayClass.RNASEQ_BULK: "transcriptomic",
    AssayClass.RNASEQ_SINGLE_CELL: "transcriptomic",
    AssayClass.SPLICING_ANALYSIS: "transcriptomic",
    AssayClass.PROTEOMICS: "proteomic",
    AssayClass.WESTERN_BLOT: "proteomic",
    AssayClass.QPCR: "transcriptomic",
    AssayClass.VIABILITY: "viability",
    AssayClass.BIOCHEMICAL_BINDING: "biochemical",
    AssayClass.OTHER: "other",
}


class PerturbationType(str, Enum):
    MUTATION_PATIENT = "mutation_patient"
    MUTATION_ENGINEERED = "mutation_engineered"
    KNOCKOUT = "knockout"
    KNOCKDOWN = "knockdown"
    OVEREXPRESSION = "overexpression"
    SMALL_MOLECULE = "small_molecule"
    ASO = "antisense_oligonucleotide"
    ISOFORM_SWITCH = "isoform_switch"
    NONE = "none"  # observational / baseline characterization
    OTHER = "other"


class SplicingEventType(str, Enum):
    """rMATS/MAJIQ-compatible event taxonomy."""

    SKIPPED_EXON = "SE"
    ALT_5_SS = "A5SS"
    ALT_3_SS = "A3SS"
    MUTUALLY_EXCLUSIVE = "MXE"
    RETAINED_INTRON = "RI"
    ALT_FIRST_EXON = "AFE"
    ALT_LAST_EXON = "ALE"
    ALT_POLYA = "APA"
    UNSPECIFIED = "unspecified"


# --------------------------------------------------------------------------
# Sub-records
# --------------------------------------------------------------------------


class Provenance(BaseModel):
    """Where in the source document this record came from. Mandatory."""

    section: str | None = None  # "methods", "results", "figure_3_legend"
    quote: str | None = Field(
        default=None,
        description="Verbatim supporting span, <=200 chars, for human audit only.",
    )
    figure_ref: str | None = None
    extractor_version: str = "unknown"
    extraction_confidence: float = Field(default=0.5, ge=0.0, le=1.0)

    @field_validator("quote")
    @classmethod
    def _truncate_quote(cls, v: str | None) -> str | None:
        if v and len(v) > 200:
            return v[:197] + "..."
        return v


class NormalizedTerm(BaseModel):
    """A surface string plus its ontology mapping. Both are retained."""

    raw: str
    ontology_id: str | None = None  # e.g. "GO:0007268", "HGNC:11138"
    ontology_label: str | None = None
    mapping_confidence: float = Field(default=0.0, ge=0.0, le=1.0)
    mapping_method: str | None = None  # "exact", "synonym", "embedding", "manual"

    def key(self) -> str:
        """Preferred identity for grouping: ontology ID if mapped, else raw."""
        return self.ontology_id or f"raw:{self.raw.strip().lower()}"


class CellModel(BaseModel):
    """The biological system. Comparability across records hinges on this."""

    model_class: CellModelClass = CellModelClass.OTHER
    raw_description: str = ""
    species: str | None = None
    differentiation_protocol: str | None = None  # "NGN2 doxycycline-inducible"
    protocol_citation: str | None = None
    days_in_vitro: int | None = None
    div_range: tuple[int, int] | None = None
    coculture: str | None = None  # "mouse astrocytes", "none"
    genetic_background: str | None = None  # cell line ID if stated
    is_isogenic_pair: bool | None = None

    def maturity_bin(self) -> str:
        """Coarse maturity bucket -- synaptic phenotypes are DIV-sensitive.

        Comparing DIV14 to DIV49 synapse counts is a category error, so this
        bin is used as a hard blocking key in experiment lookup.
        """
        d = self.days_in_vitro
        if d is None and self.div_range:
            d = sum(self.div_range) // 2
        if d is None:
            return "unknown"
        if d < 14:
            return "early"
        if d < 28:
            return "mid"
        if d < 49:
            return "late"
        return "extended"


class Perturbation(BaseModel):
    ptype: PerturbationType = PerturbationType.OTHER
    target: NormalizedTerm | None = None  # gene/protein/pathway acted upon
    agent: str | None = None  # compound name, guide RNA, ASO sequence ID
    variant_hgvs: str | None = None  # e.g. "NM_004006.2:c.1234A>G"
    variant_rsid: str | None = None
    dose: str | None = None
    duration: str | None = None
    control_description: str | None = None


class SplicingDetail(BaseModel):
    """Populated only when the record concerns a splicing outcome."""

    event_type: SplicingEventType = SplicingEventType.UNSPECIFIED
    gene: NormalizedTerm | None = None
    exon_id: str | None = None
    coordinates: str | None = None  # "chr1:1000-2000:+"
    genome_build: str | None = None
    delta_psi: float | None = Field(default=None, ge=-1.0, le=1.0)
    psi_control: float | None = Field(default=None, ge=0.0, le=1.0)
    psi_condition: float | None = Field(default=None, ge=0.0, le=1.0)
    detection_method: str | None = None  # "rMATS", "MAJIQ", "RT-PCR"
    regulator: NormalizedTerm | None = None  # RBP implicated, if any


class Outcome(BaseModel):
    """What was measured and what happened to it."""

    phenotype: NormalizedTerm
    direction: Direction = Direction.UNSPECIFIED
    effect_size: float | None = None
    effect_size_unit: str | None = None  # "percent_change", "fold", "cohens_d"
    p_value: float | None = None
    n_biological_replicates: int | None = None
    n_cells: int | None = None
    statistical_test: str | None = None
    markers_used: list[str] = Field(default_factory=list)  # SYN1, PSD95, MAP2
    is_quantified: bool = False

    @field_validator("p_value")
    @classmethod
    def _sane_p(cls, v: float | None) -> float | None:
        if v is not None and not (0.0 <= v <= 1.0):
            return None
        return v


class SourceMetadata(BaseModel):
    """Bibliographic + quality signals for the containing paper."""

    source_id: str  # "PMID:12345678", "DOI:10.1101/..."
    pmid: str | None = None
    pmcid: str | None = None
    doi: str | None = None
    title: str = ""
    journal: str | None = None
    journal_iso: str | None = None
    publication_date: date | None = None
    year: int | None = None
    is_preprint: bool = False
    is_review: bool = False
    is_retracted: bool = False
    retraction_notice: str | None = None
    authors: list[str] = Field(default_factory=list)
    corresponding_author: str | None = None
    affiliations: list[str] = Field(default_factory=list)
    lab_id: str | None = None  # assigned by normalize/labs.py
    citation_count: int | None = None
    open_access: bool = False
    full_text_available: bool = False

    # Journal quality signals -- deliberately plural, never a single number
    sjr_quartile: str | None = None  # "Q1".."Q4"
    sjr_score: float | None = None
    h_index: int | None = None

    def first_author_surname(self) -> str | None:
        if not self.authors:
            return None
        return self.authors[0].split()[-1] if self.authors[0] else None


class ExperimentRecord(BaseModel):
    """The atomic unit of the knowledge base."""

    record_id: str = ""
    source: SourceMetadata
    cell_model: CellModel
    perturbation: Perturbation
    outcome: Outcome
    assay_class: AssayClass = AssayClass.OTHER
    assay_raw: str = ""
    evidence_type: EvidenceType = EvidenceType.DIRECT_MEASUREMENT
    pathways: list[NormalizedTerm] = Field(default_factory=list)
    genes: list[NormalizedTerm] = Field(default_factory=list)
    splicing: SplicingDetail | None = None
    provenance: Provenance = Field(default_factory=Provenance)
    ingested_at: str | None = None

    def model_post_init(self, __context: Any) -> None:
        if not self.record_id:
            self.record_id = self.compute_id()

    def compute_id(self) -> str:
        """Deterministic ID so re-ingestion is idempotent."""
        basis = "|".join(
            [
                self.source.source_id,
                self.cell_model.model_class.value,
                self.cell_model.maturity_bin(),
                self.perturbation.ptype.value,
                (self.perturbation.target.key() if self.perturbation.target else ""),
                self.assay_class.value,
                self.outcome.phenotype.key(),
                self.outcome.direction.value,
            ]
        )
        return "rec_" + hashlib.sha1(basis.encode()).hexdigest()[:16]

    def orthogonality_group(self) -> str:
        return ORTHOGONALITY_GROUPS.get(self.assay_class, "other")

    def is_negative_result(self) -> bool:
        return self.outcome.direction == Direction.NO_CHANGE


class ExperimentQuery(BaseModel):
    """A description of an experiment the user is contemplating."""

    cell_model_class: CellModelClass | None = None
    differentiation_protocol: str | None = None
    days_in_vitro: int | None = None
    perturbation_type: PerturbationType | None = None
    perturbation_target: str | None = None
    variant_hgvs: str | None = None
    assay_class: AssayClass | None = None
    markers: list[str] = Field(default_factory=list)
    phenotype: str | None = None
    free_text: str | None = None

    def maturity_bin(self) -> str:
        return CellModel(days_in_vitro=self.days_in_vitro).maturity_bin()
