"""
Human Single-Cell Post-Mortem Translational Atlas Engine.

Bridges in vitro human iPSC neuron phenotypes and splicing tables to human
post-mortem single-nucleus RNA-seq (snRNA-seq) disease tissue atlases:
  * NYGC ALS Consortium / Humphrey et al. Nat Neurosci 2023 (ALS/FTD motor cortex & spinal cord)
  * Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD) / Mathys et al. Nature 2019
  * Allen Brain Cell Atlas human primary motor cortex (M1) reference cell types.

Evaluates:
  1. Human Translational Concordance: does this in vitro finding replicate in
     real human disease brain tissue?
  2. Cell-type specificity: which human brain cell type in vivo normally expresses
     and exhibits this dysregulation (Upper Motor Neurons / Betz cells, Lower Motor
     Neurons, Cortical Interneurons, Astrocytes, Microglia)?
  3. In vitro developmental maturation gap: is the gene an adult-maturation marker
     poorly expressed at early iPSC differentiation stages (DIV 21-35)?
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Any

from ..schema import Direction


class BrainCellType(str, Enum):
    """Major human brain cell types."""

    UPPER_MOTOR_NEURON = "Upper Motor Neuron (Betz cell)"
    LOWER_MOTOR_NEURON = "Lower Motor Neuron (Spinal cord)"
    CORTICAL_EXCITATORY = "Cortical Projection Neuron (Layer 2-6)"
    INTERNEURON = "Cortical Interneuron (GABAergic)"
    ASTROCYTE = "Astrocyte"
    MICROGLIA = "Microglia"
    OLIGODENDROCYTE = "Oligodendrocyte"
    OPC = "Oligodendrocyte Precursor Cell"


class DiseaseAtlas(str, Enum):
    """Human post-mortem disease atlases."""

    ALS = "Amyotrophic Lateral Sclerosis (ALS)"
    FTD = "Frontotemporal Dementia (FTD / FTLD-TDP)"
    AD = "Alzheimer's Disease (AD)"


class PostMortemStatus(str, Enum):
    """Status of an in vitro event when compared to human post-mortem brain tissue."""

    CONCORDANT = "concordant"                  # Confirmed in human disease post-mortem brain tissue
    DISCORDANT = "discordant"                  # Contradicted in human disease post-mortem brain tissue
    DEVELOPMENTAL_GAP = "developmental_gap"    # Adult-specific marker (immature in early iPSC cultures)
    CELL_TYPE_SPECIFIC = "cell_type_specific"  # Specific to a non-neuronal or non-target cell lineage
    UNCHARACTERIZED = "uncharacterized"        # Not yet mapped in existing post-mortem snRNA-seq


# --------------------------------------------------------------------------
# Curated cell-type specific marker sets (Allen Brain Cell Atlas & PsychENCODE)
# --------------------------------------------------------------------------
CELL_TYPE_MARKERS: dict[BrainCellType, set[str]] = {
    BrainCellType.UPPER_MOTOR_NEURON: {
        "NEFH", "FEZF2", "CRYM", "POU3F1", "BHLHE22", "COL5A2", "VAT1L", "NR4A2"
    },
    BrainCellType.LOWER_MOTOR_NEURON: {
        "CHAT", "ISL1", "MNX1", "SLC18A3", "PRPH", "RET", "FOXP1", "ACHE"
    },
    BrainCellType.CORTICAL_EXCITATORY: {
        "SLC17A7", "TBR1", "RORB", "CUX2", "SATB2", "NEUROD6", "CAMK2A", "GRIN1"
    },
    BrainCellType.INTERNEURON: {
        "GAD1", "GAD2", "SLC32A1", "PVALB", "SST", "VIP", "LAMP5", "CCK"
    },
    BrainCellType.ASTROCYTE: {
        "AQP4", "GFAP", "ALDH1L1", "SLC1A2", "SLC1A3", "SOX9", "GJA1", "VIM"
    },
    BrainCellType.MICROGLIA: {
        "AIF1", "TMEM119", "CX3CR1", "P2RY12", "C1QA", "C1QB", "ITGAM", "TYROBP"
    },
    BrainCellType.OLIGODENDROCYTE: {
        "MBP", "MOG", "OLIG2", "PLP1", "MAG", "CNP", "CLDN11"
    },
    BrainCellType.OPC: {
        "PDGFRA", "CSPG4", "OLIG1", "SOX10"
    },
}

# Genes that require prolonged maturation in human neurons (often adult-only)
ADULT_MATURATION_MARKERS: set[str] = {
    "NEFH", "GRIN2A", "SLC12A5", "KCC2", "PVALB", "CALB1", "SYN2", "HCN1"
}


# --------------------------------------------------------------------------
# Validated post-mortem single-cell disease signatures
# --------------------------------------------------------------------------
@dataclass
class PostMortemEvent:
    gene: str
    disease: DiseaseAtlas
    primary_cell_types: list[BrainCellType]
    direction: Direction
    mechanism: str
    evidence_source: str
    is_cryptic_inclusion: bool = False


POST_MORTEM_CATALOG: list[PostMortemEvent] = [
    # ALS post-mortem motor cortex / spinal cord
    PostMortemEvent(
        gene="STMN2",
        disease=DiseaseAtlas.ALS,
        primary_cell_types=[BrainCellType.UPPER_MOTOR_NEURON, BrainCellType.LOWER_MOTOR_NEURON],
        direction=Direction.DECREASE,
        mechanism="TDP-43 loss-induced cryptic exon inclusion triggering premature polyadenylation and NMD",
        evidence_source="Klim et al. Nat Neurosci 2019; Prudencio et al. Acta Neuropathol 2020",
        is_cryptic_inclusion=True,
    ),
    PostMortemEvent(
        gene="UNC13A",
        disease=DiseaseAtlas.ALS,
        primary_cell_types=[BrainCellType.UPPER_MOTOR_NEURON, BrainCellType.CORTICAL_EXCITATORY],
        direction=Direction.DECREASE,
        mechanism="rs12608932/rs12973192 cryptic splice donor activation and exon inclusion upon TDP-43 nuclear clearance",
        evidence_source="Brown et al. Nature 2022; Ma et al. Nature 2022",
        is_cryptic_inclusion=True,
    ),
    PostMortemEvent(
        gene="TARDBP",
        disease=DiseaseAtlas.ALS,
        primary_cell_types=[BrainCellType.UPPER_MOTOR_NEURON, BrainCellType.LOWER_MOTOR_NEURON],
        direction=Direction.DECREASE,
        mechanism="Nuclear depletion with cytoplasmic hyperphosphorylated insolubility",
        evidence_source="Neumann et al. Science 2006; NYGC ALS Consortium 2023",
    ),
    PostMortemEvent(
        gene="AARS1",
        disease=DiseaseAtlas.ALS,
        primary_cell_types=[BrainCellType.UPPER_MOTOR_NEURON, BrainCellType.LOWER_MOTOR_NEURON],
        direction=Direction.DECREASE,
        mechanism="Cryptic splicing alteration under TDP-43 nuclear loss",
        evidence_source="Humphrey et al. Nat Neurosci 2023",
        is_cryptic_inclusion=True,
    ),
    PostMortemEvent(
        gene="KALRN",
        disease=DiseaseAtlas.ALS,
        primary_cell_types=[BrainCellType.UPPER_MOTOR_NEURON, BrainCellType.CORTICAL_EXCITATORY],
        direction=Direction.DECREASE,
        mechanism="Down-regulation and exon skipping causing dendritic spine retraction",
        evidence_source="Humphrey et al. Nat Neurosci 2023",
        is_cryptic_inclusion=False,
    ),
    PostMortemEvent(
        gene="ELAVL3",
        disease=DiseaseAtlas.ALS,
        primary_cell_types=[BrainCellType.UPPER_MOTOR_NEURON, BrainCellType.LOWER_MOTOR_NEURON],
        direction=Direction.DECREASE,
        mechanism="Neuronal HuC down-regulation in degenerating motor units",
        evidence_source="NYGC ALS Consortium 2023",
    ),
    PostMortemEvent(
        gene="ATG4B",
        disease=DiseaseAtlas.ALS,
        primary_cell_types=[BrainCellType.UPPER_MOTOR_NEURON, BrainCellType.LOWER_MOTOR_NEURON],
        direction=Direction.DECREASE,
        mechanism="Autophagic flux impairment in motor neurons",
        evidence_source="NYGC ALS Consortium 2023",
    ),
    PostMortemEvent(
        gene="CELF5",
        disease=DiseaseAtlas.ALS,
        primary_cell_types=[BrainCellType.CORTICAL_EXCITATORY],
        direction=Direction.DECREASE,
        mechanism="Compensatory splicing regulator alteration",
        evidence_source="Humphrey et al. Nat Neurosci 2023",
    ),
    PostMortemEvent(
        gene="GFAP",
        disease=DiseaseAtlas.ALS,
        primary_cell_types=[BrainCellType.ASTROCYTE],
        direction=Direction.INCREASE,
        mechanism="Reactive astrogliosis in motor cortex and ventral horn",
        evidence_source="NYGC ALS Consortium 2023",
    ),
    PostMortemEvent(
        gene="AIF1",
        disease=DiseaseAtlas.ALS,
        primary_cell_types=[BrainCellType.MICROGLIA],
        direction=Direction.INCREASE,
        mechanism="Microglial proliferation and neuroinflammatory activation",
        evidence_source="NYGC ALS Consortium 2023",
    ),

    # FTD post-mortem
    PostMortemEvent(
        gene="GRN",
        disease=DiseaseAtlas.FTD,
        primary_cell_types=[BrainCellType.MICROGLIA, BrainCellType.CORTICAL_EXCITATORY],
        direction=Direction.DECREASE,
        mechanism="Haploinsufficiency leading to lysosomal lipid storage dysfunction",
        evidence_source="Baker et al. Nature 2006; Cruts et al. Nature 2006",
    ),
    PostMortemEvent(
        gene="MAPT",
        disease=DiseaseAtlas.FTD,
        primary_cell_types=[BrainCellType.CORTICAL_EXCITATORY],
        direction=Direction.INCREASE,
        mechanism="Exon 10 inclusion shift altering 3R:4R tau isoform balance",
        evidence_source="Hutton et al. Nature 1998",
        is_cryptic_inclusion=False,
    ),
    PostMortemEvent(
        gene="TBK1",
        disease=DiseaseAtlas.FTD,
        primary_cell_types=[BrainCellType.MICROGLIA, BrainCellType.CORTICAL_EXCITATORY],
        direction=Direction.DECREASE,
        mechanism="Impaired mitophagy and optineurin phosphorylation",
        evidence_source="Freischmidt et al. Nat Neurosci 2015",
    ),

    # Alzheimer's Disease post-mortem
    PostMortemEvent(
        gene="APOE",
        disease=DiseaseAtlas.AD,
        primary_cell_types=[BrainCellType.ASTROCYTE, BrainCellType.MICROGLIA],
        direction=Direction.INCREASE,
        mechanism="Lipid homeostasis disruption and plaque-associated microgliosis",
        evidence_source="Mathys et al. Nature 2019; SEA-AD Consortium 2024",
    ),
    PostMortemEvent(
        gene="TREM2",
        disease=DiseaseAtlas.AD,
        primary_cell_types=[BrainCellType.MICROGLIA],
        direction=Direction.INCREASE,
        mechanism="Disease-Associated Microglia (DAM) state transition",
        evidence_source="Keren-Shaul et al. Cell 2017; Mathys et al. Nature 2019",
    ),
    PostMortemEvent(
        gene="BIN1",
        disease=DiseaseAtlas.AD,
        primary_cell_types=[BrainCellType.CORTICAL_EXCITATORY, BrainCellType.MICROGLIA],
        direction=Direction.DECREASE,
        mechanism="Endocytic vesicle trafficking defect and tau propagation",
        evidence_source="Chapuis et al. Mol Psychiatry 2013",
    ),
    PostMortemEvent(
        gene="SYP",
        disease=DiseaseAtlas.AD,
        primary_cell_types=[BrainCellType.CORTICAL_EXCITATORY],
        direction=Direction.DECREASE,
        mechanism="Early synaptic loss in entorhinal cortex and hippocampus CA1",
        evidence_source="Terry et al. Ann Neurol 1991; SEA-AD 2024",
    ),
]


@dataclass
class GeneTranslationCall:
    gene: str
    status: PostMortemStatus
    primary_cell_type: str
    in_vitro_direction: str
    post_mortem_direction: str
    mechanism: str
    evidence_source: str
    is_cryptic_exon: bool
    is_developmental_gap: bool
    notes: str


@dataclass
class TranslationalConcordanceResult:
    disease: DiseaseAtlas
    n_genes_evaluated: int
    n_concordant: int
    n_discordant: int
    n_developmental_gap: int
    n_cell_type_mismatch: int
    n_uncharacterized: int
    concordance_score_pct: float
    enriched_cell_types: list[tuple[str, int, float]]
    gene_calls: list[GeneTranslationCall]
    translational_verdict: str
    recommendations: list[str] = field(default_factory=list)

    def summary(self) -> str:
        return (
            f"Translational Concordance against {self.disease.value}:\n"
            f"  Score: {self.concordance_score_pct:.1f}% ({self.n_concordant} concordant, "
            f"{self.n_discordant} discordant, {self.n_developmental_gap} developmental gap, "
            f"{self.n_uncharacterized} uncharacterized)\n"
            f"  Verdict: {self.translational_verdict}"
        )


class HumanTranslationalConcordance:
    """Evaluates in vitro iPSC findings against human post-mortem brain tissue."""

    def __init__(self, disease: DiseaseAtlas = DiseaseAtlas.ALS) -> None:
        self.disease = disease
        self.pm_index = {
            e.gene.upper(): e
            for e in POST_MORTEM_CATALOG
            if e.disease == disease or e.gene.upper() in ("GFAP", "AIF1", "TARDBP")
        }

    def evaluate(
        self,
        gene_entries: list[tuple[str, Direction | None, float | None]],
    ) -> TranslationalConcordanceResult:
        """Evaluate a list of (gene, direction, effect_size/delta_psi)."""
        calls: list[GeneTranslationCall] = []
        n_conc = 0
        n_disc = 0
        n_dev = 0
        n_ct_mismatch = 0
        n_unchar = 0

        cell_type_counts: dict[BrainCellType, int] = {ct: 0 for ct in BrainCellType}

        for gene_raw, in_dir, effect in gene_entries:
            gene = gene_raw.strip().upper()
            if not gene:
                continue

            # Cell type mapping
            mapped_types = [
                ct for ct, markers in CELL_TYPE_MARKERS.items() if gene in markers
            ]
            for ct in mapped_types:
                cell_type_counts[ct] += 1

            pm_event = self.pm_index.get(gene)
            is_adult_marker = gene in ADULT_MATURATION_MARKERS

            if pm_event:
                # Check direction agreement.
                # For cryptic exon splicing, an increase in cryptic exon inclusion (delta_psi > 0)
                # is the hallmark defect confirmed in human post-mortem brain tissue.
                if pm_event.is_cryptic_inclusion and in_dir == Direction.INCREASE:
                    is_dir_match = True
                    is_dir_conflict = False
                elif pm_event.is_cryptic_inclusion and in_dir == Direction.DECREASE:
                    is_dir_match = False
                    is_dir_conflict = True
                else:
                    is_dir_match = (
                        in_dir is not None
                        and pm_event.direction is not None
                        and in_dir == pm_event.direction
                    )
                    is_dir_conflict = (
                        in_dir is not None
                        and pm_event.direction is not None
                        and in_dir != Direction.NO_CHANGE
                        and pm_event.direction != Direction.NO_CHANGE
                        and in_dir != pm_event.direction
                    )

                if is_dir_match or (in_dir is None and pm_event.is_cryptic_inclusion):
                    status = PostMortemStatus.CONCORDANT
                    n_conc += 1
                    notes = f"Corroborated in post-mortem human {self.disease.value} brain tissue ({pm_event.evidence_source})."
                elif is_dir_conflict:
                    status = PostMortemStatus.DISCORDANT
                    n_disc += 1
                    notes = f"Direction conflict: in vitro={in_dir.value}, post-mortem={pm_event.direction.value} ({pm_event.evidence_source})."
                else:
                    status = PostMortemStatus.CONCORDANT
                    n_conc += 1
                    notes = f"Associated with post-mortem {self.disease.value} pathology ({pm_event.evidence_source})."

                ct_str = ", ".join(ct.value for ct in pm_event.primary_cell_types)
                pm_dir_str = pm_event.direction.value if pm_event.direction else "unknown"
                mech = pm_event.mechanism
                src = pm_event.evidence_source
                has_cryptic = pm_event.is_cryptic_inclusion
            elif is_adult_marker:
                status = PostMortemStatus.DEVELOPMENTAL_GAP
                n_dev += 1
                ct_str = "Mature adult neuron"
                pm_dir_str = "adult-dependent"
                mech = "Adult-specific neurofilament or synaptic maturation marker; poorly expressed in immature DIV21 iPSC neurons"
                src = "PsychENCODE / Allen Brain Cell Atlas developmental series"
                has_cryptic = False
                notes = "Maturation Gap: Early iPSC differentiations lack adult expression levels of this gene."
            elif mapped_types and all(ct in (BrainCellType.MICROGLIA, BrainCellType.ASTROCYTE, BrainCellType.OLIGODENDROCYTE) for ct in mapped_types):
                status = PostMortemStatus.CELL_TYPE_SPECIFIC
                n_ct_mismatch += 1
                ct_str = ", ".join(ct.value for ct in mapped_types)
                pm_dir_str = "n/a"
                mech = "Expressed predominantly in glial lineages in human brain"
                src = "Allen Brain Cell Atlas"
                has_cryptic = False
                notes = "Lineage specificity: marker of non-neuronal glial cells in human brain."
            else:
                status = PostMortemStatus.UNCHARACTERIZED
                n_unchar += 1
                ct_str = mapped_types[0].value if mapped_types else "Broad / Pan-neuronal"
                pm_dir_str = "unreported"
                mech = "Not yet cataloged in post-mortem single-cell benchmarks"
                src = "n/a"
                has_cryptic = False
                notes = "Novel in vitro finding: check PubMed before claiming human pathology."

            calls.append(
                GeneTranslationCall(
                    gene=gene,
                    status=status,
                    primary_cell_type=ct_str,
                    in_vitro_direction=in_dir.value if in_dir else "reported",
                    post_mortem_direction=pm_dir_str,
                    mechanism=mech,
                    evidence_source=src,
                    is_cryptic_exon=has_cryptic,
                    is_developmental_gap=is_adult_marker,
                    notes=notes,
                )
            )

        total = len(calls)
        if total == 0:
            score = 0.0
        else:
            # Score formula: (concordant * 1.0 + uncharacterized * 0.4 - discordant * 1.0) / total
            effective = (n_conc * 1.0) + (n_unchar * 0.4) + (n_dev * 0.5) - (n_disc * 1.0)
            score = max(0.0, min(100.0, (effective / total) * 100.0))

        # Enriched cell types
        enriched: list[tuple[str, int, float]] = []
        for ct, cnt in cell_type_counts.items():
            if cnt > 0:
                pct = round((cnt / total) * 100.0, 1)
                enriched.append((ct.value, cnt, pct))
        enriched.sort(key=lambda x: x[1], reverse=True)

        # Verdict
        if n_conc >= 3 and n_disc == 0:
            verdict = (
                f"HIGH TRANSLATIONAL FIDELITY: Multiple key dysregulations ({n_conc}) "
                f"replicate hallmark post-mortem {self.disease.value} pathology with ZERO conflicts."
            )
        elif n_disc > 0:
            verdict = (
                f"MIXED TRANSLATIONAL FIDELITY: Detected {n_disc} directional contradiction(s) "
                f"against human post-mortem brain tissue. Inspect discordant genes closely."
            )
        elif n_conc > 0:
            verdict = (
                f"MODERATE TRANSLATIONAL FIDELITY: Corroborates {n_conc} post-mortem marker(s); "
                f"majority of events are novel in vitro candidates."
            )
        else:
            verdict = (
                "UNEXPLORED TRANSLATIONAL SPACE: None of the tested genes are currently "
                f"characterized in the post-mortem {self.disease.value} snRNA-seq benchmark."
            )

        # Actionable recommendations
        recs: list[str] = []
        if n_conc > 0:
            recs.append(
                f"Highlight the {n_conc} post-mortem concordant gene(s) in grant aims and manuscript figures "
                "to demonstrate in vivo disease relevance."
            )
        if n_dev > 0:
            recs.append(
                f"{n_dev} gene(s) represent adult neurodevelopmental maturation markers. Extend iPSC culture "
                "to DIV 45-60 or introduce optogenetic/co-culture stimulation if targeting these pathways."
            )
        if n_disc > 0:
            recs.append(
                "Address discordant findings directly: investigate whether in vitro culture conditions "
                "(e.g. ambient 21% O2 vs physiological brain 3% O2) inverted the transcriptional response."
            )

        return TranslationalConcordanceResult(
            disease=self.disease,
            n_genes_evaluated=total,
            n_concordant=n_conc,
            n_discordant=n_disc,
            n_developmental_gap=n_dev,
            n_cell_type_mismatch=n_ct_mismatch,
            n_uncharacterized=n_unchar,
            concordance_score_pct=round(score, 1),
            enriched_cell_types=enriched,
            gene_calls=calls,
            translational_verdict=verdict,
            recommendations=recs,
        )


def evaluate_translational_concordance(
    gene_entries: list[tuple[str, Direction | None, float | None]],
    disease: str = "als",
) -> TranslationalConcordanceResult:
    """Helper function to run concordance against ALS, FTD, or AD."""
    dis_enum = DiseaseAtlas.ALS
    d_clean = disease.lower()
    if "ftd" in d_clean or "ftld" in d_clean:
        dis_enum = DiseaseAtlas.FTD
    elif "ad" in d_clean or "alzheimer" in d_clean:
        dis_enum = DiseaseAtlas.AD

    evaluator = HumanTranslationalConcordance(disease=dis_enum)
    return evaluator.evaluate(gene_entries)
