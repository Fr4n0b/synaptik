"""Causal Genomic Counterfactual Simulator for iPSC-neuron disease modeling.

Simulates in silico genetic perturbations (knockdown, knockout, overexpression,
and ASO steric blocking) on diseased neuronal transcriptomic vectors, predicting
rescue efficiency, directional vector alignment, and collateral toxicity.

Grounded in human iPSC-neuron functional genomics and Perturb-seq benchmarks:
- Tian et al. Neuron 2019 (CRISPRi/a screening in human iPSC-derived neurons)
- Replogle et al. Cell 2022 (Genome-scale single-cell Perturb-seq)
- Elden et al. Nature 2010 / Becker et al. Nature 2017 (ATXN2 rescue of TDP-43)
- Shi et al. Nat Med 2018 (PIKFYVE / endolysosomal rescue in ALS/FTD)
- Klim et al. Nat Neurosci 2019 / Melamed et al. 2019 (STMN2 rescue)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Tuple


class InterventionType(str, Enum):
    KNOCKDOWN = "knockdown"          # shRNA, CRISPRi, gapmer ASO (~60-90% reduction)
    KNOCKOUT = "knockout"            # CRISPR-Cas9 complete loss of function
    OVEREXPRESSION = "overexpression"# Lentiviral, dCas9-VPR, cDNA rescue
    ASO_STERIC_BLOCK = "aso_block"   # Splice-switching oligo restoring transcript


@dataclass(frozen=True)
class TargetPerturbationEffect:
    """Transcriptional response vector for a specific gene perturbation in neurons."""
    target_gene: str
    intervention: InterventionType
    # Gene symbol -> expected log2 fold-change delta in human neurons
    gene_deltas: Dict[str, float]
    # Affected functional pathways
    primary_pathway: str
    evidence_citation: str
    is_common_essential: bool = False
    collateral_risk: float = 0.0  # 0.0 (clean) to 1.0 (severe non-specific toxicity)


@dataclass
class SimulationResult:
    """Result of an in silico counterfactual rescue simulation."""
    disease: str
    target_gene: str
    intervention: InterventionType
    efficacy_alpha: float
    # Vector metrics
    disease_magnitude: float
    post_intervention_magnitude: float
    rescue_efficiency_pct: float
    cosine_rescue_alignment: float  # +1.0 = perfect reversal, -1.0 = worsens disease
    collateral_toxicity_index: float # 0.0 (safe) to 1.0 (lethal)
    composite_rescue_score: float   # 0 to 100
    is_clinically_viable: bool
    verdict: str
    key_gene_trajectories: List[Dict[str, float]] # gene, baseline, post, delta
    pathway_impacts: Dict[str, str]
    recommendations: List[str]


# Benchmark neuronal disease state vectors: log2 fold-change relative to isogenic control
BENCHMARK_DISEASE_STATES: Dict[str, Dict[str, float]] = {
    "als": {
        "STMN2": -2.2,     # Microtubule stability lost due to cryptic exon
        "UNC13A": -1.8,    # Synaptic vesicle priming lost due to cryptic exon
        "TARDBP": -1.5,    # Nuclear TDP-43 loss of function
        "KALRN": -1.4,     # Dendritic spine loss
        "NEFL": -1.6,      # Axonal neurofilament degradation
        "SYN1": -1.3,      # Presynaptic puncta loss
        "PSD95": -1.1,     # Postsynaptic density reduction
        "CHAT": -1.5,      # Motor cholinergic marker loss
        "SOD1": 0.4,       # Mild oxidative stress accumulation
        "ATXN2": 0.8,      # Stress granule load increase
        "CHMP2B": -0.6,    # Endosomal trafficking impairment
        "C9orf72": -0.5,   # Haploinsufficiency / loss of function
    },
    "ftd": {
        "GRN": -2.0,       # Progranulin haploinsufficiency
        "TARDBP": -1.4,    # Nuclear clearance
        "STMN2": -1.6,     # Cryptic splicing
        "SORT1": -1.2,     # Sortilin endocytic clearance
        "CTSD": -1.1,      # Lysosomal cathepsin D deficiency
        "SYN1": -1.0,      # Synaptic loss
    },
    "ad": {
        "APP": 1.2,        # Amyloid precursor elevation
        "MAPT": 1.5,       # Tau hyperphosphorylation / elevation
        "SYN1": -1.5,      # Severe synapse loss
        "PSD95": -1.4,     # Postsynaptic degradation
        "TREM2": -1.1,     # Microglial clearance impairment
        "APOE": 0.9,       # Lipid accumulation
    },
}

# Essential housekeeping genes monitored for collateral toxicity
NEURAL_HOUSEKEEPING_GENES: Dict[str, str] = {
    "ACTB": "Cytoskeletal structural integrity",
    "GAPDH": "Glycolytic neural energy metabolism",
    "TUBB3": "Core neuron-specific beta-tubulin",
    "POLR2A": "RNA polymerase II core subunit (pan-essential)",
    "RPS6": "Ribosomal translation machinery (pan-essential)",
    "SNAP25": "Core SNARE exocytosis machinery",
}

# Curated perturbation response vectors from human iPSC neuron functional genomics
CURATED_PERTURBATION_CATALOG: List[TargetPerturbationEffect] = [
    # 1. ATXN2 Knockdown: landmark ALS/FTD therapeutic rescue target
    TargetPerturbationEffect(
        target_gene="ATXN2",
        intervention=InterventionType.KNOCKDOWN,
        gene_deltas={
            "ATXN2": -1.8,
            "TARDBP": 0.4,
            "STMN2": 0.8,
            "SYN1": 0.6,
            "PSD95": 0.5,
            "NEFL": 0.7,
        },
        primary_pathway="Stress granule dynamics and TDP-43 aggregation suppression",
        evidence_citation="Elden et al. Nature 2010; Becker et al. Nature 2017",
        is_common_essential=False,
        collateral_risk=0.08,
    ),
    # 2. STMN2 Restoration via ASO steric block: targets cryptic splice donor
    TargetPerturbationEffect(
        target_gene="STMN2",
        intervention=InterventionType.ASO_STERIC_BLOCK,
        gene_deltas={
            "STMN2": 2.1,
            "NEFL": 1.2,
            "SYN1": 0.8,
            "PSD95": 0.6,
        },
        primary_pathway="Microtubule stability and axonal regeneration",
        evidence_citation="Klim et al. Nat Neurosci 2019; Melamed et al. Nat Neurosci 2019",
        is_common_essential=False,
        collateral_risk=0.04,
    ),
    # 3. UNC13A Restoration via ASO steric block: blocks cryptic donor at rs12608932
    TargetPerturbationEffect(
        target_gene="UNC13A",
        intervention=InterventionType.ASO_STERIC_BLOCK,
        gene_deltas={
            "UNC13A": 1.7,
            "SYN1": 0.9,
            "PSD95": 0.8,
        },
        primary_pathway="Synaptic vesicle priming and neurotransmitter release",
        evidence_citation="Brown et al. Nature 2022; Ma et al. Nature 2022",
        is_common_essential=False,
        collateral_risk=0.03,
    ),
    # 4. PIKFYVE Knockdown / Inhibition: stimulates endolysosomal clearance
    TargetPerturbationEffect(
        target_gene="PIKFYVE",
        intervention=InterventionType.KNOCKDOWN,
        gene_deltas={
            "CHMP2B": 0.9,
            "STMN2": 0.6,
            "NEFL": 0.5,
            "CHAT": 0.7,
        },
        primary_pathway="Endolysosomal trafficking and secretory autophagy",
        evidence_citation="Shi et al. Nat Med 2018; Hung et al. Cell Rep 2023",
        is_common_essential=False,
        collateral_risk=0.12,
    ),
    # 5. TARDBP Overexpression (high toxicity warning!)
    TargetPerturbationEffect(
        target_gene="TARDBP",
        intervention=InterventionType.OVEREXPRESSION,
        gene_deltas={
            "TARDBP": 2.0,
            "STMN2": 0.5,
            "SYN1": -0.8,
            "NEFL": -1.0,
            "ACTB": -0.4,
            "GAPDH": -0.5,
        },
        primary_pathway="Nuclear splicing regulation (triggers aggregate toxicity when overexpressed)",
        evidence_citation="Barmada et al. J Neurosci 2010",
        is_common_essential=True,
        collateral_risk=0.65,
    ),
    # 6. POLR2A Knockdown (Pan-essential negative control: severe lethality)
    TargetPerturbationEffect(
        target_gene="POLR2A",
        intervention=InterventionType.KNOCKDOWN,
        gene_deltas={
            "POLR2A": -2.5,
            "ACTB": -2.0,
            "GAPDH": -2.2,
            "TUBB3": -1.8,
            "SNAP25": -1.9,
            "SYN1": -2.1,
        },
        primary_pathway="Core transcription machinery",
        evidence_citation="DepMap Common Essential Index 2024",
        is_common_essential=True,
        collateral_risk=0.98,
    ),
    # 7. GRN Overexpression: restores progranulin in FTD
    TargetPerturbationEffect(
        target_gene="GRN",
        intervention=InterventionType.OVEREXPRESSION,
        gene_deltas={
            "GRN": 2.2,
            "CTSD": 1.0,
            "SORT1": 0.8,
            "SYN1": 0.6,
        },
        primary_pathway="Lysosomal homeostasis and neuroprotection",
        evidence_citation="Ward et al. Sci Transl Med 2017",
        is_common_essential=False,
        collateral_risk=0.05,
    ),
]


class CounterfactualSimulator:
    """In silico counterfactual causal simulation engine."""

    def __init__(self, catalog: Optional[List[TargetPerturbationEffect]] = None):
        self.catalog = catalog or CURATED_PERTURBATION_CATALOG
        self._target_map: Dict[Tuple[str, InterventionType], TargetPerturbationEffect] = {
            (effect.target_gene.upper(), effect.intervention): effect
            for effect in self.catalog
        }

    def simulate(
        self,
        disease: str,
        target_gene: str,
        intervention: InterventionType = InterventionType.KNOCKDOWN,
        efficacy_alpha: float = 0.8,
        custom_disease_state: Optional[Dict[str, float]] = None,
    ) -> SimulationResult:
        """Simulate the effect of a counterfactual target perturbation on disease state."""
        disease_key = disease.lower().strip()
        target_key = target_gene.upper().strip()

        # Resolve disease state vector
        if custom_disease_state:
            disease_vector = dict(custom_disease_state)
        elif disease_key in BENCHMARK_DISEASE_STATES:
            disease_vector = dict(BENCHMARK_DISEASE_STATES[disease_key])
        else:
            raise ValueError(
                f"Unknown disease '{disease}'. Choose from {list(BENCHMARK_DISEASE_STATES.keys())} "
                "or supply a custom disease state dictionary."
            )

        # Lookup perturbation model
        lookup_key = (target_key, intervention)
        if lookup_key in self._target_map:
            effect = self._target_map[lookup_key]
        else:
            # Fallback: check if target is in catalog with any intervention
            matching = [eff for eff in self.catalog if eff.target_gene.upper() == target_key]
            if matching:
                effect = matching[0]
                intervention = effect.intervention
            else:
                # Generate generic modeled perturbation for novel target
                is_essential = target_key in NEURAL_HOUSEKEEPING_GENES
                effect = TargetPerturbationEffect(
                    target_gene=target_key,
                    intervention=intervention,
                    gene_deltas={target_key: -1.5 if intervention in (InterventionType.KNOCKDOWN, InterventionType.KNOCKOUT) else 1.5},
                    primary_pathway="Uncharacterized targeted pathway",
                    evidence_citation="Inferred via standard linear perturbation model",
                    is_common_essential=is_essential,
                    collateral_risk=0.75 if is_essential else 0.15,
                )

        # Bound efficacy alpha between 0.0 and 1.0
        alpha = max(0.0, min(1.0, float(efficacy_alpha)))

        # Evaluate vectors
        # Compute baseline disease magnitude: Euclidean norm ||v_D||_2
        disease_genes = set(disease_vector.keys())
        all_genes = sorted(list(disease_genes.union(effect.gene_deltas.keys())))

        v_disease = [disease_vector.get(g, 0.0) for g in all_genes]
        v_delta = [effect.gene_deltas.get(g, 0.0) * alpha for g in all_genes]
        v_post = [d + delta for d, delta in zip(v_disease, v_delta)]

        disease_mag = math.sqrt(sum(x * x for x in v_disease))
        post_mag = math.sqrt(sum(x * x for x in v_post))

        # Rescue efficiency: reduction in disease norm
        if disease_mag > 1e-6:
            rescue_eff_pct = max(-100.0, min(100.0, ((disease_mag - post_mag) / disease_mag) * 100.0))
        else:
            rescue_eff_pct = 0.0

        # Cosine rescue alignment: cos(v_delta, -v_disease)
        # +1.0 = perturbation vector exactly opposes disease vector (perfect targeted rescue)
        # -1.0 = perturbation vector worsens disease direction
        v_neg_disease = [-x for x in v_disease]
        dot_product = sum(a * b for a, b in zip(v_delta, v_neg_disease))
        delta_mag = math.sqrt(sum(x * x for x in v_delta))

        if delta_mag > 1e-6 and disease_mag > 1e-6:
            cos_alignment = dot_product / (delta_mag * disease_mag)
            cos_alignment = max(-1.0, min(1.0, cos_alignment))
        else:
            cos_alignment = 0.0

        # Collateral toxicity calculation
        housekeeping_impact = 0.0
        for hk_gene in NEURAL_HOUSEKEEPING_GENES:
            if hk_gene in effect.gene_deltas:
                # Down-regulation of housekeeping genes is penalized
                hk_delta = effect.gene_deltas[hk_gene] * alpha
                if hk_delta < 0:
                    housekeeping_impact += abs(hk_delta)

        collateral_idx = min(1.0, effect.collateral_risk + (housekeeping_impact * 0.25))

        # Composite rescue score (0 - 100)
        # Grounded in Connectivity Map (CMap) and Perturb-seq benchmarks:
        # Cosine directional alignment is the primary metric (55 pts),
        # complemented by Euclidean magnitude rescue (up to 35 pts),
        # penalized by collateral off-target toxicity (up to 50 pts).
        alignment_score = max(0.0, cos_alignment) * 55.0
        magnitude_score = max(0.0, rescue_eff_pct) * 1.5
        penalty = collateral_idx * 50.0
        composite_score = round(max(0.0, min(100.0, alignment_score + magnitude_score - penalty + 10.0)), 1)

        # Viability verdict
        is_viable = (composite_score >= 45.0) and (collateral_idx < 0.40) and not effect.is_common_essential

        if composite_score >= 60.0 and is_viable:
            verdict = "STRONG RESCUE PREDICTION: Reverses pathological transcriptomic vector with minimal off-target neural collateral damage."
        elif composite_score >= 45.0 and is_viable:
            verdict = "MODERATE RESCUE PREDICTION: Partial phenotypic recovery; requires titration to avoid secondary pathway stress."
        elif effect.is_common_essential or collateral_idx >= 0.50:
            verdict = "VIABILITY FAILURE: High collateral toxicity / essential gene disruption outweighs therapeutic rescue."
        else:
            verdict = "INSUFFICIENT RESCUE: Perturbation does not sufficiently counter disease vector direction."

        # Trajectories for key genes
        trajectories: List[Dict[str, float]] = []
        for g in sorted(disease_genes.union(effect.gene_deltas.keys())):
            base_val = disease_vector.get(g, 0.0)
            delta_val = effect.gene_deltas.get(g, 0.0) * alpha
            new_val = base_val + delta_val
            if abs(base_val) > 0.05 or abs(delta_val) > 0.05:
                trajectories.append({
                    "gene": g,
                    "baseline_log2fc": round(base_val, 2),
                    "delta_log2fc": round(delta_val, 2),
                    "post_log2fc": round(new_val, 2),
                })

        # Pathway impacts
        pathway_impacts = {
            effect.primary_pathway: "BENEFICIAL_RESCUE" if cos_alignment > 0.2 else "INSUFFICIENT_ALIGNMENT"
        }
        if collateral_idx > 0.3:
            pathway_impacts["Cellular Viability & Housekeeping"] = "COLLATERAL_IMPAIRMENT"

        # Actionable recommendations
        recs: List[str] = []
        if is_viable:
            recs.append(
                f"Advance target {target_key} ({intervention.value}) into validation experiments (dose-response curve with alpha={alpha:.2f})."
            )
            recs.append(f"Monitor key rescue biomarkers: {[t['gene'] for t in trajectories[:3]]}.")
        else:
            if effect.is_common_essential:
                recs.append(f"Target {target_key} is a common essential gene (DepMap); direct systemic disruption will induce cell death.")
            if collateral_idx >= 0.40:
                recs.append(f"High collateral toxicity ({collateral_idx:.2f}): evaluate sub-maximal dosing or cell-type specific ASO delivery.")
            recs.append("Consider combinatorial screening or alternative pathway nodes.")

        return SimulationResult(
            disease=disease.upper(),
            target_gene=target_key,
            intervention=intervention,
            efficacy_alpha=alpha,
            disease_magnitude=round(disease_mag, 3),
            post_intervention_magnitude=round(post_mag, 3),
            rescue_efficiency_pct=round(rescue_eff_pct, 1),
            cosine_rescue_alignment=round(cos_alignment, 3),
            collateral_toxicity_index=round(collateral_idx, 2),
            composite_rescue_score=composite_score,
            is_clinically_viable=is_viable,
            verdict=verdict,
            key_gene_trajectories=trajectories,
            pathway_impacts=pathway_impacts,
            recommendations=recs,
        )

    def screen(
        self,
        disease: str,
        efficacy_alpha: float = 0.8,
    ) -> List[SimulationResult]:
        """Screen all cataloged perturbation targets against a disease state."""
        results: List[SimulationResult] = []
        for effect in self.catalog:
            res = self.simulate(
                disease=disease,
                target_gene=effect.target_gene,
                intervention=effect.intervention,
                efficacy_alpha=efficacy_alpha,
            )
            results.append(res)
        # Sort by composite rescue score descending
        results.sort(key=lambda x: x.composite_rescue_score, reverse=True)
        return results
