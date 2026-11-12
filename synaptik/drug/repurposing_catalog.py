"""Curated CNS-active Small-Molecule Benchmark Repository & Offline Tractability Engine.

Maps neurodegenerative disease targets (*STMN2*, *UNC13A*, *TARDBP*, *SOD1*,
*SMN2*, *MTOR*, *HDAC1*, *GRIN1*, etc.) to approved, investigational, and tool
compounds with physicochemical CNS MPO scoring, blood-brain barrier penetrance,
and directional mechanism-of-action compatibility.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

from ..schema import Direction
from .repurposing import CNSProfile, DrugCandidate, MoADirection, cns_mpo, direction_compatible


@dataclass
class RepurposedCompound:
    name: str
    chembl_id: str
    target_symbol: str
    target_name: str
    action_type: str
    moa_direction: MoADirection
    max_clinical_phase: int  # 4 = Approved, 3 = Phase III, 2 = Phase II, 1 = Phase I, 0 = Preclinical
    cns_mpo_score: float
    bbb_plausible: bool
    efflux_risk: bool
    mw_da: float
    clogp: float
    hbd: int
    tpsa_a2: float
    indications: list[str]
    mechanism_summary: str
    literature_precedent: str
    predicted_rescue_efficacy_pct: float = 0.0
    translational_rank_score: float = 0.0

    @property
    def phase_label(self) -> str:
        labels = {
            4: "FDA/EMA Approved",
            3: "Phase III Clinical",
            2: "Phase II Clinical",
            1: "Phase I Clinical",
            0: "Preclinical / Tool Compound",
        }
        return labels.get(self.max_clinical_phase, "Investigational")


# Curated benchmark database of brain-penetrant small-molecule drugs
BENCHMARK_CNS_DRUGS: list[dict[str, Any]] = [
    {
        "name": "Risdiplam (Evrysdi)",
        "chembl_id": "CHEMBL4297495",
        "target_symbol": "SMN2",
        "target_name": "Survival of Motor Neuron 2 (Exon 7 Splicing)",
        "action_type": "SPLICING MODULATOR",
        "moa_direction": MoADirection.ACTIVATOR,
        "max_clinical_phase": 4,
        "mw_da": 401.5,
        "clogp": 1.8,
        "hbd": 1,
        "tpsa_a2": 66.8,
        "efflux_risk": False,
        "indications": ["Spinal Muscular Atrophy (Type 1-3)"],
        "mechanism_summary": "Oral central-active small molecule that binds specifically to the pre-mRNA 5' splice site of exon 7, displacing hnRNP G and promoting U1 snRNP recruitment to restore full-length functional SMN protein.",
        "literature_precedent": "Ratni et al. J Med Chem 2018; FIREFISH & SUNFISH NEJM 2021 trials.",
    },
    {
        "name": "Branaplam (LMI070)",
        "chembl_id": "CHEMBL3989916",
        "target_symbol": "SMN2",
        "target_name": "Survival of Motor Neuron 2 / HTT Splicing",
        "action_type": "SPLICING MODULATOR",
        "moa_direction": MoADirection.ACTIVATOR,
        "max_clinical_phase": 2,
        "mw_da": 393.4,
        "clogp": 2.2,
        "hbd": 1,
        "tpsa_a2": 72.4,
        "efflux_risk": False,
        "indications": ["SMA", "Huntington's Disease (mHTT lowering)"],
        "mechanism_summary": "Promotes exon inclusion in SMN2 and targets pseudo-exons in huntingtin mRNA, demonstrating brain-wide alternative splicing modification.",
        "literature_precedent": "Palacino et al. Nat Chem Biol 2015; Novartis clinical program.",
    },
    {
        "name": "Riluzole (Rilutek)",
        "chembl_id": "CHEMBL744",
        "target_symbol": "SLC1A2",
        "target_name": "Glutamate Transporter 1 (GLT-1 / EAAT2)",
        "action_type": "MODULATOR",
        "moa_direction": MoADirection.MODULATOR,
        "max_clinical_phase": 4,
        "mw_da": 234.2,
        "clogp": 3.2,
        "hbd": 1,
        "tpsa_a2": 44.1,
        "efflux_risk": False,
        "indications": ["Amyotrophic Lateral Sclerosis"],
        "mechanism_summary": "Inhibits presynaptic glutamate release and blocks voltage-gated sodium channels, attenuating excitotoxic motor neuron degeneration.",
        "literature_precedent": "Bensimon et al. NEJM 1994; FDA approved for ALS standard of care.",
    },
    {
        "name": "Edaravone (Radicava)",
        "chembl_id": "CHEMBL183",
        "target_symbol": "GPX1",
        "target_name": "Glutathione Peroxidase / Free Radical Quenching",
        "action_type": "ANTIOXIDANT / SCAVENGER",
        "moa_direction": MoADirection.ACTIVATOR,
        "max_clinical_phase": 4,
        "mw_da": 174.2,
        "clogp": 1.4,
        "hbd": 0,
        "tpsa_a2": 32.7,
        "efflux_risk": False,
        "indications": ["Amyotrophic Lateral Sclerosis", "Acute Ischemic Stroke"],
        "mechanism_summary": "Lipophilic free radical scavenger that mitigates oxidative stress-induced mitochondrial dysfunction and motor neuron death.",
        "literature_precedent": "Writing Group, Edaravone ALS Study Group Lancet Neurol 2017.",
    },
    {
        "name": "Sodium Phenylbutyrate",
        "chembl_id": "CHEMBL1200547",
        "target_symbol": "HSPA5",
        "target_name": "BiP / Endoplasmic Reticulum Chaperone & HDAC",
        "action_type": "CHEMICAL CHAPERONE / HDAC INHIBITOR",
        "moa_direction": MoADirection.INHIBITOR,
        "max_clinical_phase": 4,
        "mw_da": 186.2,
        "clogp": 1.7,
        "hbd": 0,
        "tpsa_a2": 40.1,
        "efflux_risk": False,
        "indications": ["ALS (CENTAUR trial)", "Urea Cycle Disorders"],
        "mechanism_summary": "Alleviates endoplasmic reticulum stress and unfolded protein response (UPR) induced by cytoplasmic protein aggregates.",
        "literature_precedent": "Paganoni et al. NEJM 2020 (CENTAUR trial).",
    },
    {
        "name": "Memantine (Namenda)",
        "chembl_id": "CHEMBL703",
        "target_symbol": "GRIN2A",
        "target_name": "NMDA Receptor Subunit 2A / 1",
        "action_type": "UNCOMPETITIVE ANTAGONIST",
        "moa_direction": MoADirection.INHIBITOR,
        "max_clinical_phase": 4,
        "mw_da": 179.3,
        "clogp": 2.7,
        "hbd": 1,
        "tpsa_a2": 26.0,
        "efflux_risk": False,
        "indications": ["Alzheimer's Disease", "Frontotemporal Dementia (Investigational)"],
        "mechanism_summary": "Low-affinity, uncompetitive NMDA receptor open-channel blocker preserving physiological synaptic transmission while attenuating chronic calcium overload.",
        "literature_precedent": "Reisberg et al. NEJM 2003; Tariot et al. JAMA 2004.",
    },
    {
        "name": "Diazepam (Valium)",
        "chembl_id": "CHEMBL12",
        "target_symbol": "GABRA1",
        "target_name": "GABA-A Receptor Alpha-1 Subunit",
        "action_type": "POSITIVE ALLOSTERIC MODULATOR",
        "moa_direction": MoADirection.ACTIVATOR,
        "max_clinical_phase": 4,
        "mw_da": 284.7,
        "clogp": 2.8,
        "hbd": 0,
        "tpsa_a2": 32.7,
        "efflux_risk": False,
        "indications": ["Anxiety", "Muscle Spasticity", "Status Epilepticus"],
        "mechanism_summary": "Allosterically potentiates GABA-mediated chloride conductance, damping hyperexcitable motor neuron circuits.",
        "literature_precedent": "Classical benzodiazepine pharmacology; gold standard for CNS penetration.",
    },
    {
        "name": "Vorinostat (SAHA)",
        "chembl_id": "CHEMBL98",
        "target_symbol": "HDAC1",
        "target_name": "Histone Deacetylase 1 / 2 / 6",
        "action_type": "INHIBITOR",
        "moa_direction": MoADirection.INHIBITOR,
        "max_clinical_phase": 4,
        "mw_da": 264.3,
        "clogp": 1.9,
        "hbd": 3,
        "tpsa_a2": 78.4,
        "efflux_risk": False,
        "indications": ["Cutaneous T-Cell Lymphoma", "Neuroepigenetic Gene Reactivation"],
        "mechanism_summary": "Pan-HDAC inhibitor promoting histone acetylation, chromatin decompaction, and transcription of neuroprotective loci.",
        "literature_precedent": "Marks & Breslow Nat Biotechnol 2007; Hockly et al. PNAS 2003.",
    },
    {
        "name": "Rapamycin (Sirolimus)",
        "chembl_id": "CHEMBL410",
        "target_symbol": "MTOR",
        "target_name": "Mechanistic Target of Rapamycin Complex 1 (mTORC1)",
        "action_type": "ALLOSTERIC INHIBITOR",
        "moa_direction": MoADirection.INHIBITOR,
        "max_clinical_phase": 4,
        "mw_da": 914.2,
        "clogp": 4.3,
        "hbd": 3,
        "tpsa_a2": 195.0,
        "efflux_risk": True,
        "indications": ["Immunosuppression", "Tuberous Sclerosis", "Autophagy Induction"],
        "mechanism_summary": "Binds FKBP12 to inhibit mTORC1, stimulating macroautophagy and clearance of toxic cytoplasmic TDP-43 and tau aggregates.",
        "literature_precedent": "Bove et al. Nat Med 2011; widely used in in vitro neuronal models.",
    },
    {
        "name": "SRPIN340",
        "chembl_id": "CHEMBL519827",
        "target_symbol": "SRPK1",
        "target_name": "Serine/Arginine-Rich Protein Kinase 1",
        "action_type": "ATP-COMPETITIVE INHIBITOR",
        "moa_direction": MoADirection.INHIBITOR,
        "max_clinical_phase": 0,
        "mw_da": 349.4,
        "clogp": 2.6,
        "hbd": 1,
        "tpsa_a2": 55.8,
        "efflux_risk": False,
        "indications": ["Alternative Splicing Regulation (Preclinical)"],
        "mechanism_summary": "Selective SRPK1 kinase inhibitor that blocks phosphorylation of SR splicing factors (SRSF1/2), altering alternative splice site selection.",
        "literature_precedent": "Fukuhara et al. PNAS 2006; Amin et al. Cancer Cell 2011.",
    },
    {
        "name": "TG003",
        "chembl_id": "CHEMBL276707",
        "target_symbol": "CLK1",
        "target_name": "CDC-Like Kinase 1 / 4",
        "action_type": "ATP-COMPETITIVE INHIBITOR",
        "moa_direction": MoADirection.INHIBITOR,
        "max_clinical_phase": 0,
        "mw_da": 249.3,
        "clogp": 2.1,
        "hbd": 1,
        "tpsa_a2": 54.1,
        "efflux_risk": False,
        "indications": ["Splice Switching Chemical Tool"],
        "mechanism_summary": "Inhibits CLK family kinases to prevent SR protein hyperphosphorylation, modulating alternative exon skipping and cryptic splice activation.",
        "literature_precedent": "Muraki et al. J Biol Chem 2004; Araki et al. Science 2020.",
    },
]


def screen_repurposing_candidates(
    target_gene: str,
    required_direction: MoADirection | str = MoADirection.ACTIVATOR,
    min_clinical_phase: int = 0,
    bbb_only: bool = False,
) -> list[RepurposedCompound]:
    """Identify and rank small-molecule candidates matching target and direction."""
    gene_upper = target_gene.strip().upper()

    if isinstance(required_direction, str):
        try:
            req_dir = MoADirection(required_direction.lower())
        except ValueError:
            req_dir = MoADirection.UNKNOWN
    else:
        req_dir = required_direction

    results: list[RepurposedCompound] = []

    # Target gene synonym/pathway expansion
    relevant_targets = {gene_upper}
    if gene_upper in ["STMN2", "UNC13A", "TARDBP"]:
        relevant_targets.update(["SMN2", "SRPK1", "CLK1", "HSPA5", "SLC1A2", "GPX1"])
    elif gene_upper in ["SOD1", "FUS"]:
        relevant_targets.update(["GPX1", "HSPA5", "SLC1A2", "MTOR"])
    elif gene_upper in ["SMN2", "ELP1", "COL6A1"]:
        relevant_targets.update(["SMN2", "SRPK1", "CLK1"])

    for item in BENCHMARK_CNS_DRUGS:
        t_sym = item["target_symbol"].upper()
        if t_sym not in relevant_targets and gene_upper not in t_sym:
            continue

        phase = item["max_clinical_phase"]
        if phase < min_clinical_phase:
            continue

        # Evaluate CNS MPO
        mpo_res, bbb_res = cns_mpo(
            mw=item["mw_da"],
            logp=item["clogp"],
            tpsa=item["tpsa_a2"],
            hbd=item["hbd"],
        )
        mpo_score = mpo_res or 0.0
        bbb_plausible = (mpo_score >= 4.0 and not item["efflux_risk"]) if bbb_res is None else bbb_res

        if bbb_only and not bbb_plausible:
            continue

        moa = item["moa_direction"]
        dir_ok, rationale = direction_compatible(req_dir, moa)

        # Base rescue percentage heuristic
        base_rescue = 35.0
        if t_sym == gene_upper:
            base_rescue += 30.0
        if phase == 4:
            base_rescue += 15.0
        elif phase >= 2:
            base_rescue += 10.0
        if bbb_plausible:
            base_rescue += 10.0

        rescue_pct = round(min(88.0, base_rescue), 1)

        # Composite Translational Rank Score (0 - 100)
        rank_score = (
            (phase / 4.0) * 35.0  # Phase maturity
            + (mpo_score / 6.0) * 35.0  # CNS druggability
            + (15.0 if dir_ok else 5.0)  # Directional alignment
            + (15.0 if not item["efflux_risk"] else 0.0)  # P-gp safety
        )

        cand = RepurposedCompound(
            name=item["name"],
            chembl_id=item["chembl_id"],
            target_symbol=item["target_symbol"],
            target_name=item["target_name"],
            action_type=item["action_type"],
            moa_direction=moa,
            max_clinical_phase=phase,
            cns_mpo_score=round(mpo_score, 2),
            bbb_plausible=bbb_plausible,
            efflux_risk=item["efflux_risk"],
            mw_da=item["mw_da"],
            clogp=item["clogp"],
            hbd=item["hbd"],
            tpsa_a2=item["tpsa_a2"],
            indications=item["indications"],
            mechanism_summary=item["mechanism_summary"],
            literature_precedent=item["literature_precedent"],
            predicted_rescue_efficacy_pct=rescue_pct,
            translational_rank_score=round(rank_score, 1),
        )
        results.append(cand)

    # Sort descending by translational rank score
    results.sort(key=lambda x: -x.translational_rank_score)
    return results


def repurposed_compound_to_dict(c: RepurposedCompound) -> dict[str, Any]:
    return {
        "name": c.name,
        "chembl_id": c.chembl_id,
        "target_symbol": c.target_symbol,
        "target_name": c.target_name,
        "action_type": c.action_type,
        "moa_direction": c.moa_direction.value,
        "max_clinical_phase": c.max_clinical_phase,
        "phase_label": c.phase_label,
        "cns_mpo_score": c.cns_mpo_score,
        "bbb_plausible": c.bbb_plausible,
        "efflux_risk": c.efflux_risk,
        "mw_da": c.mw_da,
        "clogp": c.clogp,
        "hbd": c.hbd,
        "tpsa_a2": c.tpsa_a2,
        "indications": c.indications,
        "mechanism_summary": c.mechanism_summary,
        "literature_precedent": c.literature_precedent,
        "predicted_rescue_efficacy_pct": c.predicted_rescue_efficacy_pct,
        "translational_rank_score": c.translational_rank_score,
    }
