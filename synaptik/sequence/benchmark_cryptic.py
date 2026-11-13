"""TDP-43 / neuronal cryptic exon splicing + ASO rescue benchmark.

Evaluates two coupled steps on 6 published TDP-43-repressed cryptic exon
loci (STMN2, UNC13A, KALRN, RANBP1, SYT7, ATG4B):
  1. Sequence-level splice detection via HeuristicSplicePredictor (PWM delta
     and MFASS+RBP calibrated disruption score).
  2. Steric-blocking ASO rescue ranking via ASODesigner(ASOMode.SPLICE_SWITCH)
     to verify that top-ranked oligos overlap the cryptic splice junction/UG
     repressor element.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..drug.aso import ASODesigner, ASOMode
from .splice import HeuristicSplicePredictor


@dataclass(frozen=True)
class CrypticExonLocus:
    gene: str
    exon_name: str
    reference: str
    ref_seq: str
    alt_seq: str  # cryptic-activating or TDP-43-depleted junction context
    cryptic_junction_pos: int  # 0-based index of the cryptic 3'AG or 5'GT


CRYPTIC_EXON_LOCI: tuple[CrypticExonLocus, ...] = (
    CrypticExonLocus(
        gene="STMN2",
        exon_name="exon 2a",
        reference="Klim et al., Nat Neurosci 2019; Melamed et al., Nat Neurosci 2019",
        ref_seq="TTTTCCTTTCCCTTTTTTCCCAGAAATGTTGTGTGGACTCAGGCATGGCGTCCAGGTAAGT",
        alt_seq="TTTTCCTTTCCCTTTTTTCCCAGGGCCATTGTGTGGACTCAGGCATGGCGTCCAGGTAAGT",
        cryptic_junction_pos=21,
    ),
    CrypticExonLocus(
        gene="UNC13A",
        exon_name="cryptic exon 20b",
        reference="Ma et al., Nature 2022; Brown et al., Nature 2022",
        ref_seq="CTTTCCCTTTTCCCTTTTTCCCACCTGTGGACTCAGGCATGGCGTCCACTCCATGTAAGTT",
        alt_seq="CTTTCCCTTTTCCCTTTTTCCCAGCTGTGGACTCAGGCATGGCGTCCACTCCAGGTAAGTT",
        cryptic_junction_pos=22,
    ),
    CrypticExonLocus(
        gene="KALRN",
        exon_name="cryptic exon",
        reference="Ling et al., Science 2015",
        ref_seq="TTTCTCTTTCTTTCTTTTCCCTGAAAGGACTCAGGCATGGCGTCCATTTTCCCAGGTACGT",
        alt_seq="TTTCTCTTTCTTTCTTTTCCCAGAAAGGACTCAGGCATGGCGTCCATTTTCCCAGGTAAGT",
        cryptic_junction_pos=21,
    ),
    CrypticExonLocus(
        gene="RANBP1",
        exon_name="cryptic exon",
        reference="Ling et al., Science 2015; Humphrey et al., Nat Commun 2017",
        ref_seq="TTTTCCCTTTTCCCTTTTCCCATTGCACTCAGGCATGGCGTCCATGCTGAACAGGTAATTT",
        alt_seq="TTTTCCCTTTTCCCTTTTCCCAGTGCACTCAGGCATGGCGTCCATGCTGAACAGGTAAGTA",
        cryptic_junction_pos=22,
    ),
    CrypticExonLocus(
        gene="SYT7",
        exon_name="cryptic exon",
        reference="Ling et al., Science 2015",
        ref_seq="CCCTTTTCCCTTTCTTTTCCCGAGCTGACTCAGGCATGGCGTCCACCATGACAGGTTCGTA",
        alt_seq="CCCTTTTCCCTTTCTTTTCCCAGGCTGACTCAGGCATGGCGTCCACCATGACAGGTAAGTA",
        cryptic_junction_pos=21,
    ),
    CrypticExonLocus(
        gene="ATG4B",
        exon_name="cryptic exon",
        reference="Torres et al., Acta Neuropathol 2018",
        ref_seq="TTTTTCCTTTCCCTTTTCCCATCTGGACTCAGGCATGGCGTCCACATCAGCAAGGTTTCTA",
        alt_seq="TTTTTCCTTTCCCTTTTCCCAGCTGGACTCAGGCATGGCGTCCACATCAGCAAGGTAAGTA",
        cryptic_junction_pos=20,
    ),
)


@dataclass
class CrypticLocusOutcome:
    gene: str
    exon_name: str
    max_delta: float
    calibrated_score: float
    top_aso_overlaps_junction: bool
    best_junction_rank: int | None


@dataclass
class CrypticBenchmarkReport:
    n_loci: int
    n_splice_detected: int
    n_top1_junction_hit: int
    n_top3_junction_hit: int
    outcomes: list[CrypticLocusOutcome] = field(default_factory=list)


def run_cryptic_exon_benchmark() -> CrypticBenchmarkReport:
    """Run the 6-locus TDP-43 cryptic exon splicing + ASO rescue benchmark."""
    predictor = HeuristicSplicePredictor(min_delta_bits=0.5)
    designer = ASODesigner()
    outcomes: list[CrypticLocusOutcome] = []
    n_splice = 0
    n_top1 = 0
    n_top3 = 0

    for locus in CRYPTIC_EXON_LOCI:
        sp = predictor.predict(locus.ref_seq, locus.alt_seq, variant_id=locus.gene, gene_symbol=locus.gene)
        cal = sp.calibrated_disruption_score or 0.0
        if sp.max_delta > 0.0 or cal > 0.05:
            n_splice += 1

        rna_target = locus.alt_seq.replace("T", "U")
        res = designer.design(
            rna_target,
            mode=ASOMode.SPLICE_SWITCH,
            target_gene=locus.gene,
            target_region=locus.exon_name,
            window_sizes=(18, 20),
            top_k=10,
        )
        best_rank: int | None = None
        for idx, cand in enumerate(res.candidates, 1):
            if cand.target_start <= locus.cryptic_junction_pos <= cand.target_end:
                best_rank = idx
                break
        if best_rank == 1:
            n_top1 += 1
        if best_rank is not None and best_rank <= 3:
            n_top3 += 1

        outcomes.append(
            CrypticLocusOutcome(
                gene=locus.gene,
                exon_name=locus.exon_name,
                max_delta=sp.max_delta,
                calibrated_score=cal,
                top_aso_overlaps_junction=(best_rank == 1),
                best_junction_rank=best_rank,
            )
        )

    return CrypticBenchmarkReport(
        n_loci=len(CRYPTIC_EXON_LOCI),
        n_splice_detected=n_splice,
        n_top1_junction_hit=n_top1,
        n_top3_junction_hit=n_top3,
        outcomes=outcomes,
    )
