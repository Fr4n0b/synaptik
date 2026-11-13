"""
Empirical ASO Design Benchmark Against Published Wet-Lab Tiling Loci.

Evaluates ASODesigner's ability to prioritize clinically/experimentally validated
ASO target windows over flanking non-optimal windows across four published
neurogenetics loci:

  1. SMN2 intron 7 / ISS-N1 (Hua et al. AJHG 2008; Singh et al. MCB 2006)
     Steric-blocking 18-mer 2'-MOE/PS (Nusinersen / ISIS-396443) occluding
     the hnRNP A1/A2 intronic splicing silencer ISS-N1 (+10..+27).
  2. STMN2 cryptic exon 2a (Klim et al. Nat Neurosci 2019; Melamed et al. 2019;
     Baughn et al. Science 2023)
     Steric-blocking ASOs targeting the cryptic 3' splice acceptor / PPT /
     exon 2a junction to restore full-length STMN2 under TDP-43 loss.
  3. MAPT exon 10 5' splice-site stem-loop (Schoch et al. Neuron 2016)
     Structured hairpin regulating 4R/3R tau; tests whether 3D RNA folding
     accessibility (Delta G_open) distinguishes accessible vs occluded windows.
  4. SOD1 coding/3'UTR gapmer (Miller et al. Lancet Neurol 2013; Tofersen)
     5-10-5 MOE gapmer knockdown window free of G-quadruplexes and stable
     intramolecular hairpins.

Reports Hit-Rate@5 (Recall@5) and Mean Reciprocal Rank (MRR), plus an ablation
comparing sequence-only ranking vs 3D biophysics-coupled ranking.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .aso import ASODesigner, ASOMode


@dataclass
class ASOBenchmarkCase:
    gene: str
    locus: str
    citation: str
    mode: ASOMode
    target_sequence: str
    validated_start: int
    validated_end: int
    window_sizes: tuple[int, ...] = (18, 20)
    tolerance_nt: int = 8


@dataclass
class ASOLocusOutcome:
    gene: str
    locus: str
    citation: str
    best_rank_with_biophysics: int | None
    best_rank_sequence_only: int | None
    top_candidate_score: float
    top_delta_g_open: float | None
    hit_at_5_biophysics: bool
    hit_at_5_sequence_only: bool


@dataclass
class ASOBenchmarkReport:
    n_loci: int
    hit_rate_at_5_biophysics: float
    hit_rate_at_5_sequence_only: float
    mrr_biophysics: float
    mrr_sequence_only: float
    outcomes: list[ASOLocusOutcome] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "n_loci": self.n_loci,
            "hit_rate_at_5_biophysics": round(self.hit_rate_at_5_biophysics, 3),
            "hit_rate_at_5_sequence_only": round(self.hit_rate_at_5_sequence_only, 3),
            "mrr_biophysics": round(self.mrr_biophysics, 3),
            "mrr_sequence_only": round(self.mrr_sequence_only, 3),
            "outcomes": [
                {
                    "gene": o.gene,
                    "locus": o.locus,
                    "citation": o.citation,
                    "best_rank_with_biophysics": o.best_rank_with_biophysics,
                    "best_rank_sequence_only": o.best_rank_sequence_only,
                    "top_candidate_score": o.top_candidate_score,
                    "top_delta_g_open": o.top_delta_g_open,
                    "hit_at_5_biophysics": o.hit_at_5_biophysics,
                    "hit_at_5_sequence_only": o.hit_at_5_sequence_only,
                }
                for o in self.outcomes
            ],
        }


PUBLISHED_ASO_LOCI: list[ASOBenchmarkCase] = [
    ASOBenchmarkCase(
        gene="SMN2",
        locus="Exon 7 5'ss + Intron 7 ISS-N1",
        citation="Hua et al. AJHG 2008; Singh et al. MCB 2006 (Nusinersen locus)",
        mode=ASOMode.SPLICE_SWITCH,
        # Exon 7 end + Intron 7 (+1..+52) containing ISS-N1 (CCAGCAUUAUGAAAGUGAAU)
        target_sequence="GGTTTTAGACAAAATCAAAAAGAAGGAAGGTGCTCACATTCCTTAAATTAAGGAGTAAGTCTGCCAGCATTATGAAAGTGAAUCTTACTTTTGTAAAACTTTATGGTTTGT",
        validated_start=60,
        validated_end=84,
        window_sizes=(18,),
        tolerance_nt=8,
    ),
    ASOBenchmarkCase(
        gene="STMN2",
        locus="Cryptic Exon 2a 3'ss / PPT + Junction",
        citation="Klim et al. Nat Neurosci 2019; Melamed et al. Nat Neurosci 2019",
        mode=ASOMode.SPLICE_SWITCH,
        target_sequence="AAAAAGGCGCCGTTCTGTAACAGTTTTTGTTTTTTCAGGTGAGTATTCAGGTGGAGACGCACCGGAGCTCGACTTTCCAGGAAACGG",
        validated_start=20,
        validated_end=50,
        window_sizes=(18, 20),
        tolerance_nt=8,
    ),
    ASOBenchmarkCase(
        gene="MAPT",
        locus="Exon 10 5'ss Regulatory Stem-Loop",
        citation="Schoch et al. Neuron 2016",
        mode=ASOMode.SPLICE_SWITCH,
        # Exon 10 3' end + stem-loop + accessible flanking intron 10
        target_sequence="ACGTCCGGCAGTGTGGCAGTCAAGTCCAAGATCGGCTCCACTGAGAACCTGAAGCACCAGCCGGGAGGCGGGAAAGGTGAGTACCTTCACACGTCCCATGCGCCGTGCTGTG",
        validated_start=68,
        validated_end=94,
        window_sizes=(18, 20),
        tolerance_nt=10,
    ),
    ASOBenchmarkCase(
        gene="SOD1",
        locus="Exon 4 / Gapmer RNase H1 Window",
        citation="Miller et al. Lancet Neurol 2013 (Tofersen series)",
        mode=ASOMode.GAPMER,
        target_sequence="GGGGGGATGAAGAGAGGCATGTTGGAGACTTGGGCAATGTGACTGCTGACAAAGATGGTGTGGCCGATGTGTCTATTGAAGATTCTGTGATCTCACTCTCAGGAGACCATTG",
        validated_start=45,
        validated_end=85,
        window_sizes=(20,),
        tolerance_nt=10,
    ),
]


def _overlaps_validated(start: int, end: int, case: ASOBenchmarkCase) -> bool:
    """Return True if candidate [start, end) substantially overlaps the validated locus."""
    ov = min(end, case.validated_end + case.tolerance_nt) - max(start, case.validated_start - case.tolerance_nt)
    return ov >= min(12, end - start)


def run_aso_tiling_benchmark(
    cases: list[ASOBenchmarkCase] | None = None,
) -> ASOBenchmarkReport:
    """Run ASODesigner with and without 3D biophysics across published ASO loci."""
    loci = cases or PUBLISHED_ASO_LOCI
    designer = ASODesigner()
    outcomes: list[ASOLocusOutcome] = []

    for case in loci:
        res_bio = designer.design(
            target_sequence=case.target_sequence,
            mode=case.mode,
            target_gene=case.gene,
            window_sizes=case.window_sizes,
            top_k=10,
            use_biophysics=True,
        )
        res_seq = designer.design(
            target_sequence=case.target_sequence,
            mode=case.mode,
            target_gene=case.gene,
            window_sizes=case.window_sizes,
            top_k=10,
            use_biophysics=False,
        )

        rank_bio = next(
            (idx + 1 for idx, c in enumerate(res_bio.candidates) if _overlaps_validated(c.target_start, c.target_end, case)),
            None,
        )
        rank_seq = next(
            (idx + 1 for idx, c in enumerate(res_seq.candidates) if _overlaps_validated(c.target_start, c.target_end, case)),
            None,
        )

        top_c = res_bio.candidates[0] if res_bio.candidates else None
        outcomes.append(
            ASOLocusOutcome(
                gene=case.gene,
                locus=case.locus,
                citation=case.citation,
                best_rank_with_biophysics=rank_bio,
                best_rank_sequence_only=rank_seq,
                top_candidate_score=top_c.composite_score if top_c else 0.0,
                top_delta_g_open=top_c.delta_g_open_kcal_mol if top_c else None,
                hit_at_5_biophysics=(rank_bio is not None and rank_bio <= 5),
                hit_at_5_sequence_only=(rank_seq is not None and rank_seq <= 5),
            )
        )

    n = max(len(outcomes), 1)
    hit5_bio = sum(1 for o in outcomes if o.hit_at_5_biophysics) / n
    hit5_seq = sum(1 for o in outcomes if o.hit_at_5_sequence_only) / n
    mrr_bio = sum((1.0 / o.best_rank_with_biophysics) if o.best_rank_with_biophysics else 0.0 for o in outcomes) / n
    mrr_seq = sum((1.0 / o.best_rank_sequence_only) if o.best_rank_sequence_only else 0.0 for o in outcomes) / n

    return ASOBenchmarkReport(
        n_loci=len(outcomes),
        hit_rate_at_5_biophysics=hit5_bio,
        hit_rate_at_5_sequence_only=hit5_seq,
        mrr_biophysics=mrr_bio,
        mrr_sequence_only=mrr_seq,
        outcomes=outcomes,
    )
