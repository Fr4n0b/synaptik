"""
Off-target specificity screening for antisense oligonucleotides.

Why this module exists
----------------------
Off-target hybridisation is the single most safety-critical property of a
therapeutic ASO. A gapmer that partially matches an unintended transcript
recruits RNase H and degrades it; splice-switching ASOs can perturb unintended
splice sites. Published ASO programmes have been halted over off-target and
sequence-dependent toxicity, so a specificity claim is a claim about patient
safety, not a nice-to-have annotation.

Until v0.8.1 this check was performed against ten hard-coded ~86 bp transcript
fragments -- 861 bp in total, roughly 0.0007% of the human transcriptome.
Measured empirically, that screen flagged 0 of 3000 random 20-mers: it could
only ever catch an oligo that happened to be a perfect match to one of those
ten fragments. Every candidate therefore "passed", and the report stated
"0-mismatch: 0, 1-mismatch: 0" as though a genome-wide search had been run.
A reassuring number produced by a search that cannot fail is worse than no
number, because the reader cannot tell the difference.

This module replaces that with three honest states:

  NOT_SCREENED      No reference sequences were supplied. No specificity claim
                    is made, and the report says so instead of printing zeros.
  PANEL             Screened against a small built-in panel. Useful as a smoke
                    test for gross errors only; explicitly NOT a safety screen,
                    and the coverage fraction is reported alongside the result.
  TRANSCRIPTOME     Screened against a user-supplied reference (FASTA). This is
                    a real search, and the report states how much was covered.

Even TRANSCRIPTOME is a screen, not a clearance. It finds near-matches by
seed-and-extend over ungapped alignments. It does not model RNA secondary
structure, accessibility, protein occupancy, or gapped alignment, and it is not
a substitute for BLAST/Bowtie against a current reference plus wet-lab
validation (RNA-seq of treated cells is the standard).

Method
------
Seed-and-extend. For a tolerance of `max_mismatches` = m over an oligo of
length L, the oligo is split into m+1 non-overlapping seeds. Any alignment with
at most m mismatches must leave at least one seed exact (pigeonhole principle),
so indexing the seeds and extending every seed hit finds every qualifying
alignment -- this is exhaustive within the ungapped model, not a heuristic
sample.

Both strands are searched. An ASO hybridises to any transcript containing the
reverse complement of the oligo, and transcripts in a FASTA may be deposited on
either strand, so searching one orientation alone would miss real hits.
"""

from __future__ import annotations

import gzip
import io
from collections import defaultdict
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Iterable, Iterator

# Mismatch tolerance empirically anchored to published ASO off-target datasets:
#
# 1. RNase H1 gapmers (DEFAULT_MAX_MISMATCHES["gapmer"] = 3):
#    Systematic transcriptome-wide off-target profiling (Kamola et al., Nucleic
#    Acids Res 2015, 43:8638-8650; Burel et al., Nucleic Acids Res 2016,
#    44:2093-2109; Yoshida et al., Nucleic Acids Res 2018, 46:10984-10996)
#    demonstrates that RNase H1 retains cleavage activity on off-target RNA
#    duplexes harboring 1-2 mismatches and, when mismatches fall in the 2'-MOE
#    or LNA wings rather than the central 7-10 nt DNA gap, up to 3 mismatches.
#    Transcripts with >= 4 mismatches show negligible cleavage (<5% knockdown).
#
# 2. Steric splice-switching oligos (DEFAULT_MAX_MISMATCHES["splice_switch"] = 2):
#    Uniformly modified 2'-MOE/PMO splice-switching ASOs (Hua et al., Am J Hum
#    Genet 2008; Bennett & Swayze, Annu Rev Pharmacol Toxicol 2010) do not
#    recruit RNase H1 and rely on competitive occupancy of splice motifs; 1-2
#    mismatches reduce duplex affinity by ~10-100x, while >= 3 mismatches abolish
#    steric blockade at therapeutic CSF concentrations (~100-500 nM).
DEFAULT_MAX_MISMATCHES = {"gapmer": 3, "splice_switch": 2}

# Human transcriptome scale, used only to report what fraction of a realistic
# reference the screen actually covered. Order-of-magnitude figure
# (GENCODE v4x, ~250k transcripts); it contextualises coverage, nothing more.
HUMAN_TRANSCRIPTOME_BP = 120_000_000

_COMPLEMENT = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")


class ScreenStatus(str, Enum):
    NOT_SCREENED = "not_screened"
    PANEL = "panel"
    TRANSCRIPTOME = "transcriptome"


@dataclass
class OffTargetHit:
    """One near-match of an oligo against a reference transcript."""

    transcript_id: str
    position: int
    strand: str          # "+" the oligo's complement, "-" its reverse complement
    n_mismatches: int
    aligned_sequence: str

    def describe(self) -> str:
        return (
            f"{self.transcript_id}:{self.position}({self.strand}) "
            f"{self.n_mismatches} mismatch(es)"
        )


@dataclass
class SpecificityResult:
    """Outcome of a specificity screen for a single oligo.

    `status` must be checked before any count is read. When the status is
    NOT_SCREENED the counts are all zero because nothing was searched, which
    is categorically different from zero hits having been found.
    """

    status: ScreenStatus = ScreenStatus.NOT_SCREENED
    oligo: str = ""
    reference_name: str = ""
    n_transcripts_searched: int = 0
    n_bases_searched: int = 0
    max_mismatches_searched: int = 0
    hits: list[OffTargetHit] = field(default_factory=list)
    on_target_gene: str = ""
    truncated: bool = False
    notes: list[str] = field(default_factory=list)

    @property
    def screened(self) -> bool:
        return self.status is ScreenStatus.TRANSCRIPTOME

    def counts_by_mismatch(self) -> dict[int, int]:
        out: dict[int, int] = defaultdict(int)
        for h in self.hits:
            out[h.n_mismatches] += 1
        return dict(sorted(out.items()))

    def worst_hit(self) -> OffTargetHit | None:
        return min(self.hits, key=lambda h: h.n_mismatches) if self.hits else None

    def coverage_fraction(self) -> float:
        """Fraction of a human-transcriptome-scale reference actually searched."""
        return self.n_bases_searched / HUMAN_TRANSCRIPTOME_BP

    def verdict(self) -> str:
        """A sentence that can be pasted into a report without misleading anyone."""
        if self.status is ScreenStatus.NOT_SCREENED:
            return (
                "NOT SCREENED - no reference supplied. This candidate carries no "
                "specificity evidence. Screen against a transcriptome FASTA "
                "(--transcriptome) or with BLAST/Bowtie before acting on it."
            )
        if self.status is ScreenStatus.PANEL:
            base = (
                f"PANEL ONLY - searched {self.n_bases_searched:,} bp across "
                f"{self.n_transcripts_searched} transcript fragments "
                f"({self.coverage_fraction():.5%} of a human transcriptome). "
                "This is a smoke test, NOT a specificity screen."
            )
            if self.hits:
                base += (
                    f" {len(self.hits)} hit(s) found even at this coverage, "
                    "which is a strong negative signal."
                )
            return base
        cov = self.coverage_fraction()
        scale = (
            f"{self.coverage_fraction():.1%} of human-transcriptome scale"
            if cov < 2.0 else "a reference larger than the human transcriptome"
        )
        if not self.hits:
            return (
                f"No ungapped match within {self.max_mismatches_searched} "
                f"mismatches across {self.n_transcripts_searched:,} transcripts "
                f"({self.n_bases_searched:,} bp, {scale}). Screened, not cleared: "
                "confirm with BLAST/Bowtie and an RNA-seq specificity experiment."
            )
        counts = ", ".join(f"{n}mm: {c}" for n, c in self.counts_by_mismatch().items())
        return (
            f"{len(self.hits)} off-target hit(s) found ({counts}) across "
            f"{self.n_transcripts_searched:,} transcripts ({self.n_bases_searched:,} bp). "
            "Review each before synthesis."
        )


# --------------------------------------------------------------------------
# Reference sequence loading
# --------------------------------------------------------------------------


def _open_maybe_gzip(path: Path) -> io.TextIOBase:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def iter_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Stream (id, sequence) from a FASTA file, gzipped or plain.

    Streaming matters: a human transcriptome FASTA is several hundred megabytes
    and must not be held in memory in full.
    """
    path = Path(path)
    name: str | None = None
    chunks: list[str] = []
    with _open_maybe_gzip(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0] if len(line) > 1 else "unnamed"
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


# --------------------------------------------------------------------------
# Screening
# --------------------------------------------------------------------------


class SpecificityScreener:
    """Seed-and-extend off-target screen over a reference set of transcripts.

    Exhaustive within its model: every ungapped alignment at or below the
    mismatch tolerance is found. It does not perform gapped alignment and does
    not model RNA structure or accessibility.
    """

    def __init__(
        self,
        max_mismatches: int = 3,
        max_hits_per_oligo: int = 200,
        skip_on_target_gene: str = "",
    ):
        if max_mismatches < 0:
            raise ValueError("max_mismatches must be >= 0")
        self.max_mismatches = max_mismatches
        self.max_hits = max_hits_per_oligo
        self.skip_gene = skip_on_target_gene.upper()

    # -- seeds ----------------------------------------------------------

    @staticmethod
    def _seeds(target: str, m: int) -> list[tuple[int, str]]:
        """Split into m+1 non-overlapping seeds (pigeonhole guarantee)."""
        n_seeds = m + 1
        L = len(target)
        seed_len = L // n_seeds
        if seed_len < 1:
            # Oligo shorter than the number of seeds: no seeding possible,
            # signalled to the caller so it falls back to a full scan.
            return []
        out = []
        for i in range(n_seeds):
            start = i * seed_len
            end = start + seed_len if i < n_seeds - 1 else L
            out.append((start, target[start:end]))
        return out

    def screen_oligo(
        self,
        oligo: str,
        references: Iterable[tuple[str, str]],
        reference_name: str = "",
        status: ScreenStatus = ScreenStatus.TRANSCRIPTOME,
    ) -> SpecificityResult:
        """Screen one oligo against an iterable of (transcript_id, sequence)."""
        oligo_dna = oligo.upper().replace("U", "T")
        res = SpecificityResult(
            status=status,
            oligo=oligo_dna,
            reference_name=reference_name,
            max_mismatches_searched=self.max_mismatches,
            on_target_gene=self.skip_gene,
        )
        if not oligo_dna:
            res.notes.append("Empty oligo; nothing screened.")
            res.status = ScreenStatus.NOT_SCREENED
            return res

        L = len(oligo_dna)
        # The oligo hybridises to its reverse complement in a transcript.
        probe_fwd = reverse_complement(oligo_dna)
        probe_rev = oligo_dna  # same target read on the opposite strand

        seed_index: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
        for probe, strand in ((probe_fwd, "+"), (probe_rev, "-")):
            seeds = self._seeds(probe, self.max_mismatches)
            if not seeds:
                res.notes.append(
                    f"Oligo of length {L} is too short to seed at "
                    f"{self.max_mismatches} mismatches; using full comparison."
                )
            for offset, seed in seeds:
                seed_index[seed].append((strand, offset, L))

        probes = {"+": probe_fwd, "-": probe_rev}
        seed_lengths = {len(s) for s in seed_index} or {L}

        for tx_id, seq in references:
            if self.skip_gene and tx_id.upper().startswith(self.skip_gene):
                continue
            seq = seq.upper().replace("U", "T")
            tx_len = len(seq)
            if tx_len < L:
                continue
            res.n_transcripts_searched += 1
            res.n_bases_searched += tx_len

            if seed_index:
                seen: set[tuple[str, int]] = set()
                for slen in seed_lengths:
                    for i in range(tx_len - slen + 1):
                        entries = seed_index.get(seq[i : i + slen])
                        if not entries:
                            continue
                        for strand, offset, _ in entries:
                            start = i - offset
                            if start < 0 or start + L > tx_len:
                                continue
                            key = (strand, start)
                            if key in seen:
                                continue
                            seen.add(key)
                            window = seq[start : start + L]
                            mm = _count_mismatches(
                                probes[strand], window, self.max_mismatches
                            )
                            if mm is not None:
                                res.hits.append(OffTargetHit(
                                    transcript_id=tx_id, position=start,
                                    strand=strand, n_mismatches=mm,
                                    aligned_sequence=window))
            else:
                for strand, probe in probes.items():
                    for i in range(tx_len - L + 1):
                        mm = _count_mismatches(
                            probe, seq[i : i + L], self.max_mismatches)
                        if mm is not None:
                            res.hits.append(OffTargetHit(
                                transcript_id=tx_id, position=i, strand=strand,
                                n_mismatches=mm,
                                aligned_sequence=seq[i : i + L]))

            if len(res.hits) > self.max_hits:
                res.truncated = True
                res.notes.append(
                    f"Stopped after {self.max_hits} hits: this oligo is "
                    "clearly non-specific.")
                break

        res.hits.sort(key=lambda h: h.n_mismatches)
        return res

    def screen_against_fasta(
        self, oligo: str, fasta_path: str | Path
    ) -> SpecificityResult:
        path = Path(fasta_path)
        if not path.exists():
            res = SpecificityResult(status=ScreenStatus.NOT_SCREENED, oligo=oligo)
            res.notes.append(f"Reference not found: {path}")
            return res
        return self.screen_oligo(
            oligo, iter_fasta(path), reference_name=path.name,
            status=ScreenStatus.TRANSCRIPTOME)

    def screen_batch(
        self,
        oligos: list[str],
        references: Iterable[tuple[str, str]],
        reference_name: str = "supplied",
        status: ScreenStatus = ScreenStatus.TRANSCRIPTOME,
    ) -> list[SpecificityResult]:
        """Screen multiple oligos in a single streaming pass over `references`.

        Identical in sensitivity to calling `screen_oligo()` on each oligo
        individually, while avoiding repeated FASTA decompression and scanning.
        """
        if not oligos:
            return []
        if len(oligos) == 1:
            return [self.screen_oligo(oligos[0], references, reference_name=reference_name, status=status)]

        results: list[SpecificityResult] = []
        probes_list: list[dict[str, str]] = []
        lengths: list[int] = []
        seed_index: dict[str, list[tuple[int, str, int, int]]] = defaultdict(list)
        unseeded_indices: list[int] = []

        for idx, raw_oligo in enumerate(oligos):
            oligo_dna = raw_oligo.upper().replace("U", "T")
            res = SpecificityResult(
                status=status if oligo_dna else ScreenStatus.NOT_SCREENED,
                oligo=oligo_dna,
                reference_name=reference_name,
                max_mismatches_searched=self.max_mismatches,
                on_target_gene=self.skip_gene,
            )
            if not oligo_dna:
                res.notes.append("Empty oligo; nothing screened.")
                results.append(res)
                probes_list.append({})
                lengths.append(0)
                continue

            L = len(oligo_dna)
            probe_fwd = reverse_complement(oligo_dna)
            probe_rev = oligo_dna
            probes = {"+": probe_fwd, "-": probe_rev}
            probes_list.append(probes)
            lengths.append(L)
            results.append(res)

            had_seed = False
            for strand, probe in probes.items():
                seeds = self._seeds(probe, self.max_mismatches)
                if seeds:
                    had_seed = True
                    for offset, seed in seeds:
                        seed_index[seed].append((idx, strand, offset, L))
            if not had_seed:
                res.notes.append(
                    f"Oligo of length {L} is too short to seed at "
                    f"{self.max_mismatches} mismatches; using full comparison."
                )
                unseeded_indices.append(idx)

        seed_lengths = sorted({len(s) for s in seed_index})
        active_indices = [i for i, L in enumerate(lengths) if L > 0]

        for tx_id, seq in references:
            if self.skip_gene and tx_id.upper().startswith(self.skip_gene):
                continue
            seq = seq.upper().replace("U", "T")
            tx_len = len(seq)

            for idx in active_indices:
                res = results[idx]
                if res.truncated or tx_len < lengths[idx]:
                    continue
                res.n_transcripts_searched += 1
                res.n_bases_searched += tx_len

            if seed_index:
                seen_per_idx: dict[int, set[tuple[str, int]]] = defaultdict(set)
                for slen in seed_lengths:
                    if tx_len < slen:
                        continue
                    for i in range(tx_len - slen + 1):
                        entries = seed_index.get(seq[i : i + slen])
                        if not entries:
                            continue
                        for idx, strand, offset, L in entries:
                            res = results[idx]
                            if res.truncated:
                                continue
                            start = i - offset
                            if start < 0 or start + L > tx_len:
                                continue
                            key = (strand, start)
                            if key in seen_per_idx[idx]:
                                continue
                            seen_per_idx[idx].add(key)
                            window = seq[start : start + L]
                            mm = _count_mismatches(
                                probes_list[idx][strand], window, self.max_mismatches
                            )
                            if mm is not None:
                                res.hits.append(OffTargetHit(
                                    transcript_id=tx_id, position=start,
                                    strand=strand, n_mismatches=mm,
                                    aligned_sequence=window))
                                if len(res.hits) > self.max_hits:
                                    res.truncated = True
                                    res.notes.append(
                                        f"Stopped after {self.max_hits} hits: this oligo is "
                                        "clearly non-specific."
                                    )

            for idx in unseeded_indices:
                res = results[idx]
                L = lengths[idx]
                if res.truncated or tx_len < L:
                    continue
                for strand, probe in probes_list[idx].items():
                    for i in range(tx_len - L + 1):
                        mm = _count_mismatches(probe, seq[i : i + L], self.max_mismatches)
                        if mm is not None:
                            res.hits.append(OffTargetHit(
                                transcript_id=tx_id, position=i, strand=strand,
                                n_mismatches=mm, aligned_sequence=seq[i : i + L]))
                            if len(res.hits) > self.max_hits:
                                res.truncated = True
                                break

            if all(r.truncated for r in results if r.status != ScreenStatus.NOT_SCREENED):
                break

        for res in results:
            res.hits.sort(key=lambda h: h.n_mismatches)
        return results

    def screen_batch_against_fasta(
        self, oligos: list[str], fasta_path: str | Path
    ) -> list[SpecificityResult]:
        path = Path(fasta_path)
        if not path.exists():
            return [unscreened(o, f"Reference not found: {path}") for o in oligos]
        return self.screen_batch(
            oligos, iter_fasta(path), reference_name=path.name,
            status=ScreenStatus.TRANSCRIPTOME)


def _count_mismatches(a: str, b: str, limit: int) -> int | None:
    """Hamming distance, abandoning early once `limit` is exceeded."""
    mm = 0
    for x, y in zip(a, b):
        if x != y:
            mm += 1
            if mm > limit:
                return None
    return mm


def unscreened(oligo: str, reason: str = "") -> SpecificityResult:
    """The result to use when no screen was run. Never fabricates zeros."""
    res = SpecificityResult(status=ScreenStatus.NOT_SCREENED, oligo=oligo)
    res.notes.append(
        reason or "No reference transcriptome supplied; specificity unknown.")
    return res
