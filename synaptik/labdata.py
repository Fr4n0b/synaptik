"""
Bring-your-own-data: reconcile a lab's own results against the literature.

Everything else in SYNAPTIK answers questions about what is published. This
module answers the question a researcher actually has on the day the data
comes back: "which of MY results are new, which confirm what is known, which
contradict it -- and which known results did I fail to reproduce?"

Two inputs are supported:

  splicing tables   rMATS output (SE/RI/A5SS/A3SS/MXE .MATS.JC/JCEC.txt) or a
                    generic CSV/TSV with gene, delta_psi and (optionally) fdr
                    and event_type columns. MAJIQ, LeafCutter, SUPPA etc. can
                    be used through the generic format.

  phenotype tables  a CSV with one measured phenotype per row.

Design rules
------------
1. The user's data never leaves the machine and is never written to the
   knowledge base. Unpublished results are analysed in memory only.

2. The sign convention is never guessed. rMATS reports
   IncLevelDifference = IncLevel1 - IncLevel2, so whether a positive value
   means "more inclusion in the mutant" depends on which sample group was
   listed first. The caller must say which group is the perturbed condition.
   A sign error here inverts every conclusion, so there is no default.

3. "Tested and unchanged" is kept distinct from "not tested". A gene absent
   from the table says nothing; a gene present with a non-significant change
   is a (weak) negative result and is reported as such.

4. Literature matching is at gene + event-class level. Published papers rarely
   give coordinates, so a "consistent" call means "the literature reports the
   same class of event in the same gene in the same direction", not "the same
   exon". The output says so.
"""

from __future__ import annotations

import csv
import io
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

from scipy import stats

from .infer.splicing import (
    RBP_REFERENCE, NormalizedSplicingEvent, aso_strategy,
)
from .schema import Direction, SplicingEventType

RMATS_EVENT_TYPES = {
    "SE": SplicingEventType.SKIPPED_EXON,
    "RI": SplicingEventType.RETAINED_INTRON,
    "A5SS": SplicingEventType.ALT_5_SS,
    "A3SS": SplicingEventType.ALT_3_SS,
    "MXE": SplicingEventType.MUTUALLY_EXCLUSIVE,
}

REGULATOR_ALIASES = {"TDP-43": "TARDBP", "TDP43": "TARDBP", "A2BP1": "RBFOX1"}


class LabDataError(ValueError):
    """Raised when an input table cannot be interpreted safely."""


# --------------------------------------------------------------------------
# Parsed user events
# --------------------------------------------------------------------------


@dataclass
class UserSplicingEvent:
    gene: str
    event_type: SplicingEventType
    delta_psi: float | None        # positive = more inclusion in PERTURBED
    fdr: float | None = None
    p_value: float | None = None
    coordinates: str | None = None
    exon_length: int | None = None   # regulated exon, when the table gives it
    row_number: int = 0

    @property
    def direction(self) -> Direction:
        if self.delta_psi is None:
            return Direction.UNSPECIFIED
        return Direction.INCREASE if self.delta_psi > 0 else Direction.DECREASE


@dataclass
class SplicingTable:
    events: list[UserSplicingEvent]
    source_format: str
    event_type_source: str
    n_rows_read: int = 0
    n_rows_skipped: int = 0
    skip_reasons: dict[str, int] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    def genes_tested(self) -> set[str]:
        return {e.gene.upper() for e in self.events}

    def significant(
        self, max_fdr: float = 0.05, min_abs_dpsi: float = 0.1
    ) -> list[UserSplicingEvent]:
        """Events passing both thresholds.

        If the table has no FDR column at all, significance cannot be assessed
        from the table and only the effect-size filter applies; the caller is
        warned via `warnings` at parse time.
        """
        out = []
        for e in self.events:
            if e.delta_psi is None or abs(e.delta_psi) < min_abs_dpsi:
                continue
            if e.fdr is not None and e.fdr > max_fdr:
                continue
            out.append(e)
        return out


def _sniff_delimiter(text: str) -> str:
    head = text[:4096]
    return "\t" if head.count("\t") >= head.count(",") else ","


def _to_float(v: str | None) -> float | None:
    if v is None:
        return None
    s = str(v).strip()
    if s == "" or s.upper() in ("NA", "NAN", "NULL", "NONE", "N/A", "INF", "-INF"):
        return None
    if ";" in s:
        # Multiple values (e.g. MAJIQ multiple junctions per LSV) -> select max magnitude
        sub_vals = [_to_float(x) for x in s.split(";")]
        valid = [x for x in sub_vals if x is not None]
        return max(valid, key=abs) if valid else None
    try:
        f = float(s)
    except ValueError:
        return None
    return None if math.isnan(f) or math.isinf(f) else f


def _clean_gene(v: str | None) -> str:
    g = (v or "").strip().strip('"').strip("'")
    if "," in g:
        g = g.split(",")[0].strip()
    return g


def parse_splicing_table(
    source: str | Path,
    perturbed: str,
    event_type: str | None = None,
    text: str | None = None,
) -> SplicingTable:
    """Parse an rMATS table or a generic splicing CSV/TSV.

    perturbed: which side of the comparison is the perturbed/mutant condition.
        "sample1" or "sample2" for rMATS (IncLevelDifference = S1 - S2);
        for generic tables, "positive" if a positive delta_psi already means
        more inclusion in the perturbed condition, "negative" if reversed.
    event_type: SE/RI/A5SS/A3SS/MXE. Inferred from an rMATS file name or an
        event_type column when not given.
    text: file contents, when the data did not come from disk (web upload).
    """
    perturbed = (perturbed or "").strip().lower()
    name = Path(str(source)).name
    if text is None:
        try:
            text = Path(source).read_text(encoding="utf-8-sig", errors="replace")
        except OSError as e:
            raise LabDataError(f"Cannot read {source}: {e}") from e
    if not text.strip():
        raise LabDataError("The file is empty.")

    delim = _sniff_delimiter(text)
    reader = csv.DictReader(io.StringIO(text), delimiter=delim)
    if not reader.fieldnames:
        raise LabDataError("No header row found.")
    cols = {c.strip().strip('"').lower(): c for c in reader.fieldnames if c}

    is_rmats = "inclevedifference" in cols or "incleveldifference" in cols
    is_majiq = any("lsv" in c for c in cols) or any("mean_dpsi" in c for c in cols) or any("e(dpsi)" in c for c in cols)
    is_leafcutter = "cluster" in cols or any("clu_" in c for c in cols) or ("logef" in cols and any(k in cols for k in ("p.adjust", "padj", "genes", "gene")))

    if is_rmats:
        if perturbed not in ("sample1", "sample2"):
            raise LabDataError(
                "This is an rMATS table. rMATS reports IncLevelDifference as "
                "sample 1 minus sample 2, so you must say which group is the "
                "perturbed condition: perturbed='sample1' or 'sample2' "
                "(the order you gave --b1 and --b2). There is no safe default."
            )
        sign = 1.0 if perturbed == "sample1" else -1.0
        dpsi_col = cols.get("incleveldifference") or cols["inclevedifference"]
        gene_col = cols.get("genesymbol") or cols.get("geneid")
        fdr_col, p_col = cols.get("fdr"), cols.get("pvalue")
        fmt = "rMATS"
    elif is_majiq:
        if perturbed not in ("sample1", "sample2", "positive", "negative", "cond1", "cond2"):
            raise LabDataError(
                "This is a MAJIQ (VOILA) table. MAJIQ reports E(dPSI) as Condition 1 minus "
                "Condition 2. State which condition is the perturbed condition: "
                "perturbed='sample1' (or 'cond1') if condition 1 is perturbed, "
                "or 'sample2' (or 'cond2') if condition 2 is perturbed."
            )
        sign = -1.0 if perturbed in ("sample2", "cond2", "negative") else 1.0
        dpsi_col = next((cols[c] for c in (
            "e(dpsi) per lsv junction", "mean_dpsi_per_lsv_junction",
            "e(dpsi)", "mean_dpsi", "dpsi", "delta_psi"
        ) if c in cols), None)
        gene_col = next((cols[c] for c in (
            "gene name", "gene_name", "genesymbol", "gene_symbol", "gene id", "gene_id", "gene"
        ) if c in cols), None)
        fdr_col = next((cols[c] for c in (
            "p(|dpsi| > 0.20) per lsv junction", "p(|dpsi| > 0.10) per lsv junction",
            "p(|dpsi|>0.20) per lsv junction", "p(|dpsi|>0.10) per lsv junction",
            "probability_changing", "fdr", "padj", "qvalue"
        ) if c in cols), None)
        p_col = next((cols[c] for c in ("pvalue", "p_value", "p") if c in cols), None)
        fmt = "MAJIQ (VOILA)"
    elif is_leafcutter:
        if perturbed not in ("positive", "negative", "sample1", "sample2"):
            raise LabDataError(
                "This is a LeafCutter table. State the sign convention: "
                "perturbed='positive' if a positive delta_psi or log effect means more "
                "inclusion in the perturbed condition, 'negative' if it means less."
            )
        sign = -1.0 if perturbed in ("negative", "sample2") else 1.0
        dpsi_col = next((cols[c] for c in (
            "deltapsi", "delta_psi", "dpsi", "logef", "effect_size"
        ) if c in cols), None)
        gene_col = next((cols[c] for c in (
            "genes", "gene", "gene_name", "gene_symbol", "symbol"
        ) if c in cols), None)
        fdr_col = next((cols[c] for c in (
            "p.adjust", "padj", "fdr", "qvalue", "q_value"
        ) if c in cols), None)
        p_col = next((cols[c] for c in ("p", "pvalue", "p_value") if c in cols), None)
        fmt = "LeafCutter"
    else:
        if perturbed not in ("positive", "negative"):
            raise LabDataError(
                "For a generic splicing table, state the sign convention: "
                "perturbed='positive' if a positive delta_psi means more "
                "inclusion in the perturbed condition, 'negative' if it means "
                "less."
            )
        sign = 1.0 if perturbed == "positive" else -1.0
        dpsi_col = next(
            (cols[c] for c in ("delta_psi", "dpsi", "deltapsi", "delta psi",
                               "e(dpsi)", "mean_dpsi") if c in cols), None)
        gene_col = next(
            (cols[c] for c in ("gene", "gene_symbol", "genesymbol",
                               "gene_name", "symbol") if c in cols), None)
        fdr_col = next(
            (cols[c] for c in ("fdr", "padj", "qvalue", "q_value",
                               "adj_pvalue", "p_adj") if c in cols), None)
        p_col = next((cols[c] for c in ("pvalue", "p_value", "p")
                      if c in cols), None)
        fmt = "generic"

    if not dpsi_col or not gene_col:
        raise LabDataError(
            "Could not find the required columns. Need a gene column "
            "(gene / gene_symbol / geneSymbol) and a delta-PSI column "
            "(delta_psi / dpsi / IncLevelDifference). Found: "
            + ", ".join(reader.fieldnames[:20])
        )

    # Event type: explicit argument > per-row column > rMATS file name.
    et_col = cols.get("event_type") or cols.get("type")
    file_et = None
    m = re.match(r"^(SE|RI|A5SS|A3SS|MXE)\.MATS", name, re.IGNORECASE)
    if m:
        file_et = m.group(1).upper()
    fixed_et = None
    et_source = "unspecified"
    if event_type:
        key = event_type.strip().upper()
        if key not in RMATS_EVENT_TYPES:
            raise LabDataError(
                f"Unknown event type '{event_type}'. Use one of: "
                + ", ".join(RMATS_EVENT_TYPES))
        fixed_et, et_source = RMATS_EVENT_TYPES[key], "argument"
    elif et_col:
        et_source = "column"
    elif file_et:
        fixed_et, et_source = RMATS_EVENT_TYPES[file_et], "file name"

    # Coordinates: take whatever start/end pair this event type uses.
    chr_col = cols.get("chr") or cols.get("chrom") or cols.get("seqid")
    start_col = next((cols[c] for c in cols
                      if c.endswith("start_0base") or c in ("start",)), None)
    end_col = next((cols[c] for c in (
        "exonend", "riexonend", "longexonend", "1stexonend", "end")
        if c in cols), None)
    strand_col = cols.get("strand")

    table = SplicingTable(events=[], source_format=fmt,
                          event_type_source=et_source)
    skips: Counter = Counter()

    for i, row in enumerate(reader, start=2):
        table.n_rows_read += 1
        gene = _clean_gene(row.get(gene_col))
        if not gene or gene.upper() in ("NA", "."):
            skips["no gene symbol"] += 1
            continue
        raw = _to_float(row.get(dpsi_col))
        if raw is None:
            skips["delta-PSI missing or not numeric"] += 1
            continue
        if abs(raw) > 1.0:
            # PSI differences live in [-1, 1]. A value outside means the
            # column is on a percent scale or is not a dPSI at all; rescaling
            # silently would be a guess.
            skips["delta-PSI outside [-1, 1]"] += 1
            continue

        et = fixed_et
        if et is None and et_col:
            et = RMATS_EVENT_TYPES.get(
                (row.get(et_col) or "").strip().upper(),
                SplicingEventType.UNSPECIFIED)
        if et is None:
            et = SplicingEventType.UNSPECIFIED

        coords = None
        c = (row.get(chr_col) or "").strip() if chr_col else ""
        s = (row.get(start_col) or "").strip() if start_col else ""
        e = (row.get(end_col) or "").strip() if end_col else ""
        strand = (row.get(strand_col) or "").strip() if strand_col else ""

        if not (c and s and e):
            # Fallback for LeafCutter / MAJIQ cluster or junction coordinates
            clu_col = next((cols[k] for k in ("cluster", "intron", "junctions coords", "junction") if k in cols), None)
            if clu_col and row.get(clu_col):
                raw_clu = str(row[clu_col]).strip()
                first_seg = raw_clu.split(";")[0].strip()
                parts = first_seg.split(":")
                if len(parts) >= 3:
                    # LeafCutter: chr:start:end:... or chr:start:end
                    c = c or parts[0]
                    s = s or parts[1]
                    e = e or parts[2]
                    if len(parts) >= 5 and parts[4] in ("+", "-"):
                        strand = strand or parts[4]
                    elif parts[-1] in ("+", "-"):
                        strand = strand or parts[-1]
                elif len(parts) == 2 and "-" in parts[1]:
                    # MAJIQ: chr:start-end
                    c = c or parts[0]
                    se_parts = parts[1].split("-")
                    if len(se_parts) == 2:
                        s = s or se_parts[0]
                        e = e or se_parts[1]

        if c and s and e:
            coords = f"{c}:{s}-{e}" + (f":{strand}" if strand else "")

        exon_length = None
        if s and e:
            es, ee = _to_float(s), _to_float(e)
            if es is not None and ee is not None and ee > es:
                exon_length = int(ee - es)

        table.events.append(UserSplicingEvent(
            gene=gene, event_type=et, delta_psi=sign * raw,
            exon_length=exon_length,
            fdr=_to_float(row.get(fdr_col)) if fdr_col else None,
            p_value=_to_float(row.get(p_col)) if p_col else None,
            coordinates=coords, row_number=i,
        ))

    table.n_rows_skipped = sum(skips.values())
    table.skip_reasons = dict(skips)

    if not table.events:
        raise LabDataError(
            f"No usable rows ({table.n_rows_read} read). "
            + "; ".join(f"{k}: {v}" for k, v in skips.items()))
    if not fdr_col:
        table.warnings.append(
            "No FDR column found, so statistical significance cannot be "
            "assessed from this table; only the effect-size threshold is "
            "applied. Treat every call below as provisional.")
    if skips.get("delta-PSI outside [-1, 1]"):
        table.warnings.append(
            f"{skips['delta-PSI outside [-1, 1]']} row(s) had delta-PSI "
            "outside [-1, 1] and were skipped. If your table is on a 0-100 "
            "scale, divide by 100 first.")
    if et_source == "unspecified":
        table.warnings.append(
            "Event type unknown (not in the file name, no event_type column, "
            "none given). Literature matching falls back to gene level only.")
    if table.n_rows_skipped > 0.2 * max(1, table.n_rows_read):
        table.warnings.append(
            f"{table.n_rows_skipped} of {table.n_rows_read} rows were "
            "skipped. Check that the columns are what this parser assumes.")
    return table


# --------------------------------------------------------------------------
# Reconciliation against the literature corpus
# --------------------------------------------------------------------------


def frame_consequence(ue: UserSplicingEvent) -> tuple[str, str]:
    """Reading-frame consequence of a user's event, from its own coordinates.

    The generic literature heuristic ("exon inclusion goes up, so suspect a
    cryptic exon and NMD") over-calls on a lab's table, where most cassette
    exons are annotated and in-frame. When the table gives the exon's
    coordinates the question can be answered directly: an exon whose length is
    not a multiple of three shifts the frame when it is skipped or included.
    """
    if ue.event_type == SplicingEventType.RETAINED_INTRON:
        if ue.direction == Direction.INCREASE:
            return ("high", "Retained introns usually carry a premature stop "
                            "codon; the transcript is likely degraded by NMD "
                            "or held in the nucleus. Measure protein.")
        return ("low", "Less intron retention: more productive transcript.")
    if ue.event_type != SplicingEventType.SKIPPED_EXON:
        return ("unknown", "Frame impact cannot be judged for this event "
                           "class from the table alone.")
    if ue.exon_length is None:
        return ("unknown", "No exon coordinates in the table, so the frame "
                           "consequence is unknown. Check whether the exon "
                           "length is a multiple of three.")
    if ue.exon_length % 3 != 0:
        return ("high",
                f"The exon is {ue.exon_length} nt (not a multiple of three): "
                "if it is coding, changing its inclusion shifts the reading "
                "frame and one isoform is a likely NMD target. Measure "
                "protein, not only transcript.")
    return ("low",
            f"The exon is {ue.exon_length} nt (frame-preserving). An "
            "unannotated exon can still contain an in-frame stop codon.")


class LabCall(str, Enum):
    NOVEL = "novel"                    # no splicing event reported in this gene
    CONSISTENT = "consistent"          # literature reports same class + direction
    DISCORDANT = "discordant"          # literature reports the opposite direction
    SAME_GENE = "same_gene_other_event"  # gene is known, event class differs


@dataclass
class ReconciledEvent:
    event: UserSplicingEvent
    call: LabCall
    n_papers: int = 0
    n_labs: int = 0
    literature_directions: dict[str, int] = field(default_factory=dict)
    literature_event_types: list[str] = field(default_factory=list)
    regulators_reported: list[str] = field(default_factory=list)
    source_ids: list[str] = field(default_factory=list)
    nmd_risk: str = "unknown"
    nmd_note: str = ""
    aso_feasible: bool = False
    aso_approach: str = ""
    note: str = ""
    priority: float = 0.0


@dataclass
class NotReproduced:
    """A published event in a gene the user tested, with no change observed."""

    gene: str
    n_papers: int
    n_labs: int
    literature_direction: str
    regulators: list[str]
    user_best_dpsi: float | None
    user_best_fdr: float | None
    source_ids: list[str]


@dataclass
class RegulatorEnrichment:
    """Do the user's changed genes overlap an RBP's published targets?"""

    rbp: str
    n_targets_in_corpus: int
    n_targets_tested: int
    n_targets_changed: int
    overlap_genes: list[str]
    odds_ratio: float
    p_value: float
    q_value: float = 1.0
    direction_agreement: float | None = None
    n_supporting_labs: int = 0
    canonical_effect: str | None = None

    def caveats(self) -> list[str]:
        out = []
        if self.n_targets_tested < 5:
            out.append(
                f"Only {self.n_targets_tested} of this regulator's published "
                "targets were tested in your data; the estimate is unstable.")
        if self.n_supporting_labs < 2:
            out.append(
                "The published target list for this regulator comes from a "
                "single research group in the corpus.")
        if self.direction_agreement is not None and self.direction_agreement < 0.6:
            out.append(
                "Your events overlap the published targets but mostly change "
                "in the opposite direction -- consistent with the regulator "
                "being involved, but not with the same change in its activity.")
        return out


@dataclass
class SplicingDepositReport:
    table: SplicingTable
    n_events: int = 0
    n_genes_tested: int = 0
    n_significant: int = 0
    thresholds: dict = field(default_factory=dict)
    reconciled: list[ReconciledEvent] = field(default_factory=list)
    not_reproduced: list[NotReproduced] = field(default_factory=list)
    regulators: list[RegulatorEnrichment] = field(default_factory=list)
    corpus_events: int = 0
    corpus_genes: int = 0
    warnings: list[str] = field(default_factory=list)
    summary: str = ""

    def by_call(self, call: LabCall) -> list[ReconciledEvent]:
        return [r for r in self.reconciled if r.call == call]

    def counts(self) -> dict[str, int]:
        c = Counter(r.call.value for r in self.reconciled)
        return {k.value: c.get(k.value, 0) for k in LabCall}


def _canon_regulator(name: str | None) -> str | None:
    if not name:
        return None
    key = name.strip().upper()
    key = REGULATOR_ALIASES.get(key, key)
    key = re.sub(r"[^A-Z0-9]", "", key)
    return key or None


def _bh(pvals: list[float]) -> list[float]:
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [1.0] * m
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        prev = min(prev, pvals[i] * m / (rank + 1))
        q[i] = min(1.0, prev)
    return q


class LabDataReconciler:
    """Reconciles a lab's splicing results against published events."""

    def __init__(self, max_fdr: float = 0.05, min_abs_dpsi: float = 0.1):
        self.max_fdr = max_fdr
        self.min_abs_dpsi = min_abs_dpsi

    def reconcile_splicing(
        self,
        table: SplicingTable,
        corpus_events: list[NormalizedSplicingEvent],
        lab_ids: dict[str, str] | None = None,
    ) -> SplicingDepositReport:
        lab_ids = lab_ids or {}
        rep = SplicingDepositReport(
            table=table,
            n_events=len(table.events),
            n_genes_tested=len(table.genes_tested()),
            thresholds={"max_fdr": self.max_fdr,
                        "min_abs_dpsi": self.min_abs_dpsi},
            warnings=list(table.warnings),
        )

        by_gene: dict[str, list[NormalizedSplicingEvent]] = defaultdict(list)
        for ce in corpus_events:
            if ce.gene and ce.gene.lower() != "unknown":
                by_gene[ce.gene.upper()].append(ce)
        rep.corpus_events = len(corpus_events)
        rep.corpus_genes = len(by_gene)
        if not corpus_events:
            rep.warnings.append(
                "The knowledge base contains no published splicing events, so "
                "every result will be called 'novel'. That reflects the corpus, "
                "not the literature -- harvest and extract splicing papers "
                "first.")

        sig = table.significant(self.max_fdr, self.min_abs_dpsi)
        rep.n_significant = len(sig)

        for ue in sig:
            rep.reconciled.append(
                self._reconcile_one(ue, by_gene.get(ue.gene.upper(), []), lab_ids))
        rep.reconciled.sort(key=lambda r: -r.priority)

        rep.not_reproduced = self._not_reproduced(table, sig, by_gene, lab_ids)
        rep.regulators = self._regulator_enrichment(table, sig, by_gene, lab_ids)
        rep.summary = self._summarize(rep)
        return rep

    # ---------------- per-event ----------------

    def _reconcile_one(
        self, ue: UserSplicingEvent, lit: list[NormalizedSplicingEvent],
        lab_ids: dict[str, str],
    ) -> ReconciledEvent:
        level, note = frame_consequence(ue)
        strat = aso_strategy(ue.event_type, ue.direction)
        r = ReconciledEvent(
            event=ue, call=LabCall.NOVEL, nmd_risk=level, nmd_note=note,
            aso_feasible=strat.feasible, aso_approach=strat.approach,
        )
        magnitude = abs(ue.delta_psi or 0.0)

        if not lit:
            r.note = (
                "No published splicing event in this gene in the corpus. "
                "Check PubMed before claiming novelty: the corpus covers "
                "open-access full text and abstracts, not everything.")
            r.priority = 3.0 * magnitude
            return r

        papers = {e.source_id for e in lit}
        r.n_papers = len(papers)
        r.n_labs = len({lab_ids.get(s) for s in papers if lab_ids.get(s)}) \
            or len(papers)
        r.source_ids = sorted(papers)[:8]
        r.literature_event_types = sorted({
            e.event_type.value for e in lit
            if e.event_type != SplicingEventType.UNSPECIFIED})
        r.regulators_reported = sorted({
            c for c in (_canon_regulator(e.regulator) for e in lit) if c})

        typed = ue.event_type != SplicingEventType.UNSPECIFIED
        same_class = [e for e in lit if e.event_type == ue.event_type] \
            if typed else list(lit)
        # Literature events with no stated class cannot contradict the class.
        if typed and not same_class:
            same_class = [e for e in lit
                          if e.event_type == SplicingEventType.UNSPECIFIED]

        if not same_class:
            r.call = LabCall.SAME_GENE
            r.note = (
                f"The literature reports {', '.join(r.literature_event_types)} "
                f"event(s) in {ue.gene}, not {ue.event_type.value}. Your event "
                "may be a different one in the same gene.")
            r.priority = 1.8 * magnitude
            return r

        dirs = Counter(
            e.direction.value for e in same_class
            if e.direction in (Direction.INCREASE, Direction.DECREASE))
        r.literature_directions = dict(dirs)
        if not dirs:
            r.call = LabCall.CONSISTENT
            r.note = ("The literature reports this class of event in this gene "
                      "but does not state its direction.")
            r.priority = 1.0 * magnitude
            return r

        mine = ue.direction.value
        agree = dirs.get(mine, 0)
        against = sum(v for k, v in dirs.items() if k != mine)
        uncertain = sum(1 for e in same_class if e.sign_confidence != "high")

        if agree >= against:
            r.call = LabCall.CONSISTENT
            r.note = (
                f"{agree} published report(s) from {r.n_labs} group(s) agree "
                "on the direction"
                + (f"; {against} report the opposite." if against else "."))
            r.priority = 1.0 * magnitude * (1.0 + 0.2 * min(3, r.n_labs))
        else:
            r.call = LabCall.DISCORDANT
            r.note = (
                f"{against} published report(s) show the OPPOSITE direction "
                f"({agree} agree with you). Before interpreting: confirm which "
                "rMATS group was the perturbed condition, the genome build, "
                "and whether the published event is the same exon.")
            if uncertain:
                r.note += (
                    f" {uncertain} of the published events have an uncertain "
                    "dPSI sign, so the conflict may be an artefact of their "
                    "reporting.")
            r.priority = 2.5 * max(magnitude, 0.2)
        return r

    # ---------------- expected but absent ----------------

    def _not_reproduced(
        self, table: SplicingTable, sig: list[UserSplicingEvent],
        by_gene: dict[str, list[NormalizedSplicingEvent]],
        lab_ids: dict[str, str],
    ) -> list[NotReproduced]:
        """Published events in genes the user tested but saw no change in.

        Only genes PRESENT in the table qualify. Absence from the table means
        the gene was not expressed or not testable, which is not a negative
        result.
        """
        changed = {e.gene.upper() for e in sig}
        tested: dict[str, list[UserSplicingEvent]] = defaultdict(list)
        for e in table.events:
            tested[e.gene.upper()].append(e)

        out = []
        for gene, lit in by_gene.items():
            if gene not in tested or gene in changed:
                continue
            papers = {e.source_id for e in lit}
            n_labs = len({lab_ids.get(s) for s in papers if lab_ids.get(s)}) \
                or len(papers)
            dirs = Counter(e.direction.value for e in lit
                           if e.direction != Direction.UNSPECIFIED)
            best = max(tested[gene], key=lambda e: abs(e.delta_psi or 0.0))
            out.append(NotReproduced(
                gene=tested[gene][0].gene,
                n_papers=len(papers), n_labs=n_labs,
                literature_direction=(
                    dirs.most_common(1)[0][0] if dirs else "unspecified"),
                regulators=sorted({
                    c for c in (_canon_regulator(e.regulator) for e in lit) if c}),
                user_best_dpsi=best.delta_psi, user_best_fdr=best.fdr,
                source_ids=sorted(papers)[:6],
            ))
        out.sort(key=lambda n: (-n.n_labs, -n.n_papers))
        return out

    # ---------------- regulator enrichment ----------------

    def _regulator_enrichment(
        self, table: SplicingTable, sig: list[UserSplicingEvent],
        by_gene: dict[str, list[NormalizedSplicingEvent]],
        lab_ids: dict[str, str],
    ) -> list[RegulatorEnrichment]:
        """Fisher test: are published targets of an RBP over-represented
        among the user's changed genes?

        The background is the set of genes TESTED in the user's table, not the
        genome. Using the genome would make any expressed-gene list look
        enriched for any neuronal regulator.
        """
        tested = table.genes_tested()
        changed = {e.gene.upper() for e in sig}
        if not tested or not changed:
            return []
        my_dir: dict[str, Direction] = {}
        for e in sorted(sig, key=lambda e: -abs(e.delta_psi or 0.0)):
            my_dir.setdefault(e.gene.upper(), e.direction)

        targets: dict[str, dict[str, list[NormalizedSplicingEvent]]] = \
            defaultdict(lambda: defaultdict(list))
        for gene, lit in by_gene.items():
            for e in lit:
                reg = _canon_regulator(e.regulator)
                if reg:
                    targets[reg][gene].append(e)

        results = []
        for rbp, tgt in targets.items():
            t_tested = set(tgt) & tested
            if len(t_tested) < 2:
                continue
            a = len(t_tested & changed)
            b = len(t_tested) - a
            c = len(changed - t_tested)
            d = len(tested) - a - b - c
            if a == 0:
                continue
            _, p = stats.fisher_exact([[a, b], [c, max(0, d)]],
                                      alternative="greater")
            odds = ((a + 0.5) * (max(0, d) + 0.5)) / ((b + 0.5) * (c + 0.5))

            agree = total = 0
            for g in t_tested & changed:
                lit_dirs = Counter(
                    e.direction for e in tgt[g]
                    if e.direction in (Direction.INCREASE, Direction.DECREASE))
                if not lit_dirs:
                    continue
                total += 1
                if lit_dirs.most_common(1)[0][0] == my_dir.get(g):
                    agree += 1

            papers = {e.source_id for evs in tgt.values() for e in evs}
            labs = {lab_ids.get(s) for s in papers if lab_ids.get(s)}
            profile = RBP_REFERENCE.get(rbp)
            results.append(RegulatorEnrichment(
                rbp=rbp, n_targets_in_corpus=len(tgt),
                n_targets_tested=len(t_tested), n_targets_changed=a,
                overlap_genes=sorted(t_tested & changed)[:15],
                odds_ratio=round(odds, 2), p_value=float(p),
                direction_agreement=(agree / total) if total else None,
                n_supporting_labs=len(labs) or len(papers),
                canonical_effect=profile.canonical_effect if profile else None,
            ))

        if results:
            for r, q in zip(results, _bh([r.p_value for r in results])):
                r.q_value = q
            results.sort(key=lambda r: (r.q_value, -r.n_targets_changed))
        return results

    # ---------------- summary ----------------

    @staticmethod
    def _summarize(rep: SplicingDepositReport) -> str:
        c = rep.counts()
        parts = [
            f"{rep.n_significant} of {rep.n_events} events in "
            f"{rep.n_genes_tested} genes pass the thresholds."
        ]
        if rep.n_significant:
            parts.append(
                f"{c['novel']} have no published counterpart in the corpus, "
                f"{c['consistent']} agree with published reports, "
                f"{c['discordant']} contradict them, and "
                f"{c['same_gene_other_event']} fall in genes with a different "
                "published event.")
        if rep.not_reproduced:
            parts.append(
                f"{len(rep.not_reproduced)} published event(s) fall in genes "
                "you tested but saw no change in.")
        top = [r for r in rep.regulators if r.q_value < 0.1]
        if top:
            parts.append(
                "Your changed genes are enriched for published targets of "
                + ", ".join(r.rbp for r in top[:3]) + ".")
        return " ".join(parts)


# --------------------------------------------------------------------------
# Phenotype tables
# --------------------------------------------------------------------------


@dataclass
class UserPhenotype:
    phenotype: str
    direction: Direction
    effect_size: float | None = None
    gene: str | None = None
    assay: str | None = None
    days_in_vitro: int | None = None
    row_number: int = 0


_DIRECTION_WORDS = {
    "increase": Direction.INCREASE, "increased": Direction.INCREASE,
    "up": Direction.INCREASE, "higher": Direction.INCREASE,
    "decrease": Direction.DECREASE, "decreased": Direction.DECREASE,
    "reduced": Direction.DECREASE, "down": Direction.DECREASE,
    "lower": Direction.DECREASE,
    "no_change": Direction.NO_CHANGE, "no change": Direction.NO_CHANGE,
    "unchanged": Direction.NO_CHANGE, "ns": Direction.NO_CHANGE,
    "none": Direction.NO_CHANGE,
}


def parse_phenotype_table(
    source: str | Path, text: str | None = None
) -> tuple[list[UserPhenotype], list[str]]:
    """Parse a phenotype CSV/TSV.

    Required columns: phenotype, direction (increase / decrease / no_change).
    Optional: effect_size, gene, assay, div.
    """
    if text is None:
        try:
            text = Path(source).read_text(encoding="utf-8-sig", errors="replace")
        except OSError as e:
            raise LabDataError(f"Cannot read {source}: {e}") from e
    reader = csv.DictReader(io.StringIO(text), delimiter=_sniff_delimiter(text))
    if not reader.fieldnames:
        raise LabDataError("No header row found.")
    cols = {c.strip().lower(): c for c in reader.fieldnames if c}
    ph_col = cols.get("phenotype") or cols.get("readout")
    dir_col = cols.get("direction") or cols.get("change")
    if not ph_col or not dir_col:
        raise LabDataError(
            "A phenotype table needs 'phenotype' and 'direction' columns "
            "(direction: increase / decrease / no_change). Found: "
            + ", ".join(reader.fieldnames[:20]))

    rows, warnings = [], []
    for i, row in enumerate(reader, start=2):
        name = (row.get(ph_col) or "").strip()
        if not name:
            continue
        word = (row.get(dir_col) or "").strip().lower()
        direction = _DIRECTION_WORDS.get(word)
        if direction is None:
            # An unreadable direction is not "no change"; skip and say so.
            warnings.append(
                f"Row {i}: direction '{row.get(dir_col)}' not understood; "
                "row skipped.")
            continue
        div = _to_float(row.get(cols.get("div", ""))) if "div" in cols else None
        rows.append(UserPhenotype(
            phenotype=name, direction=direction,
            effect_size=_to_float(row.get(cols["effect_size"]))
            if "effect_size" in cols else None,
            gene=(row.get(cols["gene"]) or "").strip() or None
            if "gene" in cols else None,
            assay=(row.get(cols["assay"]) or "").strip() or None
            if "assay" in cols else None,
            days_in_vitro=int(div) if div is not None else None,
            row_number=i,
        ))
    if not rows:
        raise LabDataError("No usable rows in the phenotype table.")
    return rows, warnings


@dataclass
class ReconciledPhenotype:
    observation: UserPhenotype
    resolved_label: str | None
    call: str                  # novel / consistent / discordant / contested
    n_records: int = 0
    n_labs: int = 0
    literature_directions: dict[str, int] = field(default_factory=dict)
    note: str = ""
    source_ids: list[str] = field(default_factory=list)


def reconcile_phenotypes(
    observations: list[UserPhenotype], store, matcher=None
) -> list[ReconciledPhenotype]:
    """Compare each observed phenotype with what the corpus reports.

    When a gene is given, only records perturbing or mentioning that gene are
    used -- "reduced puncta" is only informative relative to the same
    perturbation.
    """
    from .normalize.matcher import PhenotypeMatcher

    matcher = matcher or PhenotypeMatcher()
    out = []
    for ob in observations:
        term = matcher.match(ob.phenotype)
        keys = [term.key()]
        recs = store.records_for_phenotype(keys)
        if ob.gene:
            g = ob.gene.upper()
            recs = [
                r for r in recs
                if (r.perturbation.target and
                    r.perturbation.target.raw.upper() == g)
                or g in {x.raw.upper() for x in r.genes}
            ]
        rp = ReconciledPhenotype(
            observation=ob, resolved_label=term.ontology_label, call="novel",
            n_records=len(recs),
            n_labs=len({r.source.lab_id for r in recs if r.source.lab_id}),
            source_ids=sorted({r.source.source_id for r in recs})[:8],
        )
        if not recs:
            rp.note = (
                "No matching record in the corpus"
                + (f" for {ob.gene}" if ob.gene else "")
                + ". Either this has not been reported, or the phenotype "
                "wording is not in the vocabulary.")
            out.append(rp)
            continue

        dirs = Counter(r.outcome.direction.value for r in recs
                       if r.outcome.direction != Direction.UNSPECIFIED)
        rp.literature_directions = dict(dirs)
        mine = ob.direction.value
        total = sum(dirs.values())
        agree = dirs.get(mine, 0)
        if total == 0:
            rp.call, rp.note = "consistent", (
                "Reported in the literature without a stated direction.")
        elif agree / total >= 0.7:
            rp.call = "consistent"
            rp.note = (f"{agree} of {total} published records "
                       f"({rp.n_labs} group(s)) report the same direction.")
        elif agree / total <= 0.3:
            rp.call = "discordant"
            rp.note = (
                f"Only {agree} of {total} published records agree with your "
                "direction. Compare maturity (DIV), cell model and controls "
                "with the published work before interpreting.")
        else:
            rp.call = "contested"
            rp.note = (
                f"The literature itself is split ({dict(dirs)}); your result "
                "adds a data point to an open question.")
        out.append(rp)
    return out
