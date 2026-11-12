"""
Surface-form -> ontology matching.

Strategy is deliberately layered, cheapest first:
  1. exact match on normalized string
  2. synonym-index match on longest matching n-gram
  3. fuzzy match (token-set ratio) above a threshold
  4. unmapped -- retained as raw, never silently dropped

Layer 4 is the important one. Tools that force every string into an ontology
produce confident garbage. We keep unmapped terms visible so they can be
reviewed and folded into the vocabulary, and we report the unmapped rate as a
first-class quality metric.
"""

from __future__ import annotations

import re
import unicodedata
from difflib import SequenceMatcher
from functools import lru_cache

from ..schema import AssayClass, NormalizedTerm
from .phenotype_ontology import SYNONYM_INDEX, BY_ID

# Greek letters and typographic variants that appear constantly in this domain
_GREEK = {
    "α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta", "κ": "kappa",
    "λ": "lambda", "μ": "mu", "σ": "sigma", "τ": "tau", "ω": "omega",
}


def normalize_string(s: str) -> str:
    """Aggressive normalization for matching only. Never stored as the raw."""
    if not s:
        return ""
    s = unicodedata.normalize("NFKD", s)
    for g, name in _GREEK.items():
        s = s.replace(g, name)
    s = s.lower().strip()
    s = re.sub(r"[\u2010-\u2015]", "-", s)  # unicode dashes
    s = re.sub(r"['\u2018\u2019\u201c\u201d]", "", s)
    s = re.sub(r"\s+", " ", s)
    s = re.sub(r"[.,;:]+$", "", s)
    return s.strip()


def _token_set_ratio(a: str, b: str) -> float:
    ta, tb = set(a.split()), set(b.split())
    if not ta or not tb:
        return 0.0
    inter = ta & tb
    # Order-insensitive comparison of the sorted token sets
    s1 = " ".join(sorted(inter)) or ""
    s2a = " ".join(sorted(inter | (ta - tb)))
    s2b = " ".join(sorted(inter | (tb - ta)))
    return max(
        SequenceMatcher(None, s1, s2a).ratio(),
        SequenceMatcher(None, s1, s2b).ratio(),
        SequenceMatcher(None, a, b).ratio(),
    )


def _contains_term(text: str, synonym: str) -> bool:
    """Substring test that respects word boundaries.

    A bare `in` test lets short synonyms match inside unrelated words: "apa"
    (alternative polyadenylation) inside "capability", "upr" inside
    "upregulation", "nmd" inside "nmda receptor current". Requiring word
    boundaries on both sides removes that whole class of false mapping, which
    otherwise sends records to the wrong phenotype silently.
    """
    if synonym not in text:
        return False
    return re.search(
        r"(?<![a-z0-9])" + re.escape(synonym) + r"(?![a-z0-9])", text
    ) is not None


class PhenotypeMatcher:
    """Maps free-text phenotype descriptions onto the controlled vocabulary."""

    def __init__(self, fuzzy_threshold: float = 0.86):
        self.fuzzy_threshold = fuzzy_threshold
        # Pre-sort synonyms longest-first so "mepsc frequency" beats "frequency"
        self._synonyms = sorted(SYNONYM_INDEX.items(), key=lambda kv: -len(kv[0]))
        self.unmapped_counter: dict[str, int] = {}

    @lru_cache(maxsize=8192)
    def _match_cached(self, text: str) -> tuple[str | None, float, str]:
        n = normalize_string(text)
        if not n:
            return (None, 0.0, "empty")

        if n in SYNONYM_INDEX:
            return (SYNONYM_INDEX[n], 1.0, "exact")

        # Substring: find the longest synonym contained in the text
        best_sub, best_len = None, 0
        for syn, tid in self._synonyms:
            if len(syn) > best_len and _contains_term(n, syn):
                best_sub, best_len = tid, len(syn)
        if best_sub:
            coverage = best_len / max(len(n), 1)
            return (best_sub, min(0.95, 0.6 + 0.35 * coverage), "synonym_substring")

        best_id, best_score = None, 0.0
        for syn, tid in self._synonyms:
            sc = _token_set_ratio(n, syn)
            if sc > best_score:
                best_id, best_score = tid, sc
        if best_score >= self.fuzzy_threshold:
            return (best_id, best_score * 0.9, "fuzzy")

        return (None, 0.0, "unmapped")

    def match(self, text: str) -> NormalizedTerm:
        tid, conf, method = self._match_cached(text)
        if tid is None:
            key = normalize_string(text)
            if key:
                self.unmapped_counter[key] = self.unmapped_counter.get(key, 0) + 1
            return NormalizedTerm(raw=text, mapping_confidence=0.0,
                                  mapping_method=method)
        term = BY_ID[tid]
        return NormalizedTerm(
            raw=text, ontology_id=tid, ontology_label=term.label,
            mapping_confidence=conf, mapping_method=method,
        )

    def unmapped_report(self, top_n: int = 50) -> list[tuple[str, int]]:
        """Frequent unmapped strings are the highest-value vocabulary additions."""
        return sorted(self.unmapped_counter.items(), key=lambda kv: -kv[1])[:top_n]


# --------------------------------------------------------------------------
# Assay classification
# --------------------------------------------------------------------------

# Method terms name the technique explicitly and outrank marker mentions.
# Without this split, "western blot for PSD95" classifies as imaging because
# PSD95 is also a puncta marker -- the antibody target is not the method.
METHOD_TERMS: dict[AssayClass, tuple[str, ...]] = {
    AssayClass.WESTERN_BLOT: ("western blot", "immunoblot", "sds-page", "blotting"),
    AssayClass.QPCR: ("qpcr", "rt-qpcr", "quantitative pcr", "real-time pcr",
                      "taqman", "ddpcr"),
    AssayClass.PROTEOMICS: ("mass spectrometry", "proteomic", "lc-ms", "tmt",
                            "silac", "label-free quant"),
    AssayClass.ELECTROPHYS_PATCH: ("patch clamp", "patch-clamp", "whole-cell",
                                   "whole cell recording", "voltage clamp",
                                   "current clamp"),
    AssayClass.ELECTROPHYS_MEA: ("multielectrode", "multi-electrode",
                                 "microelectrode array", "maestro"),
    AssayClass.CALCIUM_IMAGING: ("calcium imaging", "gcamp", "fluo-4", "fura-2",
                                 "ca2+ imaging", "calcium indicator"),
    AssayClass.SPLICING_ANALYSIS: ("rmats", "majiq", "leafcutter", "vast-tools",
                                   "suppa", "whippet", "minigene",
                                   "splicing analysis", "percent spliced"),
    AssayClass.RNASEQ_SINGLE_CELL: ("single-cell rna", "scrna-seq", "snrna-seq",
                                    "single nucleus rna", "10x genomics"),
    AssayClass.RNASEQ_BULK: ("rna-seq", "rnaseq", "rna sequencing", "bulk rna"),
    AssayClass.IMAGING_SYNAPTIC_PUNCTA: ("immunocytochem", "immunofluoresc",
                                         "confocal", "puncta analysis",
                                         "puncta quantification"),
    AssayClass.BIOCHEMICAL_BINDING: ("surface plasmon", "pulldown",
                                     "co-immunoprecipitation", "co-ip",
                                     "clip-seq", "eclip", "rip-seq"),
    AssayClass.VIABILITY: ("cell titer", "celltiter", "mtt assay", "ldh",
                           "trypan blue", "viability assay", "live/dead"),
}

# Higher rank = more specific. A downstream analysis method (rMATS) tells us
# more about what the experiment IS than the platform it ran on (RNA-seq).
SPECIFICITY_RANK: dict[AssayClass, int] = {
    AssayClass.SPLICING_ANALYSIS: 3,
    AssayClass.RNASEQ_SINGLE_CELL: 2,
    AssayClass.BIOCHEMICAL_BINDING: 2,
    AssayClass.CALCIUM_IMAGING: 2,
    AssayClass.ELECTROPHYS_PATCH: 2,
    AssayClass.ELECTROPHYS_MEA: 2,
    AssayClass.WESTERN_BLOT: 2,
    AssayClass.QPCR: 2,
    AssayClass.PROTEOMICS: 2,
    AssayClass.VIABILITY: 2,
    AssayClass.IMAGING_SYNAPTIC_PUNCTA: 1,
    AssayClass.RNASEQ_BULK: 0,
}

ASSAY_PATTERNS: list[tuple[AssayClass, tuple[str, ...]]] = [
    (AssayClass.IMAGING_SYNAPTIC_PUNCTA,
     ("puncta", "synapsin", "psd95", "psd-95", "homer", "bassoon", "vglut",
      "synaptophysin", "colocaliz", "gephyrin", "synaptic marker")),
    (AssayClass.ELECTROPHYS_PATCH,
     ("patch clamp", "patch-clamp", "whole-cell", "whole cell recording",
      "mepsc", "epsc", "ipsc recording", "voltage clamp", "current clamp")),
    (AssayClass.ELECTROPHYS_MEA,
     ("multielectrode", "multi-electrode", "mea", "maestro", "microelectrode array",
      "network burst", "spike sorting")),
    (AssayClass.CALCIUM_IMAGING,
     ("calcium imaging", "gcamp", "fluo-4", "fura-2", "ca2+ imaging",
      "calcium transient", "calcium indicator")),
    (AssayClass.SPLICING_ANALYSIS,
     ("rmats", "majiq", "leafcutter", "psi", "percent spliced", "vast-tools",
      "splicing analysis", "suppa", "whippet", "rt-pcr splicing", "minigene")),
    (AssayClass.RNASEQ_SINGLE_CELL,
     ("single-cell rna", "scrna-seq", "snrna-seq", "single nucleus rna", "10x genomics")),
    (AssayClass.RNASEQ_BULK,
     ("rna-seq", "rnaseq", "rna sequencing", "transcriptomic", "bulk rna")),
    (AssayClass.PROTEOMICS,
     ("mass spectrometry", "proteomic", "lc-ms", "tmt", "silac", "label-free quant")),
    (AssayClass.WESTERN_BLOT,
     ("western blot", "immunoblot", "sds-page")),
    (AssayClass.QPCR,
     ("qpcr", "rt-qpcr", "quantitative pcr", "real-time pcr", "taqman", "ddpcr")),
    (AssayClass.VIABILITY,
     ("cell titer", "celltiter", "mtt assay", "ldh", "propidium iodide",
      "annexin", "caspase-3", "trypan blue", "viability assay", "live/dead")),
    (AssayClass.IMAGING_MORPHOLOGY,
     ("sholl", "neurite outgrowth", "morphometr", "tracing", "imagej neurite",
      "spine density", "dendritic arbor")),
    (AssayClass.BIOCHEMICAL_BINDING,
     ("surface plasmon", "itc", "pulldown", "co-immunoprecipitation", "co-ip",
      "binding assay", "clip-seq", "eclip", "rip-seq")),
]


def classify_assay(text: str) -> tuple[AssayClass, float]:
    """Classify an assay description. Returns (class, confidence).

    Ordering matters: more specific patterns are checked first, since
    "RNA-seq used for splicing analysis" should classify as splicing.
    """
    n = normalize_string(text)
    if not n:
        return (AssayClass.OTHER, 0.0)

    # Tier 1: explicit method terms win outright.
    method_hits: list[tuple[AssayClass, int, int, int]] = []
    for cls, pats in METHOD_TERMS.items():
        matched = [p for p in pats if p in n]
        if matched:
            method_hits.append((
                cls,
                SPECIFICITY_RANK.get(cls, 0),
                len(matched),
                max(len(p) for p in matched),
            ))
    if method_hits:
        # Specificity first: "rMATS analysis of RNA-seq" is a splicing study
        # that happens to use RNA-seq as its platform, not an expression study.
        # The analysis method names the experiment; the platform does not.
        method_hits.sort(key=lambda t: (-t[1], -t[2], -t[3]))
        cls, _, count, _ = method_hits[0]
        return (cls, min(0.95, 0.7 + 0.1 * count))

    # Tier 2: fall back to marker/keyword association.
    hits: list[tuple[AssayClass, int]] = []
    for cls, pats in ASSAY_PATTERNS:
        c = sum(1 for p in pats if p in n)
        if c:
            hits.append((cls, c))
    if not hits:
        return (AssayClass.OTHER, 0.0)
    hits.sort(key=lambda kv: -kv[1])
    top, count = hits[0]
    conf = min(0.8, 0.45 + 0.12 * count)
    return (top, conf)


# --------------------------------------------------------------------------
# Gene symbol normalization
# --------------------------------------------------------------------------

# Words that look like gene symbols but are not. This list prevents a large
# class of false positives that plague naive gene NER.
GENE_STOPLIST = {
    "DNA", "RNA", "MRNA", "CDNA", "PCR", "QPCR", "WT", "KO", "KD", "OE",
    "CI", "SD", "SEM", "ANOVA", "ROI", "FOV", "DIV", "PBS", "BSA", "DMSO",
    "FBS", "EDTA", "ATP", "ADP", "GTP", "NAD", "ROS", "ER", "AD", "PD", "ALS",
    "FTD", "ASD", "SZ", "MDD", "USA", "UK", "PFA", "HRP", "GFP", "RFP", "YFP",
    "MRI", "EEG", "IQ", "II", "III", "IV", "AAV", "MOI", "IF", "ICC", "IHC",
    "FDR", "TPM", "FPKM", "RPKM", "CPM", "PSI", "UTR", "ORF", "SNP", "INDEL",
    "HEK", "PFC", "CNS", "PNS", "CSF", "BBB", "NIH", "P value", "N/A",
}

GENE_SYMBOL_RE = re.compile(r"\b([A-Z][A-Z0-9]{1,9}(?:-[A-Z0-9]{1,4})?)\b")

# Experimental-parameter prefixes that produce symbol-shaped tokens.
# DIV21, P30, E14.5, C57BL6 etc. are timepoints and strains, not genes.
_PARAM_PREFIX_RE = re.compile(
    r"^(DIV|DAY|D|P|E|W|WK|HR|H|MIN|N|FIG|TAB|SUPP|EXP|REP|LOT|CAT|RRID|"
    r"C57|BL|ICR|CD|NOD|SCID|MOI|RT|OD|MW|KDA|KB|BP|NM|UM|MM|ML|UL|NG|UG|MG)\d+$"
)


def _looks_like_parameter(token: str) -> bool:
    """Reject tokens that are experimental parameters rather than genes."""
    if _PARAM_PREFIX_RE.match(token):
        return True
    # Pure letter+number with no vowel-ish structure and a big number: PMID-like
    if re.match(r"^[A-Z]{1,3}\d{4,}$", token):
        return True
    return False


def extract_gene_candidates(text: str) -> list[str]:
    """Regex-based gene symbol candidate extraction.

    This is a recall-oriented first pass; the HGNC resolver filters it. We do
    not attempt to be precise here -- precision comes from requiring an HGNC
    hit downstream. We do filter the parameter-shaped tokens, because they are
    numerous enough to distort candidate counts before HGNC ever sees them.
    """
    if not text:
        return []
    cands = GENE_SYMBOL_RE.findall(text)
    out, seen = [], set()
    for c in cands:
        if c in GENE_STOPLIST or c in seen:
            continue
        if c.isdigit() or len(c) < 2:
            continue
        if _looks_like_parameter(c):
            continue
        seen.add(c)
        out.append(c)
    return out
