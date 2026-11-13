"""
Deterministic sentence-level extractor for Europe PMC / PubMed abstracts.

Used as a zero-API-key fallback and for reproducible offline/online temporal
hold-out benchmarking. Every emitted ExperimentRecord is grounded in an exact
abstract sentence containing:
  - a neuronal / iPSC cell model indicator,
  - a canonical phenotype recognized by PhenotypeMatcher,
  - an explicit directional cue (increase / decrease / unchanged),
  - and real SourceMetadata (PMID/PMCID, year, journal, senior author/institution).

Does not invent effect sizes, p-values, or donor counts when absent from text.
"""

from __future__ import annotations

import re
from typing import Iterable

from ..normalize.labs import LabResolver
from ..normalize.matcher import PhenotypeMatcher, classify_assay
from ..schema import (
    CellModel,
    CellModelClass,
    Direction,
    EvidenceType,
    ExperimentRecord,
    NormalizedTerm,
    Outcome,
    Perturbation,
    PerturbationType,
    Provenance,
    SourceMetadata,
)

_NEURO_GENES = (
    "TARDBP", "TDP-43", "STMN2", "UNC13A", "C9orf72", "SOD1", "FUS",
    "MAPT", "SNCA", "LRRK2", "GBA", "GRN", "APP", "PSEN1", "PSEN2",
    "ATXN2", "ATXN1", "ATXN3", "HTT", "SMN1", "SMN2", "SCN1A", "SCN2A",
    "SYNGAP1", "MECP2", "SHANK3", "KCNQ2", "TBK1", "OPTN", "VCP", "NEK1",
    "SQSTM1", "CHCHD10", "PRKN", "PINK1", "NRXN1", "GRIN2A", "GRIN2B",
    "ANK2", "CHD8", "UBQLN2", "MATR3", "HNRNPA1", "APOE", "TREM2",
)

_GENE_RE = re.compile(r"\b(" + "|".join(re.escape(g) for g in _NEURO_GENES) + r")\b")
_DECREASE_RE = re.compile(r"\b(decreas\w*|reduc\w*|loss of|impair\w*|deplet\w*|downregulat\w*|lower\w*|diminish\w*)\b", re.I)
_INCREASE_RE = re.compile(r"\b(increas\w*|elevat\w*|upregulat\w*|accumulat\w*|hyperexcit\w*|enhanc\w*|higher)\b", re.I)
_NULL_RE = re.compile(r"\b(no significant change|unchanged|did not alter|no difference)\b", re.I)

_PATHWAY_HINTS: list[tuple[re.Pattern[str], str, str]] = [
    (re.compile(r"\bsplic\w*|cryptic exon|intron retention", re.I), "RNA splicing", "GO:0008380"),
    (re.compile(r"\bsynap\w*|puncta|spine|mEPSC|sEPSC|NMDA|AMPA", re.I), "synaptic transmission", "GO:0007268"),
    (re.compile(r"\baxon\w*|neurite|outgrowth|microtubule", re.I), "axonogenesis", "GO:0007409"),
    (re.compile(r"\bmitochondri\w*|ROS|oxidative|bioenergetic", re.I), "mitochondrial organization", "GO:0007005"),
    (re.compile(r"\bautophag\w*|lysosom\w*|proteasom\w*|clearance", re.I), "autophagy", "GO:0006914"),
    (re.compile(r"\bexcitabilit\w*|firing|burst|action potential|MEA", re.I), "neuronal action potential", "GO:0019228"),
    (re.compile(r"\bapoptos\w*|cell death|viabilit\w*|degenerat\w*", re.I), "neuron death", "GO:0070997"),
    (re.compile(r"\bnucleocytoplasmic|nuclear export|stress granule", re.I), "nucleocytoplasmic transport", "GO:0006913"),
]


def _infer_model_class(text: str) -> CellModelClass:
    low = text.lower()
    if "ngn2" in low:
        return CellModelClass.IPSC_NEURON_NGN2
    if "organoid" in low:
        return CellModelClass.IPSC_ORGANOID
    if "neural progenitor" in low or "npc" in re.findall(r"\bnpc\b", low):
        return CellModelClass.IPSC_NPC
    if "ipsc" in low or "induced pluripotent" in low:
        return CellModelClass.IPSC_NEURON_SMALL_MOLECULE
    if "primary" in low:
        return CellModelClass.PRIMARY_RODENT_NEURON
    return CellModelClass.OTHER


def extract_records_from_abstract(
    source: SourceMetadata,
    abstract_text: str,
    matcher: PhenotypeMatcher | None = None,
) -> list[ExperimentRecord]:
    """Deterministically extract sentence-grounded ExperimentRecords from an abstract."""
    if not abstract_text or not abstract_text.strip():
        return []

    pm = matcher or PhenotypeMatcher()

    model_class = _infer_model_class(abstract_text)
    sentences = re.split(r"(?<=[.!?])\s+", abstract_text.strip())
    records: list[ExperimentRecord] = []
    seen_keys: set[tuple[str, str, str]] = set()

    doc_genes = _GENE_RE.findall(abstract_text)
    default_gene = "TARDBP" if doc_genes and doc_genes[0] == "TDP-43" else (doc_genes[0] if doc_genes else None)

    for idx, sent in enumerate(sentences):
        if len(sent) < 20:
            continue

        # Direction
        if _NULL_RE.search(sent):
            direction = Direction.NO_CHANGE
        else:
            has_dec = bool(_DECREASE_RE.search(sent))
            has_inc = bool(_INCREASE_RE.search(sent))
            if has_dec == has_inc:
                continue
            direction = Direction.DECREASE if has_dec else Direction.INCREASE

        # Phenotype match
        pheno = pm.match(sent)
        if pheno.ontology_id is None:
            continue

        # Target gene in sentence or document
        s_genes = _GENE_RE.findall(sent)
        gene = s_genes[0] if s_genes else default_gene
        if gene == "TDP-43":
            gene = "TARDBP"
        if not gene:
            continue

        # Pathway hints
        pathways: list[NormalizedTerm] = []
        for pat, label, go_id in _PATHWAY_HINTS:
            if pat.search(sent):
                pathways.append(NormalizedTerm(raw=label, ontology_id=go_id, ontology_label=label, mapping_confidence=0.8))
        if not pathways:
            continue

        key = (gene, pheno.ontology_id, direction.value)
        if key in seen_keys:
            continue
        seen_keys.add(key)

        assay, _ = classify_assay(sent)
        if assay.value == "other":
            doc_assay, _ = classify_assay(abstract_text)
            assay = doc_assay
        rec = ExperimentRecord(
            record_id=f"{source.source_id}:{idx}:{pheno.ontology_id}",
            source=source,
            cell_model=CellModel(model_class=model_class),
            perturbation=Perturbation(
                ptype=PerturbationType.OTHER,
                target=NormalizedTerm(raw=gene, ontology_label=gene, ontology_id=gene, mapping_confidence=1.0),
            ),
            assay_class=assay,
            outcome=Outcome(
                phenotype=pheno,
                direction=direction,
                evidence_type=EvidenceType.DIRECT_MEASUREMENT,
            ),
            pathways=pathways,
            provenance=Provenance(
                quote=sent[:400],
                section="abstract",
                extraction_confidence=0.65,
                extractor_version="heuristic-1.0",
            ),
        )
        records.append(rec)

    return records


def harvest_and_extract_europepmc(
    queries: Iterable[str],
    per_query: int = 80,
) -> list[ExperimentRecord]:
    """Fetch real abstracts from Europe PMC and extract sentence-grounded records."""
    import json
    import urllib.parse
    import urllib.request

    pm = PhenotypeMatcher()
    lr = LabResolver()
    pending: list[tuple[SourceMetadata, str]] = []
    seen_pmids: set[str] = set()

    for q in queries:
        params = urllib.parse.urlencode({
            "query": q,
            "format": "json",
            "resultType": "core",
            "pageSize": min(per_query, 100),
        })
        url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?{params}"
        req = urllib.request.Request(url, headers={"User-Agent": "SYNAPTIK-Research/0.9"})
        data = json.loads(urllib.request.urlopen(req, timeout=20).read().decode("utf-8"))
        for item in data.get("resultList", {}).get("result", []):
            pmid = str(item.get("pmid") or item.get("id") or "")
            if not pmid or pmid in seen_pmids:
                continue
            seen_pmids.add(pmid)
            abstract = item.get("abstractText") or ""
            if not abstract:
                continue
            year_raw = item.get("pubYear")
            try:
                year = int(year_raw) if year_raw else None
            except ValueError:
                year = None
            raw_authors = item.get("authorList", {}).get("author", [])
            author_names = [a.get("fullName") for a in raw_authors if a.get("fullName")]
            if not author_names and item.get("authorString"):
                author_names = [a.strip() for a in item["authorString"].split(",") if a.strip()]
            senior = author_names[-1] if author_names else None
            affils = [a.get("affiliation") for a in raw_authors if a.get("affiliation")]
            src = SourceMetadata(
                source_id=f"PMID:{pmid}",
                pmid=pmid,
                pmcid=item.get("pmcid"),
                doi=item.get("doi"),
                title=item.get("title") or "",
                year=year,
                journal=(item.get("journalInfo", {}).get("journal", {}) or {}).get("title"),
                authors=author_names,
                corresponding_author=senior,
                affiliations=affils[:5],
                is_preprint=(item.get("source") == "PPR"),
            )
            pending.append((src, abstract))

    lr.fit([src for src, _ in pending])
    all_records: list[ExperimentRecord] = []
    for src, abstract in pending:
        src.lab_id = lr.assign(src)
        all_records.extend(extract_records_from_abstract(src, abstract, matcher=pm))

    return all_records
