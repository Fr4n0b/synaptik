"""
World-scale literature harvesting.

Realistic coverage, stated honestly
-----------------------------------
"All the world's literature" is not achievable and any tool claiming it is
lying. Here is what is actually reachable, and what is not:

  REACHABLE
    PubMed abstracts             ~37M records, all biomedical literature
    Europe PMC full text         ~4.5M open-access articles
    PubMed Central OA subset     ~5M, bulk-downloadable
    bioRxiv / medRxiv            ~500k preprints, full text
    OpenAlex                     ~250M works, all disciplines, metadata
    Semantic Scholar             ~200M papers, citation graph, some full text
    CORE                         ~300M open-access documents
    Crossref                     ~150M DOIs, metadata + references
    Unpaywall                    OA location for ~50M paywalled DOIs
    Retraction Watch             ~50k retractions -- critical for integrity
    Zenodo / OSF / figshare      preprints, datasets, theses

  NOT REACHABLE without institutional licensing
    Elsevier / Springer / Wiley full text behind paywalls (~25M articles)
    Most pre-2000 literature not retrospectively digitised as full text
    Supplementary data behind publisher walls
    Non-English literature with no English abstract

  PARTIALLY REACHABLE
    Institutional subscriptions -- if your library has an API or proxy, the
    InstitutionalFullText adapter can use it. This is the legitimate route to
    paywalled content and requires YOUR credentials and YOUR compliance with
    the publisher agreement. The tool does not circumvent paywalls.

What this means in practice: near-complete coverage at abstract level, and
roughly 20-25% coverage at full-text level, heavily skewed toward recent and
well-funded work. The extraction quality difference between abstract-only and
full-text is large -- Methods sections carry the protocol details -- so the
corpus report tracks the two separately and analyses flag which they rest on.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from datetime import date, datetime
from typing import Iterator

from ..schema import SourceMetadata
from .http import PoliteSession

log = logging.getLogger(__name__)

OPENALEX = "https://api.openalex.org"
SEMANTIC_SCHOLAR = "https://api.semanticscholar.org/graph/v1"
CROSSREF = "https://api.crossref.org"
UNPAYWALL = "https://api.unpaywall.org/v2"
CORE = "https://api.core.ac.uk/v3"
RETRACTION_WATCH = (
    "https://gitlab.com/crossref/retraction-watch-data/-/raw/main/"
    "retraction_watch.csv"
)


@dataclass
class CoverageReport:
    """What fraction of the reachable literature we actually have."""

    total_records: int = 0
    with_abstract: int = 0
    with_full_text: int = 0
    open_access: int = 0
    preprints: int = 0
    retracted: int = 0
    by_source: dict[str, int] = field(default_factory=dict)
    by_decade: dict[str, int] = field(default_factory=dict)
    estimated_universe: int | None = None

    def full_text_fraction(self) -> float:
        return self.with_full_text / max(1, self.total_records)

    def narrative(self) -> str:
        ft = self.full_text_fraction()
        lines = [
            f"{self.total_records:,} records; {self.with_full_text:,} "
            f"({ft:.0%}) have full text."
        ]
        if ft < 0.3:
            lines.append(
                "Most of this corpus is abstract-only. Extraction from "
                "abstracts misses protocol details (DIV, differentiation "
                "method, antibodies) that live in Methods, so experiment "
                "lookup will be less precise for those records."
            )
        if self.estimated_universe:
            cov = self.total_records / self.estimated_universe
            lines.append(
                f"Estimated {cov:.0%} of the reachable literature on these "
                f"topics."
            )
        return " ".join(lines)


class OpenAlexHarvester:
    """OpenAlex: ~250M works, fully open, no API key required.

    The broadest single source available. Provides metadata, concepts,
    citation counts, institutional affiliations (already disambiguated via
    ROR), and open-access locations. Its author and institution
    disambiguation is substantially better than anything we could compute
    ourselves, so we prefer it for lab clustering when available.
    """

    def __init__(self, session: PoliteSession, mailto: str | None = None):
        self.s = session
        self.mailto = mailto

    def _params(self, extra: dict) -> dict:
        p = dict(extra)
        if self.mailto:
            p["mailto"] = self.mailto  # polite pool: faster and higher limits
        return p

    def search(
        self, query: str, max_results: int = 10000,
        from_year: int | None = None, per_page: int = 200,
    ) -> Iterator[SourceMetadata]:
        cursor = "*"
        fetched = 0
        filters = ["type:article"]
        if from_year:
            filters.append(f"from_publication_date:{from_year}-01-01")

        while fetched < max_results:
            params = self._params({
                "search": query,
                "per-page": min(per_page, max_results - fetched),
                "cursor": cursor,
                "filter": ",".join(filters),
            })
            data = self.s.get_json(f"{OPENALEX}/works", params)
            if not data:
                return
            results = data.get("results") or []
            if not results:
                return
            for w in results:
                meta = self._parse(w)
                if meta:
                    fetched += 1
                    yield meta
            cursor = (data.get("meta") or {}).get("next_cursor")
            if not cursor:
                return

    def _parse(self, w: dict) -> SourceMetadata | None:
        ids = w.get("ids") or {}
        doi = (ids.get("doi") or "").replace("https://doi.org/", "") or None
        pmid = (ids.get("pmid") or "").split("/")[-1] or None
        pmcid = (ids.get("pmcid") or "").split("/")[-1] or None
        source_id = (
            f"PMID:{pmid}" if pmid else f"DOI:{doi}" if doi
            else f"OPENALEX:{w.get('id', '').split('/')[-1]}"
        )
        if source_id.endswith(":"):
            return None

        host = (w.get("primary_location") or {}).get("source") or {}
        authors, affiliations, ror_ids = [], [], []
        for a in (w.get("authorships") or []):
            name = (a.get("author") or {}).get("display_name")
            if name:
                authors.append(name)
            for inst in (a.get("institutions") or []):
                if inst.get("display_name"):
                    affiliations.append(inst["display_name"])
                if inst.get("ror"):
                    ror_ids.append(inst["ror"])

        pub_date = None
        if w.get("publication_date"):
            try:
                pub_date = datetime.strptime(
                    w["publication_date"], "%Y-%m-%d"
                ).date()
            except ValueError:
                pass

        oa = w.get("open_access") or {}
        meta = SourceMetadata(
            source_id=source_id, pmid=pmid, pmcid=pmcid, doi=doi,
            title=(w.get("title") or w.get("display_name") or "").strip(),
            journal=host.get("display_name"),
            publication_date=pub_date,
            year=w.get("publication_year"),
            is_preprint=(w.get("type") == "preprint")
            or bool(host.get("type") == "repository"),
            is_review=(w.get("type") == "review"),
            is_retracted=bool(w.get("is_retracted")),
            authors=authors,
            affiliations=list(dict.fromkeys(affiliations)),
            citation_count=w.get("cited_by_count"),
            open_access=bool(oa.get("is_oa")),
            full_text_available=bool(pmcid) or bool(oa.get("is_oa")),
        )
        # ROR IDs are pre-disambiguated institution identifiers; carried
        # out-of-band so the lab resolver can prefer them over string matching.
        meta.__dict__["_ror_ids"] = list(dict.fromkeys(ror_ids))
        meta.__dict__["_oa_url"] = oa.get("oa_url")
        meta.__dict__["_concepts"] = [
            c.get("display_name") for c in (w.get("concepts") or [])[:8]
        ]
        return meta

    def count(self, query: str) -> int:
        data = self.s.get_json(
            f"{OPENALEX}/works",
            self._params({"search": query, "per-page": 1}),
        )
        return int(((data or {}).get("meta") or {}).get("count", 0))


class SemanticScholarHarvester:
    """Semantic Scholar: citation graph and influence metrics.

    Its distinctive value here is `influentialCitationCount` -- citations
    where the citing paper actually builds on the cited one, rather than
    mentioning it in passing. That is a far better relevance signal than raw
    citation count, which is dominated by reviews and methods papers.
    """

    FIELDS = (
        "paperId,externalIds,title,abstract,year,venue,authors,"
        "citationCount,influentialCitationCount,isOpenAccess,"
        "openAccessPdf,publicationTypes,publicationDate,fieldsOfStudy"
    )

    def __init__(self, session: PoliteSession, api_key: str | None = None):
        self.s = session
        self.api_key = api_key
        if api_key:
            self.s.session.headers["x-api-key"] = api_key

    def search(
        self, query: str, max_results: int = 1000, per_page: int = 100
    ) -> Iterator[SourceMetadata]:
        offset = 0
        while offset < max_results:
            data = self.s.get_json(
                f"{SEMANTIC_SCHOLAR}/paper/search",
                {"query": query, "offset": offset,
                 "limit": min(per_page, max_results - offset),
                 "fields": self.FIELDS},
            )
            if not data:
                return
            papers = data.get("data") or []
            if not papers:
                return
            for p in papers:
                meta = self._parse(p)
                if meta:
                    yield meta
            offset += len(papers)
            if data.get("next") is None:
                return

    def citations(self, doi: str, limit: int = 100) -> list[dict]:
        """Papers citing this one. Used to trace a finding's downstream fate."""
        data = self.s.get_json(
            f"{SEMANTIC_SCHOLAR}/paper/DOI:{doi}/citations",
            {"fields": "title,year,externalIds,isInfluential", "limit": limit},
        )
        return (data or {}).get("data") or []

    def references(self, doi: str, limit: int = 100) -> list[dict]:
        data = self.s.get_json(
            f"{SEMANTIC_SCHOLAR}/paper/DOI:{doi}/references",
            {"fields": "title,year,externalIds,isInfluential", "limit": limit},
        )
        return (data or {}).get("data") or []

    @staticmethod
    def _parse(p: dict) -> SourceMetadata | None:
        ext = p.get("externalIds") or {}
        pmid = str(ext["PubMed"]) if ext.get("PubMed") else None
        doi = ext.get("DOI")
        pmcid = str(ext["PubMedCentral"]) if ext.get("PubMedCentral") else None
        source_id = (
            f"PMID:{pmid}" if pmid else f"DOI:{doi}" if doi
            else f"S2:{p.get('paperId')}"
        )
        types = [t.lower() for t in (p.get("publicationTypes") or [])]
        pub_date = None
        if p.get("publicationDate"):
            try:
                pub_date = datetime.strptime(
                    p["publicationDate"], "%Y-%m-%d"
                ).date()
            except ValueError:
                pass

        meta = SourceMetadata(
            source_id=source_id, pmid=pmid, pmcid=pmcid, doi=doi,
            title=(p.get("title") or "").strip(),
            journal=p.get("venue"),
            year=p.get("year"),
            publication_date=pub_date,
            is_review="review" in types,
            is_preprint=any("preprint" in t for t in types),
            authors=[a.get("name") for a in (p.get("authors") or [])
                     if a.get("name")],
            citation_count=p.get("citationCount"),
            open_access=bool(p.get("isOpenAccess")),
            full_text_available=bool(pmcid or p.get("openAccessPdf")),
        )
        meta.__dict__["_influential_citations"] = p.get(
            "influentialCitationCount", 0
        )
        meta.__dict__["_abstract"] = p.get("abstract") or ""
        pdf = p.get("openAccessPdf") or {}
        meta.__dict__["_oa_pdf"] = pdf.get("url")
        return meta


class UnpaywallResolver:
    """Finds legal open-access copies of paywalled papers.

    Unpaywall indexes author manuscripts, repository deposits, and publisher
    OA versions. It routinely finds a legitimate free copy of a paper that
    appears paywalled -- typically raising full-text coverage by 10-20
    percentage points, entirely within publisher policy.
    """

    def __init__(self, session: PoliteSession, email: str):
        self.s = session
        self.email = email  # required by the API

    def find_oa(self, doi: str) -> dict | None:
        data = self.s.get_json(f"{UNPAYWALL}/{doi}", {"email": self.email})
        if not data or not data.get("is_oa"):
            return None
        best = data.get("best_oa_location") or {}
        return {
            "url": best.get("url_for_pdf") or best.get("url"),
            "version": best.get("version"),
            "license": best.get("license"),
            "host_type": best.get("host_type"),
            "is_oa": True,
        }


class RetractionWatchIndex:
    """Retraction Watch database, released openly via Crossref.

    This matters more than it sounds. Retracted papers continue to be cited
    for years. A literature-mining tool that does not check retraction status
    will faithfully aggregate fabricated data and present it with confidence.
    The evidence weighter zeroes retracted records, but only if it knows.
    """

    def __init__(self, session: PoliteSession):
        self.s = session
        self.by_doi: dict[str, dict] = {}
        self.by_pmid: dict[str, dict] = {}
        self.loaded = False

    def load(self, local_csv: str | None = None) -> int:
        import csv
        import io

        text = None
        if local_csv:
            try:
                text = open(local_csv, encoding="utf-8", errors="replace").read()
            except OSError as e:
                log.warning("Could not read %s: %s", local_csv, e)
        if text is None:
            text = self.s.get_text(RETRACTION_WATCH)
        if not text:
            log.warning(
                "Retraction Watch data unavailable. Retracted papers will not "
                "be down-weighted unless flagged by their source database."
            )
            return 0

        n = 0
        for row in csv.DictReader(io.StringIO(text)):
            doi = (row.get("OriginalPaperDOI") or "").strip().lower()
            pmid = (row.get("OriginalPaperPubMedID") or "").strip()
            entry = {
                "reason": row.get("Reason", ""),
                "retraction_date": row.get("RetractionDate", ""),
                "nature": row.get("RetractionNature", ""),
                "title": row.get("Title", ""),
            }
            if doi:
                self.by_doi[doi] = entry
                n += 1
            if pmid and pmid != "0":
                self.by_pmid[pmid] = entry
        self.loaded = True
        log.info("Loaded %d retraction records", n)
        return n

    def check(self, meta: SourceMetadata) -> dict | None:
        if not self.loaded:
            return None
        if meta.doi and meta.doi.lower() in self.by_doi:
            return self.by_doi[meta.doi.lower()]
        if meta.pmid and meta.pmid in self.by_pmid:
            return self.by_pmid[meta.pmid]
        return None

    def annotate(self, meta: SourceMetadata) -> bool:
        """Mark a source as retracted in place. Returns True if it was."""
        hit = self.check(meta)
        if hit:
            meta.is_retracted = True
            nature = (hit.get("nature") or "").lower()
            meta.retraction_notice = (
                f"{hit.get('nature', 'Retraction')}: "
                f"{hit.get('reason', 'reason not stated')} "
                f"({hit.get('retraction_date', 'date unknown')})"
            )
            # "Expression of concern" is not a retraction; flag but do not zero.
            if "concern" in nature:
                meta.is_retracted = False
                meta.__dict__["_expression_of_concern"] = True
            return meta.is_retracted
        return False


class InstitutionalFullText:
    """Adapter for institutional subscription access.

    Deliberately minimal and deliberately opt-in. This does NOT circumvent
    paywalls: it uses credentials the user's institution has legitimately
    licensed, subject to the terms that institution agreed to. Text and data
    mining rights vary by publisher and by contract; most major publishers
    grant TDM rights to subscribing institutions on request, often through a
    dedicated API.

    Before enabling this, confirm with your library that TDM is permitted
    under your agreements. The tool cannot verify that for you.
    """

    def __init__(
        self, session: PoliteSession, base_url: str | None = None,
        api_key: str | None = None, enabled: bool = False,
    ):
        self.s = session
        self.base_url = base_url
        self.api_key = api_key
        self.enabled = enabled and bool(base_url)
        if enabled and not base_url:
            log.warning(
                "Institutional access requested but no base URL configured."
            )

    def fetch(self, doi: str) -> str | None:
        if not self.enabled:
            return None
        headers = {"Accept": "text/xml"}
        if self.api_key:
            headers["X-API-Key"] = self.api_key
        try:
            self.s.limiter.wait(self.base_url.split("/")[2])
            resp = self.s.session.get(
                f"{self.base_url.rstrip('/')}/{doi}",
                headers=headers, timeout=self.s.timeout,
            )
            if resp.status_code == 200:
                return resp.text
            if resp.status_code in (401, 403):
                log.warning(
                    "Institutional access denied for %s (HTTP %s). Check "
                    "credentials and TDM entitlement.", doi, resp.status_code
                )
        except Exception as e:
            log.warning("Institutional fetch failed for %s: %s", doi, e)
        return None


# --------------------------------------------------------------------------
# Expanded corpus definition
# --------------------------------------------------------------------------

# Beyond the v1 iPSC-neuron queries: these broaden coverage to the surrounding
# literature that gives the association statistics a meaningful background,
# and to the disease areas where iPSC-neuron models are actually applied.
WORLD_CORPUS_QUERIES: dict[str, str] = {
    # --- core model systems ---
    "ipsc_neuron_core": (
        '"iPSC-derived neuron" OR "iPSC neuron" OR "hiPSC-derived neuron" OR '
        '"induced pluripotent stem cell-derived neuron"'
    ),
    "ngn2_induction": (
        'NGN2 OR NEUROG2 OR "Neurogenin-2" OR "induced neuron" OR iNeuron'
    ),
    "brain_organoid": (
        '"brain organoid" OR "cerebral organoid" OR "cortical organoid" OR '
        '"neural organoid" OR assembloid'
    ),
    "isogenic_ipsc": (
        '(isogenic OR "isogenic control") AND (iPSC OR "stem cell-derived") '
        'AND (neuron OR neural)'
    ),

    # --- phenotype domains ---
    "synaptic_phenotype": (
        '(synapse OR synaptic OR "synaptic puncta" OR PSD95 OR synapsin OR '
        'homer1 OR bassoon OR VGLUT1 OR gephyrin) AND '
        '(iPSC OR "stem cell-derived" OR "induced neuron")'
    ),
    "electrophysiology": (
        '(mEPSC OR "miniature excitatory postsynaptic" OR "multielectrode '
        'array" OR "patch clamp" OR "network burst" OR "calcium imaging") AND '
        '(iPSC OR "induced neuron" OR "human neuron")'
    ),
    "neurite_morphology": (
        '("neurite outgrowth" OR "dendritic arborization" OR "Sholl analysis" '
        'OR "spine density" OR "axon length") AND (neuron OR neuronal)'
    ),

    # --- splicing, deep ---
    "splicing_neuro": (
        '("alternative splicing" OR "mis-splicing" OR "cryptic exon" OR '
        '"intron retention" OR "exon skipping" OR "splice site") AND '
        '(neuron OR neuronal OR neurodegeneration OR brain)'
    ),
    "rbp_regulation": (
        '(TDP-43 OR TARDBP OR FUS OR hnRNP OR RBFOX1 OR RBFOX2 OR NOVA1 OR '
        'NOVA2 OR PTBP1 OR PTBP2 OR SRSF1 OR MATR3 OR SFPQ OR ELAVL OR QKI OR '
        'MBNL1 OR CELF) AND (splicing OR "RNA binding")'
    ),
    "splice_therapeutics": (
        '(antisense oligonucleotide OR ASO OR "splice-switching" OR '
        'nusinersen OR "exon skipping therapy" OR "splicing modulator") AND '
        '(neuron OR neurological OR neuromuscular)'
    ),
    "splicing_methods": (
        'rMATS OR MAJIQ OR LeafCutter OR SUPPA OR "vast-tools" OR Whippet OR '
        '"percent spliced in" OR SpliceAI OR Pangolin'
    ),

    # --- disease contexts where iPSC neurons are used ---
    "als_ftd": (
        '(ALS OR "amyotrophic lateral sclerosis" OR FTD OR '
        '"frontotemporal dementia" OR C9orf72) AND '
        '(iPSC OR neuron OR "motor neuron")'
    ),
    "autism_ndd": (
        '(autism OR "autism spectrum" OR "neurodevelopmental disorder" OR '
        'SHANK3 OR SCN2A OR CHD8 OR MECP2 OR FMR1) AND '
        '(neuron OR iPSC OR synapse)'
    ),
    "schizophrenia_psychiatric": (
        '(schizophrenia OR bipolar OR "major depressive") AND '
        '(iPSC OR "stem cell-derived" OR neuron)'
    ),
    "alzheimers_tau": (
        '(Alzheimer OR tau OR MAPT OR amyloid OR APP OR PSEN1) AND '
        '(iPSC OR "induced neuron" OR "human neuron")'
    ),
    "parkinsons": (
        '(Parkinson OR SNCA OR LRRK2 OR PINK1 OR PRKN OR GBA) AND '
        '(iPSC OR "dopaminergic neuron")'
    ),
    "rare_neurogenetic": (
        '("rare disease" OR "de novo variant" OR "loss of function variant") '
        'AND (neuron OR neurodevelopmental) AND (iPSC OR patient-derived)'
    ),

    # --- mechanism and pathway ---
    "mitochondria_metabolism": (
        '(mitochondrial OR "oxidative phosphorylation" OR "metabolic stress") '
        'AND (neuron OR neuronal OR iPSC)'
    ),
    "proteostasis": (
        '(autophagy OR proteasome OR "protein aggregation" OR "unfolded '
        'protein response" OR chaperone) AND (neuron OR neurodegeneration)'
    ),
    "transcriptomics_neuron": (
        '("single-cell RNA" OR scRNA-seq OR snRNA-seq OR "single nucleus") '
        'AND (neuron OR brain OR cortical)'
    ),

    # --- therapeutic ---
    "drug_rescue": (
        '(rescue OR reversal OR "phenotypic rescue" OR repurposing) AND '
        '(iPSC OR "induced neuron") AND (compound OR drug OR "small molecule")'
    ),
    "crispr_screen_neuron": (
        '(CRISPR screen OR "genome-wide screen" OR CRISPRi OR CRISPRa) AND '
        '(neuron OR iPSC OR neural)'
    ),
}


def estimate_corpus_size(
    harvester, queries: dict[str, str] | None = None
) -> dict[str, int]:
    """Count hits per query before harvesting, so scale is known in advance."""
    queries = queries or WORLD_CORPUS_QUERIES
    out = {}
    for name, q in queries.items():
        try:
            out[name] = harvester.count(q)
        except Exception as e:
            log.warning("Count failed for %s: %s", name, e)
            out[name] = -1
    return out
