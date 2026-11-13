"""
PubMed (E-utilities) and bioRxiv/medRxiv harvesters.

Europe PMC mirrors most of PubMed, but not all: some records appear in PubMed
first, and MeSH terms -- which are genuinely useful for phenotype
normalization -- are only available through E-utilities. bioRxiv matters
because iPSC-neuron work is preprint-heavy and often 12-18 months ahead of the
journal version.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import logging
from datetime import date, datetime
from typing import Iterator

from lxml import etree

from ..schema import SourceMetadata
from .http import PoliteSession

log = logging.getLogger(__name__)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BIORXIV = "https://api.biorxiv.org"


class PubMedHarvester:
    """E-utilities client. An API key raises the rate limit 3x -- get one."""

    def __init__(self, session: PoliteSession, api_key: str | None = None,
                 drift=None):
        self.s = session
        self.api_key = api_key
        if drift is None:
            from .contracts import DriftDetector
            drift = DriftDetector()
        self.drift = drift

    def _auth(self, params: dict) -> dict:
        p = dict(params)
        p["tool"] = self.s.tool_name
        if self.s.contact_email:
            p["email"] = self.s.contact_email
        if self.api_key:
            p["api_key"] = self.api_key
        return p

    def search_pmids(self, query: str, max_results: int = 1000) -> list[str]:
        pmids: list[str] = []
        retstart = 0
        while len(pmids) < max_results:
            params = self._auth({
                "db": "pubmed", "term": query, "retmode": "json",
                "retmax": min(500, max_results - len(pmids)), "retstart": retstart,
            })
            data = self.s.get_json(f"{EUTILS}/esearch.fcgi", params)
            if not data:
                break
            self.drift.check("pubmed_esearch", data)
            idlist = (data.get("esearchresult") or {}).get("idlist") or []
            if not idlist:
                break
            pmids.extend(idlist)
            retstart += len(idlist)
            total = int((data.get("esearchresult") or {}).get("count", 0))
            if retstart >= total:
                break
        return pmids[:max_results]

    def fetch_metadata(
        self,
        pmids: list[str],
        batch_size: int = 200,
        max_workers: int = 3,
    ) -> Iterator[SourceMetadata]:
        """efetch in batches with concurrent worker threads; parses MeSH terms and publication types."""
        if not pmids:
            return

        batches = [pmids[i : i + batch_size] for i in range(0, len(pmids), batch_size)]

        def _fetch_and_parse_batch(batch: list[str]) -> list[SourceMetadata]:
            params = self._auth({
                "db": "pubmed", "id": ",".join(batch), "retmode": "xml",
            })
            xml = self.s.get_text(f"{EUTILS}/efetch.fcgi", params)
            if not xml:
                return []
            try:
                return list(self._parse_efetch(xml))
            except etree.XMLSyntaxError as e:
                log.warning("efetch parse failed for batch of %d PMIDs: %s", len(batch), e)
                return []

        if len(batches) <= 1 or max_workers <= 1:
            for b in batches:
                yield from _fetch_and_parse_batch(b)
        else:
            workers = min(max_workers, len(batches))
            with ThreadPoolExecutor(max_workers=workers) as executor:
                for batch_records in executor.map(_fetch_and_parse_batch, batches):
                    yield from batch_records

    @staticmethod
    def _parse_efetch(xml: str) -> Iterator[SourceMetadata]:
        parser = etree.XMLParser(recover=True, resolve_entities=False, no_network=True)
        root = etree.fromstring(xml.encode("utf-8"), parser=parser)
        for art in root.iter("PubmedArticle"):
            pmid_el = art.find(".//PMID")
            if pmid_el is None or not pmid_el.text:
                continue
            pmid = pmid_el.text.strip()

            title = _text(art.find(".//ArticleTitle"))
            journal = _text(art.find(".//Journal/Title"))
            iso = _text(art.find(".//Journal/ISOAbbreviation"))

            doi = pmcid = None
            for aid in art.iter("ArticleId"):
                t = aid.get("IdType")
                if t == "doi":
                    doi = (aid.text or "").strip()
                elif t == "pmc":
                    pmcid = (aid.text or "").strip()

            authors, affs = [], []
            for a in art.iter("Author"):
                ln, fn = _text(a.find("LastName")), _text(a.find("ForeName"))
                if ln:
                    authors.append(f"{fn} {ln}".strip())
                for aff in a.iter("Affiliation"):
                    if aff.text:
                        affs.append(aff.text.strip())

            ptypes = [(_text(p) or "").lower() for p in art.iter("PublicationType")]
            mesh = [_text(d) for d in art.iter("DescriptorName") if _text(d)]

            retracted = any("retracted publication" in p for p in ptypes)
            pub_date = _pubmed_date(art)

            meta = SourceMetadata(
                source_id=f"PMID:{pmid}",
                pmid=pmid, pmcid=pmcid, doi=doi,
                title=title.rstrip("."), journal=journal, journal_iso=iso,
                publication_date=pub_date,
                year=pub_date.year if pub_date else None,
                is_review=any("review" in p for p in ptypes),
                is_retracted=retracted,
                authors=authors, affiliations=list(dict.fromkeys(affs)),
                full_text_available=bool(pmcid),
            )
            # MeSH is carried out-of-band for the normalizer
            meta.__dict__["_mesh_terms"] = mesh
            yield meta


class BioRxivHarvester:
    """bioRxiv/medRxiv API. Supports date-range crawling and DOI lookup."""

    def __init__(self, session: PoliteSession, server: str = "biorxiv",
                 drift=None):
        self.s = session
        self.server = server
        if drift is None:
            from .contracts import DriftDetector
            drift = DriftDetector()
        self.drift = drift

    def fetch_range(
        self, start: str, end: str, max_results: int = 2000, query: str | None = None
    ) -> Iterator[SourceMetadata]:
        """start/end as 'YYYY-MM-DD'. Paginates in blocks of 100."""
        cursor = 0
        got = 0
        while got < max_results:
            url = f"{BIORXIV}/details/{self.server}/{start}/{end}/{cursor}"
            data = self.s.get_json(url)
            if not data:
                return
            self.drift.check("biorxiv_details", data)
            items = data.get("collection") or []
            if not items:
                return
            for it in items:
                if query:
                    target_text = f"{it.get('title', '')} {it.get('abstract', '')}".lower()
                    q_lower = query.lower().strip()
                    tokens = [t for t in q_lower.replace('"', '').split() if len(t) > 2]
                    if not tokens:
                        tokens = [t for t in q_lower.split() if t]
                    if not any(t in target_text for t in tokens):
                        continue
                meta = self._parse(it)
                if meta:
                    got += 1
                    yield meta
                    if got >= max_results:
                        return
            msgs = data.get("messages") or [{}]
            total = int(msgs[0].get("total", 0) or 0)
            cursor += len(items)
            if cursor >= total:
                return

    def _parse(self, it: dict) -> SourceMetadata | None:
        doi = it.get("doi")
        if not doi:
            return None
        d = None
        if it.get("date"):
            try:
                d = datetime.strptime(it["date"], "%Y-%m-%d").date()
            except ValueError:
                pass
        authors = [a.strip() for a in (it.get("authors") or "").split(";") if a.strip()]
        return SourceMetadata(
            source_id=f"DOI:{doi}",
            doi=doi,
            title=(it.get("title") or "").strip(),
            journal=f"{self.server} (preprint)",
            publication_date=d,
            year=d.year if d else None,
            is_preprint=True,
            authors=authors,
            affiliations=[it["author_corresponding_institution"]]
            if it.get("author_corresponding_institution") else [],
            corresponding_author=it.get("author_corresponding"),
            open_access=True,
            full_text_available=False,
        )


def _text(el) -> str:
    if el is None:
        return ""
    return " ".join(el.itertext()).strip()


def _pubmed_date(art) -> date | None:
    for path in (".//ArticleDate", ".//PubDate", ".//DateCompleted"):
        el = art.find(path)
        if el is None:
            continue
        y = _text(el.find("Year"))
        m = _text(el.find("Month")) or "1"
        d = _text(el.find("Day")) or "1"
        if not y:
            medline = _text(el.find("MedlineDate"))
            if medline and medline[:4].isdigit():
                y = medline[:4]
            else:
                continue
        months = {
            "jan": 1, "feb": 2, "mar": 3, "apr": 4, "may": 5, "jun": 6,
            "jul": 7, "aug": 8, "sep": 9, "oct": 10, "nov": 11, "dec": 12,
        }
        mi = months.get(m[:3].lower(), None)
        if mi is None:
            try:
                mi = int(m)
            except ValueError:
                mi = 1
        try:
            return date(int(y), max(1, min(12, mi)), max(1, min(28, int(d))))
        except (ValueError, TypeError):
            continue
    return None
