"""
Europe PMC harvester.

Europe PMC is the primary source because, unlike the PubMed E-utilities, it
serves full text for open-access articles through a REST endpoint. Full text
matters enormously here: the Methods section is where differentiation
protocol, DIV, marker antibodies, and replicate counts live, and none of that
appears in an abstract.

Endpoints used:
  /search              -- metadata + abstract, supports rich query syntax
  /{src}/{id}/fullTextXML -- JATS XML full text for OA articles
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import logging
import re
from datetime import date, datetime
from typing import Any, Iterator

from lxml import etree

from ..schema import SourceMetadata
from .http import PoliteSession

log = logging.getLogger(__name__)

BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest"

# Section headings mapped to canonical buckets. Real papers use dozens of
# variants; this covers the common ones and falls through to "other".
SECTION_MAP = {
    "abstract": "abstract",
    "introduction": "introduction",
    "background": "introduction",
    "methods": "methods",
    "materials and methods": "methods",
    "experimental procedures": "methods",
    "star methods": "methods",
    "methods and materials": "methods",
    "results": "results",
    "results and discussion": "results",
    "discussion": "discussion",
    "conclusion": "discussion",
    "conclusions": "discussion",
    "supplementary material": "supplementary",
}


class EuropePMCHarvester:
    def __init__(self, session: PoliteSession, drift: Any = None):
        self.s = session
        # Validates the first response of the session against the endpoint
        # contract. Costs one check per run; catches an upstream field rename
        # before it silently empties thirty thousand records.
        if drift is None:
            from .contracts import DriftDetector
            drift = DriftDetector()
        self.drift = drift

    # ---------------- search ----------------

    def search(
        self,
        query: str,
        page_size: int = 100,
        max_results: int = 1000,
        open_access_only: bool = False,
    ) -> Iterator[SourceMetadata]:
        """Paginate through Europe PMC search results via cursorMark."""
        q = query
        if open_access_only:
            q = f"({query}) AND OPEN_ACCESS:y"

        cursor = "*"
        fetched = 0
        while fetched < max_results:
            params = {
                "query": q,
                "format": "json",
                "pageSize": min(page_size, max_results - fetched),
                "cursorMark": cursor,
                "resultType": "core",
            }
            data = self.s.get_json(f"{BASE}/search", params)
            if not data:
                log.warning("Search returned nothing for query: %s", q[:80])
                return
            self.drift.check("europepmc_envelope", data)
            self.drift.check("europepmc_search", data)
            results = (data.get("resultList") or {}).get("result") or []
            if not results:
                return
            for r in results:
                meta = self._parse_result(r)
                if meta:
                    fetched += 1
                    yield meta
            next_cursor = data.get("nextCursorMark")
            if not next_cursor or next_cursor == cursor:
                return
            cursor = next_cursor

    def count(self, query: str) -> int:
        data = self.s.get_json(
            f"{BASE}/search",
            {"query": query, "format": "json", "pageSize": 1, "resultType": "idlist"},
        )
        return int(data.get("hitCount", 0)) if data else 0

    # ---------------- parsing ----------------

    def _parse_result(self, r: dict) -> SourceMetadata | None:
        pmid = r.get("pmid")
        pmcid = r.get("pmcid")
        doi = r.get("doi")
        source_id = (
            f"PMID:{pmid}" if pmid else f"PMC:{pmcid}" if pmcid
            else f"DOI:{doi}" if doi else None
        )
        if not source_id:
            return None

        journal = None
        journal_iso = None
        ji = r.get("journalInfo") or {}
        j = ji.get("journal") or {}
        journal = j.get("title") or r.get("journalTitle")
        journal_iso = j.get("isoabbreviation") or j.get("medlineAbbreviation")

        pub_date = self._parse_date(r)
        authors = self._parse_authors(r)
        pub_types = [p.lower() for p in (r.get("pubTypeList") or {}).get("pubType", [])]

        return SourceMetadata(
            source_id=source_id,
            pmid=str(pmid) if pmid else None,
            pmcid=pmcid,
            doi=doi,
            title=(r.get("title") or "").strip().rstrip("."),
            journal=journal,
            journal_iso=journal_iso,
            publication_date=pub_date,
            year=pub_date.year if pub_date else self._year_only(r),
            is_preprint=r.get("source") == "PPR" or "preprint" in pub_types,
            is_review="review" in pub_types,
            authors=authors,
            affiliations=self._parse_affiliations(r),
            citation_count=r.get("citedByCount"),
            open_access=r.get("isOpenAccess") == "Y",
            full_text_available=bool(r.get("hasTextMinedTerms") == "Y" or pmcid),
        )

    @staticmethod
    def _year_only(r: dict) -> int | None:
        y = (r.get("journalInfo") or {}).get("yearOfPublication") or r.get("pubYear")
        try:
            return int(y)
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _parse_date(r: dict) -> date | None:
        for field in ("firstPublicationDate", "electronicPublicationDate"):
            v = r.get(field)
            if v:
                for fmt in ("%Y-%m-%d", "%Y-%m", "%Y"):
                    try:
                        return datetime.strptime(v, fmt).date()
                    except ValueError:
                        continue
        return None

    @staticmethod
    def _parse_authors(r: dict) -> list[str]:
        al = (r.get("authorList") or {}).get("author") or []
        out = []
        for a in al:
            name = a.get("fullName") or a.get("lastName") or ""
            if name:
                out.append(name)
        if not out and r.get("authorString"):
            out = [x.strip() for x in r["authorString"].split(",") if x.strip()]
        return out

    @staticmethod
    def _parse_affiliations(r: dict) -> list[str]:
        al = (r.get("authorList") or {}).get("author") or []
        affs: list[str] = []
        for a in al:
            aff = a.get("affiliation")
            if aff:
                affs.append(aff)
            for x in (a.get("authorAffiliationDetailsList") or {}).get(
                "authorAffiliation", []
            ):
                if x.get("affiliation"):
                    affs.append(x["affiliation"])
        seen, out = set(), []
        for a in affs:
            if a not in seen:
                seen.add(a)
                out.append(a)
        return out

    # ---------------- full text ----------------

    def fetch_full_text(self, meta: SourceMetadata) -> dict[str, str] | None:
        """Fetch and section JATS full text. Returns {section: text} or None."""
        if not meta.pmcid:
            return None
        pmcid = meta.pmcid if meta.pmcid.startswith("PMC") else f"PMC{meta.pmcid}"
        xml = self.s.get_text(f"{BASE}/{pmcid}/fullTextXML")
        if not xml or "<" not in xml:
            return None
        try:
            return self.parse_jats(xml)
        except etree.XMLSyntaxError as e:
            log.warning("XML parse failed for %s: %s", pmcid, e)
            return None

    def fetch_full_texts_concurrent(
        self,
        metas: list[SourceMetadata],
        max_workers: int = 6,
    ) -> Iterator[tuple[SourceMetadata, dict[str, str] | None]]:
        """Fetch and parse JATS full text for multiple papers concurrently.

        Submits full-text XML downloads and parsing to a ThreadPoolExecutor,
        yielding (meta, sections) as tasks complete. Papers without a PMCID
        yield (meta, None) immediately without thread pool submission.
        """
        if not metas:
            return

        to_fetch: list[SourceMetadata] = []
        for m in metas:
            if m.pmcid:
                to_fetch.append(m)
            else:
                yield (m, None)

        if not to_fetch:
            return

        if max_workers <= 1 or len(to_fetch) <= 1:
            for m in to_fetch:
                yield (m, self.fetch_full_text(m))
            return

        workers = min(max_workers, len(to_fetch))
        with ThreadPoolExecutor(max_workers=workers) as executor:
            future_to_meta = {
                executor.submit(self.fetch_full_text, meta): meta
                for meta in to_fetch
            }
            for future in as_completed(future_to_meta):
                meta = future_to_meta[future]
                try:
                    sections = future.result()
                except Exception as e:
                    log.warning("Concurrent full-text fetch failed for %s: %s", meta.pmcid, e)
                    sections = None
                yield (meta, sections)

    @staticmethod
    def parse_jats(xml: str) -> dict[str, str]:
        """Extract canonical sections from JATS XML.

        Also pulls figure captions, which are disproportionately information
        dense for quantitative phenotype claims -- effect sizes and n's often
        appear only in the legend.
        """
        parser = etree.XMLParser(recover=True, resolve_entities=False, no_network=True)
        root = etree.fromstring(xml.encode("utf-8"), parser=parser)
        sections: dict[str, list[str]] = {}

        def add(name: str, text: str) -> None:
            if text and text.strip():
                sections.setdefault(name, []).append(text.strip())

        for abst in root.iter("abstract"):
            add("abstract", _itertext(abst))

        body = root.find(".//body")
        if body is not None:
            for sec in body.findall("sec"):
                title_el = sec.find("title")
                raw_title = (
                    _itertext(title_el).strip().lower() if title_el is not None else ""
                )
                canon = SECTION_MAP.get(raw_title, "")
                if not canon:
                    for k, v in SECTION_MAP.items():
                        if k in raw_title:
                            canon = v
                            break
                add(canon or "other", _itertext(sec))

        for fig in root.iter("fig"):
            cap = fig.find(".//caption")
            if cap is not None:
                label_el = fig.find("label")
                label = _itertext(label_el) if label_el is not None else ""
                add("figures", f"{label} {_itertext(cap)}".strip())

        for tbl in root.iter("table-wrap"):
            cap = tbl.find(".//caption")
            if cap is not None:
                add("tables", _itertext(cap))

        return {k: "\n\n".join(v) for k, v in sections.items()}


def _itertext(el) -> str:
    """Flatten element text, normalizing whitespace."""
    if el is None:
        return ""
    text = " ".join(el.itertext())
    return re.sub(r"\s+", " ", text).strip()
