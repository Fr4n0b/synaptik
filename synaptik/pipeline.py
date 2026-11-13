"""
End-to-end pipeline orchestration.

Ties harvesting, extraction, normalization, and indexing into resumable
operations. Resumability matters: a full harvest of the neuronal + splicing
literature is tens of thousands of papers and many hours, and it WILL be
interrupted. Every stage checks the store before doing work, so re-running
after a crash costs only the incomplete items.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import logging
import time
from dataclasses import dataclass, field
from datetime import date, datetime, timezone
from typing import Callable, Iterable

from .extract.llm import LLMExtractor
from .ingest.contracts import DriftDetector
from .ingest.europepmc import EuropePMCHarvester
from .ingest.http import PoliteSession
from .ingest.pubmed import BioRxivHarvester, PubMedHarvester
from .infer.weighting import JournalQualityResolver
from .normalize.labs import LabResolver
from .normalize.matcher import PhenotypeMatcher, classify_assay
from .schema import ExperimentRecord, SourceMetadata
from .store import Store

log = logging.getLogger(__name__)


# --------------------------------------------------------------------------
# Corpus definition
# --------------------------------------------------------------------------

# Query set defining the corpus. Split into focused queries rather than one
# giant boolean because Europe PMC relevance-ranks per query, and a single
# broad query buries the specific work under reviews.
CORPUS_QUERIES: dict[str, str] = {
    "ipsc_neuron_core": (
        '("iPSC-derived neuron*" OR "iPSC neuron*" OR "induced pluripotent stem cell'
        '-derived neuron*" OR "hiPSC-derived neuron*")'
    ),
    "ngn2": (
        '(NGN2 OR NEUROG2 OR "Neurogenin-2" OR "Neurogenin 2") AND '
        '(neuron* OR "induced neuron*" OR iNeuron*)'
    ),
    "synapse_phenotype": (
        '(synaps* OR "synaptic puncta" OR PSD95 OR "PSD-95" OR synapsin OR '
        'SYN1 OR homer1 OR bassoon OR VGLUT1) AND '
        '(iPSC OR "stem cell-derived" OR "induced neuron*")'
    ),
    "splicing_neuro": (
        '("alternative splicing" OR "mis-splicing" OR "cryptic exon" OR '
        '"intron retention" OR "exon skipping" OR "splicing factor") AND '
        '(neuron* OR neuronal OR neurodegener* OR ALS OR FTD OR "spinal muscular '
        'atrophy" OR autism)'
    ),
    "rbp_neuro": (
        '(TDP-43 OR TARDBP OR FUS OR hnRNP OR RBFOX1 OR NOVA1 OR NOVA2 OR PTBP1 '
        'OR SRSF1 OR MATR3 OR SFPQ) AND (splicing OR neuron*)'
    ),
    "patient_mutation_ipsc": (
        '("patient-derived" OR "patient derived" OR isogenic OR CRISPR) AND '
        'iPSC AND (neuron* OR neural) AND (mutation OR variant)'
    ),
    "electrophysiology_ipsc": (
        '(mEPSC OR "miniature excitatory" OR "multielectrode array" OR MEA OR '
        '"patch clamp" OR "network burst") AND (iPSC OR "induced neuron*")'
    ),
    "drug_rescue_neuron": (
        '(rescue OR reversal OR "restored" OR repurpos*) AND '
        '(iPSC OR "induced neuron*") AND (compound OR drug OR "small molecule" '
        'OR inhibitor OR agonist)'
    ),
}


@dataclass
class PipelineStats:
    sources_seen: int = 0
    sources_new: int = 0
    fulltext_fetched: int = 0
    fulltext_missing: int = 0
    papers_extracted: int = 0
    records_created: int = 0
    records_rejected_quote: int = 0
    records_rejected_schema: int = 0
    extraction_failures: int = 0
    http_failures: int = 0
    elapsed_seconds: float = 0.0
    unmapped_phenotypes: list[tuple[str, int]] = field(default_factory=list)

    def failed(self) -> bool:
        """True when the stage did no useful work. Nothing is not success."""
        did_nothing = (
            self.sources_seen == 0
            and self.fulltext_fetched == 0
            and self.records_created == 0
        )
        return did_nothing or (self.http_failures > 0 and self.sources_seen == 0)

    def summary(self) -> str:
        if self.failed():
            return (
                f"NO DATA RETRIEVED. {self.http_failures} HTTP failure(s) in "
                f"{self.elapsed_seconds:.0f}s.\n"
                "This is a failure, not an empty result set. Check network "
                "access, then run `synaptik preflight`."
            )
        rate = (
            self.records_rejected_quote
            / max(1, self.records_created + self.records_rejected_quote)
        )
        parts = [
            f"sources: {self.sources_new} new / {self.sources_seen} seen",
            f"full text: {self.fulltext_fetched} ok, "
            f"{self.fulltext_missing} unavailable",
            f"records: {self.records_created} kept, "
            f"{self.records_rejected_quote} rejected as unsupported "
            f"({rate:.1%} hallucination rate)",
            f"{self.elapsed_seconds:.0f}s",
        ]
        out = " | ".join(parts)
        if self.http_failures:
            out += (
                f"\nWARNING: {self.http_failures} HTTP request(s) failed; "
                "this harvest is incomplete. Re-run to fill the gaps."
            )
        return out


class Pipeline:
    def __init__(
        self,
        store: Store,
        session: PoliteSession | None = None,
        extractor: LLMExtractor | None = None,
        contact_email: str | None = None,
        ncbi_api_key: str | None = None,
        scimago_csv: str | None = None,
        strict_contracts: bool = True,
    ):
        self.store = store
        self.session = session or PoliteSession(
            store=store, contact_email=contact_email
        )
        # Strict by default: a long unattended harvest that silently drops a
        # field is far more expensive than one that stops on the first page.
        self.drift = DriftDetector(strict=strict_contracts)
        self.epmc = EuropePMCHarvester(self.session, drift=self.drift)
        self.pubmed = PubMedHarvester(
            self.session, api_key=ncbi_api_key, drift=self.drift
        )
        self.biorxiv = BioRxivHarvester(self.session, server="biorxiv", drift=self.drift)
        self.medrxiv = BioRxivHarvester(self.session, server="medrxiv", drift=self.drift)
        self.extractor = extractor
        self.journals = JournalQualityResolver(scimago_csv)
        self.labs = LabResolver()
        self.matcher = PhenotypeMatcher()

    # ---------------- stage 0: preflight ----------------

    def preflight(
        self, progress: Callable[[str], None] | None = None
    ) -> tuple[bool, str]:
        """Fetch one small page from each API and validate against contracts.

        Run this before committing to a multi-hour harvest. It costs three
        requests and tells you whether the parsers still match reality.
        """
        say = progress or (lambda m: log.info(m))
        say("preflight: checking API response contracts...")

        probe = DriftDetector(strict=False)
        results = []
        n_epmc = 0
        failures_before = self.session.n_failures

        # Europe PMC: smallest possible real query
        epmc = EuropePMCHarvester(self.session, drift=probe)
        try:
            got = list(epmc.search("NGN2 neuron", max_results=3))
            n_epmc = len(got)
            say(f"  Europe PMC: {len(got)} records returned")
            if got and not got[0].title:
                results.append("Europe PMC returned records with empty titles.")
        except Exception as e:
            results.append(f"Europe PMC probe failed: {e}")

        # PubMed
        n_pubmed = 0
        try:
            pmids = self.pubmed.search_pmids("NGN2 neuron", max_results=3)
            n_pubmed = len(pmids)
            say(f"  PubMed: {len(pmids)} PMIDs returned")
        except Exception as e:
            results.append(f"PubMed probe failed: {e}")

        # An empty response is not a passing contract check -- it is the
        # absence of a check. Zero records means the contracts were never
        # exercised, which is exactly how a fully-blocked network produced a
        # "successful" harvest of nothing.
        if n_epmc == 0 and n_pubmed == 0:
            results.append(
                "No records returned from any API. Either the network is "
                "blocked, the endpoints are down, or the queries match "
                "nothing. Contracts could not be checked at all."
            )
        elif n_epmc == 0:
            results.append(
                "Europe PMC returned no records. It is the primary source; "
                "harvesting without it will produce a corpus with no full text."
            )

        new_failures = self.session.n_failures - failures_before
        if new_failures:
            results.append(
                f"{new_failures} HTTP request(s) failed during preflight. "
                "Check network access and any egress restrictions."
            )

        ok = probe.all_ok() and not results
        report = probe.summary()
        if results:
            report += "\n\nAdditional problems:\n" + "\n".join(
                f"  - {r}" for r in results
            )
        if ok:
            report += "\n\nAll contracts satisfied. Safe to harvest."
        else:
            report += (
                "\n\nDO NOT HARVEST until these are resolved -- records would "
                "be written with missing fields and no error raised.\n"
                "Fix the parsers in synaptik/ingest/, then re-run preflight."
            )
        return (ok, report)

    # ---------------- stage 1: harvest metadata ----------------

    def harvest(
        self,
        queries: dict[str, str] | None = None,
        max_per_query: int = 2000,
        concurrency: int = 5,
        include_preprints: bool = True,
        progress: Callable[[str], None] | None = None,
    ) -> PipelineStats:
        """Fetch and store source metadata for the corpus."""
        st = PipelineStats()
        t0 = time.time()
        queries = queries or CORPUS_QUERIES
        say = progress or (lambda m: log.info(m))
        http_failures_before = self.session.n_failures

        collected: list[SourceMetadata] = []

        def _flush_batch() -> None:
            nonlocal collected
            if collected:
                self._flush_sources(collected)
                collected = []

        def _harvest_single_query(name: str, q: str) -> tuple[str, list[SourceMetadata]]:
            results = []
            for meta in self.epmc.search(q, max_results=max_per_query):
                results.append(meta)
            return (name, results)

        if concurrency <= 1 or len(queries) <= 1:
            for name, q in queries.items():
                say(f"harvesting '{name}'...")
                n_before = st.sources_new
                for meta in self.epmc.search(q, max_results=max_per_query):
                    st.sources_seen += 1
                    if self.store.has_source(meta.source_id):
                        continue
                    if meta.is_preprint and not include_preprints:
                        continue
                    collected.append(meta)
                    st.sources_new += 1
                    if len(collected) >= 500:
                        _flush_batch()
                say(f"  '{name}': +{st.sources_new - n_before} new")
            if include_preprints:
                start_date = "2020-01-01"
                end_date = date.today().isoformat()
                for srv_name, harvester in [("biorxiv", self.biorxiv), ("medrxiv", self.medrxiv)]:
                    say(f"harvesting {srv_name} preprints...")
                    for name, q in queries.items():
                        n_before = st.sources_new
                        for meta in harvester.fetch_range(start_date, end_date, max_results=max_per_query, query=q):
                            st.sources_seen += 1
                            if self.store.has_source(meta.source_id):
                                continue
                            collected.append(meta)
                            st.sources_new += 1
                            if len(collected) >= 500:
                                _flush_batch()
                        say(f"  '{srv_name}:{name}': +{st.sources_new - n_before} new")
        else:
            say(f"harvesting {len(queries)} query sets (concurrency={concurrency})...")
            tasks = {}
            with ThreadPoolExecutor(max_workers=concurrency) as executor:
                for name, q in queries.items():
                    f = executor.submit(_harvest_single_query, name, q)
                    tasks[f] = f"epmc:{name}"
                if include_preprints:
                    start_date = "2020-01-01"
                    end_date = date.today().isoformat()
                    for srv_name, harvester in [("biorxiv", self.biorxiv), ("medrxiv", self.medrxiv)]:
                        for name, q in queries.items():
                            def _make_fn(h=harvester, q_str=q):
                                return lambda: list(h.fetch_range(start_date, end_date, max_results=max_per_query, query=q_str))
                            f = executor.submit(_make_fn(harvester, q))
                            tasks[f] = f"{srv_name}:{name}"

                for future in as_completed(tasks):
                    task_name = tasks[future]
                    try:
                        res = future.result()
                        metas = res[1] if isinstance(res, tuple) else res
                    except Exception as e:
                        log.warning("Harvest task failed for %s: %s", task_name, e)
                        metas = []
                    n_before = st.sources_new
                    for meta in metas:
                        st.sources_seen += 1
                        if self.store.has_source(meta.source_id):
                            continue
                        if meta.is_preprint and not include_preprints:
                            continue
                        collected.append(meta)
                        st.sources_new += 1
                        if len(collected) >= 500:
                            _flush_batch()
                    say(f"  '{task_name}': +{st.sources_new - n_before} new (seen {len(metas)})")

        _flush_batch()

        st.http_failures = self.session.n_failures - http_failures_before
        st.elapsed_seconds = time.time() - t0
        if st.failed():
            return st

        # Lab clustering needs the whole corpus, so it runs after harvesting.
        self._recluster_labs(say)
        st.elapsed_seconds = time.time() - t0
        return st

    def _flush_sources(self, metas: list[SourceMetadata]) -> None:
        with self.store.tx():
            for m in metas:
                self.journals.annotate(m)
                self.store.upsert_source(m)
                self.store.index_source_text(m.source_id, m.title, "")

    def _recluster_labs(self, say: Callable[[str], None]) -> None:
        """Re-run lab clustering over all known sources and persist the IDs."""
        say("clustering research groups...")
        rows = self.store.conn.execute("SELECT payload FROM sources").fetchall()
        sources = [SourceMetadata.model_validate_json(r["payload"]) for r in rows]
        self.labs = LabResolver()
        self.labs.fit(sources)
        with self.store.tx():
            for s in sources:
                s.lab_id = self.labs.assign(s)
                self.store.upsert_source(s)
        n_labs = len({s.lab_id for s in sources})
        say(f"  {len(sources)} sources -> {n_labs} distinct research groups")

    # ---------------- stage 2: fetch full text ----------------

    def fetch_full_texts(
        self,
        limit: int | None = None,
        max_workers: int = 6,
        progress: Callable[[str], None] | None = None,
    ) -> PipelineStats:
        st = PipelineStats()
        t0 = time.time()
        say = progress or (lambda m: log.info(m))
        http_failures_before = self.session.n_failures

        rows = self.store.conn.execute(
            """SELECT s.payload FROM sources s
               LEFT JOIN documents d ON d.source_id = s.source_id
               WHERE d.source_id IS NULL AND s.pmcid IS NOT NULL
               ORDER BY s.year DESC"""
        ).fetchall()
        if limit:
            rows = rows[:limit]

        if not rows:
            say("no papers pending full-text fetch.")
            st.elapsed_seconds = time.time() - t0
            return st

        say(f"fetching full text for {len(rows)} papers (max_workers={max_workers})...")
        metas = [SourceMetadata.model_validate_json(r["payload"]) for r in rows]

        count = 0
        for meta, sections in self.epmc.fetch_full_texts_concurrent(metas, max_workers=max_workers):
            count += 1
            if not sections:
                st.fulltext_missing += 1
                continue
            self.store.store_document(
                meta.source_id,
                abstract=sections.get("abstract", ""),
                full_text="\n\n".join(sections.values()),
                methods=sections.get("methods", ""),
                results=sections.get("results", ""),
                sections=sections,
            )
            self.store.index_source_text(
                meta.source_id, meta.title, sections.get("abstract", "")
            )
            st.fulltext_fetched += 1
            if count % 50 == 0 or count == len(rows):
                self.store.conn.commit()
                say(f"  {count}/{len(rows)} ({st.fulltext_fetched} retrieved)")

        self.store.conn.commit()
        st.http_failures = self.session.n_failures - http_failures_before
        st.elapsed_seconds = time.time() - t0
        return st

    # ---------------- stage 3: extract records ----------------

    def extract_records(
        self,
        limit: int | None = None,
        skip_reviews: bool = True,
        progress: Callable[[str], None] | None = None,
    ) -> PipelineStats:
        """Run LLM extraction over documents that have no records yet."""
        st = PipelineStats()
        t0 = time.time()
        say = progress or (lambda m: log.info(m))

        if not self.extractor:
            raise RuntimeError(
                "No extractor configured. Construct Pipeline with "
                "extractor=LLMExtractor(...) to run extraction."
            )

        sql = """SELECT s.payload AS src, d.sections_json AS sections
                 FROM sources s JOIN documents d ON d.source_id = s.source_id
                 WHERE s.source_id NOT IN (SELECT DISTINCT source_id FROM records)"""
        if skip_reviews:
            sql += " AND s.is_review = 0"
        sql += " ORDER BY s.year DESC"

        rows = self.store.conn.execute(sql).fetchall()
        if limit:
            rows = rows[:limit]
        say(f"extracting from {len(rows)} papers...")

        import json as _json

        for i, row in enumerate(rows, 1):
            meta = SourceMetadata.model_validate_json(row["src"])
            try:
                sections = _json.loads(row["sections"] or "{}")
            except _json.JSONDecodeError:
                sections = {}
            if not sections:
                continue

            result = self.extractor.extract(meta, sections)
            st.papers_extracted += 1
            st.records_rejected_quote += result.n_quote_failed
            st.records_rejected_schema += result.n_schema_failed
            if not result.records and result.n_proposed == 0:
                st.extraction_failures += 1

            normalized = [self.normalize_record(r) for r in result.records]
            if normalized:
                self.store.insert_records(normalized)
                st.records_created += len(normalized)

            if i % 20 == 0:
                say(f"  {i}/{len(rows)} papers, {st.records_created} records")

        st.unmapped_phenotypes = self.matcher.unmapped_report(30)
        st.elapsed_seconds = time.time() - t0
        return st

    # ---------------- normalization ----------------

    def normalize_record(self, rec: ExperimentRecord) -> ExperimentRecord:
        """Apply ontology mapping, journal tier, and lab ID to a raw record."""
        rec.outcome.phenotype = self.matcher.match(rec.outcome.phenotype.raw)

        if rec.assay_raw:
            cls, conf = classify_assay(rec.assay_raw)
            if conf > 0.5:
                rec.assay_class = cls

        self.journals.annotate(rec.source)
        if not rec.source.lab_id:
            rec.source.lab_id = self.labs.assign(rec.source)

        # Pathway terms are kept as raw text plus a normalized key. We do not
        # force them into an ontology here -- pathway naming in papers is far
        # too heterogeneous, and a wrong mapping silently merges distinct
        # biology. Grouping happens on the normalized key in the store.
        for p in rec.pathways:
            if not p.ontology_id:
                p.ontology_label = _tidy_pathway(p.raw)

        rec.ingested_at = datetime.now(timezone.utc).isoformat()
        rec.record_id = rec.compute_id()
        return rec

    # ---------------- convenience ----------------

    def run_all(
        self,
        max_per_query: int = 500,
        fulltext_limit: int | None = None,
        extract_limit: int | None = None,
        progress: Callable[[str], None] | None = None,
    ) -> dict[str, PipelineStats]:
        say = progress or (lambda m: log.info(m))
        out = {}
        say("=== stage 1: harvest ===")
        out["harvest"] = self.harvest(max_per_query=max_per_query, progress=say)
        say(out["harvest"].summary())
        say("=== stage 2: full text ===")
        out["fulltext"] = self.fetch_full_texts(limit=fulltext_limit, progress=say)
        say(out["fulltext"].summary())
        if self.extractor:
            say("=== stage 3: extraction ===")
            out["extract"] = self.extract_records(limit=extract_limit, progress=say)
            say(out["extract"].summary())
        else:
            say("(skipping extraction: no extractor configured)")
        return out


def _tidy_pathway(raw: str) -> str:
    """Light canonicalization of pathway names for grouping.

    Deliberately conservative: lowercase, strip punctuation and the word
    'pathway'/'signaling' when trailing. Anything more aggressive risks merging
    genuinely different pathways.
    """
    import re

    s = raw.strip().lower()
    s = re.sub(r"[^\w\s/-]", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    s = re.sub(r"\s+(pathway|signalling|signaling|cascade)$", "", s)
    return s.strip()
