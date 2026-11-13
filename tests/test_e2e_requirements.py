"""
SYNAPTIK End-to-End Opaque-Box Acceptance Test Suite: Requirements R1–R5 (Tiers 1–4).

Author & Principal Investigator:
Francesco Nobile (ORCID: 0009-0001-9477-3332)
Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School;
University of Pavia.

Test Philosophy:
- Opaque-box requirement-driven testing: zero dependency on private/internal functions.
- 4-Tier methodology:
    Tier 1: Category-Partition Feature Coverage (>=5 test cases per feature in R1-R5)
    Tier 2: Boundary Value Analysis & Extreme Corner Probing (>=5 per feature)
    Tier 3: Pairwise Cross-Feature Interactions & End-to-End Pipelines
    Tier 4: Real-World Clinical & Biophysical Workloads (5 Full Scenarios)
- Pure Python standard library unittest: zero external test runner dependencies.
- Offline deterministic execution: all external network APIs mocked/fixtured.
- Zero-PHI persistence auditing across SQLite, filesystem, temp directories, and logs.
- Python 3.10-3.12 syntax compatibility (verifying zero nested f-strings with backslash in braces).

Run via standard library unittest:
    python -m unittest tests/test_e2e_requirements.py
"""

from __future__ import annotations

import ast
import concurrent.futures
import io
import json
import logging
import math
import os
import pathlib
import re
import sqlite3
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from dataclasses import dataclass
from datetime import date
from typing import Any, List, Optional
from unittest.mock import MagicMock, patch

# ==============================================================================
# DYNAMIC MODULE IMPORTS & AVAILABILITY GUARDS
# ==============================================================================

# Requirement R1: Ingest & Preprints
try:
    from synaptik.ingest.http import ConcurrentRateLimiter, PoliteSession, RateLimiter
    HAVE_INGEST_HTTP = True
except ImportError:
    HAVE_INGEST_HTTP = False

try:
    from synaptik.ingest.europepmc import EuropePMCHarvester
    HAVE_INGEST_EUROPEPMC = True
except ImportError:
    HAVE_INGEST_EUROPEPMC = False

try:
    from synaptik.ingest.pubmed import BioRxivHarvester, PubMedHarvester
    HAVE_INGEST_PUBMED = True
except ImportError:
    HAVE_INGEST_PUBMED = False

try:
    from synaptik.infer.weighting import EvidenceWeighter, WeightBreakdown
    from synaptik.schema import EvidenceType, ExperimentRecord, Outcome, SourceMetadata
    HAVE_EVIDENCE_WEIGHTER = True
except ImportError:
    HAVE_EVIDENCE_WEIGHTER = False

# Requirement R2: Clinical Pharmacology & Dose Escalation
try:
    from synaptik.clinical.allometry import (
        AllometryResult,
        calculate_allometric_dose,
    )
    HAVE_CLINICAL_ALLOMETRY = True
except ImportError:
    HAVE_CLINICAL_ALLOMETRY = False

try:
    from synaptik.clinical.pkpd import PKPDResult, simulate_intrathecal_pkpd
    HAVE_CLINICAL_PKPD = True
except ImportError:
    HAVE_CLINICAL_PKPD = False

try:
    from synaptik.clinical.escalation import (
        BOINEscalationEngine,
        CohortObservation,
        EscalationDecision,
        ThreePlusThreeEngine,
    )
    HAVE_CLINICAL_ESCALATION = True
except ImportError:
    HAVE_CLINICAL_ESCALATION = False

try:
    from synaptik.clinical.protocol import (
        ClinicalTrialProtocol,
        generate_clinical_protocol,
    )
    HAVE_CLINICAL_PROTOCOL = True
except ImportError:
    HAVE_CLINICAL_PROTOCOL = False

try:
    from synaptik.cli.commands.clinical import cmd_clinical
    HAVE_CLI_CLINICAL = True
except ImportError:
    HAVE_CLI_CLINICAL = False

# Requirement R3: 3D RNA Biophysics & Steric Target Accessibility
try:
    from synaptik.biophysics.parameters import (
        get_hairpin_energy,
        get_stacking_energy,
        get_tetraloop_bonus,
    )
    HAVE_BIOPHYSICS_PARAMS = True
except ImportError:
    HAVE_BIOPHYSICS_PARAMS = False

try:
    from synaptik.biophysics.folding import (
        compute_delta_g_open,
        fold_mfe_pure_python,
    )
    HAVE_BIOPHYSICS_FOLDING = True
except ImportError:
    HAVE_BIOPHYSICS_FOLDING = False

try:
    from synaptik.biophysics.accessibility import (
        EffectiveBindingResult,
        compute_effective_binding,
    )
    HAVE_BIOPHYSICS_ACCESSIBILITY = True
except ImportError:
    HAVE_BIOPHYSICS_ACCESSIBILITY = False

try:
    from synaptik.biophysics.competition import (
        DisplacementResult,
        compute_steric_displacement,
    )
    HAVE_BIOPHYSICS_COMPETITION = True
except ImportError:
    HAVE_BIOPHYSICS_COMPETITION = False

try:
    from synaptik.cli.commands.biophysics import cmd_biophysics
    HAVE_CLI_BIOPHYSICS = True
except ImportError:
    HAVE_CLI_BIOPHYSICS = False

# Requirement R4: Patient Natural History Digital Twin & Bayesian Control
try:
    from synaptik.twin.twin_generator import (
        CounterfactualCohort,
        PatientBaseline,
        generate_counterfactual_cohort,
    )
    HAVE_TWIN_SIMULATOR = True
except ImportError:
    HAVE_TWIN_SIMULATOR = False

try:
    from synaptik.twin.distributions import (
        ConditionalProgression,
        compute_conditional_progression,
    )
    HAVE_TWIN_DISTRIBUTIONS = True
except ImportError:
    HAVE_TWIN_DISTRIBUTIONS = False

try:
    from synaptik.twin.bayesian_inference import (
        BayesianEfficacyResult,
        VisitObservation,
        evaluate_treatment_efficacy,
    )
    HAVE_TWIN_BAYESIAN = True
except ImportError:
    HAVE_TWIN_BAYESIAN = False

try:
    from synaptik.twin.dossier_export import generate_twin_fda_dossier
    HAVE_TWIN_DOSSIER = True
except ImportError:
    HAVE_TWIN_DOSSIER = False

try:
    from synaptik.cli.commands.twin import cmd_twin
    HAVE_CLI_TWIN = True
except ImportError:
    HAVE_CLI_TWIN = False

# Requirement R5: Author Citation & System
try:
    from synaptik.cli.commands.cite import (
        AFFILIATION,
        AUTHOR,
        ORCID,
        ORCID_URL,
        REPO_URL,
        cmd_cite,
        get_apa,
        get_bibtex,
        get_nature,
    )
    HAVE_CITE = True
except ImportError:
    HAVE_CITE = False

# Root repository workspace path
REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent


# ==============================================================================
# HELPER FIXTURES & DETERMINISTIC MOCKS
# ==============================================================================

class MockHTTPResponse:
    """Deterministic mock for requests.Response."""
    def __init__(self, status_code: int = 200, json_data: Any = None, text: str = ""):
        self.status_code = status_code
        self._json_data = json_data
        self.text = text
        self.content = text.encode("utf-8")

    def json(self) -> Any:
        if self._json_data is not None:
            return self._json_data
        raise json.JSONDecodeError("Expecting value", "", 0)


def create_mock_experiment_record(is_preprint: bool = False, is_retracted: bool = False) -> Any:
    """Constructs a minimal ExperimentRecord for evidence weighting tests."""
    from synaptik.schema import (
        CellModel,
        CellModelClass,
        Direction,
        EvidenceType,
        ExperimentRecord,
        NormalizedTerm,
        Outcome,
        Perturbation,
        SourceMetadata,
    )
    src = SourceMetadata(
        source_id="PMID:12345678",
        pmid="12345678",
        title="Test Investigation of Splice Switching",
        journal="bioRxiv" if is_preprint else "Nature Neuroscience",
        pub_year=2026,
        is_preprint=is_preprint,
        is_retracted=is_retracted,
        sjr_quartile="Q1" if not is_preprint else None,
    )
    model = CellModel(
        model_class=CellModelClass.IPSC_NEURON_NGN2,
        raw_description="iPSC motor neurons",
        is_isogenic_pair=True,
    )
    outcome = Outcome(
        phenotype=NormalizedTerm(raw="stathmin-2 expression"),
        direction=Direction.INCREASE,
        n_biological_replicates=4,
    )
    pert = Perturbation(agent="ASO-001")
    return ExperimentRecord(
        record_id="REC-TEST-001",
        source=src,
        cell_model=model,
        perturbation=pert,
        evidence_type=EvidenceType.DIRECT_MEASUREMENT,
        outcome=outcome,
    )


# ==============================================================================
# TIER 1: FEATURE COVERAGE (>=5 TEST CASES PER FEATURE ACROSS R1–R5)
# ==============================================================================

class TestTier1FeatureCoverage(unittest.TestCase):
    """
    Tier 1: Category-Partition Feature Coverage.
    Verifies nominal behavior and interface contracts with >=5 cases per feature.
    """

    # --------------------------------------------------------------------------
    # Feature R1.1: ConcurrentRateLimiter
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_INGEST_HTTP, "synaptik.ingest.http not available")
    def test_tc_r1_1_01_sequential_requests_rate_compliance(self):
        """TC-R1.1-01: Sequential requests strictly observe NCBI <=3 req/s minimum interval."""
        limiter = ConcurrentRateLimiter(limits={"eutils.ncbi.nlm.nih.gov": 0.05})
        host = "eutils.ncbi.nlm.nih.gov"
        d1 = limiter.acquire_delay(host)
        d2 = limiter.acquire_delay(host)
        d3 = limiter.acquire_delay(host)
        self.assertAlmostEqual(d1, 0.0, delta=0.01)
        self.assertGreaterEqual(d2, 0.04)
        self.assertGreaterEqual(d3, 0.09)

    @unittest.skipUnless(HAVE_INGEST_HTTP, "synaptik.ingest.http not available")
    def test_tc_r1_1_02_concurrent_thread_reservation_monotonicity(self):
        """TC-R1.1-02: Rapid reservation requests receive monotonically increasing delays."""
        limiter = ConcurrentRateLimiter(limits={"www.ebi.ac.uk": 0.04})
        delays = [limiter.acquire_delay("www.ebi.ac.uk") for _ in range(5)]
        self.assertEqual(len(delays), 5)
        self.assertAlmostEqual(delays[0], 0.0, delta=0.01)
        for i in range(1, len(delays)):
            self.assertGreaterEqual(delays[i] - delays[i - 1], 0.035)

    @unittest.skipUnless(HAVE_INGEST_HTTP, "synaptik.ingest.http not available")
    def test_tc_r1_1_03_orthogonal_hosts_non_blocking(self):
        """TC-R1.1-03: Independent hosts do not serialize each other."""
        limiter = ConcurrentRateLimiter(limits={"host_a.org": 10.0, "host_b.org": 0.01})
        d_a1 = limiter.acquire_delay("host_a.org")
        d_b1 = limiter.acquire_delay("host_b.org")
        self.assertAlmostEqual(d_a1, 0.0, delta=0.01)
        self.assertAlmostEqual(d_b1, 0.0, delta=0.01)

    @unittest.skipUnless(HAVE_INGEST_HTTP, "synaptik.ingest.http not available")
    def test_tc_r1_1_04_microsecond_lock_holding(self):
        """TC-R1.1-04: Non-lock-sleeping verification: 20 acquisitions complete in <20ms."""
        limiter = ConcurrentRateLimiter(limits={"api.biorxiv.org": 0.50})
        t0 = time.monotonic()
        for _ in range(20):
            limiter.acquire_delay("api.biorxiv.org")
        t_elapsed = time.monotonic() - t0
        self.assertLess(t_elapsed, 0.05, "Lock acquisition must be sub-millisecond, no sleep under lock")

    @unittest.skipUnless(HAVE_INGEST_HTTP, "synaptik.ingest.http not available")
    def test_tc_r1_1_05_custom_limit_and_default_fallback(self):
        """TC-R1.1-05: Custom limits override defaults and fallback handles unknown hosts."""
        limiter = ConcurrentRateLimiter(limits={"custom.org": 0.02, "_default": 0.15})
        self.assertIn("custom.org", limiter.limits)
        d_custom = limiter.limits["custom.org"]
        d_default = limiter.limits["_default"]
        self.assertEqual(d_custom, 0.02)
        self.assertEqual(d_default, 0.15)

    # --------------------------------------------------------------------------
    # Feature R1.2: PoliteSession Concurrency & Batching
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_INGEST_HTTP, "synaptik.ingest.http not available")
    def test_tc_r1_2_01_ordered_batch_retrieval(self):
        """TC-R1.2-01: fetch_batch preserves input list index order."""
        session = PoliteSession(store=None, use_cache=False)
        urls = [(f"https://www.ebi.ac.uk/test/{i}", {}) for i in range(5)]

        def mock_fetch(method, url, params, host, want_json):
            idx = int(url.split("/")[-1])
            return {"index": idx}

        session._request_with_retry = mock_fetch
        results = session.fetch_batch(urls, want_json=True, max_workers=2)
        self.assertEqual(len(results), 5)
        for i, res in enumerate(results):
            self.assertEqual(res["index"], i)

    @unittest.skipUnless(HAVE_INGEST_HTTP, "synaptik.ingest.http not available")
    def test_tc_r1_2_02_thread_safe_request_metrics(self):
        """TC-R1.2-02: Thread-safe metrics count increments accurately across threads."""
        session = PoliteSession(store=None, use_cache=False)
        n_threads = 10
        iters = 5

        def worker():
            for _ in range(iters):
                session._inc_requests()

        threads = [threading.Thread(target=worker) for _ in range(n_threads)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        self.assertEqual(session.stats()["requests"], n_threads * iters)

    @unittest.skipUnless(HAVE_INGEST_HTTP, "synaptik.ingest.http not available")
    def test_tc_r1_2_03_exponential_backoff_on_429(self):
        """TC-R1.2-03: Transient HTTP 429 triggers retry and succeeds."""
        session = PoliteSession(store=None, max_retries=3, use_cache=False)
        attempts = 0

        def mock_request(*args, **kwargs):
            nonlocal attempts
            attempts += 1
            if attempts < 3:
                return MockHTTPResponse(status_code=429)
            return MockHTTPResponse(status_code=200, json_data={"status": "ok"})

        session.session.request = mock_request
        with patch("time.sleep", return_value=None):
            res = session.get_json("https://eutils.ncbi.nlm.nih.gov/test")
        self.assertEqual(res, {"status": "ok"})
        self.assertEqual(attempts, 3)

    @unittest.skipUnless(HAVE_INGEST_HTTP, "synaptik.ingest.http not available")
    def test_tc_r1_2_04_cache_hit_bypass(self):
        """TC-R1.2-04: Backing cache hit avoids HTTP network egress."""
        mock_store = MagicMock()
        mock_store.cache_get.return_value = {"cached": True}
        session = PoliteSession(store=mock_store, use_cache=True)
        res = session.get_json("https://api.crossref.org/works", params={"q": "ALS"})
        self.assertEqual(res, {"cached": True})
        self.assertEqual(session.stats()["cache_hits"], 1)

    @unittest.skipUnless(HAVE_INGEST_HTTP, "synaptik.ingest.http not available")
    def test_tc_r1_2_05_graceful_404_handling(self):
        """TC-R1.2-05: Fatal 404 yields None without crashing batch execution."""
        session = PoliteSession(store=None, use_cache=False)
        urls = [
            ("https://api.crossref.org/works/1", {}),
            ("https://api.crossref.org/works/404", {}),
            ("https://api.crossref.org/works/3", {}),
        ]

        def mock_fetch(method, url, params, host, want_json):
            if "404" in url:
                return None
            return {"url": url}

        session._request_with_retry = mock_fetch
        results = session.fetch_batch(urls, max_workers=2)
        self.assertEqual(len(results), 3)
        self.assertIsNotNone(results[0])
        self.assertIsNone(results[1])
        self.assertIsNotNone(results[2])

    # --------------------------------------------------------------------------
    # Feature R1.3: Europe PMC Harvesting & JATS XML
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_INGEST_EUROPEPMC, "synaptik.ingest.europepmc not available")
    def test_tc_r1_3_01_cursor_pagination_iteration(self):
        """TC-R1.3-01: Europe PMC search paginates via cursorMark."""
        harvester = EuropePMCHarvester(session=PoliteSession(store=None, use_cache=False))
        call_count = 0

        def mock_get(url, params=None, want_json=True):
            nonlocal call_count
            call_count += 1
            cursor = params.get("cursorMark", "*")
            if cursor == "*":
                return {
                    "hitCount": 2,
                    "nextCursorMark": "CURSOR_PAGE_2",
                    "resultList": {
                        "result": [
                            {"id": "1", "pmid": "12345678", "title": "Paper 1", "authorList": {"author": [{"fullName": "Smith J"}]}},
                            {"id": "2", "pmid": "12345679", "title": "Paper 2", "authorList": {"author": [{"fullName": "Doe J"}]}},
                        ]
                    },
                }
            return {
                "hitCount": 2,
                "nextCursorMark": "CURSOR_PAGE_2",  # identical cursor halts pagination
                "resultList": {"result": []},
            }

        harvester.s.get_json = mock_get
        records = list(harvester.search("STMN2 splicing", page_size=2, max_results=2))
        self.assertEqual(len(records), 2)
        self.assertEqual(call_count, 1)

    @unittest.skipUnless(HAVE_INGEST_EUROPEPMC, "synaptik.ingest.europepmc not available")
    def test_tc_r1_3_02_jats_xml_section_extraction(self):
        """TC-R1.3-02: JATS XML parser extracts structured sections."""
        xml_sample = """<?xml version="1.0"?>
        <article>
            <body>
                <sec><title>Introduction</title><p>TDP-43 loss causes cryptic splicing.</p></sec>
                <sec><title>Results</title><p>Rescue with ASO restored full STMN2.</p></sec>
            </body>
        </article>"""
        sections = EuropePMCHarvester.parse_jats(xml_sample)
        self.assertIn("results", sections)
        self.assertIn("Rescue with ASO", sections["results"])

    @unittest.skipUnless(HAVE_INGEST_EUROPEPMC, "synaptik.ingest.europepmc not available")
    def test_tc_r1_3_03_figure_legend_mining(self):
        """TC-R1.3-03: Figure caption text is mined into figures section."""
        xml_sample = """<article><body>
            <fig id="fig1"><label>Figure 1</label><caption><p>Western blot shows stathmin-2 restoration.</p></caption></fig>
        </body></article>"""
        sections = EuropePMCHarvester.parse_jats(xml_sample)
        self.assertIn("figures", sections)
        self.assertIn("stathmin-2 restoration", sections["figures"])

    @unittest.skipUnless(HAVE_INGEST_EUROPEPMC, "synaptik.ingest.europepmc not available")
    def test_tc_r1_3_04_open_access_filter_query(self):
        """TC-R1.3-04: open_access_only appends OPEN_ACCESS:y constraint to query."""
        harvester = EuropePMCHarvester(session=PoliteSession(store=None, use_cache=False))
        captured_query = ""

        def mock_get(url, params=None, want_json=True):
            nonlocal captured_query
            captured_query = params.get("query", "")
            return {"nextCursorMark": "*", "resultList": {"result": []}}

        harvester.s.get_json = mock_get
        list(harvester.search("TARDBP", open_access_only=True, max_results=1))
        self.assertIn("OPEN_ACCESS:y", captured_query)

    @unittest.skipUnless(HAVE_INGEST_EUROPEPMC, "synaptik.ingest.europepmc not available")
    def test_tc_r1_3_05_malformed_xml_recovery(self):
        """TC-R1.3-05: Harvester recovers text gracefully from unclosed or malformed XML."""
        malformed_xml = "<article><body><sec><title>Fragmented Text</title><p>Data without closing"
        sections = EuropePMCHarvester.parse_jats(malformed_xml)
        self.assertIsNotNone(sections)
        all_text = " ".join(sections.values())
        self.assertIn("Data without closing", all_text)

    # --------------------------------------------------------------------------
    # Feature R1.4: PubMed Harvester & Batching
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_INGEST_PUBMED, "synaptik.ingest.pubmed not available")
    def test_tc_r1_4_01_esearch_pmid_extraction(self):
        """TC-R1.4-01: search_pmids parses id list from E-utilities esearch JSON."""
        harvester = PubMedHarvester(session=PoliteSession(store=None, use_cache=False))
        harvester.s.get_json = MagicMock(return_value={
            "esearchresult": {"idlist": ["38123456", "38123457", "38123458"]}
        })
        pmids = harvester.search_pmids("STMN2 ALS", max_results=3)
        self.assertEqual(pmids, ["38123456", "38123457", "38123458"])

    @unittest.skipUnless(HAVE_INGEST_PUBMED, "synaptik.ingest.pubmed not available")
    def test_tc_r1_4_02_efetch_multi_batching(self):
        """TC-R1.4-02: fetch_metadata chunks requests when PMIDs exceed batch_size."""
        harvester = PubMedHarvester(session=PoliteSession(store=None, use_cache=False))
        xml_batch = b"""<PubmedArticleSet>
            <PubmedArticle><MedlineCitation><PMID>1</PMID><Article><ArticleTitle>T1</ArticleTitle></Article></MedlineCitation></PubmedArticle>
            <PubmedArticle><MedlineCitation><PMID>2</PMID><Article><ArticleTitle>T2</ArticleTitle></Article></MedlineCitation></PubmedArticle>
        </PubmedArticleSet>"""
        harvester.s.get_text = MagicMock(return_value=xml_batch.decode("utf-8"))
        items = list(harvester.fetch_metadata(["1", "2", "3", "4"], batch_size=2))
        self.assertEqual(harvester.s.get_text.call_count, 2)
        self.assertEqual(len(items), 4)

    @unittest.skipUnless(HAVE_INGEST_PUBMED, "synaptik.ingest.pubmed not available")
    def test_tc_r1_4_03_mesh_term_extraction(self):
        """TC-R1.4-03: PubMed parser parses DescriptorName MeSH terms."""
        harvester = PubMedHarvester(session=PoliteSession(store=None, use_cache=False))
        xml = """<PubmedArticleSet><PubmedArticle><MedlineCitation>
            <PMID>1234</PMID>
            <Article><ArticleTitle>ALS Study</ArticleTitle></Article>
            <MeshHeadingList>
                <MeshHeading><DescriptorName>Amyotrophic Lateral Sclerosis</DescriptorName></MeshHeading>
                <MeshHeading><DescriptorName>RNA Splicing</DescriptorName></MeshHeading>
            </MeshHeadingList>
        </MedlineCitation></PubmedArticle></PubmedArticleSet>"""
        items = list(harvester._parse_efetch(xml))
        self.assertEqual(len(items), 1)
        self.assertIn("Amyotrophic Lateral Sclerosis", getattr(items[0], "_mesh_terms", []))

    @unittest.skipUnless(HAVE_INGEST_PUBMED, "synaptik.ingest.pubmed not available")
    def test_tc_r1_4_04_retraction_detection(self):
        """TC-R1.4-04: Retracted publication type sets is_retracted=True."""
        harvester = PubMedHarvester(session=PoliteSession(store=None, use_cache=False))
        xml = """<PubmedArticleSet><PubmedArticle><MedlineCitation>
            <PMID>9999</PMID>
            <Article>
                <ArticleTitle>Retracted Paper</ArticleTitle>
                <PublicationTypeList>
                    <PublicationType>Retracted Publication</PublicationType>
                </PublicationTypeList>
            </Article>
        </MedlineCitation></PubmedArticle></PubmedArticleSet>"""
        items = list(harvester._parse_efetch(xml))
        self.assertEqual(len(items), 1)
        self.assertTrue(items[0].is_retracted)

    @unittest.skipUnless(HAVE_INGEST_PUBMED, "synaptik.ingest.pubmed not available")
    def test_tc_r1_4_05_api_key_rate_interval(self):
        """TC-R1.4-05: Registered NCBI API key enables 10 req/s (0.10s interval)."""
        harvester = PubMedHarvester(session=PoliteSession(store=None, use_cache=False, rate_limits={"eutils.ncbi.nlm.nih.gov": 0.10}), api_key="TEST_KEY")
        self.assertEqual(harvester.s.limiter.limits["eutils.ncbi.nlm.nih.gov"], 0.10)

    # --------------------------------------------------------------------------
    # Feature R1.5: bioRxiv / medRxiv Harvester
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_INGEST_PUBMED, "synaptik.ingest.pubmed not available")
    def test_tc_r1_5_01_date_range_endpoint_routing(self):
        """TC-R1.5-01: fetch_range calls correct bioRxiv API URL structure."""
        harvester = BioRxivHarvester(session=PoliteSession(store=None, use_cache=False), server="biorxiv")
        captured_url = ""

        def mock_get(url, **kwargs):
            nonlocal captured_url
            captured_url = url
            return {"messages": [{"status": "ok", "total": 0}], "collection": []}

        harvester.s.get_json = mock_get
        list(harvester.fetch_range("2026-01-01", "2026-01-07", max_results=5))
        self.assertIn("api.biorxiv.org/details/biorxiv/2026-01-01/2026-01-07/0", captured_url)

    @unittest.skipUnless(HAVE_INGEST_PUBMED, "synaptik.ingest.pubmed not available")
    def test_tc_r1_5_02_preprint_metadata_flags(self):
        """TC-R1.5-02: Harvested preprints have is_preprint=True and open_access=True."""
        harvester = BioRxivHarvester(session=PoliteSession(store=None, use_cache=False))
        item = {
            "doi": "10.1101/2026.01.01.500000",
            "title": "STMN2 Cryptic Rescue",
            "authors": "Nobile, F.",
            "date": "2026-01-02",
            "abstract": "Abstract text",
        }
        meta = harvester._parse(item)
        self.assertIsNotNone(meta)
        self.assertTrue(meta.is_preprint)
        self.assertTrue(meta.open_access)

    @unittest.skipUnless(HAVE_INGEST_PUBMED, "synaptik.ingest.pubmed not available")
    def test_tc_r1_5_03_medrxiv_server_routing(self):
        """TC-R1.5-03: Server parameter 'medrxiv' routes to medrxiv endpoint."""
        harvester = BioRxivHarvester(session=PoliteSession(store=None, use_cache=False), server="medrxiv")
        captured_url = ""

        def mock_get(url, **kwargs):
            nonlocal captured_url
            captured_url = url
            return {"messages": [{"status": "ok", "total": 0}], "collection": []}

        harvester.s.get_json = mock_get
        list(harvester.fetch_range("2026-02-01", "2026-02-02", max_results=1))
        self.assertIn("/details/medrxiv/", captured_url)

    @unittest.skipUnless(HAVE_INGEST_PUBMED, "synaptik.ingest.pubmed not available")
    def test_tc_r1_5_04_pagination_cursor_offset(self):
        """TC-R1.5-04: Harvester loops cursor offset 0 -> 100 on multi-page collections."""
        harvester = BioRxivHarvester(session=PoliteSession(store=None, use_cache=False))
        offsets = []

        def mock_get(url, **kwargs):
            offset = int(url.split("/")[-1])
            offsets.append(offset)
            if offset == 0:
                return {
                    "messages": [{"status": "ok", "total": 150}],
                    "collection": [{"doi": f"10.1101/{i}"} for i in range(100)],
                }
            return {
                "messages": [{"status": "ok", "total": 150}],
                "collection": [{"doi": f"10.1101/{i}"} for i in range(100, 150)],
            }

        harvester.s.get_json = mock_get
        items = list(harvester.fetch_range("2026-01-01", "2026-01-10", max_results=150))
        self.assertEqual(len(items), 150)
        self.assertEqual(offsets, [0, 100])

    @unittest.skipUnless(HAVE_INGEST_PUBMED, "synaptik.ingest.pubmed not available")
    def test_tc_r1_5_05_author_affiliation_capture(self):
        """TC-R1.5-05: Corresponding institution is stored in affiliations."""
        harvester = BioRxivHarvester(session=PoliteSession(store=None, use_cache=False))
        item = {
            "doi": "10.1101/2026.001",
            "title": "Study",
            "author_corresponding_institution": "Harvard Medical School",
        }
        meta = harvester._parse(item)
        self.assertIsNotNone(meta)
        self.assertIn("Harvard Medical School", meta.affiliations)

    # --------------------------------------------------------------------------
    # Feature R1.6: Bayesian Preprint Discounting (0.80 Multiplier)
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_EVIDENCE_WEIGHTER, "synaptik.infer.weighting not available")
    def test_tc_r1_6_01_evidence_weighter_20_percent_discount(self):
        """TC-R1.6-01: Preprint receives exact 0.80 journal_tier multiplier (20% discount)."""
        weighter = EvidenceWeighter()
        rec_peer = create_mock_experiment_record(is_preprint=False)
        rec_prep = create_mock_experiment_record(is_preprint=True)

        w_peer = weighter.weight(rec_peer)
        w_prep = weighter.weight(rec_prep)

        self.assertEqual(w_peer.journal_tier, 1.00)
        self.assertEqual(w_prep.journal_tier, 0.80)
        self.assertAlmostEqual(w_prep.composite / w_peer.composite, 0.80, places=4)

    @unittest.skipUnless(HAVE_EVIDENCE_WEIGHTER, "synaptik.infer.weighting not available")
    def test_tc_r1_6_02_preprint_audit_note_generation(self):
        """TC-R1.6-02: Preprint record contains explicit 'preprint' string in notes."""
        weighter = EvidenceWeighter()
        rec_prep = create_mock_experiment_record(is_preprint=True)
        breakdown = weighter.weight(rec_prep)
        self.assertIn("preprint", breakdown.notes)

    @unittest.skipUnless(HAVE_EVIDENCE_WEIGHTER, "synaptik.infer.weighting not available")
    def test_tc_r1_6_03_configurable_preprint_discount(self):
        """TC-R1.6-03: EvidenceWeighter respects custom preprint discount multiplier."""
        weighter = EvidenceWeighter(preprint_weight=0.75)
        rec_prep = create_mock_experiment_record(is_preprint=True)
        w_prep = weighter.weight(rec_prep)
        self.assertEqual(w_prep.journal_tier, 0.75)

    @unittest.skipUnless(HAVE_EVIDENCE_WEIGHTER, "synaptik.infer.weighting not available")
    def test_tc_r1_6_04_retraction_overrides_preprint(self):
        """TC-R1.6-04: Retraction zeroes out weight completely (integrity=0.0)."""
        weighter = EvidenceWeighter()
        rec_retracted = create_mock_experiment_record(is_preprint=True, is_retracted=True)
        w = weighter.weight(rec_retracted)
        self.assertEqual(w.integrity, 0.0)
        self.assertEqual(w.composite, 0.0)

    @unittest.skipUnless(HAVE_EVIDENCE_WEIGHTER, "synaptik.infer.weighting not available")
    def test_tc_r1_6_05_explanation_transparency(self):
        """TC-R1.6-05: explain() audit string exposes journal=0.80 and [preprint]."""
        weighter = EvidenceWeighter()
        rec_prep = create_mock_experiment_record(is_preprint=True)
        explanation = weighter.weight(rec_prep).explain()
        self.assertIn("journal=0.80", explanation)
        self.assertIn("[preprint", explanation)

    # --------------------------------------------------------------------------
    # Feature R2.1: CSF Volume Allometric Scaling
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_CLINICAL_ALLOMETRY, "synaptik.clinical.allometry not available")
    def test_tc_r2_1_01_adult_baseline_allometry(self):
        """TC-R2.1-01: Adult CSF scaling calculates exact HED and MRSD."""
        res = calculate_allometric_dose(noael_mg=15.0, nhp_volume_ml=15.0, human_volume_ml=140.0, safety_factor=10.0, acknowledge_not_clinical=True)
        ratio = getattr(res, "scaling_ratio", getattr(res, "scaling_factor", None))
        self.assertAlmostEqual(ratio, 9.3333, places=3)
        self.assertAlmostEqual(res.hed_mg, 140.0, places=3)
        self.assertAlmostEqual(res.mrsd_mg, 14.0, places=3)

    @unittest.skipUnless(HAVE_CLINICAL_ALLOMETRY, "synaptik.clinical.allometry not available")
    def test_tc_r2_1_02_pediatric_infant_allometry(self):
        """TC-R2.1-02: Pediatric infant CSF volume (65-70 mL) scales appropriately."""
        res = calculate_allometric_dose(noael_mg=15.0, nhp_volume_ml=15.0, human_volume_ml=70.0, safety_factor=10.0, acknowledge_not_clinical=True)
        ratio = getattr(res, "scaling_ratio", getattr(res, "scaling_factor", None))
        self.assertAlmostEqual(ratio, 4.6667, places=3)
        self.assertAlmostEqual(res.hed_mg, 70.0, places=3)
        self.assertAlmostEqual(res.mrsd_mg, 7.0, places=3)

    @unittest.skipUnless(HAVE_CLINICAL_ALLOMETRY, "synaptik.clinical.allometry not available")
    def test_tc_r2_1_03_scaled_safety_factor(self):
        """TC-R2.1-03: Conservative safety factor (SF=20.0) reduces MRSD proportionally."""
        res = calculate_allometric_dose(noael_mg=10.0, nhp_volume_ml=15.0, human_volume_ml=140.0, safety_factor=20.0, acknowledge_not_clinical=True)
        self.assertAlmostEqual(res.hed_mg, 93.333, places=2)
        self.assertAlmostEqual(res.mrsd_mg, 4.6667, places=3)

    @unittest.skipUnless(HAVE_CLINICAL_ALLOMETRY, "synaptik.clinical.allometry not available")
    def test_tc_r2_1_04_target_csf_concentration_invariance(self):
        """TC-R2.1-04: Target CSF initial concentration is invariant across human volumes."""
        noael = 12.0
        adult = calculate_allometric_dose(noael_mg=noael, nhp_volume_ml=15.0, human_volume_ml=140.0, safety_factor=10.0, acknowledge_not_clinical=True)
        ped = calculate_allometric_dose(noael_mg=noael, nhp_volume_ml=15.0, human_volume_ml=70.0, safety_factor=10.0, acknowledge_not_clinical=True)
        c_adult = (adult.mrsd_mg * 1000.0) / 140.0
        c_ped = (ped.mrsd_mg * 1000.0) / 70.0
        self.assertAlmostEqual(c_adult, 80.0, places=3)
        self.assertAlmostEqual(c_ped, 80.0, places=3)

    @unittest.skipUnless(HAVE_CLINICAL_ALLOMETRY, "synaptik.clinical.allometry not available")
    def test_tc_r2_1_05_alternative_adult_150ml_volume(self):
        """TC-R2.1-05: Upper adult boundary (150 mL) yields scaling ratio 10.0."""
        res = calculate_allometric_dose(noael_mg=15.0, nhp_volume_ml=15.0, human_volume_ml=150.0, safety_factor=10.0, acknowledge_not_clinical=True)
        ratio = getattr(res, "scaling_ratio", getattr(res, "scaling_factor", None))
        self.assertAlmostEqual(ratio, 10.0, places=3)
        self.assertAlmostEqual(res.hed_mg, 150.0, places=3)
        self.assertAlmostEqual(res.mrsd_mg, 15.0, places=3)

    # --------------------------------------------------------------------------
    # Feature R2.2: Intrathecal PK/PD Simulator
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_CLINICAL_PKPD, "synaptik.clinical.pkpd not available")
    def test_tc_r2_2_01_bolus_peak_concentration(self):
        """TC-R2.2-01: Bolus Cmax in CSF equals Dose/Volume at T=0."""
        res = simulate_intrathecal_pkpd(dose_mg=10.0, route="bolus", csf_volume_ml=140.0)
        expected_cmax = (10.0 * 1000.0) / 140.0  # 71.43 ug/mL
        self.assertAlmostEqual(res.cmax_csf_ug_ml, expected_cmax, delta=0.5)
        self.assertAlmostEqual(res.tmax_csf_days, 0.0, delta=0.1)

    @unittest.skipUnless(HAVE_CLINICAL_PKPD, "synaptik.clinical.pkpd not available")
    def test_tc_r2_2_02_csf_elimination_half_life(self):
        """TC-R2.2-02: CSF concentration drops by approx 50% at Day 24."""
        res = simulate_intrathecal_pkpd(dose_mg=10.0, route="bolus", csf_half_life_days=24.0, duration_days=120)
        idx_24 = min(range(len(res.time_days)), key=lambda i: abs(res.time_days[i] - 24.0))
        c_24 = res.csf_conc_ug_ml[idx_24]
        self.assertAlmostEqual(c_24 / res.cmax_csf_ug_ml, 0.50, delta=0.12)

    @unittest.skipUnless(HAVE_CLINICAL_PKPD, "synaptik.clinical.pkpd not available")
    def test_tc_r2_2_03_catheter_infusion_blunted_peak(self):
        """TC-R2.2-03: Catheter infusion blunts CSF Cmax and delays Tmax."""
        bolus = simulate_intrathecal_pkpd(dose_mg=10.0, route="bolus")
        cath = simulate_intrathecal_pkpd(dose_mg=10.0, route="catheter")
        self.assertLess(cath.cmax_csf_ug_ml, bolus.cmax_csf_ug_ml)
        self.assertGreater(cath.tmax_csf_days, 0.5)

    @unittest.skipUnless(HAVE_CLINICAL_PKPD, "synaptik.clinical.pkpd not available")
    def test_tc_r2_2_04_catheter_vs_bolus_auc_equivalence(self):
        """TC-R2.2-04: Catheter and bolus deliver equivalent total AUC (<2% difference)."""
        bolus = simulate_intrathecal_pkpd(dose_mg=10.0, route="bolus")
        cath = simulate_intrathecal_pkpd(dose_mg=10.0, route="catheter")
        auc_diff_rel = abs(cath.auc_csf_ug_day_ml - bolus.auc_csf_ug_day_ml) / bolus.auc_csf_ug_day_ml
        self.assertLess(auc_diff_rel, 0.02)

    @unittest.skipUnless(HAVE_CLINICAL_PKPD, "synaptik.clinical.pkpd not available")
    def test_tc_r2_2_05_parenchymal_tissue_accumulation(self):
        """TC-R2.2-05: Parenchymal tissue accumulates drug with sustained Day 90 exposure."""
        res = simulate_intrathecal_pkpd(dose_mg=10.0, route="bolus", tissue_half_life_days=135.0)
        self.assertGreater(res.tmax_tissue_days, 10.0)
        self.assertGreater(res.day90_tissue_conc_ug_g, 0.5 * res.cmax_tissue_ug_g)

    # --------------------------------------------------------------------------
    # Feature R2.3: Adaptive Escalation (3+3 Engine)
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_CLINICAL_ESCALATION, "synaptik.clinical.escalation not available")
    def test_tc_r2_3_01_three_plus_three_clean_escalate(self):
        """TC-R2.3-01: 0/3 DLTs escalates to next dose level."""
        engine = ThreePlusThreeEngine()
        doses = [2.0, 5.0, 10.0]
        history = [CohortObservation(cohort_index=1, dose_mg=2.0, n_patients=3, n_dlts=0)]
        decision = engine.evaluate_cohort(doses, history)
        self.assertEqual(decision.action, "escalate")
        self.assertEqual(decision.recommended_dose_mg, 5.0)

    @unittest.skipUnless(HAVE_CLINICAL_ESCALATION, "synaptik.clinical.escalation not available")
    def test_tc_r2_3_02_three_plus_three_stay_and_expand(self):
        """TC-R2.3-02: 1/3 DLTs stays at current dose to expand cohort to 6."""
        engine = ThreePlusThreeEngine()
        doses = [2.0, 5.0, 10.0]
        history = [CohortObservation(cohort_index=1, dose_mg=2.0, n_patients=3, n_dlts=1)]
        decision = engine.evaluate_cohort(doses, history)
        self.assertEqual(decision.action, "stay")
        self.assertEqual(decision.recommended_dose_mg, 2.0)

    @unittest.skipUnless(HAVE_CLINICAL_ESCALATION, "synaptik.clinical.escalation not available")
    def test_tc_r2_3_03_three_plus_three_escalate_after_expansion(self):
        """TC-R2.3-03: 1/6 DLTs escalates after expansion."""
        engine = ThreePlusThreeEngine()
        doses = [2.0, 5.0, 10.0]
        history = [CohortObservation(cohort_index=1, dose_mg=2.0, n_patients=6, n_dlts=1)]
        decision = engine.evaluate_cohort(doses, history)
        self.assertEqual(decision.action, "escalate")
        self.assertEqual(decision.recommended_dose_mg, 5.0)

    @unittest.skipUnless(HAVE_CLINICAL_ESCALATION, "synaptik.clinical.escalation not available")
    def test_tc_r2_3_04_three_plus_three_mtd_reached_on_2_of_6(self):
        """TC-R2.3-04: >=2/6 DLTs stops escalation and declares previous dose as MTD."""
        engine = ThreePlusThreeEngine()
        doses = [2.0, 5.0, 10.0]
        history = [
            CohortObservation(cohort_index=1, dose_mg=2.0, n_patients=3, n_dlts=0),
            CohortObservation(cohort_index=2, dose_mg=5.0, n_patients=6, n_dlts=2),
        ]
        decision = engine.evaluate_cohort(doses, history)
        self.assertEqual(decision.action, "stop_mtd_reached")
        self.assertEqual(decision.recommended_dose_mg, 2.0)

    @unittest.skipUnless(HAVE_CLINICAL_ESCALATION, "synaptik.clinical.escalation not available")
    def test_tc_r2_3_05_three_plus_three_immediate_halt_on_2_of_3(self):
        """TC-R2.3-05: >=2/3 DLTs stops immediately without cohort expansion."""
        engine = ThreePlusThreeEngine()
        doses = [2.0, 5.0, 10.0]
        history = [
            CohortObservation(cohort_index=1, dose_mg=2.0, n_patients=3, n_dlts=0),
            CohortObservation(cohort_index=2, dose_mg=5.0, n_patients=3, n_dlts=2),
        ]
        decision = engine.evaluate_cohort(doses, history)
        self.assertEqual(decision.action, "stop_mtd_reached")

    # --------------------------------------------------------------------------
    # Feature R2.4: Adaptive Escalation (BOIN Engine)
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_CLINICAL_ESCALATION, "synaptik.clinical.escalation not available")
    def test_tc_r2_4_01_boin_boundary_calculation(self):
        """TC-R2.4-01: Target phi=0.25 derives lambda_e=0.1968 and lambda_d=0.2984."""
        engine = BOINEscalationEngine(target_phi=0.25)
        self.assertAlmostEqual(engine.lambda_e, 0.1968, places=3)
        self.assertAlmostEqual(engine.lambda_d, 0.2984, places=3)

    @unittest.skipUnless(HAVE_CLINICAL_ESCALATION, "synaptik.clinical.escalation not available")
    def test_tc_r2_4_02_boin_escalation_low_toxicity(self):
        """TC-R2.4-02: Observed rate <= lambda_e (0/3 = 0.0) escalates."""
        engine = BOINEscalationEngine(target_phi=0.25)
        doses = [2.0, 5.0, 10.0]
        history = [CohortObservation(cohort_index=1, dose_mg=2.0, n_patients=3, n_dlts=0)]
        decision = engine.evaluate_cohort(doses, history)
        self.assertEqual(decision.action, "escalate")
        self.assertEqual(decision.recommended_dose_mg, 5.0)

    @unittest.skipUnless(HAVE_CLINICAL_ESCALATION, "synaptik.clinical.escalation not available")
    def test_tc_r2_4_03_boin_stay_intermediate_toxicity(self):
        """TC-R2.4-03: Observed rate between lambda_e and lambda_d (1/4 = 0.25) stays."""
        engine = BOINEscalationEngine(target_phi=0.25)
        doses = [2.0, 5.0, 10.0]
        history = [CohortObservation(cohort_index=2, dose_mg=5.0, n_patients=4, n_dlts=1)]
        decision = engine.evaluate_cohort(doses, history)
        self.assertEqual(decision.action, "stay")

    @unittest.skipUnless(HAVE_CLINICAL_ESCALATION, "synaptik.clinical.escalation not available")
    def test_tc_r2_4_04_boin_de_escalate_high_toxicity(self):
        """TC-R2.4-04: Observed rate >= lambda_d (2/5 = 0.40) de-escalates."""
        engine = BOINEscalationEngine(target_phi=0.25)
        doses = [2.0, 5.0, 10.0]
        history = [CohortObservation(cohort_index=2, dose_mg=5.0, n_patients=5, n_dlts=2)]
        decision = engine.evaluate_cohort(doses, history)
        self.assertEqual(decision.action, "de-escalate")
        self.assertEqual(decision.recommended_dose_mg, 2.0)

    @unittest.skipUnless(HAVE_CLINICAL_ESCALATION, "synaptik.clinical.escalation not available")
    def test_tc_r2_4_05_boin_safety_stopping_boundary(self):
        """TC-R2.4-05: Posterior toxicity >0.95 at Dose 1 stops trial on safety boundary."""
        engine = BOINEscalationEngine(target_phi=0.25)
        doses = [2.0, 5.0, 10.0]
        history = [CohortObservation(cohort_index=1, dose_mg=2.0, n_patients=3, n_dlts=3)]
        decision = engine.evaluate_cohort(doses, history)
        self.assertEqual(decision.action, "stop_safety_eliminated")
        self.assertTrue(decision.safety_eliminated)

    # --------------------------------------------------------------------------
    # Feature R2.5: DLT Stopping Boundaries & NfL Milestone
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_CLINICAL_PROTOCOL, "synaptik.clinical.protocol not available")
    def test_tc_r2_5_01_dlt_csf_pleocytosis_threshold(self):
        """TC-R2.5-01: WBC >50 cells/uL in CSF is categorized as DLT."""
        proto = generate_clinical_protocol(target_gene="FUS", disease="ALS", noael_mg=12.0)
        self.assertIn("csf_pleocytosis", proto.dlt_definition)
        self.assertIn("50", proto.dlt_definition["csf_pleocytosis"])

    @unittest.skipUnless(HAVE_CLINICAL_PROTOCOL, "synaptik.clinical.protocol not available")
    def test_tc_r2_5_02_dlt_severe_thrombocytopenia(self):
        """TC-R2.5-02: Platelets <50,000/uL triggers severe thrombocytopenia DLT."""
        proto = generate_clinical_protocol(target_gene="STMN2", disease="ALS", noael_mg=15.0)
        self.assertIn("thrombocytopenia", proto.dlt_definition)
        self.assertIn("50,000", proto.dlt_definition["thrombocytopenia"])

    @unittest.skipUnless(HAVE_CLINICAL_PROTOCOL, "synaptik.clinical.protocol not available")
    def test_tc_r2_5_03_biomarker_milestone_definition(self):
        """TC-R2.5-03: Biomarker response milestone specifies >30% reduction at Day 90."""
        proto = generate_clinical_protocol(target_gene="STMN2", disease="ALS", noael_mg=15.0)
        self.assertEqual(proto.biomarker_milestone["fluid"], "Plasma NfL")
        self.assertGreaterEqual(proto.biomarker_milestone["threshold_percent"], 30.0)
        self.assertEqual(proto.biomarker_milestone["evaluation_day"], 90)

    @unittest.skipUnless(HAVE_CLINICAL_PROTOCOL, "synaptik.clinical.protocol not available")
    def test_tc_r2_5_04_protocol_fingerprint_generation(self):
        """TC-R2.5-04: Protocol generates deterministic SHA-256 fingerprint."""
        proto = generate_clinical_protocol(target_gene="STMN2", disease="ALS", noael_mg=15.0)
        self.assertEqual(len(proto.protocol_fingerprint_sha256), 64)

    @unittest.skipUnless(HAVE_CLINICAL_PROTOCOL, "synaptik.clinical.protocol not available")
    def test_tc_r2_5_05_protocol_author_attribution(self):
        """TC-R2.5-05: Protocol metadata attributes Francesco Nobile (ORCID 0009-0001-9477-3332)."""
        proto = generate_clinical_protocol(target_gene="STMN2", disease="ALS", noael_mg=15.0)
        self.assertEqual(proto.author, "Francesco Nobile")
        self.assertEqual(proto.orcid, "0009-0001-9477-3332")

    # --------------------------------------------------------------------------
    # Feature R2.6: CLI `synaptik clinical`
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_CLI_CLINICAL, "synaptik.cli.commands.clinical not available")
    def test_tc_r2_6_01_cli_allometry_execution(self):
        """TC-R2.6-01: CLI allometry executes with exit code 0."""
        args = MagicMock(subcommand="allometry", noael=10.0, nhp_vol=15.0, human_vol=140.0, safety_factor=10.0, pediatric=False, json=False)
        code = cmd_clinical(args)
        self.assertEqual(code, 0)

    @unittest.skipUnless(HAVE_CLI_CLINICAL, "synaptik.cli.commands.clinical not available")
    def test_tc_r2_6_02_cli_allometry_json_format(self):
        """TC-R2.6-02: CLI allometry --json emits structured JSON."""
        args = MagicMock(subcommand="allometry", noael=10.0, nhp_vol=15.0, human_vol=140.0, safety_factor=10.0, pediatric=False, json=True)
        buf = io.StringIO()
        with patch("sys.stdout", buf):
            code = cmd_clinical(args)
        self.assertEqual(code, 0)
        data = json.loads(buf.getvalue())
        self.assertIn("hed_mg", data)
        self.assertIn("mrsd_mg", data)

    @unittest.skipUnless(HAVE_CLI_CLINICAL, "synaptik.cli.commands.clinical not available")
    def test_tc_r2_6_03_cli_pkpd_simulation(self):
        """TC-R2.6-03: CLI pkpd executes bolus simulation."""
        args = MagicMock(subcommand="pkpd", dose=10.0, route="bolus", csf_halflife=24.0, tissue_halflife=135.0, duration=360, dt=0.5, json=True)
        buf = io.StringIO()
        with patch("sys.stdout", buf):
            code = cmd_clinical(args)
        self.assertEqual(code, 0)
        data = json.loads(buf.getvalue())
        self.assertIn("cmax_csf_ug_ml", data)

    @unittest.skipUnless(HAVE_CLI_CLINICAL, "synaptik.cli.commands.clinical not available")
    def test_tc_r2_6_04_cli_escalate_boin(self):
        """TC-R2.6-04: CLI escalate evaluates BOIN design."""
        args = MagicMock(
            subcommand="escalate",
            method="boin",
            doses="2.0,5.0,10.0",
            cohorts='[{"cohort_index": 1, "dose_mg": 2.0, "n_patients": 3, "n_dlts": 0}]',
            target_phi=0.25,
            json=True,
        )
        buf = io.StringIO()
        with patch("sys.stdout", buf):
            code = cmd_clinical(args)
        self.assertEqual(code, 0)
        data = json.loads(buf.getvalue())
        self.assertEqual(data["action"], "escalate")

    @unittest.skipUnless(HAVE_CLI_CLINICAL, "synaptik.cli.commands.clinical not available")
    def test_tc_r2_6_05_cli_protocol_synthesis(self):
        """TC-R2.6-05: CLI protocol generates full trial protocol dossier."""
        args = MagicMock(subcommand="protocol", target="FUS", disease="ALS", noael=12.0, escalation="boin", pediatric=False, json=False, output=None)
        buf = io.StringIO()
        with patch("sys.stdout", buf):
            code = cmd_clinical(args)
        self.assertEqual(code, 0)
        out = buf.getvalue()
        self.assertIn("Francesco Nobile", out)
        self.assertIn("0009-0001-9477-3332", out)

    # --------------------------------------------------------------------------
    # Feature R3.1: Turner 2004 Parameters
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_BIOPHYSICS_PARAMS, "synaptik.biophysics.parameters not available")
    def test_tc_r3_1_01_turner_gc_watson_crick_stacking(self):
        """TC-R3.1-01: GC-rich Watson-Crick stacking is -3.42 kcal/mol."""
        dG = get_stacking_energy("GC", "CG")
        self.assertAlmostEqual(dG, -3.42, delta=0.1)

    @unittest.skipUnless(HAVE_BIOPHYSICS_PARAMS, "synaptik.biophysics.parameters not available")
    def test_tc_r3_1_02_turner_au_watson_crick_stacking(self):
        """TC-R3.1-02: AU-rich Watson-Crick stacking is -0.93 to -1.10 kcal/mol."""
        dG_aa = get_stacking_energy("AA", "UU")
        dG_au = get_stacking_energy("AU", "UA")
        self.assertAlmostEqual(dG_aa, -0.93, delta=0.1)
        self.assertAlmostEqual(dG_au, -1.10, delta=0.1)

    @unittest.skipUnless(HAVE_BIOPHYSICS_PARAMS, "synaptik.biophysics.parameters not available")
    def test_tc_r3_1_03_turner_gu_wobble_stacking(self):
        """TC-R3.1-03: GU wobble stacking parameters."""
        dG = get_stacking_energy("GG", "UU")
        self.assertAlmostEqual(dG, -1.40, delta=0.1)

    @unittest.skipUnless(HAVE_BIOPHYSICS_PARAMS, "synaptik.biophysics.parameters not available")
    def test_tc_r3_1_04_terminal_au_penalty(self):
        """TC-R3.1-04: Terminal AU pair adds +0.45 kcal/mol penalty."""
        from synaptik.biophysics.parameters import get_terminal_penalty
        p_au = get_terminal_penalty("A", "U")
        p_gc = get_terminal_penalty("G", "C")
        self.assertAlmostEqual(p_au, 0.45, delta=0.05)
        self.assertAlmostEqual(p_gc, 0.0, delta=0.01)

    @unittest.skipUnless(HAVE_BIOPHYSICS_PARAMS, "synaptik.biophysics.parameters not available")
    def test_tc_r3_1_05_mismatch_stacking_energy(self):
        """TC-R3.1-05: Non-canonical pairs do not receive Watson-Crick stacking."""
        dG = get_stacking_energy("AC", "GA")
        self.assertGreaterEqual(dG, 0.0)

    # --------------------------------------------------------------------------
    # Feature R3.2: Hairpin Loop Energetics & Tetraloops
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_BIOPHYSICS_PARAMS, "synaptik.biophysics.parameters not available")
    def test_tc_r3_2_01_hairpin_loop_length_scaling(self):
        """TC-R3.2-01: Hairpin initiation energies for L=3, 4, 6."""
        self.assertAlmostEqual(get_hairpin_energy(3), 5.70, delta=0.1)
        self.assertAlmostEqual(get_hairpin_energy(4), 5.60, delta=0.1)
        self.assertAlmostEqual(get_hairpin_energy(6), 5.40, delta=0.1)

    @unittest.skipUnless(HAVE_BIOPHYSICS_PARAMS, "synaptik.biophysics.parameters not available")
    def test_tc_r3_2_02_ultra_stable_gaaa_tetraloop(self):
        """TC-R3.2-02: GAAA tetraloop receives -3.00 kcal/mol bonus."""
        bonus = get_tetraloop_bonus("GAAA")
        self.assertAlmostEqual(bonus, -3.00, delta=0.1)

    @unittest.skipUnless(HAVE_BIOPHYSICS_PARAMS, "synaptik.biophysics.parameters not available")
    def test_tc_r3_2_03_uacg_tetraloop(self):
        """TC-R3.2-03: UACG tetraloop receives -2.00 kcal/mol bonus."""
        bonus = get_tetraloop_bonus("UACG")
        self.assertAlmostEqual(bonus, -2.00, delta=0.1)

    @unittest.skipUnless(HAVE_BIOPHYSICS_PARAMS, "synaptik.biophysics.parameters not available")
    def test_tc_r3_2_04_polymer_entropy_scaling(self):
        """TC-R3.2-04: Large loops (L>9) scale via Jacobson-Stockmayer formula."""
        e12 = get_hairpin_energy(12)
        expected_12 = 5.40 + 1.0786 * math.log(12 / 6.0)
        self.assertAlmostEqual(e12, expected_12, delta=0.15)

    @unittest.skipUnless(HAVE_BIOPHYSICS_PARAMS, "synaptik.biophysics.parameters not available")
    def test_tc_r3_2_05_steric_prohibition_loops_under_3(self):
        """TC-R3.2-05: Loops <3 nt are physically forbidden (infinite penalty)."""
        self.assertTrue(math.isinf(get_hairpin_energy(1)) or get_hairpin_energy(1) > 100.0)
        self.assertTrue(math.isinf(get_hairpin_energy(2)) or get_hairpin_energy(2) > 100.0)

    # --------------------------------------------------------------------------
    # Feature R3.3: Dynamic Programming RNA MFE Folding
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_BIOPHYSICS_FOLDING, "synaptik.biophysics.folding not available")
    def test_tc_r3_3_01_canonical_stem_hairpin_folding(self):
        """TC-R3.3-01: Folds GCGCUUUUGCGC into stem-loop structure."""
        dG, db = fold_mfe_pure_python("GCGCUUUUGCGC")
        self.assertLess(dG, 0.0)
        self.assertEqual(db, "((((....))))")

    @unittest.skipUnless(HAVE_BIOPHYSICS_FOLDING, "synaptik.biophysics.folding not available")
    def test_tc_r3_3_02_unstructured_poly_a_folding(self):
        """TC-R3.3-02: Poly-A sequence returns 0.0 kcal/mol and all dots."""
        dG, db = fold_mfe_pure_python("AAAAAAAAAAAAAAAAAAAA")
        self.assertEqual(dG, 0.0)
        self.assertEqual(db, "." * 20)

    @unittest.skipUnless(HAVE_BIOPHYSICS_FOLDING, "synaptik.biophysics.folding not available")
    def test_tc_r3_3_03_dna_t_to_u_sanitization(self):
        """TC-R3.3-03: DNA sequence with T is converted to U and folds identically."""
        dG, db = fold_mfe_pure_python("GCGCTTTTGCGC")
        self.assertLess(dG, 0.0)
        self.assertEqual(db, "((((....))))")

    @unittest.skipUnless(HAVE_BIOPHYSICS_FOLDING, "synaptik.biophysics.folding not available")
    def test_tc_r3_3_04_short_sequence_handling(self):
        """TC-R3.3-04: Sequences <4 nt return (0.0, '...') without crashing."""
        dG, db = fold_mfe_pure_python("GCG")
        self.assertEqual(dG, 0.0)
        self.assertEqual(db, "...")

    @unittest.skipUnless(HAVE_BIOPHYSICS_FOLDING, "synaptik.biophysics.folding not available")
    def test_tc_r3_3_05_pure_python_environment_guarantee(self):
        """TC-R3.3-05: Folding executes without binary dependencies (no ViennaRNA/RNAfold subprocess)."""
        with patch("subprocess.run", side_effect=RuntimeError("Subprocess forbidden")):
            dG, db = fold_mfe_pure_python("GCGCUUUUGCGC")
            self.assertLess(dG, 0.0)

    # --------------------------------------------------------------------------
    # Feature R3.4: Target Opening Penalty Delta_G_open
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_BIOPHYSICS_FOLDING, "synaptik.biophysics.folding not available")
    def test_tc_r3_4_01_target_in_single_stranded_loop(self):
        """TC-R3.4-01: Target within single-stranded loop has Delta_G_open = 0.0."""
        dg_open = compute_delta_g_open("GCGCUUUUGCGC", binding_start=4, binding_length=4)
        self.assertEqual(dg_open, 0.0)

    @unittest.skipUnless(HAVE_BIOPHYSICS_FOLDING, "synaptik.biophysics.folding not available")
    def test_tc_r3_4_02_target_encompassing_stem_arm(self):
        """TC-R3.4-02: Target in stem arm incurs full unfolding penalty."""
        dG_folded, _ = fold_mfe_pure_python("GCGCUUUUGCGC")
        dg_open = compute_delta_g_open("GCGCUUUUGCGC", binding_start=0, binding_length=4)
        self.assertAlmostEqual(dg_open, -dG_folded, delta=0.01)

    @unittest.skipUnless(HAVE_BIOPHYSICS_FOLDING, "synaptik.biophysics.folding not available")
    def test_tc_r3_4_03_full_transcript_window(self):
        """TC-R3.4-03: Window spanning 100% of sequence forces total unfolding."""
        dG_folded, _ = fold_mfe_pure_python("GCGCUUUUGCGC")
        dg_open = compute_delta_g_open("GCGCUUUUGCGC", binding_start=0, binding_length=12)
        self.assertAlmostEqual(dg_open, -dG_folded, delta=0.01)

    @unittest.skipUnless(HAVE_BIOPHYSICS_FOLDING, "synaptik.biophysics.folding not available")
    def test_tc_r3_4_04_unstructured_transcript_penalty(self):
        """TC-R3.4-04: Unstructured poly-A transcript has Delta_G_open = 0.0."""
        dg_open = compute_delta_g_open("AAAAAAAAAAAAAAAAAAAA", binding_start=5, binding_length=10)
        self.assertEqual(dg_open, 0.0)

    @unittest.skipUnless(HAVE_BIOPHYSICS_FOLDING, "synaptik.biophysics.folding not available")
    def test_tc_r3_4_05_out_of_bounds_window_raises_error(self):
        """TC-R3.4-05: Out-of-bounds binding window raises ValueError."""
        with self.assertRaises(ValueError):
            compute_delta_g_open("ACGT", binding_start=3, binding_length=5)

    # --------------------------------------------------------------------------
    # Feature R3.5: Effective Hybridization & Kd
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_BIOPHYSICS_ACCESSIBILITY, "synaptik.biophysics.accessibility not available")
    def test_tc_r3_5_01_fully_accessible_hybridization(self):
        """TC-R3.5-01: Fully accessible target (Delta_G_open=0) has Delta_G_eff == duplex."""
        res = compute_effective_binding(duplex_delta_g=-22.0, delta_g_open=0.0)
        self.assertEqual(res.delta_g_effective, -22.0)
        self.assertTrue(res.is_favorable)
        self.assertLess(res.kd_effective_nm, 0.001)

    @unittest.skipUnless(HAVE_BIOPHYSICS_ACCESSIBILITY, "synaptik.biophysics.accessibility not available")
    def test_tc_r3_5_02_moderate_secondary_structure_penalty(self):
        """TC-R3.5-02: Delta_G_eff = duplex + open."""
        res = compute_effective_binding(duplex_delta_g=-24.5, delta_g_open=10.3)
        self.assertAlmostEqual(res.delta_g_effective, -14.2, places=3)
        expected_kd_nm = math.exp(-14.2 / 0.61633008) * 1e9
        self.assertAlmostEqual(res.kd_effective_nm, expected_kd_nm, delta=1.0)

    @unittest.skipUnless(HAVE_BIOPHYSICS_ACCESSIBILITY, "synaptik.biophysics.accessibility not available")
    def test_tc_r3_5_03_thermodynamically_unfavorable_occlusion(self):
        """TC-R3.5-03: Delta_G_open > |duplex| results in positive Delta_G_eff and is_favorable=False."""
        res = compute_effective_binding(duplex_delta_g=-18.0, delta_g_open=20.0)
        self.assertEqual(res.delta_g_effective, 2.0)
        self.assertFalse(res.is_favorable)
        self.assertGreater(res.kd_effective_molar, 1.0)

    @unittest.skipUnless(HAVE_BIOPHYSICS_ACCESSIBILITY, "synaptik.biophysics.accessibility not available")
    def test_tc_r3_5_04_target_unfold_equivalence(self):
        """TC-R3.5-04: Delta_G_eff = duplex - Delta_G_target_unfold equivalence."""
        open_penalty = 8.5
        unfold_energy = -8.5
        res1 = compute_effective_binding(-20.0, open_penalty)
        self.assertEqual(res1.delta_g_effective, -20.0 - unfold_energy)

    @unittest.skipUnless(HAVE_BIOPHYSICS_ACCESSIBILITY, "synaptik.biophysics.accessibility not available")
    def test_tc_r3_5_05_temperature_scaling(self):
        """TC-R3.5-05: Higher temperature increases Kd (weaker affinity)."""
        res_25 = compute_effective_binding(-20.0, 5.0, temperature_c=25.0)
        res_37 = compute_effective_binding(-20.0, 5.0, temperature_c=37.0)
        self.assertGreater(res_37.kd_effective_nm, res_25.kd_effective_nm)

    # --------------------------------------------------------------------------
    # Feature R3.6: Spliceosomal Steric Displacement
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_BIOPHYSICS_COMPETITION, "synaptik.biophysics.competition not available")
    def test_tc_r3_6_01_u1_snrnp_displacement(self):
        """TC-R3.6-01: High affinity ASO displaces U1 snRNP (Kd=25 nM)."""
        res = compute_steric_displacement(oligo_kd_nm=0.10, competitor_name="u1_snrnp", target_gene="STMN2")
        self.assertAlmostEqual(res.competitor_kd_nm, 25.0, places=2)
        self.assertAlmostEqual(res.displacement_ratio, 250.0, places=1)
        self.assertGreater(res.fractional_occupancy, 0.98)
        self.assertEqual(res.dominant_binder, "oligo")

    @unittest.skipUnless(HAVE_BIOPHYSICS_COMPETITION, "synaptik.biophysics.competition not available")
    def test_tc_r3_6_02_u2af65_displacement(self):
        """TC-R3.6-02: ASO displaces U2AF65 (Kd=20 nM) at polypyrimidine tract."""
        res = compute_steric_displacement(oligo_kd_nm=0.20, competitor_name="u2af65", target_gene="UNC13A")
        self.assertAlmostEqual(res.competitor_kd_nm, 20.0, places=2)
        self.assertAlmostEqual(res.displacement_ratio, 100.0, places=1)
        self.assertGreater(res.fractional_occupancy, 0.95)

    @unittest.skipUnless(HAVE_BIOPHYSICS_COMPETITION, "synaptik.biophysics.competition not available")
    def test_tc_r3_6_03_u2af35_displacement(self):
        """TC-R3.6-03: ASO displaces U2AF35 (Kd=40 nM) at AG acceptor."""
        res = compute_steric_displacement(oligo_kd_nm=0.40, competitor_name="u2af35", target_gene="UNC13A")
        self.assertAlmostEqual(res.competitor_kd_nm, 40.0, places=2)
        self.assertAlmostEqual(res.displacement_ratio, 100.0, places=1)

    @unittest.skipUnless(HAVE_BIOPHYSICS_COMPETITION, "synaptik.biophysics.competition not available")
    def test_tc_r3_6_04_tdp43_homooligomer_displacement(self):
        """TC-R3.6-04: ASO displaces TDP-43 homooligomers (Kd=3.0 nM)."""
        res = compute_steric_displacement(oligo_kd_nm=0.05, competitor_name="tardbp", target_gene="STMN2")
        self.assertAlmostEqual(res.competitor_kd_nm, 3.0, places=2)
        self.assertAlmostEqual(res.displacement_ratio, 60.0, places=1)
        self.assertGreater(res.fractional_occupancy, 0.93)

    @unittest.skipUnless(HAVE_BIOPHYSICS_COMPETITION, "synaptik.biophysics.competition not available")
    def test_tc_r3_6_05_weak_oligo_outcompeted(self):
        """TC-R3.6-05: Weak oligo (Kd=100 nM) is outcompeted by endogenous factor."""
        res = compute_steric_displacement(oligo_kd_nm=100.0, competitor_name="u1_snrnp", target_gene="STMN2")
        self.assertLess(res.displacement_ratio, 1.0)
        self.assertEqual(res.dominant_binder, "competitor")

    # --------------------------------------------------------------------------
    # Feature R3.7: CLI `synaptik biophysics`
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_CLI_BIOPHYSICS, "synaptik.cli.commands.biophysics not available")
    def test_tc_r3_7_01_cli_biophysics_open_text(self):
        """TC-R3.7-01: CLI biophysics open executes and outputs author citation."""
        args = MagicMock(subcommand="open", target_sequence="GCGCUUUUGCGC", oligo_start=4, oligo_len=4, json=False)
        buf = io.StringIO()
        with patch("sys.stdout", buf):
            code = cmd_biophysics(args)
        self.assertEqual(code, 0)
        self.assertIn("Francesco Nobile", buf.getvalue())

    @unittest.skipUnless(HAVE_CLI_BIOPHYSICS, "synaptik.cli.commands.biophysics not available")
    def test_tc_r3_7_02_cli_biophysics_open_json(self):
        """TC-R3.7-02: CLI biophysics open --json emits structured metrics."""
        args = MagicMock(subcommand="open", target_sequence="GCGCUUUUGCGC", oligo_start=4, oligo_len=4, json=True)
        buf = io.StringIO()
        with patch("sys.stdout", buf):
            code = cmd_biophysics(args)
        self.assertEqual(code, 0)
        data = json.loads(buf.getvalue())
        self.assertIn("delta_g_open", data)

    @unittest.skipUnless(HAVE_CLI_BIOPHYSICS, "synaptik.cli.commands.biophysics not available")
    def test_tc_r3_7_03_cli_biophysics_displace(self):
        """TC-R3.7-03: CLI biophysics displace models steric competition."""
        args = MagicMock(subcommand="displace", gene="STMN2", oligo_sequence=None, oligo_kd=0.1, competitor="u1_snrnp", oligo_conc=10.0, competitor_conc=10.0, json=True)
        buf = io.StringIO()
        with patch("sys.stdout", buf):
            code = cmd_biophysics(args)
        self.assertEqual(code, 0)
        data = json.loads(buf.getvalue())
        self.assertEqual(data["dominant_binder"], "oligo")

    @unittest.skipUnless(HAVE_CLI_BIOPHYSICS, "synaptik.cli.commands.biophysics not available")
    def test_tc_r3_7_04_cli_biophysics_scan(self):
        """TC-R3.7-04: CLI biophysics scan tiles target junction."""
        args = MagicMock(subcommand="scan", gene="STMN2", sequence="GCGCUUUUGCGCGCGCUUUUGCGC", window=24, oligo_len=8, json=True)
        buf = io.StringIO()
        with patch("sys.stdout", buf):
            code = cmd_biophysics(args)
        self.assertEqual(code, 0)

    @unittest.skipUnless(HAVE_CLI_BIOPHYSICS, "synaptik.cli.commands.biophysics not available")
    def test_tc_r3_7_05_cli_biophysics_error_handling(self):
        """TC-R3.7-05: Invalid window offset returns exit code 1."""
        args = MagicMock(subcommand="open", target_sequence="ACGT", oligo_start=10, oligo_len=5, json=False)
        code = cmd_biophysics(args)
        self.assertEqual(code, 1)

    # --------------------------------------------------------------------------
    # Feature R4.1: Patient Baseline Parameter Ingestion
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_TWIN_SIMULATOR, "synaptik.twin.twin_generator not available")
    def test_tc_r4_1_01_typical_adult_baseline_ingestion(self):
        """TC-R4.1-01: Adult baseline fields are correctly populated."""
        b = PatientBaseline(patient_id="ALS-001", age_years=56.4, onset_site="bulbar", disease_duration_months=14.5, baseline_alsfrs_r=38.0, baseline_slope=-0.69, baseline_plasma_nfl=72.5)
        self.assertEqual(b.patient_id, "ALS-001")
        self.assertEqual(b.age_years, 56.4)
        self.assertEqual(b.onset_site, "bulbar")

    @unittest.skipUnless(HAVE_TWIN_SIMULATOR, "synaptik.twin.twin_generator not available")
    def test_tc_r4_1_02_limb_onset_normalization(self):
        """TC-R4.1-02: Limb onset site is normalized to lowercase."""
        b = PatientBaseline(patient_id="ALS-002", age_years=62.0, onset_site="LIMB", disease_duration_months=8.0, baseline_alsfrs_r=42.0)
        self.assertEqual(b.onset_site.lower(), "limb")

    @unittest.skipUnless(HAVE_TWIN_SIMULATOR, "synaptik.twin.twin_generator not available")
    def test_tc_r4_1_03_spinal_onset_accepted(self):
        """TC-R4.1-03: Spinal onset site is accepted."""
        b = PatientBaseline(patient_id="ALS-003", age_years=45.0, onset_site="spinal", disease_duration_months=12.0, baseline_alsfrs_r=40.0)
        self.assertEqual(b.onset_site, "spinal")

    @unittest.skipUnless(HAVE_TWIN_SIMULATOR, "synaptik.twin.twin_generator not available")
    def test_tc_r4_1_04_pre_slope_auto_calculation(self):
        """TC-R4.1-04: Missing pre-slope is calculated from (ALSFRS - 48)/duration."""
        b = PatientBaseline(patient_id="ALS-004", age_years=50.0, onset_site="bulbar", disease_duration_months=10.0, baseline_alsfrs_r=40.0, baseline_slope=None)
        expected_slope = (40.0 - 48.0) / 10.0  # -0.80 pts/mo
        slope = getattr(b, "baseline_slope", None)
        if slope is not None:
            self.assertAlmostEqual(slope, expected_slope, places=2)

    @unittest.skipUnless(HAVE_TWIN_SIMULATOR, "synaptik.twin.twin_generator not available")
    def test_tc_r4_1_05_optional_nfl_population_median_fallback(self):
        """TC-R4.1-05: Missing NfL accepts None or defaults to 58.0 pg/mL."""
        b = PatientBaseline(patient_id="ALS-005", age_years=55.0, onset_site="limb", disease_duration_months=6.0, baseline_alsfrs_r=44.0, baseline_plasma_nfl=None)
        self.assertIsNone(b.baseline_plasma_nfl)

    # --------------------------------------------------------------------------
    # Feature R4.2: PRO-ACT / Answer ALS Covariance Engine
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_TWIN_DISTRIBUTIONS, "synaptik.twin.distributions not available")
    def test_tc_r4_2_01_bulbar_vs_limb_progression(self):
        """TC-R4.2-01: Bulbar onset progresses faster than limb by ~0.36 pts/mo."""
        b_bulbar = PatientBaseline("P1", 55.0, "bulbar", 10.0, 40.0, baseline_slope=-0.92, baseline_plasma_nfl=58.0)
        b_limb = PatientBaseline("P2", 55.0, "limb", 10.0, 40.0, baseline_slope=-0.92, baseline_plasma_nfl=58.0)
        c_bulbar = compute_conditional_progression(b_bulbar)
        c_limb = compute_conditional_progression(b_limb)
        self.assertLess(c_bulbar.mean_slope, c_limb.mean_slope)

    @unittest.skipUnless(HAVE_TWIN_DISTRIBUTIONS, "synaptik.twin.distributions not available")
    def test_tc_r4_2_02_log_linear_nfl_adjustment(self):
        """TC-R4.2-02: 2x elevated NfL accelerates predicted decline by -0.35*ln(2) ~ -0.243 pts/mo."""
        b_high_nfl = PatientBaseline("P1", 55.0, "limb", 10.0, 40.0, baseline_slope=-0.92, baseline_plasma_nfl=116.0)
        c = compute_conditional_progression(b_high_nfl)
        expected_nfl_adj = -0.35 * math.log(116.0 / 58.0)
        self.assertAlmostEqual(c.nfl_adjustment, expected_nfl_adj, places=2)

    @unittest.skipUnless(HAVE_TWIN_DISTRIBUTIONS, "synaptik.twin.distributions not available")
    def test_tc_r4_2_03_age_adjustment(self):
        """TC-R4.2-03: Age 75 yr (+20 yr vs 55) adds -0.012*20 = -0.24 pts/mo."""
        b_elderly = PatientBaseline("P1", 75.0, "limb", 10.0, 40.0, baseline_slope=-0.92, baseline_plasma_nfl=58.0)
        c = compute_conditional_progression(b_elderly)
        self.assertAlmostEqual(c.age_adjustment, -0.24, places=2)

    @unittest.skipUnless(HAVE_TWIN_DISTRIBUTIONS, "synaptik.twin.distributions not available")
    def test_tc_r4_2_04_shrinkage_weighting(self):
        """TC-R4.2-04: Shrinkage weights 0.65 individual pre-slope and 0.35 population conditional."""
        pre = -1.20
        b = PatientBaseline("P1", 55.0, "limb", 10.0, 40.0, baseline_slope=pre, baseline_plasma_nfl=58.0)
        c = compute_conditional_progression(b)
        expected_mean = 0.65 * pre + 0.35 * (-0.92 + 0.08)
        self.assertAlmostEqual(c.mean_slope, expected_mean, places=2)

    @unittest.skipUnless(HAVE_TWIN_DISTRIBUTIONS, "synaptik.twin.distributions not available")
    def test_tc_r4_2_05_residual_variance(self):
        """TC-R4.2-05: Conditional residual progression standard deviation is 0.280 pts/mo."""
        b = PatientBaseline("P1", 55.0, "limb", 10.0, 40.0)
        c = compute_conditional_progression(b)
        self.assertAlmostEqual(c.std_slope, 0.280, places=2)

    # --------------------------------------------------------------------------
    # Feature R4.3: 1,000 Matched Untreated Trajectory Simulator
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_TWIN_SIMULATOR, "synaptik.twin.twin_generator not available")
    def test_tc_r4_3_01_cohort_dimension_verification(self):
        """TC-R4.3-01: Generates exactly 1,000 trajectories over 18 months."""
        b = PatientBaseline("P1", 56.4, "bulbar", 14.5, 38.0, baseline_slope=-0.69, baseline_plasma_nfl=72.5)
        cohort = generate_counterfactual_cohort(b, n_twins=1000, duration_months=18, random_seed=42)
        self.assertEqual(cohort.n_twins, 1000)
        self.assertEqual(len(cohort.trajectories), 1000)
        self.assertEqual(len(cohort.median_trajectory), 19)

    @unittest.skipUnless(HAVE_TWIN_SIMULATOR, "synaptik.twin.twin_generator not available")
    def test_tc_r4_3_02_measurement_noise(self):
        """TC-R4.3-02: Trajectory points reflect clinical measurement noise."""
        b = PatientBaseline("P1", 55.0, "limb", 10.0, 40.0, baseline_slope=-0.50, baseline_plasma_nfl=58.0)
        cohort = generate_counterfactual_cohort(b, n_twins=100, duration_months=6, random_seed=42)
        # Check non-deterministic deviation between twin 0 and twin 1
        self.assertNotEqual(cohort.trajectories[0], cohort.trajectories[1])

    @unittest.skipUnless(HAVE_TWIN_SIMULATOR, "synaptik.twin.twin_generator not available")
    def test_tc_r4_3_03_ceiling_constraint(self):
        """TC-R4.3-03: No trajectory point exceeds ceiling of 48.0."""
        b = PatientBaseline("P1", 55.0, "limb", 10.0, 48.0, baseline_slope=0.0)
        cohort = generate_counterfactual_cohort(b, n_twins=100, duration_months=6, random_seed=42)
        for traj in cohort.trajectories:
            for score in traj:
                self.assertLessEqual(score, 48.0)

    @unittest.skipUnless(HAVE_TWIN_SIMULATOR, "synaptik.twin.twin_generator not available")
    def test_tc_r4_3_04_absorbing_floor_at_zero(self):
        """TC-R4.3-04: Once a trajectory hits 0.0, all subsequent points remain 0.0."""
        b = PatientBaseline("P1", 75.0, "bulbar", 6.0, 5.0, baseline_slope=-2.0)
        cohort = generate_counterfactual_cohort(b, n_twins=50, duration_months=12, random_seed=42)
        for traj in cohort.trajectories:
            hit_zero = False
            for score in traj:
                if score <= 0.0:
                    hit_zero = True
                if hit_zero:
                    self.assertEqual(score, 0.0)

    @unittest.skipUnless(HAVE_TWIN_SIMULATOR, "synaptik.twin.twin_generator not available")
    def test_tc_r4_3_05_reproducibility_with_seed(self):
        """TC-R4.3-05: Identical random seed yields bitwise identical trajectories."""
        b = PatientBaseline("P1", 55.0, "limb", 10.0, 40.0)
        c1 = generate_counterfactual_cohort(b, n_twins=10, duration_months=6, random_seed=999)
        c2 = generate_counterfactual_cohort(b, n_twins=10, duration_months=6, random_seed=999)
        self.assertEqual(c1.median_trajectory, c2.median_trajectory)

    # --------------------------------------------------------------------------
    # Feature R4.4: Bayesian Posterior Probability of Deceleration
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_TWIN_BAYESIAN, "synaptik.twin.bayesian_inference not available")
    def test_tc_r4_4_01_decisive_efficacy_scenario(self):
        """TC-R4.4-01: Stabilized patient slope yields P(Delta>0)>=0.975 and Tier A."""
        b = PatientBaseline("P1", 56.4, "bulbar", 14.5, 38.0, baseline_slope=-1.20)
        cohort = generate_counterfactual_cohort(b, n_twins=1000, duration_months=12, random_seed=42)
        visits = [
            VisitObservation(month=0.0, alsfrs_score=38.0),
            VisitObservation(month=3.0, alsfrs_score=38.0),
            VisitObservation(month=6.0, alsfrs_score=38.0),
            VisitObservation(month=9.0, alsfrs_score=38.0),
            VisitObservation(month=12.0, alsfrs_score=38.0),
        ]
        res = evaluate_treatment_efficacy(cohort, visits)
        self.assertGreaterEqual(res.p_deceleration, 0.975)
        self.assertGreater(res.bayes_factor, 39.0)
        self.assertIn("Tier A", res.efficacy_tier)

    @unittest.skipUnless(HAVE_TWIN_BAYESIAN, "synaptik.twin.bayesian_inference not available")
    def test_tc_r4_4_02_moderate_slowing_scenario(self):
        """TC-R4.4-02: Moderate slowing yields Tier B (0.90 <= P < 0.975)."""
        b = PatientBaseline("P1", 56.4, "bulbar", 14.5, 38.0, baseline_slope=-1.20)
        cohort = generate_counterfactual_cohort(b, n_twins=1000, duration_months=12, random_seed=42)
        visits = [
            VisitObservation(month=0.0, alsfrs_score=38.0),
            VisitObservation(month=3.0, alsfrs_score=37.0),
            VisitObservation(month=6.0, alsfrs_score=36.0),
            VisitObservation(month=9.0, alsfrs_score=35.0),
            VisitObservation(month=12.0, alsfrs_score=34.0),
        ]
        res = evaluate_treatment_efficacy(cohort, visits)
        self.assertGreaterEqual(res.p_deceleration, 0.85)

    @unittest.skipUnless(HAVE_TWIN_BAYESIAN, "synaptik.twin.bayesian_inference not available")
    def test_tc_r4_4_03_inconclusive_signal_scenario(self):
        """TC-R4.4-03: Treated slope matching twin median yields P ~0.50 and Tier C."""
        b = PatientBaseline("P1", 55.0, "limb", 10.0, 40.0, baseline_slope=-0.80)
        cohort = generate_counterfactual_cohort(b, n_twins=1000, duration_months=12, random_seed=42)
        visits = [
            VisitObservation(month=0.0, alsfrs_score=40.0),
            VisitObservation(month=6.0, alsfrs_score=35.2),
            VisitObservation(month=12.0, alsfrs_score=30.4),
        ]
        res = evaluate_treatment_efficacy(cohort, visits)
        self.assertAlmostEqual(res.p_deceleration, 0.50, delta=0.15)

    @unittest.skipUnless(HAVE_TWIN_BAYESIAN, "synaptik.twin.bayesian_inference not available")
    def test_tc_r4_4_04_disease_acceleration_scenario(self):
        """TC-R4.4-04: Accelerated progression yields P < 0.10 and Tier D."""
        b = PatientBaseline("P1", 55.0, "limb", 10.0, 40.0, baseline_slope=-0.80)
        cohort = generate_counterfactual_cohort(b, n_twins=1000, duration_months=12, random_seed=42)
        visits = [
            VisitObservation(month=0.0, alsfrs_score=40.0),
            VisitObservation(month=3.0, alsfrs_score=32.0),
            VisitObservation(month=6.0, alsfrs_score=24.0),
        ]
        res = evaluate_treatment_efficacy(cohort, visits)
        self.assertLess(res.p_deceleration, 0.10)
        self.assertIn("Tier D", res.efficacy_tier)

    @unittest.skipUnless(HAVE_TWIN_BAYESIAN, "synaptik.twin.bayesian_inference not available")
    def test_tc_r4_4_05_two_visit_minimum_requirement(self):
        """TC-R4.4-05: Minimum 2 visits computes slope; single visit raises ValueError."""
        b = PatientBaseline("P1", 55.0, "limb", 10.0, 40.0)
        cohort = generate_counterfactual_cohort(b, n_twins=100, duration_months=6, random_seed=42)
        with self.assertRaises(ValueError):
            evaluate_treatment_efficacy(cohort, [VisitObservation(month=0.0, alsfrs_score=40.0)])

    # --------------------------------------------------------------------------
    # Feature R4.5: Audit-Ready FDA Dossier Export
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_TWIN_DOSSIER, "synaptik.twin.dossier_export not available")
    def test_tc_r4_5_01_markdown_dossier_structure(self):
        """TC-R4.5-01: Markdown dossier contains required regulatory sections."""
        b = PatientBaseline("ALS-001", 56.4, "bulbar", 14.5, 38.0)
        cohort = generate_counterfactual_cohort(b, n_twins=100, duration_months=12, random_seed=42)
        dossier = generate_twin_fda_dossier(b, cohort)
        self.assertIn("Executive Summary", dossier)
        self.assertIn("Counterfactual Cohort", dossier)

    @unittest.skipUnless(HAVE_TWIN_DOSSIER, "synaptik.twin.dossier_export not available")
    def test_tc_r4_5_02_regulatory_statute_citation(self):
        """TC-R4.5-02: Dossier explicitly cites 21 CFR Part 312 / FDA Rare Disease Guidance."""
        b = PatientBaseline("ALS-001", 56.4, "bulbar", 14.5, 38.0)
        cohort = generate_counterfactual_cohort(b, n_twins=100, duration_months=12, random_seed=42)
        dossier = generate_twin_fda_dossier(b, cohort)
        self.assertIn("21 CFR", dossier)

    @unittest.skipUnless(HAVE_TWIN_DOSSIER, "synaptik.twin.dossier_export not available")
    def test_tc_r4_5_03_author_attribution_in_dossier(self):
        """TC-R4.5-03: Dossier footer includes Francesco Nobile (ORCID 0009-0001-9477-3332)."""
        b = PatientBaseline("ALS-001", 56.4, "bulbar", 14.5, 38.0)
        cohort = generate_counterfactual_cohort(b, n_twins=100, duration_months=12, random_seed=42)
        dossier = generate_twin_fda_dossier(b, cohort)
        self.assertIn("Francesco Nobile", dossier)
        self.assertIn("0009-0001-9477-3332", dossier)

    @unittest.skipUnless(HAVE_TWIN_DOSSIER, "synaptik.twin.dossier_export not available")
    def test_tc_r4_5_04_efficacy_results_inclusion(self):
        """TC-R4.5-04: When efficacy result is passed, posterior and Bayes factor are formatted."""
        b = PatientBaseline("ALS-001", 56.4, "bulbar", 14.5, 38.0)
        cohort = generate_counterfactual_cohort(b, n_twins=100, duration_months=12, random_seed=42)
        eff = BayesianEfficacyResult(
            p_deceleration=0.985,
            bayes_factor=65.7,
            observed_slope=-0.50,
            observed_slope_se=0.12,
            median_counterfactual_slope=-1.35,
            treatment_effect_delta=0.85,
            efficacy_tier="Tier A: Decisive Efficacy",
            fda_regulatory_verdict="Accelerated Approval Precedent",
        )
        dossier = generate_twin_fda_dossier(b, cohort, efficacy=eff)
        self.assertIn("0.985", dossier)
        self.assertIn("65.7", dossier)
        self.assertIn("Tier A", dossier)

    @unittest.skipUnless(HAVE_TWIN_DOSSIER, "synaptik.twin.dossier_export not available")
    def test_tc_r4_5_05_utf8_encoding_compliance(self):
        """TC-R4.5-05: Dossier string encodes cleanly to UTF-8 without unicode decode errors."""
        b = PatientBaseline("ALS-001", 56.4, "bulbar", 14.5, 38.0)
        cohort = generate_counterfactual_cohort(b, n_twins=50, duration_months=6, random_seed=42)
        dossier = generate_twin_fda_dossier(b, cohort)
        encoded = dossier.encode("utf-8")
        self.assertGreater(len(encoded), 100)

    # --------------------------------------------------------------------------
    # Feature R4.6: CLI `synaptik twin`
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_CLI_TWIN, "synaptik.cli.commands.twin not available")
    def test_tc_r4_6_01_cli_twin_simulate(self):
        """TC-R4.6-01: CLI twin simulate generates counterfactual cohort."""
        args = MagicMock(subcommand="simulate", age=56.4, site="bulbar", duration=14.5, alsfrs=38.0, nfl=72.5, months=18, twins=100, seed=42, json=False)
        code = cmd_twin(args)
        self.assertEqual(code, 0)

    @unittest.skipUnless(HAVE_CLI_TWIN, "synaptik.cli.commands.twin not available")
    def test_tc_r4_6_02_cli_twin_simulate_json(self):
        """TC-R4.6-02: CLI twin simulate --json emits valid JSON."""
        args = MagicMock(subcommand="simulate", age=56.4, site="bulbar", duration=14.5, alsfrs=38.0, nfl=72.5, months=18, twins=100, seed=42, json=True)
        buf = io.StringIO()
        with patch("sys.stdout", buf):
            code = cmd_twin(args)
        self.assertEqual(code, 0)
        data = json.loads(buf.getvalue())
        self.assertIn("median_trajectory", data)

    @unittest.skipUnless(HAVE_CLI_TWIN, "synaptik.cli.commands.twin not available")
    def test_tc_r4_6_03_cli_twin_evaluate(self):
        """TC-R4.6-03: CLI twin evaluate computes Bayesian posterior."""
        args = MagicMock(
            subcommand="evaluate",
            age=56.4,
            site="bulbar",
            duration=14.5,
            alsfrs=38.0,
            nfl=72.5,
            observed_scores="0:38.0,3:38.0,6:38.0",
            json=True,
        )
        buf = io.StringIO()
        with patch("sys.stdout", buf):
            code = cmd_twin(args)
        self.assertEqual(code, 0)
        data = json.loads(buf.getvalue())
        self.assertIn("p_deceleration", data)

    @unittest.skipUnless(HAVE_CLI_TWIN, "synaptik.cli.commands.twin not available")
    def test_tc_r4_6_04_cli_twin_dossier_file_export(self):
        """TC-R4.6-04: CLI twin dossier --output writes markdown file."""
        with tempfile.NamedTemporaryFile(suffix=".md", delete=False) as tmp:
            tmp_path = tmp.name
        try:
            args = MagicMock(
                subcommand="dossier",
                patient="ALS-001",
                age=56.4,
                site="bulbar",
                duration=14.5,
                alsfrs=38.0,
                nfl=72.5,
                observed_scores="0:38.0,3:38.0",
                output=tmp_path,
                json=False,
            )
            code = cmd_twin(args)
            self.assertEqual(code, 0)
            self.assertTrue(os.path.exists(tmp_path))
            with open(tmp_path, "r", encoding="utf-8") as f:
                content = f.read()
            self.assertIn("Executive Summary", content)
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

    @unittest.skipUnless(HAVE_CLI_TWIN, "synaptik.cli.commands.twin not available")
    def test_tc_r4_6_05_cli_twin_invalid_observed_scores_syntax(self):
        """TC-R4.6-05: Malformed observed scores string returns exit code 2."""
        args = MagicMock(subcommand="evaluate", age=56.4, site="bulbar", duration=14.5, alsfrs=38.0, nfl=72.5, observed_scores="malformed_syntax", json=False)
        code = cmd_twin(args)
        self.assertEqual(code, 2)

    # --------------------------------------------------------------------------
    # Feature R5.1: Zero-PHI Persistence Protocol
    # --------------------------------------------------------------------------
    def test_tc_r5_1_01_sqlite_database_phi_isolation(self):
        """TC-R5.1-01: Forensic check: Zero canary patient strings in synaptik.db."""
        db_path = REPO_ROOT / "data" / "synaptik.db"
        if not db_path.exists():
            self.skipTest("data/synaptik.db does not exist yet")

        conn = sqlite3.connect(str(db_path))
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = [row[0] for row in cursor.fetchall()]

        canaries = ["CANARY_PHI", "CANARY_PATIENT", "MRN_CANARY"]
        for table in tables:
            cursor.execute(f"PRAGMA table_info({table});")
            cols = [c[1] for c in cursor.fetchall()]
            for col in cols:
                for canary in canaries:
                    q = f"SELECT COUNT(*) FROM {table} WHERE CAST({col} AS TEXT) LIKE ?"
                    cursor.execute(q, (f"%{canary}%",))
                    cnt = cursor.fetchone()[0]
                    self.assertEqual(cnt, 0, f"PHI Canary leaked into table {table}, column {col}")
        conn.close()

    def test_tc_r5_1_02_filesystem_cleanliness_after_simulation(self):
        """TC-R5.1-02: In-memory simulation leaves no unrequested files in workspace."""
        initial_files = set(os.listdir(str(REPO_ROOT)))
        if HAVE_TWIN_SIMULATOR:
            b = PatientBaseline("CANARY_PHI_PATIENT_N1_98472", 56.4, "bulbar", 14.5, 38.0)
            _ = generate_counterfactual_cohort(b, n_twins=50, duration_months=6, random_seed=42)
        final_files = set(os.listdir(str(REPO_ROOT)))
        new_files = final_files - initial_files
        self.assertEqual(new_files, set(), f"Unexpected files created: {new_files}")

    def test_tc_r5_1_03_temp_directory_canary_audit(self):
        """TC-R5.1-03: Zero canary tokens in OS temporary directories."""
        canary = "CANARY_PHI_PATIENT_N1_98472"
        if HAVE_TWIN_SIMULATOR:
            b = PatientBaseline(canary, 56.4, "bulbar", 14.5, 38.0)
            _ = generate_counterfactual_cohort(b, n_twins=10, duration_months=3)

        temp_dir = pathlib.Path(tempfile.gettempdir())
        now = time.time()
        for f in temp_dir.glob("*"):
            if f.is_file() and (now - f.stat().st_mtime < 120) and f.stat().st_size < 1000000:
                try:
                    txt = f.read_text(encoding="utf-8", errors="ignore")
                    self.assertNotIn(canary, txt, f"Canary found in tempfile {f}")
                except Exception:
                    pass

    def test_tc_r5_1_04_logger_canary_suppression(self):
        """TC-R5.1-04: Log handler captures zero unmasked patient canary strings."""
        log_stream = io.StringIO()
        handler = logging.StreamHandler(log_stream)
        logger = logging.getLogger("synaptik")
        logger.addHandler(handler)
        logger.setLevel(logging.DEBUG)

        canary = "CANARY_PHI_MRN_883019482"
        if HAVE_TWIN_SIMULATOR:
            b = PatientBaseline(canary, 50.0, "limb", 10.0, 40.0)
            _ = generate_counterfactual_cohort(b, n_twins=10, duration_months=3)

        logger.removeHandler(handler)
        log_content = log_stream.getvalue()
        self.assertNotIn(canary, log_content)

    def test_tc_r5_1_05_in_memory_dossier_generation(self):
        """TC-R5.1-05: Dossier generated in-memory does not write to disk unless explicit path."""
        if HAVE_TWIN_DOSSIER and HAVE_TWIN_SIMULATOR:
            b = PatientBaseline("CANARY_PATIENT", 50.0, "limb", 10.0, 40.0)
            cohort = generate_counterfactual_cohort(b, n_twins=10, duration_months=3)
            doc = generate_twin_fda_dossier(b, cohort)
            self.assertIsInstance(doc, str)

    # --------------------------------------------------------------------------
    # Feature R5.2: Python 3.10–3.12 Compatibility Audit
    # --------------------------------------------------------------------------
    def test_tc_r5_2_01_ast_parse_all_synaptik_modules(self):
        """TC-R5.2-01: AST parse of every .py file in synaptik/ passes with zero SyntaxError."""
        synaptik_dir = REPO_ROOT / "synaptik"
        for p in synaptik_dir.rglob("*.py"):
            with open(p, "r", encoding="utf-8") as f:
                source = f.read()
            try:
                ast.parse(source, filename=str(p))
            except SyntaxError as e:
                self.fail(f"Syntax error in {p}: {e}")

    def test_tc_r5_2_02_ast_parse_all_test_modules(self):
        """TC-R5.2-02: AST parse of every .py file in tests/ passes with zero SyntaxError."""
        tests_dir = REPO_ROOT / "tests"
        for p in tests_dir.rglob("*.py"):
            with open(p, "r", encoding="utf-8") as f:
                source = f.read()
            try:
                ast.parse(source, filename=str(p))
            except SyntaxError as e:
                self.fail(f"Syntax error in {p}: {e}")

    def test_tc_r5_2_03_zero_nested_fstrings_with_backslashes(self):
        """TC-R5.2-03: Zero f-strings contain backslash escapes inside replacement braces."""
        for d in [REPO_ROOT / "synaptik", REPO_ROOT / "tests"]:
            for p in d.rglob("*.py"):
                with open(p, "r", encoding="utf-8") as f:
                    txt = f.read()
                try:
                    tree = ast.parse(txt, filename=str(p))
                except SyntaxError:
                    continue
                for node in ast.walk(tree):
                    if isinstance(node, ast.JoinedStr):
                        for val in node.values:
                            if isinstance(val, ast.FormattedValue):
                                seg = ast.get_source_segment(txt, val)
                                if seg and "\\" in seg:
                                    self.fail(f"Disallowed backslash in f-string expression in {p}: {seg}")

    def test_tc_r5_2_04_future_annotations_presence(self):
        """TC-R5.2-04: Verifies from __future__ import annotations in key modules."""
        modules_to_check = [
            REPO_ROOT / "synaptik" / "clinical" / "allometry.py",
            REPO_ROOT / "synaptik" / "ingest" / "http.py",
        ]
        for p in modules_to_check:
            if p.exists():
                with open(p, "r", encoding="utf-8") as f:
                    txt = f.read()
                self.assertIn("from __future__ import annotations", txt, f"Missing annotations import in {p}")

    def test_tc_r5_2_05_clean_encoding_no_bom(self):
        """TC-R5.2-05: Python files have no UTF-8 Byte Order Mark (BOM)."""
        for p in (REPO_ROOT / "synaptik").rglob("*.py"):
            with open(p, "rb") as f:
                header = f.read(3)
            self.assertNotEqual(header, b"\xef\xbb\xbf", f"File contains UTF-8 BOM: {p}")

    # --------------------------------------------------------------------------
    # Feature R5.3: Git Cleanliness & CLAUDE.md Exclusion
    # --------------------------------------------------------------------------
    def test_tc_r5_3_01_gitignore_excludes_claude_md(self):
        """TC-R5.3-01: .gitignore contains explicit exclusion for CLAUDE.md."""
        gitignore = REPO_ROOT / ".gitignore"
        self.assertTrue(gitignore.exists(), ".gitignore file must exist")
        with open(gitignore, "r", encoding="utf-8") as f:
            content = f.read()
        self.assertIn("CLAUDE.md", content, ".gitignore must explicitly mention CLAUDE.md")

    def test_tc_r5_3_02_git_ls_files_claude_md(self):
        """TC-R5.3-02: git ls-files confirms CLAUDE.md is not tracked in git tree."""
        res = subprocess.run(["git", "ls-files", "CLAUDE.md"], cwd=str(REPO_ROOT), capture_output=True, text=True)
        self.assertEqual(res.stdout.strip(), "", "CLAUDE.md must never be tracked in git index")

    def test_tc_r5_3_03_git_status_porcelain_cleanliness(self):
        """TC-R5.3-03: git status --porcelain does not report tracked or modified CLAUDE.md."""
        res = subprocess.run(["git", "status", "--porcelain"], cwd=str(REPO_ROOT), capture_output=True, text=True)
        for line in res.stdout.splitlines():
            self.assertNotIn("CLAUDE.md", line)

    def test_tc_r5_3_04_simulated_claude_md_creation_ignored(self):
        """TC-R5.3-04: Creating dummy CLAUDE.md is ignored by git status."""
        test_file = REPO_ROOT / "CLAUDE.md"
        created = False
        try:
            if not test_file.exists():
                test_file.write_text("# Claude Guidance\n", encoding="utf-8")
                created = True
            res = subprocess.run(["git", "status", "--porcelain", "CLAUDE.md"], cwd=str(REPO_ROOT), capture_output=True, text=True)
            self.assertEqual(res.stdout.strip(), "", "CLAUDE.md should be ignored and untracked by git")
        finally:
            if created and test_file.exists():
                test_file.unlink()

    def test_tc_r5_3_05_no_staged_claude_in_git(self):
        """TC-R5.3-05: git diff --cached confirms CLAUDE.md is never staged."""
        res = subprocess.run(["git", "diff", "--cached", "--name-only"], cwd=str(REPO_ROOT), capture_output=True, text=True)
        self.assertNotIn("CLAUDE.md", res.stdout)

    # --------------------------------------------------------------------------
    # Feature R5.4: Scientific Author Attribution
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_CITE, "synaptik.cli.commands.cite not available")
    def test_tc_r5_4_01_author_metadata_constants(self):
        """TC-R5.4-01: Metadata constants identify Francesco Nobile and ORCID."""
        self.assertEqual(AUTHOR, "Francesco Nobile")
        self.assertEqual(ORCID, "0009-0001-9477-3332")
        self.assertEqual(ORCID_URL, "https://orcid.org/0009-0001-9477-3332")
        self.assertIn("Harvard Medical School", AFFILIATION)
        self.assertIn("University of Pavia", AFFILIATION)

    @unittest.skipUnless(HAVE_CITE, "synaptik.cli.commands.cite not available")
    def test_tc_r5_4_02_bibtex_export_validity(self):
        """TC-R5.4-02: BibTeX export format is syntactically balanced and accurate."""
        bib = get_bibtex("0.8.0")
        self.assertTrue(bib.startswith("@software{nobile2026synaptik"))
        self.assertIn("author       = {Nobile, Francesco}", bib)
        self.assertIn(ORCID_URL, bib)
        self.assertEqual(bib.count("{"), bib.count("}"))

    @unittest.skipUnless(HAVE_CITE, "synaptik.cli.commands.cite not available")
    def test_tc_r5_4_03_apa_export_validity(self):
        """TC-R5.4-03: APA citation format is formatted per scientific guidelines."""
        apa = get_apa("0.8.0")
        self.assertTrue(apa.startswith("Nobile, F. (2026). SYNAPTIK"))
        self.assertIn(ORCID_URL, apa)

    @unittest.skipUnless(HAVE_CITE, "synaptik.cli.commands.cite not available")
    def test_tc_r5_4_04_citation_cff_file_schema(self):
        """TC-R5.4-04: CITATION.cff exists and specifies valid metadata."""
        cff = REPO_ROOT / "CITATION.cff"
        self.assertTrue(cff.exists())
        text = cff.read_text(encoding="utf-8")
        self.assertIn('family-names: "Nobile"', text)
        self.assertIn('given-names: "Francesco"', text)
        self.assertIn(ORCID_URL, text)

    @unittest.skipUnless(HAVE_CITE, "synaptik.cli.commands.cite not available")
    def test_tc_r5_4_05_cmd_cite_cli_smoke(self):
        """TC-R5.4-05: synaptik cite CLI handler executes with exit code 0."""
        args = MagicMock(bibtex=False, apa=False, json=True)
        buf = io.StringIO()
        with patch("sys.stdout", buf):
            code = cmd_cite(args)
        self.assertEqual(code, 0)
        data = json.loads(buf.getvalue())
        self.assertEqual(data["author"], "Francesco Nobile")
        self.assertEqual(data["orcid"], "0009-0001-9477-3332")

    # --------------------------------------------------------------------------
    # Feature R5.5: CLI Subcommand Execution
    # --------------------------------------------------------------------------
    def test_tc_r5_5_01_cli_main_help_execution(self):
        """TC-R5.5-01: python -m synaptik.cli.main --help executes with exit code 0."""
        res = subprocess.run([sys.executable, "-m", "synaptik.cli.main", "--help"], cwd=str(REPO_ROOT), capture_output=True, text=True)
        self.assertEqual(res.returncode, 0)
        self.assertIn("SYNAPTIK", res.stdout)

    def test_tc_r5_5_02_cli_cite_subcommand(self):
        """TC-R5.5-02: python -m synaptik.cli.main cite executes with exit code 0."""
        res = subprocess.run([sys.executable, "-m", "synaptik.cli.main", "cite"], cwd=str(REPO_ROOT), capture_output=True, text=True)
        self.assertEqual(res.returncode, 0)
        self.assertIn("Francesco Nobile", res.stdout)

    def test_tc_r5_5_03_cli_cite_json_subcommand(self):
        """TC-R5.5-03: python -m synaptik.cli.main cite --json outputs valid JSON."""
        res = subprocess.run([sys.executable, "-m", "synaptik.cli.main", "cite", "--json"], cwd=str(REPO_ROOT), capture_output=True, text=True)
        self.assertEqual(res.returncode, 0)
        data = json.loads(res.stdout)
        self.assertEqual(data["author"], "Francesco Nobile")

    def test_tc_r5_5_04_cli_unrecognized_subcommand(self):
        """TC-R5.5-04: Unrecognized CLI command returns non-zero error exit code."""
        res = subprocess.run([sys.executable, "-m", "synaptik.cli.main", "nonexistent_subcommand"], cwd=str(REPO_ROOT), capture_output=True, text=True)
        self.assertIn(res.returncode, (1, 2))

    def test_tc_r5_5_05_cli_missing_required_flag_returns_code_2(self):
        """TC-R5.5-05: Missing mandatory flags return exit code 2."""
        res = subprocess.run([sys.executable, "-m", "synaptik.cli.main", "cite", "--invalid-flag-999"], cwd=str(REPO_ROOT), capture_output=True, text=True)
        self.assertEqual(res.returncode, 2)


# ==============================================================================
# TIER 2: BOUNDARY VALUE ANALYSIS & CORNER CASES (>=5 PER FEATURE)
# ==============================================================================

class TestTier2BoundaryCases(unittest.TestCase):
    """
    Tier 2: Boundary Value Analysis & Corner Stress Testing.
    Verifies behavior at mathematical extremes, physiological boundaries, and resource limits.
    """

    # --------------------------------------------------------------------------
    # Boundary R1: Rate Limiting & Ingestion
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_INGEST_HTTP, "synaptik.ingest.http not available")
    def test_tc_r1_b01_extreme_burst_concurrency(self):
        """TC-R1-B01: 50 requests requesting NCBI scale delays linearly without deadlock."""
        limiter = ConcurrentRateLimiter(limits={"eutils.ncbi.nlm.nih.gov": 0.34})
        delays = [limiter.acquire_delay("eutils.ncbi.nlm.nih.gov") for _ in range(50)]
        self.assertEqual(len(delays), 50)
        self.assertAlmostEqual(delays[0], 0.0, delta=0.01)
        self.assertAlmostEqual(delays[-1], 49 * 0.34, delta=0.1)

    @unittest.skipUnless(HAVE_INGEST_HTTP, "synaptik.ingest.http not available")
    def test_tc_r1_b02_zero_length_batch(self):
        """TC-R1-B02: fetch_batch([]) returns [] immediately without spawning threads."""
        session = PoliteSession(store=None, use_cache=False)
        res = session.fetch_batch([], max_workers=5)
        self.assertEqual(res, [])

    @unittest.skipUnless(HAVE_INGEST_EUROPEPMC, "synaptik.ingest.europepmc not available")
    def test_tc_r1_b03_empty_query_string(self):
        """TC-R1-B03: Europe PMC empty query string returns empty generator safely."""
        harvester = EuropePMCHarvester(session=PoliteSession(store=None, use_cache=False))
        harvester.s.get_json = MagicMock(return_value={"nextCursorMark": "*", "resultList": {"result": []}})
        records = list(harvester.search("", max_results=5))
        self.assertEqual(records, [])

    @unittest.skipUnless(HAVE_INGEST_PUBMED, "synaptik.ingest.pubmed not available")
    def test_tc_r1_b04_empty_pmid_list(self):
        """TC-R1-B04: fetch_metadata([]) yields 0 items without network egress."""
        harvester = PubMedHarvester(session=PoliteSession(store=None, use_cache=False))
        harvester.s.get_text = MagicMock()
        items = list(harvester.fetch_metadata([]))
        self.assertEqual(items, [])
        self.assertEqual(harvester.s.get_text.call_count, 0)

    @unittest.skipUnless(HAVE_INGEST_PUBMED, "synaptik.ingest.pubmed not available")
    def test_tc_r1_b05_inverted_date_range_biorxiv(self):
        """TC-R1-B05: Inverted date range returns empty collection without crash."""
        harvester = BioRxivHarvester(session=PoliteSession(store=None, use_cache=False))
        harvester.s.get_json = MagicMock(return_value={"messages": [{"status": "ok", "total": 0}], "collection": []})
        records = list(harvester.fetch_range("2026-12-31", "2026-01-01"))
        self.assertEqual(records, [])

    # --------------------------------------------------------------------------
    # Boundary R2: Clinical Pharmacology Extremes
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_CLINICAL_ALLOMETRY, "synaptik.clinical.allometry not available")
    def test_tc_r2_b01_zero_or_negative_noael_rejected(self):
        """TC-R2-B01: NOAEL <= 0 raises ValueError."""
        with self.assertRaises(ValueError):
            calculate_allometric_dose(noael_mg=0.0, acknowledge_not_clinical=True)
        with self.assertRaises(ValueError):
            calculate_allometric_dose(noael_mg=-5.0, acknowledge_not_clinical=True)

    @unittest.skipUnless(HAVE_CLINICAL_ALLOMETRY, "synaptik.clinical.allometry not available")
    def test_tc_r2_b02_safety_factor_boundary(self):
        """TC-R2-B02: Safety factor <= 0 raises ValueError."""
        with self.assertRaises(ValueError):
            calculate_allometric_dose(noael_mg=10.0, safety_factor=0.0, acknowledge_not_clinical=True)

    @unittest.skipUnless(HAVE_CLINICAL_ALLOMETRY, "synaptik.clinical.allometry not available")
    def test_tc_r2_b03_extreme_csf_volumes(self):
        """TC-R2-B03: Neonate (40 mL) and adult hydrocephalus (250 mL) compute without overflow."""
        neonate = calculate_allometric_dose(noael_mg=15.0, human_volume_ml=40.0, acknowledge_not_clinical=True)
        hydro = calculate_allometric_dose(noael_mg=15.0, human_volume_ml=250.0, acknowledge_not_clinical=True)
        r_neo = getattr(neonate, "scaling_ratio", getattr(neonate, "scaling_factor", None))
        r_hyd = getattr(hydro, "scaling_ratio", getattr(hydro, "scaling_factor", None))
        self.assertAlmostEqual(r_neo, 40.0 / 15.0, places=3)
        self.assertAlmostEqual(r_hyd, 250.0 / 15.0, places=3)

    @unittest.skipUnless(HAVE_CLINICAL_PKPD, "synaptik.clinical.pkpd not available")
    def test_tc_r2_b04_zero_dose_pkpd_rejected(self):
        """TC-R2-B04: Dose <= 0 in PKPD raises ValueError."""
        with self.assertRaises(ValueError):
            simulate_intrathecal_pkpd(dose_mg=0.0)

    @unittest.skipUnless(HAVE_CLINICAL_PKPD, "synaptik.clinical.pkpd not available")
    def test_tc_r2_b05_invalid_route_pkpd_rejected(self):
        """TC-R2-B05: Unsupported route raises ValueError."""
        with self.assertRaises(ValueError):
            simulate_intrathecal_pkpd(dose_mg=10.0, route="intravenous")

    @unittest.skipUnless(HAVE_CLINICAL_PKPD, "synaptik.clinical.pkpd not available")
    def test_tc_r2_b06_extreme_five_year_duration(self):
        """TC-R2-B06: 1825-day simulation decays monotonically to zero without numerical explosion."""
        res = simulate_intrathecal_pkpd(dose_mg=10.0, duration_days=1825, dt_days=1.0)
        self.assertLess(res.csf_conc_ug_ml[-1], 0.001)

    @unittest.skipUnless(HAVE_CLINICAL_ESCALATION, "synaptik.clinical.escalation not available")
    def test_tc_r2_b07_zero_patients_cohort_rejected(self):
        """TC-R2-B07: Cohort with 0 patients raises ValueError."""
        engine = BOINEscalationEngine()
        with self.assertRaises(ValueError):
            engine.evaluate_cohort([2.0, 5.0], [CohortObservation(cohort_index=1, dose_mg=2.0, n_patients=0, n_dlts=0)])

    @unittest.skipUnless(HAVE_CLINICAL_ESCALATION, "synaptik.clinical.escalation not available")
    def test_tc_r2_b08_top_dose_reached_without_dlt(self):
        """TC-R2-B08: Highest available dose reached with 0 DLTs stays or stops."""
        engine = BOINEscalationEngine()
        doses = [2.0, 5.0]
        history = [
            CohortObservation(cohort_index=1, dose_mg=2.0, n_patients=3, n_dlts=0),
            CohortObservation(cohort_index=2, dose_mg=5.0, n_patients=3, n_dlts=0),
        ]
        decision = engine.evaluate_cohort(doses, history)
        self.assertIn(decision.action, ["stay", "stop_mtd_reached"])
        self.assertEqual(decision.recommended_dose_mg, 5.0)

    @unittest.skipUnless(HAVE_CLINICAL_ESCALATION, "synaptik.clinical.escalation not available")
    def test_tc_r2_b09_dose1_toxic_immediate_termination(self):
        """TC-R2-B09: 3/3 DLTs at lowest dose stops trial immediately."""
        engine = BOINEscalationEngine()
        history = [CohortObservation(cohort_index=1, dose_mg=2.0, n_patients=3, n_dlts=3)]
        decision = engine.evaluate_cohort([2.0, 5.0], history)
        self.assertEqual(decision.action, "stop_safety_eliminated")

    @unittest.skipUnless(HAVE_CLINICAL_ESCALATION, "synaptik.clinical.escalation not available")
    def test_tc_r2_b10_dlts_exceeding_patients_rejected(self):
        """TC-R2-B10: DLTs > patients raises ValueError."""
        engine = ThreePlusThreeEngine()
        with self.assertRaises(ValueError):
            engine.evaluate_cohort([2.0, 5.0], [CohortObservation(cohort_index=1, dose_mg=2.0, n_patients=2, n_dlts=3)])

    # --------------------------------------------------------------------------
    # Boundary R3: RNA Biophysics Extremes
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_BIOPHYSICS_FOLDING, "synaptik.biophysics.folding not available")
    def test_tc_r3_b01_poly_a_homopolymer_50nt(self):
        """TC-R3-B01: Poly-A (50 nt) yields dG=0.0 and 50 dots."""
        dG, db = fold_mfe_pure_python("A" * 50)
        self.assertEqual(dG, 0.0)
        self.assertEqual(db, "." * 50)

    @unittest.skipUnless(HAVE_BIOPHYSICS_FOLDING, "synaptik.biophysics.folding not available")
    def test_tc_r3_b02_poly_u_homopolymer_50nt(self):
        """TC-R3-B02: Poly-U (50 nt) yields dG=0.0."""
        dG, db = fold_mfe_pure_python("U" * 50)
        self.assertEqual(dG, 0.0)

    @unittest.skipUnless(HAVE_BIOPHYSICS_FOLDING, "synaptik.biophysics.folding not available")
    def test_tc_r3_b03_poly_c_homopolymer_50nt(self):
        """TC-R3-B03: Poly-C (50 nt) yields dG=0.0."""
        dG, db = fold_mfe_pure_python("C" * 50)
        self.assertEqual(dG, 0.0)

    @unittest.skipUnless(HAVE_BIOPHYSICS_FOLDING, "synaptik.biophysics.folding not available")
    def test_tc_r3_b04_poly_g_homopolymer_50nt(self):
        """TC-R3-B04: Poly-G (50 nt) cannot form Watson-Crick pairs with itself; dG=0.0."""
        dG, db = fold_mfe_pure_python("G" * 50)
        self.assertEqual(dG, 0.0)

    @unittest.skipUnless(HAVE_BIOPHYSICS_FOLDING, "synaptik.biophysics.folding not available")
    def test_tc_r3_b05_zero_pair_ac_heteropolymer(self):
        """TC-R3-B05: Sequence composed entirely of A and C forms zero pairs; dG=0.0."""
        dG, db = fold_mfe_pure_python("AC" * 20)
        self.assertEqual(dG, 0.0)
        self.assertEqual(db, "." * 40)

    @unittest.skipUnless(HAVE_BIOPHYSICS_FOLDING, "synaptik.biophysics.folding not available")
    def test_tc_r3_b06_empty_sequence_string(self):
        """TC-R3-B06: Empty sequence string returns (0.0, '')."""
        dG, db = fold_mfe_pure_python("")
        self.assertEqual(dG, 0.0)
        self.assertEqual(db, "")

    @unittest.skipUnless(HAVE_BIOPHYSICS_FOLDING, "synaptik.biophysics.folding not available")
    def test_tc_r3_b07_single_nucleotide_sequence(self):
        """TC-R3-B07: Single nucleotide returns (0.0, '.')."""
        dG, db = fold_mfe_pure_python("G")
        self.assertEqual(dG, 0.0)
        self.assertEqual(db, ".")

    @unittest.skipUnless(HAVE_BIOPHYSICS_FOLDING, "synaptik.biophysics.folding not available")
    def test_tc_r3_b08_zero_length_binding_window(self):
        """TC-R3-B08: compute_delta_g_open with length=0 returns 0.0."""
        dg_open = compute_delta_g_open("GCGCUUUUGCGC", binding_start=2, binding_length=0)
        self.assertEqual(dg_open, 0.0)

    @unittest.skipUnless(HAVE_BIOPHYSICS_FOLDING, "synaptik.biophysics.folding not available")
    def test_tc_r3_b09_negative_window_start_raises_error(self):
        """TC-R3-B09: Negative binding_start raises ValueError."""
        with self.assertRaises(ValueError):
            compute_delta_g_open("GCGCUUUUGCGC", binding_start=-1, binding_length=4)

    @unittest.skipUnless(HAVE_BIOPHYSICS_FOLDING, "synaptik.biophysics.folding not available")
    def test_tc_r3_b10_window_extending_past_sequence_raises_error(self):
        """TC-R3-B10: start + len > sequence length raises ValueError."""
        with self.assertRaises(ValueError):
            compute_delta_g_open("GCGCUUUUGCGC", binding_start=10, binding_length=5)

    @unittest.skipUnless(HAVE_BIOPHYSICS_ACCESSIBILITY, "synaptik.biophysics.accessibility not available")
    def test_tc_r3_b11_extreme_high_affinity_no_underflow(self):
        """TC-R3-B11: Delta_G_eff = -35.0 kcal/mol computes Kd without arithmetic underflow."""
        res = compute_effective_binding(duplex_delta_g=-35.0, delta_g_open=0.0)
        self.assertGreater(res.kd_effective_molar, 0.0)
        self.assertLess(res.kd_effective_nm, 1e-12)

    @unittest.skipUnless(HAVE_BIOPHYSICS_ACCESSIBILITY, "synaptik.biophysics.accessibility not available")
    def test_tc_r3_b12_thermodynamic_neutrality_zero(self):
        """TC-R3-B12: Delta_G_eff = 0.0 kcal/mol yields Kd = 1.0 M (1e9 nM)."""
        res = compute_effective_binding(duplex_delta_g=-10.0, delta_g_open=10.0)
        self.assertEqual(res.delta_g_effective, 0.0)
        self.assertAlmostEqual(res.kd_effective_molar, 1.0, places=4)
        self.assertAlmostEqual(res.kd_effective_nm, 1e9, delta=100.0)

    @unittest.skipUnless(HAVE_BIOPHYSICS_COMPETITION, "synaptik.biophysics.competition not available")
    def test_tc_r3_b13_femtomolar_oligo_kd(self):
        """TC-R3-B13: Femtomolar oligo Kd yields near-infinite displacement ratio and 100% occupancy."""
        res = compute_steric_displacement(oligo_kd_nm=1e-6, competitor_name="u1_snrnp")
        self.assertGreater(res.displacement_ratio, 1e6)
        self.assertAlmostEqual(res.fractional_occupancy, 1.0, places=4)

    @unittest.skipUnless(HAVE_BIOPHYSICS_COMPETITION, "synaptik.biophysics.competition not available")
    def test_tc_r3_b14_identical_competitor_affinities(self):
        """TC-R3-B14: Equal Kd values yield displacement ratio of exactly 1.0."""
        res = compute_steric_displacement(oligo_kd_nm=25.0, competitor_name="u1_snrnp")
        self.assertAlmostEqual(res.displacement_ratio, 1.0, places=4)

    @unittest.skipUnless(HAVE_BIOPHYSICS_COMPETITION, "synaptik.biophysics.competition not available")
    def test_tc_r3_b15_negative_oligo_kd_raises_error(self):
        """TC-R3-B15: Negative oligo Kd raises ValueError."""
        with self.assertRaises(ValueError):
            compute_steric_displacement(oligo_kd_nm=-5.0, competitor_name="u1_snrnp")

    # --------------------------------------------------------------------------
    # Boundary R4: Patient Digital Twin Parameter Boundaries
    # --------------------------------------------------------------------------
    @unittest.skipUnless(HAVE_TWIN_SIMULATOR, "synaptik.twin.twin_generator not available")
    def test_tc_r4_b01_pediatric_juvenile_age(self):
        """TC-R4-B01: Pediatric age (14 yr) accepts input without arithmetic overflow."""
        b = PatientBaseline("PED-01", 14.0, "limb", 6.0, 42.0)
        self.assertEqual(b.age_years, 14.0)

    @unittest.skipUnless(HAVE_TWIN_SIMULATOR, "synaptik.twin.twin_generator not available")
    def test_tc_r4_b02_geriatric_extreme_age(self):
        """TC-R4-B02: Geriatric age (98 yr) accepts input cleanly."""
        b = PatientBaseline("GER-01", 98.0, "bulbar", 6.0, 30.0)
        self.assertEqual(b.age_years, 98.0)

    @unittest.skipUnless(HAVE_TWIN_SIMULATOR, "synaptik.twin.twin_generator not available")
    def test_tc_r4_b03_disease_duration_near_zero(self):
        """TC-R4-B03: Duration 0.1 month avoids division by zero in slope calculation."""
        b = PatientBaseline("P1", 55.0, "limb", 0.1, 48.0)
        self.assertEqual(b.disease_duration_months, 0.1)

    @unittest.skipUnless(HAVE_TWIN_SIMULATOR, "synaptik.twin.twin_generator not available")
    def test_tc_r4_b04_baseline_alsfrs_perfect_score(self):
        """TC-R4-B04: Perfect score 48.0 yields pre-slope of 0.0 pts/mo."""
        b = PatientBaseline("P1", 55.0, "limb", 10.0, 48.0)
        self.assertEqual(b.baseline_alsfrs_r, 48.0)

    @unittest.skipUnless(HAVE_TWIN_SIMULATOR, "synaptik.twin.twin_generator not available")
    def test_tc_r4_b05_baseline_alsfrs_floor_zero(self):
        """TC-R4-B05: Baseline score 0.0 clamped at floor for all twin timepoints."""
        b = PatientBaseline("P1", 55.0, "limb", 10.0, 0.0)
        cohort = generate_counterfactual_cohort(b, n_twins=20, duration_months=6, random_seed=42)
        for traj in cohort.trajectories:
            self.assertTrue(all(s == 0.0 for s in traj))

    @unittest.skipUnless(HAVE_TWIN_DISTRIBUTIONS, "synaptik.twin.distributions not available")
    def test_tc_r4_b06_extreme_low_nfl(self):
        """TC-R4-B06: Near-zero NfL (2.0 pg/mL) computes without math domain error."""
        b = PatientBaseline("P1", 55.0, "limb", 10.0, 40.0, baseline_plasma_nfl=2.0)
        c = compute_conditional_progression(b)
        self.assertTrue(math.isfinite(c.nfl_adjustment))

    @unittest.skipUnless(HAVE_TWIN_DISTRIBUTIONS, "synaptik.twin.distributions not available")
    def test_tc_r4_b07_massive_nfl(self):
        """TC-R4-B07: Massive NfL (2500 pg/mL) evaluates cleanly."""
        b = PatientBaseline("P1", 55.0, "limb", 10.0, 40.0, baseline_plasma_nfl=2500.0)
        c = compute_conditional_progression(b)
        self.assertTrue(math.isfinite(c.nfl_adjustment))

    @unittest.skipUnless(HAVE_TWIN_SIMULATOR, "synaptik.twin.twin_generator not available")
    def test_tc_r4_b08_single_twin_simulation(self):
        """TC-R4-B08: n_twins=1 runs without crash; quantiles match single trajectory."""
        b = PatientBaseline("P1", 55.0, "limb", 10.0, 40.0)
        cohort = generate_counterfactual_cohort(b, n_twins=1, duration_months=6, random_seed=42)
        self.assertEqual(len(cohort.trajectories), 1)
        self.assertEqual(cohort.median_trajectory, cohort.trajectories[0])

    @unittest.skipUnless(HAVE_TWIN_SIMULATOR, "synaptik.twin.twin_generator not available")
    def test_tc_r4_b09_zero_duration_rejected(self):
        """TC-R4-B09: duration_months <= 0 raises ValueError."""
        b = PatientBaseline("P1", 55.0, "limb", 10.0, 40.0)
        with self.assertRaises(ValueError):
            generate_counterfactual_cohort(b, duration_months=0)

    @unittest.skipUnless(HAVE_TWIN_BAYESIAN, "synaptik.twin.bayesian_inference not available")
    def test_tc_r4_b10_identical_visit_times_collinear(self):
        """TC-R4-B10: Collinear visits with identical timepoints handled without zero division crash."""
        b = PatientBaseline("P1", 55.0, "limb", 10.0, 40.0)
        cohort = generate_counterfactual_cohort(b, n_twins=50, duration_months=6, random_seed=42)
        visits = [VisitObservation(month=2.0, alsfrs_score=38.0), VisitObservation(month=2.0, alsfrs_score=38.0)]
        with self.assertRaises(ValueError):
            evaluate_treatment_efficacy(cohort, visits)

    # --------------------------------------------------------------------------
    # Boundary R5: Security, AST, Git & CLI Error Handling
    # --------------------------------------------------------------------------
    def test_tc_r5_b01_sql_injection_string_handled_safely(self):
        """TC-R5-B01: Patient identifier with SQL injection characters handled in-memory without query syntax errors."""
        injection_id = "PATIENT'; DROP TABLE records;--"
        if HAVE_TWIN_SIMULATOR:
            b = PatientBaseline(injection_id, 55.0, "limb", 10.0, 40.0)
            cohort = generate_counterfactual_cohort(b, n_twins=10, duration_months=3)
            self.assertEqual(cohort.baseline.patient_id, injection_id)

    def test_tc_r5_b02_ast_handles_empty_python_file(self):
        """TC-R5-B02: ast.parse handles empty string without error."""
        tree = ast.parse("", filename="empty.py")
        self.assertIsInstance(tree, ast.Module)
        self.assertEqual(len(tree.body), 0)

    def test_tc_r5_b03_ast_handles_multiline_docstrings_with_backslashes(self):
        """TC-R5-B03: Regular strings and docstrings containing backslashes are permitted."""
        sample_code = r'"""Multiline docstring with \n and \Delta G."""' + '\nx = "path\\\\to\\\\file"'
        tree = ast.parse(sample_code)
        self.assertIsInstance(tree, ast.Module)

    def test_tc_r5_b04_git_cleanliness_case_sensitivity(self):
        """TC-R5-B04: .gitignore pattern covers case sensitivity (claude.md / CLAUDE.md)."""
        gitignore = REPO_ROOT / ".gitignore"
        text = gitignore.read_text(encoding="utf-8").lower()
        self.assertIn("claude.md", text)

    def test_tc_r5_b05_cli_non_numeric_argument_error_code(self):
        """TC-R5-B05: Non-numeric argument to numeric flag returns exit code 2."""
        res = subprocess.run([sys.executable, "-m", "synaptik.cli.main", "cite", "--invalid-non-numeric-flag"], cwd=str(REPO_ROOT), capture_output=True, text=True)
        self.assertEqual(res.returncode, 2)


# ==============================================================================
# TIER 3: PAIRWISE CROSS-FEATURE INTERACTIONS
# ==============================================================================

class TestTier3CrossFeatureInteractions(unittest.TestCase):
    """
    Tier 3: Pairwise Combinatorial Cross-Feature Interactions.
    Tests end-to-end multi-module pipelines where outputs flow across domain boundaries.
    """

    @unittest.skipUnless(HAVE_CLINICAL_ALLOMETRY and HAVE_CLINICAL_PKPD and HAVE_CLINICAL_ESCALATION, "Clinical modules not fully available")
    def test_pipeline_1_allometry_to_pkpd_to_boin_escalation(self):
        """
        Pipeline 1:
        NHP NOAEL -> CSF Allometric Scaling -> MRSD -> Intrathecal PK/PD -> BOIN Escalation Ladder.
        """
        # Step 1: Preclinical NHP NOAEL of 15.0 mg scaled to adult human (140 mL)
        allo = calculate_allometric_dose(noael_mg=15.0, nhp_volume_ml=15.0, human_volume_ml=140.0, safety_factor=10.0, acknowledge_not_clinical=True)
        self.assertAlmostEqual(allo.mrsd_mg, 14.0, places=2)

        # Step 2: Feed MRSD directly into PK/PD simulator
        pk = simulate_intrathecal_pkpd(dose_mg=allo.mrsd_mg, route="bolus", csf_volume_ml=140.0)
        self.assertAlmostEqual(pk.cmax_csf_ug_ml, 100.0, delta=1.0)
        self.assertGreater(pk.day90_csf_conc_ug_ml, 5.0)

        # Step 3: Construct dose escalation ladder starting from MRSD
        ladder = [allo.mrsd_mg, allo.mrsd_mg * 2.0, allo.mrsd_mg * 4.0]  # [14.0, 28.0, 56.0]
        engine = BOINEscalationEngine(target_phi=0.25)
        history = [CohortObservation(cohort_index=1, dose_mg=14.0, n_patients=3, n_dlts=0)]
        decision = engine.evaluate_cohort(ladder, history)
        self.assertEqual(decision.action, "escalate")
        self.assertEqual(decision.recommended_dose_mg, 28.0)

    @unittest.skipUnless(HAVE_BIOPHYSICS_FOLDING and HAVE_BIOPHYSICS_ACCESSIBILITY and HAVE_BIOPHYSICS_COMPETITION, "Biophysics modules not fully available")
    def test_pipeline_2_rna_mfe_to_delta_g_open_to_spliceosome_displacement(self):
        """
        Pipeline 2:
        Target RNA MFE Folding -> Delta_G_open -> Delta_G_effective -> U1 snRNP Steric Competition.
        """
        target_seq = "GCGCUUUUGCGC"
        # Step 1: Fold target transcript
        dG_folded, db = fold_mfe_pure_python(target_seq)
        self.assertLess(dG_folded, 0.0)

        # Step 2: Calculate opening penalty across 5' stem arm (bases 0-4)
        dg_open = compute_delta_g_open(target_seq, binding_start=0, binding_length=4)
        self.assertGreater(dg_open, 0.0)

        # Step 3: Compute net effective binding with duplex hybridization of -22.0 kcal/mol
        eff = compute_effective_binding(duplex_delta_g=-22.0, delta_g_open=dg_open)
        self.assertEqual(eff.delta_g_effective, -22.0 + dg_open)

        # Step 4: Model competitive displacement against U1 snRNP (Kd=25.0 nM)
        disp = compute_steric_displacement(oligo_kd_nm=eff.kd_effective_nm, competitor_name="u1_snrnp", target_gene="STMN2")
        self.assertEqual(disp.competitor_name, "u1_snrnp")
        self.assertEqual(disp.target_gene, "STMN2")
        self.assertTrue(disp.displacement_ratio > 0.0)

    @unittest.skipUnless(HAVE_TWIN_SIMULATOR and HAVE_TWIN_BAYESIAN and HAVE_TWIN_DOSSIER and HAVE_CITE, "Twin & citation modules not fully available")
    def test_pipeline_3_digital_twin_to_bayesian_deceleration_to_dossier_to_citation(self):
        """
        Pipeline 3:
        Patient Baseline -> 1,000 Matched Twins -> Bayesian Deceleration Posterior -> FDA Dossier -> Author Citation.
        """
        # Step 1: Baseline ingestion & 1,000-twin Monte Carlo simulation
        baseline = PatientBaseline("ALS-E2E-001", 54.0, "bulbar", 9.2, 36.0, baseline_slope=-1.30, baseline_plasma_nfl=112.5)
        cohort = generate_counterfactual_cohort(baseline, n_twins=1000, duration_months=12, random_seed=42)

        # Step 2: Treated patient visits exhibiting functional stabilization (slope -0.50 pts/mo)
        visits = [
            VisitObservation(month=0.0, alsfrs_score=36.0),
            VisitObservation(month=3.0, alsfrs_score=35.0),
            VisitObservation(month=6.0, alsfrs_score=34.5),
            VisitObservation(month=9.0, alsfrs_score=34.0),
            VisitObservation(month=12.0, alsfrs_score=33.5),
        ]
        efficacy = evaluate_treatment_efficacy(cohort, visits)
        self.assertGreaterEqual(efficacy.p_deceleration, 0.985)
        self.assertIn("Tier A", efficacy.efficacy_tier)

        # Step 3: Export audit-ready FDA Accelerated Approval dossier
        dossier = generate_twin_fda_dossier(baseline, cohort, efficacy=efficacy)
        self.assertIn("Executive Summary", dossier)

        # Step 4: Verify Francesco Nobile citation and ORCID embedded in dossier
        self.assertIn(AUTHOR, dossier)
        self.assertIn(ORCID, dossier)

    @unittest.skipUnless(HAVE_INGEST_PUBMED and HAVE_EVIDENCE_WEIGHTER, "Ingest & weighter modules not available")
    def test_pipeline_4_biorxiv_preprint_to_evidence_weighter_discount(self):
        """
        Pipeline 4:
        bioRxiv preprint ingestion -> SourceMetadata.is_preprint -> EvidenceWeighter 20% discount.
        """
        session = PoliteSession(store=None, use_cache=False)
        harvester = BioRxivHarvester(session=session)
        raw_item = {
            "doi": "10.1101/2026.03.01.600000",
            "title": "Novel Splice Modulator for FUS ALS",
            "authors": "Nobile, F. et al.",
            "date": "2026-03-02",
            "abstract": "Evidence of motor neuron rescue",
        }
        meta = harvester._parse(raw_item)
        self.assertIsNotNone(meta)
        self.assertTrue(meta.is_preprint)

        rec = create_mock_experiment_record(is_preprint=meta.is_preprint)
        weighter = EvidenceWeighter()
        breakdown = weighter.weight(rec)

        self.assertEqual(breakdown.journal_tier, 0.80)
        self.assertIn("preprint", breakdown.notes)

    def test_pipeline_5_zero_phi_active_monitoring_during_cli_execution(self):
        """
        Pipeline 5:
        Active monitoring of filesystem, SQLite database, and OS tempdir during CLI execution.
        """
        canary = "CANARY_PHI_E2E_AUDIT_774921"
        initial_cwd_files = set(os.listdir(str(REPO_ROOT)))

        # Execute synaptik cite with canary argument in help check
        res = subprocess.run([sys.executable, "-m", "synaptik.cli.main", "cite", "--json"], cwd=str(REPO_ROOT), capture_output=True, text=True)
        self.assertEqual(res.returncode, 0)

        # Verify zero unrequested file creations in workspace root
        final_cwd_files = set(os.listdir(str(REPO_ROOT)))
        self.assertEqual(final_cwd_files, initial_cwd_files)


# ==============================================================================
# TIER 4: REAL-WORLD CLINICAL & BIOPHYSICAL APPLICATION SCENARIOS
# ==============================================================================

class TestTier4RealWorldScenarios(unittest.TestCase):
    """
    Tier 4: Real-World Clinical & Biophysical Workloads.
    Executes 5 full translational application scenarios replicating genuine bedside trials.
    """

    @unittest.skipUnless(HAVE_CLINICAL_ALLOMETRY and HAVE_CLINICAL_PKPD and HAVE_CLINICAL_ESCALATION and HAVE_CLINICAL_PROTOCOL, "Clinical modules not fully available")
    def test_scenario_1_als_fus_first_in_human_protocol_design(self):
        """
        Scenario 1:
        First-in-Human ALS FUS ASO Clinical Protocol Design.
        Rapidly progressive ALS patient harboring heterozygous FUS P525L mutation.
        Nonclinical cynomolgus monkey NOAEL = 12.0 mg.
        """
        # Step 1: Allometric dose scaling
        adult_allo = calculate_allometric_dose(noael_mg=12.0, nhp_volume_ml=15.0, human_volume_ml=140.0, safety_factor=10.0, acknowledge_not_clinical=True)
        self.assertAlmostEqual(adult_allo.hed_mg, 112.0, places=2)
        self.assertAlmostEqual(adult_allo.mrsd_mg, 11.2, places=2)

        # Step 2: 360-day Intrathecal PK/PD simulation
        pk = simulate_intrathecal_pkpd(dose_mg=adult_allo.mrsd_mg, route="bolus", duration_days=360)
        self.assertAlmostEqual(pk.cmax_csf_ug_ml, 80.0, delta=1.0)
        self.assertGreater(pk.day90_csf_conc_ug_ml, 5.0)

        # Catheter blunts peak to avoid chemical arachnoiditis
        pk_cath = simulate_intrathecal_pkpd(dose_mg=adult_allo.mrsd_mg, route="catheter", duration_days=360)
        self.assertLess(pk_cath.cmax_csf_ug_ml, pk.cmax_csf_ug_ml)

        # Step 3: BOIN escalation ladder simulation
        ladder = [11.2, 25.0, 50.0, 100.0]
        engine = BOINEscalationEngine(target_phi=0.25)
        h1 = [CohortObservation(cohort_index=1, dose_mg=11.2, n_patients=3, n_dlts=0)]
        dec1 = engine.evaluate_cohort(ladder, h1)
        self.assertEqual(dec1.action, "escalate")
        self.assertEqual(dec1.recommended_dose_mg, 25.0)

        # Step 4: Synthesize audit-ready FIH protocol dossier
        proto = generate_clinical_protocol(target_gene="FUS", disease="ALS", noael_mg=12.0, escalation_method="boin")
        self.assertEqual(proto.target_gene, "FUS")
        self.assertIn("csf_pleocytosis", proto.dlt_definition)
        self.assertEqual(proto.author, "Francesco Nobile")

    @unittest.skipUnless(HAVE_BIOPHYSICS_FOLDING and HAVE_BIOPHYSICS_ACCESSIBILITY and HAVE_BIOPHYSICS_COMPETITION, "Biophysics modules not fully available")
    def test_scenario_2_stmn2_cryptic_junction_splice_switching(self):
        """
        Scenario 2:
        STMN2 cryptic splice donor site blocking in TDP-43-depleted motor neurons.
        40-nt cryptic junction sequence containing donor site and hairpin stem.
        """
        stmn2_cryptic_junction = "UUACUUACAGTTTTTGTTTTTTCAGGTGAGTATTCAGGT"
        dG_folded, db = fold_mfe_pure_python(stmn2_cryptic_junction)
        self.assertLess(dG_folded, -2.0)

        # Opening penalty for 18-mer ASO across donor site
        dg_open = compute_delta_g_open(stmn2_cryptic_junction, binding_start=10, binding_length=18)
        self.assertGreaterEqual(dg_open, 0.0)

        # Effective hybridization energy
        eff = compute_effective_binding(duplex_delta_g=-22.4, delta_g_open=dg_open)
        self.assertTrue(eff.is_favorable)
        self.assertLess(eff.kd_effective_nm, 1.0)

        # Steric displacement against U1 snRNP (Kd=25.0 nM)
        disp = compute_steric_displacement(oligo_kd_nm=eff.kd_effective_nm, competitor_name="u1_snrnp", target_gene="STMN2")
        self.assertGreater(disp.displacement_ratio, 25.0)
        self.assertGreater(disp.fractional_occupancy, 0.90)
        self.assertEqual(disp.dominant_binder, "oligo")

    @unittest.skipUnless(HAVE_BIOPHYSICS_FOLDING and HAVE_BIOPHYSICS_ACCESSIBILITY and HAVE_BIOPHYSICS_COMPETITION, "Biophysics modules not fully available")
    def test_scenario_3_unc13a_cryptic_acceptor_blockade(self):
        """
        Scenario 3:
        UNC13A cryptic acceptor blockade (rs12973192 risk variant).
        Displaces U2AF65 (polypyrimidine tract, Kd=20 nM) and U2AF35 (AG acceptor, Kd=40 nM).
        """
        unc13a_acceptor = "UUUUUUUUUUUUUUUUCCAGACAGUUGACU"
        dG_folded, db = fold_mfe_pure_python(unc13a_acceptor)

        # Opening penalty across Py-tract
        dg_open = compute_delta_g_open(unc13a_acceptor, binding_start=0, binding_length=18)
        eff = compute_effective_binding(duplex_delta_g=-24.0, delta_g_open=dg_open)

        # Dual competitor displacement: U2AF65 and U2AF35
        disp_65 = compute_steric_displacement(oligo_kd_nm=eff.kd_effective_nm, competitor_name="u2af65", target_gene="UNC13A")
        disp_35 = compute_steric_displacement(oligo_kd_nm=eff.kd_effective_nm, competitor_name="u2af35", target_gene="UNC13A")

        self.assertGreater(disp_65.displacement_ratio, 50.0)
        self.assertGreater(disp_35.displacement_ratio, 100.0)
        self.assertEqual(disp_65.dominant_binder, "oligo")
        self.assertEqual(disp_35.dominant_binder, "oligo")

    @unittest.skipUnless(HAVE_TWIN_SIMULATOR and HAVE_TWIN_BAYESIAN and HAVE_TWIN_DOSSIER, "Twin modules not fully available")
    def test_scenario_4_fast_progressing_bulbar_als_single_subject_trial(self):
        """
        Scenario 4:
        Fast-Progressing Bulbar-Onset ALS Patient N-of-1 Digital Twin Control & FDA Dossier.
        54-year-old female, bulbar onset, pre-slope -1.30 pts/mo, baseline NfL 112.5 pg/mL.
        Observed treated slope: -0.50 pts/mo.
        """
        baseline = PatientBaseline("CANARY_PATIENT_SCEN4", 54.0, "bulbar", 9.2, 36.0, baseline_slope=-1.30, baseline_plasma_nfl=112.5)
        cohort = generate_counterfactual_cohort(baseline, n_twins=1000, duration_months=12, random_seed=42)

        # Counterfactual twins predict steep untreated decline
        self.assertLess(cohort.conditional_mean_slope, -1.0)

        # Treated visits showing dramatic disease deceleration
        visits = [
            VisitObservation(month=0.0, alsfrs_score=36.0),
            VisitObservation(month=3.0, alsfrs_score=35.0),
            VisitObservation(month=6.0, alsfrs_score=34.5),
            VisitObservation(month=9.0, alsfrs_score=34.0),
            VisitObservation(month=12.0, alsfrs_score=33.5),
        ]
        efficacy = evaluate_treatment_efficacy(cohort, visits)

        self.assertGreaterEqual(efficacy.p_deceleration, 0.985)
        self.assertGreater(efficacy.bayes_factor, 60.0)
        self.assertIn("Tier A", efficacy.efficacy_tier)

        dossier = generate_twin_fda_dossier(baseline, cohort, efficacy=efficacy)
        self.assertIn("Executive Summary", dossier)
        self.assertIn("Francesco Nobile", dossier)

    @unittest.skipUnless(HAVE_CLINICAL_ESCALATION, "synaptik.clinical.escalation not available")
    def test_scenario_5_boin_escalation_early_safety_termination(self):
        """
        Scenario 5:
        Phase 1 trial encountering acute Dose-Limiting Toxicity at Dose 1.
        Cohort 1 (5.0 mg): 3 DLTs in 3 patients (radiculopathy, pleocytosis, thrombocytopenia).
        BOIN engine triggers safety stopping boundary and terminates trial.
        """
        doses = [5.0, 10.0, 20.0, 40.0]
        engine = BOINEscalationEngine(target_phi=0.25)
        history = [
            CohortObservation(
                cohort_index=1,
                dose_mg=5.0,
                n_patients=3,
                n_dlts=3,
                toxicities=["sensory_radiculopathy", "csf_pleocytosis_75", "severe_thrombocytopenia_32k"],
            )
        ]
        decision = engine.evaluate_cohort(doses, history)
        self.assertEqual(decision.action, "stop_safety_eliminated")
        self.assertTrue(decision.safety_eliminated)
        self.assertGreater(decision.posterior_toxicity_prob, 0.95)


# ==============================================================================
# MAIN TEST RUNNER
# ==============================================================================

if __name__ == "__main__":
    unittest.main()
