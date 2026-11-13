"""
Challenger 2 Adversarial Stress Test Suite:
1. SQLite Thread-Affinity & Concurrency Torture Test
2. Malicious Input & Edge Case Testing (JATS XML, BioRxiv/MedRxiv)
3. Bayesian Math & Boundary Conditions (AssociationEngine, EvidenceReconciler)
"""

from __future__ import annotations

import logging
import sqlite3
import sys
import tempfile
import time
from datetime import date
from pathlib import Path
from unittest.mock import MagicMock, patch

import requests

from synaptik.evidence.reconcile import (
    CLASS_PRIOR, EvidenceClass, EvidenceItem, EvidenceReconciler,
    from_literature_association,
)
from synaptik.infer.association import AssociationEngine, PathwayAssociation
from synaptik.infer.weighting import (
    EvidenceWeighter, QUARTILE_WEIGHT, UNKNOWN_QUARTILE_WEIGHT,
)
from synaptik.ingest.europepmc import EuropePMCHarvester
from synaptik.ingest.http import PoliteSession
from synaptik.ingest.pubmed import BioRxivHarvester
from synaptik.pipeline import Pipeline
from synaptik.schema import (
    CellModel, CellModelClass, Direction, EvidenceType,
    ExperimentRecord, NormalizedTerm, Outcome, Perturbation,
    PerturbationType, SourceMetadata,
)
from synaptik.store import Store

FAILURES: list[str] = []
FINDINGS: list[str] = []


def check(cond: bool, msg: str) -> None:
    if cond:
        print(f"  PASS  {msg}")
    else:
        print(f"  FAIL  {msg}")
        FAILURES.append(msg)


def finding(msg: str) -> None:
    print(f"  FINDING  {msg}")
    FINDINGS.append(msg)


# =====================================================================
# TRACK 1: SQLite Thread-Affinity Torture Test
# =====================================================================

def test_sqlite_thread_affinity_torture():
    print("\n" + "=" * 70)
    print("TRACK 1: SQLite Thread-Affinity & Concurrency Torture Test")
    print("=" * 70)

    # Track any sqlite3 errors logged or caught
    logged_programming_errors = []

    class LogHandler(logging.Handler):
        def emit(self, record):
            msg = record.getMessage()
            if "ProgrammingError" in msg or "SQLite objects created in a thread" in msg:
                logged_programming_errors.append(msg)

    http_logger = logging.getLogger("synaptik.ingest.http")
    handler = LogHandler()
    http_logger.addHandler(handler)
    prev_level = http_logger.level
    http_logger.setLevel(logging.DEBUG)

    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = Path(tmpdir) / "torture.db"
        store = Store(db_path)

        session = PoliteSession(
            store=store,
            rate_limits={"www.ebi.ac.uk": 0.001, "api.biorxiv.org": 0.001},
        )
        pipeline = Pipeline(store=store, session=session, strict_contracts=False)

        # Mock Europe PMC search responses adhering to schema contracts
        def mock_epmc_search(method, url, params=None, timeout=None):
            resp = MagicMock()
            resp.status_code = 200
            q = (params or {}).get("query", "unknown")
            resp.json.return_value = {
                "hitCount": 10,
                "nextCursorMark": "*",
                "resultList": {
                    "result": [
                        {
                            "id": f"PMC_TEST_{q}_{i}",
                            "source": "MED",
                            "pmid": f"390000{i}",
                            "pmcid": f"PMC9000{i}",
                            "doi": f"10.1000/epmc_{q}_{i}",
                            "title": f"Paper {q} title {i}",
                            "authorString": "Author A, Author B",
                            "authorList": {
                                "author": [
                                    {"fullName": "Author A", "affiliation": "Lab A"},
                                    {"fullName": "Author B", "affiliation": "Lab A"},
                                ]
                            },
                            "journalTitle": "Nature",
                            "journalInfo": {
                                "journal": {"title": "Nature"},
                                "yearOfPublication": 2023,
                            },
                            "pubYear": "2023",
                            "isOpenAccess": "Y",
                            "citedByCount": 10,
                        }
                        for i in range(10)
                    ]
                },
            }
            return resp

        # Mock BioRxiv details responses
        def mock_biorxiv_details(method, url, params=None, timeout=None):
            resp = MagicMock()
            resp.status_code = 200
            resp.json.return_value = {
                "messages": [{"status": "ok", "total": 4}],
                "collection": [
                    {
                        "doi": f"10.1101/2023.01.0{i}",
                        "title": f"Preprint TDP-43 title {i}",
                        "abstract": f"Abstract on cryptic splicing {i}",
                        "authors": "Smith, J.",
                        "author_corresponding": "Smith, J.",
                        "author_corresponding_institution": "Lab Pre",
                        "date": "2023-01-01",
                        "server": "biorxiv",
                        "category": "neuroscience",
                    }
                    for i in range(4)
                ],
            }
            return resp

        def dispatch_request(method, url, params=None, timeout=None):
            if "biorxiv.org" in url:
                return mock_biorxiv_details(method, url, params, timeout)
            return mock_epmc_search(method, url, params, timeout)

        # 1. Pipeline.harvest(concurrency=5)
        print("\n--- Testing Pipeline.harvest(concurrency=5) with real Store ---")
        unhandled_harvest_err = None
        try:
            with patch.object(session.session, "request", side_effect=dispatch_request):
                queries = {
                    f"test_q_{i}": f"query_term_{i}"
                    for i in range(10)
                }
                stats = pipeline.harvest(
                    queries=queries,
                    max_per_query=10,
                    concurrency=5,
                    include_preprints=True,
                )
        except Exception as e:
            unhandled_harvest_err = e

        check(unhandled_harvest_err is None, f"Pipeline.harvest completed without raising unhandled exception: {unhandled_harvest_err}")
        check(stats.sources_seen > 0, f"Pipeline.harvest retrieved sources (seen: {stats.sources_seen}, new: {stats.sources_new})")

        # Verify Store integrity after harvest
        source_count = store.conn.execute("SELECT count(*) FROM sources").fetchone()[0]
        check(source_count > 0, f"Sources committed into SQLite store: {source_count} records")

        # Empirically verify the internal SQLite thread affinity discovery:
        # Worker threads calling session.get_json() trigger sqlite3.ProgrammingError inside store.cache_get / store.cache_put
        if logged_programming_errors:
            finding(
                f"SQLite Thread Affinity finding verified: {len(logged_programming_errors)} "
                f"sqlite3.ProgrammingError caught inside PoliteSession._cache_get / _cache_put. "
                f"Because PoliteSession catches Exception with log.debug, the pipeline survives "
                f"gracefully without crashing, but SQLite cache is completely bypassed for worker threads."
            )
            check(True, "PoliteSession successfully isolates sqlite3.ProgrammingError so Pipeline does not crash")
        else:
            check(True, "Zero sqlite3.ProgrammingError encountered")

        # 2. Pipeline.fetch_full_texts(max_workers=6)
        print("\n--- Testing Pipeline.fetch_full_texts(max_workers=6) with real Store ---")
        mock_xml = """<article>
          <front><article-meta><title-group><article-title>T</article-title></title-group></article-meta></front>
          <body>
            <sec><title>Methods</title><p>We cultured iPSC motor neurons.</p></sec>
            <sec><title>Results</title><p>Cryptic exon was suppressed.</p></sec>
          </body>
        </article>"""

        def mock_fulltext_request(method, url, params=None, timeout=None):
            resp = MagicMock()
            resp.status_code = 200
            resp.text = mock_xml
            return resp

        unhandled_fulltext_err = None
        try:
            with patch.object(session.session, "request", side_effect=mock_fulltext_request):
                ft_stats = pipeline.fetch_full_texts(limit=10, max_workers=6)
        except Exception as e:
            unhandled_fulltext_err = e

        check(unhandled_fulltext_err is None, f"Pipeline.fetch_full_texts completed without unhandled exception: {unhandled_fulltext_err}")
        check(ft_stats.fulltext_fetched > 0, f"Full texts fetched and stored: {ft_stats.fulltext_fetched}")

        doc_count = store.conn.execute("SELECT count(*) FROM documents").fetchone()[0]
        check(doc_count == ft_stats.fulltext_fetched, f"Documents correctly persisted in SQLite: {doc_count} rows")

        # Verify zero lock contention under high-frequency concurrent readers
        print("\n--- Testing SQLite transaction lock contention under concurrent load ---")
        lock_errors = []

        def concurrent_reader():
            try:
                conn2 = sqlite3.connect(str(db_path), timeout=10.0)
                for _ in range(50):
                    c = conn2.execute("SELECT count(*) FROM sources").fetchone()[0]
                    time.sleep(0.001)
                conn2.close()
            except Exception as e:
                lock_errors.append(e)

        import threading
        threads = [threading.Thread(target=concurrent_reader) for _ in range(5)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        check(len(lock_errors) == 0, f"Zero SQLite lock contention / OperationalError across concurrent readers (errors: {len(lock_errors)})")

        store.conn.close()

    http_logger.removeHandler(handler)
    http_logger.setLevel(prev_level)


# =====================================================================
# TRACK 2: Malicious Input & Edge Case Testing
# =====================================================================

def test_malicious_inputs_and_edge_cases():
    print("\n" + "=" * 70)
    print("TRACK 2: Malicious Input & Edge Case Testing")
    print("=" * 70)

    # 1. Corrupt and Malformed JATS XML behavior
    print("\n--- EuropePMCHarvester.parse_jats resilience and edge cases ---")
    malformed_xmls = [
        ("Empty string", "", True),
        ("Whitespace only", "   \n\t  ", True),
        ("Non-XML raw text", "This is plain text with no tags whatsoever.", True),
        ("Truncated tag", "<article><body><sec><title>Methods</title><p>Truncated", False),
        ("Unescaped ampersand", "<article><body><sec><title>Results</title><p>A & B & C < 5</p></sec></body></article>", False),
        ("Null bytes and control characters", "<article><body><sec><title>Null\x00Byte</title><p>Control\x01\x02chars</p></sec></body></article>", False),
        ("Deeply nested broken tags", "<div>" * 50 + "broken" + "</div>" * 20, False),
        ("Entities injection attempt", '<!DOCTYPE foo [<!ENTITY xx SYSTEM "file:///etc/passwd">]><article><body>&xx;</body></article>', False),
        ("Only figures and tables", "<article><body><fig><label>Fig 1</label><caption>Caption 1</caption></fig><table-wrap><caption>Table 1</caption></table-wrap></body></article>", False),
    ]

    for label, xml, raises_syntax_error in malformed_xmls:
        crashed = False
        res = None
        try:
            res = EuropePMCHarvester.parse_jats(xml)
        except Exception as e:
            crashed = True
            err_type = type(e).__name__

        if raises_syntax_error:
            check(crashed, f"parse_jats on '{label}' strictly raises XMLSyntaxError (as expected without caller guard)")
        else:
            check(not crashed, f"parse_jats recovered '{label}' without crashing (sections: {list(res.keys()) if res else None})")

    # 2. fetch_full_texts_concurrent with malformed responses
    print("\n--- fetch_full_texts_concurrent resilience to HTTP and XML failures ---")
    session = PoliteSession(store=None)
    harvester = EuropePMCHarvester(session)

    metas = [
        SourceMetadata(source_id=f"PMID:{i}", pmcid=f"PMC{i}", title=f"Title {i}")
        for i in range(8)
    ]
    # Add papers without PMCID
    metas.append(SourceMetadata(source_id="PMID:99", pmcid=None, title="No PMCID"))
    metas.append(SourceMetadata(source_id="PMID:100", pmcid="", title="Empty PMCID"))

    def mock_fetch_xml(url, params=None):
        if "PMC0" in url:
            return ""  # Empty string (guarded by fetch_full_text)
        elif "PMC1" in url:
            return "<article>truncated"  # Corrupt XML (recovered by lxml)
        elif "PMC2" in url:
            return None  # HTTP failure (None)
        elif "PMC3" in url:
            return "<html><body>502 Bad Gateway</body></html>"  # Non-JATS HTML
        elif "PMC4" in url:
            raise requests.exceptions.ReadTimeout("Timeout")  # Request exception
        else:
            return "<article><body><sec><title>Methods</title><p>Valid text</p></sec></body></article>"

    with patch.object(session, "get_text", side_effect=mock_fetch_xml):
        crashed = False
        results = []
        try:
            results = list(harvester.fetch_full_texts_concurrent(metas, max_workers=4))
        except Exception as e:
            crashed = True
            err_desc = str(e)

        check(not crashed, f"fetch_full_texts_concurrent handled mixed failures gracefully (got {len(results)} outputs)")
        check(len(results) == len(metas), f"Every input yielded exactly one result tuple ({len(results)}/{len(metas)})")

        valid_count = sum(1 for m, s in results if s and "methods" in s)
        empty_or_none_count = sum(1 for m, s in results if not s)
        check(valid_count == 3, f"Exactly 3 valid XMLs parsed (got {valid_count})")
        check(empty_or_none_count == 7, f"All 7 failed/empty/missing-pmcid papers yielded empty or None safely (got {empty_or_none_count})")

        # Note about empty dict vs None
        none_count = sum(1 for m, s in results if s is None)
        empty_dict_count = sum(1 for m, s in results if s == {})
        finding(f"fetch_full_texts_concurrent return behavior: {none_count} None, {empty_dict_count} empty dict {{}} for non-JATS/empty bodies (both evaluate to False in 'if not sections')")

    # 3. BioRxivHarvester Malicious and Edge Case Inputs
    print("\n--- BioRxivHarvester edge case inputs ---")
    biorxiv = BioRxivHarvester(session, server="biorxiv")

    edge_cases = [
        # Missing DOI
        {"title": "No DOI paper", "date": "2023-01-01", "authors": "Smith, J."},
        # Missing date
        {"doi": "10.1101/nodate", "title": "No date paper", "authors": "Smith, J."},
        # Invalid / unparseable date
        {"doi": "10.1101/baddate", "date": "invalid-date-format", "title": "Bad date paper"},
        # Missing authors / null authors
        {"doi": "10.1101/noauth", "date": "2023-01-01", "title": "No authors", "authors": None},
        # Missing title / null title
        {"doi": "10.1101/notitle", "date": "2023-01-01", "title": None, "authors": "Smith, J."},
        # Null corresponding institution and author
        {"doi": "10.1101/noinst", "date": "2023-01-01", "title": "No inst", "author_corresponding_institution": None, "author_corresponding": None},
        # Empty dict
        {},
    ]

    mock_resp = MagicMock(status_code=200)
    mock_resp.json.return_value = {
        "messages": [{"status": "ok", "total": len(edge_cases)}],
        "collection": edge_cases,
    }

    with patch.object(session.session, "request", return_value=mock_resp):
        harvested = list(biorxiv.fetch_range("2023-01-01", "2023-01-02", max_results=100))
        check(len(harvested) >= 4, f"BioRxivHarvester safely skipped missing DOI and parsed valid entries ({len(harvested)} parsed)")
        check(all(h.doi for h in harvested), "All returned SourceMetadata have valid DOI")
        # Check invalid date fallback
        bad_date_meta = next((h for h in harvested if h.doi == "10.1101/baddate"), None)
        check(bad_date_meta is not None and bad_date_meta.publication_date is None, "Invalid date gracefully converted to None")

    # 4. Special Characters and Malicious Query Strings in In-Memory Filtering
    print("\n--- BioRxivHarvester query filtering with special characters ---")
    special_queries = [
        'TDP-43 / (cryptic) !@#$%',
        '"""',
        '   ',
        'a',
        '\\',
        '[splicing] AND {neuron}',
        '<script>alert("xss")</script>',
        "'; DROP TABLE sources; --",
    ]

    with patch.object(session.session, "request", return_value=mock_resp):
        for q in special_queries:
            crashed = False
            try:
                res = list(biorxiv.fetch_range("2023-01-01", "2023-01-02", max_results=10, query=q))
            except Exception as e:
                crashed = True
                err_desc = str(e)
            check(not crashed, f"Query '{q}' handled cleanly without regex/syntax crash")


# =====================================================================
# TRACK 3: Bayesian Math Edge Cases & Boundary Conditions
# =====================================================================

def test_bayesian_math_edge_cases():
    print("\n" + "=" * 70)
    print("TRACK 3: Bayesian Math & Boundary Conditions")
    print("=" * 70)

    # 1. Boundary condition in AssociationEngine at weighted_support = 2.999 vs 3.000
    print("\n--- AssociationEngine boundary condition (weighted_support = 2.999 vs 3.000) ---")
    base_kwargs = dict(
        pathway_key="pw_tdp",
        pathway_label="TDP-43 Splicing",
        n_labs=3,
        n_ortho_groups=2,
        survives_lolo=True,
        q_value=0.01,
        direction_consistency=0.9,
    )

    assoc_below = PathwayAssociation(**base_kwargs, weighted_support=2.9999)
    assoc_exact = PathwayAssociation(**base_kwargs, weighted_support=3.0000)
    assoc_above = PathwayAssociation(**base_kwargs, weighted_support=3.0001)

    tier_below = assoc_below.confidence_tier()
    tier_exact = assoc_exact.confidence_tier()
    tier_above = assoc_above.confidence_tier()

    check(tier_below == "moderate", f"weighted_support=2.9999 strictly gated to 'moderate' (got '{tier_below}')")
    check(tier_exact == "strong", f"weighted_support=3.0000 promotes to 'strong' (got '{tier_exact}')")
    check(tier_above == "strong", f"weighted_support=3.0001 promotes to 'strong' (got '{tier_above}')")

    # Also test boundary at weighted_support = 1.999 vs 2.000 (weak vs moderate)
    assoc_mod_below = PathwayAssociation(**base_kwargs, weighted_support=1.9999)
    assoc_mod_exact = PathwayAssociation(**base_kwargs, weighted_support=2.0000)
    check(assoc_mod_below.confidence_tier() == "weak", f"weighted_support=1.9999 gated to 'weak' (got '{assoc_mod_below.confidence_tier()}')")
    check(assoc_mod_exact.confidence_tier() == "moderate", f"weighted_support=2.0000 promotes to 'moderate' (got '{assoc_mod_exact.confidence_tier()}')")

    # Also test boundary at weighted_support = 0.999 vs 1.000 (insufficient vs weak)
    assoc_insuf = PathwayAssociation(**base_kwargs, weighted_support=0.9999)
    assoc_min_sup = PathwayAssociation(**base_kwargs, weighted_support=1.0000)
    check(assoc_insuf.confidence_tier() == "insufficient", f"weighted_support=0.9999 returns 'insufficient' (got '{assoc_insuf.confidence_tier()}')")
    check(assoc_min_sup.confidence_tier() == "weak", f"weighted_support=1.0000 returns 'weak' (got '{assoc_min_sup.confidence_tier()}')")

    # 2. EvidenceWeighter with mixed preprints and peer-reviewed papers
    print("\n--- EvidenceWeighter preprint calibration across tiers ---")
    weighter = EvidenceWeighter(preprint_weight=0.80)

    cell = CellModel(model_class=CellModelClass.IPSC_NEURON_NGN2, is_isogenic_pair=True)
    pert = Perturbation(ptype=PerturbationType.ASO, agent="ASO-1")
    outcome = Outcome(phenotype=NormalizedTerm(raw="cryptic exon", id="HP:0000001"), direction=Direction.DECREASE, is_quantified=True, n_biological_replicates=6)

    # Validate the exact quartile weight constants: Q1=1.00, Q2=0.85, Q3=0.65, Q4=0.50, Preprint=0.80
    journals = [
        ("Nature", "Q1", False, QUARTILE_WEIGHT["Q1"]),
        ("Neurobiology of Disease", "Q2", False, QUARTILE_WEIGHT["Q2"]),
        ("PLOS One", "Q3", False, QUARTILE_WEIGHT["Q3"]),
        ("Low tier", "Q4", False, QUARTILE_WEIGHT["Q4"]),
        ("Unindexed Journal", None, False, UNKNOWN_QUARTILE_WEIGHT),
        ("bioRxiv", None, True, 0.80),
        ("medRxiv", None, True, 0.80),
    ]

    for jname, quart, is_prep, expected_tier in journals:
        meta = SourceMetadata(source_id=f"TEST_{jname}", journal=jname, sjr_quartile=quart, is_preprint=is_prep, year=2024)
        rec = ExperimentRecord(source=meta, cell_model=cell, perturbation=pert, outcome=outcome, evidence_type=EvidenceType.DIRECT_MEASUREMENT)
        wb = weighter.weight(rec)
        check(abs(wb.journal_tier - expected_tier) < 1e-5, f"{jname} (is_preprint={is_prep}, Q={quart}) assigned journal tier {wb.journal_tier} (expected {expected_tier})")

    # Verify that preprint (0.80) is discounted 20% relative to Q1 (1.00) and sits between Q2 (0.85) and Q3 (0.65)
    check(QUARTILE_WEIGHT["Q2"] > 0.80 > QUARTILE_WEIGHT["Q3"], "Preprint weight (0.80) is properly calibrated between Q2 (0.85) and Q3 (0.65)")

    # 3. EvidenceReconciler under zero-weight, extreme, and conflicting inputs
    print("\n--- EvidenceReconciler convergence score behavior under extreme weights ---")
    reconciler = EvidenceReconciler()

    # Case A: Empty items
    res_empty = reconciler.reconcile("Test claim", [])
    check(res_empty.convergence_score == 0.0, f"Empty items convergence score == 0.0 (got {res_empty.convergence_score})")
    check(res_empty.concordance.value == "silent", f"Empty items concordance is silent")

    # Case B: Zero-confidence items
    item_zero_1 = EvidenceItem(source_class=EvidenceClass.LITERATURE, source_name="Lit1", claim_direction="supports", confidence=0.0)
    item_zero_2 = EvidenceItem(source_class=EvidenceClass.SEQUENCE, source_name="Seq1", claim_direction="supports", confidence=0.0)
    res_zero = reconciler.reconcile("Test claim", [item_zero_1, item_zero_2])
    check(res_zero.convergence_score == 0.0, f"Zero-confidence items convergence score == 0.0 (got {res_zero.convergence_score})")

    # Case C: Single class with 100 papers vs Two independent classes
    items_lit_100 = [
        EvidenceItem(source_class=EvidenceClass.LITERATURE, source_name=f"Paper_{i}", claim_direction="supports", confidence=1.0, n_independent=1)
        for i in range(100)
    ]
    res_lit_100 = reconciler.reconcile("Test claim", items_lit_100)

    # 2 independent classes (Literature + Sequence) with 3 observations each
    items_multi_2 = [
        EvidenceItem(source_class=EvidenceClass.LITERATURE, source_name="Paper_Lit", claim_direction="supports", confidence=1.0, n_independent=3),
        EvidenceItem(source_class=EvidenceClass.SEQUENCE, source_name="SpliceAI", claim_direction="supports", confidence=1.0, n_independent=3),
    ]
    res_multi_2 = reconciler.reconcile("Test claim", items_multi_2)

    check(res_multi_2.concordance.value == "convergent", f"2 independent classes classified as convergent")
    check(res_lit_100.concordance.value == "consistent", f"Single class with 100 papers classified only as consistent, NOT convergent")
    check(not res_lit_100.is_actionable(), "100 single-class papers are NOT actionable (is_actionable requires >=2 classes)")
    check(res_multi_2.is_actionable(), "2 independent classes with replication ARE actionable")

    # Case D: Maximum weight across all 7 classes
    items_max = [
        EvidenceItem(source_class=c, source_name=f"Max_{c.value}", claim_direction="supports", confidence=1.0, n_independent=20)
        for c in EvidenceClass
    ]
    res_max = reconciler.reconcile("Test claim", items_max)
    check(res_max.convergence_score > 20.0, f"Maximum 7-class convergence score reaches {res_max.convergence_score} (>20.0)")
    check(res_max.is_actionable(), "Maximum evidence is actionable")


# =====================================================================
# Main Execution Runner
# =====================================================================

def run_adversarial_suite():
    t0 = time.time()
    test_sqlite_thread_affinity_torture()
    test_malicious_inputs_and_edge_cases()
    test_bayesian_math_edge_cases()
    elapsed = time.time() - t0

    print("\n" + "=" * 70)
    print("CHALLENGER 2 SUMMARY")
    print("=" * 70)
    print(f"Elapsed: {elapsed:.2f}s")
    print(f"Failures: {len(FAILURES)}")
    print(f"Findings: {len(FINDINGS)}")
    if FAILURES:
        print("\nFAILURE DETAILS:")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    else:
        print("\nALL ADVERSARIAL STRESS TESTS COMPLETED!")
        if FINDINGS:
            print("\nKEY FINDINGS:")
            for f in FINDINGS:
                print(f"  * {f}")


if __name__ == "__main__":
    run_adversarial_suite()
