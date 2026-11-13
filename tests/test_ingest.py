"""
Comprehensive test suite for high-throughput literature harvesting, rate limiting,
bioRxiv/medRxiv preprints, and Bayesian preprint discounting.

100% offline, deterministic, mock-based execution (<3 seconds total).

Run: python -m tests.test_ingest
"""

from __future__ import annotations

import asyncio
import json
import math
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date
from pathlib import Path
from unittest.mock import MagicMock, patch

from synaptik.evidence.reconcile import (
    CLASS_PRIOR, EvidenceClass, EvidenceItem, EvidenceReconciler,
    from_literature_association,
)
from synaptik.infer.association import AssociationEngine, PathwayAssociation
from synaptik.infer.weighting import EvidenceWeighter
from synaptik.ingest.europepmc import EuropePMCHarvester
from synaptik.ingest.http import (
    ConcurrentRateLimiter, PoliteSession, RateLimiter, DEFAULT_RATE_LIMITS,
)
from synaptik.ingest.pubmed import BioRxivHarvester, PubMedHarvester
from synaptik.schema import (
    AssayClass, CellModel, CellModelClass, Direction, EvidenceType,
    ExperimentRecord, Outcome, Perturbation, PerturbationType,
    Provenance, SourceMetadata, NormalizedTerm,
)
from synaptik.store import Store

FAILURES: list[str] = []


def check(cond: bool, msg: str) -> None:
    if cond:
        print(f"  PASS  {msg}")
    else:
        print(f"  FAIL  {msg}")
        FAILURES.append(msg)


# =====================================================================
# 1. ConcurrentRateLimiter Unit Tests
# =====================================================================

def test_rate_limiter_instant_reservation() -> None:
    print("\n[rate limiter: instant reservation & multi-host isolation]")
    limiter = ConcurrentRateLimiter(limits={"eutils.ncbi.nlm.nih.gov": 0.34, "www.ebi.ac.uk": 0.15})

    # Test 10 concurrent threads calling acquire_delay instantly (zero sleeping)
    delays = []
    sched_times = []
    lock = threading.Lock()

    def acquire(host):
        t_call = time.monotonic()
        d = limiter.acquire_delay(host)
        with lock:
            delays.append(d)
            sched_times.append(t_call + d)

    threads = [threading.Thread(target=acquire, args=("eutils.ncbi.nlm.nih.gov",)) for _ in range(10)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    check(len(delays) == 10, "10 threads acquired reservations successfully")
    sorted_d = sorted(delays)
    check(sorted_d[0] == 0.0, f"first request has 0.0 delay (got {sorted_d[0]})")
    # Verify spacing between consecutive scheduled transmission slots is exactly 0.34s
    sorted_sched = sorted(sched_times)
    sched_spacings = [sorted_sched[i + 1] - sorted_sched[i] for i in range(len(sorted_sched) - 1)]
    valid_sched = all(abs(s - 0.34) < 1e-3 for s in sched_spacings)
    check(valid_sched, f"all reservation scheduled slots are exactly 0.34s (min {min(sched_spacings):.3f}, max {max(sched_spacings):.3f})")

    # Verify delays spacing allows slight thread lock dispatch jitter
    spacings = [sorted_d[i + 1] - sorted_d[i] for i in range(len(sorted_d) - 1)]
    valid_spacing = all(abs(s - 0.34) < 0.01 for s in spacings)
    check(valid_spacing, f"all reservation intervals are within 10ms of 0.34s (min {min(spacings):.3f}, max {max(spacings):.3f})")

    # Host isolation: EBI request must not be blocked or delayed by NCBI reservations
    ebi_delay = limiter.acquire_delay("www.ebi.ac.uk")
    check(ebi_delay == 0.0, f"independent host EBI acquires immediately with 0.0 delay (got {ebi_delay})")

    # Backward compatibility checks
    check(issubclass(RateLimiter, ConcurrentRateLimiter), "RateLimiter is subclass of ConcurrentRateLimiter")
    rl = RateLimiter()
    check(isinstance(rl, ConcurrentRateLimiter), "RateLimiter instance is ConcurrentRateLimiter")


def test_rate_limiter_concurrent_execution() -> None:
    print("\n[rate limiter: multi-threaded concurrent timing & host isolation]")
    interval = 0.02
    limiter = ConcurrentRateLimiter(limits={"fast.test": interval, "slow.host": 0.40, "quick.host": 0.01})

    timestamps = []
    lock = threading.Lock()

    def worker():
        limiter.wait("fast.test")
        with lock:
            timestamps.append(time.monotonic())

    threads = [threading.Thread(target=worker) for _ in range(10)]
    t0 = time.monotonic()
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    total_elapsed = time.monotonic() - t0

    check(len(timestamps) == 10, "10 concurrent worker threads completed wait()")
    sorted_ts = sorted(timestamps)
    diffs = [sorted_ts[i + 1] - sorted_ts[i] for i in range(len(sorted_ts) - 1)]
    check(min(diffs) >= interval * 0.75, f"minimum inter-request interval satisfied ({min(diffs):.4f}s >= {interval * 0.75:.4f}s)")
    check(total_elapsed < 1.0, f"10-thread test completed fast ({total_elapsed:.2f}s < 1.0s)")

    # Cross-host isolation under real wait()
    quick_times = []

    def quick_worker():
        limiter.wait("quick.host")
        with lock:
            quick_times.append(time.monotonic())

    def slow_worker():
        limiter.wait("slow.host")

    t_slow = threading.Thread(target=slow_worker)
    t_slow.start()
    t_quick_start = time.monotonic()
    q_threads = [threading.Thread(target=quick_worker) for _ in range(4)]
    for qt in q_threads:
        qt.start()
    for qt in q_threads:
        qt.join()
    quick_duration = time.monotonic() - t_quick_start
    t_slow.join()

    check(len(quick_times) == 4, "4 quick host requests completed")
    check(quick_duration < 0.25, f"quick host requests completed in {quick_duration:.3f}s without slow host stall")


def test_rate_limiter_async_wait() -> None:
    print("\n[rate limiter: async_wait compatibility]")
    limiter = ConcurrentRateLimiter(limits={"async.test": 0.02})

    async def async_worker(results: list[float]):
        await limiter.async_wait("async.test")
        results.append(time.monotonic())

    async def main_async():
        results: list[float] = []
        await asyncio.gather(
            async_worker(results),
            async_worker(results),
            async_worker(results),
        )
        return results

    t0 = time.monotonic()
    res = asyncio.run(main_async())
    elapsed = time.monotonic() - t0

    check(len(res) == 3, "3 async tasks completed async_wait")
    check(elapsed >= 0.035, f"async tasks spaced properly ({elapsed:.3f}s >= 0.035s)")


# =====================================================================
# 2. PoliteSession Thread Safety, Backoff & Batching
# =====================================================================

def test_polite_session_thread_safety_and_batch() -> None:
    print("\n[polite session: thread safety and metrics]")
    session = PoliteSession(store=None, rate_limits={"api.test": 0.001})

    # Mock response
    mock_resp = MagicMock()
    mock_resp.status_code = 200
    mock_resp.json.return_value = {"status": "ok"}

    with patch.object(session.session, "request", return_value=mock_resp):
        urls_with_params = [(f"https://api.test/{i}", {"p": i}) for i in range(20)]
        results = session.fetch_batch(urls_with_params, want_json=True, max_workers=5)

        check(len(results) == 20, f"fetch_batch returned all 20 responses (got {len(results)})")
        check(session.n_requests == 20, f"n_requests counter atomically reached 20 (got {session.n_requests})")
        check(session.n_failures == 0, f"zero failures recorded (got {session.n_failures})")

    # Metric thread-safety under concurrent increment
    def inc_workers():
        for _ in range(50):
            session._inc_requests()
            session._inc_cache_hits()
            session._inc_failures()

    threads = [threading.Thread(target=inc_workers) for _ in range(10)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    st = session.stats()
    check(st["requests"] == 20 + 500, f"atomic requests metric exact ({st['requests']} == 520)")
    check(st["cache_hits"] == 500, f"atomic cache_hits metric exact ({st['cache_hits']} == 500)")
    check(st["failures"] == 500, f"atomic failures metric exact ({st['failures']} == 500)")


def test_polite_session_caching_and_batch() -> None:
    print("\n[polite session: two-stage caching in fetch_batch]")
    store = Store(":memory:")
    session = PoliteSession(store=store, rate_limits={"api.test": 0.001})

    # Pre-populate cache for 2 items
    key0 = session._cache_key("GET", "https://api.test/item0", None)
    key1 = session._cache_key("GET", "https://api.test/item1", None)
    store.cache_put(key0, {"id": 0, "source": "cache"})
    store.cache_put(key1, {"id": 1, "source": "cache"})

    batch_items = [
        "https://api.test/item0",  # cache hit
        "https://api.test/item1",  # cache hit
        "https://api.test/item2",  # network fetch
        "https://api.test/item3",  # network fetch
        "https://api.test/item4",  # network fetch
    ]

    def mock_request(method, url, params=None, timeout=None):
        resp = MagicMock()
        resp.status_code = 200
        idx = int(url.split("item")[-1])
        resp.json.return_value = {"id": idx, "source": "network"}
        return resp

    with patch.object(session.session, "request", side_effect=mock_request):
        results = session.fetch_batch(batch_items, want_json=True, max_workers=3)

        check(len(results) == 5, f"fetch_batch returned all 5 items (got {len(results)})")
        check(results[0] == {"id": 0, "source": "cache"}, "item 0 served from cache")
        check(results[1] == {"id": 1, "source": "cache"}, "item 1 served from cache")
        check(results[2] == {"id": 2, "source": "network"}, "item 2 fetched from network")
        check(results[3] == {"id": 3, "source": "network"}, "item 3 fetched from network")
        check(results[4] == {"id": 4, "source": "network"}, "item 4 fetched from network")

        check(session.n_cache_hits == 2, f"n_cache_hits == 2 (got {session.n_cache_hits})")
        check(session.n_requests == 3, f"n_requests == 3 (got {session.n_requests})")

        # Verify newly fetched items were saved into cache
        key2 = session._cache_key("GET", "https://api.test/item2", None)
        cached2 = store.cache_get(key2)
        check(cached2 == {"id": 2, "source": "network"}, "item 2 was persisted into cache")


def test_polite_session_backoff_and_retry() -> None:
    print("\n[polite session: retries and exponential backoff]")
    session = PoliteSession(store=None, max_retries=3, rate_limits={"api.test": 0.001, "blocked.test": 0.001})

    resp_429 = MagicMock(status_code=429)
    resp_503 = MagicMock(status_code=503)
    resp_200 = MagicMock(status_code=200)
    resp_200.json.return_value = {"recovered": True}

    with patch.object(session.session, "request", side_effect=[resp_429, resp_503, resp_200]), \
         patch("synaptik.ingest.http.time.sleep") as mock_sleep:
        data = session.get_json("https://api.test/retry")
        check(data == {"recovered": True}, "session recovered after 429 and 503 retries")
        check(mock_sleep.call_count >= 2, f"exponential backoff slept between retries (count {mock_sleep.call_count})")

    # 403 egress denial test
    resp_403 = MagicMock(status_code=403)
    resp_403.headers = {"x-deny-reason": "Sandbox Blocked"}
    with patch.object(session.session, "request", return_value=resp_403):
        res403 = session.get_json("https://blocked.test/egress")
        check(res403 is None, "403 sandbox block returns None cleanly")
        check(session.n_failures >= 1, "failure recorded for sandbox block")


# =====================================================================
# 3. Europe PMC Concurrent Full-Text & JATS XML Parsing
# =====================================================================

SAMPLE_JATS_XML = """<?xml version="1.0" encoding="UTF-8"?>
<article>
  <front>
    <article-meta>
      <abstract>
        <p>Loss of TDP-43 induces cryptic exon inclusion in STMN2 and impaired axonal regrowth.</p>
      </abstract>
    </article-meta>
  </front>
  <body>
    <sec>
      <title>Introduction</title>
      <p>Amyotrophic lateral sclerosis is a fatal neurodegenerative disorder.</p>
    </sec>
    <sec>
      <title>Materials and Methods</title>
      <p>Human iPSC-derived motor neurons were differentiated for 21 DIV using NGN2 induction.</p>
    </sec>
    <sec>
      <title>Results</title>
      <p>Cryptic exon inclusion in STMN2 was rescued by splice-switching ASO treatment.</p>
      <fig id="fig1">
        <label>Figure 1</label>
        <caption><p>Dose-dependent rescue of STMN2 protein levels (n=6 biological replicates, p &lt; 0.001).</p></caption>
      </fig>
    </sec>
    <sec>
      <title>Discussion</title>
      <p>These data provide strong rationale for clinical translation.</p>
    </sec>
  </body>
</article>
"""

def test_europepmc_jats_xml_and_concurrent_fetch() -> None:
    print("\n[europe pmc: JATS XML parsing & concurrent extraction]")
    sections = EuropePMCHarvester.parse_jats(SAMPLE_JATS_XML)
    check("abstract" in sections, "JATS XML parser extracted abstract")
    check("methods" in sections, "JATS XML parser extracted methods")
    check("results" in sections, "JATS XML parser extracted results")
    check("figures" in sections, "JATS XML parser extracted figure captions")
    check("n=6 biological replicates" in sections["figures"], "figure caption retained replicate count")

    # Concurrent fetch test
    session = PoliteSession(store=None, rate_limits={"www.ebi.ac.uk": 0.001})
    harvester = EuropePMCHarvester(session)
    mock_resp = MagicMock(status_code=200, text=SAMPLE_JATS_XML)

    metas = [
        SourceMetadata(source_id=f"PMID:100{i}", pmcid=f"PMC100{i}", title=f"Paper {i}")
        for i in range(8)
    ]
    # Add 2 papers without PMCID
    metas.append(SourceMetadata(source_id="PMID:2001", pmcid=None, title="Paper No PMC 1"))
    metas.append(SourceMetadata(source_id="PMID:2002", pmcid=None, title="Paper No PMC 2"))

    with patch.object(session.session, "request", return_value=mock_resp):
        res = list(harvester.fetch_full_texts_concurrent(metas, max_workers=4))
        check(len(res) == 10, f"concurrent fetch retrieved all 10 papers (got {len(res)})")
        with_text = [r for r in res if r[1] is not None]
        without_text = [r for r in res if r[1] is None]
        check(len(with_text) == 8, f"8 papers with PMCIDs extracted full text (got {len(with_text)})")
        check(len(without_text) == 2, f"2 papers without PMCID returned None immediately (got {len(without_text)})")
        check(all("methods" in r[1] for r in with_text), "all extracted full texts contain methods section")


# =====================================================================
# 4. bioRxiv and medRxiv Harvester & Keyword Filtering
# =====================================================================

SAMPLE_CSHL_RESPONSE = {
    "messages": [{"status": "ok", "total": 2}],
    "collection": [
        {
            "doi": "10.1101/2023.05.01.123456",
            "title": "TDP-43 cryptic exon inclusion in human motor neurons",
            "authors": "Smith, J.; Doe, A.",
            "author_corresponding": "Jane Smith",
            "author_corresponding_institution": "Harvard Medical School",
            "date": "2023-05-01",
            "version": "1",
            "category": "neuroscience",
            "abstract": "We describe cryptic splicing events in ALS patient iPSC models.",
        },
        {
            "doi": "10.1101/2023.05.02.654321",
            "title": "Cryo-EM structure of human cardiac sodium channel",
            "authors": "Brown, C.; White, D.",
            "author_corresponding": "Charlie Brown",
            "author_corresponding_institution": "Stanford University",
            "date": "2023-05-02",
            "version": "1",
            "category": "biophysics",
            "abstract": "Structure of Nav1.5 in lipid nanodiscs reveals channel inactivation.",
        },
    ],
}

def test_biorxiv_medrxiv_parsing_and_filtering() -> None:
    print("\n[biorxiv / medrxiv: parsing and in-memory keyword filtering]")
    session = PoliteSession(store=None, rate_limits={"api.biorxiv.org": 0.001})
    biorxiv = BioRxivHarvester(session, server="biorxiv")
    medrxiv = BioRxivHarvester(session, server="medrxiv")

    check(biorxiv.server == "biorxiv", "biorxiv server configured")
    check(medrxiv.server == "medrxiv", "medrxiv server configured")

    mock_resp = MagicMock(status_code=200)
    mock_resp.json.return_value = SAMPLE_CSHL_RESPONSE

    with patch.object(session.session, "request", return_value=mock_resp):
        # 1. Unfiltered query
        all_results = list(biorxiv.fetch_range("2023-05-01", "2023-05-02", max_results=10))
        check(len(all_results) == 2, f"unfiltered harvest fetched all 2 records (got {len(all_results)})")
        check(all_results[0].is_preprint is True, "is_preprint flag is True")
        check(all_results[0].journal == "biorxiv (preprint)", "journal set to 'biorxiv (preprint)'")

        # 2. medRxiv test
        med_results = list(medrxiv.fetch_range("2023-05-01", "2023-05-02", max_results=10))
        check(len(med_results) == 2, "medrxiv fetched 2 records")
        check(med_results[0].journal == "medrxiv (preprint)", "journal set to 'medrxiv (preprint)'")

        # 3. In-memory keyword filtering (filtering for "TDP-43 cryptic")
        filtered = list(biorxiv.fetch_range("2023-05-01", "2023-05-02", max_results=10, query="TDP-43 cryptic"))
        check(len(filtered) == 1, f"keyword filtered harvest returned exactly 1 matching record (got {len(filtered)})")
        check("TDP-43" in filtered[0].title, "returned record matches title keyword")

        # 4. Filter with no match
        empty_res = list(biorxiv.fetch_range("2023-05-01", "2023-05-02", max_results=10, query="nonexistent_target_xyz"))
        check(len(empty_res) == 0, "non-matching query returns 0 records")


# =====================================================================
# 5. PubMed Concurrent Metadata Batching
# =====================================================================

SAMPLE_EFETCH_XML = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//DTD PubMedArticle, 1st January 2024//EN" "https://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_240101.dtd">
<PubmedArticleSet>
  <PubmedArticle>
    <MedlineCitation>
      <PMID>31000001</PMID>
      <Article>
        <ArticleTitle>Rescue of STMN2 splicing in human motor neurons.</ArticleTitle>
        <Journal>
          <Title>Neuron</Title>
          <ISOAbbreviation>Neuron</ISOAbbreviation>
        </Journal>
        <AuthorList>
          <Author>
            <LastName>Smith</LastName>
            <ForeName>Jane</ForeName>
            <AffiliationInfo><Affiliation>Harvard</Affiliation></AffiliationInfo>
          </Author>
        </AuthorList>
        <PublicationTypeList>
          <PublicationType>Journal Article</PublicationType>
        </PublicationTypeList>
      </Article>
    </MedlineCitation>
  </PubmedArticle>
</PubmedArticleSet>
"""

def test_pubmed_concurrent_metadata() -> None:
    print("\n[pubmed: concurrent metadata batching]")
    session = PoliteSession(store=None, rate_limits={"eutils.ncbi.nlm.nih.gov": 0.001})
    pubmed = PubMedHarvester(session)

    pmids = [str(31000000 + i) for i in range(15)]

    def mock_request(method, url, params=None, timeout=None):
        resp = MagicMock()
        resp.status_code = 200
        # return sample XML with requested IDs
        ids_str = (params or {}).get("id", "")
        requested_ids = ids_str.split(",") if ids_str else []
        articles = []
        for pid in requested_ids:
            articles.append(f"""
            <PubmedArticle>
              <MedlineCitation>
                <PMID>{pid}</PMID>
                <Article>
                  <ArticleTitle>Article title for PMID {pid}</ArticleTitle>
                  <Journal><Title>Neuron</Title><ISOAbbreviation>Neuron</ISOAbbreviation></Journal>
                </Article>
              </MedlineCitation>
            </PubmedArticle>
            """)
        xml = f"<PubmedArticleSet>{''.join(articles)}</PubmedArticleSet>"
        resp.text = xml
        return resp

    with patch.object(session.session, "request", side_effect=mock_request):
        # 15 PMIDs in batches of 5 across 3 workers
        results = list(pubmed.fetch_metadata(pmids, batch_size=5, max_workers=3))
        check(len(results) == 15, f"concurrent fetch_metadata returned 15 records (got {len(results)})")
        check(results[0].pmid == "31000000", f"first PMID preserved order ({results[0].pmid})")
        check(results[-1].pmid == "31000014", f"last PMID preserved order ({results[-1].pmid})")


# =====================================================================
# 6. Bayesian Preprint Discount Calibration Verification
# =====================================================================

def test_preprint_discount_bayesian_pipeline() -> None:
    print("\n[bayesian calibration: preprint discount propagation]")
    weighter = EvidenceWeighter(preprint_weight=0.80)

    # 1. EvidenceWeighter test: 20% discount
    meta_q1 = SourceMetadata(source_id="PMID:111", journal="Neuron", sjr_quartile="Q1", is_preprint=False, year=2024)
    meta_pre = SourceMetadata(source_id="DOI:10.1101/222", journal="biorxiv (preprint)", is_preprint=True, year=2024)

    cell = CellModel(model_class=CellModelClass.IPSC_NEURON_NGN2, is_isogenic_pair=True)
    pert = Perturbation(ptype=PerturbationType.ASO, agent="ASO-1")
    outcome = Outcome(phenotype=NormalizedTerm(raw="cryptic exon", id="HP:0000001"), direction=Direction.DECREASE, is_quantified=True, n_biological_replicates=6)

    rec_q1 = ExperimentRecord(source=meta_q1, cell_model=cell, perturbation=pert, outcome=outcome, evidence_type=EvidenceType.DIRECT_MEASUREMENT)
    rec_pre = ExperimentRecord(source=meta_pre, cell_model=cell, perturbation=pert, outcome=outcome, evidence_type=EvidenceType.DIRECT_MEASUREMENT)

    wb_q1 = weighter.weight(rec_q1)
    wb_pre = weighter.weight(rec_pre)

    check(wb_q1.journal_tier == 1.00, f"Q1 journal tier weight is 1.00 (got {wb_q1.journal_tier})")
    check(wb_pre.journal_tier == 0.80, f"preprint journal tier weight is 0.80 (got {wb_pre.journal_tier})")
    check(abs(wb_pre.composite - 0.80 * wb_q1.composite) < 1e-4, "composite weight has exact 20% preprint discount")

    # 2. AssociationEngine tier gating test
    # 3 Q1 papers (total weight ~ 3.0) reach 'strong'; 3 preprints (total weight ~ 2.4) remain 'moderate'
    assoc_q1 = PathwayAssociation(pathway_key="pw_1", pathway_label="Path 1", n_labs=3, n_ortho_groups=2, survives_lolo=True, q_value=0.01, direction_consistency=0.9, weighted_support=3.0)
    assoc_pre = PathwayAssociation(pathway_key="pw_1", pathway_label="Path 1", n_labs=3, n_ortho_groups=2, survives_lolo=True, q_value=0.01, direction_consistency=0.9, weighted_support=2.4)

    check(assoc_q1.confidence_tier() == "strong", f"Q1 association tiered 'strong' (support {assoc_q1.weighted_support})")
    check(assoc_pre.confidence_tier() == "moderate", f"preprint association gated to 'moderate' (support {assoc_pre.weighted_support})")

    # 3. EvidenceReconciler convergence score test
    reconciler = EvidenceReconciler()
    item_q1 = from_literature_association(assoc_q1)
    item_pre = from_literature_association(assoc_pre)

    check(item_q1.confidence == 0.90, f"strong tier produces confidence 0.90 (got {item_q1.confidence})")
    check(item_pre.confidence == 0.70, f"moderate tier produces confidence 0.70 (got {item_pre.confidence})")

    # Pair with orthogonal evidence (e.g. sequence prediction)
    item_seq = EvidenceItem(source_class=EvidenceClass.SEQUENCE, source_name="SpliceAI", claim_direction="supports", confidence=0.80)

    res_q1 = reconciler.reconcile("STMN2 rescue", [item_q1, item_seq])
    res_pre = reconciler.reconcile("STMN2 rescue", [item_pre, item_seq])

    check(res_q1.convergence_score > res_pre.convergence_score, f"peer-reviewed convergence score ({res_q1.convergence_score}) exceeds preprint score ({res_pre.convergence_score})")


# =====================================================================
# 7. Git Cleanliness: .gitignore Test
# =====================================================================

def test_git_cleanliness_claude_md() -> None:
    print("\n[git cleanliness: CLAUDE.md exclusion]")
    gitignore_path = Path(__file__).parent.parent / ".gitignore"
    check(gitignore_path.exists(), ".gitignore file exists")
    lines = gitignore_path.read_text(encoding="utf-8").splitlines()
    has_claude = any(line.strip() == "CLAUDE.md" for line in lines)
    check(has_claude, ".gitignore explicitly contains 'CLAUDE.md'")


# =====================================================================
# Main Test Runner
# =====================================================================

def run_tests() -> None:
    print("=" * 70)
    print("SYNAPTIK Milestone 1 Ingest, Rate Limiter & Preprint Test Suite")
    print("=" * 70)
    t0 = time.time()

    test_rate_limiter_instant_reservation()
    test_rate_limiter_concurrent_execution()
    test_rate_limiter_async_wait()
    test_polite_session_thread_safety_and_batch()
    test_polite_session_caching_and_batch()
    test_polite_session_backoff_and_retry()
    test_europepmc_jats_xml_and_concurrent_fetch()
    test_biorxiv_medrxiv_parsing_and_filtering()
    test_pubmed_concurrent_metadata()
    test_preprint_discount_bayesian_pipeline()
    test_git_cleanliness_claude_md()

    elapsed = time.time() - t0
    print("\n" + "=" * 70)
    if FAILURES:
        print(f"FAILED: {len(FAILURES)} check(s) failed in {elapsed:.2f}s")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    else:
        print(f"ALL INGEST AND PREPRINT TESTS PASSED ({elapsed:.2f}s)")
        print("=" * 70)


if __name__ == "__main__":
    run_tests()
