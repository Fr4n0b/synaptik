"""
Shared HTTP plumbing for all harvesters.

Every external API this system touches has rate limits and intermittent
failures. Rather than sprinkle retry logic through each harvester, all network
egress goes through PoliteSession, which enforces per-host rate limits,
retries with exponential backoff on transient errors, and optionally caches
responses in the Store so a re-run costs nothing.
"""

from __future__ import annotations

import asyncio
import concurrent.futures
import hashlib
import logging
import random
import threading
import time
from typing import Any

import requests
from requests.adapters import HTTPAdapter

log = logging.getLogger(__name__)

# NCBI asks for <=3 req/s without a key (0.34s interval), 10/s with one (0.10s).
# Europe PMC is conservative at 0.15s (~6.67 req/s). bioRxiv is 0.5s (2 req/s).
DEFAULT_RATE_LIMITS: dict[str, float] = {
    "eutils.ncbi.nlm.nih.gov": 0.34,
    "www.ebi.ac.uk": 0.15,
    "api.biorxiv.org": 0.5,
    "api.crossref.org": 0.1,
    "www.ebi.ac.uk/chembl": 0.2,
    "api.platform.opentargets.org": 0.2,
    "reactome.org": 0.2,
    "rest.uniprot.org": 0.2,
    "gnomad.broadinstitute.org": 0.25,
    "_default": 0.25,
}

TRANSIENT_STATUS = {429, 500, 502, 503, 504}


class ConcurrentRateLimiter:
    """Thread-safe and async-compatible per-host rate limiter.

    Uses scheduled timestamp reservations so threads hold the lock for
    microseconds, sleeping independently outside the lock without blocking
    other hosts or threads.
    """

    def __init__(self, limits: dict[str, float] | None = None):
        self.limits: dict[str, float] = dict(limits or DEFAULT_RATE_LIMITS)
        self._next_allowed: dict[str, float] = {}
        self._lock = threading.Lock()

    def _normalize_host(self, host: str) -> str:
        if "://" in host:
            host = host.split("/")[2]
        if ":" in host:
            host = host.split(":")[0]
        return host.strip().lower()

    def acquire_delay(self, host: str) -> float:
        """Reserve a time slot for the given host and return delay in seconds."""
        h = self._normalize_host(host)
        interval = self.limits.get(h, self.limits.get("_default", 0.25))
        with self._lock:
            now = time.monotonic()
            scheduled = max(now, self._next_allowed.get(h, 0.0))
            self._next_allowed[h] = scheduled + interval
            return max(0.0, scheduled - now)

    def wait(self, host: str) -> None:
        """Synchronously wait for the host rate limit delay outside the lock."""
        delay = self.acquire_delay(host)
        if delay > 0:
            time.sleep(delay)

    async def async_wait(self, host: str) -> None:
        """Asynchronously wait for the host rate limit delay outside the lock."""
        delay = self.acquire_delay(host)
        if delay > 0:
            await asyncio.sleep(delay)


class RateLimiter(ConcurrentRateLimiter):
    """Backward compatibility alias for ConcurrentRateLimiter."""
    pass


class PoliteSession:
    """Thread-safe requests.Session with rate limiting, retries, and optional caching."""

    def __init__(
        self,
        store: Any = None,
        contact_email: str | None = None,
        tool_name: str = "SYNAPTIK",
        max_retries: int = 4,
        timeout: float = 45.0,
        use_cache: bool = True,
        rate_limits: dict[str, float] | None = None,
    ):
        self.store = store
        self.contact_email = contact_email
        self.tool_name = tool_name
        self.max_retries = max_retries
        self.timeout = timeout
        self.use_cache = use_cache and store is not None
        self.limiter = ConcurrentRateLimiter(rate_limits)
        self.session = requests.Session()

        # Mount robust connection pool adapter for multi-threaded harvesting
        adapter = HTTPAdapter(pool_connections=25, pool_maxsize=25, max_retries=0)
        self.session.mount("https://", adapter)
        self.session.mount("http://", adapter)

        ua = f"{tool_name}/1.0 (research tool"
        if contact_email:
            ua += f"; mailto:{contact_email}"
        ua += ")"
        self.session.headers.update({"User-Agent": ua, "Accept": "application/json"})

        self._stats_lock = threading.Lock()
        self._cache_lock = threading.Lock()
        self.n_requests = 0
        self.n_cache_hits = 0
        self.n_failures = 0

    @staticmethod
    def _cache_key(method: str, url: str, params: dict | None) -> str:
        basis = f"{method}|{url}|{sorted((params or {}).items())}"
        return "http:" + hashlib.sha1(basis.encode()).hexdigest()

    def _inc_requests(self) -> None:
        with self._stats_lock:
            self.n_requests += 1

    def _inc_cache_hits(self) -> None:
        with self._stats_lock:
            self.n_cache_hits += 1

    def _inc_failures(self) -> None:
        with self._stats_lock:
            self.n_failures += 1

    def _cache_get(self, key: str) -> dict | None:
        if not (self.use_cache and self.store is not None):
            return None
        try:
            with self._cache_lock:
                return self.store.cache_get(key)
        except Exception as e:
            log.debug("Cache get bypassed or failed: %s", e)
            return None

    def _cache_put(self, key: str, payload: dict) -> None:
        if not (self.use_cache and self.store is not None):
            return
        try:
            with self._cache_lock:
                self.store.cache_put(key, payload)
        except Exception as e:
            log.debug("Cache put bypassed or failed: %s", e)

    def get_json(
        self, url: str, params: dict | None = None, cacheable: bool = True
    ) -> dict | None:
        """GET returning parsed JSON, or None on unrecoverable failure.

        Returning None rather than raising is deliberate: a harvest over
        thousands of IDs must not abort because one record 404s.
        """
        key = self._cache_key("GET", url, params)
        if self.use_cache and cacheable:
            cached = self._cache_get(key)
            if cached is not None:
                self._inc_cache_hits()
                return cached

        host = url.split("/")[2] if "://" in url else "_default"
        payload = self._request_with_retry("GET", url, params, host, want_json=True)
        if payload is not None and self.use_cache and cacheable:
            self._cache_put(key, payload)
        return payload

    def get_text(self, url: str, params: dict | None = None) -> str | None:
        """GET returning raw text response, or None on failure."""
        host = url.split("/")[2] if "://" in url else "_default"
        return self._request_with_retry("GET", url, params, host, want_json=False)

    def fetch_batch(
        self,
        items: list[tuple[str, dict] | str],
        want_json: bool = True,
        max_workers: int = 5,
        cacheable: bool = True,
    ) -> list[Any]:
        """Fetch multiple URLs concurrently using a ThreadPoolExecutor.

        Args:
            items: List of (url, params_dict) tuples or URL strings.
            want_json: If True, returns parsed JSON; if False, returns raw text.
            max_workers: Max concurrent worker threads (default 5).
            cacheable: Whether to check and populate cache.

        Returns:
            List of responses in the exact order of the input items.
        """
        if not items:
            return []

        norm: list[tuple[str, dict | None]] = []
        for item in items:
            if isinstance(item, str):
                norm.append((item, None))
            elif isinstance(item, (tuple, list)):
                u = item[0]
                p = item[1] if len(item) > 1 else None
                norm.append((u, p))
            else:
                raise TypeError(f"Invalid batch item format: {item!r}")

        results: list[Any] = [None] * len(norm)
        pending: list[tuple[int, str, dict | None]] = []

        if self.use_cache and cacheable and want_json:
            for idx, (url, params) in enumerate(norm):
                key = self._cache_key("GET", url, params)
                cached = self._cache_get(key)
                if cached is not None:
                    self._inc_cache_hits()
                    results[idx] = cached
                else:
                    pending.append((idx, url, params))
        else:
            for idx, (url, params) in enumerate(norm):
                pending.append((idx, url, params))

        if not pending:
            return results

        def _fetch_worker(entry: tuple[int, str, dict | None]) -> tuple[int, Any]:
            idx, url, params = entry
            host = url.split("/")[2] if "://" in url else "_default"
            res = self._request_with_retry(
                "GET", url, params, host, want_json=want_json
            )
            return idx, res

        workers = min(max(1, max_workers), len(pending))
        if workers == 1:
            fetched = [_fetch_worker(entry) for entry in pending]
        else:
            with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
                fetched = list(pool.map(_fetch_worker, pending))

        for idx, res in fetched:
            results[idx] = res
            if res is not None and self.use_cache and cacheable and want_json:
                url, params = norm[idx]
                key = self._cache_key("GET", url, params)
                self._cache_put(key, res)

        return results

    def _request_with_retry(
        self, method: str, url: str, params: dict | None, host: str, want_json: bool
    ) -> Any:
        last_err: str | None = None
        for attempt in range(self.max_retries):
            self.limiter.wait(host)
            try:
                self._inc_requests()
                resp = self.session.request(
                    method, url, params=params, timeout=self.timeout
                )
                if resp.status_code == 200:
                    if not want_json:
                        return resp.text
                    try:
                        return resp.json()
                    except ValueError:
                        log.warning("Non-JSON 200 from %s", url)
                        return None
                if resp.status_code == 404:
                    return None
                if resp.status_code == 403:
                    deny = resp.headers.get("x-deny-reason")
                    if deny:
                        log.error(
                            "Network egress blocked for %s (%s). If running in a "
                            "sandbox, allow this host in network settings.", host, deny
                        )
                        self._inc_failures()
                        return None
                    last_err = "403 Forbidden"
                elif resp.status_code in TRANSIENT_STATUS:
                    last_err = f"HTTP {resp.status_code}"
                else:
                    log.warning("HTTP %s from %s", resp.status_code, url)
                    self._inc_failures()
                    return None
            except requests.RequestException as e:
                last_err = str(e)

            backoff = (2**attempt) + random.uniform(0, 0.5)
            log.debug("Retry %d/%d for %s after %s (sleep %.1fs)",
                      attempt + 1, self.max_retries, url, last_err, backoff)
            time.sleep(backoff)

        log.warning("Giving up on %s after %d attempts: %s",
                    url, self.max_retries, last_err)
        self._inc_failures()
        return None

    def stats(self) -> dict:
        with self._stats_lock:
            return {
                "requests": self.n_requests,
                "cache_hits": self.n_cache_hits,
                "failures": self.n_failures,
            }
