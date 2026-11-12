"""
Shared HTTP plumbing for all harvesters.

Every external API this system touches has rate limits and intermittent
failures. Rather than sprinkle retry logic through each harvester, all network
egress goes through PoliteSession, which enforces per-host rate limits,
retries with exponential backoff on transient errors, and optionally caches
responses in the Store so a re-run costs nothing.
"""

from __future__ import annotations

import hashlib
import logging
import random
import threading
import time
from typing import Any

import requests

log = logging.getLogger(__name__)

# NCBI asks for <=3 req/s without a key, 10/s with one. Europe PMC is more
# permissive but we stay conservative -- this pipeline runs overnight anyway.
DEFAULT_RATE_LIMITS: dict[str, float] = {
    "eutils.ncbi.nlm.nih.gov": 0.34,
    "www.ebi.ac.uk": 0.15,
    "api.biorxiv.org": 0.5,
    "api.crossref.org": 0.1,
    "www.ebi.ac.uk/chembl": 0.2,
    "api.platform.opentargets.org": 0.2,
    "reactome.org": 0.2,
    "rest.uniprot.org": 0.2,
    "_default": 0.25,
}

TRANSIENT_STATUS = {429, 500, 502, 503, 504}


class RateLimiter:
    """Per-host minimum interval between requests. Thread-safe."""

    def __init__(self, limits: dict[str, float] | None = None):
        self.limits = dict(limits or DEFAULT_RATE_LIMITS)
        self._last: dict[str, float] = {}
        self._lock = threading.Lock()

    def wait(self, host: str) -> None:
        interval = self.limits.get(host, self.limits.get("_default", 0.25))
        with self._lock:
            last = self._last.get(host, 0.0)
            delta = time.monotonic() - last
            if delta < interval:
                time.sleep(interval - delta)
            self._last[host] = time.monotonic()


class PoliteSession:
    """requests.Session with rate limiting, retries, and optional caching."""

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
        self.limiter = RateLimiter(rate_limits)
        self.session = requests.Session()
        ua = f"{tool_name}/1.0 (research tool"
        if contact_email:
            ua += f"; mailto:{contact_email}"
        ua += ")"
        self.session.headers.update({"User-Agent": ua, "Accept": "application/json"})
        self.n_requests = 0
        self.n_cache_hits = 0
        self.n_failures = 0

    @staticmethod
    def _cache_key(method: str, url: str, params: dict | None) -> str:
        basis = f"{method}|{url}|{sorted((params or {}).items())}"
        return "http:" + hashlib.sha1(basis.encode()).hexdigest()

    def get_json(
        self, url: str, params: dict | None = None, cacheable: bool = True
    ) -> dict | None:
        """GET returning parsed JSON, or None on unrecoverable failure.

        Returning None rather than raising is deliberate: a harvest over
        thousands of IDs must not abort because one record 404s.
        """
        key = self._cache_key("GET", url, params)
        if self.use_cache and cacheable:
            cached = self.store.cache_get(key)
            if cached is not None:
                self.n_cache_hits += 1
                return cached

        host = url.split("/")[2] if "://" in url else "_default"
        payload = self._request_with_retry("GET", url, params, host, want_json=True)
        if payload is not None and self.use_cache and cacheable:
            self.store.cache_put(key, payload)
        return payload

    def get_text(self, url: str, params: dict | None = None) -> str | None:
        host = url.split("/")[2] if "://" in url else "_default"
        return self._request_with_retry("GET", url, params, host, want_json=False)

    def _request_with_retry(
        self, method: str, url: str, params: dict | None, host: str, want_json: bool
    ) -> Any:
        last_err: str | None = None
        for attempt in range(self.max_retries):
            self.limiter.wait(host)
            try:
                self.n_requests += 1
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
                        self.n_failures += 1
                        return None
                    last_err = "403 Forbidden"
                elif resp.status_code in TRANSIENT_STATUS:
                    last_err = f"HTTP {resp.status_code}"
                else:
                    log.warning("HTTP %s from %s", resp.status_code, url)
                    self.n_failures += 1
                    return None
            except requests.RequestException as e:
                last_err = str(e)

            backoff = (2**attempt) + random.uniform(0, 0.5)
            log.debug("Retry %d/%d for %s after %s (sleep %.1fs)",
                      attempt + 1, self.max_retries, url, last_err, backoff)
            time.sleep(backoff)

        log.warning("Giving up on %s after %d attempts: %s",
                    url, self.max_retries, last_err)
        self.n_failures += 1
        return None

    def stats(self) -> dict:
        return {
            "requests": self.n_requests,
            "cache_hits": self.n_cache_hits,
            "failures": self.n_failures,
        }
