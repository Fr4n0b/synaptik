"""
Research-group identity resolution.

This module exists because of the single largest failure mode in literature
meta-analysis: counting one lab's five papers as five independent
confirmations. If a group publishes a finding in 2019, extends it in 2021, and
reviews it in 2023, a naive count reports "3 studies agree" when the true
independent-replication count is one.

We resolve group identity by clustering on (a) shared senior authors and
(b) shared institutional affiliation, then expose the cluster ID so all
downstream statistics can be computed per-lab rather than per-paper.

This is heuristic and will make mistakes -- two different Smith labs at
different institutions, or one PI who moved institutions. We therefore expose
the clustering for inspection and let the user override it via a curation file.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

from ..schema import SourceMetadata

# Institution name variants -> canonical form. Affiliation strings are
# notoriously inconsistent; this normalizes the highest-frequency offenders in
# the neuroscience/stem-cell literature.
INSTITUTION_CANONICAL = {
    "massachusetts general hospital": "MGH",
    "mass general": "MGH",
    "mgh": "MGH",
    "brigham and womens hospital": "BWH",
    "harvard medical school": "HMS",
    "harvard university": "Harvard",
    "broad institute": "Broad",
    "stanley center": "Broad",
    "massachusetts institute of technology": "MIT",
    "mit": "MIT",
    "whitehead institute": "Whitehead",
    "stanford university": "Stanford",
    "university of california san francisco": "UCSF",
    "ucsf": "UCSF",
    "university of california san diego": "UCSD",
    "ucsd": "UCSD",
    "columbia university": "Columbia",
    "new york stem cell foundation": "NYSCF",
    "johns hopkins": "JHU",
    "university of pennsylvania": "UPenn",
    "washington university": "WashU",
    "northwestern university": "Northwestern",
    "university of michigan": "UMich",
    "yale university": "Yale",
    "national institutes of health": "NIH",
    "cedars-sinai": "Cedars-Sinai",
    "salk institute": "Salk",
    "scripps research": "Scripps",
    "university college london": "UCL",
    "university of cambridge": "Cambridge",
    "university of oxford": "Oxford",
    "francis crick institute": "Crick",
    "karolinska institutet": "Karolinska",
    "max planck": "MaxPlanck",
    "riken": "RIKEN",
    "university of tokyo": "UTokyo",
    "kyoto university": "Kyoto",
}

_DEPT_NOISE = re.compile(
    r"\b(department|dept|division|institute|center|centre|program|programme|"
    r"school|faculty|laboratory|lab|group|unit|of|for|the|and|graduate)\b"
)


def canonical_institution(affiliation: str) -> str | None:
    """Map a raw affiliation string to a canonical institution token."""
    if not affiliation:
        return None
    a = affiliation.lower()
    a = re.sub(r"[^\w\s-]", " ", a)
    a = re.sub(r"\s+", " ", a).strip()
    for variant, canon in INSTITUTION_CANONICAL.items():
        if variant in a:
            return canon
    # Fall back: strip departmental noise, take the longest remaining phrase
    cleaned = _DEPT_NOISE.sub(" ", a)
    cleaned = re.sub(r"\s+", " ", cleaned).strip()
    tokens = [t for t in cleaned.split() if len(t) > 3]
    if not tokens:
        return None
    return "".join(t.capitalize() for t in tokens[:2])[:40]


def normalize_author(name: str) -> str:
    """Reduce an author name to 'surname_firstinitial'.

    Full-name matching fails on 'J Smith' vs 'John Smith' vs 'Smith J'.
    Surname + first initial is the standard compromise; it over-merges common
    names, which we accept because affiliation is required to co-occur.
    """
    if not name:
        return ""

    # "Smith, John A." -- the comma is an explicit surname-first marker and is
    # unambiguous, so handle it before any heuristic.
    if "," in name:
        surname_part, _, given_part = name.partition(",")
        surname_tokens = [t for t in re.split(r"[\s.]+", surname_part) if t]
        given_tokens = [t for t in re.split(r"[\s.]+", given_part) if t]
        if surname_tokens and given_tokens:
            return f"{surname_tokens[0].lower()}_{given_tokens[0][0].lower()}"

    n = re.sub(r"[^\w\s.-]", " ", name).strip()
    parts = [p for p in re.split(r"[\s.]+", n) if p]
    if not parts:
        return ""
    if len(parts) == 1:
        return parts[0].lower() + "_"
    # "Smith J" / "Smith JA" -- trailing initials mean surname came first
    if len(parts[-1]) <= 2 and parts[-1].isupper():
        surname, initial = parts[0], parts[-1][0]
    else:
        surname, initial = parts[-1], parts[0][0]
    return f"{surname.lower()}_{initial.lower()}"


class LabResolver:
    """Assigns a stable lab_id to each source via senior-author clustering."""

    def __init__(self, curation_file: str | Path | None = None):
        self.curation: dict[str, str] = {}
        if curation_file and Path(curation_file).exists():
            self.curation = json.loads(Path(curation_file).read_text())
        # union-find over (author, institution) signatures
        self._parent: dict[str, str] = {}
        self._signature_to_sources: dict[str, list[str]] = defaultdict(list)

    def _find(self, x: str) -> str:
        self._parent.setdefault(x, x)
        while self._parent[x] != x:
            self._parent[x] = self._parent[self._parent[x]]
            x = self._parent[x]
        return x

    def _union(self, a: str, b: str) -> None:
        ra, rb = self._find(a), self._find(b)
        if ra != rb:
            self._parent[rb] = ra

    @staticmethod
    def senior_author_signatures(meta: SourceMetadata) -> list[str]:
        """Senior authors are last and second-to-last; also the corresponding author.

        In biology the last author is the PI. We include second-to-last because
        co-senior authorship is now common, and the corresponding author
        because it is explicit when available.
        """
        sigs: list[str] = []
        authors = [a for a in meta.authors if a.strip()]
        candidates: list[str] = []
        if authors:
            candidates.append(authors[-1])
            if len(authors) >= 3:
                candidates.append(authors[-2])
        if meta.corresponding_author:
            candidates.append(meta.corresponding_author)

        institutions = [
            canonical_institution(a) for a in meta.affiliations
        ]
        institutions = [i for i in institutions if i] or [None]

        for c in candidates:
            an = normalize_author(c)
            if not an:
                continue
            for inst in institutions:
                sigs.append(f"{an}@{inst or 'unknown'}")
        return sigs

    def fit(self, sources: list[SourceMetadata]) -> None:
        """Build the clustering. Call once over the corpus before assigning."""
        for meta in sources:
            sigs = self.senior_author_signatures(meta)
            if not sigs:
                continue
            # All senior-author signatures on one paper belong to one lab-cluster
            base = sigs[0]
            for s in sigs[1:]:
                self._union(base, s)
            for s in sigs:
                self._signature_to_sources[s].append(meta.source_id)

    def assign(self, meta: SourceMetadata) -> str:
        """Return a stable lab_id for this source."""
        if meta.source_id in self.curation:
            return self.curation[meta.source_id]
        sigs = self.senior_author_signatures(meta)
        if not sigs:
            # No author info: treat as its own singleton lab, never merged
            return "lab_unknown_" + hashlib.sha1(
                meta.source_id.encode()
            ).hexdigest()[:8]
        root = self._find(sigs[0])
        return "lab_" + hashlib.sha1(root.encode()).hexdigest()[:10]

    def cluster_report(self) -> dict[str, list[str]]:
        """lab_id -> source_ids, for human inspection of the clustering."""
        out: dict[str, list[str]] = defaultdict(list)
        for sig, srcs in self._signature_to_sources.items():
            lab = "lab_" + hashlib.sha1(self._find(sig).encode()).hexdigest()[:10]
            for s in srcs:
                if s not in out[lab]:
                    out[lab].append(s)
        return dict(out)


def independent_lab_count(sources: list[SourceMetadata]) -> int:
    """Number of distinct labs represented. The honest denominator."""
    return len({s.lab_id for s in sources if s.lab_id})
