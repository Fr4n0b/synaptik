"""`synaptik cite`: display citation metadata and export BibTeX / APA formats."""

from __future__ import annotations

import json
from ..main import Fmt, header

AUTHOR = "Francesco Nobile"
ORCID = "0009-0001-9477-3332"
ORCID_URL = "https://orcid.org/0009-0001-9477-3332"
AFFILIATION = "Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School; University of Pavia"
REPO_URL = "https://github.com/Fr4n0b/synaptik"
YEAR = "2026"


def _resolve_version(version: str | None) -> str:
    if version is not None:
        return version
    from ... import __version__
    return __version__


def get_bibtex(version: str | None = None) -> str:
    v = _resolve_version(version)
    return f"""@software{{nobile{YEAR}synaptik,
  author       = {{Nobile, Francesco}},
  title        = {{{{SYNAPTIK: Evidence-structure analysis for iPSC-neuron research}}}},
  year         = {{{YEAR}}},
  version      = {{{v}}},
  url          = {{{REPO_URL}}},
  note         = {{ORCID: {ORCID_URL}}},
  organization = {{Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School}}
}}"""


def get_apa(version: str | None = None) -> str:
    v = _resolve_version(version)
    return (
        f"Nobile, F. ({YEAR}). SYNAPTIK: Evidence-structure analysis for iPSC-neuron research "
        f"(Version {v}). Massachusetts General Hospital / Harvard Medical School; University of Pavia. "
        f"{REPO_URL}. ORCID: {ORCID_URL}."
    )


def get_nature(version: str | None = None) -> str:
    v = _resolve_version(version)
    return (
        f"Nobile F. SYNAPTIK: Evidence-structure analysis for iPSC-neuron research. "
        f"v{v} ({YEAR}). {REPO_URL} (ORCID: {ORCID})."
    )


def cmd_cite(args) -> int:
    from ... import __version__ as version

    if getattr(args, "bibtex", False):
        print(get_bibtex(version))
        return 0

    if getattr(args, "apa", False):
        print(get_apa(version))
        return 0

    if getattr(args, "json", False):
        data = {
            "author": AUTHOR,
            "orcid": ORCID,
            "orcid_url": ORCID_URL,
            "affiliation": AFFILIATION,
            "title": "SYNAPTIK: Evidence-structure analysis for iPSC-neuron research",
            "version": version,
            "year": YEAR,
            "repository": REPO_URL,
            "citations": {
                "bibtex": get_bibtex(version),
                "apa": get_apa(version),
                "nature": get_nature(version),
            },
        }
        print(json.dumps(data, indent=2))
        return 0

    header("How to Cite SYNAPTIK")
    print("If you use SYNAPTIK in published research, grant applications, or preclinical studies,")
    print("please cite Francesco Nobile and this tool:\n")

    print(Fmt.bold("APA:"))
    print(f"  {get_apa(version)}\n")

    print(Fmt.bold("Nature / Vancouver:"))
    print(f"  {get_nature(version)}\n")

    print(Fmt.bold("BibTeX (for LaTeX / Overleaf):"))
    for line in get_bibtex(version).splitlines():
        print(f"  {line}")

    print("\n" + "-" * 78)
    print(f"Author ORCID Record: {Fmt.green(ORCID_URL)}")
    print(f"GitHub Repository  : {REPO_URL}")
    print("Tip: Run 'synaptik cite --bibtex >> references.bib' to append directly.")
    return 0


def register(sub, common) -> None:
    s = sub.add_parser(
        "cite", parents=[common],
        help="display citation instructions and export BibTeX / APA formats with author ORCID",
    )
    s.add_argument("--bibtex", action="store_true", help="output only the BibTeX entry")
    s.add_argument("--apa", action="store_true", help="output only the APA citation string")
    s.add_argument("--json", action="store_true", help="output citation metadata as JSON")
    s.set_defaults(func=cmd_cite)
