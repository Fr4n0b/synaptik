"""Tests for scientific citation and author recognition (`synaptik cite`)."""

import json
from pathlib import Path
from synaptik.cli.commands.cite import (
    AUTHOR,
    ORCID,
    ORCID_URL,
    AFFILIATION,
    REPO_URL,
    get_bibtex,
    get_apa,
    get_nature,
    cmd_cite,
)


def check(desc, cond):
    if cond:
        print(f"  PASS  {desc}")
    else:
        print(f"  FAIL  {desc}")
        raise AssertionError(desc)


class MockArgs:
    def __init__(self, bibtex=False, apa=False, json=False):
        self.bibtex = bibtex
        self.apa = apa
        self.json = json


def run_tests():
    print("=" * 70)
    print("SYNAPTIK author recognition & citation test suite")
    print("=" * 70)

    # ------------------------------------------------------------------
    # 1. Author and ORCID metadata constants
    # ------------------------------------------------------------------
    print("\n[author metadata constants]")
    check("author is Francesco Nobile", AUTHOR == "Francesco Nobile")
    check("ORCID matches verified identifier", ORCID == "0009-0001-9477-3332")
    check("ORCID URL is HTTPS and valid", ORCID_URL == "https://orcid.org/0009-0001-9477-3332")
    check("affiliation includes Harvard Medical School", "Harvard Medical School" in AFFILIATION)
    check("affiliation includes University of Pavia", "University of Pavia" in AFFILIATION)
    check("repository URL points to Fr4n0b/synaptik", REPO_URL == "https://github.com/Fr4n0b/synaptik")

    # ------------------------------------------------------------------
    # 2. BibTeX format generation
    # ------------------------------------------------------------------
    print("\n[bibtex citation export]")
    bib = get_bibtex("0.6.1")
    check("bibtex contains @software entry", bib.startswith("@software{nobile"))
    check("bibtex contains author Nobile, Francesco", "author       = {Nobile, Francesco}" in bib)
    check("bibtex contains SYNAPTIK title", "SYNAPTIK" in bib)
    check("bibtex contains version 0.6.1", "version      = {0.6.1}" in bib)
    check("bibtex contains ORCID in note field", f"ORCID: {ORCID_URL}" in bib)
    check("bibtex contains repository URL", f"url          = {{{REPO_URL}}}" in bib)
    check("bibtex curly braces balanced", bib.count("{") == bib.count("}"))

    # ------------------------------------------------------------------
    # 3. APA and Nature citation generation
    # ------------------------------------------------------------------
    print("\n[apa and nature citation export]")
    apa = get_apa("0.6.1")
    check("apa format begins with Nobile, F.", apa.startswith("Nobile, F."))
    check("apa includes version 0.6.1", "(Version 0.6.1)" in apa)
    check("apa includes author ORCID", f"ORCID: {ORCID_URL}" in apa)
    check("apa includes repository URL", REPO_URL in apa)

    nature = get_nature("0.6.1")
    check("nature format begins with Nobile F.", nature.startswith("Nobile F."))
    check("nature format includes v0.6.1", "v0.6.1" in nature)
    check("nature format includes ORCID", f"(ORCID: {ORCID})" in nature)

    # ------------------------------------------------------------------
    # 4. CITATION.cff file validation
    # ------------------------------------------------------------------
    print("\n[CITATION.cff schema validation]")
    cff_path = Path(__file__).parent.parent / "CITATION.cff"
    check("CITATION.cff file exists", cff_path.exists())
    cff_text = cff_path.read_text(encoding="utf-8")
    check("cff specifies cff-version 1.2.0", 'cff-version: 1.2.0' in cff_text)
    check("cff specifies author family-names Nobile", 'family-names: "Nobile"' in cff_text)
    check("cff specifies author given-names Francesco", 'given-names: "Francesco"' in cff_text)
    check("cff specifies verified ORCID URL", f'orcid: "{ORCID_URL}"' in cff_text)
    check("cff specifies repository-code", f'repository-code: "{REPO_URL}"' in cff_text)

    # ------------------------------------------------------------------
    # 5. CLI command execution smoke tests
    # ------------------------------------------------------------------
    print("\n[cli cmd_cite handler]")
    res_default = cmd_cite(MockArgs())
    check("default cmd_cite returns 0", res_default == 0)

    res_bib = cmd_cite(MockArgs(bibtex=True))
    check("cmd_cite --bibtex returns 0", res_bib == 0)

    res_apa = cmd_cite(MockArgs(apa=True))
    check("cmd_cite --apa returns 0", res_apa == 0)

    res_json = cmd_cite(MockArgs(json=True))
    check("cmd_cite --json returns 0", res_json == 0)

    print("\n" + "=" * 70)
    print("ALL CITATION TESTS PASSED")
    print("=" * 70)


if __name__ == "__main__":
    run_tests()
