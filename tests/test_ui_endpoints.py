"""Tests for web interface HTTP endpoints (`synaptik/ui/server.py`)."""

import json
from pathlib import Path
from synaptik.ui.server import create_app
from synaptik.store import Store


def check(desc, cond):
    if cond:
        print(f"  PASS  {desc}")
    else:
        print(f"  FAIL  {desc}")
        raise AssertionError(desc)


def run_tests():
    print("=" * 70)
    print("SYNAPTIK UI & API endpoints test suite")
    print("=" * 70)

    # Initialize a clean temporary database
    import tempfile
    td = tempfile.TemporaryDirectory()
    db_path = str(Path(td.name) / "test_synaptik.db")
    Store(db_path).close()

    app = create_app(db_path)
    client = app.test_client()

    # 1. /api/cite endpoint
    print("\n[/api/cite]")
    res = client.get("/api/cite")
    check("/api/cite returns 200 OK", res.status_code == 200)
    data = res.get_json()
    check("/api/cite returns Francesco Nobile as author", data["author"] == "Francesco Nobile")
    check("/api/cite returns verified ORCID", "0009-0001-9477-3332" in data["orcid"])
    check("/api/cite includes bibtex and apa citations", "bibtex" in data["citations"] and "apa" in data["citations"])

    # 2. /api/nof1 endpoint
    print("\n[/api/nof1]")
    res = client.post("/api/nof1", json={
        "gene": "STMN2",
        "patient": "PATIENT-UI-01",
        "disease": "als",
        "variant": "cryptic_donor_activation",
        "nfl": 72.0,
    })
    check("/api/nof1 returns 200 OK", res.status_code == 200)
    data = res.get_json()
    check("/api/nof1 returns dossier and markdown", "dossier" in data and "markdown" in data)
    d = data["dossier"]
    check("/api/nof1 dossier contains Tier 1 actionability", "Tier 1" in d["actionability_tier"])
    check("/api/nof1 includes lead ASO sequence", len(d["stage3_therapeutic_aso"]["lead_oligo_sequence"]) >= 15)
    check("/api/nof1 includes protocol SHA-256 digest", len(d["stage6_preclinical_protocol"]["protocol_sha256"]) == 64)

    # 3. /api/aso endpoint
    print("\n[/api/aso]")
    res = client.post("/api/aso", json={
        "sequence": "TTCTGTAACAGTTTTTGTTTTTTCAGGTGAGTATTCAGGT",
        "gene": "STMN2",
        "mode": "splice_switch",
    })
    check("/api/aso returns 200 OK", res.status_code == 200)
    data = res.get_json()
    check("/api/aso returns candidates", len(data["candidates"]) > 0)
    check("/api/aso lead candidate has vendor IDT format", "vendor_idt_format" in data["candidates"][0])

    # 4. /api/translate endpoint
    print("\n[/api/translate]")
    res = client.post("/api/translate", json={
        "genes": ["STMN2", "UNC13A"],
        "disease": "als",
    })
    check("/api/translate returns 200 OK", res.status_code == 200)
    data = res.get_json()
    check("/api/translate returns gene calls", len(data["gene_calls"]) >= 2)

    # 5. /api/simulate endpoint
    print("\n[/api/simulate]")
    res = client.post("/api/simulate", json={
        "target": "STMN2",
        "disease": "als",
        "intervention": "aso_steric_block",
        "alpha": 0.85,
    })
    check("/api/simulate returns 200 OK", res.status_code == 200)
    data = res.get_json()
    check("/api/simulate returns rescue efficiency", data["rescue_efficiency_pct"] > 15.0)

    # 6. /api/biomarker endpoint
    print("\n[/api/biomarker]")
    res = client.post("/api/biomarker", json={
        "target": "STMN2",
        "disease": "als",
        "rescue_pct": 75.0,
        "baseline_nfl": 65.0,
    })
    check("/api/biomarker returns 200 OK", res.status_code == 200)
    data = res.get_json()
    check("/api/biomarker returns stratification", "stratification" in data)

    # 7. /api/protocol endpoint
    print("\n[/api/protocol]")
    res = client.post("/api/protocol", json={
        "target": "TARDBP",
        "disease": "als",
        "pairs": 8,
        "rho": 0.75,
    })
    check("/api/protocol returns 200 OK", res.status_code == 200)
    data = res.get_json()
    check("/api/protocol returns SHA-256 fingerprint", len(data["protocol_fingerprint_sha256"]) == 64)

    td.cleanup()
    print("\n" + "=" * 70)
    print("ALL UI ENDPOINT TESTS PASSED")
    print("=" * 70)


if __name__ == "__main__":
    run_tests()
