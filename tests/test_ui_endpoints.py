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

    # 8. /api/vcf endpoint
    print("\n[/api/vcf]")
    sample_vcf_str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr8\t79616822\trs1\tC\tT\t99\tPASS\tGENE=STMN2;CLNSIG=Pathogenic;SpliceAI=T|STMN2|0.85|0.0|0.0|0.0|0|0|0|0"
    res = client.post("/api/vcf", json={"vcf_text": sample_vcf_str})
    check("/api/vcf returns 200 OK", res.status_code == 200)
    data = res.get_json()
    check("/api/vcf returns candidates list", "candidates" in data and len(data["candidates"]) == 1)
    check("/api/vcf candidate is STMN2", data["candidates"][0]["gene_symbol"] == "STMN2")

    # 9. /api/repurpose endpoint
    print("\n[/api/repurpose]")
    res = client.post("/api/repurpose", json={"target": "STMN2", "direction": "activator", "bbb_only": True})
    check("/api/repurpose returns 200 OK", res.status_code == 200)
    data = res.get_json()
    check("/api/repurpose returns candidates", "candidates" in data and len(data["candidates"]) > 0)
    check("/api/repurpose candidates are BBB penetrant", data["candidates"][0]["bbb_plausible"] is True)

    # 10. /api/visualize/splicing endpoint
    print("\n[/api/visualize/splicing]")
    res = client.get("/api/visualize/splicing?gene=STMN2&aso=CTGAATACTCACCTGAAA")
    check("/api/visualize/splicing returns 200 OK", res.status_code == 200)
    check("/api/visualize/splicing mimetype is image/svg+xml", res.mimetype == "image/svg+xml")
    check("/api/visualize/splicing contains SVG root", b"<svg" in res.data)

    # 11. /api/nof1/preind endpoint
    print("\n[/api/nof1/preind]")
    res = client.post("/api/nof1/preind", json={"gene": "STMN2", "patient": "PATIENT-01"})
    check("/api/nof1/preind returns 200 OK", res.status_code == 200)
    data = res.get_json()
    check("/api/nof1/preind returns markdown", "markdown" in data and "21 CFR Part 312" in data["markdown"])

    td.cleanup()
    print("\n" + "=" * 70)
    print("ALL UI ENDPOINT TESTS PASSED")
    print("=" * 70)


if __name__ == "__main__":
    run_tests()
