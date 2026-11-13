"""
Safety-claim regression suite.

Every check here encodes a way SYNAPTIK could mislead someone into treating an
unvalidated sequence as a therapeutic. They are not style rules. They exist
because each one was a real defect found in v0.8.0:

  * The off-target screen searched ten hard-coded transcript fragments
    totalling 861 bp (~0.0007% of the human transcriptome) and reported
    "0-mismatch: 0, 1-mismatch: 0". Measured empirically it flagged 0 of 3000
    random 20-mers, so every candidate passed and the report looked like a
    clean genome-wide result.
  * "Predicted CSF Half-Life: ~135 days" was a constant, identical for every
    sequence, presented as a per-candidate prediction.
  * The N-of-1 output was titled "Translational Preclinical IND Dossier",
    asserted that in vitro target engagement had been demonstrated and that a
    "whole-transcriptome interrogation" had been run, and stated that the
    sequence "lacks unmethylated CpG motifs" without ever computing CpG
    content.
  * The Pre-IND output declared "21 CFR Part 312 compliant" and named a real
    person as "Sponsor / Lead Investigator".

The common failure is a confident claim about something that was never
computed. A reader cannot distinguish a screen that found nothing from a
screen that could not find anything, so the tool must never print the first
when it did the second.

Run: python3 -m tests.test_safety_claims
"""

from __future__ import annotations

import random
import sys

from synaptik.drug.aso import (
    ASODesigner, ASOMode, count_cpg_motifs, screen_off_targets,
)
from synaptik.drug.specificity import (
    ScreenStatus, SpecificityScreener, reverse_complement, unscreened,
)
from synaptik.nof1.pipeline import (
    Nof1PatientCase, generate_fda_pre_ind_package, nof1_dossier_to_markdown,
    run_nof1_pipeline,
)

FAILURES: list[str] = []

TARGET = ("GGCAGCAGTAAGCAGGTAAGTCCCTTTGCAATGCTGGATTTTGTGGAACTTT"
          "ATCAGCCAGTAGCTGAAGAA")


def check(cond: bool, msg: str) -> None:
    print(f"  {'PASS' if cond else 'FAIL'}  {msg}")
    if not cond:
        FAILURES.append(msg)


# --------------------------------------------------------------------------


def test_screen_is_exhaustive() -> None:
    print("\n[specificity: the search must not miss hits]")

    def brute(oligo: str, seq: str, m: int):
        L = len(oligo)
        out = []
        for probe, strand in ((reverse_complement(oligo), "+"), (oligo, "-")):
            for i in range(len(seq) - L + 1):
                d = sum(1 for a, b in zip(probe, seq[i:i + L]) if a != b)
                if d <= m:
                    out.append((strand, i, d))
        return sorted(out)

    rng = random.Random(42)
    agree = True
    for _ in range(40):
        m = rng.choice([0, 1, 2, 3])
        oligo = "".join(rng.choice("ACGT") for _ in range(rng.choice([16, 18, 20, 22])))
        seq = "".join(rng.choice("ACGT") for _ in range(3000))
        planted = list(reverse_complement(oligo))
        for pos in rng.sample(range(len(planted)), m):
            planted[pos] = rng.choice("ACGT")
        at = rng.randrange(0, 2500)
        seq = seq[:at] + "".join(planted) + seq[at + len(planted):]
        sc = SpecificityScreener(max_mismatches=m, max_hits_per_oligo=10 ** 6)
        got = sorted((h.strand, h.position, h.n_mismatches)
                     for h in sc.screen_oligo(oligo, [("tx", seq)]).hits)
        if got != brute(oligo, seq, m):
            agree = False
            break
    check(agree, "seed-and-extend finds exactly what brute force finds "
                 "(40 randomised trials, 0-3 mismatches)")

    # Both strands, because a transcript may be deposited either way.
    oligo = "ACGTACGTACGTACGTACGT"
    fwd = reverse_complement(oligo)
    sc = SpecificityScreener(max_mismatches=0)
    check(bool(sc.screen_oligo(oligo, [("t", "TTTT" + fwd + "TTTT")]).hits),
          "a perfect target match is found on the forward strand")
    check(bool(sc.screen_oligo(oligo, [("t", "TTTT" + oligo + "TTTT")]).hits),
          "a perfect target match is found on the reverse strand")


def test_no_fabricated_zeros() -> None:
    print("\n[specificity: silence is not a clean result]")
    r = unscreened("ACGTACGTACGTACGTACGT")
    check(r.status is ScreenStatus.NOT_SCREENED,
          "an unscreened oligo is marked NOT_SCREENED")
    check("NOT SCREENED" in r.verdict(),
          "the verdict says NOT SCREENED rather than reporting zero hits")
    check("no specificity evidence" in r.verdict().lower(),
          "the verdict states that no specificity evidence exists")

    res = ASODesigner().design(TARGET, target_gene="STMN2", top_k=2)
    lead = res.candidates[0]
    check(not lead.specificity_screened,
          "a candidate designed without a reference is not marked screened")
    check(lead.specificity.status is ScreenStatus.NOT_SCREENED,
          "its specificity status is NOT_SCREENED")
    check(any("NO OFF-TARGET SCREEN WAS RUN" in c for c in res.caveats),
          "the result carries a prominent unscreened caveat")
    check(any("not screened" in w.lower() for w in lead.penalties_and_warnings),
          "the candidate itself records that specificity was not screened")


def test_panel_is_not_a_screen() -> None:
    print("\n[specificity: the legacy panel must not be mistaken for a screen]")
    rng = random.Random(7)
    flagged = sum(
        1 for _ in range(500)
        if any(screen_off_targets(
            "".join(rng.choice("ACGT") for _ in range(20)))[:2]))
    check(flagged == 0,
          f"the 861 bp panel flags nothing on random 20-mers ({flagged}/500) "
          "- documented as a smoke test, never cited as a safety result")
    check("NOT a specificity screen" in (screen_off_targets.__doc__ or ""),
          "its docstring states plainly that it is not a specificity screen")


def test_real_screen_catches_what_panel_missed() -> None:
    print("\n[specificity: a real reference changes the answer]")
    rng = random.Random(3)
    res = ASODesigner().design(TARGET, target_gene="STMN2", top_k=1)
    oligo = res.candidates[0].oligo_sequence_5_to_3

    # The panel sees nothing; a reference containing a perfect off-target does.
    check(screen_off_targets(oligo)[:2] == (0, 0),
          "the panel reports zero hits for this oligo")
    decoys = [(f"TX{i}", "".join(rng.choice("ACGT") for _ in range(2000)))
              for i in range(50)]
    danger = ("".join(rng.choice("ACGT") for _ in range(500))
              + reverse_complement(oligo)
              + "".join(rng.choice("ACGT") for _ in range(500)))
    hits = SpecificityScreener(max_mismatches=2).screen_oligo(
        oligo, decoys + [("DANGER", danger)]).hits
    check(any(h.n_mismatches == 0 for h in hits),
          "a real screen finds the perfect off-target the panel missed")


def test_half_life_not_presented_as_prediction() -> None:
    print("\n[half-life: a constant must not be labelled a prediction]")
    rng = random.Random(5)
    seen = set()
    for _ in range(5):
        seq = "".join(rng.choice("ACGT") for _ in range(72))
        r = ASODesigner().design(seq, top_k=1)
        if r.candidates:
            seen.add(r.candidates[0].chemistry_class_csf_half_life_days)
    check(len(seen) == 1,
          "the value is identical across unrelated sequences, confirming it "
          "is a class constant")
    check(not hasattr(ASODesigner().design(TARGET, top_k=1).candidates[0],
                      "predicted_cns_half_life_days"),
          "it is no longer exposed as 'predicted_cns_half_life_days'")


def test_cpg_is_computed() -> None:
    print("\n[immunogenicity: CpG content must be counted, not asserted]")
    check(count_cpg_motifs("ACGTACGTACGT")[0] == 3,
          "CpG dinucleotides are counted correctly")
    check(count_cpg_motifs("AAATTTAAA")[0] == 0,
          "a CpG-free sequence counts zero")
    check(count_cpg_motifs("acgtacgt")[0] == 2, "counting is case-insensitive")
    check(count_cpg_motifs("ACGUACGU")[0] == 2, "RNA input is handled")


# --------------------------------------------------------------------------


def _documents():
    case = Nof1PatientCase(
        patient_id="TEST-001", target_gene="STMN2",
        variant_id="chr8:79611214:G>A", target_sequence=TARGET,
    )
    d = run_nof1_pipeline(case)
    return d, nof1_dossier_to_markdown(d), generate_fda_pre_ind_package(d)


def test_documents_are_not_regulatory_submissions() -> None:
    print("\n[documents: must not read as completed regulatory packages]")
    _, md, preind = _documents()

    for name, doc in (("dossier", md), ("pre-IND outline", preind)):
        check("COMPUTATIONAL DRAFT - NOT A REGULATORY SUBMISSION" in doc,
              f"{name} opens with the not-a-submission banner")
        check("No experiment described here has been performed" in doc,
              f"{name} states that no experiment was performed")
        check("Do not submit this document to a regulator" in doc,
              f"{name} warns against submitting it to a regulator")

    check("IND Dossier" not in md,
          "the dossier is not titled as an IND dossier")
    check("21 CFR Part 312 compliant" not in preind,
          "the outline does not claim 21 CFR Part 312 compliance")
    check("**Sponsor / Lead Investigator**" not in preind,
          "the outline does not name a regulatory sponsor")
    check("DRAFT PLANNING OUTLINE" in preind,
          "the outline is titled as a draft planning document")


def test_documents_make_no_false_experimental_claims() -> None:
    print("\n[documents: no claim that an experiment happened]")
    _, md, preind = _documents()
    both = md + preind

    banned = {
        "Demonstrates complete steric exclusion":
            "asserts in vitro target engagement was shown",
        "Whole-transcriptome interrogation":
            "claims a whole-transcriptome screen was run",
        "Sequence lacks unmethylated CpG motifs":
            "asserts absence of CpG motifs without computing them",
        "Orderable IDT Vendor Sequence":
            "describes an unvalidated oligo as orderable",
        "Longitudinal In Vivo Suppression":
            "labels simulated values as in vivo",
    }
    for phrase, why in banned.items():
        check(phrase not in both, f"no text that {why}")

    check("STUDY NOT PERFORMED" in preind,
          "the unperformed in vitro study is explicitly flagged")
    check("CpG / TLR9 Motif Count" in preind,
          "CpG content is reported as a computed count")
    check("not measured" in both.lower() or "chemistry-class" in both.lower(),
          "the half-life is qualified as unmeasured / class-level")


def test_specificity_status_reaches_the_documents() -> None:
    print("\n[documents: specificity status must propagate]")
    d, md, preind = _documents()
    s3 = d.stage3_therapeutic_aso
    check("specificity_status" in s3,
          "the dossier payload carries a specificity status")
    check(s3["specificity_status"] == ScreenStatus.NOT_SCREENED.value,
          "an unscreened run reports NOT_SCREENED, not zero hits")
    check("NOT PERFORMED" in md or "NOT SCREENED" in md.upper(),
          "the dossier surfaces that no screen was performed")
    check(s3["cpg_motif_count"] >= 0, "a CpG count is present in the payload")


def test_clinical_allometry_guardrail() -> None:
    print("\n[clinical: intrathecal dose calculations must require explicit safety acknowledgment]")
    from synaptik.clinical.allometry import calculate_allometric_dose

    # Casual invocation without acknowledgment MUST raise ValueError
    blocked = False
    try:
        calculate_allometric_dose(noael_mg=30.0)
    except ValueError as e:
        blocked = "RESEARCH PLANNING ONLY" in str(e)
    check(blocked, "calling calculate_allometric_dose without acknowledgment raises ValueError")

    # When acknowledged, outputs must carry explicit research planning disclaimers
    res = calculate_allometric_dose(noael_mg=30.0, acknowledge_not_clinical=True)
    check(res.is_clinical_recommendation is False,
          "AllometryResult marks is_clinical_recommendation as False")
    check("RESEARCH PLANNING ONLY" in res.disclaimer,
          "AllometryResult carries RESEARCH PLANNING ONLY disclaimer")
    d = res.to_dict()
    check(d.get("is_clinical_recommendation") is False,
          "serialized dictionary exposes is_clinical_recommendation as False")
    check("RESEARCH PLANNING ONLY" in d.get("disclaimer", ""),
          "serialized dictionary exposes research planning disclaimer")


def test_specificity_screened_property_and_gapmer_thermo() -> None:
    print("\n[specificity & chemistry: .screened property and gapmer 2'-MOE wing count]")
    from synaptik.drug.aso import ASOChemistry, compute_thermodynamics
    from synaptik.drug.specificity import ScreenStatus, SpecificityResult, unscreened

    u = unscreened("TCACTTTCATAATGCTGG")
    check(u.screened is False, "unscreened SpecificityResult.screened must be False")
    p = SpecificityResult(status=ScreenStatus.PANEL, oligo="TCACTTTCATAATGCTGG")
    check(p.screened is False, "PANEL-only SpecificityResult.screened must be False")
    t = SpecificityResult(status=ScreenStatus.TRANSCRIPTOME, oligo="TCACTTTCATAATGCTGG")
    check(t.screened is True, "TRANSCRIPTOME SpecificityResult.screened must be True")

    # 20-mer 5-10-5 gapmer has 10 MOE wing nt, so Tm < uniform 20-mer MOE/PS
    tofersen = "CAGGATACATTTCTACAGCT"
    _, _, dg_gap, tm_gap = compute_thermodynamics(tofersen, chemistry=ASOChemistry.GAPMER_5_10_5)
    _, _, dg_moe, tm_moe = compute_thermodynamics(tofersen, chemistry=ASOChemistry.MOE_PS)
    check(tm_gap < tm_moe, f"5-10-5 gapmer Tm ({tm_gap}C) < uniform 2'-MOE Tm ({tm_moe}C)")
    check(dg_gap > dg_moe, f"5-10-5 gapmer dG37 ({dg_gap}) less negative than uniform 2'-MOE ({dg_moe})")

    from synaptik.drug.specificity import SpecificityScreener
    refs = [
        ("TX1", "AAAAATTTTTTCACTTTCATAATGCTGGGGGGGG"),
        ("TX2", "CCCCCCCCAGCTGTAGAAATGTATCCTGGGGGGG"),
    ]
    oligos = ["CCAGCATTATGAAAGTGA", "CAGGATACATTTCTACAGCT"]
    sc = SpecificityScreener(max_mismatches=2)
    batch = sc.screen_batch(oligos, refs)
    singles = [sc.screen_oligo(o, refs) for o in oligos]
    check(
        all(
            [(h.transcript_id, h.position, h.strand, h.n_mismatches) for h in b.hits]
            == [(h.transcript_id, h.position, h.strand, h.n_mismatches) for h in s.hits]
            for b, s in zip(batch, singles)
        ),
        "screen_batch matches single-oligo screen_oligo exactly",
    )


def main() -> int:
    print("=" * 70)
    print("SYNAPTIK safety-claim regression suite")
    print("=" * 70)
    test_screen_is_exhaustive()
    test_no_fabricated_zeros()
    test_panel_is_not_a_screen()
    test_real_screen_catches_what_panel_missed()
    test_half_life_not_presented_as_prediction()
    test_cpg_is_computed()
    test_documents_are_not_regulatory_submissions()
    test_documents_make_no_false_experimental_claims()
    test_specificity_status_reaches_the_documents()
    test_clinical_allometry_guardrail()
    test_specificity_screened_property_and_gapmer_thermo()
    print("\n" + "=" * 70)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S):")
        for f in FAILURES:
            print("  -", f)
        return 1
    print("All safety-claim checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
