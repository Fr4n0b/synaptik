"""Tests for Small-Molecule Target Repurposing and CNS Druggability (`synaptik repurpose`)."""

import json
from synaptik.drug.repurposing import MoADirection, cns_mpo, direction_compatible
from synaptik.drug.repurposing_catalog import (
    BENCHMARK_CNS_DRUGS,
    repurposed_compound_to_dict,
    screen_repurposing_candidates,
)
from synaptik.cli.commands.repurpose import cmd_repurpose


def check(desc, cond):
    if cond:
        print(f"  PASS  {desc}")
    else:
        print(f"  FAIL  {desc}")
        raise AssertionError(desc)


class MockArgs:
    def __init__(self, target="STMN2", direction="activator", min_phase=0, bbb_only=False, json=False):
        self.target = target
        self.direction = direction
        self.min_phase = min_phase
        self.bbb_only = bbb_only
        self.json = json


def run_tests():
    print("=" * 70)
    print("SYNAPTIK Small-Molecule Drug Repurposing test suite")
    print("=" * 70)

    # ------------------------------------------------------------------
    # 1. Benchmark CNS drugs repository validation
    # ------------------------------------------------------------------
    print("\n[benchmark drugs repository]")
    check("benchmark drug database populated", len(BENCHMARK_CNS_DRUGS) >= 10)
    names = [d["name"] for d in BENCHMARK_CNS_DRUGS]
    check("Risdiplam present in benchmark repository", any("Risdiplam" in n for n in names))
    check("Riluzole present in benchmark repository", any("Riluzole" in n for n in names))
    check("Edaravone present in benchmark repository", any("Edaravone" in n for n in names))
    check("Memantine present in benchmark repository", any("Memantine" in n for n in names))

    # ------------------------------------------------------------------
    # 2. CNS MPO scoring calculation
    # ------------------------------------------------------------------
    print("\n[cns mpo physicochemical scoring]")
    # Diazepam: ideal CNS profile (MW 284, clogP 2.8, TPSA 32.7, HBD 0)
    mpo_diazepam, bbb_diazepam = cns_mpo(mw=284.7, logp=2.8, tpsa=32.7, hbd=0)
    check("diazepam CNS MPO >= 5.0", mpo_diazepam is not None and mpo_diazepam >= 5.0)
    check("diazepam BBB is plausible", bbb_diazepam is True)

    # Rapamycin: massive macrolide (MW 914, TPSA 195, HBD 3)
    mpo_rapamycin, bbb_rapamycin = cns_mpo(mw=914.2, logp=4.3, tpsa=195.0, hbd=3)
    check("rapamycin CNS MPO < 3.0", mpo_rapamycin is not None and mpo_rapamycin < 3.0)
    check("rapamycin BBB is unlikely", bbb_rapamycin is False)

    # ------------------------------------------------------------------
    # 3. Directional compatibility evaluation
    # ------------------------------------------------------------------
    print("\n[mechanism directional compatibility]")
    ok, _ = direction_compatible(MoADirection.ACTIVATOR, MoADirection.ACTIVATOR)
    check("activator matches activator", ok)

    ok_opp, _ = direction_compatible(MoADirection.ACTIVATOR, MoADirection.INHIBITOR)
    check("activator rejects inhibitor", not ok_opp)

    ok_mod, _ = direction_compatible(MoADirection.ACTIVATOR, MoADirection.MODULATOR)
    check("activator accepts modulator with caveat", ok_mod)

    # ------------------------------------------------------------------
    # 4. Small-molecule screening for neurodegenerative targets
    # ------------------------------------------------------------------
    print("\n[small-molecule target screening]")
    cands_stmn2 = screen_repurposing_candidates(target_gene="STMN2", required_direction="activator")
    check("STMN2 screening returns candidate drugs", len(cands_stmn2) > 0)
    top_stmn2 = cands_stmn2[0]
    check("top candidate has translational rank score > 50", top_stmn2.translational_rank_score > 50.0)
    check("top candidate has predicted rescue efficacy > 30%", top_stmn2.predicted_rescue_efficacy_pct >= 30.0)

    # Test BBB filter
    cands_bbb = screen_repurposing_candidates(target_gene="STMN2", bbb_only=True)
    check("bbb_only filters out non-BBB drugs", all(c.bbb_plausible for c in cands_bbb))

    # Test clinical phase filter
    cands_app = screen_repurposing_candidates(target_gene="STMN2", min_clinical_phase=4)
    check("min_phase=4 returns only approved drugs", all(c.max_clinical_phase == 4 for c in cands_app))

    # ------------------------------------------------------------------
    # 5. Serialization to dict
    # ------------------------------------------------------------------
    print("\n[serialization to dict]")
    d = repurposed_compound_to_dict(top_stmn2)
    check("dict contains name", "name" in d)
    check("dict contains chembl_id", "chembl_id" in d)
    check("dict contains cns_mpo_score", "cns_mpo_score" in d)
    check("dict contains bbb_plausible", "bbb_plausible" in d)
    check("dict contains predicted_rescue_efficacy_pct", "predicted_rescue_efficacy_pct" in d)

    # ------------------------------------------------------------------
    # 6. CLI smoke tests
    # ------------------------------------------------------------------
    print("\n[cli cmd_repurpose execution]")
    ret_default = cmd_repurpose(MockArgs(target="STMN2"))
    check("cmd_repurpose default returns 0", ret_default == 0)

    ret_json = cmd_repurpose(MockArgs(target="SMN2", json=True))
    check("cmd_repurpose --json returns 0", ret_json == 0)

    ret_bbb = cmd_repurpose(MockArgs(target="STMN2", bbb_only=True, min_phase=2))
    check("cmd_repurpose --bbb-only --min-phase 2 returns 0", ret_bbb == 0)

    print("\n" + "=" * 70)
    print("ALL DRUG REPURPOSING TESTS PASSED")
    print("=" * 70)


if __name__ == "__main__":
    run_tests()
