"""Tests for Patient Genomic VCF Parsing and Splice Prioritization (`synaptik vcf`)."""

import json
from pathlib import Path
from synaptik.sequence.vcf import (
    KNOWN_NEURO_TARGETS,
    GenomicVariant,
    PatientVCFParser,
    variant_to_dict,
)
from synaptik.cli.commands.vcf import cmd_vcf


def check(desc, cond):
    if cond:
        print(f"  PASS  {desc}")
    else:
        print(f"  FAIL  {desc}")
        raise AssertionError(desc)


class MockArgs:
    def __init__(self, vcf_file, target=None, min_splice_score=0.20, json=False):
        self.vcf_file = vcf_file
        self.target = target
        self.min_splice_score = min_splice_score
        self.json = json


SAMPLE_VCF = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr8\t79616822\trs183204964\tC\tT\t99.0\tPASS\tGENE=STMN2;CLNSIG=Pathogenic;gnomAD_AF=0.00002;CADD_PHRED=26.4;SpliceAI=T|STMN2|0.88|0.01|0.02|0.00|12|-4|0|0
chr19\t17642414\trs12608932\tA\tC\t99.0\tPASS\tGENE=UNC13A;CLNSIG=Pathogenic;gnomAD_AF=0.00015;CADD_PHRED=24.1;SpliceAI=C|UNC13A|0.12|0.00|0.82|0.01|-5|0|18|0
chr5\t70951950\trs121909192\tC\tT\t99.0\tPASS\tGENE=SMN2;CLNSIG=Pathogenic;gnomAD_AF=0.00001;CADD_PHRED=28.2;SpliceAI=T|SMN2|0.05|0.92|0.01|0.00|0|-22|0|0
chr2\t10000000\trs99999999\tG\tA\t45.0\tPASS\tGENE=NOT_TARGET;CLNSIG=Benign;gnomAD_AF=0.45;CADD_PHRED=2.1;SpliceAI=A|NOT_TARGET|0.01|0.01|0.00|0.00|0|0|0|0
"""


def run_tests():
    print("=" * 70)
    print("SYNAPTIK Patient Genomic VCF Splice Prioritization test suite")
    print("=" * 70)

    # ------------------------------------------------------------------
    # 1. Target loci catalog validation
    # ------------------------------------------------------------------
    print("\n[target loci catalog]")
    check("STMN2 in target loci", "STMN2" in KNOWN_NEURO_TARGETS)
    check("UNC13A in target loci", "UNC13A" in KNOWN_NEURO_TARGETS)
    check("SMN2 in target loci", "SMN2" in KNOWN_NEURO_TARGETS)
    check("STMN2 has cryptic junction defined", "cryptic_junction" in KNOWN_NEURO_TARGETS["STMN2"])

    # ------------------------------------------------------------------
    # 2. VCF string content parsing
    # ------------------------------------------------------------------
    print("\n[vcf string content parsing]")
    parser = PatientVCFParser()
    variants = parser.parse_vcf_content(SAMPLE_VCF)
    check("parsed 4 variants from sample VCF", len(variants) == 4)

    v1 = variants[0]
    check("variant 1 is on chr8", v1.chrom == "chr8")
    check("variant 1 pos is 79616822", v1.pos == 79616822)
    check("variant 1 ID is rs183204964", v1.id == "rs183204964")
    check("variant 1 gene is STMN2", v1.gene_symbol == "STMN2")
    check("variant 1 ClinVar is Pathogenic", "Pathogenic" in v1.clinical_significance)
    check("variant 1 gnomAD AF parsed correctly", v1.gnomad_af == 0.00002)
    check("variant 1 CADD parsed correctly", v1.cadd_phred == 26.4)
    check("variant 1 SpliceAI delta is 0.88", v1.spliceai_max_delta == 0.88)
    check("variant 1 is cryptic splice candidate", v1.is_cryptic_splice_candidate)

    # ------------------------------------------------------------------
    # 3. Actionable variant prioritization
    # ------------------------------------------------------------------
    print("\n[actionable variant prioritization]")
    prioritized = parser.prioritize_actionable_variants(variants, min_splice_score=0.20)
    check("prioritized filters out non-target benign variant", len(prioritized) == 3)
    check("top candidate is STMN2 or SMN2", prioritized[0].gene_symbol in ["STMN2", "SMN2"])
    check("top candidate priority score > 80", prioritized[0].priority_score >= 80.0)

    # Filtering by target gene
    parser_stmn2 = PatientVCFParser(target_genes=["STMN2"])
    p_stmn2 = parser_stmn2.prioritize_actionable_variants(variants)
    check("filtering by target STMN2 returns 1 candidate", len(p_stmn2) == 1)
    check("filtered candidate is STMN2", p_stmn2[0].gene_symbol == "STMN2")

    # ------------------------------------------------------------------
    # 4. Parsing from file on disk
    # ------------------------------------------------------------------
    print("\n[file on disk ingestion]")
    example_vcf = Path(__file__).parent.parent / "examples" / "patient_stmn2_wgs.vcf"
    check("examples/patient_stmn2_wgs.vcf exists", example_vcf.exists())
    file_vars = parser.parse_file(example_vcf)
    check("parsed variants from example file", len(file_vars) == 4)

    # ------------------------------------------------------------------
    # 5. Serialization to dict
    # ------------------------------------------------------------------
    print("\n[serialization to dict]")
    d = variant_to_dict(v1)
    check("dict contains variant_id", "variant_id" in d)
    check("dict contains therapeutic_modality", "therapeutic_modality" in d)
    check("dict contains priority_score", "priority_score" in d)

    # ------------------------------------------------------------------
    # 6. CLI smoke tests
    # ------------------------------------------------------------------
    print("\n[cli cmd_vcf execution]")
    ret_default = cmd_vcf(MockArgs(str(example_vcf)))
    check("cmd_vcf default returns 0", ret_default == 0)

    ret_json = cmd_vcf(MockArgs(str(example_vcf), json=True))
    check("cmd_vcf --json returns 0", ret_json == 0)

    ret_filtered = cmd_vcf(MockArgs(str(example_vcf), target="UNC13A"))
    check("cmd_vcf --target UNC13A returns 0", ret_filtered == 0)

    print("\n" + "=" * 70)
    print("ALL VCF PRIORITIZATION TESTS PASSED")
    print("=" * 70)


if __name__ == "__main__":
    run_tests()
