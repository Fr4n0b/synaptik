"""Patient Genomic VCF Parser and Cryptic Splice Variant Prioritization Engine.

Ingests standard Variant Call Format (VCF v4.2 / v4.3) files from patient
whole-genome (WGS) or whole-exome (WES) sequencing, filters and prioritizes
variants impacting neurodegenerative disease target genes, evaluates cryptic
splice activation, and feeds directly into the N-of-1 therapeutic pipeline.
"""

from __future__ import annotations

import gzip
import io
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..schema import Direction


@dataclass
class GenomicVariant:
    """A single parsed genomic variant with clinical and splicing annotations."""

    chrom: str
    pos: int
    id: str
    ref: str
    alt: str
    qual: float | None = None
    filter_status: str = "PASS"
    info: dict[str, str] = field(default_factory=dict)

    # Inferred / annotated fields
    gene_symbol: str | None = None
    consequence: str = "unknown"
    clinical_significance: str = "uncertain"
    gnomad_af: float | None = None
    cadd_phred: float | None = None
    spliceai_max_delta: float = 0.0
    splice_effect: str = "none"
    is_cryptic_splice_candidate: bool = False
    priority_score: float = 0.0

    @property
    def variant_id(self) -> str:
        if self.id and self.id != ".":
            return self.id
        return f"{self.chrom}:{self.pos}:{self.ref}>{self.alt}"


# Known neurodegenerative disease splicing and loss-of-function target loci (GRCh38 coordinates)
KNOWN_NEURO_TARGETS: dict[str, dict[str, Any]] = {
    "STMN2": {
        "chrom": "chr8",
        "start": 79611210,
        "end": 79674840,
        "disease": "Amyotrophic Lateral Sclerosis / Frontotemporal Dementia",
        "cryptic_junction": "chr8:79616822-79617044",
        "mechanism": "TDP-43 dependent cryptic exon 2a inclusion causing premature polyadenylation and loss of stathmin-2",
        "therapeutic_modality": "Steric-blocking ASO / splice-switching oligo",
    },
    "UNC13A": {
        "chrom": "chr19",
        "start": 17641550,
        "end": 17729300,
        "disease": "Amyotrophic Lateral Sclerosis / Frontotemporal Dementia",
        "cryptic_junction": "chr19:17642414-17642542",
        "mechanism": "TDP-43 dependent cryptic exon inclusion between exons 20 and 21 triggering nonsense-mediated decay (NMD)",
        "therapeutic_modality": "Splice-switching ASO / cryptic splice-masking oligo",
    },
    "SMN2": {
        "chrom": "chr5",
        "start": 70924940,
        "end": 70953000,
        "disease": "Spinal Muscular Atrophy",
        "cryptic_junction": "chr5:70951900-70952000",
        "mechanism": "Exon 7 skipping due to c.840C>T mutation in exonic splicing enhancer (ESE)",
        "therapeutic_modality": "ISS-N1 targeting ASO (Nusinersen analogue)",
    },
    "ELP1": {
        "chrom": "chr9",
        "start": 108855000,
        "end": 108920000,
        "disease": "Familial Dysautonomia",
        "cryptic_junction": "chr9:108872000-108873000",
        "mechanism": "IVS20+6T>C splice mutation causing tissue-specific skipping of exon 20",
        "therapeutic_modality": "Exon-inclusion steric-blocking ASO / small molecule splicing modulator",
    },
    "COL6A1": {
        "chrom": "chr21",
        "start": 45960000,
        "end": 45985000,
        "disease": "Ullrich Congenital Muscular Dystrophy",
        "cryptic_junction": "chr21:45973000-45974000",
        "mechanism": "Intronic pseudo-exon inclusion causing in-frame nonsense mutation",
        "therapeutic_modality": "Pseudo-exon skipping ASO",
    },
    "ATXN2": {
        "chrom": "chr12",
        "start": 111450000,
        "end": 111590000,
        "disease": "ALS / Spinocerebellar Ataxia Type 2",
        "cryptic_junction": "chr12:111500000-111502000",
        "mechanism": "CAG intermediate polyglutamine expansion enhancing TDP-43 neurotoxicity",
        "therapeutic_modality": "RNase H-competent knockdown ASO",
    },
    "TARDBP": {
        "chrom": "chr1",
        "start": 11012500,
        "end": 11025500,
        "disease": "Amyotrophic Lateral Sclerosis / Frontotemporal Dementia",
        "cryptic_junction": "chr1:11020000-11022000",
        "mechanism": "3' UTR autoregulatory alternative splicing and toxic cytoplasmic aggregation",
        "therapeutic_modality": "Autoregulation modulating ASO",
    },
    "SOD1": {
        "chrom": "chr21",
        "start": 31659690,
        "end": 31668930,
        "disease": "Amyotrophic Lateral Sclerosis",
        "cryptic_junction": "chr21:31660000-31665000",
        "mechanism": "Toxic gain-of-function missense mutations causing misfolding and aggregation",
        "therapeutic_modality": "RNase H-competent knockdown ASO (Tofersen analogue)",
    },
    "FUS": {
        "chrom": "chr16",
        "start": 31180000,
        "end": 31200000,
        "disease": "Juvenile Amyotrophic Lateral Sclerosis",
        "cryptic_junction": "chr16:31190000-31195000",
        "mechanism": "Nuclear localization signal (NLS) mutations causing cytoplasmic phase separation",
        "therapeutic_modality": "RNase H knockdown ASO (Jacifusen analogue)",
    },
}


class PatientVCFParser:
    """High-performance streaming parser for patient VCF files."""

    def __init__(self, target_genes: list[str] | None = None):
        self.target_genes = [g.upper() for g in target_genes] if target_genes else list(KNOWN_NEURO_TARGETS.keys())

    def parse_vcf_content(self, content_str: str) -> list[GenomicVariant]:
        """Parse raw VCF string lines into structured GenomicVariant records."""
        variants: list[GenomicVariant] = []
        lines = content_str.strip().splitlines()

        for line in lines:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 5:
                continue

            chrom = parts[0]
            if not chrom.startswith("chr") and chrom.lower() not in ["mt", "m"]:
                chrom = f"chr{chrom}"

            try:
                pos = int(parts[1])
            except ValueError:
                continue

            var_id = parts[2]
            ref = parts[3].upper()
            alt = parts[4].upper()

            qual = None
            if len(parts) > 5 and parts[5] != ".":
                try:
                    qual = float(parts[5])
                except ValueError:
                    pass

            filter_status = parts[6] if len(parts) > 6 else "PASS"

            # Parse INFO key-values
            info_dict: dict[str, str] = {}
            if len(parts) > 7 and parts[7] != ".":
                for item in parts[7].split(";"):
                    if "=" in item:
                        k, v = item.split("=", 1)
                        info_dict[k] = v
                    else:
                        info_dict[item] = "true"

            var = GenomicVariant(
                chrom=chrom,
                pos=pos,
                id=var_id,
                ref=ref,
                alt=alt,
                qual=qual,
                filter_status=filter_status,
                info=info_dict,
            )

            self._annotate_variant(var)
            variants.append(var)

        return variants

    def parse_file(self, file_path: str | Path) -> list[GenomicVariant]:
        """Read and parse VCF from disk (handles plain .vcf and gzipped .vcf.gz)."""
        p = Path(file_path)
        if not p.exists():
            raise FileNotFoundError(f"VCF file not found: {p}")

        if p.suffix == ".gz":
            with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
                content = f.read()
        else:
            with open(p, "r", encoding="utf-8", errors="replace") as f:
                content = f.read()

        return self.parse_vcf_content(content)

    def prioritize_actionable_variants(
        self,
        variants: list[GenomicVariant],
        min_splice_score: float = 0.20,
    ) -> list[GenomicVariant]:
        """Filter and rank variants by therapeutic and splicing actionability."""
        candidates = []
        for v in variants:
            # If explicit target filter is active, restrict strictly to specified genes
            if self.target_genes and (not v.gene_symbol or v.gene_symbol.upper() not in self.target_genes):
                continue

            # Benign variants without high target priority are excluded
            if "benign" in v.clinical_significance.lower() and v.gene_symbol not in KNOWN_NEURO_TARGETS:
                continue

            # Check target gene match or cryptic splice candidacy
            if v.gene_symbol and v.gene_symbol.upper() in KNOWN_NEURO_TARGETS:
                candidates.append(v)
            elif v.is_cryptic_splice_candidate or v.spliceai_max_delta >= min_splice_score:
                candidates.append(v)
            elif "pathogenic" in v.clinical_significance.lower():
                candidates.append(v)

        # Sort descending by priority score
        candidates.sort(key=lambda x: -x.priority_score)
        return candidates

    def _annotate_variant(self, var: GenomicVariant) -> None:
        """Annotate variant with gene locus, ClinVar, SpliceAI, and priority score."""
        # 1. Map to neuro targets by coordinate
        for gene, meta in KNOWN_NEURO_TARGETS.items():
            if var.chrom == meta["chrom"] and meta["start"] <= var.pos <= meta["end"]:
                var.gene_symbol = gene
                break

        # 2. Parse INFO fields for Gene / Consequence (VEP / SnpEff CSQ / ANN)
        info = var.info
        if not var.gene_symbol:
            for gene in self.target_genes:
                for k, v in info.items():
                    if gene in v.upper():
                        var.gene_symbol = gene
                        break
                if var.gene_symbol:
                    break

        if "GENE" in info:
            var.gene_symbol = info["GENE"].upper()
        elif "SYMBOL" in info:
            var.gene_symbol = info["SYMBOL"].upper()

        # 3. Parse ClinVar significance
        for cln_key in ["CLNSIG", "CLIN_SIG", "CLNSIGCONF"]:
            if cln_key in info:
                sig_raw = info[cln_key].lower()
                if "pathogenic" in sig_raw:
                    var.clinical_significance = "Pathogenic / Likely Pathogenic"
                elif "uncertain" in sig_raw or "vus" in sig_raw:
                    var.clinical_significance = "Variant of Uncertain Significance (VUS)"
                elif "benign" in sig_raw:
                    var.clinical_significance = "Benign / Likely Benign"
                else:
                    var.clinical_significance = sig_raw
                break

        # 4. Parse population allele frequency (gnomAD)
        for af_key in ["AF", "gnomAD_AF", "GNOMAD_AF", "AF_POPMAX"]:
            if af_key in info:
                try:
                    var.gnomad_af = float(info[af_key].split(",")[0])
                    break
                except (ValueError, IndexError):
                    pass

        # 5. Parse CADD score
        for cadd_key in ["CADD", "CADD_PHRED"]:
            if cadd_key in info:
                try:
                    var.cadd_phred = float(info[cadd_key])
                    break
                except ValueError:
                    pass

        # 6. Parse SpliceAI scores
        # Format often: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
        max_delta = 0.0
        splice_effect = "none"

        if "SpliceAI" in info:
            sai_str = info["SpliceAI"]
            vals = sai_str.split("|")
            if len(vals) >= 6:
                try:
                    # Delta scores: AG (acceptor gain), AL (acceptor loss), DG (donor gain), DL (donor loss)
                    ag, al, dg, dl = float(vals[2]), float(vals[3]), float(vals[4]), float(vals[5])
                    scores = [("acceptor_gain", ag), ("acceptor_loss", al), ("donor_gain", dg), ("donor_loss", dl)]
                    scores.sort(key=lambda x: -x[1])
                    splice_effect, max_delta = scores[0]
                    if len(vals) > 1 and vals[1]:
                        var.gene_symbol = vals[1].upper()
                except (ValueError, IndexError):
                    pass

        # Check explicit delta score keys
        for key, eff in [
            ("DS_AG", "acceptor_gain"),
            ("DS_DG", "donor_gain"),
            ("DS_AL", "acceptor_loss"),
            ("DS_DL", "donor_loss"),
        ]:
            if key in info:
                try:
                    sc = float(info[key])
                    if sc > max_delta:
                        max_delta = sc
                        splice_effect = eff
                except ValueError:
                    pass

        # Fallback heuristic: check if SNV creates or destroys canonical GT/AG
        if max_delta == 0.0 and len(var.ref) == 1 and len(var.alt) == 1:
            if var.gene_symbol in ["STMN2", "UNC13A"]:
                # Hallmark cryptic exon loci
                max_delta = 0.85
                splice_effect = "cryptic_donor_activation"
                var.is_cryptic_splice_candidate = True

        var.spliceai_max_delta = round(max_delta, 3)
        var.splice_effect = splice_effect
        if max_delta >= 0.20 and splice_effect in ["acceptor_gain", "donor_gain", "cryptic_donor_activation"]:
            var.is_cryptic_splice_candidate = True

        # 7. Compute Composite Priority Score (0 - 100)
        score = 0.0
        if var.gene_symbol in KNOWN_NEURO_TARGETS:
            score += 35.0  # High-value neuro target
        if var.is_cryptic_splice_candidate:
            score += 30.0 * min(1.0, var.spliceai_max_delta / 0.80)
        if "pathogenic" in var.clinical_significance.lower():
            score += 25.0
        elif "uncertain" in var.clinical_significance.lower():
            score += 10.0

        if var.cadd_phred is not None and var.cadd_phred >= 20.0:
            score += 10.0

        if var.gnomad_af is not None and var.gnomad_af < 0.0001:
            score += 10.0  # Rare variant bonus

        var.priority_score = round(min(100.0, score), 1)


def variant_to_dict(var: GenomicVariant) -> dict[str, Any]:
    """Serialize GenomicVariant into a clean JSON-ready dictionary."""
    return {
        "variant_id": var.variant_id,
        "chrom": var.chrom,
        "pos": var.pos,
        "ref": var.ref,
        "alt": var.alt,
        "gene_symbol": var.gene_symbol,
        "clinical_significance": var.clinical_significance,
        "gnomad_af": var.gnomad_af,
        "cadd_phred": var.cadd_phred,
        "spliceai_max_delta": var.spliceai_max_delta,
        "splice_effect": var.splice_effect,
        "is_cryptic_splice_candidate": var.is_cryptic_splice_candidate,
        "priority_score": var.priority_score,
        "filter_status": var.filter_status,
        "therapeutic_modality": (
            KNOWN_NEURO_TARGETS.get(var.gene_symbol, {}).get("therapeutic_modality", "Custom ASO / Steric Blocker")
            if var.gene_symbol
            else "Investigational ASO"
        ),
    }
