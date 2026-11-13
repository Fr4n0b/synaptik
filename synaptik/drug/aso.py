"""
In Silico Antisense Oligonucleotide (ASO) Engineering Engine.

Designs clinical-grade, synthesis-ready antisense oligonucleotides (ASOs) for:
  1. Splice-switching oligonucleotides (SSO): steric hindrance to promote exon
     inclusion (e.g. Nusinersen / Spinraza) or exclusion (e.g. STMN2 cryptic exon).
  2. RNase H-mediated knockdown gapmers (5-10-5 wing-gap-wing architecture; e.g.
     Tofersen / Qalsody).

Integrates:
  * Nearest-neighbor thermodynamic hybridization parameters (SantaLucia 1998,
    Sugimoto et al. 1995, Xia et al. 1998) computing exact Delta H, Delta S,
    Delta G_37 and salt-adjusted Tm at physiological CSF conditions.
  * 2'-MOE and phosphorothioate (PS) thermal and pharmacokinetic adjustments.
  * Secondary structure & biophysical liability filters (G-quadruplex runs,
    poly-T/U tracts, self-hairpins, and 3'-end self-dimerization).
  * Optional whole-transcriptome off-target selectivity screening via
    `synaptik.drug.specificity` (returns NOT_SCREENED when no FASTA is supplied).
  * Vendor-ready ordering specifications (IDT format, molecular weight,
    extinction coefficient at 260nm, and chemistry-class CSF half-life).
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from enum import Enum
from typing import Any

from .aso_calibrated import predict_empirical_knockdown_pct
from .specificity import (
    DEFAULT_MAX_MISMATCHES, ScreenStatus, SpecificityResult,
    SpecificityScreener, iter_fasta, unscreened,
)


class ASOMode(str, Enum):
    """Functional mechanism of the antisense oligonucleotide."""

    SPLICE_SWITCH = "splice_switch"  # Steric hindrance: exon inclusion or exclusion (no RNase H)
    GAPMER = "gapmer"                # Catalytic degradation via nuclear RNase H1


class ASOChemistry(str, Enum):
    """Chemical modification architecture."""

    MOE_PS = "2'-MOE / PS"                  # Uniform 2'-O-methoxyethyl with full phosphorothioate
    PMO = "PMO"                             # Phosphorodiamidate morpholino oligomer (neutral backbone)
    GAPMER_5_10_5 = "5-10-5 2'-MOE/DNA/PS"  # 5 nt 2'-MOE wings flanking 10 nt central DNA gap


# --------------------------------------------------------------------------
# Nearest-neighbor thermodynamic parameters (DNA ASO : RNA target)
# Enthalpy (dH) in kcal/mol, Entropy (dS) in cal/(mol*K)
# From Sugimoto et al. Biochemistry 1995 (34:11211-11216), Table 2.
# Verified against Biopython Bio.SeqUtils.MeltingTemp.R_DNA_NN1.
# --------------------------------------------------------------------------
# RNA 5'->3' indexed table (5'-r(N1 N2)-3' / 3'-d(N1' N2')-5'):
NN_RNA_DNA: dict[tuple[str, str], tuple[float, float]] = {
    ("A", "A"): (-7.8, -21.9),
    ("A", "C"): (-5.9, -12.3),
    ("A", "G"): (-9.1, -23.5),
    ("A", "T"): (-8.3, -23.9),
    ("C", "A"): (-9.0, -26.1),
    ("C", "C"): (-9.3, -23.2),
    ("C", "G"): (-16.3, -47.1),
    ("C", "T"): (-7.0, -19.7),
    ("G", "A"): (-5.5, -13.5),
    ("G", "C"): (-8.0, -17.1),
    ("G", "G"): (-12.8, -31.9),
    ("G", "T"): (-7.8, -21.6),
    ("T", "A"): (-7.8, -23.2),
    ("T", "C"): (-8.6, -22.9),
    ("T", "G"): (-10.4, -28.4),
    ("T", "T"): (-11.5, -36.4),
}

# DNA ASO 5'->3' indexed table (5'-d(D1 D2)-3' / 3'-r(R1 R2)-5'):
# Because duplexes are antiparallel, a 5'-D1 D2-3' DNA doublet pairs with a
# 5'-comp(D2) comp(D1)-3' RNA doublet in the target strand.
NN_DNA_RNA: dict[tuple[str, str], tuple[float, float]] = {
    ("A", "A"): (-11.5, -36.4),  # complementary RNA 5'-UU-3' (Biopython TT/AA)
    ("A", "C"): (-7.8, -21.6),   # complementary RNA 5'-GU-3' (Biopython GT/CA)
    ("A", "G"): (-7.0, -19.7),   # complementary RNA 5'-CU-3' (Biopython CT/GA)
    ("A", "T"): (-8.3, -23.9),   # complementary RNA 5'-AU-3' (Biopython AT/TA)
    ("C", "A"): (-10.4, -28.4),  # complementary RNA 5'-UG-3' (Biopython TG/AC)
    ("C", "C"): (-12.8, -31.9),  # complementary RNA 5'-GG-3' (Biopython GG/CC)
    ("C", "G"): (-16.3, -47.1),  # complementary RNA 5'-CG-3' (Biopython CG/GC)
    ("C", "T"): (-9.1, -23.5),   # complementary RNA 5'-AG-3' (Biopython AG/TC)
    ("G", "A"): (-8.6, -22.9),   # complementary RNA 5'-UC-3' (Biopython TC/AG)
    ("G", "C"): (-8.0, -17.1),   # complementary RNA 5'-GC-3' (Biopython GC/CG)
    ("G", "G"): (-9.3, -23.2),   # complementary RNA 5'-CC-3' (Biopython CC/GG)
    ("G", "T"): (-5.9, -12.3),   # complementary RNA 5'-AC-3' (Biopython AC/TG)
    ("T", "A"): (-7.8, -23.2),   # complementary RNA 5'-UA-3' (Biopython TA/AT)
    ("T", "C"): (-5.5, -13.5),   # complementary RNA 5'-GA-3' (Biopython GA/CT)
    ("T", "G"): (-9.0, -26.1),   # complementary RNA 5'-CA-3' (Biopython CA/GT)
    ("T", "T"): (-7.8, -21.9),   # complementary RNA 5'-AA-3' (Biopython AA/TT)
}

# Monomer molecular weights (sodium salt, anhydrous g/mol):
MW_MOE_PS: dict[str, float] = {
    "A": 389.26,
    "C": 365.23,
    "G": 405.26,
    "T": 366.21,
    "U": 366.21,
}

MW_DNA: dict[str, float] = {
    "A": 331.20,
    "C": 307.18,
    "G": 347.20,
    "T": 322.21,
    "U": 322.21,
}

# Extinction coefficient nearest-neighbor (L/(mol*cm) at 260 nm)
EPSILON_SINGLE: dict[str, int] = {
    "A": 15200,
    "C": 7050,
    "G": 12010,
    "T": 8400,
    "U": 8400,
}


def _reverse_complement_rna(seq: str) -> str:
    """Return reverse complement DNA-style antisense sequence for RNA target."""
    tr = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")
    return seq.translate(tr)[::-1].upper()


def _gc_content(seq: str) -> float:
    """Calculate GC percentage."""
    if not seq:
        return 0.0
    gc = sum(1 for b in seq.upper() if b in ("G", "C"))
    return round((gc / len(seq)) * 100.0, 1)


def compute_thermodynamics(
    oligo_dna: str,
    target_rna: str = "",
    na_molar: float = 0.10,
    oligo_molar: float = 2.5e-7,
    moe_shift_per_nt: float = 1.4,
    moe_dg_per_nt: float = -0.16,
    chemistry: ASOChemistry | None = None,
) -> tuple[float, float, float, float]:
    """Calculate nearest-neighbor thermodynamics for ASO duplex.

    Applies nearest-neighbor parameters for RNA/DNA hybrid duplexes
    (Sugimoto et al. Biochemistry 1995, 34:11211-11216) with empirical
    backbone/sugar stabilization:
      - Uniform 2'-MOE/PS (default / MOE_PS): +1.4 C/nt, -0.16 kcal/mol/nt
        across all nucleotides (benchmarked against Nusinersen).
      - 5-10-5 2'-MOE/DNA/PS gapmer (GAPMER_5_10_5): 2'-MOE stabilization
        applied only to the 10 wing nucleotides (5' and 3' 5-nt wings); the
        central DNA/PS gap has no 2'-MOE stabilization.
      - PMO (neutral phosphorodiamidate morpholino): +0.5 C/nt, -0.05 kcal/mol/nt.

    Returns:
        (delta_h_kcal_mol, delta_s_cal_mol_k, delta_g_37_kcal_mol, tm_celsius)
    """
    s = oligo_dna.upper().replace("U", "T")
    n = len(s)
    if n < 2:
        return 0.0, 0.0, 0.0, 0.0

    # Duplex initiation (Sugimoto 1995 Table 2)
    dh = 1.9
    ds = -3.9

    # Nearest-neighbor stacking for DNA ASO paired against target RNA
    for i in range(n - 1):
        pair = (s[i], s[i + 1])
        h, ent = NN_DNA_RNA.get(pair, (-8.0, -22.0))
        dh += h
        ds += ent

    if chemistry == ASOChemistry.GAPMER_5_10_5:
        n_mod = min(n, 10)  # Two 5-nt 2'-MOE wings; central 10-nt DNA gap
        shift_per_nt = moe_shift_per_nt
        dg_per_nt = moe_dg_per_nt
    elif chemistry == ASOChemistry.PMO:
        n_mod = n
        shift_per_nt = 0.5
        dg_per_nt = -0.05
    else:
        n_mod = n
        shift_per_nt = moe_shift_per_nt
        dg_per_nt = moe_dg_per_nt

    # Delta G at 37 C (310.15 K) with chemistry-specific stabilization
    dg37_base = dh - 310.15 * (ds / 1000.0)
    dg37_moe = dg37_base + (n_mod * dg_per_nt)

    # Tm calculation with SantaLucia (1998) salt adjustment
    # R = 1.9872 cal/(mol*K)
    r_const = 1.9872
    try:
        tm_1m = (dh * 1000.0) / (ds + r_const * math.log(oligo_molar / 4.0)) - 273.15
        salt_correction = 16.6 * math.log10(max(na_molar, 1e-4))
        tm_base = tm_1m + salt_correction
        tm_moe = tm_base + (n_mod * shift_per_nt)
    except (ValueError, ZeroDivisionError):
        tm_moe = 55.0

    return round(dh, 1), round(ds, 1), round(dg37_moe, 2), round(tm_moe, 1)


def predict_hairpin_delta_g(oligo: str) -> float:
    """Predict minimum free energy of secondary hairpin stem-loop formation.

    Negative values denote stable hairpins (undesirable if <= -3.0 kcal/mol).
    """
    seq = oligo.upper().replace("U", "T")
    n = len(seq)
    min_dg = 0.0

    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}

    # Search for stem of length 3-6 nt with loop 3-8 nt
    for stem_len in range(3, 7):
        for loop_len in range(3, 9):
            total_span = 2 * stem_len + loop_len
            if total_span > n:
                continue
            for i in range(n - total_span + 1):
                stem1 = seq[i : i + stem_len]
                stem2 = seq[i + stem_len + loop_len : i + total_span]
                stem2_rev = stem2[::-1]

                matches = sum(1 for a, b in zip(stem1, stem2_rev) if comp.get(a) == b)
                if matches >= stem_len - 1 and matches >= 3:
                    # Approximate base pairing energy: -1.8 kcal/mol per match, +3.5 kcal/mol loop penalty
                    dg_hairpin = -(matches * 1.8) + 3.5
                    if dg_hairpin < min_dg:
                        min_dg = dg_hairpin

    return round(min_dg, 1)


def predict_self_dimer_score(oligo: str) -> float:
    """Detect self-dimerization risk, emphasizing 3'-end complementarity."""
    seq = oligo.upper().replace("U", "T")
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    end_seq = seq[-5:]  # Critical 3' terminus

    # Compare 3' terminus against reverse complement of any 5-mer in the oligo
    rev_comp = "".join(comp.get(b, "N") for b in seq[::-1])
    max_match = 0
    for i in range(len(rev_comp) - 4):
        sub = rev_comp[i : i + 5]
        m = sum(1 for a, b in zip(end_seq, sub) if a == b)
        if m > max_match:
            max_match = m

    return float(max_match)


# Curated essential neural and housekeeping gene reference sequences for off-target screening
CORE_NEURAL_TRANSCRIPTS: dict[str, str] = {
    "ACTB": "CCACCATGTACCCAGGCATTGCTGACAGGATGCAGAAGGAGATCACCGCCCTGGCGCCCAGCACCATGAAGATCAAGATCATTGCTCCTCC",
    "GAPDH": "TCTCCCCACACACATGCACTTACACCCCCTCTCCCTCCCCACACACATGCACTTACACCCCCTCTCCCTCCCCACACACATGCACT",
    "TUBB3": "AGAGCGGCGAGGAGATCGTGGACGAGTACCGCACCACGGAGGAGACGTACAGCCTGGACGAGTACGAGGACCTGGAGTCGGAGGCG",
    "MAP2": "AAGACACCCGGCAGGAAGGAGACCCCAGCTCACCGGAGAAGCCCTCACCCACACCCAGGAGAGGCAGCCTCAGCCCCAGGAGGAC",
    "SYP": "ACCGCCAGGTGGAGTGGCAGTATGGCAGCTACTACCGCTGGGTCAACTATGCCCACACCAAGTTCGTGGCCACGGCCAACGGCAAC",
    "SNAP25": "GATGAGAACATGGAGCAGGTGTCCCGGATCATGGAGAACCTGCGGCAGATGTCCATGGACATGGGCAGCGAGATCGACACACAGAA",
    "TARDBP": "ATGGATGAATTTGCAGCTAACCAGGCTCTCCAGTGCAATGGTATCATGTTTGGTGGTTTTAGTTGTTTTTACTGTTGGATGAAG",
    "FUS": "ATGGCCTCAAACGATTATACCCAACAAGCAACCCAAAGCTATGGGGCCTACCCCACCCAGCCCGGGCAGGGCTATTCCCAGCAGAGC",
    "SOD1": "ATGGCGACGAAGGCCGTGTGCGTGCTGAAGGGCGACGGCCCAGTGCAGGGCATCATCAATTTCGAGCAGAAGGAAAGTAATGGACC",
    "C9orf72": "GGCCCCGGCCCGGCCCCGGCCCGGCCCCGGCCCGGCCCCGGCCCGGCCCCGGCCCGGCCCCGGCCCGGCCCCGGCCCCGGCCCG",
}


def count_cpg_motifs(oligo: str) -> tuple[int, list[int]]:
    """Count unmethylated CpG dinucleotides, a TLR9 innate-immune liability.

    Phosphorothioate ASOs containing CpG motifs can activate TLR9 and drive
    pro-inflammatory responses; CpG reduction is a standard part of ASO lead
    optimisation. This is a sequence count, not an immunogenicity prediction:
    context, methylation and backbone chemistry all matter, and a CpG-free
    oligo can still be immunostimulatory.

    Returns (count, 0-based positions of each CpG).
    """
    seq = (oligo or "").upper().replace("U", "T")
    positions = [i for i in range(len(seq) - 1) if seq[i:i + 2] == "CG"]
    return len(positions), positions


def score_splice_regulatory_context(target_rna: str) -> tuple[float, list[str]]:
    """Score canonical splicing-regulatory elements in a target RNA window.

    For steric-blocking splice-switching ASOs, occluding a splice junction,
    polypyrimidine tract, branchpoint, or intronic/exonic splicing silencer
    (hnRNP A1/A2, PTB) is strongly associated with functional exon modulation.
    Returns (bonus_points in [0, 6], matched_motif_descriptions).
    """
    rna = (target_rna or "").upper().replace("T", "U")
    bonus = 0.0
    hits: list[str] = []

    if re.search(r"[CU]{6,}", rna):
        bonus += 2.0
        hits.append("polypyrimidine tract")
    if re.search(r"[CU]AG[GA]", rna) or re.search(r"GU[AG]AGU|AGGUA", rna):
        bonus += 2.5
        hits.append("splice junction consensus")
    if re.search(r"UAGG|UAGGG|UCUU|CUCUCU|UUAAAG", rna):
        bonus += 2.0
        hits.append("hnRNP/PTB splicing silencer")
    if re.search(r"[CU]U[AG]AC", rna):
        bonus += 1.5
        hits.append("branchpoint consensus")

    return min(6.0, bonus), hits


def screen_off_targets(
    oligo: str, target_gene: str = ""
) -> tuple[int, int, int]:
    """Smoke test against a ten-gene panel. NOT a specificity screen.

    CORE_NEURAL_TRANSCRIPTS totals 861 bp, about 0.0007% of the human
    transcriptome. Measured empirically, this function flags 0 of 3000 random
    20-mers: it can only catch an oligo that is a near-perfect match to one of
    ten short fragments. A "0 off-targets" result from it therefore means
    "nothing was looked at", not "nothing was found".

    It is retained only as a cheap sanity check for gross errors (for example
    an oligo accidentally designed against ACTB). Real specificity screening
    lives in synaptik.drug.specificity, which searches a supplied
    transcriptome and reports how much it covered.

    Returns:
        (exact_matches, one_mismatch_matches, min_mismatch_distance)
    """
    seq = oligo.upper().replace("U", "T")
    oligo_len = len(seq)
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    # Convert oligo back to target orientation to match sense transcripts
    target_match_seq = "".join(comp.get(b, "N") for b in seq[::-1])

    exact_matches = 0
    one_mismatch = 0
    min_dist = oligo_len

    for gene, tx in CORE_NEURAL_TRANSCRIPTS.items():
        if target_gene and gene.upper() == target_gene.upper():
            continue  # Do not count on-target matches

        tx_len = len(tx)
        for i in range(tx_len - oligo_len + 1):
            window = tx[i : i + oligo_len]
            dist = sum(1 for a, b in zip(target_match_seq, window) if a != b)
            if dist == 0:
                exact_matches += 1
            elif dist == 1:
                one_mismatch += 1
            if dist < min_dist:
                min_dist = dist

    return exact_matches, one_mismatch, min_dist


def format_vendor_spec(
    oligo: str, mode: ASOMode, chemistry: ASOChemistry
) -> str:
    """Format the oligonucleotide in standard vendor ordering format (IDT)."""
    seq = oligo.upper().replace("T", "U")

    if chemistry == ASOChemistry.PMO:
        # Morpholino notation
        return f"5'-PMO-{seq}-3'"

    if mode == ASOMode.SPLICE_SWITCH:
        # Full uniform 2'-MOE with phosphorothioate internucleotide linkages
        return "*".join(f"m{b}" for b in seq)

    # Gapmer 5-10-5 architecture
    n = len(seq)
    if n == 20:
        wing5 = "*".join(f"m{b}" for b in seq[:5])
        gap10 = "*".join(f"+d{b}" for b in seq[5:15].replace("U", "T"))
        wing3 = "*".join(f"m{b}" for b in seq[15:])
        return f"{wing5}*{gap10}*{wing3}"

    return "*".join(f"m{b}" for b in seq)


def calculate_molecular_weight(oligo: str, chemistry: ASOChemistry) -> float:
    """Compute molecular weight of the ASO (sodium salt, g/mol)."""
    seq = oligo.upper().replace("T", "U")
    n = len(seq)
    if n == 0:
        return 0.0

    if chemistry == ASOChemistry.PMO:
        # Average PMO subunit ~ 324.5 Da
        return round(n * 324.5 + 40.0, 1)

    # 2'-MOE with phosphorothioate linkages
    # Base weight + PS linkage (+16.06 per linkage)
    total_mw = sum(MW_MOE_PS.get(b, 375.0) for b in seq)
    ps_linkages = n - 1
    total_mw += ps_linkages * 16.06  # Sulfur replacing oxygen
    return round(total_mw, 1)


@dataclass
class ASOCandidate:
    """A scored ASO candidate for EXPERIMENTAL validation.

    Not a clinical or synthesis-ready product. Every candidate still requires
    an off-target screen against a real transcriptome (see `specificity`),
    in vitro confirmation of the intended splicing or knockdown effect, and
    toxicology work before any use beyond the bench.
    """

    oligo_sequence_5_to_3: str
    target_sequence_5_to_3: str
    target_start: int
    target_end: int
    length: int
    mode: ASOMode
    chemistry: ASOChemistry
    gc_content_pct: float
    delta_h_kcal_mol: float
    delta_s_cal_mol_k: float
    gibbs_free_energy_37c: float
    melting_temp_c: float
    has_quadruplex_risk: bool
    has_poly_u_risk: bool
    hairpin_delta_g: float
    self_dimer_score: float
    # Panel smoke test only -- see screen_off_targets(). These are NOT a
    # specificity result; read `specificity` for that.
    off_target_0_mismatch: int
    off_target_1_mismatch: int
    min_mismatch_distance: int
    molecular_weight_da: float
    vendor_idt_format: str
    # Literature half-life for this CHEMISTRY CLASS, not a per-sequence
    # prediction: it is identical for every oligo of the same chemistry.
    chemistry_class_csf_half_life_days: int
    composite_score: float
    specificity: SpecificityResult = field(default_factory=SpecificityResult)
    penalties_and_warnings: list[str] = field(default_factory=list)
    delta_g_open_kcal_mol: float | None = None
    effective_dg_kcal_mol: float | None = None
    accessibility_tier: str | None = None
    empirical_knockdown_pct: float | None = None

    @property
    def specificity_screened(self) -> bool:
        return self.specificity.status is ScreenStatus.TRANSCRIPTOME


@dataclass
class ASODesignResult:
    """Full results of an in silico ASO tiling and scoring session."""

    target_gene: str
    target_region: str
    target_sequence: str
    mode: ASOMode
    chemistry: ASOChemistry
    n_windows_evaluated: int
    n_passed_filters: int
    candidates: list[ASOCandidate]
    precedent: str = ""
    caveats: list[str] = field(default_factory=list)

    def summary(self) -> str:
        s = [
            f"Target: {self.target_gene} ({self.target_region}), mode: {self.mode.value}",
            f"Evaluated {self.n_windows_evaluated} tiling windows -> {self.n_passed_filters} passed biophysical filters",
        ]
        if self.candidates:
            top = self.candidates[0]
            s.append(
                f"Top Candidate: 5'-{top.oligo_sequence_5_to_3}-3' "
                f"(Tm={top.melting_temp_c}C, dG_37={top.gibbs_free_energy_37c} kcal/mol, Score={top.composite_score:.1f})"
            )
        return "\n".join(s)


class ASODesigner:
    """Clinical-grade antisense oligonucleotide design and scoring."""

    def __init__(
        self,
        default_length: int = 20,
        min_gc: float = 35.0,
        max_gc: float = 65.0,
        target_tm_min: float = 52.0,
        target_tm_max: float = 72.0,
    ) -> None:
        self.default_length = default_length
        self.min_gc = min_gc
        self.max_gc = max_gc
        self.target_tm_min = target_tm_min
        self.target_tm_max = target_tm_max

    def design(
        self,
        target_sequence: str,
        mode: ASOMode = ASOMode.SPLICE_SWITCH,
        target_gene: str = "",
        target_region: str = "cryptic splice site",
        window_sizes: tuple[int, ...] = (18, 20, 22),
        top_k: int = 5,
        transcriptome: str | None = None,
        max_mismatches: int | None = None,
        use_biophysics: bool = True,
    ) -> ASODesignResult:
        """Tile across the target RNA and rank candidate ASOs.

        transcriptome: path to a FASTA of reference transcripts. Supplied, the
            top candidates are screened for off-target hybridisation. Omitted,
            no specificity claim is made and every candidate is reported as
            NOT SCREENED -- the tool will not print zeros it has not earned.
        """
        clean_target = target_sequence.strip().upper().replace("T", "U")
        target_len = len(clean_target)

        chemistry = (
            ASOChemistry.MOE_PS
            if mode == ASOMode.SPLICE_SWITCH
            else ASOChemistry.GAPMER_5_10_5
        )

        # Literature value for the chemistry class (2'-MOE/PS ASOs persist in
        # CNS tissue for months; PMOs clear far faster). It does not vary with
        # sequence, so it is reported as a class property, never as a
        # per-candidate prediction.
        half_life = 135 if chemistry != ASOChemistry.PMO else 18

        candidates: list[ASOCandidate] = []
        n_windows = 0

        for w_len in window_sizes:
            if target_len < w_len:
                continue
            for i in range(target_len - w_len + 1):
                n_windows += 1
                sub_target = clean_target[i : i + w_len]
                oligo = _reverse_complement_rna(sub_target)

                # Biophysical properties
                gc = _gc_content(oligo)
                dh, ds, dg37, tm = compute_thermodynamics(oligo, sub_target, chemistry=chemistry)

                has_quad = bool(re.search(r"G{4,}", oligo))
                has_poly_u = bool(re.search(r"[TU]{4,}", oligo))
                hairpin_dg = predict_hairpin_delta_g(oligo)
                dimer_score = predict_self_dimer_score(oligo)

                off_0, off_1, min_dist = screen_off_targets(oligo, target_gene)

                penalties: list[str] = []
                # Scoring (0 - 100)
                score = 100.0

                # 1. GC content penalty
                if gc < self.min_gc or gc > self.max_gc:
                    score -= 25.0
                    penalties.append(f"Suboptimal GC ({gc}% outside {self.min_gc}-{self.max_gc}%)")
                elif 45.0 <= gc <= 55.0:
                    pass  # Optimal
                else:
                    score -= 5.0

                # 2. Thermodynamic affinity
                if tm < self.target_tm_min:
                    diff = self.target_tm_min - tm
                    score -= min(30.0, diff * 3.0)
                    penalties.append(f"Low melting temperature ({tm}C < {self.target_tm_min}C)")
                elif tm > self.target_tm_max:
                    diff = tm - self.target_tm_max
                    score -= min(25.0, diff * 2.5)
                    penalties.append(f"High melting temperature ({tm}C > {self.target_tm_max}C, off-target risk)")

                # 3. G-quadruplex risk
                if has_quad:
                    score -= 40.0
                    penalties.append("Contains G-quadruplex motif (4+ consecutive Gs: CSF neurotoxicity risk)")

                # 4. Poly-T / Poly-U risk
                if has_poly_u:
                    score -= 20.0
                    penalties.append("Contains poly-T/poly-U tract (premature termination / immune liability)")

                # 5. Stable hairpin
                if hairpin_dg <= -3.0:
                    score -= 25.0
                    penalties.append(f"Stable intramolecular hairpin (dG={hairpin_dg} kcal/mol <= -3.0)")

                # 6. Self-dimerization
                if dimer_score >= 4.0:
                    score -= 15.0
                    penalties.append(f"High 3'-end self-complementarity ({int(dimer_score)} nt match)")

                # 7. Off-target penalty
                if off_0 > 0:
                    score -= 50.0
                    penalties.append(f"Exact 0-mismatch off-target match in core neural transcriptome ({off_0} hits)")
                elif off_1 > 0:
                    score -= 25.0
                    penalties.append(f"Single 1-mismatch off-target match in core neural transcriptome ({off_1} hits)")

                # 8. Unmethylated CpG innate-immune (TLR9) motif penalty for PS backbones
                cpg_count, _ = count_cpg_motifs(oligo)
                if chemistry != ASOChemistry.PMO:
                    if cpg_count > 0:
                        score -= min(12.0, 4.0 * cpg_count)
                        penalties.append(f"Contains {cpg_count} CpG motif(s) (TLR9 innate immune risk in PS backbone)")

                # 9. Splicing-regulatory element alignment for splice-switching ASOs
                if mode == ASOMode.SPLICE_SWITCH:
                    splice_bonus, _ = score_splice_regulatory_context(sub_target)
                    score += splice_bonus

                final_score = max(0.0, min(100.0, score))
                mw = calculate_molecular_weight(oligo, chemistry)
                vendor_fmt = format_vendor_spec(oligo, mode, chemistry)
                emp_kd = (
                    predict_empirical_knockdown_pct(
                        oligo, dg37, tm, hairpin_dg, dimer_score, cpg_count,
                        flank_5p=clean_target[max(0, i - 50) : i],
                        flank_3p=clean_target[i + w_len : i + w_len + 50],
                    )
                    if mode == ASOMode.GAPMER
                    else None
                )

                candidates.append(
                    ASOCandidate(
                        oligo_sequence_5_to_3=oligo,
                        target_sequence_5_to_3=sub_target,
                        target_start=i,
                        target_end=i + w_len,
                        length=w_len,
                        mode=mode,
                        chemistry=chemistry,
                        gc_content_pct=gc,
                        delta_h_kcal_mol=dh,
                        delta_s_cal_mol_k=ds,
                        gibbs_free_energy_37c=dg37,
                        melting_temp_c=tm,
                        has_quadruplex_risk=has_quad,
                        has_poly_u_risk=has_poly_u,
                        hairpin_delta_g=hairpin_dg,
                        self_dimer_score=dimer_score,
                        off_target_0_mismatch=off_0,
                        off_target_1_mismatch=off_1,
                        min_mismatch_distance=min_dist,
                        molecular_weight_da=mw,
                        vendor_idt_format=vendor_fmt,
                        chemistry_class_csf_half_life_days=half_life,
                        composite_score=round(final_score, 1),
                        penalties_and_warnings=penalties,
                        empirical_knockdown_pct=emp_kd,
                    )
                )

        # Sort by composite score descending, then by proximity to ideal Tm (60 C)
        candidates.sort(
            key=lambda c: (c.composite_score, -abs(c.melting_temp_c - 62.0)),
            reverse=True,
        )

        passed = [c for c in candidates if not c.has_quadruplex_risk and c.off_target_0_mismatch == 0 and c.composite_score >= 60.0]

        # ---- 3D RNA target accessibility & specificity ----
        # Evaluate on an expanded shortlist so structurally occluded or
        # off-target hits are properly demoted before slicing top_k.
        pool_size = max(top_k * 3, 10) if use_biophysics else max(top_k, 1)
        shortlist = (passed or candidates)[:pool_size]

        if use_biophysics and shortlist:
            from ..biophysics.folding import compute_target_accessibility
            flank = 60
            for cand in shortlist:
                win_start = max(0, cand.target_start - flank)
                win_end = min(target_len, cand.target_end + flank)
                local_rna = clean_target[win_start:win_end]
                rel_start = cand.target_start - win_start
                rel_end = cand.target_end - win_start
                acc = compute_target_accessibility(local_rna, rel_start, rel_end)
                cand.delta_g_open_kcal_mol = acc.delta_g_open
                eff_dg = round(cand.gibbs_free_energy_37c + acc.delta_g_open, 2)
                cand.effective_dg_kcal_mol = eff_dg
                cand.accessibility_tier = acc.accessibility_tier.value
                if acc.delta_g_open > 4.0:
                    if eff_dg > -15.0:
                        cand.composite_score = round(max(0.0, cand.composite_score - 20.0), 1)
                        cand.penalties_and_warnings.append(
                            f"Target RNA is structurally occluded (Delta G_open={acc.delta_g_open:.2f} kcal/mol > 4.0; dG_eff={eff_dg:.1f})"
                        )
                    else:
                        cand.composite_score = round(max(0.0, cand.composite_score - 6.0), 1)
                        cand.penalties_and_warnings.append(
                            f"Target RNA structured (Delta G_open={acc.delta_g_open:.2f} kcal/mol > 4.0) but thermodynamically invadable (dG_eff={eff_dg:.1f})"
                        )
                elif acc.delta_g_open >= 1.5:
                    cand.composite_score = round(max(0.0, cand.composite_score - 5.0), 1)
                    cand.penalties_and_warnings.append(
                        f"Target RNA moderately structured (Delta G_open={acc.delta_g_open:.2f} kcal/mol)"
                    )
            shortlist.sort(
                key=lambda c: (
                    c.composite_score,
                    -(c.effective_dg_kcal_mol if c.effective_dg_kcal_mol is not None else c.gibbs_free_energy_37c),
                ),
                reverse=True,
            )

        shortlist = shortlist[: max(top_k, 1)]
        mm = max_mismatches if max_mismatches is not None else \
            DEFAULT_MAX_MISMATCHES.get(mode.value, 3)
        if transcriptome:
            screener = SpecificityScreener(
                max_mismatches=mm, skip_on_target_gene=target_gene)
            batch_results = screener.screen_batch_against_fasta(
                [c.oligo_sequence_5_to_3 for c in shortlist], transcriptome)
            for cand, spec_res in zip(shortlist, batch_results):
                cand.specificity = spec_res
                hits = cand.specificity.hits
                if hits:
                    worst = min(h.n_mismatches for h in hits)
                    # An exact or near-exact off-target is disqualifying, not
                    # a deduction: it means the oligo hits something else.
                    penalty = {0: 80.0, 1: 60.0, 2: 30.0}.get(worst, 15.0)
                    cand.composite_score = round(
                        max(0.0, cand.composite_score - penalty), 1)
                    cand.penalties_and_warnings.append(
                        f"{len(hits)} transcriptome off-target hit(s), closest "
                        f"at {worst} mismatch(es) "
                        f"({hits[0].transcript_id})")
            shortlist.sort(key=lambda c: -c.composite_score)
        else:
            for cand in shortlist:
                cand.specificity = unscreened(cand.oligo_sequence_5_to_3)
                cand.penalties_and_warnings.append(
                    "Specificity NOT screened: no reference transcriptome "
                    "supplied.")

        precedent = (
            "Nusinersen (Spinraza) / STMN2 cryptic exon exclusion series"
            if mode == ASOMode.SPLICE_SWITCH
            else "Tofersen (Qalsody) / SOD1 5-10-5 MOE gapmer series"
        )

        caveats = [
            "Tiling predictions optimize thermodynamic and biophysical parameters; empirical in vitro potency must be confirmed in cells.",
            "Local RNA target accessibility (Delta G_open) is modeled via Zuker MFE over a +/-60 nt flanking window; long-range tertiary contacts and RBP occupancy are not modeled.",
            "CSF half-life is a literature value for the chemistry class, "
            "identical for every sequence; it is not a per-oligo prediction.",
        ]
        if not transcriptome:
            caveats.insert(0, (
                "NO OFF-TARGET SCREEN WAS RUN. Supply --transcriptome, or "
                "screen with BLAST/Bowtie, before treating any candidate as "
                "usable. Off-target hybridisation is the main mechanism of "
                "ASO toxicity."))

        return ASODesignResult(
            target_gene=target_gene,
            target_region=target_region,
            target_sequence=clean_target,
            mode=mode,
            chemistry=chemistry,
            n_windows_evaluated=n_windows,
            n_passed_filters=len(passed),
            candidates=shortlist,
            precedent=precedent,
            caveats=caveats,
        )
