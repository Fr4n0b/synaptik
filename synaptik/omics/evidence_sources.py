"""
Population-genetics and functional-genomics evidence.

These supply the evidence classes that literature and sequence models cannot:

  POPULATION  gnomAD constraint. Whether a gene tolerates loss of function in
              ~800,000 sequenced humans is a natural experiment nobody had to
              fund. A gene with LOEUF < 0.35 is under strong selection; a
              phenotype attributed to its partial loss is a priori more
              plausible. Conversely, a gene that tolerates homozygous LoF in
              healthy people is a hard constraint on any haploinsufficiency
              story, and one that literature mining alone will never surface.

  FUNCTIONAL  DepMap CRISPR screens (~1,100 cell lines), and where available
              neuronal CRISPR screens. Direct perturbation is the most causal
              evidence class available without doing the experiment yourself.

  EXPRESSION  GTEx tissue specificity, Allen Brain Atlas, Human Protein Atlas.
              A gene not expressed in neurons is not causing a neuronal
              phenotype through a cell-autonomous mechanism, whatever the
              literature says.

  CLINICAL    ClinVar variant interpretations, OMIM gene-disease links. What
              clinical geneticists have concluded from patients.

All of these are free, open, and API-accessible. Together they turn a
literature-only tool into something that can say "the field believes X, but
population genetics is inconsistent with X" -- which is the kind of statement
that actually moves a field.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

from ..evidence.reconcile import EvidenceClass, EvidenceItem
from ..ingest.http import PoliteSession

log = logging.getLogger(__name__)

GNOMAD_GRAPHQL = "https://gnomad.broadinstitute.org/api"
DEPMAP_API = "https://depmap.org/portal/api"
CLINVAR_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
HPA_API = "https://www.proteinatlas.org"
MYGENE = "https://mygene.info/v3"


# --------------------------------------------------------------------------
# gnomAD constraint
# --------------------------------------------------------------------------


@dataclass
class GeneConstraint:
    """gnomAD constraint metrics for one gene."""

    gene_symbol: str
    ensembl_id: str | None = None
    pli: float | None = None            # P(loss-of-function intolerant)
    loeuf: float | None = None          # LoF observed/expected upper bound
    oe_lof: float | None = None
    oe_mis: float | None = None
    mis_z: float | None = None
    n_lof_observed: int | None = None
    n_lof_expected: float | None = None
    n_homozygous_lof: int | None = None
    error: str | None = None

    def constraint_class(self) -> str:
        """LOEUF decile interpretation, per gnomAD's own guidance."""
        if self.loeuf is None:
            return "unknown"
        if self.loeuf < 0.35:
            return "highly_constrained"
        if self.loeuf < 0.60:
            return "constrained"
        if self.loeuf < 1.00:
            return "moderately_tolerant"
        return "tolerant"

    def interpretation(self) -> str:
        c = self.constraint_class()
        if c == "unknown":
            return "No constraint data available for this gene."
        if c == "highly_constrained":
            return (
                f"Strongly depleted of loss-of-function variants in the "
                f"population (LOEUF {self.loeuf:.2f}). Heterozygous loss is "
                f"likely to be consequential, consistent with a "
                f"haploinsufficiency mechanism."
            )
        if c == "constrained":
            return (
                f"Moderately constrained (LOEUF {self.loeuf:.2f}). Partial "
                f"loss of function is plausible as a disease mechanism."
            )
        if c == "moderately_tolerant":
            return (
                f"Weakly constrained (LOEUF {self.loeuf:.2f}). A "
                f"haploinsufficiency mechanism is possible but not supported "
                f"by population data."
            )
        note = (
            f" {self.n_homozygous_lof} homozygous LoF individuals are present "
            f"in gnomAD."
            if self.n_homozygous_lof else ""
        )
        return (
            f"Tolerant of loss of function in healthy populations "
            f"(LOEUF {self.loeuf:.2f}).{note} This argues against simple "
            f"haploinsufficiency and should be reconciled with any published "
            f"claim that heterozygous loss causes disease."
        )

    def to_evidence(self, claim_is_lof_mechanism: bool = True) -> EvidenceItem:
        c = self.constraint_class()
        if c == "unknown":
            direction, conf = "neutral", 0.1
        elif c in ("highly_constrained", "constrained"):
            direction = "supports" if claim_is_lof_mechanism else "neutral"
            conf = 0.8 if c == "highly_constrained" else 0.6
        elif c == "tolerant":
            direction = "refutes" if claim_is_lof_mechanism else "neutral"
            conf = 0.75
        else:
            direction, conf = "neutral", 0.4

        caveats = []
        if self.n_lof_expected is not None and self.n_lof_expected < 10:
            caveats.append(
                "Few LoF variants expected for this gene, so the constraint "
                "estimate is imprecise."
            )
        caveats.append(
            "Constraint reflects selection on the whole gene, not on the "
            "specific variant or isoform in question."
        )
        return EvidenceItem(
            source_class=EvidenceClass.POPULATION,
            source_name="gnomAD constraint",
            claim_direction=direction,
            magnitude=self.loeuf,
            magnitude_unit="LOEUF",
            confidence=conf,
            provenance=f"gnomAD:{self.gene_symbol}",
            caveats=caveats,
            raw={"pLI": self.pli, "LOEUF": self.loeuf,
                 "constraint_class": c},
        )


GNOMAD_QUERY = """
query GeneConstraint($symbol: String!) {
  gene(gene_symbol: $symbol, reference_genome: GRCh38) {
    gene_id
    symbol
    gnomad_constraint {
      pLI
      oe_lof
      oe_lof_upper
      oe_mis
      mis_z
      obs_lof
      exp_lof
    }
  }
}
"""


class GnomADClient:
    """gnomAD GraphQL client. No key required; be polite."""

    def __init__(self, session: PoliteSession):
        self.s = session

    def constraint(self, gene_symbol: str) -> GeneConstraint:
        out = GeneConstraint(gene_symbol=gene_symbol)
        payload = {
            "query": GNOMAD_QUERY,
            "variables": {"symbol": gene_symbol.upper()},
        }
        data = self._post(payload)
        if not data:
            out.error = "gnomAD request failed or returned nothing"
            return out
        gene = ((data.get("data") or {}).get("gene")) or {}
        if not gene:
            out.error = f"gene {gene_symbol} not found in gnomAD"
            return out
        out.ensembl_id = gene.get("gene_id")
        c = gene.get("gnomad_constraint") or {}
        out.pli = c.get("pLI")
        out.loeuf = c.get("oe_lof_upper")
        out.oe_lof = c.get("oe_lof")
        out.oe_mis = c.get("oe_mis")
        out.mis_z = c.get("mis_z")
        out.n_lof_observed = c.get("obs_lof")
        out.n_lof_expected = c.get("exp_lof")
        return out

    def _post(self, payload: dict) -> dict | None:
        import json as _json

        key = "gnomad:" + _json.dumps(payload, sort_keys=True)[:300]
        if self.s.use_cache and self.s.store:
            cached = self.s.store.cache_get(key)
            if cached is not None:
                self.s.n_cache_hits += 1
                return cached
        try:
            self.s.limiter.wait("gnomad.broadinstitute.org")
            resp = self.s.session.post(
                GNOMAD_GRAPHQL, json=payload, timeout=self.s.timeout
            )
            if resp.status_code != 200:
                log.warning("gnomAD returned HTTP %s", resp.status_code)
                return None
            data = resp.json()
            if self.s.use_cache and self.s.store:
                self.s.store.cache_put(key, data)
            return data
        except Exception as e:
            log.warning("gnomAD request failed: %s", e)
            return None


# --------------------------------------------------------------------------
# Expression specificity
# --------------------------------------------------------------------------


@dataclass
class ExpressionProfile:
    """Where a gene is expressed. Constrains cell-autonomous mechanisms."""

    gene_symbol: str
    brain_enriched: bool | None = None
    neuron_expressed: bool | None = None
    max_tissue: str | None = None
    brain_tpm: float | None = None
    tissue_specificity: str | None = None
    ipsc_neuron_expressed: bool | None = None
    error: str | None = None

    def interpretation(self) -> str:
        if self.error:
            return f"No expression data: {self.error}"
        if self.neuron_expressed is False:
            return (
                f"{self.gene_symbol} is not appreciably expressed in neurons. "
                f"A cell-autonomous neuronal mechanism is implausible; any "
                f"published neuronal phenotype likely acts through another "
                f"cell type or an indirect route."
            )
        if self.brain_enriched:
            return (
                f"{self.gene_symbol} is brain-enriched"
                + (f" (highest in {self.max_tissue})" if self.max_tissue else "")
                + ". Consistent with a neuronal mechanism."
            )
        return (
            f"{self.gene_symbol} is expressed in neurons but is not "
            f"brain-specific, so a neuronal phenotype is possible but the "
            f"expression pattern does not by itself support it."
        )

    def to_evidence(self) -> EvidenceItem:
        if self.error or self.neuron_expressed is None:
            direction, conf = "neutral", 0.1
        elif self.neuron_expressed is False:
            direction, conf = "refutes", 0.8
        elif self.brain_enriched:
            direction, conf = "supports", 0.6
        else:
            direction, conf = "neutral", 0.35
        return EvidenceItem(
            source_class=EvidenceClass.EXPRESSION,
            source_name="tissue expression",
            claim_direction=direction,
            magnitude=self.brain_tpm,
            magnitude_unit="TPM",
            confidence=conf,
            provenance=self.gene_symbol,
            caveats=[
                "Bulk tissue expression does not resolve cell type; a gene "
                "with low bulk brain TPM may still be highly expressed in a "
                "rare neuronal subtype.",
                "iPSC-derived neurons are developmentally immature and their "
                "expression profile differs from adult brain.",
            ],
        )


class ExpressionClient:
    """Expression lookup via mygene.info, which aggregates several sources."""

    def __init__(self, session: PoliteSession):
        self.s = session

    def profile(self, gene_symbol: str) -> ExpressionProfile:
        out = ExpressionProfile(gene_symbol=gene_symbol)
        data = self.s.get_json(
            f"{MYGENE}/query",
            {"q": f"symbol:{gene_symbol}", "species": "human",
             "fields": "symbol,name,summary,type_of_gene,genomic_pos", "size": 1},
        )
        if not data or not data.get("hits"):
            out.error = f"gene {gene_symbol} not found"
            return out
        hit = data["hits"][0]
        summary = (hit.get("summary") or "").lower()

        # Heuristic from the RefSeq summary. Coarse, and labelled as such --
        # a proper implementation would query GTEx or HPA directly, which
        # requires bulk downloads rather than a light API call.
        neural_terms = ("neuron", "neuronal", "brain", "synap", "axon",
                        "dendrit", "cortex", "cortical", "cns",
                        "nervous system")
        out.brain_enriched = any(t in summary for t in neural_terms)
        out.neuron_expressed = out.brain_enriched or None
        if not summary:
            out.error = "no functional summary available"
        return out


# --------------------------------------------------------------------------
# ClinVar
# --------------------------------------------------------------------------


@dataclass
class ClinicalEvidence:
    """What clinical genetics says about a gene or variant."""

    gene_symbol: str
    variant_id: str | None = None
    n_pathogenic: int = 0
    n_likely_pathogenic: int = 0
    n_vus: int = 0
    n_benign: int = 0
    associated_conditions: list[str] = field(default_factory=list)
    variant_classification: str | None = None
    error: str | None = None

    def has_clinical_support(self) -> bool:
        return (self.n_pathogenic + self.n_likely_pathogenic) > 0

    def interpretation(self) -> str:
        if self.error:
            return f"No ClinVar data: {self.error}"
        if self.variant_classification:
            return (
                f"This variant is classified as "
                f"{self.variant_classification} in ClinVar."
            )
        if self.has_clinical_support():
            conds = ", ".join(self.associated_conditions[:3]) or "reported conditions"
            return (
                f"{self.n_pathogenic} pathogenic and "
                f"{self.n_likely_pathogenic} likely-pathogenic variants are "
                f"reported in {self.gene_symbol}, associated with {conds}. "
                f"The gene has established clinical relevance."
            )
        if self.n_vus:
            return (
                f"{self.n_vus} variants of uncertain significance are reported "
                f"in {self.gene_symbol}, but none classified pathogenic. "
                f"Clinical relevance is not established."
            )
        return f"No ClinVar entries for {self.gene_symbol}."

    def to_evidence(self) -> EvidenceItem:
        if self.error:
            direction, conf = "neutral", 0.1
        elif self.variant_classification:
            vc = self.variant_classification.lower()
            if "pathogenic" in vc and "benign" not in vc:
                direction, conf = "supports", 0.9
            elif "benign" in vc:
                direction, conf = "refutes", 0.85
            else:
                direction, conf = "neutral", 0.3
        elif self.has_clinical_support():
            direction, conf = "supports", 0.7
        elif self.n_vus > 5 and self.n_benign > self.n_vus:
            direction, conf = "refutes", 0.4
        else:
            direction, conf = "neutral", 0.2

        return EvidenceItem(
            source_class=EvidenceClass.CLINICAL,
            source_name="ClinVar",
            claim_direction=direction,
            confidence=conf,
            n_independent=max(1, self.n_pathogenic + self.n_likely_pathogenic),
            provenance=self.variant_id or self.gene_symbol,
            caveats=[
                "ClinVar classifications vary in evidence quality; check the "
                "review status of individual submissions.",
                "Gene-level clinical relevance does not establish that a "
                "particular variant or mechanism is causal.",
            ],
        )


class ClinVarClient:
    """ClinVar via NCBI E-utilities."""

    def __init__(self, session: PoliteSession, api_key: str | None = None):
        self.s = session
        self.api_key = api_key

    def gene_summary(self, gene_symbol: str) -> ClinicalEvidence:
        out = ClinicalEvidence(gene_symbol=gene_symbol)
        counts = {}
        for label, term in (
            ("pathogenic", "pathogenic[CLNSIG]"),
            ("likely_pathogenic", '"likely pathogenic"[CLNSIG]'),
            ("vus", '"uncertain significance"[CLNSIG]'),
            ("benign", "benign[CLNSIG]"),
        ):
            params = {
                "db": "clinvar",
                "term": f"{gene_symbol}[gene] AND {term}",
                "retmode": "json", "retmax": 0,
            }
            if self.api_key:
                params["api_key"] = self.api_key
            data = self.s.get_json(f"{CLINVAR_EUTILS}/esearch.fcgi", params)
            try:
                counts[label] = int(
                    (data or {}).get("esearchresult", {}).get("count", 0)
                )
            except (TypeError, ValueError):
                counts[label] = 0

        out.n_pathogenic = counts.get("pathogenic", 0)
        out.n_likely_pathogenic = counts.get("likely_pathogenic", 0)
        out.n_vus = counts.get("vus", 0)
        out.n_benign = counts.get("benign", 0)
        if sum(counts.values()) == 0:
            out.error = "no ClinVar records found"
        return out


# --------------------------------------------------------------------------
# DepMap functional screens
# --------------------------------------------------------------------------


@dataclass
class FunctionalEvidence:
    """CRISPR screen data. The most causal evidence class available."""

    gene_symbol: str
    is_common_essential: bool | None = None
    mean_gene_effect: float | None = None
    n_dependent_lines: int | None = None
    n_lines_tested: int | None = None
    neuronal_screen_hit: bool | None = None
    screen_context: str | None = None
    error: str | None = None

    def interpretation(self) -> str:
        if self.error:
            return f"No functional screen data: {self.error}"
        if self.is_common_essential:
            return (
                f"{self.gene_symbol} is commonly essential across cell lines. "
                f"Any phenotype from its loss may reflect general cell "
                f"toxicity rather than a neuron-specific mechanism -- a "
                f"critical confound for viability and morphology readouts."
            )
        if self.mean_gene_effect is not None and self.mean_gene_effect > -0.2:
            return (
                f"{self.gene_symbol} shows little fitness effect in CRISPR "
                f"screens (mean gene effect {self.mean_gene_effect:.2f}), so "
                f"observed phenotypes are unlikely to be secondary to "
                f"reduced viability."
            )
        return (
            f"{self.gene_symbol} shows a moderate fitness effect "
            f"(mean {self.mean_gene_effect:.2f})."
            if self.mean_gene_effect is not None else
            f"Functional screen data available for {self.gene_symbol}."
        )

    def to_evidence(self) -> EvidenceItem:
        caveats = [
            "DepMap screens are in cancer cell lines, not neurons. Essentiality "
            "transfers imperfectly to post-mitotic cells.",
        ]
        if self.is_common_essential:
            caveats.append(
                "Common essentiality is a major confound: phenotypes may be "
                "downstream of impaired viability rather than the mechanism "
                "of interest."
            )
        direction = "neutral"
        conf = 0.3
        if self.neuronal_screen_hit:
            direction, conf = "supports", 0.85
        elif self.is_common_essential:
            direction, conf = "neutral", 0.5
        return EvidenceItem(
            source_class=EvidenceClass.FUNCTIONAL,
            source_name="CRISPR screen",
            claim_direction=direction,
            magnitude=self.mean_gene_effect,
            magnitude_unit="gene_effect",
            confidence=conf,
            provenance=self.screen_context or "DepMap",
            caveats=caveats,
        )


class DepMapClient:
    """DepMap gene-effect lookup.

    DepMap's public API is limited; the full dataset is distributed as CSV
    downloads. This client supports a local CSV, which is the realistic
    deployment: download CRISPRGeneEffect.csv once (~400MB) and query it.
    """

    def __init__(self, session: PoliteSession,
                 local_gene_effect_csv: str | None = None):
        self.s = session
        self.csv_path = local_gene_effect_csv
        self._cache: dict[str, FunctionalEvidence] = {}
        self._common_essentials: set[str] = set()

    def load_common_essentials(self, path: str) -> int:
        """Load DepMap's CommonEssentials.csv."""
        try:
            import csv
            with open(path, newline="", encoding="utf-8") as fh:
                for row in csv.DictReader(fh):
                    g = (row.get("Gene") or "").split(" ")[0].strip()
                    if g:
                        self._common_essentials.add(g.upper())
        except OSError as e:
            log.warning("Could not load common essentials: %s", e)
            return 0
        return len(self._common_essentials)

    def gene_effect(self, gene_symbol: str) -> FunctionalEvidence:
        g = gene_symbol.upper()
        if g in self._cache:
            return self._cache[g]
        out = FunctionalEvidence(gene_symbol=gene_symbol)
        if self._common_essentials:
            out.is_common_essential = g in self._common_essentials
        if not self.csv_path:
            out.error = (
                "No local DepMap data configured. Download "
                "CRISPRGeneEffect.csv from depmap.org and pass its path."
            )
        self._cache[g] = out
        return out
