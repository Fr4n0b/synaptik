# CLAUDE.md — SYNAPTIK Project Context & Engineering Guide

> **Note for AI Assistant / Claude Code:**  
> This document provides the complete context, architectural blueprint, development history, test commands, and roadmap for **SYNAPTIK**. Read this first before making code modifications, proposing refactors, or adding new evidence adapters.

---

## 1. Project Overview & Core Philosophy

- **Project Name:** SYNAPTIK (`v0.2.0 Alpha`)
- **Repository:** [`https://github.com/Fr4n0b/synaptik`](https://github.com/Fr4n0b/synaptik)
- **Documentation Website:** [`https://fr4n0b.github.io/synaptik/`](https://fr4n0b.github.io/synaptik/)
- **License:** MIT
- **Target Biological Domain:** Translational neurobiology and cellular physiology, specifically **patient-derived human iPSC-differentiated neurons** (NGN2 induction, dual-SMAD protocols).
- **Maintainer:** Francesco Nobile (MD & Biomedical Engineering Candidate, University of Pavia & IUSS Pavia; Armenise-Harvard Fellow at Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School).

### The Fundamental Problem SYNAPTIK Solves
When bench researchers search PubMed for a cellular readout (e.g., reduced SYN1/PSD95 puncta density, altered mEPSC frequency), standard tools rank associations by p-value or enrichment. In practice, high-scoring associations often trace back to a **single research group using one assay**, with no independent replication or published null results across other labs.

SYNAPTIK changes the unit of analysis from the *paper* to the **experiment record** (`paper × cell model × maturity [DIV] × perturbation × assay × outcome`). It reports **evidence structure**:
1. Independent-lab count (research group clustering).
2. Orthogonal assay diversity.
3. Leave-one-lab-out statistical robustness.
4. Direction consistency (detecting conflicting up/down findings).
5. Cross-layer reconciliation against population genetics (gnomAD constraint, ClinVar, DepMap) and sequence models (SpliceAI) to highlight contradictions (e.g., a published haploinsufficiency phenotype in a gene that tolerates loss-of-function in healthy cohorts).

---

## 2. Codebase Architecture & File Tree

The codebase contains ~15,500 lines of modular Python with strict scientific contracts:

```
synaptik/
├── ingest/                # Harvesters & API contract validation
│   ├── europepmc.py       # Full-text XML and abstract harvester
│   ├── pubmed.py          # NCBI E-utilities harvester
│   ├── biorxiv.py         # Preprint details harvester
│   ├── contracts.py       # Proactive API schema drift detector (offline & live)
│   ├── http.py            # HTTP client utilities and rate limiting
│   └── world.py           # Multi-source coordinator
├── extract/               # Extraction layer
│   └── llm.py             # LLM extraction requiring verbatim quotes from source text
├── normalize/             # Normalization & entity resolution
│   ├── phenotype_ontology.py # Controlled vocabulary for synaptic/neuronal phenotypes
│   ├── matcher.py         # Fuzzy matching & synonym resolution
│   └── labs.py            # Author affiliation parsing and heuristic research-group clustering
├── infer/                 # Statistical inference & evidence scoring
│   ├── association.py     # Background-conditioned Fisher exact test, FDR, leave-one-out
│   ├── lookup.py          # Graded experiment lookup (exact / close / related / blocked)
│   ├── splicing.py        # Alternative splicing inference and dPSI direction reconciliation
│   └── weighting.py       # 6-component evidence weighting (journal, design, recency, etc.)
├── evidence/              # Multi-class reconciliation
│   └── reconcile.py       # Cross-layer reconciliation (literature vs sequence vs genetics)
├── sequence/              # Sequence-level models
│   ├── splice.py          # SpliceAI adapter & PWM-based splice site prediction
│   └── reconcile_splice.py# Sequence-to-literature reconciliation
├── omics/                 # Functional & population genetics evidence
│   └── evidence_sources.py# gnomAD constraint (pLI, LOEUF), ClinVar pathogenicity, DepMap
├── drug/                  # Therapeutic candidate routing
│   └── repurposing.py     # ChEMBL / Open Targets with mechanism-direction filtering
├── cli/                   # Command-line interface
│   └── main.py            # Central CLI entrypoint (1,050 lines; target for modularization)
└── ui/                    # Flask web interface
    ├── server.py          # Local Flask backend (serves at http://127.0.0.1:8420)
    └── static/index.html  # Dark-theme "laboratory instrument" frontend
```

### Supporting Infrastructure
- `tests/test_inference.py`: 114 regression tests encoding biological correctness constraints.
- `tests/synthetic.py` & `demo_data.py`: Synthetic test corpus mimicking real literature traps.
- `notebooks/quickstart.ipynb`: Interactive 1-click Google Colab notebook.
- `docs/`: Full documentation (`RUNBOOK.md`, `DESIGN.md`, `DEPLOYMENT.md`, `INTERPRETING.md`, `index.html`).
- `.github/workflows/tests.yml`: Matrix CI on GitHub Actions across Python 3.10, 3.11, and 3.12.

---

## 3. Current Project State & Accomplishments

As of release `v0.2.0 (Alpha)`:
1. **GitHub Repository:** Publicly live at `Fr4n0b/synaptik` with 16 scientific SEO topics.
2. **CI / Automation:** 100% Green on GitHub Actions across Python 3.10, 3.11, and 3.12.
3. **Community Health:** **100% Score** on GitHub Community Standards (Code of Conduct v2.1, Security Policy, Contributing Guide, Pull Request Template, MIT License).
4. **Distribution:** Built wheels (`.whl`) and source archives (`.tar.gz`) attached to the `v0.2.0` release.
5. **Interactive Docs:** Live GitHub Pages site at `https://fr4n0b.github.io/synaptik/` with interactive Google Colab badge.
6. **Community PRs:** Active upstream Pull Requests on `Awesome-Bioinformatics` (PR #192), `Awesome-Neuroscience` (PR #59), and `Awesome-Genomic-Skills` (PR #14).
7. **Active Issues:** Seeded with newcomer contributor tags (`good first issue`, `help wanted`).

---

## 4. Development Rules & Scientific Conventions

When editing or extending SYNAPTIK, follow these core disciplines:

### 1. Test Execution
Before committing any changes, run the full validation suite:
```bash
./run_tests.sh
# or: bash run_tests.sh
```
`run_tests.sh` executes 4 distinct stages:
1. Unit & Regression Suite (`python3 -m tests.test_inference`)
2. API Contract Self-Test (`python3 -m synaptik.ingest.contracts`)
3. CLI Smoke Test on temporary SQLite database (`demo`, `phenotype`, `splicing`, `lookup`, `compare`, `report`)
4. CLI JSON Output Validity (`--json` output must strictly parse as valid JSON)

### 2. Scientific Properties Over Mere Code Paths
Tests must assert **scientific and statistical correctness properties**:
- A single-lab association with high enrichment and small p-value **must not** outrank an independently replicated association across multiple labs.
- Absence of evidence is **not** evidence of absence: missing data must remain `None`, never silently coerced to `0` or `"no effect"`.
- If an upstream API removes or renames an essential field (e.g. ChEMBL `action_type`), contract validation must fail immediately rather than silently writing empty records.

### 3. Syntax & Compatibility
- Maintain strict compatibility across **Python 3.10, 3.11, and 3.12**.
- Avoid nested quotes in f-strings like `f'{d['key']}'` which break in Python 3.10 and 3.11.

---

## 5. Technical Roadmap for Claude to Execute

These are the primary milestones where Claude Code can provide substantial improvements:

### Milestone 1: CLI Architecture Modularization
- **File:** `synaptik/cli/main.py` (currently ~1,050 lines in a single file).
- **Task:** Split into dedicated sub-command modules under `synaptik/cli/commands/`:
  - `demo.py`, `phenotype.py`, `lookup.py`, `gene.py`, `variant.py`, `compare.py`, `report.py`, `ui.py`.
- **Requirement:** Ensure all CLI commands and JSON output flags continue to pass `run_tests.sh`.

### Milestone 2: Async & Batched Harvesting
- **Files:** `synaptik/ingest/europepmc.py`, `pubmed.py`, `http.py`.
- **Task:** Migrate sequential HTTP GET calls to async batched requests using `httpx.AsyncClient` or `asyncio.gather`.
- **Requirement:** Respect API rate limits (NCBI max 3 req/sec unauthenticated, 10 req/sec with API key; Europe PMC concurrency guidelines).

### Milestone 3: SpliceAI Adapter Benchmarking & Pangolin Integration
- **Files:** `synaptik/sequence/splice.py`, `reconcile_splice.py`.
- **Task 3A:** Add a benchmark dataset of published, experimentally validated splice-altering variants and test the SpliceAI adapter recovery rate.
- **Task 3B:** Introduce a Pangolin adapter for tissue-aware alternative splicing prediction, prioritizing brain and neuronal transcripts.

### Milestone 4: Phenotype Ontology Expansion
- **File:** `synaptik/normalize/phenotype_ontology.py`.
- **Task:** Add synonym mappings for:
  - Microelectrode array (MEA) metrics: Mean firing rate, burst duration, synchrony index, network burst frequency.
  - Patch-clamp electrophysiology: Resting membrane potential, input resistance, action potential rheobase.
  - High-content confocal imaging: Synaptic colocalization coefficients (Mander's/Pearson's), dendritic spine morphometry (mushroom, stubby, thin).

### Milestone 5: Web UI Data Deposit & User Tables Ingestion
- **Files:** `synaptik/ui/server.py`, `synaptik/ui/static/index.html`.
- **Task:** Allow wet-lab researchers to drag-and-drop user-generated differential splicing tables (rMATS `SE.MATS.JC.txt`, MAJIQ `.voila.tsv`) or puncta quantification CSVs to immediately reconcile their lab's unpublished results against the global literature corpus.
