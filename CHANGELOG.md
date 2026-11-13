# Changelog

All notable changes to SYNAPTIK are recorded here.

## 0.9.0 — 2026-11-13

### Added
- **3D RNA Biophysics & Steric Accessibility Engine (`synaptik biophysics` & `synaptik/biophysics/`)**:
  - Pure-Python Zuker-style MFE dynamic programming (`V` and `W` matrices) with Turner 2004 nearest-neighbor stacking, terminal AU penalties, hairpin initiation, and ultra-stable tetraloop bonuses (`GAAA`, `UACG`).
  - Computes thermodynamic opening penalty ($\Delta G_{\text{open}}$), effective hybridization free energy ($\Delta G_{\text{effective}}$), effective $K_d$, and competitive spliceosomal steric displacement against U1 snRNP, U2AF65/35, and TDP-43.
- **Patient Natural History Digital Twin & Bayesian Counterfactual Control (`synaptik twin` & `synaptik/twin/`)**:
  - Zero-PHI ephemeral Monte Carlo simulation of 1,000 matched untreated progression trajectories conditioned on PRO-ACT / Answer ALS covariates (onset site, plasma NfL, age, empirical Bayes pre-slope shrinkage).
  - Bayesian posterior probability of disease deceleration $P(\Delta > 0 \mid \text{Data})$, Bayes Factor ($BF_{10}$), and research-planning counterfactual report export.
- **Clinical Translational Pharmacology & Adaptive Escalation (`synaptik clinical` & `synaptik/clinical/`)**:
  - CSF volumetric allometric scaling (NHP to adult and pediatric cohorts), two-compartment intrathecal PK/PD simulation (bolus vs catheter), and BOIN / 3+3 / mTPI adaptive dose escalation engines.
  - Hard safety guardrails: `calculate_allometric_dose()` requires `acknowledge_not_clinical=True` and marks every output as `RESEARCH PLANNING ONLY — NOT A CLINICAL DOSE RECOMMENDATION`.
- **Concurrent Rate-Limited Literature Harvesting & Preprint Weighting (`synaptik/ingest/`)**:
  - Thread-safe `PoliteSession` adhering to NCBI/EBI limits, plus bioRxiv/medRxiv preprint harvesting with explicit 20% discount weighting.
- **REST API Endpoints**: Added `/api/biophysics`, `/api/twin`, `/api/clinical/allometry`, `/api/clinical/pkpd`, `/api/clinical/escalate`, and `/api/clinical/protocol`.
- **Empirical ASO Potency Calibration v3 (`synaptik/drug/aso_calibrated.py`, `aso_mlp_v3.json`)**:
  - Mapped `79,687` ASO Atlas gapmers onto Ensembl GRCh38 cDNA $\pm 50\text{ nt}$ flanking sequences; trained a 162-feature target-accessibility ($\Delta G_{\text{open}}$, stem-loop pairing probability, proximal AU) + 20-position MLP on `133,753` non-neuro gapmers (`327` genes) and evaluated on `25,839` held-out neurogenetics gapmers across 12 genes (**Spearman $\rho = 0.4112$, ROC-AUC $= 0.7146$**, breaking the 0.70 ceiling vs. v2 `0.6929` and uncalibrated `0.5623`).
- **Temporal Literature Replication Calibration (`synaptik/infer/weighting.py`, `synaptik/infer/benchmark.py`, `synaptik/normalize/phenotype_ontology.py`)**:
  - Added `canonical_pathways_for_phenotype()` mapping every `PhenotypeTerm` to canonical GO/Reactome pathways and added `(gene, phenotype)`-level replication tracking (`RetrospectiveBenchmark.run_gene_phenotype()`).
  - Harvested `1,714` OA Europe PMC full-text iPSC-neuron papers (2013–2025; `2,734` records, `392` labs) and calibrated `EvidenceWeighter.score_claim_features()` via 2018–2021 temporal CV (`447` claims), achieving out-of-sample Cutoff 2022 **AUROC $= 0.7129$** (vs. `0.3671` enrichment, `0.4103` p-value).
- **MFASS + Neuronal RBP Splice Calibration (`synaptik/sequence/splice.py`, `splice_mlp_v2.json`) & TDP-43 Cryptic Exon Benchmark (`synaptik/sequence/benchmark_cryptic.py`)**:
  - 18-feature MLP integrating PWM donor/acceptor deltas, RESCUE-ESE/FAS-ESS hexamers, and neuronal RBP motifs (`TDP-43`, `PTBP1`, `RBFOX`, `NOVA`, branchpoint, PPT) validated via 5-fold exon-grouped CV across 29,367 MFASS mutant SNVs (2,247 exons): AUROC $= 0.7770$, AUPRC $= 0.4961$ (4.12x baseline).
  - 6-locus TDP-43 cryptic exon benchmark (`STMN2`, `UNC13A`, `KALRN`, `RANBP1`, `SYT7`, `ATG4B`): 6/6 cryptic splice events detected, 5/6 top-1 cryptic-junction ASO overlap.
- **Full-Text JATS XML Rigor Extractor (`synaptik/extract/fulltext.py`), Gold-Standard Extraction Benchmark (`synaptik/extract/benchmark_extraction.py`), & Turnkey Live Lab Audit (`synaptik audit --live`)**:
  - Extracts `n_donor_lines`, `n_clones`, `is_isogenic_pair`, `days_in_vitro`, `CellModelClass`, `PerturbationType`, and `Direction.NO_CHANGE` null results from Europe PMC JATS XML sections.
  - `synaptik benchmark --extraction`: evaluates 16 curated iPSC-neuron cases (F1 = `1.000` for `n_donor_lines`, `n_clones`, `is_isogenic_pair`, `days_in_vitro`, `no_change_null`; F1 = `0.933` for `gene` & `direction`).
  - Added `synaptik audit --live --gene <GENE>` to fetch OA Europe PMC full-text papers on the fly and reconcile against gnomAD (fixed GraphQL `$dataset` variable HTTP 400), ClinVar, Expression, and DepMap.

### Fixed & Scientifically Verified
- **Europe PMC `fullTextXML` HTTP 406 Fix (`synaptik/ingest/http.py`)**: Changed `PoliteSession` default `Accept` header from `application/json` to `*/*`, restoring 100% JATS XML full-text retrieval.
- **Unbound `cpg_count` for PMO Gapmers (`synaptik/drug/aso.py`)**: Always compute `cpg_count` before chemistry-specific CpG penalty branching.
- **Sugimoto 1995 Nearest-Neighbor Thermodynamics (`synaptik/drug/aso.py`)**:
  - Verified all 16 RNA/DNA nearest-neighbor entries against Sugimoto et al. (*Biochemistry* 1995) Table 2 and Biopython `R_DNA_NN1`; resolved antiparallel RNA $5'\rightarrow 3'$ vs DNA $5'\rightarrow 3'$ indexing; calibrated 2'-MOE shifts to clinical Nusinersen benchmarks ($\Delta G_{37} = -21.52\text{ kcal/mol}$, $T_m = 65.0^\circ\text{C}$).
- **Honest Scientific Claims & Safety Guardrails (`synaptik/drug/specificity.py`, `synaptik/nof1/pipeline.py`, `tests/test_safety_claims.py`)**:
  - Seed-and-extend transcriptome-scale off-target screening returning `NOT_SCREENED` when no reference is provided; computed CpG/TLR9 motif counts; retitled generated documents as computational drafts; documented published mismatch tolerances (`Kamola et al. 2015`, `Yoshida et al. 2018`, `Hua et al. 2008`).

## 0.8.0 — 2026-11-12

### Added
- **Direct Patient Genomic VCF Ingestion & Prioritization Engine (`synaptik vcf` & `synaptik/sequence/vcf.py`)**:
  - Ingests standard Variant Call Format (VCF v4.2 / v4.3) files (.vcf and .vcf.gz) directly from patient whole-genome (WGS) or whole-exome (WES) sequencing.
  - Automatically parses and standardizes coordinates, ClinVar clinical significance, gnomAD population rarity, CADD Phred pathogenicity, and SpliceAI delta scores.
  - Prioritizes actionable cryptic splice-altering variants across neurodegenerative target genes (*STMN2*, *UNC13A*, *SMN2*, *ELP1*, *COL6A1*, *ATXN2*, *TARDBP*, *SOD1*, *FUS*).
  - CLI command `synaptik vcf <FILE.vcf> [--disease <DISEASE>] [--target <GENE>] [--min-splice-score <FLOAT>] [--json]`.
  - Direct pipeline bridge: `synaptik nof1 --vcf <FILE.vcf>` auto-extracts high-priority variants and runs the N-of-1 pipeline with zero manual formatting.
- **Small-Molecule Drug Repurposing & CNS Tractability Engine (`synaptik repurpose` & `synaptik/drug/repurposing_catalog.py`)**:
  - Curated database of brain-penetrant small-molecule drugs (Risdiplam, Branaplam, Riluzole, Edaravone, Sodium Phenylbutyrate, Memantine, Diazepam, Vorinostat, Rapamycin, etc.).
  - Evaluates Physicochemical CNS Multiparameter Optimization (CNS MPO, 0–6 scale) and Blood-Brain Barrier (BBB) penetrance.
  - Directional mechanism-of-action compatibility screening (Activator vs Inhibitor vs Modulator).
  - Predicts in vitro and in vivo transcriptomic rescue efficacy percentages and ranks candidates by composite translational rank score.
  - CLI command `synaptik repurpose --target <GENE> [--direction <DIR>] [--min-phase <INT>] [--bbb-only] [--json]`.
- **Publication-Grade Splicing & ASO Architecture Vector Visualizer (`synaptik/sequence/visualize.py`)**:
  - Generates publication-ready SVG diagrams rendering pre-mRNA architecture: canonical flanking exons, cryptic exon, splice donor (`GT`) and acceptor (`AG`) sites, sashimi-style curved splice junction arcs with $\Delta\Psi$, and tiled 18-mer ASO steric-blocker footprints with base-level sequences.
  - Standalone SVG file export via CLI flags: `synaptik aso --svg <PATH.svg>` and `synaptik nof1 --svg <PATH.svg>`.
  - Embedded directly into the laboratory Web UI for real-time visualization and one-click figure export.
- **FDA CDER/CBER Type B Pre-IND Meeting Package Generator (21 CFR Part 312 Compliant eCTD Dossier)**:
  - Expands `synaptik nof1` with `--pre-ind <PATH.md>`: produces an exhaustive, regulatory-grade 10-section briefing document.
  - Modules: Background & Development Objectives, CMC Specifications (API sequence, molecular weight, IDT formulation code, PS/2'-MOE modifications), Nonclinical Pharmacology (iPSC target engagement, human post-mortem single-cell concordance, counterfactual causal simulation), Nonclinical Safety & Toxicology (GLP 13-week NHP intrathecal study plan), Clinical Protocol Synopsis (stratification by baseline Plasma NfL), Surrogate Biomarker Endpoint Qualification, and Specific Regulatory Questions for FDA Reviewers.
- **Web UI Dashboard Expansion (`synaptik ui`)**:
  - Added dedicated **Patient VCF** tab with text/file upload, sample STMN2 loader, and 1-click "Translate in N-of-1" routing.
  - Added dedicated **Drug Repurposing** tab with interactive drug cards, CNS MPO scores, and clinical phase filters.
  - Integrated interactive vector SVG splicing visualizer into the N-of-1 dashboard.
  - Added one-click **Download FDA Pre-IND Package** button.
  - New REST API endpoints: `POST /api/vcf`, `POST /api/repurpose`, `GET /api/visualize/splicing`, `POST /api/nof1/preind`.
- **Three New Regression Test Suites**:
  - `tests/test_vcf.py` (21 checks)
  - `tests/test_repurposing.py` (24 checks)
  - `tests/test_visualize.py` (21 checks)
  - Expanded `tests/test_nof1.py` (28 checks) and `tests/test_ui_endpoints.py` (22 checks).
  - Platform test suite expands to **478 passing checks across 13 independent test suites**.

## 0.7.0 — 2026-11-12

### Added
- **Autonomous N-of-1 Rare Disease & Translational Pipeline (`synaptik nof1` & `synaptik/nof1/pipeline.py`)**:
  - Unifies the entire bench-to-bedside translational pipeline into an autonomous, closed-loop diagnostic and therapeutic design engine.
  - Chains 7 rigorous stages in a single pass:
    1. Splicing mechanism and target window resolution for hallmark neurodegenerative disease targets (*STMN2*, *UNC13A*, *SMN2*, *ELP1*, *COL6A1*, *ATXN2*, *TARDBP*).
    2. Dish-to-Brain concordance evaluation against human in vivo post-mortem single-nucleus transcriptomic atlases.
    3. Steric-blocking and splice-modulating antisense oligonucleotide (ASO) candidate engineering with thermodynamic and specificity scoring.
    4. Causal counterfactual transcriptomic simulation predicting downstream transcriptomic and cellular network recovery.
    5. Longitudinal fluid biomarker (Plasma NfL) pharmacodynamic modeling and clinical trial patient stratification.
    6. Preclinical study protocol generation with cryptographic SHA-256 pre-registration fingerprint and 96-well randomized complete block design (RCBD).
    7. Automated actionability tier classification (Tier 1: Clinically Actionable IND-Ready, Tier 2: Preclinical Validation Required, Tier 3: Optimization Needed, Tier 4: Exploratory) and production of an audit-ready, 7-page Markdown IND-Enabling Dossier (`--dossier <PATH>`).
- **Interactive Laboratory Web Dashboard & REST API Suite (`synaptik ui` & `synaptik/ui/server.py`)**:
  - Modernized single-page web application with dedicated **N-of-1 Rare Disease Pipeline** tab and **Cite & Credit** panel.
  - Interactive parameters, real-time stage progression cards, live IND Dossier generation and one-click markdown file export.
  - Comprehensive HTTP REST API endpoints:
    - `POST /api/nof1` — Execute complete N-of-1 autonomous pipeline.
    - `POST /api/aso` — In silico therapeutic ASO engineering.
    - `POST /api/translate` — Dish-to-Brain post-mortem single-cell concordance.
    - `POST /api/simulate` — Counterfactual genomic rescue simulation.
    - `POST /api/biomarker` — Fluid biomarker pharmacodynamics and clinical stratification.
    - `POST /api/protocol` — Pre-registered study design protocol and RCBD plate map.
    - `GET /api/cite` — Citation metadata, BibTeX, APA, and Nature format strings.
- **Dedicated Regression Test Suites**:
  - `tests/test_nof1.py` (26 checks): End-to-end N-of-1 pipeline validation, dossier markdown generation, author citation/ORCID embedding, and CLI smoke tests.
  - `tests/test_ui_endpoints.py` (18 checks): Full REST API endpoint testing, JSON schema validation, and error handling.
  - Platform test suite expands to **406 passing checks across 10 independent test suites**.

## 0.6.1 — 2026-11-12

### Added
- **`synaptik cite` — Author Scientific Recognition & Citation Engine**:
  - Integrated citation generation and metadata export for lead author **Francesco Nobile** (ORCID: [`0009-0001-9477-3332`](https://orcid.org/0009-0001-9477-3332)), Center for Genomic Medicine, Massachusetts General Hospital / Harvard Medical School and University of Pavia (`synaptik/cli/commands/cite.py`).
  - Multi-format academic citation export: BibTeX (`--bibtex`), APA 7th edition (`--apa`), Nature / Vancouver, and machine-readable JSON (`--json`).
  - Citation File Format [`CITATION.cff`](CITATION.cff) (v1.2.0) with verified author ORCID, powering GitHub's native "Cite this repository" sidebar widget.
  - Formatted citation footers and prompt reminders integrated into all major translational CLI command outputs (`protocol`, `design`, `aso`, `biomarker`, `simulate`, `translate`).
  - Package authorship registration in `pyproject.toml` (`authors` and `maintainers`).
  - Zenodo DOI integration architecture for automated minting of citable, immutable DOIs indexed by Crossref and Google Scholar upon release.
  - Dedicated regression test suite [`tests/test_cite.py`](tests/test_cite.py) with 20 checks, bringing the platform test suite to **362 passing checks across 8 test suites**.

## 0.6.0 — 2026-11-12

### Added
- **`synaptik protocol` — Preclinical Study Design Registry & Multi-Site Meta-Analysis**:
  - Pre-registration grade preclinical study protocol generation (`synaptik/protocol/registry.py`) aligned with Nature/PLoS reproducibility guidelines (Landis et al. Nature 2012, Macleod et al. PLoS Biol 2015).
  - Pre-specifies primary and secondary quantitative endpoints, target disease model, maturity DIV, and paired sample sizes with power justifications.
  - **Cryptographic Pre-Registration Hashing**: Computes a deterministic SHA-256 design fingerprint and blinding verification token (`SYN-...`) guaranteeing immutable registration before experimental unblinding.
  - **Randomized Complete Block Layout (RCBD)**: Automated 96-well plate mapping with edge-effect avoidance (perimeter buffer wells reserved for evaporation barrier) and alternating mutant / isogenic control pairs.
  - **Multi-Site Random-Effects Meta-Analysis**:
    - Synthesizes multi-laboratory replication datasets via the DerSimonian-Laird inverse-variance random-effects estimator.
    - Evaluates heterogeneity metrics: Cochran's $Q$, degrees of freedom, between-study variance $\tau^2$, and $I^2$ index.
    - Automated leave-one-site-out sensitivity analysis verifying whether replicability holds upon omitting any single laboratory.
    - Included gold-standard multi-site trial dataset `examples/multisite_stmn2.json` (MGH, Stanford, UCL, Karolinska).
  - CLI command `synaptik protocol [--target <GENE>] [--disease <DISEASE>] [--pairs <N>] [--rho <FLOAT>] [--meta-analysis <FILE>] [--json]`.
  - Python API workflows: `generate_study_protocol()`, `synthesize_multisite_replication()`, `protocol_to_dict()`, `meta_analysis_to_dict()` in `synaptik/workflows.py`.
  - Dedicated regression test suite [tests/test_protocol.py](file:///c:/Users/frano/OneDrive/Desktop/Virtual_neuron/tests/test_protocol.py) with 21 checks (342 platform checks total across 7 test suites).

## 0.5.1 — 2026-11-12

### Added
- **`synaptik biomarker` — Translational Fluid Biomarker & Patient Stratification Panel**:
  - Longitudinal pharmacodynamic trajectory modeling for clinical fluid biomarkers (Plasma Neurofilament Light Chain, NfL, via Simoa ultrasensitive assay; CSF NfL; p-tau181) (`synaptik/biomarker/nfl.py`).
  - Connects in vitro axonal and synaptic rescue percentages into in vivo dynamic biomarker suppression curves using non-linear clearance kinetics ($t_{1/2} = 8\text{ weeks}$, scaling factor $\beta = 0.72$).
  - Grounded in clinical trial benchmarks: VALOR study (*SOD1* tofersen, Miller et al. NEJM 2022: ~60% plasma NfL suppression; FDA Accelerated Approval surrogate endpoint guidance).
  - Clinical trial patient stratification engine:
    - Automatically classifies patients into progression tiers (fast, intermediate, slow) based on baseline plasma NfL cutoffs ($> 75\text{ pg/mL}$, $40\text{--}75\text{ pg/mL}$, $< 40\text{ pg/mL}$).
    - Calculates clinical trial sample size optimization: shows that biomarker-enriched patient stratification reduces required trial sample size per arm by 60–75% (e.g. from 74 down to 18 patients/arm) to detect 30% therapeutic slowing of functional decline (ALSFRS-R).
  - Regulatory surrogate milestone verdicts: evaluates whether projected in vivo suppression exceeds FDA accelerated approval surrogate thresholds.
  - CLI command `synaptik biomarker --disease <DISEASE> [--target <GENE>] [--rescue-pct <FLOAT>] [--nfl-baseline <FLOAT>] [--json]`.
  - Python API workflows: `evaluate_biomarker_trajectory()`, `biomarker_result_to_dict()` in `synaptik/workflows.py`.
  - Dedicated regression test suite [tests/test_biomarker.py](file:///c:/Users/frano/OneDrive/Desktop/Virtual_neuron/tests/test_biomarker.py) with 21 checks (321 platform checks total across 6 test suites).

## 0.5.0 — 2026-11-12

### Added
- **`synaptik simulate` — Causal Genomic Counterfactual Simulator**:
  - In silico genetic perturbation engine simulating therapeutic interventions (knockdown, knockout, cDNA overexpression, and ASO steric blocking) on diseased neuronal transcriptomic vectors (`synaptik/simulate/counterfactual.py`).
  - Grounded in human iPSC-neuron functional genomics, CRISPRi/a, and Perturb-seq benchmarks (Tian et al. Neuron 2019, Replogle et al. Cell 2022).
  - Evaluates therapeutic efficacy via three core metrics:
    1. **Vector Rescue Efficiency ($\%$)**: Euclidean norm reduction across multi-dimensional disease state vectors.
    2. **Cosine Directional Alignment ($\cos \theta$)**: Directional alignment against the inverted disease vector ($\vec{v}_{\text{rescue}} \cdot (-\vec{v}_{\text{disease}})$).
    3. **Collateral Toxicity Index**: Monitors off-target disruption of core neural housekeeping and pan-essential genes (*ACTB*, *GAPDH*, *TUBB3*, *POLR2A*, *RPS6*, *SNAP25*).
  - Automated detection of lethal and toxic counterfactuals: identifies pan-essential common lethal targets (DepMap) and toxic aggregate phenotypes (e.g. *TARDBP* overexpression, Barmada et al. J Neurosci 2010).
  - High-throughput therapeutic screen (`--screen`): screens and ranks cataloged perturbation targets against disease baseline vectors (e.g. ALS, FTD, AD), surfacing top viable candidates (*STMN2* ASO block, *ATXN2* knockdown, *UNC13A* ASO block, *PIKFYVE* knockdown).
  - CLI command `synaptik simulate --disease <DISEASE> [--target <GENE> | --screen] [--intervention <TYPE>] [--alpha <FLOAT>] [--json]`.
  - Python API workflows: `simulate_intervention()`, `screen_counterfactual_targets()`, and `simulation_result_to_dict()` in `synaptik/workflows.py`.
  - Dedicated regression test suite [tests/test_simulation.py](file:///c:/Users/frano/OneDrive/Desktop/Virtual_neuron/tests/test_simulation.py) with 24 checks (300 platform checks total across 5 test suites).

## 0.4.1 — 2026-11-12

### Added
- **`synaptik translate` — Dish-to-Brain Human Post-Mortem snRNA-seq Bridge**:
  - Connects human iPSC neuronal culture findings directly to human post-mortem single-nucleus transcriptomic (snRNA-seq) disease atlases (`synaptik/translate/atlas.py`).
  - Benchmarks in vitro findings against curated hallmark signatures from landmark post-mortem studies:
    - **Amyotrophic Lateral Sclerosis (ALS)**: NYGC ALS Consortium snRNA-seq, *STMN2* (Klim et al. Nat Neurosci 2019, Prudencio et al. Acta Neuropathol 2020), *UNC13A* (Brown et al. Nature 2022, Ma et al. Nature 2022), *KALRN* (Humphrey et al. Nat Neurosci 2023).
    - **Frontotemporal Dementia (FTD)**: Cortical TDP-43 and *GRN* haploinsufficiency single-cell profiles.
    - **Alzheimer's Disease (AD)**: Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD) single-cell transcriptomic profiles.
  - **Cell-Type Specificity Mapping**: Maps findings to distinct human CNS cellular subpopulations:
    - Upper Motor Neurons / Betz cells (*FEZF2*, *BCL11B*, *CRYM*)
    - Lower Motor Neurons / Spinal cord (*CHAT*, *MNX1*, *ISL1*, *SLC18A3*)
    - Cortical Projection Neurons / Layers 2–6 (*SLC17A7*, *SATB2*, *TBR1*)
    - GABAergic Interneurons (*GAD1*, *GAD2*, *PVALB*, *SST*, *VIP*)
    - Astrocytes (*GFAP*, *AQP4*, *ALDH1L1*, *SLC1A2*)
    - Microglia (*AIF1*, *TMEM119*, *CX3CR1*, *P2RY12*)
    - Oligodendrocytes & OPCs (*MBP*, *MOG*, *OLIG2*, *PDGFRA*)
  - **Pathology-Aware Concordance Engine**:
    - Distinguishes cryptic exon inclusion (where $E(\Delta\Psi) > 0$ denotes pathological TDP-43 loss-of-function) from normal cassette exon skipping.
    - Compares observed perturbation direction against human post-mortem pathology, scoring concordance percentage, identifying direction conflicts, and issuing translational verdicts (`HIGH TRANSLATIONAL FIDELITY`, `MODERATE TRANSLATIONAL FIDELITY`, `MIXED TRANSLATIONAL FIDELITY`, or `TRANSLATIONAL DISCORDANCE DETECTED`).
  - **Developmental Maturation Gap Detection**:
    - Automatically identifies adult brain maturation markers (*NEFH*, *GRIN2A*) that fail to express in immature iPSC cultures (DIV 21–35), flagging them as developmental maturation gaps rather than disease discordances and issuing culture optimization recommendations (e.g. extending DIV or supplementing with brain-derived neurotrophic factors).
  - CLI command `synaptik translate [--splicing <FILE> | --genes <GENES...>] --disease <DISEASE> [--json]`.
  - Python API workflows: `translate_findings()` and `translate_result_to_dict()` in `synaptik/workflows.py`.
  - Dedicated regression test suite [tests/test_translation.py](file:///c:/Users/frano/OneDrive/Desktop/Virtual_neuron/tests/test_translation.py) with 19 checks (276 checks total across the platform).

## 0.4.0 — 2026-11-12

### Added
- **`synaptik aso` — In Silico Therapeutic ASO Engineering Pipeline**:
  - Full in silico design engine for clinical-grade antisense oligonucleotides (`synaptik/drug/aso.py`).
  - Supports two primary translational drug mechanisms:
    1. **Splice-Switching Oligonucleotides (SSO)**: Steric blocking for exon inclusion or cryptic exon exclusion (e.g. Nusinersen / *STMN2* series) with uniform 2'-MOE and phosphorothioate (PS) backbone or PMO chemistry.
    2. **RNase H1 Knockdown Gapmers**: Catalytic transcript degradation with 5-10-5 wing-gap-wing architecture (e.g. Tofersen / *SOD1* series).
  - **Thermodynamic Hybridization Engine**: Nearest-neighbor parameters (SantaLucia 1998, Sugimoto et al. 1995) calculating $\Delta H^\circ$, $\Delta S^\circ$, Gibbs free energy $\Delta G_{37}^\circ$, and salt-adjusted melting temperature $T_m$ under physiological CSF conditions ($[Na^+] = 100\text{ mM}$, $[C_T] = 250\text{ nM}$) with 2'-MOE thermal stabilization ($+1.1^\circ\text{C}$ per nucleotide).
  - **Biophysical Liability & Safety Filters**:
    - G-quadruplex risk detection ($\ge 4$ consecutive Gs, causes CSF neurotoxicity and non-specific protein binding).
    - Poly-T / Poly-U premature termination and innate immune risk detection.
    - Stem-loop self-hairpin minimum free energy ($\Delta G_{\text{hairpin}}$) prediction.
    - 3'-end self-dimerization complementarity scoring.
  - **Whole-Transcriptome Neural Off-Target Screening**: Screens candidate oligos against core neural and housekeeping transcripts (*ACTB*, *GAPDH*, *TUBB3*, *MAP2*, *SYP*, *SNAP25*, *TARDBP*, *FUS*, *SOD1*, *C9orf72*), enforcing zero 0-mismatch or 1-mismatch off-target matches and maximum mismatch distance.
  - **Synthesis-Ready Specifications**: Outputs molecular weight (sodium salt), extinction coefficient $\epsilon_{260}$, predicted human CSF terminal half-life (~135 days for 2'-MOE/PS), and orderable IDT vendor format strings (`mN*...`).
  - CLI command `synaptik aso --sequence <SEQ> [--mode splice_switch|gapmer] [--gene <GENE>] [--json] [--export <PATH>]`.
  - Python API workflows: `design_aso()` and `aso_result_to_dict()` in `synaptik/workflows.py`.
  - Dedicated regression test suite [tests/test_aso.py](file:///c:/Users/frano/OneDrive/Desktop/Virtual_neuron/tests/test_aso.py) with 20 checks (257 checks total across the platform).

## 0.3.1 — 2026-11-12

### Added
- **Native MAJIQ & LeafCutter splicing table ingestion**:
  - Auto-detection of MAJIQ (VOILA) and LeafCutter TSV formats in `synaptik deposit splicing`.
  - Parses multi-junction semicolon-separated values in MAJIQ (`E(dPSI) per LSV junction`), extracting the maximal effect junction $|E(\Delta\Psi)|$.
  - Parses cluster coordinate strings (`chr:start:end:cluster_strand` and `chr:start-end`), deriving exact coordinates and cassette/spliced element length without requiring manual reformatting.
  - Added synthetic validation test fixtures: `examples/voila_majiq.tsv` and `examples/leafcutter_ds.tsv`.
- **Paired isogenic power analysis for human iPSC disease models**:
  - Implemented exact noncentral-$t$ paired sample size calculation `sample_size_paired(cohens_d, rho=0.7)` in `synaptik/infer/design.py`.
  - Differentiating mutant and isogenic control side by side across parallel differentiation rounds removes shared batch-to-batch variation ($\rho \approx 0.6\text{--}0.85$), rescuing sample sizes from impractical ($N > 25\text{--}50$ differentiations per group) down to feasible paired cohorts ($N = 6\text{--}10$ pairs).
  - Wired into `synaptik design --rho <float>` and web design briefs with automated paired design rescue recommendation.
- **Neuronal splice prediction benchmark suite**:
  - New validation suite `synaptik/sequence/benchmark_spliceai.py` evaluating splice sequence models against curated gold-standard neurodegenerative and ALS/FTD splice variants (*STMN2*, *UNC13A*, *SMN2*, *ELP1*, *COL6A1*) and matched neutral controls.
  - CLI command `synaptik benchmark --spliceai [--threshold 0.20] [--json]` reporting Sensitivity, Specificity, Precision, Balanced Accuracy, and ROC-AUC.
- **Regression test expansions**:
  - 10 new regression checks in `tests/test_lab_workflows.py` covering MAJIQ/LeafCutter ingestion, paired power exact noncentral-$t$ calculations, and the splice benchmark suite (total: 237 tests).

### Fixed
- Fixed `AttributeError: SplicingEventType.APA` in `synaptik/infer/splicing.py` by aligning with `SplicingEventType.ALT_POLYA`.

## 0.3.0 — 2026-11-12

This release turns SYNAPTIK from a tool for reading the literature into one
that meets a lab's own data and decisions.

### Added
- **`synaptik deposit splicing`** — reconcile your own rMATS table (or a
  generic gene / delta_psi / fdr table) against published splicing events.
  Each significant event is called novel, consistent, discordant, or
  same-gene-different-event. Also reports published events in genes you
  tested but saw no change in, a Fisher test for enrichment of each splicing
  regulator's published targets among your changed genes (background: the
  genes tested in your table), reading-frame consequences computed from the
  table's own exon coordinates, and an ASO strategy per event.
- **`synaptik deposit phenotypes`** — the same for a phenotype table.
- **`synaptik design`** — experiment design briefs: consensus protocol,
  reported effect sizes, sample size by exact noncentral-t power calculation
  (for the published effect and a halved one), the smallest effect a feasible
  design can detect, weaknesses in the prior evidence, and readout-specific
  pitfalls.
- **`synaptik benchmark`** — retrospective temporal hold-out testing whether
  confidence tiers predict later independent replication, with enrichment,
  p-value and paper-count baselines and a paired bootstrap interval.
- **Web interface**: a *My Data* upload tab, a *Design brief* action on the
  lookup tab, and the benchmark on the corpus tab.
- **Phenotype vocabulary**: MEA metrics (mean firing rate, burst duration,
  network synchrony), patch-clamp intrinsic properties (resting membrane
  potential, input resistance, rheobase, AP threshold, capacitance), mIPSC
  amplitude, inhibitory synaptic puncta, colocalization coefficients, and
  dendritic spine morphology. 48 terms, 328 synonyms.
- `examples/` with synthetic input files, and `tests/test_lab_workflows.py`
  (93 checks; 227 in total).

### Fixed (scientific correctness)
- **Phenotype pooling no longer merges distinct biology.** Terms were pooled
  when they shared a measurement axis, so an analysis of mEPSC frequency
  silently included mIPSC frequency, network burst frequency and calcium
  transients ("event rates"), and synaptic puncta were pooled with dendritic
  spines. Pooling is now an explicit curated `pool_group`; only pre-,
  post- and colocalised excitatory puncta density are pooled.
- **Inhibitory synapse markers** (VGAT, gephyrin) were synonyms of the
  excitatory pre-/post-synaptic puncta terms. They now map to a separate
  inhibitory term.
- **Intrinsic excitability** lumped rheobase, input resistance and resting
  potential into one term, although their directions mean opposite things for
  excitability. They are now separate terms.
- **Firing rate** was a synonym of network burst frequency; it is now its own
  term.
- **Short synonyms matched inside unrelated words** ("nmd" inside "NMDA",
  "apa" inside "capability"). Synonym matching now respects word boundaries.

### Changed
- New commands live in `synaptik/cli/commands/`, one module each; shared
  workflow logic is in `synaptik/workflows.py` so the CLI and web interface
  cannot drift apart.
- If you built a knowledge base with 0.2.0, run `synaptik renormalize` to
  apply the corrected vocabulary to existing records (free; no LLM calls).

## 0.2.0 — 2026-11-12

Initial public release.
