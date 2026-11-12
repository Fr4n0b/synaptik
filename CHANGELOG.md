# Changelog

All notable changes to SYNAPTIK are recorded here.

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
  network synchrony; including synonyms contributed by @HarshRajSinghania), patch-clamp intrinsic properties (resting membrane
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
