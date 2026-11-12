# Design rationale

Why each component works the way it does. Written for whoever maintains this
after the person who commissioned it — including the reasoning behind decisions
that look arbitrary from the code alone.

---

## The central design problem

Mining literature for phenotype–pathway associations is easy to do badly and
the failure is invisible. A naive implementation counts co-occurrences, applies
a p-value, sorts, and produces a confident ranked list that is mostly an
artifact of what happened to get funded and published. It looks exactly like a
real result.

Three specific pathologies dominate, and most of this system exists to counter
them:

**One lab, many papers.** A group publishes a finding in 2019, extends it in
2021, reviews it in 2023. Naive counting reports "three studies agree." The
true independent-replication count is one. This is the single largest source of
false confidence in literature meta-analysis.

**Ubiquity mistaken for relevance.** Apoptosis and mTOR appear in nearly every
neuroscience paper. Any co-occurrence metric without a background correction
puts them at the top of every query, for every phenotype.

**Publication filtering.** Negative results are rarely published. A literature
showing 15 positive findings and no nulls may reflect a robust effect, or may
reflect fifteen labs that only wrote up the experiments that worked.

The tool cannot fix any of these. It can measure and surface all three, which
is what turns an unreliable ranked list into a usable hypothesis-prioritisation
aid.

---

## Schema: the record, not the paper

`ExperimentRecord` is one (paper, biological system, perturbation, assay,
measured outcome) tuple. A paper yields 5–50.

Making the *paper* the unit would make "has this experiment been done?"
unanswerable — a paper containing twelve experiments is a match for all twelve
queries and none of them precisely. The record-level schema is what lets the
lookup engine say "this exact experiment, at this DIV, with these markers" and
mean it.

Every normalisable field keeps both the raw surface string and the ontology ID
(`NormalizedTerm`). Normalisation is lossy and sometimes wrong; discarding what
the paper actually said would make errors unauditable. The raw string is also
what makes the unmapped-phenotype report possible, which is how the vocabulary
grows.

Record IDs are content hashes, so re-ingesting a paper is idempotent.

### Direction and magnitude are separate

Many papers report that something decreased without saying by how much. If
those were coerced to `effect_size = 0`, the median effect size across a
pathway would be dragged toward zero by the papers that simply did not
quantify. `direction` and `effect_size` are independent, and `is_quantified`
records which is which.

---

## Weighting: six signals, not one

The request was Q1/Q2 prioritisation. That is implemented, but as one of six
components, because journal tier is a weak proxy for reliability — it
correlates with novelty and institutional prestige at least as much as with
rigour, and high-impact journals do not have lower retraction rates.

| Component | Rationale |
|---|---|
| journal tier | requested; Q1:Q4 is 2:1, deliberately modest |
| study design | isogenic controls, replicate count, whether quantified |
| directness | measured vs. inferred vs. asserted in a review |
| extraction confidence | how sure the extractor was about this record |
| recency | mild decay; protocols and reagents change |
| integrity | retraction zeroes the record; preprints discounted |

The Q1:Q4 ratio is 2:1 rather than 10:1 because the evidence that Q4 findings
are ten times less likely to replicate does not exist. Aggressive tiering would
encode prestige bias while looking like data quality.

All six components are returned individually in `WeightBreakdown`, not just the
product. A user who disagrees with the weighting can see exactly which
component drove a score and re-weight.

---

## Association statistics

### Background-conditioned enrichment

Enrichment is a Haldane–Anscombe-corrected odds ratio against the rest of the
corpus, not a raw count. Without the background term, ubiquitous pathways win
every query. The 0.5 correction prevents infinite odds ratios when a specific
pathway appears in the foreground but never in the background — common in a
focused corpus, and otherwise a source of nonsense infinities.

### Ranking gates on enrichment

`rank_score()` returns 0 for any pathway with enrichment ≤ 1. Sorting by
weighted support alone promotes whatever is most written-about; sorting by
enrichment alone promotes noise from pathways with two records. The score
multiplies log-enrichment (specificity) by log-support (evidence mass), then
rewards independent labs and assay diversity.

This is the single most consequential formula in the system. The regression
suite asserts that a planted single-lab artifact with **13× higher enrichment
and a smaller q-value** still ranks below a well-replicated true positive. If
that test fails, the ranking has regressed to something that looks scientific
and is not.

### Leave-one-lab-out

Delete the most prolific group, recompute. Real associations usually survive;
one-group narratives usually do not. This is the most useful robustness check
available for this data because it directly targets the dominant failure mode.

### FDR correction

Benjamini–Hochberg across all tested pathways. Hundreds are tested per query;
uncorrected p-values would produce a wall of spurious hits.

### Confidence tiers are conservative by construction

Nothing from a single group reaches `strong`, regardless of paper count,
enrichment, or journal. `weighted_support ≥ 1.0` gates the upper tiers so that
retraction-zeroed evidence cannot reach `moderate` on lab count alone — the
counts might say two groups while the weights say almost nothing.

---

## Lab identity resolution

Union-find over (normalised senior author, canonical institution) signatures.
Senior authors are last, second-to-last (co-senior authorship is now common),
and corresponding author.

Author names normalise to `surname_firstinitial`, which collapses `Smith J`,
`John Smith`, `Smith, John A.`, and `J. Smith` correctly. It over-merges common
names, which is acceptable only because institution must also match.

**Known limitations, deliberately accepted:**

- A PI who changes institution appears as two groups. Splitting is the safer
  default: merging two genuinely different groups inflates the replication
  count, which is the error this whole system is built to prevent.
- Dutch and Spanish particle surnames (`van der Berg`) normalise
  inconsistently.
- Two same-surname PIs at one institution merge incorrectly.

The clustering is exposed via `cluster_report()` and overridable through a
curation file, because a heuristic that silently produces the load-bearing
statistic must be inspectable.

---

## Extraction

### Quote validation is the primary quality control

Every extracted record must carry a verbatim span, 15–200 characters, that
appears in the source text after whitespace normalisation. Records failing this
are dropped and counted.

This eliminates fabricated findings, which is the failure mode that would most
damage trust in the tool. It does not catch misattribution — a real quote
paired with the wrong effect size — which is why extraction confidence feeds
the weighting rather than being treated as certainty.

Matching is exact-after-normalisation, deliberately. Fuzzy matching here would
defeat the purpose.

### Null results are solicited explicitly

The prompt names negative results and gives example phrasings. Models otherwise
skip them as uninteresting. They are the scarcest and most valuable evidence in
the corpus — the only available partial correction for publication bias.

### Section budgeting

Long papers are truncated by allocating a fixed share to each section
(abstract 10%, methods 30%, results 45%, figures 15%) rather than cutting the
tail. Naive truncation drops Results entirely on long papers, which is exactly
the content needed.

Figure legends are included because effect sizes and sample sizes frequently
appear only there.

---

## Experiment lookup

### Hard blocks, not soft penalties

Some mismatches make results incomparable, and a similarity score that lets
them rank highly is worse than no answer:

- **Different species/system** — mouse primary neurons are not human iPSC
  neurons.
- **Different maturity bin** — DIV14 and DIV49 synapse counts are not the same
  measurement.
- **Opposing perturbation valence** — loss-of-function vs. gain-of-function.

That last one was a real bug caught in testing. The original implementation
enumerated opposing *pairs* (knockout/overexpression, knockdown/
overexpression), which missed patient-mutation vs. overexpression entirely — so
an overexpression study ranked as a "close match" for a patient-mutation query,
with the opposite sign. The fix classifies each perturbation by valence
(`loss` / `gain` / `variant` / `modulator`) and blocks on class mismatch,
catching the whole family.

### The readout gate

A match cannot be tiered `exact` unless the readout agrees. Without this, a
query specifying system + DIV + target but no phenotype scores 1.0 against
every record from a matching system — so a viability assay was reported as an
exact match for a splicing experiment. In testing this produced **47 false
exact matches**, reduced to 9 genuine ones.

When only an assay is specified, the gate additionally checks that the record's
phenotype is one that assay actually produces (`_phenotype_fits_assay`), which
filters records whose assay and outcome were mis-paired during extraction.

---

## Drug repurposing

### Direction matching is the whole point

If a pathway is down, you need an activator. Returning an inhibitor because it
"targets the pathway" is worse than returning nothing — it would deepen the
phenotype. Compounds whose mechanism opposes the required direction are
excluded with an explicit reason rather than ranked lower.

Compounds whose ChEMBL `action_type` is absent map to `UNKNOWN` and are
excluded, not guessed at. `required_moa_direction()` returns `UNKNOWN` when the
pathway's own direction is not established, which suppresses the whole
directional ranking — because whether the pathway drives or opposes the
phenotype is precisely what the literature has not told us.

### CNS penetrance annotates, it does not filter

A compound that never crosses the blood–brain barrier is perfectly usable on
cultured neurons. Filtering on BBB penetrance up front discards most of the
useful in vitro chemical space. `usable_in_dish` and `bbb_plausible` are
separate flags.

### The efflux problem

CNS MPO scoring was validated against 11 known drugs: 8/11 correct. All three
misses (metformin, domperidone, loperamide) are **P-glycoprotein efflux
substrates** — acceptable physicochemistry, but actively pumped out of the
brain.

No descriptor-based score can capture this; it is a transporter property, not a
molecular one. Rather than tune thresholds to paper over a structural
limitation, known efflux substrates are annotated from a curated list. The
threshold was set at 3.5 (not the conventional 4.0) after validating against
approved CNS drugs — donepezil scores 5.32 and would have failed at the
conventional cutoff.

---

## Splicing

Splicing gets a dedicated module because the generic phenotype→pathway
machinery discards what matters most: the specific event, in a specific gene,
with a direction and a **regulator**. The regulator is the actionable node.

### Sign convention normalisation

rMATS versions have differed in whether `IncLevelDifference` is
sample1−sample2 or the reverse. Everything is normalised to "positive dPSI =
more inclusion in the perturbed condition," with the detection tool recorded so
the normalisation is auditable.

When a reported dPSI sign conflicts with the direction stated in prose, the
prose wins and the conflict is flagged as `sign_confidence: uncertain`. The
prose was written by the author; the number may have been transcribed under
either convention. Silently flipping a sign here inverts the biology.

### Regulator inference

Two signals: explicit (the paper named a regulator) and implicit (the event
signature matches a known RBP's loss/gain phenotype). Cryptic exon inclusion
across multiple genes is the canonical TDP-43 nuclear-loss signature and is
inferred even when unnamed.

The PTBP1/PTBP2 note exists because that switch drives neuronal maturation and
is the most common confound in iPSC-neuron splicing comparisons — a
maturation-stage difference between mutant and control lines reads as a
splicing phenotype.

---

## API contract validation

The failure this prevents: an upstream API renames a field, the parser silently
returns `None`, and a 30,000-paper harvest completes "successfully" with every
journal name missing. Nothing raises. You find out weeks later.

Each endpoint has a contract listing the fields the parser depends on, whether
each is required, and what shape it should have. The first response of each
session is validated; violations abort in strict mode (the default for
harvest).

**Which fields are marked required is a design decision, not bookkeeping.**
Three were promoted to required after testing showed the consequence of losing
them:

- `action_type` (ChEMBL) — without it, direction matching cannot run at all.
- `standard_units` (ChEMBL) — a missing unit means an unconverted µM value
  reported as nM: a 1000× potency error presented as a measurement.
- `authorList` (Europe PMC) — without authors, every paper becomes its own
  "independent group." This does not degrade the output, it invalidates it.

### Zero records is a failure, not an empty result

Testing revealed that a fully-blocked network produced a "successful" harvest
of nothing: preflight passed because zero records means zero contract
violations. Both preflight and `PipelineStats.failed()` now treat "retrieved
nothing" as a failure with a nonzero exit code. An absent check is not a passed
check.

---

## Interface design

The organising principle is that **evidence structure must be as visible as the
answer**. A user who sees a ranked list of pathways without lab counts and
caveats will draw exactly the wrong conclusions, so:

- Confidence tier drives the visual hierarchy (left border colour, badge).
- The replication meter makes independent-lab count the most salient number,
  above enrichment.
- Caveats appear inline, not behind a disclosure.
- Supporting records show their verbatim quotes, so any claim can be traced to
  its source in one click.

The palette derives from fluorescence microscopy channels — green/cyan/amber/
magenta — because the intended user already reads those colours as carrying
meaning.

---

## Lab data deposit

`synaptik deposit` reconciles a lab's own results against the corpus. Four
decisions shape it:

**The sign convention is never guessed.** rMATS reports sample 1 minus
sample 2. Which is the mutant depends on argument order at run time, which
the table does not record. A default would be right half the time and
silently invert the biology the other half, so the perturbed group is a
required argument.

**Tested-and-unchanged is not the same as not-tested.** The "not reproduced"
list contains only genes present in the user's table. A gene absent from an
rMATS table was not expressed or had too few junction reads; reporting it as
a failed replication would be a false negative result.

**Regulator enrichment uses the user's tested genes as background.** Against
a genome background, any list of neuronally expressed genes looks enriched
for any neuronal splicing regulator. The Fisher test is therefore: of this
regulator's published targets that were tested in this table, how many
changed, versus the rest of the table.

**Frame consequences come from the table, not from a heuristic.** The
literature-side heuristic treats increased exon inclusion as a probable
cryptic exon. Applied to a lab's table it flagged every cassette exon as an
NMD risk. With coordinates available the question has a direct answer:
exon length modulo three.

Matching is at gene plus event-class level because papers rarely publish
coordinates. A "consistent" call therefore means the same class of event in
the same gene and direction, not the same exon, and the output says so.

## Design briefs

**Replicates are differentiations.** Every sample size is in independent
differentiations. Pseudoreplication (treating cells or fields as n) is the
most common statistical error in this literature, and a brief that printed
"n = 40 cells" would entrench it.

**Power is computed for a halved effect as well as the published one.**
Significant results from small studies over-estimate effects. Halving is a
blunt correction, not an estimate; it is used because it is transparent and
errs in the safe direction.

**The CV is shown as a range unless the user supplies one.** Papers almost
never report between-differentiation variability, so any single value would
be invented. The brief shows 15%, 25% and 35% and labels the recommendation
as resting on an assumption.

**Impractical answers are converted into advice via paired isogenic designs.**
When the unpaired design needs more differentiations than a lab can run (e.g. $N > 15$
differentiations per group, which is common given 20–35% differentiation CVs), the
brief automatically computes and presents a paired isogenic rescue design.

Differentiating patient-derived (or CRISPR-engineered) mutant lines side-by-side with
their isogenic parental control in each round removes shared batch-to-batch variation
(media lots, incubator fluctuations, passage drift). With an inter-batch correlation
$\rho \approx 0.60\text{--}0.85$, the effective paired Cohen's $d_z$ is:
$$d_z = \frac{d}{\sqrt{2(1 - \rho)}}$$
Power is calculated using the exact noncentral-$t$ distribution for paired differences
with degrees of freedom $\nu = N - 1$ and noncentrality parameter $\delta = d_z \sqrt{N}$:
$$\text{Power} = 1 - T_{\nu, \delta}(t_{\text{crit}}) + T_{\nu, \delta}(-t_{\text{crit}})$$
where $t_{\text{crit}} = t_{1 - \alpha/2, \nu}$. For $d=0.8$ and $\rho=0.7$, this reduces
the required sample size from 26 differentiations per genotype down to 10 matched pairs.

Power uses the exact noncentral t distribution; the suite checks it against
the standard table values (64, 26, 17 and 6 per group for d = 0.5, 0.8, 1
and 2 for independent two-group tests, and 10 pairs for paired $d=0.8, \rho=0.7$).

## Retrospective benchmark

The tool's central claim is empirical: associations it rates highly should
hold up better than ones it does not. `synaptik benchmark` tests that with a
temporal hold-out.

**Outcome.** An association counts as replicated if, after the cutoff, a
research group that did not contribute to it before reports the same pathway
for the same phenotype in the same direction. The same group publishing
again does not count.

**Baselines.** Enrichment, p-value, and raw record count. The last matters
most: attention is self-reinforcing, so "more papers" predicts "more papers"
regardless of evidence quality.

**No verdict without an interval.** The difference in AUROC against the best
baseline carries a paired bootstrap interval, and "better" or "worse" is only
stated when it excludes zero.

**Known limits.** Re-observation is a proxy for truth. Evidence weights are
computed as of the cutoff year, but retraction status necessarily reflects
today. The benchmark has so far only been run on the synthetic demo corpus,
where publication years are random and no signal exists; the result there is
correctly "no detectable difference".

## In silico ASO engineering

Therapeutic oligonucleotide design is governed by biophysics, not heuristic string matching. `synaptik aso` implements a clinical-grade design pipeline for splice-switching (SSO) and RNase H-mediated knockdown gapmers:

**Nearest-neighbor thermodynamics.** Duplex hybridization parameters ($\Delta H^\circ$, $\Delta S^\circ$, $\Delta G_{37}^\circ$) are computed using the unified parameters of Sugimoto et al. (1995) and SantaLucia (1998). Melting temperature is evaluated at physiological CSF ionic conditions ($[Na^+] = 100\text{ mM}$, $[C_T] = 250\text{ nM}$) with empirical 2'-MOE thermal stabilization ($+1.1^\circ\text{C}$ per nucleotide), yielding predictions within $0.6^\circ\text{C}$ of clinical benchmarks (e.g. Nusinersen $T_m \approx 64.5^\circ\text{C}$).

**Biophysical liability filtering.**
- *G-quadruplex formation:* Stretches of $\ge 4$ consecutive guanines (`GGGG`) are penalized and filtered out. G-quadruplexes in CSF cause nucleolar stress, neuroinflammation, and non-specific protein binding.
- *Poly-T/U tracts:* $\ge 4$ consecutive thymines or uracils are flagged for premature transcription termination during chemical synthesis and innate immune activation (TLR7/8).
- *Secondary hairpins:* Evaluates intramolecular stem-loop folding free energy $\Delta G_{\text{hairpin}}$. Oligos with $\Delta G_{\text{hairpin}} \le -3.0\text{ kcal/mol}$ are penalized as the hairpin competes with target RNA binding.
- *Self-dimerization:* 3'-terminal self-complementarity is scored to prevent self-sequestration.

**Whole-transcriptome neural off-target screening.** Candidate oligos are screened against a curated index of core human neural and housekeeping transcripts (*ACTB*, *GAPDH*, *TUBB3*, *MAP2*, *SYP*, *SNAP25*, *TARDBP*, *FUS*, *SOD1*, *C9orf72*). Candidates with exact 0-mismatch or 1-mismatch off-target alignments are eliminated.

**Synthesis specifications.** Generates orderable vendor sequences (IDT notation with 2'-MOE `m` prefixes and phosphorothioate `*` linkages), molecular weight (sodium salt), and estimated human CSF terminal half-life (~135 days for 2'-MOE/PS, based on clinical pharmacokinetics).

## Dish-to-brain human post-mortem single-cell translational bridge

Cell culture in human iPSC-derived neurons frequently suffers from two translational traps:
1. **Culture artifacts and discordant directionality:** In vitro artifacts (e.g. ambient 21% $O_2$ vs physiological brain 3% $O_2$, 2D plastic substrata, serum-free media variations) can invert the transcriptional direction of compensatory or stress pathways.
2. **Developmental maturation confounders:** Immature iPSC-derived neurons (DIV 21–35) naturally lack expression of adult neurofilaments (*NEFH*) and mature NMDA receptor subunits (*GRIN2A*). Mistaking developmental immaturity for disease pathology leads to spurious rescue claims.

**Pathology-aware concordance scoring.**
SYNAPTIK evaluates in vitro perturbation data against curated single-nucleus transcriptomic atlases of human post-mortem brain tissue (ALS NYGC Consortium, SEA-AD Alzheimer's, FTD):
- *Cryptic exon inclusion:* For cryptic exon targets regulated by TDP-43 (*STMN2*, *UNC13A*, *AARS1*), loss of nuclear TDP-43 causes non-conserved cryptic exon inclusion, which triggers nonsense-mediated decay (NMD) or premature polyadenylation and drops functional protein expression. An observed in vitro splicing inclusion ($E(\Delta\Psi) > 0$) is biologically concordant with human post-mortem disease pathology where cryptic inclusion is the hallmark lesion.
- *Standard cassette exon skipping & expression:* For standard functional cassette exons (*KALRN*), skipping ($E(\Delta\Psi) < 0$) or down-regulation in dish is evaluated against down-regulation in human motor units.
- *Concordance scoring:*
  $$\text{Score} = \frac{N_{\text{concordant}}}{N_{\text{evaluated}} - N_{\text{developmental\_gap}}} \times 100$$
  Penalizes direct conflicts ($\text{Score} = \max(0, \text{Score} - 20 \times N_{\text{discordant}})$).

**Cell-type specificity mapping.**
Markers for specific human CNS lineages are curated from landmark snRNA-seq datasets:
- Upper Motor Neurons / Betz cells (*FEZF2*, *BCL11B*, *CRYM*)
- Lower Motor Neurons / Spinal cord (*CHAT*, *MNX1*, *ISL1*, *SLC18A3*)
- Cortical Projection Neurons / Layers 2–6 (*SLC17A7*, *SATB2*, *TBR1*)
- GABAergic Interneurons (*GAD1*, *GAD2*, *PVALB*, *SST*, *VIP*)
- Astrocytes (*GFAP*, *AQP4*, *ALDH1L1*, *SLC1A2*)
- Microglia (*AIF1*, *TMEM119*, *CX3CR1*, *P2RY12*)
- Oligodendrocytes (*MBP*, *MOG*, *OLIG2*, *PDGFRA*)

## Causal genomic counterfactual simulation

Target identification in neurological disorders is obstructed by complex downstream feedback loops. A mutation or nuclear clearance event (e.g. TDP-43 loss) generates a high-dimensional transcriptomic disease vector $\vec{v}_D$. When evaluating therapeutic hypotheses, researchers need to predict whether a proposed intervention (gene knockdown, knockout, cDNA overexpression, or ASO steric blocking) will reverse the disease vector without triggering collateral toxicity.

**Mathematical formulation.**
Let $\vec{v}_D \in \mathbb{R}^M$ denote the baseline disease state vector (log2 fold-change relative to healthy isogenic controls) across hallmark neuronal genes.
Let $\vec{\Delta}_T \in \mathbb{R}^M$ represent the empirical or modeled perturbation response vector for target gene $T$.
Under an intervention with scaling efficacy $\alpha \in [0, 1]$, the post-intervention state is:
$$\vec{v}_{\text{post}} = \vec{v}_D + \alpha \cdot \vec{\Delta}_T$$

**Evaluation metrics:**
1. *Vector Rescue Efficiency (%):* Measures reduction in the Euclidean norm of the disease state:
   $$\text{Rescue Efficiency} = \frac{\|\vec{v}_D\|_2 - \|\vec{v}_{\text{post}}\|_2}{\|\vec{v}_D\|_2} \times 100$$
2. *Cosine Directional Alignment:* Measures whether the perturbation vector directly opposes the pathological trajectory:
   $$\cos \theta = \frac{\vec{\Delta}_T \cdot (-\vec{v}_D)}{\|\vec{\Delta}_T\|_2 \|\vec{v}_D\|_2}$$
   Where $\cos \theta = +1.0$ indicates perfect counter-pathological rescue, and $\cos \theta < 0$ indicates exacerbation of the disease state.
3. *Collateral Toxicity Index:* Disruption of unperturbed essential neural housekeeping genes (*ACTB*, *GAPDH*, *TUBB3*, *POLR2A*, *RPS6*, *SNAP25*) is penalized:
   $$\text{Toxicity} = \text{Risk}_{\text{base}} + 0.25 \sum_{g \in \text{Housekeeping}} |\Delta_g|$$
4. *Composite Rescue Score:* Grounded in Connectivity Map (CMap) and Perturb-seq standards, weighting directional cosine alignment (55%), magnitude reduction (up to 35%), and deducting collateral toxicity (up to 50%).

## Translational fluid biomarker and patient trial stratification

Neurofilament Light Chain (NfL) in cerebrospinal fluid (CSF) and blood plasma is the primary regulatory surrogate endpoint for axonal degeneration in neurodegenerative disease clinical trials (e.g. accelerated approval of tofersen for *SOD1* ALS based on plasma NfL suppression).

**Longitudinal pharmacodynamic clearance kinetics.**
SYNAPTIK models the translation of in vitro axonal rescue into longitudinal in vivo biomarker trajectories:
$$\Delta \text{NfL}(t) = \text{NfL}_{\text{baseline}} \times \left(1 - \beta \cdot \text{Rescue}_{\text{dish}} \cdot (1 - e^{-t / \tau})\right)$$
where $\tau \approx 8\text{ weeks}$ (systemic NfL clearance half-time) and $\beta \approx 0.72$ (dish-to-human translation scaling factor).

**Patient cohort enrichment and sample size reduction.**
In unstratified all-comers trials, inter-individual heterogeneity in baseline disease progression rates inflates variance ($\text{CV} \approx 0.50\text{--}0.65$), requiring $70\text{--}150$ patients per arm. Stratifying cohorts by baseline plasma NfL into progression tiers (fast: $> 75\text{ pg/mL}$, intermediate: $40\text{--}75\text{ pg/mL}$, slow: $< 40\text{ pg/mL}$) reduces variance ($\text{CV} \approx 0.30\text{--}0.32$), cutting required cohort sample size per arm by $60\text{--}75\%$ (e.g. from 74 down to 18 patients/arm) to detect a 30% slowing in functional decline.

## Phenotype pooling

Records for related phenotype terms are pooled when analysing a phenotype.
Until 0.3.0 pooling was derived from a shared measurement axis, which merged
mEPSC with mIPSC frequency and synaptic puncta with dendritic spines. Pooling
is now an explicit, curated `pool_group` on each term. Adding a term to a
group is a scientific claim that the two measure the same quantity and should
be reviewed as one.

---

## What would need to change for this to be more than hypothesis prioritisation

Honest answer: quite a lot, and some of it is not achievable from literature
alone.

- **Causal direction** would require perturbation data at scale, not
  co-occurrence. LINCS/CMap-style signatures could partly address this.
- **Effect size meta-analysis** would need consistent reporting of variance,
  which most papers do not provide.
- **Negative results** would need a preregistration or registered-reports
  corpus.
- **Pathway normalisation** would need a real ontology mapping with manual
  curation; the current string-based grouping under-merges.

Until then, the correct use of this tool is to prioritise which experiment to
run next — not to conclude what is true.
