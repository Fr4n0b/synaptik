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
