# Changelog

All notable changes to SYNAPTIK are recorded here.

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
