# SYNAPTIK: Complete Runbook

From nothing installed to answering your own research questions.

Every command here has been run. Expected output is shown so you can tell
whether a step worked. Where something commonly goes wrong, the fix is inline
rather than in a troubleshooting appendix.

**Time budget**
| Part | What you get | Time |
|---|---|---|
| 1 | Installed, exploring demo data | 10 min |
| 2 | Accounts and keys | 20 min |
| 3 | Real corpus for your topic | 2–4 h (mostly unattended) |
| 4 | Structured records extracted | 1–3 h (unattended) |
| 5 | Answering research questions | ongoing |
| 6 | Sequence + population evidence | 1–2 h setup |

Parts 1–2 are required. Part 3 onward is where it becomes useful to you.

---

# PART 1 — Install and verify (10 minutes)

## 1.1 Check Python

```bash
python3 --version
```

Need **3.10 or higher**. If you have 3.9 or older:

```bash
# macOS
brew install python@3.12
# Ubuntu/Debian
sudo apt update && sudo apt install python3.12 python3.12-venv
```

## 1.2 Unpack and install

```bash
tar xzf synaptik-2.0.0.tar.gz
cd synaptik

python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate

pip install --upgrade pip
pip install -e .
```

**Verify:**
```bash
synaptik --help
```

You should see 14 commands listed (`init`, `demo`, `preflight`, `contracts`,
`harvest`, `extract`, `lookup`, `phenotype`, `splicing`, `compare`, `gene`,
`variant`, `ui`, `report`).

> **`synaptik: command not found`** — your venv isn't active. Run
> `source .venv/bin/activate` again. You need this in every new terminal.

## 1.3 Run the test suite

```bash
./run_tests.sh
```

Expect `All tests passed.` — 114 checks.

This matters: the suite asserts the properties that make output trustworthy
(a single-lab artifact with 12× higher enrichment must not outrank a
well-replicated finding; retracted evidence must weigh zero; an overexpression
study must never be called a close match for a patient-mutation query). If any
check fails, **stop and tell me** — the tool would be producing misleading
output.

## 1.4 Explore the demo

```bash
synaptik demo --db data/demo.db
```

```
Building demo knowledge base at data/demo.db ...
  sources: 310
  records: 310
  splicing_events: 11
  distinct_labs: 114
```

The data is **synthetic** — ignore the biology entirely. It is deliberately
structured like the real literature: one over-productive lab, a retracted
paper, conflicting reports, a publication-bias signature. What you are
evaluating is whether the *reasoning* is something you would trust.

```bash
synaptik phenotype "SYN1 puncta density" --db data/demo.db -v
```

```
pathway                        tier          enrich         q  labs  assays
neurotransmitter release       strong          16.9  1.58e-06     9       2
autophagy                      weak           111.4  2.33e-05     5       1
synaptic vesicle cycle         strong           7.1  7.26e-05    10       2
exotic kinase                  single-lab     225.2  3.46e-09     1       1
```

**Read row 4 carefully.** `exotic kinase` has **13× the enrichment** of the
top hit and a smaller q-value, and it ranks fourth, tiered `single-lab`. One
group, one assay, no null results across nine papers. A tool that sorted by
enrichment or p-value — which is what most literature mining does — would put
it first and you would chase it.

That inversion is the entire point of the tool. If it doesn't strike you as
right, stop here and tell me, because everything downstream depends on it.

Two more:

```bash
synaptik splicing --db data/demo.db
synaptik gene TARDBP --db data/demo.db --offline
synaptik ui --db data/demo.db        # then open http://127.0.0.1:8420
```

Note what `gene` says: 301 papers is still only **"CONSISTENT"**, not
convergent, because it is all one epistemic class. That is deliberate.

---

# PART 2 — Accounts and keys (20 minutes)

Three of these are free and take two minutes. The fourth costs money.

## 2.1 NCBI API key — free, required in practice

Triples your PubMed rate limit (3 → 10 requests/sec). Without it, harvesting
takes 3× longer.

1. Go to https://www.ncbi.nlm.nih.gov/account/
2. Register or sign in (use your MGH address)
3. Account settings → **API Key Management** → Create an API Key
4. Copy the key

## 2.2 SCImago journal table — free, recommended

Gives real journal quartiles instead of the small built-in table.

1. Go to https://www.scimagojr.com/journalrank.php
2. Click **Download data** (top right) — semicolon-delimited CSV
3. Save it:

```bash
mkdir -p ~/synaptik-data
mv ~/Downloads/scimagojr*.csv ~/synaptik-data/scimago.csv
```

## 2.3 Anthropic API key — costs money, required for extraction

Extraction turns papers into structured records. Without it you have a
searchable paper index but no analysis.

1. Go to https://console.anthropic.com/
2. Sign up → **API Keys** → Create Key
3. Add credit — **$20 is plenty to start**

Cost is roughly **$0.01–0.03 per full-text paper**. A 1,000-paper corpus is
$10–30, one time. You will see an exact estimate before spending anything.

> If this needs a grant code, do Part 3 first — harvesting is free, and you
> can decide about extraction once you see the corpus.

## 2.4 Save your configuration

```bash
cat > ~/.synaptik.env <<'EOF'
export SYNAPTIK_DB=$HOME/synaptik-data/synaptik.db
export SYNAPTIK_EMAIL=you@mgh.harvard.edu
export NCBI_API_KEY=paste_ncbi_key_here
export ANTHROPIC_API_KEY=sk-ant-paste_here
export SCIMAGO_CSV=$HOME/synaptik-data/scimago.csv
EOF

chmod 600 ~/.synaptik.env
```

`chmod 600` matters — this file holds API keys.

**Load it in every new terminal:**
```bash
source ~/.synaptik.env
```

Or add `source ~/.synaptik.env` to your `~/.bashrc` / `~/.zshrc`.

**Replace `you@mgh.harvard.edu` with your real address.** NCBI and Europe PMC
require a contact address; a fake one can get your IP blocked.

---

# PART 3 — Build a real corpus (2–4 hours, mostly unattended)

## 3.1 Preflight — never skip this

```bash
cd ~/path/to/synaptik
source .venv/bin/activate
source ~/.synaptik.env

synaptik preflight --email $SYNAPTIK_EMAIL --ncbi-key $NCBI_API_KEY
```

**Expected:**
```
preflight: checking API response contracts...
  Europe PMC: 3 records returned
  PubMed: 3 PMIDs returned
[Europe PMC /search] OK (3 samples checked)
All contracts satisfied. Safe to harvest.
PASS - safe to harvest.
```

This costs three requests and prevents the worst failure mode in the system:
an upstream API renames a field, the parser silently returns nothing, and a
30,000-paper harvest completes "successfully" with every journal name missing.
You would find out weeks later.

> **"No records returned from any API"** — network blocked. If you are on an
> institutional network, ask IT to allow `www.ebi.ac.uk` and
> `eutils.ncbi.nlm.nih.gov`. Full list in `docs/DEPLOYMENT.md`.
>
> **"CONTRACT VIOLATION"** — an API changed. The message names the exact
> field. Send it to me; it is a one-line fix.

## 3.2 Small probe first — do not skip to the full harvest

```bash
synaptik harvest --email $SYNAPTIK_EMAIL --ncbi-key $NCBI_API_KEY \
    --scimago $SCIMAGO_CSV \
    --only ipsc_neuron_core --max-per-query 20
```

Takes about a minute. Then **look at what actually landed**:

```bash
synaptik report
```

Check that **titles, journals, and years are populated**. A parser mismatch
preflight didn't catch shows up here as systematically empty fields — far
cheaper to find on 20 papers than on 30,000.

## 3.3 Choose your scope

Do **not** harvest everything. Pick the question you are stuck on now.

See all available query sets:
```bash
synaptik harvest --help | grep -A2 "only"
```

The 22 sets include: `ipsc_neuron_core`, `ngn2_induction`, `brain_organoid`,
`isogenic_ipsc`, `synaptic_phenotype`, `electrophysiology`,
`neurite_morphology`, `splicing_neuro`, `rbp_regulation`,
`splice_therapeutics`, `splicing_methods`, `als_ftd`, `autism_ndd`,
`schizophrenia_psychiatric`, `alzheimers_tau`, `parkinsons`,
`rare_neurogenetic`, `mitochondria_metabolism`, `proteostasis`,
`transcriptomics_neuron`, `drug_rescue`, `crispr_screen_neuron`.

**Recommended starting scope** for synaptic phenotypes + splicing in
NGN2 neurons:

```bash
nohup synaptik harvest \
    --email $SYNAPTIK_EMAIL \
    --ncbi-key $NCBI_API_KEY \
    --scimago $SCIMAGO_CSV \
    --only ipsc_neuron_core ngn2_induction synaptic_phenotype \
           splicing_neuro rbp_regulation isogenic_ipsc \
    --max-per-query 800 \
    --fulltext \
    > ~/harvest.log 2>&1 &
```

`nohup ... &` runs it in the background so closing your laptop won't kill it.

**Monitor:**
```bash
tail -f ~/harvest.log     # Ctrl-C to stop watching (harvest continues)
```

Expect **2–4 hours** for ~4,000 papers. It is **resumable** — if it dies,
re-run the identical command and it skips everything already stored.

## 3.4 Check the corpus

```bash
synaptik report
```

```
  sources          3,847
  documents        1,102        <- papers with full text
  records          0            <- nothing extracted yet, expected
  distinct labs    1,455
```

**The number that matters is `documents` vs `sources`.** That ratio is your
full-text coverage — typically 25–35%. Only full-text papers give good
extraction, because Methods sections carry DIV, differentiation protocol, and
antibodies. Abstract-only papers still contribute, but with less detail.

---

# PART 4 — Extract structured records (1–3 hours, unattended)

## 4.1 Estimate the cost

```bash
synaptik extract --dry-run
```

```
1102 papers queued for extraction.
Estimated cost at ~$0.02/paper: $22.04
```

No money is spent by `--dry-run`.

## 4.2 Test on 50 papers first

```bash
synaptik extract --limit 50
```

Costs about $1. **Read the summary carefully:**

```
sources: 0 new / 0 seen | records: 312 kept, 24 rejected as unsupported
(7.1% hallucination rate) | 284s
```

**The hallucination rate is the number to watch.** It is the fraction of
proposed records rejected because their supporting quote did not appear
verbatim in the source.

| Rate | Meaning |
|---|---|
| **< 10%** | Normal. Proceed. |
| **10–25%** | Acceptable but watch it. |
| **> 25%** | Something is wrong — usually extraction from abstracts rather than full text. Check `documents` in `synaptik report`. |

It also prints **unmapped phenotypes** — phrasings the vocabulary didn't
recognise. Note these; you'll fix them in 4.4.

## 4.3 Run the rest

```bash
nohup synaptik extract > ~/extract.log 2>&1 &
tail -f ~/extract.log
```

Resumable, like harvest.

## 4.4 Extend the vocabulary — this is the highest-leverage 10 minutes

```bash
synaptik report
```

Look at **`unmapped pheno`**. Above 25% means analyses are blind to a quarter
of your records.

The extraction log lists the most frequent unmapped strings. Open
`synaptik/normalize/phenotype_ontology.py`, find the term closest to what you
saw, and add the phrasing to its `synonyms` tuple:

```python
PhenotypeTerm(
    "SYN:0003", "colocalized synaptic puncta density",
    (
        "colocalized puncta", "synapse density", "synapse number",
        "puncta per 10 um dendrite",        # <-- your addition
        "syn1/psd95 overlap",               # <-- your addition
    ),
    measurement_axis="puncta_count_per_length",
    ...
)
```

Then re-normalize. This does **not** re-run the LLM, so it costs nothing:

```bash
synaptik renormalize
```

```
310 records; 61 currently unmapped (20%)
This does not call the LLM and costs nothing.

  61 record(s) changed mapping
  unmapped: 61 -> 0 (0%)
  recovered 61 record(s) for analysis
```

If nothing changes, it prints the still-unmapped phrasings so you can see
exactly what to add next. Repeat until unmapped is below ~15%.

Re-run your analyses afterwards — records that were previously invisible now
contribute.

## 4.5 Read the lab clustering — do not skip this

```bash
synaptik report
```

Scroll to **"Most prolific groups"**.

Independent-lab count is the denominator for **every** confidence claim the
system makes. The clustering that produces it is heuristic — senior author
plus institution.

- A "group" with an implausible number of papers → it **over-merged**
- A PI you know appearing twice → it **under-merged** (they moved institutions)

Everything else in the system degrades gracefully when wrong. This number does
not. To correct it, create `~/synaptik-data/lab_curation.json`:

```json
{
  "PMID:12345678": "lab_smith_mgh",
  "PMID:23456789": "lab_smith_mgh"
}
```

---

# PART 5 — Answer your research questions

Everything below is free and instant — no more API calls.

## 5.1 "Has this experiment been done?"

Describe the experiment you are planning:

```bash
synaptik lookup \
    --model ipsc_neuron_ngn2 \
    --div 21 \
    --perturbation mutation_patient \
    --target SCN2A \
    --assay imaging_synaptic_puncta \
    --markers SYN1 PSD95 \
    --phenotype "colocalized puncta density"
```

**Reading the output:**

| Tier | Meaning |
|---|---|
| **EXACT** | Same system, maturity, perturbation class, readout. **Read these before running anything.** |
| **CLOSE** | Differs in one dimension — check the `differs:` line, since *which* dimension matters enormously |
| **RELATED** | Useful context, not a direct answer |
| **ADJACENT** | Background |

Each match shows the **actual protocol**: differentiation method, DIV,
coculture, markers, replicate count, statistical test, whether isogenic.

Two verdict lines to take seriously:

- **"Direction of effect is disputed"** — do not average. Read the conflicting
  reports and work out what differs between them. That difference is often the
  real finding.
- **"All matches come from a single group"** — done once, by one lab. Your
  replication has genuine value.

> Not all flags are required. `--target SCN2A --assay imaging_synaptic_puncta`
> alone works fine. More flags = more precise tiering.

## 5.2 "What does my phenotype mean?"

```bash
synaptik phenotype "reduced colocalized SYN1/PSD95 puncta density" \
    --direction decrease -v
```

**Read the tier, not the ranking:**

| Tier | What to do |
|---|---|
| **strong** | ≥3 labs, ≥2 assay families, survives leave-one-lab-out. Reasonable to build on. |
| **moderate** | Real but thinner. Read the primary sources. |
| **weak** | Hypothesis. Usually directions conflict or one method dominates. |
| **single-lab** | One group. **High value to test yourself precisely because nobody else has.** |

Four numbers that should change your behaviour:

- **labs = 1** → everything else is noise until this is ≥2
- **assays = 1** → may be a measurement artifact, not biology
- **leave-one-lab-out FAILS** → one group's narrative, not a consensus
- **direction consistency < 70%** → the field disagrees about the *sign*; the
  discrepancy usually points at an uncontrolled variable

## 5.3 Cross-check two readouts

If you have both an imaging and an ephys phenotype:

```bash
synaptik compare "SYN1 puncta density" "mEPSC frequency"
```

Convergent pathways across orthogonal readouts are genuine support. Divergence
means one may be an artifact.

## 5.4 Splicing analysis

```bash
synaptik splicing --gene STMN2
synaptik splicing --regulator TDP-43
synaptik splicing                    # everything in the corpus
```

Three things to act on:

- **Regulator inference** — cryptic exon inclusion across multiple genes is the
  canonical TDP-43 nuclear-loss signature, inferred even when no paper names it
- **`NMD risk: high`** — changes the therapeutic question entirely. If the
  aberrant transcript is degraded, your defect is *loss of protein*. **Measure
  protein, not just transcript.**
- **`⚠` on ΔPSI** — sign is uncertain; verify against the source figure before
  believing the direction

## 5.5 Drug candidates

```bash
synaptik phenotype "reduced synaptic puncta density" \
    --direction decrease --drugs --email $SYNAPTIK_EMAIL
```

Only **direction-matched** compounds appear: a compound is listed because its
mechanism *opposes* the observed change. An inhibitor for a down-regulated
pathway is excluded with an explicit reason, not ranked lower.

- **`bbb_plausible: false` does not disqualify a compound for dish work** — it
  only matters for in vivo follow-up
- **Efflux warnings override the CNS score** — loperamide is the textbook case
- **"already tested" annotations are the highest-value line in the output** —
  somebody put this compound on neurons and reported what happened

## 5.6 Web interface

```bash
synaptik ui
```

Open http://127.0.0.1:8420. Five tabs covering everything above. There is
**no authentication** — to share with the lab, use an SSH tunnel:

```bash
# from a labmate's machine
ssh -L 8420:localhost:8420 you@yourserver
```

---

# PART 6 — Sequence and population evidence

This is where the tool does something a literature search cannot.

## 6.1 Gene dossier — works immediately, no setup

```bash
synaptik gene SCN2A --email $SYNAPTIK_EMAIL
```

Pulls gnomAD constraint, ClinVar, expression, and your literature corpus, then
reports where they **disagree**.

**The section to read first is "Cross-class contradictions."** Examples it
surfaces:

- Literature reports haploinsufficiency phenotypes, but gnomAD shows the gene
  is LoF-tolerant with homozygous LoF carriers in healthy people
- Literature reports neuronal phenotypes, but the gene isn't expressed in
  neurons
- The gene is commonly essential in CRISPR screens, so published phenotypes may
  be secondary to impaired viability

These are statements a literature-only tool structurally cannot make.

Run offline (literature only, no external APIs):
```bash
synaptik gene SCN2A --offline
```

## 6.2 Install SpliceAI

Required for variant analysis. Two options.

### Option A — precomputed scores (recommended)

Fast, no GPU, no TensorFlow. **SNVs only.**

```bash
pip install pysam
```

Download the precomputed VCF (free BaseSpace account required):
https://basespace.illumina.com/s/otSPW8hnhaZR

Get `spliceai_scores.raw.snv.hg38.vcf.gz` **and** its `.tbi` index (~30 GB).

```bash
mv ~/Downloads/spliceai_scores.raw.snv.hg38.vcf.gz* ~/synaptik-data/
echo 'export SPLICEAI_VCF=$HOME/synaptik-data/spliceai_scores.raw.snv.hg38.vcf.gz' \
    >> ~/.synaptik.env
source ~/.synaptik.env
```

### Option B — local inference

Handles indels. Needs TensorFlow and a reference genome.

```bash
pip install spliceai tensorflow
# GRCh38 FASTA from https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/
echo 'export SPLICEAI_REF=$HOME/synaptik-data/GRCh38.fa' >> ~/.synaptik.env
```

## 6.3 Validate SpliceAI before trusting it

**Do not skip this.** I could not test the adapter against real model weights.

Run 3–5 published pathogenic splice variants where you know the expected
answer:

```bash
synaptik variant "chr2:165991235:G>A" --gene SCN2A \
    --spliceai-vcf $SPLICEAI_VCF
```

Check the delta scores match published values. If they are all zero or all
identical, the adapter is misreading the file — **tell me the output** and it
is a quick fix.

## 6.4 Analyze your patient variants

```bash
synaptik variant \
    "chr8:79611214:G>A" \
    "chr2:165991235:C>T" \
    --gene STMN2 \
    --spliceai-vcf $SPLICEAI_VCF
```

**Four states, each meaning something different:**

| State | Meaning | What to do |
|---|---|---|
| **`predicted_only`** | Model is confident; nobody has published it | **Highest experimental value.** A minigene assay tests a specific, sequence-grounded hypothesis with no prior art. |
| **`corroborated`** | Prediction and literature agree | Two independent epistemic classes. Reasonable to build on. |
| **`contradicted`** | They disagree | Check genome build and transcript isoform first — coordinate mismatches explain many apparent contradictions. |
| **`observed_only`** | Literature reports it; sequence doesn't predict it | Mechanism is probably indirect — acting through a regulator, not a splice site. |

Results are ranked by **experimental priority**, with untested predictions
first. That ordering is deliberate: the tool is most useful pointing at what
nobody has done.

> **Caveat that ships with every `predicted_only` result:** absence from your
> corpus is not absence from the literature. The paper may be paywalled or not
> harvested. **Check PubMed manually before claiming novelty.**

---

# Daily workflow, once set up

```bash
cd ~/path/to/synaptik
source .venv/bin/activate
source ~/.synaptik.env

# before designing an experiment
synaptik lookup --target YOURGENE --assay imaging_synaptic_puncta --div 21

# after seeing a phenotype
synaptik phenotype "what you observed" --direction decrease -v

# for a candidate gene
synaptik gene YOURGENE

# for patient variants
synaptik variant "chr:pos:REF>ALT" --gene YOURGENE --spliceai-vcf $SPLICEAI_VCF

# or just
synaptik ui
```

**Monthly:** refresh the corpus and check for API drift.

```bash
synaptik harvest --email $SYNAPTIK_EMAIL --ncbi-key $NCBI_API_KEY \
    --scimago $SCIMAGO_CSV --fulltext
synaptik extract
```

**Backup** (single SQLite file):
```bash
sqlite3 $SYNAPTIK_DB ".backup ~/synaptik-data/backup-$(date +%F).db"
```

---

# Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `synaptik: command not found` | venv not active | `source .venv/bin/activate` |
| `ANTHROPIC_API_KEY not set` | env not loaded | `source ~/.synaptik.env` |
| Preflight: "No records returned" | network blocked | Ask IT to allow `www.ebi.ac.uk`, `eutils.ncbi.nlm.nih.gov` |
| Preflight: "CONTRACT VIOLATION" | API changed a field | Message names the field — send it to me |
| Harvest very slow | rate limits | Add `--ncbi-key`; do not raise limits in code |
| `database is locked` | two writers | Only one harvest/extract per database at a time |
| Hallucination rate > 25% | extracting from abstracts | Check `documents` vs `sources` in `report` |
| "No records match" in phenotype | vocabulary gap | Add synonyms (step 4.4) |
| All SpliceAI deltas zero | adapter misreading VCF | Send me the output |
| Lab counts look wrong | clustering heuristic | Curation file (step 4.5) |

---

# Three things that will bite you if you forget them

**1. Run `preflight` before every real harvest.** Every API integration was
validated against fixtures, not live responses — the build sandbox blocked all
network egress. Preflight is what confirms the parsers still match reality.

**2. Read the lab clustering before trusting any replication count.**
Independent-lab count is the denominator for every confidence claim. It is
heuristic. Everything else degrades gracefully when wrong; that number doesn't.

**3. Validate SpliceAI on known variants first.** The adapter is written
against the documented output format and tested against fixtures, not real
weights.

---

# What this tool does and does not claim

It computes **evidence structure**, not causal inference.

When it says a pathway is associated with your phenotype, the claim is: *among
published experiments that measured P, pathway W is implicated more often than
background, supported by N independent groups across K orthogonal assay
families, with these robustness properties and these failure modes.*

It is **not** claiming W causes P. Literature co-occurrence cannot establish
that.

**It cannot tell you:**
- Whether pathway W *causes* phenotype P
- The true effect size (most papers don't report variance)
- Whether anyone tried and failed (unpublished nulls are unrecoverable)
- Whether a finding holds in *your* line and background

**Coverage is real but bounded:** near-complete at abstract level (~37M PubMed
records reachable), roughly 20–25% at full-text level. The remainder is
paywalled and mining it would violate publisher agreements MGH has signed.

If a result matters to your next six months of work, **read the primary
sources.** This tool tells you which ones, and how sceptical to be.
