# Deployment

Setting SYNAPTIK up on a lab server, and running it without babysitting it.

## Why a lab server

A laptop cannot comfortably hold the corpus or run an overnight harvest without
sleeping. Full cloud deployment adds procurement friction and a security review
you do not need for a read-only literature tool. A lab server sits in the
middle: persistent storage, cron jobs that survive you closing your laptop, and
a database file you can `scp` to a collaborator.

**Sizing.** For ~30,000 papers with full text and extracted records:

| Resource | Requirement |
|---|---|
| Disk | 20–40 GB (full text dominates) |
| RAM | 4 GB (8 GB if you run extraction in parallel) |
| CPU | 2 cores; the pipeline is I/O-bound, not compute-bound |
| Network | Outbound HTTPS to the APIs listed below |

## Install

```bash
sudo useradd -m -s /bin/bash synaptik
sudo -u synaptik -i

git clone https://github.com/Fr4n0b/synaptik.git ~/synaptik && cd ~/synaptik
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
pip install anthropic      # only if using LLM extraction

mkdir -p ~/data ~/logs
```

### Required network access

If your institution restricts egress, request these hosts:

| Host | Purpose | Required |
|---|---|---|
| `www.ebi.ac.uk` | Europe PMC — metadata and full text | **yes** |
| `eutils.ncbi.nlm.nih.gov` | PubMed — MeSH terms, PubMed-only records | recommended |
| `api.biorxiv.org` | bioRxiv/medRxiv preprints | optional |
| `www.ebi.ac.uk/chembl` | ChEMBL — compounds and mechanisms | drug layer only |
| `api.platform.opentargets.org` | Open Targets — target tractability | drug layer only |
| `api.anthropic.com` | LLM extraction | extraction only |

Europe PMC is the only hard requirement — it is the sole source of full text,
and without full text you get abstract-level extraction, which loses the
Methods section where DIV, protocol, and antibodies live.

### Get an NCBI API key

Free, immediate, and triples your PubMed rate limit (3→10 req/s):
<https://www.ncbi.nlm.nih.gov/account/> → Settings → API Key Management.

### Get the journal quality table

Download the SCImago Journal Rank CSV from
<https://www.scimagojr.com/journalrank.php> ("Download data" — semicolon
delimited) and save it as `~/data/scimago.csv`. Without it, journal tiering
falls back to a small built-in table covering the major neuroscience and
stem-cell journals.

## Configuration

```bash
cat > ~/.synaptik.env <<'EOF'
export SYNAPTIK_DB=/home/synaptik/data/synaptik.db
export SYNAPTIK_EMAIL=you@mgh.harvard.edu
export NCBI_API_KEY=your_key_here
export ANTHROPIC_API_KEY=sk-ant-...
export SCIMAGO_CSV=/home/synaptik/data/scimago.csv
EOF
chmod 600 ~/.synaptik.env
echo 'source ~/.synaptik.env' >> ~/.bashrc
```

`chmod 600` matters — the file holds API keys.

## First harvest

Do not start a full harvest blind. Work up in three steps.

```bash
cd ~/synaptik && source .venv/bin/activate && source ~/.synaptik.env

# 1. Confirm the APIs still return what the parsers expect
synaptik preflight --email $SYNAPTIK_EMAIL --ncbi-key $NCBI_API_KEY

# 2. Small probe — inspect what actually lands before scaling up
synaptik harvest --email $SYNAPTIK_EMAIL --only ipsc_neuron_core \
    --max-per-query 20 --scimago $SCIMAGO_CSV
synaptik report

# 3. Full harvest once the probe looks right
nohup synaptik harvest --email $SYNAPTIK_EMAIL --ncbi-key $NCBI_API_KEY \
    --scimago $SCIMAGO_CSV --fulltext > ~/logs/harvest.log 2>&1 &
```

**Inspect the probe output before step 3.** Check that titles, journals, years,
and author lists are populated. A parser mismatch that preflight did not catch
shows up here as systematically empty fields — much cheaper to find on 20
papers than on 30,000.

Expect the full harvest to take 4–12 hours depending on rate limits and how
much full text is available. It is resumable: if it dies, re-run the same
command and it skips everything already stored.

## Extraction

Extraction is the expensive stage. Estimate before committing:

```bash
synaptik extract --dry-run
```

Roughly $0.01–0.03 per full-text paper with a small model, so a 30,000-paper
corpus is a few hundred dollars once, with incremental updates nearly free.

```bash
# Test on a small batch and check the hallucination rate first
synaptik extract --limit 50
# Then the rest
nohup synaptik extract > ~/logs/extract.log 2>&1 &
```

The run reports a **hallucination rate** — the fraction of proposed records
rejected because their supporting quote did not appear in the source. Under 10%
is normal. Above 25% suggests the model is struggling with your corpus; try a
larger model or check whether full text is being retrieved correctly.

It also reports the most frequent **unmapped phenotypes**. These are the
highest-value additions to `normalize/phenotype_ontology.py` — each one you add
recovers records that are currently invisible to analysis.

## Scheduled updates

```cron
# Weekly incremental harvest, Sunday 02:00
0 2 * * 0 cd /home/synaptik/synaptik && . .venv/bin/activate && \
  . ~/.synaptik.env && synaptik harvest --email $SYNAPTIK_EMAIL \
  --ncbi-key $NCBI_API_KEY --scimago $SCIMAGO_CSV --fulltext \
  >> ~/logs/cron-harvest.log 2>&1

# Extraction, Sunday 08:00 (after the harvest has settled)
0 8 * * 0 cd /home/synaptik/synaptik && . .venv/bin/activate && \
  . ~/.synaptik.env && synaptik extract >> ~/logs/cron-extract.log 2>&1

# Monthly contract check — catches API drift before it corrupts a harvest
0 6 1 * * cd /home/synaptik/synaptik && . .venv/bin/activate && \
  . ~/.synaptik.env && synaptik preflight --email $SYNAPTIK_EMAIL \
  >> ~/logs/preflight.log 2>&1 || \
  mail -s "SYNAPTIK: API contract drift detected" you@mgh.harvard.edu \
  < ~/logs/preflight.log
```

That last job is the important one. It emails you when an upstream API changes,
rather than letting a scheduled harvest quietly write records with missing
fields for months.

## Web interface

```bash
synaptik ui --db $SYNAPTIK_DB
```

Binds to `127.0.0.1:8420` by default. **There is no authentication.** To share
it with the lab, use an SSH tunnel rather than binding to `0.0.0.0`:

```bash
# from a lab member's machine
ssh -L 8420:localhost:8420 synaptik@labserver
# then open http://localhost:8420
```

If you genuinely need it on the network, put it behind an authenticating
reverse proxy and a production WSGI server:

```bash
pip install gunicorn
gunicorn -w 4 -b 127.0.0.1:8420 \
  "synaptik.ui.server:create_app('/home/synaptik/data/synaptik.db')"
```

Flask's built-in server is for development only.

## Backups

The knowledge base is a single SQLite file. Back it up with SQLite's own
command so you get a consistent snapshot even while the server is reading:

```bash
sqlite3 $SYNAPTIK_DB ".backup /backup/synaptik-$(date +%F).db"
```

Weekly is plenty. A rebuild from scratch costs a harvest plus an extraction run,
so the backup mainly saves you the extraction spend.

## Health checks

```bash
synaptik report
```

Watch three numbers:

- **quantified fraction** — below 40% means most records have no effect size,
  so effect-size statistics will be thin.
- **null-result fraction** — below 5% means publication bias is severe in your
  corpus, or extraction is missing negative findings. Check a few papers by
  hand.
- **unmapped phenotype fraction** — above 25% means the vocabulary needs
  extending; those records are invisible to phenotype analysis.

The report also lists the most prolific research groups. **Read this list.**
Lab clustering is heuristic, and it is the denominator for every replication
claim the system makes. If one "group" has an implausible number of papers, it
has probably over-merged; if a PI you know appears twice, it has under-merged.

### Correcting lab clustering

Create `~/data/lab_curation.json` mapping source IDs to lab IDs:

```json
{
  "PMID:12345678": "lab_smith_mgh",
  "PMID:23456789": "lab_smith_mgh"
}
```

Then pass it to `LabResolver(curation_file=...)`. Curated assignments override
the automatic clustering.

## Troubleshooting

**`host_not_allowed` / 403 from every request.** Institutional egress
filtering. The hosts to request are listed above.

**Preflight fails with "No records returned from any API."** Network is blocked
or the endpoints are down. This is treated as a failure, not an empty result,
precisely so it cannot pass silently.

**Preflight reports a contract violation.** An upstream API changed a field
name. The message names the field and shows the keys actually present. Fix the
parser in `synaptik/ingest/`, add the new field name to the contract in
`contracts.py`, and re-run.

**Harvest is very slow.** Expected — rate limits are deliberately conservative.
An NCBI API key triples the PubMed limit. Do not raise the limits in
`ingest/http.py` unless you enjoy being blocked.

**`database is locked`.** Another process is writing. SQLite is in WAL mode so
concurrent reads are fine, but only one writer at a time. Do not run two
harvests against one database.

**Extraction reports a high hallucination rate.** Check that full text is
actually being retrieved (`synaptik report`, `documents` count). Extraction from
abstracts alone is much harder and produces more rejected records.
