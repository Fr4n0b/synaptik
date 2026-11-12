# Contributing to SYNAPTIK

Thank you for considering a contribution. SYNAPTIK is most useful when it is
shaped by the people who run iPSC-neuron experiments, so contributions from wet
labs are as valuable as code.

## Ways to contribute

### 1. Phenotype vocabulary (no coding experience needed)

The phenotype ontology lives in `synaptik/normalize/phenotype_ontology.py`. If
your lab describes a readout in a way the tool does not recognise, add the
phrasing to the closest term's `synonyms` tuple:

```python
PhenotypeTerm(
    "SYN:0003", "colocalized synaptic puncta density",
    ("colocalized puncta", "synapse density", "your lab's phrasing here"),
    ...
)
```

Run `synaptik renormalize` against your own knowledge base to confirm the change
recovers records, then open a pull request. Or open a **Phenotype vocabulary**
issue and someone will add it.

### 2. Validation cases

Known results the tool should recover: a published phenotype with a well-
established mechanism, a validated splice-altering variant, a confirmed drug
rescue. Each becomes a regression test. Open an issue with the PMIDs and the
expected outcome.

### 3. Correctness reports

If the tool reaches a conclusion you know to be wrong, please open a **Wrong
conclusion** issue. These are the most useful reports we receive: a confident
wrong answer is the failure mode this project exists to prevent.

### 4. Code

```bash
git clone https://github.com/Fr4n0b/synaptik.git
cd synaptik
python3 -m venv .venv && source .venv/bin/activate
pip install -e .
./run_tests.sh
```

All three stages of `run_tests.sh` must pass before a pull request is merged.

## Conventions

- **Tests encode scientific properties, not just code paths.** If you change
  ranking, weighting, tiering, or matching logic, add a test that asserts the
  property you care about (for example: "a single-lab association with higher
  enrichment must not outrank a replicated one"). See `tests/test_inference.py`.
- **Never silently coerce missing data.** A missing effect size is `None`, not
  `0`. A missing prediction is an error state, not "no effect". Absence of
  evidence and evidence of absence must stay distinguishable.
- **Caveats travel with results.** If a method has a known limitation, it belongs
  in the output object, not only in the documentation.
- **External APIs get a contract.** New harvesters or evidence sources must
  register the fields they depend on in `synaptik/ingest/contracts.py`, with a
  fixture.
- **No paywall circumvention.** Full text is fetched only from open-access
  sources or through an institution's licensed text-and-data-mining access.

## Larger changes

For new evidence classes, new data sources, or changes to the statistical core
(`infer/association.py`, `evidence/reconcile.py`), please open an issue first to
discuss the approach. These modules carry most of the scientific judgement in
the project and changes to them benefit from review by someone who will use the
output.

## Collaborating labs

If your lab is using SYNAPTIK and would like to coordinate (shared vocabulary,
shared validation sets, or a shared knowledge-base snapshot), open an issue
titled "Collaboration: <lab or institution>" and we will follow up.

## Code of conduct

Be constructive and assume good faith. Scientific disagreement is welcome;
personal criticism is not.
