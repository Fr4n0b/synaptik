#!/usr/bin/env bash
# Full test suite. Run before any commit.
set -e
cd "$(dirname "$0")"
export PYTHONPATH="$PWD"

echo "=== unit + regression suite ==="
python3 -m tests.test_inference

echo
echo "=== lab-workflow suite (deposit, design, benchmark, pooling) ==="
python3 -m tests.test_lab_workflows

echo
echo "=== antisense oligonucleotide (ASO) engineering suite ==="
python3 -m tests.test_aso

echo
echo "=== human post-mortem single-cell translational bridge ==="
python3 -m tests.test_translation

echo
echo "=== causal genomic counterfactual simulation suite ==="
python3 -m tests.test_simulation

echo
echo "=== translational fluid biomarker & patient stratification suite ==="
python3 -m tests.test_biomarker

echo
echo "=== API contract self-test (offline) ==="
python3 -m synaptik.ingest.contracts

echo
echo "=== CLI smoke test ==="
TMPDB=$(mktemp -u /tmp/synaptik-smoke-XXXX.db)
python3 -m synaptik.cli.main demo --db "$TMPDB" > /dev/null
for cmd in \
  "phenotype 'SYN1 puncta density' --limit 3" \
  "splicing" \
  "lookup --target TARDBP --assay splicing_analysis --div 21" \
  "compare 'SYN1 puncta density' 'mEPSC frequency'" \
  "report" \
  "design --target TARDBP --assay imaging_synaptic_puncta --div 21" \
  "design --target TARDBP --assay imaging_synaptic_puncta --cv 0.2" \
  "deposit splicing examples/SE.MATS.JC.txt --perturbed sample1" \
  "deposit phenotypes examples/phenotypes.csv" \
  "benchmark" \
  "benchmark --spliceai" \
  "aso --sequence TTCTGTAACAGTTTTTGTTTTTTCAGGTGAGTATTCAGGT --gene STMN2" \
  "translate --splicing examples/voila_majiq.tsv --disease als" \
  "translate --genes STMN2 UNC13A TARDBP --disease als" \
  "simulate --disease als --target ATXN2" \
  "simulate --disease als --screen" \
  "biomarker --disease als --target STMN2 --rescue-pct 75" \
  "contracts --self-test"
do
  eval "python3 -m synaptik.cli.main $cmd --db '$TMPDB'" > /dev/null 2>&1 \
    && echo "  PASS  synaptik $cmd" \
    || { echo "  FAIL  synaptik $cmd"; exit 1; }
done

echo
echo "=== JSON output validity ==="
for cmd in \
  "phenotype 'SYN1 puncta density'" \
  "splicing" \
  "report" \
  "design --target TARDBP --assay imaging_synaptic_puncta" \
  "deposit splicing examples/SE.MATS.JC.txt --perturbed sample1" \
  "deposit phenotypes examples/phenotypes.csv" \
  "benchmark" \
  "benchmark --spliceai" \
  "aso --sequence TTCTGTAACAGTTTTTGTTTTTTCAGGTGAGTATTCAGGT --gene STMN2" \
  "translate --splicing examples/voila_majiq.tsv --disease als" \
  "translate --genes STMN2 UNC13A TARDBP --disease als" \
  "simulate --disease als --target ATXN2" \
  "simulate --disease als --screen" \
  "biomarker --disease als --target STMN2 --rescue-pct 75"
do
  eval "python3 -m synaptik.cli.main $cmd --db '$TMPDB' --json" \
    | python3 -c "import json,sys; json.load(sys.stdin)" \
    && echo "  PASS  $cmd --json" \
    || { echo "  FAIL  $cmd --json"; exit 1; }
done

echo
echo "=== refusing unsafe input ==="
# An rMATS table without a stated perturbed group must be refused.
if python3 -m synaptik.cli.main deposit splicing examples/SE.MATS.JC.txt \
     --db "$TMPDB" > /dev/null 2>&1; then
  echo "  FAIL  deposit accepted an rMATS table with no --perturbed"; exit 1
else
  echo "  PASS  deposit refuses an rMATS table with no --perturbed"
fi

rm -f "$TMPDB"*
echo
echo "All tests passed."
