#!/usr/bin/env bash
# Full test suite. Run before any commit.
set -e
cd "$(dirname "$0")"
export PYTHONPATH="$PWD"

echo "=== unit + regression suite ==="
python3 -m tests.test_inference

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
  "contracts --self-test"
do
  eval "python3 -m synaptik.cli.main $cmd --db '$TMPDB'" > /dev/null 2>&1 \
    && echo "  PASS  synaptik $cmd" \
    || { echo "  FAIL  synaptik $cmd"; exit 1; }
done

echo
echo "=== JSON output validity ==="
for cmd in phenotype splicing report; do
  case $cmd in
    phenotype) a="'SYN1 puncta density'";;
    *) a="";;
  esac
  eval "python3 -m synaptik.cli.main $cmd $a --db '$TMPDB' --json" \
    | python3 -c "import json,sys; json.load(sys.stdin)" \
    && echo "  PASS  $cmd --json"
done

rm -f "$TMPDB"*
echo
echo "All tests passed."
