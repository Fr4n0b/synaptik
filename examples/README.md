# Example input files

**Synthetic data** for trying `synaptik deposit` against the demo knowledge
base. Nothing here is a real measurement.

- `SE.MATS.JC.txt`: an rMATS skipped-exon table (sample 1 = perturbed).
- `phenotypes.csv`: a phenotype table.

```bash
synaptik demo --db data/demo.db
synaptik deposit splicing examples/SE.MATS.JC.txt --perturbed sample1 --db data/demo.db
synaptik deposit phenotypes examples/phenotypes.csv --db data/demo.db
```
