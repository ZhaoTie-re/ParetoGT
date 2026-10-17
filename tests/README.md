# Tests

## Synthetic data (`-profile test`)

`data/` holds synthetic test data written by `make_test_data.py` (deterministic; no real individuals):
- 40 samples, 28 at 15x and 12 at 30x;
- two chromosomes, with 300 WGS sites each, 260 of them also on the array;
- a WGS VCF per chromosome (GT, DP, GQ, AF), an array VCF (GT) and the sample sheet.

The WGS calls are simulated from read counts, so low depth gives missing and wrong calls as in real data.

```bash
nextflow run main.nf -profile test     # all steps, 84 settings, local; results in tests/results/
```

To rewrite the data: `python3 tests/make_test_data.py --bcftools <bcftools>`.

## Comparing two runs

`compare_pickles.py` checks that every result pickle of one run equals its namesake in another: the same keys, and equal tables (values, data types, index and column order).

```bash
python3 tests/compare_pickles.py --reference <earlier>/02.concordance_vmiss_calculator/chr22 \
    --test <new>/02.concordance_vmiss_calculator/chr22
```

Checks made for version 0.2.0 (DSL2, the matrices read once per task, several settings per task):
- **Synthetic data.** The 84 settings of `-profile test`, chromosome 22, against the evaluation script of version 0.1.0 run on the same matrices: identical.
- **First application.** Chromosome 22, 16 settings, run from the prepared matrices of the original run (`runs/real_chr22/run.sh`, `--matrixDir`), against that run's results: the evaluation tables (`02`) and the per-chromosome metrics (`03`) are identical.
