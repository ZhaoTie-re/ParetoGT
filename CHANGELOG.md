# Changelog

## Unreleased

**Docs**
- README:
  - states the purpose of the software, following the slides: depth-dependent genotypes, bias in mixed-depth case–control cohorts, what each threshold removes, and the accuracy-retention trade-off;
  - shows how it works with the grid-search figure of the slides (`docs/grid_search.png`), in place of a flowchart.
- `docs/METHODS.md` opens with the purpose.
- The slides as PDF (`docs/qc.gt_tuning.v6.rev1.pdf`), rendered with metric-compatible substitutes for the slide fonts.

## 0.2.0 (2026-10-18)

Renamed **ParetoGT** (repository and folder `ParetoGT`, formerly `wgs_gt_qc_calibration`).

**Pipeline**
- Rewritten in Nextflow DSL2 (`main.nf`, Nextflow 25.04 or later), lint-clean in strict syntax.
- The threshold grid, chromosomes and batch sizes are parameters.
- Default resources for SLURM (`--rsc`), with `local` and `test` profiles.

**Faster evaluation, same results**
- Each `EVALUATE` task reads its chromosome's matrices once and evaluates a batch of settings (`--combosPerTask`, default 36). The first version re-read the matrices for every block of 1,000 variants and ran one setting per task.
- The per-chromosome metrics are computed in the same task.
- Merging and the genome-wide metrics share one task per batch.
- Task count for the full grid: 2,095, down from 149,062.
- The metric code is unchanged. A failing block now stops the task; it used to be dropped with an empty result.

**New**
- `--matrixDir`: start from the matrices of an earlier run, without reading the VCFs again.
- `SELECT_SETTING` (`scripts/select_setting.py`): the Pareto front and weighted Kneedle of `notebooks/02`, as a step. The weighting runs on the released `kneed`, and gives the same knees as the patched one.
- `FIGURES` (`scripts/plot_figures.py`): threshold scan, Pareto front with knee points, allele-fraction window, depth groups.
- `EXPORT_TABLES` reads the settings from the file names, so any grid works.

**Inputs**
- The WGS and array VCFs are now required parameters, with no defaults; the VCFs of the first application no longer exist at their old paths.

**Tests**
- `tests/`: synthetic data and its generator, `-profile test`, and `compare_pickles.py`.
- Checks:
  - the evaluation tables of 84 synthetic settings equal those of version 0.1.0's script;
  - those of chromosome 22 of the first application (16 settings, from its matrices) equal the original run's.

**Docs**
- README with the figures of the first application; methods updated.
- `results/selection/` and `results/figures/` added.
- The slides no longer name individual samples: two sample IDs in footnotes, on slides 4 and 6, became "two samples".

## 0.1.0 (2026-10-18)

First version in its own repository, imported from the analysis folder `cteph_agp3k.v5/tuning.concordance.v5`. Nothing was re-run.
- The pipeline, scripts and notebooks are those of the v5 run. Code comments and docstrings were translated into English, and the logic is unchanged.
- Paths:
  - the pipeline locates its scripts, sample sheet and results through `projectDir`;
  - the input VCF names and the sample sheet are parameters;
  - the notebooks read the results relative to the repository.
- Script names lost their revision suffixes:
  - `*.rev3.py` became `evaluate_genotype_concordance_and_vmiss.py`;
  - `*.summary.rev1.py` became `genotype_concordance_vmiss_summary.py`;
  - `merge_concordance_results.rev1.py` became `merge_concordance_results.py`.
- Notebooks renamed `01_`–`04_`, each opening with what it shows and where its saved outputs come from.
- New:
  - `scripts/export_summary_tables.py`: one table per sample group from the genome-wide summaries; the notebook export it replaces is commented out in `02`;
  - `patches/kneed-0.8.5-weighted.patch`: the weighted Kneedle the notebooks use, until now a local edit of `kneed`;
  - `docs/METHODS.md`;
  - `docs/qc.gt_tuning.v6.rev1.pptx`, the slides.
- Results:
  - `results/summary/` re-exported from the v5 run's genome-wide summaries;
  - the summary tables previously kept next to the notebooks came from an earlier run, and were not carried over.
