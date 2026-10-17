# Changelog

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
