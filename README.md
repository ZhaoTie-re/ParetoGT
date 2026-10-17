# wgs_gt_qc_calibration

Calibrates genotype-level QC thresholds for whole-genome sequencing (WGS) against array genotypes of the same samples.

The tool scans 3,240 combinations of:
- read depth (DP);
- genotype quality (GQ);
- the heterozygote allele-fraction window [LAF, HAF].

For each combination it masks the WGS calls that fail and measures, against the array genotypes:
- concordance;
- miss rate;
- false positive and false negative rates.

It does this in all samples and per sequencing-depth group. A setting is then chosen on the Pareto front of concordance against miss rate, at the knee point found by a weighted Kneedle.

**Status: under development (0.1.0).** The pipeline and the analysis notebooks are those of the first application, now in their own repository (see [CHANGELOG.md](CHANGELOG.md)). Interfaces and file layout may still change.

## How it works

| Step | What it does | Where |
|---|---|---|
| 1 | Shared samples and variants of the WGS and array VCFs; GT, DP, GQ and AF matrices per chromosome | `FormatMatrixPrepare`, `scripts/extract_vcf_format.py` |
| 2 | For each combination: WGS genotypes failing DP or GQ, and heterozygotes outside [LAF, HAF], set to missing; the confusion matrix against the array, VMISS and SMISS | `ConcordanceVmissCalculator`, `scripts/evaluate_genotype_concordance_and_vmiss.py` |
| 3 | Metrics per chromosome | `ConcordanceVmissSummary`, `scripts/genotype_concordance_vmiss_summary.py` |
| 4–5 | Genome-wide merge and metrics, per combination | `MergeConcordanceVmiss`, `MergeConcordanceVmissSummary`, `scripts/merge_concordance_results.py` |
| 6 | One table per sample group (ALL, 15X, 30X) | `scripts/export_summary_tables.py` → `results/summary/` |
| 7 | Pareto front and weighted knee point; threshold scans; depth-group comparison | `notebooks/01`–`04` |

The metric definitions, the grid and the selection rule are in [docs/METHODS.md](docs/METHODS.md). The slides of the first application are in [docs/qc.gt_tuning.v6.rev1.pptx](docs/qc.gt_tuning.v6.rev1.pptx).

## First application

**Data.**
- 1,750 samples with both WGS and array genotypes, sequenced at 15x or 30x;
- 382,521 shared autosomal variants.

**Knee point by weight ratio.** W is the weight on concordance relative to the weight on miss rate. The metrics are for all samples.

| W | DP | GQ | LAF | HAF | Concordance | Miss rate |
|---|---:|---:|---:|---:|---:|---:|
| 1/7 – 1/4 | 8 | 20 | 0.20 | 0.80 | 99.975% | 3.534% |
| 1/3 | 9 | 20 | 0.20 | 0.80 | 99.978% | 4.970% |
| 1/2 | 9 | 20 | 0.25 | 0.80 | 99.979% | 5.259% |
| 1 | 9 | 20 | 0.25 | 0.75 | 99.979% | 5.551% |
| 2 – 6 | 11 | 30 | 0.25 | 0.75 | 99.984% | 13.232% |

**Adopted setting: DP ≥ 8, GQ ≥ 20, heterozygote AF in [0.20, 0.80].** The same setting is the knee for every W from 1/7 to 1/4. At this setting:

| Group | Concordance | Miss rate | False positive rate | False negative rate |
|---|---:|---:|---:|---:|
| All | 99.975% | 3.534% | 0.013% | 0.011% |
| 15x | 99.977% | 3.842% | 0.013% | 0.010% |
| 30x | 99.969% | 2.095% | 0.017% | 0.015% |

The 15x and 30x samples agree on accuracy; they differ in missingness.

All the metrics of every combination are in [`results/summary/`](results/summary/), one table per group.

## Requirements

- Nextflow 22.10.x. The pipeline is written in DSL1, which later releases no longer support. Here it is the conda environment `dsl1`.
- A SLURM cluster. Resources are set as `--rsc` cluster options in `nextflow.config`, and the partition in the `.nf`.
- bcftools (tested with 1.21).
- Python 3 with the packages in [requirements.txt](requirements.txt). The pipeline steps need pandas and numpy; the notebooks need the rest.
- For the notebooks: `kneed` 0.8.5 with the weighted-Kneedle patch:

  ```bash
  pip install kneed==0.8.5
  SITE=$(python -c 'import kneed, os; print(os.path.dirname(os.path.dirname(kneed.__file__)))')
  patch -p1 -d "$SITE" < patches/kneed-0.8.5-weighted.patch
  ```

## Run

```bash
sbatch run_nextflow.sh          # resumes the last run; `sbatch run_nextflow.sh clean` starts afresh
./monitor_status.sh -l          # progress, by process, from trace.txt and the launcher log
```

**Inputs.** Set them on the command line or in the `.nf`:

| Parameter | Meaning |
|---|---|
| `--wgsDir`, `--callVcfSuffix` | WGS VCFs, one per chromosome: `<chr><suffix>` |
| `--trueVcf` | array VCF |
| `--sampleInfo` | sample sheet with the depth group of each sample ([info/README.md](info/README.md)) |

**Outputs.** They go to `results/01`–`05`, about 0.5 TB per run for the full grid. Then:

```bash
python3 scripts/export_summary_tables.py \
    --summary_dir results/05.merge_concordance_vmiss_summary --out_dir results/summary
```

and run the notebooks from `notebooks/`.

## Layout

| Path | Content |
|---|---|
| `wgs_gt_qc_calibration.nf`, `nextflow.config` | the pipeline (DSL1) and its SLURM resources |
| `run_nextflow.sh`, `monitor_status.sh` | launcher (SLURM job) and progress monitor |
| `scripts/` | the pipeline steps and the summary export |
| `notebooks/` | `01_threshold_scans`, `02_pareto_knee`, `03_cross_platform`, `04_cohort_weighting_check` |
| `patches/` | the weighted Kneedle for `kneed` |
| `docs/` | methods and slides |
| `info/` | the sample sheet (not under version control) |
| `results/summary/` | the metrics of every combination, per sample group |

Two kinds of data are kept locally and are not under version control:
- per-sample data: the sample sheet, and the full results of the adopted and other knee settings (`results/selected_settings/`);
- run records (`results/_run_info/`).

## Roadmap

- Port the pipeline to DSL2.
- Make the threshold grid a parameter.
- Move the Pareto and weighted-Kneedle selection from the notebooks into a script that writes the chosen setting and its metrics.
- Keep the per-chromosome intermediates smaller.
- Allow sample groups other than 15x and 30x.

## Author

Zhao Tie ([@ZhaoTie-re](https://github.com/ZhaoTie-re)).
