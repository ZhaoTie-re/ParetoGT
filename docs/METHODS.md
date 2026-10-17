# Methods

## 1. Truth set and comparison set

The array genotypes of the same individuals serve as the truth. Before the comparison they pass array QC and a frequency filter against a large reference panel of the same population, so that the remaining array sites are high-confidence anchors.

The comparison set is the samples and variants present in both the WGS calls and the array genotypes (`FormatMatrixPrepare`, `scripts/extract_vcf_format.py`):
- the shared samples;
- the shared variant IDs (`chr:pos:ref:alt`), one chromosome at a time.

For each chromosome the step writes sample × variant matrices of the WGS FORMAT fields GT, DP, GQ and AF, and of the array GT. Genotypes are coded 0 (reference homozygote), 1 (heterozygote), 2 (alternate homozygote) or missing.

## 2. Threshold grid

Each combination of four thresholds defines one mask:

| Threshold | Values |
|---|---|
| DP: minimum read depth | 1, 2, …, 30 |
| GQ: minimum genotype quality | 10, 20, 30 |
| LAF: lower bound of the heterozygote allele fraction | 0.0, 0.1, 0.15, 0.2, 0.25, 0.3 |
| HAF: upper bound of the heterozygote allele fraction | 1.0, 0.9, 0.85, 0.8, 0.75, 0.7 |

That makes 30 × 3 × 6 × 6 = 3,240 combinations. The grid is set in `wgs_gt_qc_calibration.nf`.

A WGS genotype is set to missing when any of the following holds (`scripts/evaluate_genotype_concordance_and_vmiss.py`):
- DP < DP threshold;
- GQ < GQ threshold;
- the call is heterozygous and AF < LAF;
- the call is heterozygous and AF > HAF.

A genotype that is missing in the WGS calls stays missing. A missing DP, GQ or AF value does not mask a call, since a comparison with a missing value is false.

## 3. Confusion matrix and metrics

For each combination, chromosome and sample group, the masked WGS genotypes are cross-tabulated against the array genotypes:
- rows: array genotype t ∈ {0, 1, 2};
- columns: WGS genotype c ∈ {0, 1, 2, missing};
- cells: the counts n<sub>tc</sub>.

Sites where the array genotype is missing are left out. The metrics follow Kishikawa et al. (2019) (`scripts/genotype_concordance_vmiss_summary.py`):

| Metric | Definition |
|---|---|
| Genotype concordance | Σ<sub>g</sub> n<sub>gg</sub> / N, where N counts the genotypes called on both platforms (t, c ∈ {0, 1, 2}) |
| Genotype miss rate | Σ<sub>t</sub> n<sub>t,missing</sub> / (N + Σ<sub>t</sub> n<sub>t,missing</sub>) |
| False positive rate | (n<sub>01</sub> + n<sub>02</sub> + n<sub>12</sub>) / N: the WGS call carries more alternate alleles than the array |
| False negative rate | (n<sub>10</sub> + n<sub>20</sub> + n<sub>21</sub>) / N: the WGS call carries fewer alternate alleles than the array |
| Class counts and rates | each n<sub>tc</sub> (for example HOMREF>HET = n<sub>01</sub>) and n<sub>tc</sub> / N |
| VMISS, SMISS | the missing rate of each variant across samples and of each sample across variants; their mean, SD and median |

**Sample groups.** Each metric is computed in all samples (ALL) and separately by sequencing depth (15X, 30X). The depth groups come from the sample sheet (`info/`, column `Target_DP`).

**Merging.** The per-chromosome counts of one combination are added up (`scripts/merge_concordance_results.py`), and the genome-wide metrics are computed from the sums (`MergeConcordanceVmissSummary`). `scripts/export_summary_tables.py` collects the genome-wide metrics of all combinations into one table per sample group.

## 4. Choosing a setting

Stricter thresholds raise concordance and lower the error rates, but they also set more genotypes to missing. No single combination is best on both counts, so the choice is a trade-off between two objectives: maximize concordance and minimize the miss rate (`notebooks/02_pareto_knee.ipynb`).

**Pareto front.** A combination is on the front when no other combination has both:
- a concordance at least as high and a miss rate at least as low;
- and is strictly better on one of the two.

Only combinations on the front are considered.

**Knee point.** The front is a monotonic, concave curve: x is the miss rate and y the concordance. Its knee point is found with Kneedle (Satopää et al., 2011), as implemented in `kneed` (`curve='concave'`, `direction='increasing'`, S = 1):
- x and y are min–max normalized;
- the knee lies where the difference curve stops rising, that is, where further gains in concordance cost more miss rate than they return.

**Weights.** A weight ratio W = w<sub>y</sub> : w<sub>x</sub> (weight on concordance : weight on miss rate) is applied to the difference curve:

  D<sub>W</sub>(x) = w<sub>y</sub> · y<sub>norm</sub> − w<sub>x</sub> · x<sub>norm</sub>

- W = 1 is the original Kneedle.
- W > 1 values concordance more, and the knee moves to stricter settings.
- W < 1 values retention more, and the knee moves to looser settings.

The weights are not part of the released `kneed`. They are added by `patches/kneed-0.8.5-weighted.patch` (three changed lines).

The notebooks scan W from 1/7 to 6. A setting chosen for a whole interval of W is preferred, because it does not depend on one weighting.

**Checks.**
- `notebooks/03_cross_platform.ipynb` compares the metrics of the knee settings between depth groups.
- `notebooks/04_cohort_weighting_check.ipynb` re-weights the 15x and 30x metrics to the depth composition of the whole cohort.

## References

- Kishikawa, T., Momozawa, Y., Ozeki, T., Mushiroda, T., Inohara, H., Kamatani, Y., Kubo, M., & Okada, Y. (2019). Empirical evaluation of variant calling accuracy using ultra-deep whole-genome sequencing data. *Scientific Reports, 9*(1), 1784.
- Satopää, V., Albrecht, J., Irwin, D., & Raghavan, B. (2011). Finding a "Kneedle" in a haystack: Detecting knee points in system behavior. *31st International Conference on Distributed Computing Systems Workshops*, 166–171.
