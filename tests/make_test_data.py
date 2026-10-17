#!/usr/bin/env python3
"""
make_test_data.py
Writes the synthetic test data of `-profile test` to tests/data/ (deterministic; no real individuals):

  wgs/chr21.test.vcf.gz, wgs/chr22.test.vcf.gz   WGS calls, FORMAT GT:DP:GQ:AF, 40 samples (28 at 15x, 12 at 30x)
  array/array.test.vcf.gz                        array genotypes (truth), FORMAT GT, 38 samples (36 shared)
  sample_info.csv                                ID,Target_DP

Per chromosome, 300 WGS sites, 260 of them also on the array. True genotypes follow Hardy-Weinberg at a
random allele frequency. A WGS call is simulated from read counts: depth ~ Poisson(target depth x site
effect), alternate reads ~ Binomial(depth, allele fraction), genotype = the most likely of 0/1/2, GQ = the
Phred-scaled ratio of the two most likely genotypes. Low depth therefore gives missing or wrong calls, and
the thresholds matter as in real data. The array has 0.5% missing and 0.05% wrong genotypes.

Usage: python3 make_test_data.py [--out tests/data] [--bcftools bcftools]   (needs numpy and bcftools)
"""
import argparse
import os
import subprocess

import numpy as np

SEED = 2026
N_SITES, N_SHARED = 300, 260
CHROMS = {'chr21': 46709983, 'chr22': 50818468}


def genotype_likelihoods(alt, dp):
    """Binomial likelihoods of 0/1/2 given alternate reads and depth (error 1%)."""
    p = np.array([0.01, 0.5, 0.99])
    with np.errstate(divide='ignore'):
        logl = alt[..., None] * np.log(p) + (dp - alt)[..., None] * np.log(1 - p)
    return logl


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data'))
    ap.add_argument('--bcftools', default='bcftools')
    a = ap.parse_args()
    rng = np.random.default_rng(SEED)
    os.makedirs(os.path.join(a.out, 'wgs'), exist_ok=True)
    os.makedirs(os.path.join(a.out, 'array'), exist_ok=True)

    wgs_samples = [f'S{i:03d}' for i in range(1, 41)]
    depth = {s: (15 if i < 28 else 30) for i, s in enumerate(wgs_samples)}
    array_samples = wgs_samples[:36] + ['A001', 'A002']
    with open(os.path.join(a.out, 'sample_info.csv'), 'w') as f:
        f.write('ID,Target_DP\n')
        for s in wgs_samples:
            f.write(f'{s},{depth[s]}x\n')

    array_records = []
    for chrom, length in CHROMS.items():
        pos = np.sort(rng.choice(np.arange(10_000_000, length - 1000), N_SITES + 20, replace=False))
        bases = np.array(list('ACGT'))
        ref = rng.choice(bases, len(pos))
        alt = np.array([rng.choice(bases[bases != r]) for r in ref])
        ids = [f'{chrom}:{p}:{r}:{x}' for p, r, x in zip(pos, ref, alt)]
        maf = rng.uniform(0.05, 0.5, len(pos))
        n_all = len(wgs_samples) + 2                         # WGS samples, then A001, A002
        truth = rng.binomial(2, maf[:, None], size=(len(pos), n_all))

        # WGS: the first N_SITES sites; the array: N_SHARED of them plus 20 array-only sites
        wgs_idx = np.arange(N_SITES)
        site_effect = rng.gamma(8.0, 1 / 8.0, size=(N_SITES, 1))
        mean_dp = np.array([depth[s] for s in wgs_samples])[None, :] * site_effect
        dp = rng.poisson(mean_dp)
        af_true = np.where(truth[wgs_idx, :40] == 1, rng.beta(20, 20, size=(N_SITES, 40)),
                           np.where(truth[wgs_idx, :40] == 0, 0.005, 0.995))
        alt_reads = rng.binomial(dp, af_true)
        logl = genotype_likelihoods(alt_reads, dp)
        order = np.argsort(-logl, axis=-1)
        gt = order[..., 0]
        best = np.take_along_axis(logl, order[..., :1], -1)[..., 0]
        second = np.take_along_axis(logl, order[..., 1:2], -1)[..., 0]
        gq = np.minimum(99, np.round(10 * (best - second) / np.log(10))).astype(int)

        vcf = os.path.join(a.out, 'wgs', f'{chrom}.test.vcf')
        with open(vcf, 'w') as f:
            f.write('##fileformat=VCFv4.2\n')
            for c, ln in CHROMS.items():
                f.write(f'##contig=<ID={c},length={ln}>\n')
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n'
                    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
                    '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Alternate allele fraction">\n')
            f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(wgs_samples) + '\n')
            for i in wgs_idx:
                cells = []
                for j in range(40):
                    if dp[i, j] == 0:
                        cells.append('./.:0:.:.')
                    else:
                        g = ('0/0', '0/1', '1/1')[gt[i, j]]
                        cells.append(f'{g}:{dp[i, j]}:{gq[i, j]}:{alt_reads[i, j] / dp[i, j]:.3f}')
                f.write(f'{chrom}\t{pos[i]}\t{ids[i]}\t{ref[i]}\t{alt[i]}\t.\tPASS\t.\tGT:DP:GQ:AF\t' +
                        '\t'.join(cells) + '\n')
        subprocess.run([a.bcftools, 'view', '-Oz', '-o', vcf + '.gz', vcf], check=True)
        subprocess.run([a.bcftools, 'index', '-t', '-f', vcf + '.gz'], check=True)
        os.remove(vcf)

        cols = list(range(36)) + [40, 41]                     # S001-S036, A001, A002
        for i in list(rng.choice(N_SITES, N_SHARED, replace=False)) + list(range(N_SITES, N_SITES + 20)):
            g = truth[i, cols].copy()
            wrong = rng.random(len(g)) < 0.0005
            g[wrong] = (g[wrong] + 1) % 3
            miss = rng.random(len(g)) < 0.005
            cells = ['./.' if m else ('0/0', '0/1', '1/1')[x] for x, m in zip(g, miss)]
            array_records.append((chrom, pos[i], ids[i], ref[i], alt[i], cells))

    vcf = os.path.join(a.out, 'array', 'array.test.vcf')
    with open(vcf, 'w') as f:
        f.write('##fileformat=VCFv4.2\n')
        for c, ln in CHROMS.items():
            f.write(f'##contig=<ID={c},length={ln}>\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(array_samples) + '\n')
        for chrom, p, vid, r, x, cells in sorted(array_records, key=lambda t: (t[0], t[1])):
            f.write(f'{chrom}\t{p}\t{vid}\t{r}\t{x}\t.\tPASS\t.\tGT\t' + '\t'.join(cells) + '\n')
    subprocess.run([a.bcftools, 'view', '-Oz', '-o', vcf + '.gz', vcf], check=True)
    subprocess.run([a.bcftools, 'index', '-t', '-f', vcf + '.gz'], check=True)
    os.remove(vcf)
    print(f'test data written to {a.out}')


if __name__ == '__main__':
    main()
