#!/usr/bin/env python3
"""
export_summary_tables.py
Collects the one-row summary of every threshold combination from the genome-wide summary pickles
(results/05.merge_concordance_vmiss_summary/merged._DP*_GQ*_LAF*_HAF*_.summary.pkl) into three tables, one per
sample group: merged_summary_all.csv, merged_summary_15x.csv and merged_summary_30x.csv. The tables hold the
metrics only (concordance, miss rate, false positive and false negative rates, the discordant-call counts and
rates, and the VMISS / SMISS distribution statistics); no per-sample or per-variant rows.

This is the export the notebooks start from (notebooks/02_pareto_knee.ipynb loads the same summaries). It reads
the pickles only; nothing is recomputed.

Usage:
    python3 export_summary_tables.py --summary_dir <results/05.merge_concordance_vmiss_summary> \
        --out_dir <results/summary> [--workers 8]
"""
import argparse
import os
import pickle
from concurrent.futures import ProcessPoolExecutor
from itertools import product

import pandas as pd

# The grid of the pipeline (wgs_gt_qc_calibration.nf)
DP_RANGE = range(1, 31)
GQ_VALUES = [10, 20, 30]
LAF_VALUES = [0.0, 0.1, 0.15, 0.2, 0.25, 0.3]
HAF_VALUES = [1.0, 0.9, 0.85, 0.8, 0.75, 0.7]
GROUPS = {'ALL': 'all', '15X': '15x', '30X': '30x'}


def load_summary(job):
    summary_dir, dp, gq, laf, haf = job
    path = os.path.join(summary_dir, f'merged._DP{dp}_GQ{gq}_LAF{laf}_HAF{haf}_.summary.pkl')
    with open(path, 'rb') as f:
        result_dict = pickle.load(f)
    rows = {}
    for group in GROUPS:
        key = frozenset([f'DP{dp}', f'GQ{gq}', f'LAF{laf}', f'HAF{haf}', group])
        summary = result_dict[key][3].copy()
        summary['PLATFORM'] = group
        rows[group] = summary
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--summary_dir', required=True, help='folder of the merged.*.summary.pkl files')
    ap.add_argument('--out_dir', required=True, help='folder for merged_summary_{all,15x,30x}.csv')
    ap.add_argument('--workers', type=int, default=8)
    a = ap.parse_args()

    jobs = [(a.summary_dir, dp, gq, laf, haf) for dp, gq, laf, haf in product(DP_RANGE, GQ_VALUES, LAF_VALUES, HAF_VALUES)]
    missing = [j for j in jobs if not os.path.exists(os.path.join(
        a.summary_dir, f'merged._DP{j[1]}_GQ{j[2]}_LAF{j[3]}_HAF{j[4]}_.summary.pkl'))]
    if missing:
        raise SystemExit(f'{len(missing)} of {len(jobs)} summary files are missing, e.g. {missing[0][1:]}')

    tables = {group: [] for group in GROUPS}
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for rows in ex.map(load_summary, jobs, chunksize=4):
            for group, row in rows.items():
                tables[group].append(row)

    os.makedirs(a.out_dir, exist_ok=True)
    for group, suffix in GROUPS.items():
        d = pd.concat(tables[group], ignore_index=True).sort_values(['DP', 'GQ', 'LAF', 'HAF']).reset_index(drop=True)
        assert len(d) == len(jobs), (group, len(d))
        out = os.path.join(a.out_dir, f'merged_summary_{suffix}.csv')
        d.to_csv(out, index=False)
        print(f'{out}: {len(d)} rows, {d.shape[1]} columns')


if __name__ == '__main__':
    main()
