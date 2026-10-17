#!/usr/bin/env python3
"""
export_summary_tables.py
Collects the one-row summary of every threshold combination from the genome-wide summary pickles
(merged._DP*_GQ*_LAF*_HAF*_.summary.pkl) into three tables, one per sample group: merged_summary_all.csv,
merged_summary_15x.csv and merged_summary_30x.csv, sorted by DP, GQ, LAF and HAF. The tables hold the metrics
only (concordance, miss rate, false positive and false negative rates, the discordant-call counts and rates,
and the VMISS / SMISS distribution statistics); no per-sample or per-variant rows.

The combinations are read from the file names, so any grid works; --expect checks their number. This is the
export the notebooks start from. It reads the pickles only; nothing is recomputed.

Usage:
    python3 export_summary_tables.py --summary_dir results/05.merge_concordance_vmiss_summary \
        --out_dir results/summary [--expect 3240] [--workers 8]
    python3 export_summary_tables.py --summary_files merged.*.summary.pkl --out_dir . [--workers 8]
"""
import argparse
import glob
import os
import pickle
import re
from concurrent.futures import ProcessPoolExecutor

import pandas as pd

GROUPS = {'ALL': 'all', '15X': '15x', '30X': '30x'}
NAME = re.compile(r'^merged\._DP(\d+)_GQ(\d+)_LAF([0-9.]+)_HAF([0-9.]+)_\.summary\.pkl$')


def load_summary(path):
    m = NAME.match(os.path.basename(path))
    dp, gq, laf, haf = m.group(1), m.group(2), m.group(3), m.group(4)
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
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument('--summary_dir', help='folder of the merged.*.summary.pkl files')
    src.add_argument('--summary_files', nargs='+', help='the merged.*.summary.pkl files')
    ap.add_argument('--out_dir', required=True, help='folder for merged_summary_{all,15x,30x}.csv')
    ap.add_argument('--expect', type=int, default=0, help='number of combinations expected (0: no check)')
    ap.add_argument('--workers', type=int, default=8)
    a = ap.parse_args()

    files = a.summary_files or glob.glob(os.path.join(a.summary_dir, 'merged.*.summary.pkl'))
    files = sorted(f for f in files if NAME.match(os.path.basename(f)))
    if not files:
        raise SystemExit('no merged.*.summary.pkl files found')
    if a.expect and len(files) != a.expect:
        raise SystemExit(f'{len(files)} summary files found, {a.expect} expected')

    tables = {group: [] for group in GROUPS}
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for rows in ex.map(load_summary, files, chunksize=4):
            for group, row in rows.items():
                tables[group].append(row)

    os.makedirs(a.out_dir, exist_ok=True)
    for group, suffix in GROUPS.items():
        d = pd.concat(tables[group], ignore_index=True).sort_values(['DP', 'GQ', 'LAF', 'HAF']).reset_index(drop=True)
        out = os.path.join(a.out_dir, f'merged_summary_{suffix}.csv')
        d.to_csv(out, index=False)
        print(f'{out}: {len(d)} rows, {d.shape[1]} columns')


if __name__ == '__main__':
    main()
