#!/usr/bin/env python3
"""
compare_pickles.py
Checks that two runs give the same results: for every result pickle in --test with a namesake in --reference
(e.g. chr22._DP8_GQ20_LAF0.2_HAF0.8_.pkl or its .summary.pkl), the same keys and, under each key, equal
tables (pandas.testing.assert_frame_equal: values, dtypes, index and column order).

Usage: python3 compare_pickles.py --reference <dir of the earlier run> --test <dir of the new run> [--min 1]
Exit status 1 when a table differs or fewer than --min pairs were compared.
"""
import argparse
import glob
import os
import pickle
import sys

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--reference', required=True)
    ap.add_argument('--test', required=True)
    ap.add_argument('--min', type=int, default=1, help='minimum number of pairs to compare')
    a = ap.parse_args()
    pairs, bad = 0, 0
    for test in sorted(glob.glob(os.path.join(a.test, '*.pkl'))):
        ref = os.path.join(a.reference, os.path.basename(test))
        if not os.path.exists(ref):
            continue
        r, t = pickle.load(open(ref, 'rb')), pickle.load(open(test, 'rb'))
        pairs += 1
        if set(r) != set(t):
            print(f'DIFF  {os.path.basename(test)}: keys differ')
            bad += 1
            continue
        for key in r:
            for i, (x, y) in enumerate(zip(r[key], t[key])):
                try:
                    pd.testing.assert_frame_equal(x, y)
                except AssertionError as e:
                    print(f'DIFF  {os.path.basename(test)} {sorted(key)} table {i}: {str(e).splitlines()[0]}')
                    bad += 1
    print(f'{pairs} pickles compared, {bad} differences')
    sys.exit(1 if bad or pairs < a.min else 0)


if __name__ == '__main__':
    main()
