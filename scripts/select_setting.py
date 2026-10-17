#!/usr/bin/env python3
"""
select_setting.py
Chooses a QC setting from the summary tables: the Pareto front of genotype miss rate against genotype
concordance, and its knee point by weighted Kneedle for each weight ratio W (weight on concordance : weight
on miss rate). The algorithm is that of notebooks/02_pareto_knee.ipynb, moved into a script.

The weighted Kneedle uses the difference curve  D_W = w_y * y_norm - w_x * x_norm  (w_x = 1, w_y = W); with
W = 1 it is the original Kneedle (Satopaa et al. 2011, as implemented in kneed). WeightedKneeLocator below
gives the same knee as kneed patched with patches/kneed-0.8.5-weighted.patch, on a stock kneed.

Inputs : merged_summary_all.csv, merged_summary_15x.csv, merged_summary_30x.csv (scripts/export_summary_tables.py)
Outputs: pareto_front_<group>.tsv  the non-dominated settings of each group, by miss rate
         knee_by_weight.tsv        the knee of each group and weight, with the metrics of that setting in
                                   all samples, 15x and 30x
Usage  : python3 select_setting.py --tables_dir results/summary --out_dir results/selection
"""
import argparse
import os
from fractions import Fraction

import numpy as np
import pandas as pd
from kneed import KneeLocator
from scipy.signal import argrelextrema

GROUPS = {'ALL': 'all', '15X': '15x', '30X': '30x'}
KEYS = ['DP', 'GQ', 'LAF', 'HAF']
X_COL, Y_COL = 'GENOTYPE_MISS_RATE', 'GENOTYPE_CONCORDANCE'      # cost (minimize), benefit (maximize)
METRICS = ['GENOTYPE_CONCORDANCE', 'GENOTYPE_MISS_RATE', 'FALSE_POSITIVE_RATE', 'FALSE_NEGATIVE_RATE']
DEFAULT_WEIGHTS = '1/7,1/6,1/5,1/4,1/3,1/2,1,2,3,4,5,6'


class WeightedKneeLocator(KneeLocator):
    """Kneedle on the weighted difference curve weight_y * y_norm - weight_x * x_norm.

    kneed computes the difference curve, its local extrema and thresholds in __init__ and then calls
    find_knee(); this class recomputes those three with the weights just before the search, so the knee is
    the one the patched kneed returns."""

    def __init__(self, x, y, weight_x=1.0, weight_y=1.0, **kwargs):
        self._weight_x = weight_x
        self._weight_y = weight_y
        super().__init__(x, y, **kwargs)

    def find_knee(self):
        self.y_difference = self._weight_y * self.y_normalized - self._weight_x * self.x_normalized
        self.maxima_indices = argrelextrema(self.y_difference, np.greater_equal)[0]
        self.x_difference_maxima = self.x_difference[self.maxima_indices]
        self.y_difference_maxima = self.y_difference[self.maxima_indices]
        self.minima_indices = argrelextrema(self.y_difference, np.less_equal)[0]
        self.x_difference_minima = self.x_difference[self.minima_indices]
        self.y_difference_minima = self.y_difference[self.minima_indices]
        self.Tmx = self.y_difference_maxima - (self.S * np.abs(np.diff(self.x_normalized).mean()))
        return super().find_knee()


def find_pareto_front(df, metric1_col, metric2_col, maximize_metric1=True, maximize_metric2=False):
    """True for the settings no other setting beats on both metrics (as in notebooks/02_pareto_knee)."""
    m1 = df[metric1_col].values
    m2 = df[metric2_col].values
    is_efficient = np.ones(len(df), dtype=bool)
    for i in range(len(df)):
        if is_efficient[i]:
            better_m1 = m1 >= m1[i] if maximize_metric1 else m1 <= m1[i]
            better_m2 = m2 >= m2[i] if maximize_metric2 else m2 <= m2[i]
            strictly_better_m1 = m1 > m1[i] if maximize_metric1 else m1 < m1[i]
            strictly_better_m2 = m2 > m2[i] if maximize_metric2 else m2 < m2[i]
            is_dominated = (better_m1 & better_m2) & (strictly_better_m1 | strictly_better_m2)
            is_dominated[i] = False
            if np.any(is_dominated):
                is_efficient[i] = False
    return pd.Series(is_efficient, index=df.index)


def knee_point(front, weight_y, weight_x=1.0, S=1.0):
    """The knee of the front (sorted by miss rate) for one weight; None when Kneedle finds none."""
    d = front[[X_COL, Y_COL] + KEYS].drop_duplicates().sort_values(by=X_COL)
    x = d[X_COL].values
    y = d[Y_COL].values
    kl = WeightedKneeLocator(x, y, weight_x=weight_x, weight_y=weight_y, S=S, curve='concave', direction='increasing')
    if kl.knee is None:
        return None
    hits = np.where(x == kl.knee)[0]
    return d.iloc[hits[0]] if len(hits) else None


def parse_weights(text):
    return [Fraction(w.strip()) for w in text.split(',') if w.strip()]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--tables_dir', required=True, help='folder of merged_summary_{all,15x,30x}.csv')
    ap.add_argument('--out_dir', required=True)
    ap.add_argument('--weights', default=DEFAULT_WEIGHTS, help=f'weight ratios W (default {DEFAULT_WEIGHTS})')
    a = ap.parse_args()

    tables = {g: pd.read_csv(os.path.join(a.tables_dir, f'merged_summary_{s}.csv')) for g, s in GROUPS.items()}
    weights = parse_weights(a.weights)
    os.makedirs(a.out_dir, exist_ok=True)

    rows = []
    for group, summary in tables.items():
        front = summary[find_pareto_front(summary, Y_COL, X_COL, True, False)].sort_values(by=X_COL)
        front[KEYS + METRICS].to_csv(os.path.join(a.out_dir, f'pareto_front_{GROUPS[group]}.tsv'), sep='\t',
                                     index=False)
        for w in weights:
            k = knee_point(front, weight_y=float(w))
            row = {'GROUP': group, 'W': str(w), 'W_VALUE': float(w)}
            if k is None:
                rows.append(row)
                continue
            row.update({c: k[c] for c in KEYS})
            for g2, t in tables.items():                         # the setting's metrics in every group
                hit = t[(t.DP == k.DP) & (t.GQ == k.GQ) & np.isclose(t.LAF, k.LAF) & np.isclose(t.HAF, k.HAF)]
                for m in METRICS:
                    row[f'{m}_{g2}'] = hit[m].iloc[0] if len(hit) else np.nan
            rows.append(row)
        print(f'{group}: {len(front)} settings on the Pareto front')

    out = pd.DataFrame(rows)
    out[['DP', 'GQ']] = out[['DP', 'GQ']].astype('Int64')
    out.to_csv(os.path.join(a.out_dir, 'knee_by_weight.tsv'), sep='\t', index=False)
    show = out[out.GROUP == 'ALL'][['W'] + KEYS + ['GENOTYPE_CONCORDANCE_ALL', 'GENOTYPE_MISS_RATE_ALL']]
    print(show.to_string(index=False))


if __name__ == '__main__':
    main()
