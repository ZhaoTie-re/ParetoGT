#!/usr/bin/env python3
"""
plot_figures.py
Draws the summary figures from the summary tables (scripts/export_summary_tables.py) and the knee table
(scripts/select_setting.py). Every figure is written as PNG and PDF:

  threshold_scan   concordance and miss rate against the DP threshold, one line per GQ threshold, at the
                   allele-fraction window of the adopted setting (all samples)
  pareto_knee      every setting by miss rate and concordance; the Pareto front; the knee points by weight W
  af_window        false heterozygous calls (array HOMREF, WGS HET) and retained true heterozygotes over the
                   allele-fraction window, at the DP and GQ of the adopted setting, 15x and 30x
  depth_groups     miss rate and concordance against the DP threshold in 15x and 30x samples

The adopted setting is marked: the knee point (all samples) at the weight --adopt_weight (default 1/4, the
choice of the first application: W from 1/7 to 1/4 all give DP 8, GQ 20, AF 0.20-0.80), or --setting.

Usage: python3 plot_figures.py --tables_dir results/summary --knee_table results/selection/knee_by_weight.tsv
           --out_dir results/figures [--adopt_weight 1/4 | --setting 8,20,0.2,0.8]
"""
import argparse
import os
from fractions import Fraction

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.ticker import FuncFormatter, MultipleLocator  # noqa: E402

GROUPS = {'ALL': 'all', '15X': '15x', '30X': '30x'}
GQ_COLORS = {10: '#9ECAE1', 20: '#4292C6', 30: '#08306B'}
GROUP_COLORS = {'ALL': '#525252', '15X': '#D95F02', '30X': '#1B9E77'}
FRONT_COLOR = '#B2182B'
INK = '#222222'

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'], 'font.size': 9,
    'axes.titlesize': 10, 'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8,
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': INK, 'axes.labelcolor': INK,
    'xtick.color': INK, 'ytick.color': INK, 'savefig.bbox': 'tight', 'savefig.dpi': 200, 'pdf.fonttype': 42,
})


def pct(decimals):
    return FuncFormatter(lambda v, _: f'{100 * v:.{decimals}f}%')


def select(t, dp=None, gq=None, laf=None, haf=None):
    m = np.ones(len(t), dtype=bool)
    for col, v in (('DP', dp), ('GQ', gq), ('LAF', laf), ('HAF', haf)):
        if v is not None:
            m &= np.isclose(t[col], v)
    return t[m]


def setting_label(s):
    return f'DP {s[0]:g}, GQ {s[1]:g}, AF {s[2]:.2f}–{s[3]:.2f}'


def adopted_setting(knees, weight, front):
    """The knee point (all samples) at the given weight ratio W; else the knee at the nearest W; else the
    setting with the lowest miss rate on the Pareto front (both with a warning)."""
    k = knees[(knees.GROUP == 'ALL') & knees.DP.notna()]
    if k.empty:
        print('warning: no knee point was found; marking the first setting of the Pareto front')
        return tuple(float(v) for v in front[['DP', 'GQ', 'LAF', 'HAF']].iloc[0])
    w = float(Fraction(weight))
    nearest = k.iloc[[int(np.argmin(np.abs(k.W_VALUE.values - w)))]]
    if not np.isclose(nearest.W_VALUE.iloc[0], w):
        print(f'warning: no knee point at W = {weight}; marking the knee at W = {nearest.W.iloc[0]}')
    return tuple(float(v) for v in nearest[['DP', 'GQ', 'LAF', 'HAF']].iloc[0])


def save(fig, out_dir, name):
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(out_dir, f'{name}.{ext}'))
    plt.close(fig)


def fig_threshold_scan(T, s, out_dir):
    t = select(T['ALL'], laf=s[2], haf=s[3])
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8))
    for ax, col, title, dec in ((axes[0], 'GENOTYPE_CONCORDANCE', 'Genotype concordance', 3),
                                (axes[1], 'GENOTYPE_MISS_RATE', 'Genotype miss rate', 0)):
        for gq, color in GQ_COLORS.items():
            d = t[t.GQ == gq].sort_values('DP')
            ax.plot(d.DP, d[col], color=color, lw=1.6, label=f'GQ ≥ {gq}')
        p = select(t, dp=s[0], gq=s[1])
        ax.scatter(p.DP, p[col], s=36, facecolor='white', edgecolor=FRONT_COLOR, lw=1.5, zorder=5,
                   label='adopted' if ax is axes[0] else None)
        ax.set_title(title, loc='left', fontweight='bold')
        ax.set_xlabel('DP threshold')
        ax.yaxis.set_major_formatter(pct(dec))
        ax.set_xlim(0, 31)
        ax.grid(axis='y', color='#E5E5E5', lw=0.6)
    axes[0].legend(frameon=False, loc='lower right')
    fig.suptitle(f'All samples; heterozygote AF window {s[2]:.2f}–{s[3]:.2f}', x=0.01, ha='left', fontsize=8,
                 color='#555555')
    fig.tight_layout()
    save(fig, out_dir, 'threshold_scan')


def fig_pareto_knee(T, fronts, knees, s, out_dir):
    t = T['ALL']
    front = fronts['ALL']
    k = knees[(knees.GROUP == 'ALL') & knees.DP.notna()].copy()
    k['key'] = list(zip(k.DP, k.GQ, k.LAF, k.HAF))
    xmax = max(0.2, 1.3 * k.GENOTYPE_MISS_RATE_ALL.max())
    near = t[t.GENOTYPE_MISS_RATE <= xmax]
    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    ax.scatter(near.GENOTYPE_MISS_RATE, near.GENOTYPE_CONCORDANCE, s=6, color='#C8C8C8', lw=0,
               label='dominated settings')
    f = front[front.GENOTYPE_MISS_RATE <= xmax]
    ax.plot(f.GENOTYPE_MISS_RATE, f.GENOTYPE_CONCORDANCE, color=FRONT_COLOR, lw=1.2, marker='o', ms=2.2,
            label='Pareto front')
    groups = sorted(k.groupby('key', sort=False), key=lambda kv: kv[1].GENOTYPE_MISS_RATE_ALL.iloc[0])
    key_lines = []
    for n, (key, g) in enumerate(groups, start=1):
        ws = list(g.W)
        wtxt = ws[0] if len(ws) == 1 else f'{ws[0]} to {ws[-1]}'
        x, y = g.GENOTYPE_MISS_RATE_ALL.iloc[0], g.GENOTYPE_CONCORDANCE_ALL.iloc[0]
        chosen = np.allclose(key, s)
        color = FRONT_COLOR if chosen else INK
        ax.scatter([x], [y], marker='D', s=40 if chosen else 26, color=color, zorder=6, edgecolor='white', lw=0.6)
        ax.annotate(str(n), (x, y), xytext=(-3, 5), textcoords='offset points', ha='right', fontsize=8,
                    color=color, fontweight='bold')
        key_lines.append((f'{n}  W = {wtxt}: {setting_label(key)}' + ('  (adopted)' if chosen else ''), color))
    ylo = near.GENOTYPE_CONCORDANCE.quantile(0.02)
    yhi = t.GENOTYPE_CONCORDANCE.max()
    ax.set_ylim(ylo, yhi + 0.12 * (yhi - ylo))
    ax.set_xlim(0, xmax)
    ax.xaxis.set_major_locator(MultipleLocator(0.02))
    ax.xaxis.set_major_formatter(pct(0))
    ax.yaxis.set_major_formatter(pct(3))
    ax.set_xlabel('Genotype miss rate')
    ax.set_ylabel('Genotype concordance')
    ax.set_title('Pareto front and knee points by weight W (all samples)', loc='left', fontweight='bold')
    ax.legend(frameon=False, loc='upper left')
    ax.grid(color='#EFEFEF', lw=0.6)
    for i, (line, color) in enumerate(key_lines):
        ax.text(0.98, 0.06 + 0.065 * (len(key_lines) - 1 - i), line, transform=ax.transAxes, ha='right',
                fontsize=7.5, color=color, fontweight='bold' if color == FRONT_COLOR else 'normal')
    fig.tight_layout()
    save(fig, out_dir, 'pareto_knee')


def fig_af_window(T, s, out_dir):
    lafs = sorted(T['ALL'].LAF.unique())
    hafs = sorted(T['ALL'].HAF.unique(), reverse=True)
    panels = [('HOMREF>HET_RATE', 'False heterozygous calls (array HOMREF, WGS HET)', 'Reds', 4),
              ('HET>HET_RATE', 'Heterozygotes called correctly (array HET, WGS HET)', 'Blues', 2)]
    fig = plt.figure(figsize=(6.8, 6.6))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.06], left=0.1, right=0.9, top=0.86, bottom=0.08,
                          wspace=0.18, hspace=0.42)
    for r, (col, title, cmap, dec) in enumerate(panels):
        mats = {}
        for g in ('15X', '30X'):
            d = select(T[g], dp=s[0], gq=s[1])
            mats[g] = d.pivot_table(index='LAF', columns='HAF', values=col).reindex(index=lafs, columns=hafs)
        vmin = min(m.values.min() for m in mats.values())
        vmax = max(m.values.max() for m in mats.values())
        row_axes = []
        for c, g in enumerate(('15X', '30X')):
            ax = fig.add_subplot(gs[r, c])
            row_axes.append(ax)
            m = mats[g]
            im = ax.imshow(m.values, cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')
            for i in range(m.shape[0]):
                for j in range(m.shape[1]):
                    v = m.values[i, j]
                    dark = (v - vmin) / (vmax - vmin + 1e-30) > 0.6
                    ax.text(j, i, f'{100 * v:.{dec - 1}f}', ha='center', va='center', fontsize=6.5,
                            color='white' if dark else INK)
            ax.set_xticks(range(len(hafs)), [f'{h:g}' for h in hafs])
            ax.set_yticks(range(len(lafs)), [f'{lv:g}' for lv in lafs] if c == 0 else [])
            if r == 1:
                ax.set_xlabel('HAF (upper AF bound)')
            if c == 0:
                ax.set_ylabel('LAF (lower AF bound)')
            ax.set_title(f'{g.lower()} samples', loc='left', fontsize=9)
            ax.tick_params(length=0)
            ax.spines[:].set_visible(False)
            if s[2] in lafs and s[3] in hafs:
                ax.add_patch(plt.Rectangle((hafs.index(s[3]) - 0.5, lafs.index(s[2]) - 0.5), 1, 1, fill=False,
                                           ec=FRONT_COLOR, lw=1.6))
        cb = fig.colorbar(im, cax=fig.add_subplot(gs[r, 2]))
        cb.ax.yaxis.set_major_formatter(pct(dec))
        cb.outline.set_visible(False)
        pos = row_axes[0].get_position()
        fig.text(pos.x0, pos.y1 + 0.045, title, fontsize=9.5, fontweight='bold', ha='left')
    fig.text(0.1, 0.965, f'DP {s[0]:g}, GQ {s[1]:g}. Cells: percent of the genotypes called on both platforms; '
             'red box: adopted window', fontsize=8, color='#555555', ha='left')
    save(fig, out_dir, 'af_window')


def fig_depth_groups(T, s, out_dir):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8))
    for ax, col, title, dec in ((axes[0], 'GENOTYPE_MISS_RATE', 'Genotype miss rate', 0),
                                (axes[1], 'GENOTYPE_CONCORDANCE', 'Genotype concordance', 3)):
        for g in ('15X', '30X', 'ALL'):
            d = select(T[g], gq=s[1], laf=s[2], haf=s[3]).sort_values('DP')
            d = d[d.DP <= 15]
            ax.plot(d.DP, d[col], color=GROUP_COLORS[g], lw=1.6 if g != 'ALL' else 1.0,
                    ls='-' if g != 'ALL' else '--', label={'ALL': 'all', '15X': '15x', '30X': '30x'}[g])
            p = select(d, dp=s[0])
            ax.scatter(p.DP, p[col], s=22, color=GROUP_COLORS[g], zorder=5)
        ax.axvline(s[0], color='#BBBBBB', lw=0.8, ls=':')
        ax.set_title(title, loc='left', fontweight='bold')
        ax.set_xlabel('DP threshold')
        ax.yaxis.set_major_formatter(pct(dec))
        ax.set_xlim(0.5, 15.5)
        ax.grid(axis='y', color='#E5E5E5', lw=0.6)
    axes[0].legend(frameon=False, loc='upper left')
    fig.suptitle(f'GQ {s[1]:g}, heterozygote AF window {s[2]:.2f}–{s[3]:.2f}; dotted line: adopted DP',
                 x=0.01, ha='left', fontsize=8, color='#555555')
    fig.tight_layout()
    save(fig, out_dir, 'depth_groups')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--tables_dir', required=True)
    ap.add_argument('--knee_table', required=True)
    ap.add_argument('--out_dir', required=True)
    ap.add_argument('--adopt_weight', default='1/4', help='W whose knee point is the adopted setting (default 1/4)')
    ap.add_argument('--setting', default='', help='DP,GQ,LAF,HAF to mark instead (overrides --adopt_weight)')
    a = ap.parse_args()
    T = {g: pd.read_csv(os.path.join(a.tables_dir, f'merged_summary_{s}.csv')) for g, s in GROUPS.items()}
    knees = pd.read_csv(a.knee_table, sep='\t')
    fronts = {g: pd.read_csv(os.path.join(os.path.dirname(a.knee_table), f'pareto_front_{s}.tsv'), sep='\t')
              for g, s in GROUPS.items()}
    s = (tuple(float(v) for v in a.setting.split(',')) if a.setting
         else adopted_setting(knees, a.adopt_weight, fronts['ALL']))
    os.makedirs(a.out_dir, exist_ok=True)
    fig_threshold_scan(T, s, a.out_dir)
    fig_pareto_knee(T, fronts, knees, s, a.out_dir)
    fig_af_window(T, s, a.out_dir)
    fig_depth_groups(T, s, a.out_dir)
    print(f'figures written to {a.out_dir} (adopted setting: {setting_label(s)})')


if __name__ == '__main__':
    main()
