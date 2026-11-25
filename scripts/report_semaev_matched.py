#!/usr/bin/env python3
"""Render the matched benchmark figure and descriptive table from saved medians."""
import json
import math
from pathlib import Path
import statistics

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parents[1]/'experiments/semaev-matched-20261124'


def main():
    summary = json.loads((OUT/'summary.json').read_text())
    cases = summary['cases']
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    for cat, label, color, marker in [('yes', 'PDP solution', '#1976a3', 'o'),
                                       ('no', 'No PDP solution', '#ca6b24', 's')]:
        cs = [c for c in cases if c['category']==cat and 'gb_speedup' in c]
        ax.scatter([c['original']['gb_seconds'] for c in cs],
                   [c['split']['gb_seconds'] for c in cs], c=color, marker=marker,
                   s=35, alpha=.8, label=label)
    values=[c[m]['gb_seconds'] for c in cases if 'gb_speedup' in c for m in ('original','split')]
    lo, hi = min(values)/1.4, max(values)*1.4
    ax.plot([lo,hi], [lo,hi], '--', color='#555555', linewidth=1, label='Equal time')
    ax.set(xscale='log', yscale='log', xlim=(lo,hi), ylim=(lo,hi),
           xlabel='Original GB time (seconds)', ylabel='Split GB time (seconds)',
           title='All 33 paired GB medians favor original')
    ax.legend(frameon=False, fontsize=8)
    grid=[(3,1),(3,2),(4,2),(5,2),(5,3),(7,2)]
    table=[]
    for i,(n,ell) in enumerate(grid):
        cs=[c for c in cases if (c['n'],c['ell'])==(n,ell)]
        ratios=[c['gb_speedup'] for c in cs if 'gb_speedup' in c]
        gm=math.exp(statistics.mean(map(math.log,ratios)))
        bx.scatter([i]*len(ratios), ratios, s=25, alpha=.6, color='#1976a3')
        bx.scatter([i], [gm], s=70, marker='_', color='black', zorder=3)
        table.append(f'| {n} | {ell} | {len(ratios)} | {gm:.3f} | {min(ratios):.3f}–{max(ratios):.3f} |')
    bx.axhline(1,color='#555555',linestyle='--',linewidth=1)
    bx.set(xticks=range(len(grid)),xticklabels=[f'{n}, {ell}' for n,ell in grid],
           xlabel='Extension degree n, subspace dimension ell',
           ylabel='Original / split GB time (below 1 favors original)',
           title='No crossover on the frozen grid',yscale='log',ylim=(.06,1.2))
    bx.text(.02,.04,'Points: paired medians; black marks: geometric means',
            transform=bx.transAxes,fontsize=8)
    fig.suptitle('Corrected PDP formulations · SymPy 1.14.0 Buchberger · 3 repetitions',fontsize=12)
    fig.savefig(OUT/'comparison.svg',metadata={'Date':None})
    fig.savefig(OUT/'comparison.png',dpi=150)
    plt.close(fig)
    (OUT/'parameter-table.md').write_text(
        '| n | ell | Complete pairs | Geometric mean original/split GB time | Range |\n'
        '| ---: | ---: | ---: | ---: | ---: |\n'+'\n'.join(table)+'\n')


if __name__=='__main__':
    main()
