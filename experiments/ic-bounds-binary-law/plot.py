#!/usr/bin/env python3
"""Figure: LFD against k + e (e = r_eff - m_eff) and against log2|Z|, affine class, existing data
plus the preregistered systems (score.json, if present).  Writes lfd_vs_kplus_e.{svg,png}."""
import json, math, os, sys, collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from analysis import load  # noqa: E402


def main():
    S, _ = load()
    A = [s for s in S if s['cls'] == 'affine']
    pre = []
    p = os.path.join(HERE, 'score.json')
    if os.path.exists(p):
        pre = [d for d in json.load(open(p))['systems'] if d['cls'] == 'affine']
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, xf, xl in ((axes[0], lambda s: s['k'] + s['e'], 'k + e   (e = r_eff - m_eff, log2 of expected |Z|)'),
                       (axes[1], lambda s: s['e'], 'e alone')):
        cnt = collections.Counter((xf(s), s['lfd']) for s in A)
        xs, ys, sz = zip(*[(x, y, c) for (x, y), c in cnt.items()])
        ax.scatter(xs, ys, s=[8 + 3 * math.sqrt(c) * 4 for c in sz], color='#4a6fa5', alpha=0.55,
                   label='existing certified systems (area ~ count)')
        if pre:
            c2 = collections.Counter((xf(d), d['lfd']) for d in pre)
            ax.scatter([x for x, _ in c2], [y for _, y in c2], marker='x', s=70, color='#c0392b',
                       label='preregistered check')
        ax.set_xlabel(xl)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel('certified exact LFD')
    for t in (3, 9, 19):
        axes[0].axvline(t - 0.5, color='grey', ls=':', lw=1)
    axes[0].set_title('law B thresholds k+e >= 3, 9, 19 (dotted)')
    axes[1].set_title('law A: e only')
    axes[0].legend(loc='upper left', fontsize=8)
    fig.suptitle('Two-summand binary/Koblitz systems, affine class (x_R != 0)', fontsize=11)
    fig.tight_layout()
    for ext in ('svg', 'png'):
        fig.savefig(os.path.join(HERE, 'lfd_vs_kplus_e.' + ext), dpi=110)


if __name__ == '__main__':
    main()
