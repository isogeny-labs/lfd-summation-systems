#!/usr/bin/env python3
"""The three line and scatter figures of paper/curve-lfd-bounds.tex.

  python3 scripts/plot_ic_lfd_figures.py          # writes paper/figures/lfd-*.pdf and .png

lfd-depth-profile   exact two-summand LFD against the depth, with the proved bounds and the conjectured law
lfd-closure-share   share of the degree-D space reached by the closure V_D, D = 4, 5, 6
lfd-variables       exact LFD against the number mk of chain variables, m = 2, ..., 8

Closure dimensions and the more-summand values are read from the committed records; the two-summand
profile and the three-summand values are the ones of Table `tab:profile` and of the text (sources named
next to each constant).  Colours: the first three slots of a palette checked for colour-vision
deficiency (all pairs); every series also has its own marker shape and a direct label.
"""
import json
from math import comb
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / 'paper' / 'figures'
RUNS = ROOT / 'experiments' / 'ic-bounds-k14' / 'runs'

BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'
INK, INK2, GRID = '#0b0b0b', '#52514e', '#dddcd6'

plt.rcParams.update({
    'font.family': 'cmr10', 'mathtext.fontset': 'cm', 'axes.formatter.use_mathtext': True, 'axes.unicode_minus': False,
    'font.size': 9.5, 'axes.labelsize': 9.5,
    'axes.edgecolor': INK2, 'axes.linewidth': 0.6, 'axes.labelcolor': INK, 'text.color': INK,
    'xtick.color': INK2, 'ytick.color': INK2, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True, 'grid.color': GRID,
    'grid.linewidth': 0.5, 'grid.linestyle': '-', 'legend.frameon': False, 'legend.fontsize': 8,
    'lines.linewidth': 1.6, 'lines.markersize': 5.5, 'pdf.fonttype': 42, 'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.03,
})
RING = dict(markeredgecolor='white', markeredgewidth=1.0)      # surface ring around every marker


def save(fig, name):
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / (name + '.pdf'), metadata={'CreationDate': None, 'ModDate': None})
    fig.savefig(FIG / (name + '.png'), dpi=200)
    plt.close(fig)


def low(D, k):
    """|M_<=D| at depth k for two summands: squarefree monomials of weight <= D in 2k variables."""
    return sum(comb(2 * k, w) for w in range(D + 1))


# -------------------------------------------------------------------------------------------------
# 1. depth profile
# -------------------------------------------------------------------------------------------------
# Exact two-summand values on generic instances (Table tab:profile): 4 for k = 2..5, 5 for 6..9,
# 6 for 10..13, 7 at k = 14 (upper bound at 14 by the degree lift).
EXACT = {k: 4 for k in range(2, 6)} | {k: 5 for k in range(6, 10)} | {k: 6 for k in range(10, 14)} | {14: 7}
# Averaged degree of regularity of Petit, Kosters and Messeng, Table 4: 4.00 (n1 <= 5), 5.00 (6..9), 6.00 (10, 11).
PKM = {k: 4 for k in range(2, 6)} | {k: 5 for k in range(6, 10)} | {10: 6, 11: 6}


def depth_profile():
    fig, ax = plt.subplots(figsize=(6.1, 3.3))
    ks = list(range(2, 17))
    ax.plot(ks, [2 * k + 1 for k in ks], color=ORANGE, marker='s', markersize=4.5, **RING, label=r'$2k+1$: kernel-checked upper bound')
    ax.plot(ks, [k + 2 for k in ks], color=AQUA, marker='^', markersize=5, **RING, label=r'$k+2$: upper bound with a written proof')
    ax.step([k - 0.5 for k in ks] + [ks[-1] + 0.5], [(k + 14) // 4 for k in ks] + [(ks[-1] + 14) // 4], where='post',
            color=INK2, linewidth=1.0, label=r'$\lfloor (k+14)/4\rfloor$: conjectured law')
    kp = sorted(PKM)
    ax.plot(kp, [PKM[k] for k in kp], linestyle='none', marker='o', markersize=9.5, markerfacecolor='none', markeredgecolor=INK2,
            markeredgewidth=0.8, label='degree of regularity reported in [3]')
    ke = sorted(EXACT)
    ax.plot(ke, [EXACT[k] for k in ke], linestyle='none', color=BLUE, marker='o', **RING, label='exact last fall degree')
    ax.annotate(r'$2k+1$', (16, 33), xytext=(6, 0), textcoords='offset points', va='center', fontsize=8, color=INK)
    ax.annotate(r'$k+2$', (16, 18), xytext=(6, 0), textcoords='offset points', va='center', fontsize=8, color=INK)
    ax.text(15.25, 8.3, r'conjectured beyond $k=14$', ha='center', va='bottom', fontsize=8, color=INK2)
    ax.set_xlabel(r'depth $k$ (the factor base is the subgroup of order $2^k$)')
    ax.set_ylabel('last fall degree')
    ax.set_xticks(range(2, 17, 2))
    ax.set_xlim(1.4, 17.6)
    ax.set_ylim(0, 35)
    ax.legend(loc='upper left', handlelength=1.8)
    save(fig, 'lfd-depth-profile')


# -------------------------------------------------------------------------------------------------
# 2. share of the degree-D space reached by the closure
# -------------------------------------------------------------------------------------------------
def rank(name):
    return json.loads((RUNS / name).read_text())['rank']


def closure_share():
    # D = 4: dim V_4 = 4k^2 + 172k - 571 for k = 6..13 (Section "Results"); the whole space for k <= 5.
    d4 = {5: low(4, 5)} | {k: 4 * k * k + 172 * k - 571 for k in range(6, 14)}
    # D = 5: the whole space through k = 9; records of the unit instance over p = 65537, k = 16 from the random curve.
    d5 = {9: low(5, 9)} | {k: rank('unit-p65537-k%d-D5.json' % k) for k in range(10, 16)} | {16: rank('proper-p65537-k16-D5.json')}
    # D = 6: the whole space through k = 13; one deficient depth computed.
    d6 = {13: rank('unit-p65537-k13-D6.json'), 14: rank('unit-p65537-k14-D6.json')}
    assert d5[10] == 18569 and d5[16] == 66073 and d6[13] == low(6, 13) and d6[14] == 448922 and d4[6] == 605
    fig, ax = plt.subplots(figsize=(6.1, 3.0))
    for D, d, col, mk in ((4, d4, BLUE, 'o'), (5, d5, ORANGE, 's'), (6, d6, AQUA, '^')):
        ks = sorted(d)
        ys = [d[k] / low(D, k) for k in ks]
        ax.plot(ks, ys, color=col, marker=mk, **RING, label=r'$D=%d$' % D)
        ax.annotate(r'$D=%d$' % D, (ks[-1], ys[-1]), xytext=(7, 0), textcoords='offset points', va='center', fontsize=8, color=INK)
    for k0 in (6, 10, 14):
        ax.annotate(r'$k=%d$' % k0, (k0, 1.0), xytext=(0, 9), textcoords='offset points', ha='center', fontsize=8, color=INK2)
        ax.plot([k0, k0], [1.0, 1.035], color=INK2, linewidth=0.6, clip_on=False)
    ax.set_xlabel(r'depth $k$')
    ax.set_ylabel(r'$\dim V_D\,/\,|M_{\leq D}|$')
    ax.set_xticks(range(5, 17))
    ax.set_xlim(4.5, 17.3)
    ax.set_ylim(0, 1.0)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.legend(loc='lower left', ncol=3, handlelength=1.8, title='degree of the closure', title_fontsize=8)
    save(fig, 'lfd-closure-share')


# -------------------------------------------------------------------------------------------------
# 3. exact LFD against the number of chain variables
# -------------------------------------------------------------------------------------------------
def variables():
    rows = [json.loads(l) for l in (ROOT / 'experiments/ic-bounds-msum/results.jsonl').read_text().splitlines()]
    cells = sorted({(r['m'], r['k'], r['lfd_exact']) for r in rows if r['m'] >= 4 and 2 <= r['k'] <= r['m'] - 1 and r['lfd_exact']})
    assert all(v == m * k + 1 for m, k, v in cells) and len(cells) == 12
    three = [(6, 7), (9, 7), (12, 7), (15, 7)]             # m = 3: exact 7 at k = 2 (on the bound) and at k = 3, 4, 5
    two = [(2 * k, v) for k, v in sorted(EXACT.items())]
    fig, ax = plt.subplots(figsize=(6.1, 3.5))
    xs = list(range(3, 30))
    ax.plot(xs, [x + 1 for x in xs], color=INK2, linewidth=1.0, label=r'$mk+1$: kernel-checked upper bound')
    ax.plot([m * k for m, k, _ in cells], [v for _, _, v in cells], linestyle='none', color=AQUA, marker='^', markersize=6.5, **RING,
            label=r'$m=4,\ldots,8$ summands, $k\leq m-1$')
    ax.plot([x for x, _ in three], [y for _, y in three], linestyle='none', color=ORANGE, marker='s', **RING, label=r'$m=3$ summands, $k=2,\ldots,5$')
    ax.plot([x for x, _ in two], [y for _, y in two], linestyle='none', color=BLUE, marker='o', **RING, label=r'$m=2$ summands, $k=2,\ldots,14$')
    note = dict(textcoords='offset points', fontsize=8, color=INK2, arrowprops=dict(arrowstyle='-', color=INK2, linewidth=0.5))
    ax.annotate('$m=8,\\ k=3$\n$m=6,\\ k=4$', (24, 25), xytext=(16, -22), ha='left', va='center', **note)
    ax.annotate(r'$m=3,\ k=5$', (15, 7), xytext=(10, 14), ha='left', va='center', **note)
    ax.annotate(r'$m=2,\ k=14$', (28, 7), xytext=(-6, 16), ha='right', va='center', **note)
    ax.annotate(r'$mk+1$', (29, 30), xytext=(4, 0), textcoords='offset points', va='center', fontsize=8, color=INK)
    ax.set_xlabel(r'number of chain variables $mk$')
    ax.set_ylabel('exact last fall degree')
    ax.set_xticks(range(4, 30, 4))
    ax.set_xlim(2.5, 31)
    ax.set_ylim(0, 31)
    ax.legend(loc='upper left', handlelength=1.8)
    save(fig, 'lfd-variables')


if __name__ == '__main__':
    depth_profile()
    closure_share()
    variables()
    print('wrote', ', '.join(sorted(p.name for p in FIG.glob('lfd-*.pdf'))))
