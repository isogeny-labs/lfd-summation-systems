"""Figures of paper/curve-lfd-bounds.tex built from repository records only.

Writes five vector PDFs into paper/figures/ (the pipeline diagram is TikZ in the .tex):
  curve-sweep.pdf    exact LFD histograms of the seeded sweep, with the bound 2k+1 on its own axis
                     (redraws the provisional plot of scripts/build_ic_curve_paper.py, which runs first)
  lfd-depth.pdf      two-summand LFD against the depth k, with bounds and the conjecture
  closure-dims.pdf   closure dimensions d_D(k) against |M_<=D|, and the codimension at depth 4D-10
  msummands.pdf      LFD against deg g for m = 3, 4, 5 summands, with the kernel bound mk+1
  binary-deep.pdf    certified binary LFD at k = 11..16 by coordinate-space class

Deterministic: no randomness except a fixed-seed jitter, and PDF metadata dates removed.
Every plotted number is read from the files named in SOURCES, and the values the paper quotes are
asserted, so the build stops if a record and the text disagree.
"""
from math import comb
from pathlib import Path
import json, re

ROOT = Path(__file__).resolve().parents[1]
EXP = ROOT / 'experiments'
SOURCES = {
    'sweep': 'experiments/ic-curve-sweep/summary.json',
    'depth8': 'experiments/ic-chain-depth8/summary.json',
    'depth9': 'experiments/ic-chain-depth9/summary.json',
    'k14runs': 'experiments/ic-bounds-k14/runs/',
    'k14results': 'experiments/ic-bounds-k14/RESULTS.md',
    'depthstep': 'experiments/ic-depth-step/growth.json',
    'd5runs': 'experiments/ic-bounds-d5/runs/',
    'msum': 'experiments/ic-bounds-msummands/records/ and pod-records/*/out/',
    'binary': 'experiments/ic-bounds-binary-deep/outcomes-2026-11-12.json',
}

# Okabe-Ito colour-blind-safe palette
BLACK, ORANGE, SKY, GREEN, BLUE, VERMIL, PURPLE, GREY = (
    '#000000', '#E69F00', '#56B4E9', '#009E73', '#0072B2', '#D55E00', '#CC79A7', '#999999')


def load(rel):
    return json.loads((ROOT / rel).read_text())


def run(name):
    return load(f'experiments/ic-bounds-k14/runs/{name}.json')


def m_le(k, D):
    """Number of squarefree monomials of degree <= D in 2k chain variables."""
    return sum(comb(2 * k, j) for j in range(D + 1))


def style(plt):
    plt.rcParams.update({
        'font.family': 'serif', 'font.serif': ['DejaVu Serif'], 'mathtext.fontset': 'cm',
        'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 8.5, 'legend.fontsize': 7,
        'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5, 'axes.linewidth': 0.6,
        'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'lines.linewidth': 1.1,
        'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
        'svg.hashsalt': 'lfd'})


def save(fig, path):
    fig.savefig(path, metadata={'CreationDate': None, 'ModDate': None, 'Creator': None})


# ---------------------------------------------------------------- data

def two_summand_profile():
    """Exact LFD of the generic two-summand prime chain systems per depth, from the records."""
    exact = {}
    sweep = load(SOURCES['sweep'])
    for g in sweep['groups']:
        if g['kind'] == 'prime':
            exact[g['k']] = max(exact.get(g['k'], 0), g['max_exact_lfd'])
    for key, k in (('depth8', '8'), ('depth9', '9')):
        lfd = load(SOURCES[key])['by_depth'][k]['lfd']
        assert list(lfd) == ['5'], lfd
        exact[int(k)] = 5
    for k in (10, 11, 12, 13):
        d5, d6 = run(f'unit-p65537-k{k}-D5'), run(f'unit-p65537-k{k}-D6')
        assert d5['closed'] and not d5['contains_one'] and d6['contains_one'] and d5['points'] == 0
        exact[k] = 6
    d6 = run('unit-p65537-k14-D6')
    assert d6['closed'] and not d6['contains_one']          # LFD >= 7 at k = 14 (certified)
    lift = run('lift-unit-p65537-k14-D6to7')
    assert lift['contains_one'] and lift['dim_VD1'] == 1683218  # LFD = 7 by the degree lift
    assert [exact[k] for k in range(2, 14)] == [4] * 4 + [5] * 4 + [6] * 4, exact
    return exact


def closure_dims():
    rows = [r for r in load(SOURCES['depthstep'])['rows'] if (r['prime'], r['a'], r['b']) == (12289, 1, 0)]
    d4 = {r['k']: r['dim_V4'] for r in rows}
    assert all(r['solutions'] == 0 and r['monomials_le4'] == m_le(r['k'], 4) for r in rows)
    for k in range(6, 11):
        assert d4[k] == 4 * k * k + 172 * k - 571
    d5 = {}
    for name in ('k6-p12289-s10-D5-split', 'k8-p12289-s10-D5-split-full', 'k9-p12289-s10-D5-split-full'):
        k = int(name[1])
        prof = load(f'experiments/ic-bounds-d5/runs/{name}.json')['profile']
        d5[k] = prof[-1]
        assert d5[k] == m_le(k, 5)
    prof5 = {}
    for k in range(10, 16):
        r = run(f'unit-p65537-k{k}-D5')
        assert r['closed'] and r['columns'] == m_le(k, 5)
        prof5[k] = r['profile']
    prof5[16] = run('proper-p65537-k16-D5')['profile']       # random curve, seed 12
    for k, p in prof5.items():
        d5[k] = p[-1]
    d6 = {}
    for k in range(10, 15):
        r = run(f'unit-p65537-k{k}-D6')
        assert r['closed'] and r['columns'] == m_le(k, 6)
        d6[k] = r['rank']
    assert [d5[k] for k in range(10, 17)] == [18569, 25033, 32057, 39657, 47849, 56649, 66073]
    assert d6[14] == 448922
    # recursion p_j(k+1) = p_j(k) + 2 p_{j-1}(k) + p_{j-2}(k), predicted top value from the previous profile
    rec = {}
    for k in range(10, 16):
        p = prof5[k]
        pred = [p[j] + 2 * (p[j - 1] if j >= 1 else 0) + (p[j - 2] if j >= 2 else 0) for j in range(6)]
        assert pred == prof5[k + 1], (k, pred, prof5[k + 1])
        rec[k + 1] = pred[-1]
    # codimension of V_D in M_<=D at the first deficient depth 4D - 10
    codim = {4: m_le(6, 4) - d4[6], 5: m_le(10, 5) - d5[10], 6: m_le(14, 6) - d6[14]}
    text = (ROOT / SOURCES['k14results']).read_text()
    stated = re.search(r'first deficient depth 4D - 10: ([0-9, ]+) for D = 3, 4, 5, 6', text).group(1)
    stated = [int(s) for s in stated.replace(', ', ';').replace(',', '').split(';')]
    codim[3] = stated[0]
    assert codim[3] == m_le(2, 3) - (2 * 2 + 1)              # dim V_3(2) = 2k + 1 = 5
    assert [codim[D] for D in (3, 4, 5, 6)] == stated == [10, 189, 3131, 50256], (codim, stated)
    return d4, d5, d6, rec, codim


def msummand_points():
    pts = []
    for f in sorted((EXP / 'ic-bounds-msummands/records').glob('P-m*.json')):
        r = json.loads(f.read_text())
        if r.get('exact_lfd') is not None:
            pts.append((r['m'], r['k'], r['generator_degree'], r['exact_lfd'], r['case']))
    for f in sorted((EXP / 'ic-bounds-msummands/pod-records').glob('*/out/GP-m*.json')):
        r = json.loads(f.read_text())
        m, k = map(int, re.match(r'GP-m(\d)k(\d)', f.name).groups())
        if r.get('exact_lfd') is not None:
            pts.append((m, k, r['generator_degree'], r['exact_lfd'], f.stem))
    assert all(lfd == deg + 1 for _, _, deg, lfd, _ in pts), pts
    cells = sorted({(m, k, deg, lfd) for m, k, deg, lfd, _ in pts})
    assert cells == [(3, 2, 6, 7), (3, 3, 6, 7), (3, 4, 6, 7), (4, 2, 8, 9), (4, 3, 12, 13), (4, 4, 12, 13),
                     (5, 2, 10, 11), (5, 3, 15, 16)], cells
    count = {c: sum(1 for p in pts if p[:4] == c) for c in cells}
    return cells, count


def binary_rows():
    rows = load(SOURCES['binary'])['rows']
    assert len(rows) == 57 and all(r['status'].startswith('exact') for r in rows)

    def cls(r):
        if r['base'] == 'subfield':
            return 'subfield'
        if r['base'] == 'trace-zero':                 # hyperplane k = n - 1
            assert r['k'] == r['n'] - 1
            return 'trace-zero'
        if r['n'] == 2 * r['k'] - 1:
            return 'n=2k-1'
        return 'other'
    out = [(r['k'], r['lfd'], cls(r), r['Z']) for r in rows]
    tz = [(k, l) for k, l, c, _ in out if c == 'trace-zero']
    assert sorted(tz) == [(11, 5), (11, 5), (11, 5), (12, 6), (12, 6), (13, 6), (13, 6)], tz
    assert all(l == 4 and z <= 4 for k, l, c, z in out if c == 'n=2k-1') and sum(c == 'n=2k-1' for _, _, c, _ in out) == 11
    assert sorted((k, l) for k, l, c, _ in out if c == 'subfield') == [(10, 6), (10, 6), (12, 7)]
    return out


# ---------------------------------------------------------------- figures

def fig_lfd_depth(plt, out):
    exact = two_summand_profile()
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.1, 2.95), layout='constrained',
                                 gridspec_kw={'width_ratios': [1.25, 1]})
    ks = list(range(2, 21))
    conj = [(k + 14) // 4 for k in ks]
    ax.step(ks, conj, where='mid', color=ORANGE, lw=1.2, ls='--', label=r'conjecture $\lfloor(k+14)/4\rfloor$')
    ax.plot(ks, [k + 2 for k in ks], color=BLUE, lw=1.0, label=r'kernel upper bound $k+2$')
    ax.axhline(4, color=GREY, lw=0.8, ls=':', label=r'kernel lower bound $4$ (unit, $3t^2+a\neq0$)')
    kx = sorted(k for k in exact if k <= 13)
    ax.plot(kx, [exact[k] for k in kx], 'o', color=BLACK, ms=4, label='exact (certified)', zorder=5)
    ax.plot([14], [7], 's', mfc='white', mec=VERMIL, mew=1.3, ms=5.5, zorder=6,
            label=r'$k=14$: $\geq7$ certified, $=7$ by the lift')
    ax.annotate('', xy=(15, 7.9), xytext=(15, 7.0), arrowprops=dict(arrowstyle='->', color=VERMIL, lw=1))
    ax.plot([15], [7], '_', color=VERMIL, ms=7, mew=1.5)
    ax.text(15.3, 7.05, r'$\geq7$', color=VERMIL, fontsize=7, va='bottom')
    ax.set(xlabel='depth $k$', ylabel='last fall degree', xlim=(1.5, 20.5), ylim=(2.5, 12.5),
           xticks=range(2, 21, 2), yticks=range(3, 13))
    ax.grid(alpha=0.15, lw=0.5)
    ax.set_title('(a) measured range')
    kk = list(range(2, 81))
    bx.plot(kk, [k + 2 for k in kk], color=BLUE, lw=1.0, label=r'$k+2$ (kernel)')
    bx.step(kk, [(k + 14) // 4 for k in kk], where='mid', color=ORANGE, ls='--', lw=1.1, label='conjecture')
    kg = [k for k in kk if k >= 38]
    bx.step(kg, [(k - 10) // 14 + 1 for k in kg], where='post', color=GREEN, lw=1.4,
            label=r'written growth bound $\lfloor(k-10)/14\rfloor+1$, $k\geq38$')
    bx.axhline(4, color=GREY, lw=0.8, ls=':')
    bx.plot(kx, [exact[k] for k in kx], 'o', color=BLACK, ms=2.5)
    bx.plot([14], [7], 's', mfc='white', mec=VERMIL, mew=1, ms=3.5)
    bx.set(xlabel='depth $k$', xlim=(0, 81), ylim=(0, 30), xticks=range(0, 81, 20))
    bx.grid(alpha=0.15, lw=0.5)
    bx.set_title('(b) proved bounds at large depth')
    handles, labels = ax.get_legend_handles_labels()
    hb, lb = bx.get_legend_handles_labels()
    handles.append(hb[2]); labels.append(lb[2])
    fig.legend(handles, labels, loc='outside lower center', ncol=3, fontsize=6.3, handlelength=1.8,
               columnspacing=1.2)
    save(fig, out / 'lfd-depth.pdf')
    plt.close(fig)


def fig_closure_dims(plt, out):
    d4, d5, d6, rec, codim = closure_dims()
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.1, 2.6), layout='constrained',
                                 gridspec_kw={'width_ratios': [1.45, 1]})
    for D, data, col, mk in ((4, d4, SKY, 'o'), (5, d5, BLUE, 's'), (6, d6, VERMIL, 'D')):
        kr = range(2, 17)
        ax.plot(list(kr), [m_le(k, D) for k in kr], color=col, lw=0.9, alpha=0.9,
                label=rf'$|M_{{\leq{D}}}|=\binom{{2k}}{{\leq{D}}}$')
        ks = sorted(data)
        ax.plot(ks, [data[k] for k in ks], mk, color=col, mfc=col, ms=3.6, ls='none',
                label=rf'$\dim V_{D}(k)$')
    ax.plot(sorted(rec), [rec[k] for k in sorted(rec)], '+', color=BLACK, ms=7, mew=0.9,
            label=r'recursion from $V_5(k-1)$')
    ax.set_yscale('log')
    ax.set(xlabel='depth $k$', ylabel='dimension', xlim=(2.5, 16.5), xticks=range(3, 17), ylim=(30, 4e7))
    ax.grid(alpha=0.15, lw=0.5, which='major')
    ax.legend(loc='upper left', fontsize=5.9, ncol=2, handlelength=1.4, columnspacing=0.8)
    ax.set_title(r'(a) closures against the degree-$D$ space')
    Ds = [3, 4, 5, 6]
    bx.semilogy(Ds, [codim[D] for D in Ds], 'o-', color=BLACK, ms=4, lw=0.9)
    for D in Ds:
        bx.annotate(f'{codim[D]:,}\n$k={4 * D - 10}$', (D, codim[D]), textcoords='offset points',
                    xytext=(-4, 6), ha='right', fontsize=6.5)
    bx.set(xlabel='degree $D$', ylabel=r'codim of $V_D$ in $M_{\leq D}$ at $k=4D-10$',
           xticks=Ds, xlim=(2.4, 6.4), ylim=(3, 3e5))
    bx.grid(alpha=0.15, lw=0.5)
    bx.set_title('(b) first deficient depth')
    save(fig, out / 'closure-dims.pdf')
    plt.close(fig)


def fig_msummands(plt, out):
    from matplotlib.lines import Line2D
    cells, count = msummand_points()
    fig, ax = plt.subplots(figsize=(4.2, 2.6), layout='constrained')
    xs = range(5, 18)
    ax.plot(list(xs), [x + 1 for x in xs], color=GREY, lw=0.8, ls='--')
    ax.text(13.4, 14.7, r'$\deg g+1$', color=GREY, fontsize=6.5, ha='right', va='bottom')
    colors = {3: GREEN, 4: BLUE, 5: VERMIL}
    marks = {3: 'o', 4: 's', 5: '^'}
    for m in (3, 4, 5):
        sub = [c for c in cells if c[0] == m]
        exact = sorted({(c[2], c[3]) for c in sub})
        for x, y in exact:
            ks = ','.join(str(c[1]) for c in sub if (c[2], c[3]) == (x, y))
            ax.plot([x], [y], marks[m], color=colors[m], ms=5.5, ls='none', zorder=4)
            ax.annotate(f'$k={ks}$', (x, y), textcoords='offset points', xytext=(5, -9), fontsize=6,
                        color=colors[m])
        for c in sub:
            bound = m * c[1] + 1
            if bound != c[3]:
                ax.plot([c[2]], [bound], marks[m], mfc='white', mec=colors[m], ms=5.5, mew=1, zorder=4)
                ax.annotate(f'$k={c[1]}$', (c[2], bound), textcoords='offset points', xytext=(5, -2.5),
                            fontsize=6, color=colors[m])
                ax.plot([c[2], c[2]], [c[3] + 0.3, bound - 0.3], color=colors[m], lw=0.6, ls=':')
    handles = [Line2D([], [], marker=marks[m], color=colors[m], ls='none', ms=5, label=f'$m={m}$') for m in (3, 4, 5)]
    handles += [Line2D([], [], marker='o', color=BLACK, ls='none', ms=5, label='exact LFD'),
                Line2D([], [], marker='o', mfc='white', mec=BLACK, ls='none', ms=5,
                       label=r'kernel bound $mk+1$ (where larger)')]
    ax.legend(handles=handles, loc='upper left', fontsize=6, handlelength=1.0)
    ax.set(xlabel=r'chain degree of the generator $g=\mathrm{NF}(S_{m+1})$', ylabel='last fall degree',
           xlim=(5, 16.5), ylim=(5.5, 17.8), xticks=range(6, 17, 2))
    ax.grid(alpha=0.15, lw=0.5)
    assert count[(4, 4, 12, 13)] == 6
    save(fig, out / 'msummands.pdf')
    plt.close(fig)


def fig_binary(plt, out):
    import random
    rows = binary_rows()
    rng = random.Random(20261112)
    fig, ax = plt.subplots(figsize=(4.4, 2.5), layout='constrained')
    spec = {'n=2k-1': (SKY, 'o', r'random, $n=2k-1$ ($|Z|\leq4$)'),
            'other': (GREY, 'o', 'random or invariant, other $(n,k)$'),
            'trace-zero': (VERMIL, 's', 'trace-zero hyperplane'),
            'subfield': (GREEN, 'D', 'scaled subfield')}
    for c, (col, mk, lab) in spec.items():
        sub = [r for r in rows if r[2] == c]
        ax.plot([r[0] + rng.uniform(-0.18, 0.18) for r in sub], [r[1] + rng.uniform(-0.12, 0.12) for r in sub],
                mk, color=col, ms=4, alpha=0.85, ls='none', label=f'{lab} ({len(sub)})')
    ax.set(xlabel='dimension $k$ of the coordinate space', ylabel='exact LFD', xticks=range(10, 17),
           yticks=range(3, 9), xlim=(9.4, 16.6), ylim=(3.3, 8.2))
    ax.grid(alpha=0.15, lw=0.5)
    ax.legend(loc='upper right', fontsize=6, handlelength=1.2)
    save(fig, out / 'binary-deep.pdf')
    plt.close(fig)


SWEEP_KINDS = [('prime', 'Prime'), ('binary', 'Binary ordinary'), ('koblitz', 'Koblitz')]


def sweep_groups():
    s = load(SOURCES['sweep'])
    assert s['capped_systems'] == 0 and all(g['capped_systems'] == 0 for g in s['groups'])  # nothing capped to show
    assert sum(sum(g['exact_lfd_histogram'].values()) for g in s['groups']) == s['unique_systems'] == 3005
    return s['groups']


def fig_sweep(plt, out):
    groups_all = sweep_groups()
    colors = {'0': '#c9ced6', '1': '#aab3c0', '2': '#8d99ae', '3': '#90be6d', '4': '#277da8', '5': '#f8961e',
              '6': '#d1495b'}
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6), layout='constrained')
    for ax, (kind, label) in zip(axes, SWEEP_KINDS):
        groups = [g for g in groups_all if g['kind'] == kind]
        ks = sorted({g['k'] for g in groups})
        bases = sorted({g['base_kind'] for g in groups})
        width = 0.8 / max(1, len(bases))
        top = 0
        for bi, base in enumerate(bases):
            xs = [k + (bi - (len(bases) - 1) / 2) * width for k in ks]
            bottom = [0] * len(ks)
            for degree in sorted(colors):
                vals = [next((g['exact_lfd_histogram'].get(degree, 0) for g in groups
                              if g['base_kind'] == base and g['k'] == k), 0) for k in ks]
                ax.bar(xs, vals, width=width, bottom=bottom, color=colors[degree], edgecolor='#444444',
                       linewidth=0.3, hatch='//' if base.startswith('random') else None,
                       label=f'LFD {degree}' if bi == 0 and any(vals) else None)
                bottom = [a + b for a, b in zip(bottom, vals)]
            top = max(top, max(bottom))
        ax.set(xlabel='$k$', ylabel='distinct systems', title=f'{label} (hatched = random base)', xticks=ks)
        ax.set_ylim(0, 1.6 * top)            # headroom so that no bar is clipped and the legend sits above the bars
        ax.grid(alpha=.2)
        rx = ax.twinx()                      # the bound is a degree, so it gets its own labelled axis
        rx.spines['right'].set_visible(True)
        rx.plot(ks, [2 * k + 1 for k in ks], color='#176b87', linewidth=1, linestyle='--', marker='.',
                label='kernel bound $2k+1$')
        rx.set_ylabel('degree (bound $2k+1$)', color='#176b87')
        rx.set_ylim(0, 2 * max(ks) + 3)
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = rx.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=6, frameon=False, loc='upper left', ncol=2)
    save(fig, out / 'curve-sweep.pdf')
    fig.savefig(out / 'curve-sweep.png', dpi=160, metadata={'Software': None})
    plt.close(fig)


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    style(plt)
    out = ROOT / 'paper/figures'
    out.mkdir(exist_ok=True)
    fig_sweep(plt, out)
    fig_lfd_depth(plt, out)
    fig_closure_dims(plt, out)
    fig_msummands(plt, out)
    fig_binary(plt, out)
    return ['paper/figures/curve-sweep.pdf', 'paper/figures/lfd-depth.pdf', 'paper/figures/closure-dims.pdf',
            'paper/figures/msummands.pdf', 'paper/figures/binary-deep.pdf']


if __name__ == '__main__':
    print('\n'.join(main()))
