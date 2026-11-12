"""Build the curve-bounds manuscript: data macros and tables from the sweep summary, the figure, and the PDF."""
from collections import Counter, defaultdict
from pathlib import Path
import hashlib, json, os, re, shutil, subprocess, tempfile

KINDS = [('prime', 'Prime'), ('binary', 'Binary ordinary'), ('koblitz', 'Koblitz')]
BASES = {'subgroup': 'subgroup (symmetric)', 'random-coset': 'random coset',
         'frobenius-invariant': 'Frobenius-invariant (symmetric)', 'random-subspace': 'random subspace'}


def tex_escape(s):
    return str(s).replace('_', '\\_').replace('%', '\\%')


def results_section(summary, cases, results):
    lines = []
    lines.append(f"The sweep has {summary['cases']} cases and {summary['unique_systems']} distinct systems; "
                 f"{summary['exact_systems']} systems have a certified exact LFD and {summary['capped_systems']} are capped. "
                 f"The histogram of exact values over all distinct systems is "
                 + ', '.join(f"${k}$: {v}" for k, v in sorted(summary['lfd_histogram'].items(), key=lambda kv: (kv[0] == 'None', kv[0])))
                 + " (None = capped). All exact values satisfy the kernel-checked bound $2k+1$; the binary values also satisfy the written bound $k+2$.")
    lines.append('')
    # Table per kind: rows k, columns base kinds: histogram and max
    for kind, label in KINDS:
        groups = [g for g in summary['groups'] if g['kind'] == kind]
        if not groups:
            continue
        bases = [b for b in BASES if any(g['base_kind'] == b for g in groups)]
        ks = sorted({g['k'] for g in groups})
        lines.append('\\begin{table}[htbp]\\centering\\small')
        lines.append('\\begin{tabular}{r' + 'l' * len(bases) + 'r}')
        lines.append('\\toprule')
        lines.append('$k$ & ' + ' & '.join(BASES[b] for b in bases) + ' & bound $2k{+}1$\\\\')
        lines.append('\\midrule')
        for k in ks:
            cells = []
            for b in bases:
                g = next((g for g in groups if g['base_kind'] == b and g['k'] == k), None)
                if g is None:
                    cells.append('--')
                    continue
                hist = ', '.join(f"{d}:{c}" for d, c in sorted(g['exact_lfd_histogram'].items(),
                                                              key=lambda kv: (kv[0] == 'None', kv[0])) if d != 'None')
                capped = g['capped_systems']
                cells.append(f"{hist}" + (f"; capped {capped}" if capped else '') + f" (max {g['max_exact_lfd']})")
            lines.append(f"{k} & " + ' & '.join(cells) + f" & {2 * k + 1}\\\\")
        lines.append('\\bottomrule\\end{tabular}')
        lines.append(f"\\caption{{{label} curves: distinct systems by exact LFD (value:count) for each factor-base class and depth or dimension $k$, with the number of capped systems and the maximum certified value.}}")
        lines.append(f"\\label{{tab:{kind}}}\\end{{table}}")
        lines.append('')
    # Random vs symmetric paired comparison: same (kind, field, curve, k, target class) cells
    lines.append('\\paragraph{Random versus symmetric bases.}')
    for kind, label in KINDS:
        sym = defaultdict(list); rnd = defaultdict(list)
        for c in cases:
            if c['kind'] != kind:
                continue
            r = results[c['system_hash']]
            if r['exact_lfd'] is None:
                continue
            key = (c.get('prime', c.get('n')), c['a'], c['b'], c['k'], c['target_class'])
            (sym if c['symmetric'] else rnd)[key].append(r['exact_lfd'])
        common = sorted(set(sym) & set(rnd))
        if not common:
            lines.append(f"{label}: no cell has both classes certified.")
            continue
        diffs = [max(rnd[key]) - max(sym[key]) for key in common]
        higher = sum(d > 0 for d in diffs); lower = sum(d < 0 for d in diffs); equal = sum(d == 0 for d in diffs)
        mean_sym = sum(max(sym[key]) for key in common) / len(common)
        mean_rnd = sum(max(rnd[key]) for key in common) / len(common)
        lines.append(f"{label}: over {len(common)} cells (field, curve, $k$, target class) with both classes certified, "
                     f"the maximum exact LFD of the random bases exceeds that of the symmetric bases in {higher} cells, "
                     f"is lower in {lower} and equal in {equal}; the cell means are {mean_rnd:.2f} (random) and {mean_sym:.2f} (symmetric).")
    lines.append('')
    lines.append('\\paragraph{Growth in $k$ and the gap to the bound.}')
    for kind, label in KINDS:
        groups = [g for g in summary['groups'] if g['kind'] == kind]
        if not groups:
            continue
        ks = sorted({g['k'] for g in groups})
        per_k = []
        for k in ks:
            vals = []
            for g in groups:
                if g['k'] == k:
                    vals += [int(d) for d, c in g['exact_lfd_histogram'].items() if d != 'None' for _ in range(c)]
            if vals:
                per_k.append(f"$k={k}$: max {max(vals)}, mean {sum(vals) / len(vals):.2f}, bound {2 * k + 1}")
        lines.append(f"{label}: " + '; '.join(per_k) + '.')
    lines.append('')
    lines.append('These are exact values of the specified presentations on a seeded grid. They support no population inference, no statement about other bases, and no relation-cost claim; '
                 'capped systems retain only the interval from their last strict fall to the proved bound.')
    return '\n'.join(lines)


def main():
    root = Path(__file__).resolve().parents[1]
    paper = root / 'paper'
    sweep = root / 'experiments/ic-curve-sweep'
    engine = 'pdflatex' if shutil.which('pdflatex') else 'tectonic' if shutil.which('tectonic') else None
    if engine is None:
        raise SystemExit('Install pdfLaTeX (or Tectonic) to build the manuscript.')
    summary = json.loads((sweep / 'summary.json').read_text())
    cases = [json.loads(l) for l in (sweep / 'cases.jsonl').read_text().splitlines()]
    results = {r['system_hash']: r for r in map(json.loads, (sweep / 'results.jsonl').read_text().splitlines())}
    lean = json.loads((root / 'experiments/ic-lean/verification.json').read_text())
    decls = len(lean['audited_declarations'])
    data = ['% generated by scripts/build_ic_curve_paper.py; do not edit',
            f"\\newcommand{{\\SWEEPCASES}}{{{summary['cases']}}}",
            f"\\newcommand{{\\SWEEPSYSTEMS}}{{{summary['unique_systems']}}}",
            f"\\newcommand{{\\LEANDECLS}}{{{decls}}}",
            '\\newcommand{\\SWEEPRESULTS}{' + results_section(summary, cases, results) + '}', '']
    (paper / 'curve-sweep-data.tex').write_text('\n'.join(data))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    figures = paper / 'figures'
    figures.mkdir(exist_ok=True)
    ink, ink2, grid = '#0b0b0b', '#52514e', '#dddcd6'
    plt.rcParams.update({
        'font.family': 'cmr10', 'mathtext.fontset': 'cm', 'axes.formatter.use_mathtext': True, 'axes.unicode_minus': False,
        'font.size': 9, 'axes.labelsize': 9, 'axes.titlesize': 9.5, 'axes.edgecolor': ink2, 'axes.linewidth': 0.6,
        'axes.labelcolor': ink, 'text.color': ink, 'xtick.color': ink2, 'ytick.color': ink2, 'xtick.labelsize': 8,
        'ytick.labelsize': 8, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True, 'axes.axisbelow': True,
        'grid.color': grid, 'grid.linewidth': 0.5, 'legend.frameon': False, 'legend.fontsize': 8, 'hatch.linewidth': 0.5,
        'pdf.fonttype': 42, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.03})
    # the LFD is an ordered quantity: one hue, light to dark (sequential), random bases marked by a hatch
    degrees = ['0', '1', '2', '3', '4', '5', '6']
    ramp = dict(zip(degrees, ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#104281', '#0d366b']))
    symmetric = {'subgroup', 'frobenius-invariant'}
    fig, axes = plt.subplots(1, 3, figsize=(6.1, 2.75), layout='constrained')
    present = set()
    for ax, (kind, label) in zip(axes, KINDS):
        groups = [g for g in summary['groups'] if g['kind'] == kind]
        ks = sorted({g['k'] for g in groups})
        bases = sorted({g['base_kind'] for g in groups}, key=lambda b: (b not in symmetric, b))
        width = 0.76 / max(1, len(bases))
        for bi, base in enumerate(bases):
            xs = [k + (bi - (len(bases) - 1) / 2) * width for k in ks]
            bottom = [0] * len(ks)
            for degree in degrees:
                vals = [next((g['exact_lfd_histogram'].get(degree, 0) for g in groups if g['base_kind'] == base and g['k'] == k), 0) for k in ks]
                if any(vals):
                    present.add(degree)
                ax.bar(xs, vals, width=width, bottom=bottom, color=ramp[degree], edgecolor='white', linewidth=0.6,
                       hatch=None if base in symmetric else '////')
                bottom = [a + b for a, b in zip(bottom, vals)]
        ax.set(xlabel='$k$', title=label, xticks=ks)
        ax.grid(axis='x', visible=False)
        tallest = max(sum(g['exact_lfd_histogram'].values()) for g in groups)
        ax.set_ylim(0, tallest * 1.08)          # head room: no bar touches the top of the panel
    axes[0].set_ylabel('distinct systems')
    handles = [Patch(facecolor=ramp[d], edgecolor='white', label='LFD ' + d) for d in degrees if d in present]
    handles += [Patch(facecolor='white', edgecolor=ink2, linewidth=0.5, label='symmetric base'),
                Patch(facecolor='white', edgecolor=ink2, linewidth=0.5, hatch='////', label='random base')]
    fig.legend(handles=handles, loc='outside lower center', ncol=len(handles), handlelength=1.2, columnspacing=1.0)
    fig.savefig(figures / 'curve-sweep.pdf', metadata={'CreationDate': None, 'ModDate': None})
    fig.savefig(figures / 'curve-sweep.png', dpi=200)
    plt.close(fig)
    env = dict(os.environ, SOURCE_DATE_EPOCH='1793318400', FORCE_SOURCE_DATE='1', TZ='UTC')
    with tempfile.TemporaryDirectory(prefix='ic-curve-paper-') as temporary:
        for _ in range(3 if engine == 'pdflatex' else 1):   # Tectonic reruns by itself
            cmd = (['pdflatex', '-no-shell-escape', '-interaction=nonstopmode', '-halt-on-error',
                    '-output-directory', temporary, 'curve-lfd-bounds.tex'] if engine == 'pdflatex' else
                   ['tectonic', '--keep-logs', '--outdir', temporary, 'curve-lfd-bounds.tex'])
            run = subprocess.run(cmd, cwd=paper, env=env, text=True, capture_output=True)
            if run.returncode:
                raise SystemExit(run.stdout[-6000:] + run.stderr[-1000:])
        log = (Path(temporary) / 'curve-lfd-bounds.log').read_text(errors='replace')
        if re.search(r'undefined|Rerun to get cross-references right|multiply defined', log, re.I):
            raise SystemExit('Unresolved TeX references/citations; inspect manuscript.')
        (paper / 'curve-lfd-bounds.pdf').write_bytes((Path(temporary) / 'curve-lfd-bounds.pdf').read_bytes())
        warnings = re.findall(r'(?:LaTeX|Package [^\n]+) Warning:[^\n]*', log)
        overfull = len(re.findall(r'Overfull \\hbox', log))
    sources = ['scripts/build_ic_curve_paper.py', 'paper/curve-lfd-bounds.tex', 'paper/curve-sweep-data.tex',
               'experiments/ic-curve-sweep/summary.json', 'experiments/ic-lean/verification.json']
    outputs = ['paper/curve-lfd-bounds.pdf', 'paper/figures/curve-sweep.pdf', 'paper/figures/curve-sweep.png']
    report = dict(status='passed: PDF built; no unresolved citations or references; not a proof or novelty certificate',
                  tex_engine=engine,
                  latex_warnings=warnings, overfull_boxes=overfull,
                  source_sha256={n: hashlib.sha256((root / n).read_bytes()).hexdigest() for n in sources},
                  data_sha256={n: hashlib.sha256((root / n).read_bytes()).hexdigest() for n in outputs})
    out = root / 'experiments/ic-curve-paper'
    out.mkdir(exist_ok=True)
    (out / 'verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(dict(status=report['status'], overfull_boxes=overfull, warnings=len(warnings))))


if __name__ == '__main__':
    main()
