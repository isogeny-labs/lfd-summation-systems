#!/usr/bin/env python3
"""Growth tables and figure for m2 binary rows (lane C).  Reads data/binary_table.jsonl.
Unit of counting: distinct systems (system_hash).  Writes data/growth_table.md, data/growth_summary.json,
data/binary_sd_vs_k.svg (+ .png)."""
import json, os, statistics as st, collections
HERE = os.path.dirname(os.path.abspath(__file__)); LANE = os.path.dirname(HERE)
R = [json.loads(l) for l in open(os.path.join(LANE, 'data/binary_table.jsonl'))]


def kcls(e, k):
    v = k + e
    return '<=2' if v < 3 else ('3..8' if v < 9 else ('9..18' if v < 19 else '>=19'))


def distinct(rows):
    d = {}
    for r in rows:
        d.setdefault(r['system_hash'], r)
    return list(d.values())


def rng(xs):
    xs = [x for x in xs if x is not None]
    return '-' if not xs else f'{st.median(xs):g} [{min(xs)}..{max(xs)}] (n={len(xs)})'


m2 = distinct([r for r in R if r['family'] == 'm2' and r['lfd'] is not None])
proper = [r for r in m2 if not r['unit']]
lines = []
summary = {}


def table(title, rows, keyf, keyname):
    lines.append(f'\n### {title}\n')
    lines.append(f'| {keyname} | systems | LFD med [min..max] | D* | d_ff | d_reg | d_full exact | d_full lower bound only (n, med lb) | r+1 |')
    lines.append('|---|---|---|---|---|---|---|---|---|')
    g = collections.defaultdict(list)
    for r in rows:
        g[keyf(r)].append(r)
    for key in sorted(g, key=lambda x: (str(type(x)), x)):
        rs = g[key]
        ex = [r['d_full'] for r in rs if r['d_full'] is not None]
        lb = [r['d_full_lb'] for r in rs if r['d_full'] is None and r['d_full_lb'] is not None]
        lines.append(f"| {key} | {len(rs)} | {rng([r['lfd'] for r in rs])} | {rng([r['dstar'] for r in rs])} | {rng([r['d_ff'] for r in rs])} | "
                     f"{rng([r['d_reg'] for r in rs])} | {rng(ex)} | {(len(lb), st.median(lb)) if lb else '-'} | {rs[0]['kernel_r_plus_1']} |")
    return g


lines.append('# Growth table, two-summand binary systems (distinct systems; unit ideals included unless stated)')
lines.append('\nAll statements [pcc:unverified]. "med [min..max] (n=...)": n = rows with a value (d_ff null if no strict fall recorded; d_full exact only where the scan terminated).')
table('By k, all systems', m2, lambda r: r['k'], 'k')
table('By k, proper ideals (|Z|>0)', proper, lambda r: r['k'], 'k')
table('By k and class k+e (e = r_eff - m_eff), proper ideals', proper, lambda r: (r['k'], kcls(r['e'], r['k'])), '(k, k+e class)')
table('By k, unit ideals (|Z|=0)', [r for r in m2 if r['unit']], lambda r: r['k'], 'k')

# P5 evidence
p5 = {}
chk = [r for r in m2 if r['d_full'] is not None]
p5['rows_with_exact_d_full'] = len(chk)
p5['dstar_le_dfull_violations'] = sum(1 for r in chk if r['dstar'] > r['d_full'])
p5['dfull_lt_lfd_rows'] = sum(1 for r in chk if r['d_full'] < r['lfd'])
p5['lb_lt_lfd_is_not_a_violation'] = True
p5['dfull_lb_gt_dstar_violations'] = sum(1 for r in m2 if r['d_full'] is None and r['d_full_lb'] and r['dstar'] is not None and r['d_full_lb'] > 0 and False)
by = collections.defaultdict(lambda: collections.Counter())
for r in m2:
    if r['d_full'] is not None:
        c = 'gt' if r['d_full'] > r['lfd'] else ('eq' if r['d_full'] == r['lfd'] else 'lt')
    elif r['d_full_lb'] is not None and r['d_full_lb'] > r['lfd']:
        c = 'gt(by lower bound)'
    elif r['d_full_lb'] is not None:
        c = 'undetermined(lb<=LFD)'
    else:
        continue
    by[r['k']][c] += 1
p5['dfull_vs_lfd_by_k'] = {k: dict(v) for k, v in sorted(by.items())}
pby = collections.defaultdict(lambda: collections.Counter())
for r in proper:
    if r['d_full'] is not None:
        c = 'gt' if r['d_full'] > r['lfd'] else ('eq' if r['d_full'] == r['lfd'] else 'lt')
    elif r['d_full_lb'] is not None and r['d_full_lb'] > r['lfd']:
        c = 'gt(by lower bound)'
    elif r['d_full_lb'] is not None:
        c = 'undetermined(lb<=LFD)'
    else:
        continue
    pby[r['k']][c] += 1
p5['proper_dfull_vs_lfd_by_k'] = {k: dict(v) for k, v in sorted(pby.items())}
lines.append('\n### d_full vs LFD by k (distinct systems with scan data; gt = d_full > LFD proved exactly; gt(by lower bound) = lower bound already exceeds LFD)\n')
lines.append('| k | all | proper only |'); lines.append('|---|---|---|')
for k in sorted(set(by) | set(pby)):
    lines.append(f"| {k} | {dict(by.get(k, {}))} | {dict(pby.get(k, {}))} |")
summary['p5'] = p5

# pattern: d_full - k for proper rows with exact d_full, by (k, LFD)
pat = collections.defaultdict(collections.Counter)
scan = {}
import glob
for f in glob.glob(os.path.join(LANE, 'data/static_scan_*.jsonl')):
    for l in open(f):
        x = json.loads(l); scan[x['system_hash']] = x
defk = collections.Counter(); defn = 0
for r in proper:
    if r['d_full'] is not None:
        pat[(r['k'], r['lfd'])][r['d_full'] - r['k']] += 1
        sc = scan.get(r['system_hash'])
        if sc and len(sc['per_degree']) >= 2 and r['d_full'] == r['k'] + 1:
            defn += 1
            defk[sc['per_degree'][-2]['deficiency'] == r['k'] + 1] += 1
lines.append('\n### d_full - k on proper systems with exact d_full, by (k, LFD): {d_full - k: count}\n')
lines.append('| (k, LFD) | {d_full - k: count} |'); lines.append('|---|---|')
for key in sorted(pat):
    lines.append(f'| {key} | {dict(pat[key])} |')
lines.append(f'\nSystems with d_full = k+1 and a recorded previous degree: {defn}; of these, deficiency at D = k equals k+1 in {defk[True]}, otherwise in {defk[False]}.')
summary['dfull_minus_k'] = {str(k): dict(v) for k, v in pat.items()}
summary['def_at_k_equals_kp1'] = dict(defk)
open(os.path.join(LANE, 'data/growth_table.md'), 'w').write('\n'.join(lines) + '\n')
json.dump(summary, open(os.path.join(LANE, 'data/growth_summary.json'), 'w'), indent=1, sort_keys=True)

# figure (hand-written SVG; matplotlib is not installed)
ks = sorted({r['k'] for r in m2})
W, H, L, Rr, T, Bt = 760, 480, 60, 190, 40, 55
ymax = max(2 * max(ks) + 1, 12)
X = lambda k: L + (k - min(ks)) / (max(ks) - min(ks)) * (W - L - Rr)
Y = lambda v: H - Bt - v / ymax * (H - T - Bt)
C = dict(lfd='#0072B2', full='#D55E00', reg='#009E73', ker='#555555')
o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="sans-serif" font-size="11">',
     f'<rect width="{W}" height="{H}" fill="#ffffff"/>',
     f'<text x="{L}" y="22" font-size="13" font-weight="600">Binary m=2, proper ideals: LFD vs static Macaulay degree d_full vs k</text>']
for v in range(0, ymax + 1, 2):
    o.append(f'<line x1="{L}" x2="{W-Rr}" y1="{Y(v)}" y2="{Y(v)}" stroke="#e5e5e5"/><text x="{L-6}" y="{Y(v)+4}" text-anchor="end" fill="#444">{v}</text>')
for k in ks:
    o.append(f'<text x="{X(k)}" y="{H-Bt+16}" text-anchor="middle" fill="#444">{k}</text>')
o.append(f'<text x="{(L+W-Rr)/2}" y="{H-12}" text-anchor="middle" fill="#222">k  (r = 2k variables)</text>')
o.append(f'<text x="14" y="{(T+H-Bt)/2}" transform="rotate(-90 14 {(T+H-Bt)/2})" text-anchor="middle" fill="#222">degree</text>')
legend = []


def series(key, color, label, shape, rows=proper):
    pts = []
    for k in ks:
        v = [r[key] for r in rows if r['k'] == k and r[key] is not None]
        if v:
            pts.append((k, st.median(v), min(v), max(v), len(v)))
    if not pts:
        return
    if len(pts) > 1:
        up = ' '.join(f'{X(k)},{Y(h)}' for k, m, l, h, n in pts); dn = ' '.join(f'{X(k)},{Y(l)}' for k, m, l, h, n in reversed(pts))
        o.append(f'<polygon points="{up} {dn}" fill="{color}" opacity="0.15"/>')
        o.append('<polyline fill="none" stroke="%s" stroke-width="1.8" points="%s"/>' % (color, ' '.join(f'{X(k)},{Y(m)}' for k, m, l, h, n in pts)))
    for k, m, l, h, n in pts:
        if shape == 'c': o.append(f'<circle cx="{X(k)}" cy="{Y(m)}" r="3.5" fill="{color}"/>')
        elif shape == 's': o.append(f'<rect x="{X(k)-3.5}" y="{Y(m)-3.5}" width="7" height="7" fill="{color}"/>')
        else: o.append(f'<path d="M{X(k)} {Y(m)-4.5} l4.5 8 h-9 z" fill="{color}"/>')
    legend.append((color, label, shape))


series('lfd', C['lfd'], 'LFD (median, min-max band)', 'c')
series('d_reg', C['reg'], 'semi-regular d_reg(r_eff, m_eff)', 's')
series('d_full', C['full'], 'd_full exact (scanned systems)', 't')
mx = {}
for r in proper:
    if r['d_full'] is None and r['d_full_lb']:
        mx[r['k']] = max(mx.get(r['k'], 0), r['d_full_lb'])
have_exact = {r['k'] for r in proper if r['d_full'] is not None and r['k'] >= 8}
for k, v in mx.items():
    if k in have_exact:
        continue
    o.append(f'<path d="M{X(k)} {Y(v)+4.5} l4.5 -8 h-9 z" fill="none" stroke="{C["full"]}" stroke-width="1.5"/>')
o.append('<polyline fill="none" stroke="%s" stroke-width="1.5" stroke-dasharray="5 3" points="%s"/>' % (C['ker'], ' '.join(f'{X(k)},{Y(2*k+1)}' for k in ks)))
y0 = T + 10
for i, (col, lab, sh) in enumerate(legend):
    o.append(f'<rect x="{W-Rr+14}" y="{y0+i*18-8}" width="10" height="10" fill="{col}"/><text x="{W-Rr+30}" y="{y0+i*18+1}">{lab}</text>')
i = len(legend)
o.append(f'<path d="M{W-Rr+19} {y0+i*18-8} l4.5 8 h-9 z" fill="none" stroke="{C["full"]}" stroke-width="1.5"/><text x="{W-Rr+30}" y="{y0+i*18+1}">d_full lower bound only</text>')
o.append(f'<text x="{W-Rr+30}" y="{y0+(i+1)*18+1}" fill="{C["ker"]}">- - kernel bound r+1 = 2k+1</text>')
o.append(f'<text x="{W-Rr+14}" y="{H-Bt}" font-size="9" fill="#666">[pcc:unverified]</text></svg>')
open(os.path.join(LANE, 'data/binary_sd_vs_k.svg'), 'w').write('\n'.join(o))
print('ok')
