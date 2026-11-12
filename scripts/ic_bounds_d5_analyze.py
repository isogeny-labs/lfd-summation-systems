"""Collect the d_5(k) runs, fit the cost model, and project memory/disk/time to k = 12, 13 (and 14).

Reads experiments/ic-bounds-d5/runs/*.json (degree 5, split engine) and d4-check.json, writes
experiments/ic-bounds-d5/summary.json, table.md and d5.svg/png.  Every projection is labelled as a
model output (empirical scaling), never as a measurement.

Cost model (split engine, see csd/ic_bounds_d5.py).  Reducing one row against a part of rank r on
c_part columns costs r * (c_part - r) <= c_part^2 / 4 multiply-adds.  Each of the lambda = dim(V_5 cap
M_<=4) low rows is shifted k times in the closure phase and k times in the verification sweep, and each
shifted row is reduced in both parts, so
    madds(verify)  <= k * lambda * 2 * c_part^2 / 4,
    madds(closure) <= k * lambda * 2 * c_part^2 / 4 (clearing included; measured 0.8 of the bound at k = 11),
with c_part ~ |M_<=5| / 2.  Storage: r * (c_part - r) <= c_part^2 / 4 uint16 entries per part.
lambda(k) is unknown beyond the measured depths: the projection uses the measured lambda(10), lambda(11)
extrapolated linearly (scenario A) and the hard bound lambda <= |M_<=4| (scenario B).
"""
import glob, json, math, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'experiments/ic-bounds-d5'
sys.path.insert(0, str(ROOT))

C = lambda k, D: sum(math.comb(2 * k, d) for d in range(D + 1))


def part_cols(k):
    """(sym, anti) column counts of M_<=5 under the summand swap (orbits; fixed points have a1 = a2)."""
    fixed = sum(math.comb(k, w) for w in range(0, 3))  # a1 = a2: weight 2*wt(a1) <= 5
    total = C(k, 5)
    return (total + fixed) // 2, (total - fixed) // 2


runs = []
for f in sorted(glob.glob(str(OUT / 'runs/*.json'))):
    r = json.loads(Path(f).read_text())
    if r.get('degree') != 5 or r.get('engine') != 'split':
        continue
    runs.append(dict(file=Path(f).name, k=r['k'], prime=r['prime'], a=r['a'], b=r['b'], seed=r['seed'],
                     target=r['target'], points=r['points'], columns=r['columns'], status=r['result']['status'],
                     dim_V5=r['result'].get('rank', r.get('ranks') and sum(r['ranks'])),
                     ideal_dim_lower=r['columns'] - r['points'], low_dimension=r['low_dimension'],
                     profile=r['profile'], sym_profile=r['sym_profile'], anti_profile=r['anti_profile'],
                     ranks=r['ranks'], dim_V4=(r['prev'] or {}).get('dim'), prev_engine=(r['prev'] or {}).get('engine'),
                     work=r.get('work'), total_seconds=r.get('total_seconds'), maxrss_mb=r['maxrss_mb'],
                     peak_free_entries=r.get('peak_free_entries'), basis_sha256=r['basis_sha256'],
                     no_early_stop='full' in Path(f).name))

# measured d_5 per depth (deficient-closed or closed-full runs)
by_k = {}
for r in runs:
    if r['status'] in ('deficient-closed', 'closed-full') or r['status'] == 'reached':
        by_k.setdefault(r['k'], []).append(r)
measured = {}
for k, rs in sorted(by_k.items()):
    dims = sorted({r['dim_V5'] for r in rs if r['status'] == 'deficient-closed'})
    measured[k] = dict(instances=len(rs), statuses=sorted({r['status'] for r in rs}),
                       deficient_dims=dims, low_dims=sorted({r['low_dimension'] for r in rs}),
                       profiles=sorted({tuple(r['profile']) for r in rs}), columns=C(k, 5))

# throughput from the deficient runs at k >= 10 with work counters
thr = {}
for phase in ('closure', 'verify'):
    m = sum(r['work'][phase][0] for r in runs if r['work'] and r['k'] >= 10 and r['status'] == 'deficient-closed')
    s = sum(r['work'][phase][1] for r in runs if r['work'] and r['k'] >= 10 and r['status'] == 'deficient-closed')
    thr[phase] = m / s if s else None
bound_check = []
for r in runs:
    if r['work'] and r['k'] >= 10 and r['status'] == 'deficient-closed':
        cs, ca = part_cols(r['k'])
        bound = r['k'] * r['low_dimension'] * (cs * cs + ca * ca) / 4
        bound_check.append(dict(k=r['k'], file=r['file'], verify_madds=r['work']['verify'][0],
                                closure_madds=r['work']['closure'][0], per_phase_bound=bound,
                                verify_over_bound=r['work']['verify'][0] / bound,
                                closure_over_bound=r['work']['closure'][0] / bound,
                                peak_free_entries=r['peak_free_entries'],
                                storage_bound_entries=(cs * cs + ca * ca) // 4, maxrss_mb=r['maxrss_mb']))
lam = {k: v['low_dims'][0] for k, v in measured.items() if v['low_dims'] and k >= 10
       and any(r['status'] == 'deficient-closed' for r in by_k[k])}
proj = []
ks = sorted(lam)
for k in (12, 13, 14):
    cs, ca = part_cols(k)
    if len(ks) >= 2:
        k0, k1 = ks[-2], ks[-1]
        lamA = lam[k1] + (lam[k1] - lam[k0]) * (k - k1) / (k1 - k0)
    else:
        lamA = None
    lamB = C(k, 4)
    row = dict(k=k, columns=C(k, 5), part_columns=[cs, ca], low_columns=C(k, 4),
               prime=12289 if k <= 12 else (40961 if k == 13 else 65537),
               storage='uint16' if k <= 13 else 'uint32')
    bytes_per = 2 if k <= 13 else 4
    store = (cs * cs + ca * ca) / 4 * bytes_per
    row['storage_bound_GB'] = round(store / 1e9, 2)
    # conservative: 1.5 x storage (one part copied during an update) + row chunks (2 x 2 x 256 x |M| int64)
    # + product tiles (2 x 1024 x c_part float64) + 0.2 GB; overestimates the measured k = 11 peak (0.86 GB) by 1.45x
    row['peak_memory_estimate_GB'] = round((store * 1.5 + 4 * 256 * C(k, 5) * 8 + 2 * 1024 * cs * 8) / 1e9 + 0.2, 2)
    row['checkpoint_bound_GB'] = round(store / 1e9, 2)
    for name, lm in (('A_linear_lambda', lamA), ('B_lambda_le_M4', lamB)):
        if lm is None:
            continue
        per = k * lm * (cs * cs + ca * ca) / 4
        # measured madds/bound: 0.80-0.85 per phase at k = 11, so the bound is used as is (conservative)
        sec = per / thr['verify'] + per / thr['closure'] if thr['verify'] and thr['closure'] else None
        row[name] = dict(lambda_=round(lm), madds_per_phase_bound=per,
                         single_core_hours=round(sec / 3600, 2) if sec else None)
    proj.append(row)

summary = dict(runs=runs, measured=measured, throughput_madds_per_second=thr, bound_check=bound_check,
               projection=proj,
               notes=['throughput measured on one contended core (4-core host shared with three CPU-bound jobs); '
                      'projections are model outputs, not measurements',
                      'storage bound: r*(c_part - r) <= c_part^2/4 entries per part, attained near r = c_part/2'])
d4 = OUT / 'd4-check.json'
if d4.exists():
    summary['d4_check'] = json.loads(d4.read_text())
(OUT / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')

# table
lines = ['| k | p | a | b | seed | |Z| | |M_<=5| | dim V_4 | dim V_5 | dim(V_5 cap M_<=w), w=0..5 | status | s |',
         '|---|---|---|---|---|---|---|---|---|---|---|---|']
for r in sorted(runs, key=lambda r: (r['k'], r['prime'], r['seed'])):
    lines.append(f"| {r['k']} | {r['prime']} | {r['a']} | {r['b']} | {r['seed']} | {r['points']} | {r['columns']} | "
                 f"{r['dim_V4']} | {r['dim_V5']} | {r['profile']} | {r['status']} | {r['total_seconds']} |")
lines += ['', '| k | |M_<=5| | part columns | storage bound (GB) | peak RAM est. (GB) | lambda (A) | hours, 1 core (A) | hours, 1 core (B) |',
          '|---|---|---|---|---|---|---|---|']
for r in proj:
    A, Bs = r.get('A_linear_lambda', {}), r.get('B_lambda_le_M4', {})
    lines.append(f"| {r['k']} | {r['columns']} | {r['part_columns']} | {r['storage_bound_GB']} | {r['peak_memory_estimate_GB']} | "
                 f"{A.get('lambda_')} | {A.get('single_core_hours')} | {Bs.get('single_core_hours')} |")
(OUT / 'table.md').write_text('\n'.join(lines) + '\n')

# figure
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    kk = list(range(6, 14))
    ax[0].semilogy(kk, [C(k, 5) for k in kk], 'k--', label='|M_<=5| (ideal dim. bound)')
    ax[0].semilogy(kk, [C(k, 4) for k in kk], 'k:', label='|M_<=4|')
    pts = [(r['k'], r['dim_V5']) for r in runs if r['status'] == 'deficient-closed']
    full = [(r['k'], r['dim_V5']) for r in runs if r['status'] in ('closed-full', 'reached')]
    lows = [(r['k'], r['low_dimension']) for r in runs if r['status'] == 'deficient-closed']
    if full:
        ax[0].semilogy(*zip(*full), 'o', mfc='none', color='C0', label='dim V_5 = ideal (k<=9)')
    if pts:
        ax[0].semilogy(*zip(*pts), 'o', color='C3', label='dim V_5, verified closed, deficient')
    if lows:
        ax[0].semilogy(*zip(*lows), 's', color='C2', label='lambda = dim(V_5 cap M_<=4)')
    ax[0].set_xlabel('depth k'); ax[0].set_ylabel('dimension'); ax[0].legend(fontsize=7)
    ax[0].set_title('degree-5 closure, exact (measured)')
    pk = [r['k'] for r in proj if 'A_linear_lambda' in r]
    ax[1].semilogy(pk, [r['A_linear_lambda']['single_core_hours'] for r in proj if 'A_linear_lambda' in r], 'o-', label='scenario A (linear lambda)')
    ax[1].semilogy([r['k'] for r in proj], [r['B_lambda_le_M4']['single_core_hours'] for r in proj], 's--', label='scenario B (lambda <= |M_<=4|)')
    mk = [(b['k'], (b['verify_madds'] / thr['verify'] + b['closure_madds'] / thr['closure']) / 3600) for b in bound_check]
    if mk:
        ax[1].semilogy(*zip(*mk), 'k*', ms=10, label='measured (k=10, 11)')
    ax[1].set_xlabel('depth k'); ax[1].set_ylabel('single-core hours (model)'); ax[1].legend(fontsize=7)
    ax[1].set_title('cost model (projection, not measurement)')
    fig.tight_layout()
    fig.savefig(OUT / 'd5.svg'); fig.savefig(OUT / 'd5.png', dpi=110)
except Exception as e:  # the figure is optional
    print('figure skipped:', e)
print(json.dumps(dict(measured={k: dict(v, profiles=[list(p) for p in v['profiles']]) for k, v in measured.items()},
                      throughput=thr, projection=proj), indent=1, default=str))
