#!/usr/bin/env python3
"""Build the planned cells: instances (equations, Z by field enumeration, structure), cells.json (run order, degree
caps, estimates) and the frozen predictions file.  NO closure is run here (only equation construction, the linear
algebra of the structure e = r_eff - m_eff, and the field enumeration of Z).

  python3 make_cells.py          -> instances/*.json, cells.json, predictions-<date>.json, plan-numbers.json
The predictions file is write-once (the script refuses to overwrite it).
"""
import json
import hashlib
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import ic_bounds_bin_big_run as R  # noqa: E402
from csd import ic_bounds_bin_big as EB  # noqa: E402

DATE = '2026-11-12'
SEED = 20261112
RATE = 25e12          # effective float32 flops/s assumed on the H200 (conservative; see PLAN.md)
RATIO = 0.28          # measured total flops / N^3 (update 0.19 + residual 0.09), see PLAN.md calibration

# cell: name, purpose, systems [(n, k, base, curve, target, count, extra)], max_degree per k rule
CELLS = [
    dict(cell='C1', purpose='discriminating cell (n,k)=(14,11): e=8, k+e=19: law B 5, post-hoc 4', n=14, k=11, maxD=6,
         groups=[('random', 'ordinary', 'random-rational', 3), ('random', 'koblitz0', 'random-rational', 1),
                 ('random', 'koblitz1', 'random-rational', 1), ('invariant', 'koblitz0', 'random-rational', 2),
                 ('random', 'ordinary', 'order-two', 1)]),
    dict(cell='C2a', purpose='trace-zero hyperplane n=k+1, k=11 (e=11)', n=12, k=11, maxD=6,
         groups=[('trace-zero', 'koblitz0', 'random-rational', 2), ('trace-zero', 'koblitz1', 'order-two', 1)]),
    dict(cell='C3a', purpose='small |Z| n=2k-1, k=11 (e=1): does the floor rise to 5?', n=21, k=11, maxD=6,
         groups=[('random', 'ordinary', 'random-rational', 2), ('random', 'koblitz0', 'random-rational', 1)]),
    dict(cell='C3b', purpose='small |Z| n=2k-1, k=12', n=23, k=12, maxD=6, groups=[('random', 'ordinary', 'random-rational', 2)]),
    dict(cell='C2b', purpose='trace-zero hyperplane n=k+1, k=12 (e=12)', n=13, k=12, maxD=6,
         groups=[('trace-zero', 'koblitz0', 'random-rational', 2)]),
    dict(cell='C5a', purpose='subfield F_(2^10) in F_(2^20), order-two, written value k/2+1 = 6', n=20, k=10, maxD=7,
         groups=[('subfield', 'koblitz0', 'order-two', 1), ('subfield', 'koblitz1', 'order-two', 1)]),
    dict(cell='C6a', purpose='discriminating (16,12): e=8, k+e=20: law B 5, post-hoc 4', n=16, k=12, maxD=6,
         groups=[('random', 'ordinary', 'random-rational', 3), ('random', 'koblitz0', 'random-rational', 1)]),
    dict(cell='C3c', purpose='small |Z| n=2k-1, k=13', n=25, k=13, maxD=6, groups=[('random', 'ordinary', 'random-rational', 2)]),
    dict(cell='C4a', purpose='Koblitz (15,11): random and Frobenius-invariant bases (e=7)', n=15, k=11, maxD=6,
         groups=[('random', 'koblitz0', 'random-rational', 2), ('invariant', 'koblitz0', 'random-rational', 2)]),
    dict(cell='C2c', purpose='trace-zero hyperplane n=k+1, k=13 (e=13)', n=14, k=13, maxD=6,
         groups=[('trace-zero', 'koblitz1', 'random-rational', 2)]),
    dict(cell='C4b', purpose='Koblitz (21,12): random and Frobenius-invariant bases (e=3..4)', n=21, k=12, maxD=6,
         groups=[('random', 'koblitz1', 'random-rational', 2), ('invariant', 'koblitz1', 'random-rational', 2)]),
    dict(cell='C6b', purpose='discriminating (18,13): e=8, k+e=21: law B 5, post-hoc 4', n=18, k=13, maxD=6,
         groups=[('random', 'ordinary', 'random-rational', 3)]),
    dict(cell='C3d', purpose='small |Z| n=2k-1, k=14', n=27, k=14, maxD=6, groups=[('random', 'ordinary', 'random-rational', 2)]),
    dict(cell='C5b', purpose='subfield F_(2^12) in F_(2^24), order-two, written value 7 (needs D = 7: 536,155 columns)', n=24, k=12, maxD=7,
         groups=[('subfield', 'koblitz0', 'order-two', 1)]),
    dict(cell='C3e', purpose='small |Z| n=2k-1, k=15 (D=6 is 768,212 columns: cap 5 unless budget remains)', n=29, k=15, maxD=5,
         groups=[('random', 'ordinary', 'random-rational', 1)]),
    dict(cell='C3f', purpose='small |Z| n=2k-1, k=16 (D=6 is 1,149,017 columns: cap 5)', n=31, k=16, maxD=5,
         groups=[('random', 'ordinary', 'random-rational', 1)]),
    # reserve cells (only if the budget allows; ordered after everything else)
    dict(cell='R1', purpose='reserve: more systems at the discriminating cell (14,11) (within-cell variation)', n=14, k=11, maxD=6, reserve=True,
         groups=[('random', 'ordinary', 'random-rational', 6), ('random', 'koblitz0', 'random-rational', 3)]),
    dict(cell='R2', purpose='reserve: (17,12): e=7, k+e=19: law B 5, post-hoc 4', n=17, k=12, maxD=6, reserve=True,
         groups=[('random', 'ordinary', 'random-rational', 4)]),
]


def A_law(e):
    return 2 if e <= -1 else 3 if e <= 3 else 4 if e <= 8 else 5


def B_law(x):
    return 2 + (x >= 3) + (x >= 9) + (x >= 19)


def posthoc(k, e):
    return max(2 + ((k + e) >= 3) + ((k + e) >= 9), A_law(e))


def build(spec, want_class):
    for t in range(60):
        s = dict(spec, seed=SEED + t)
        inst = R.build_instance(s)
        if inst['class'] == want_class or want_class is None:
            return inst, t
    raise RuntimeError('no instance of class %s for %s' % (want_class, spec['name']))


def main():
    pred_path = os.path.join(HERE, 'predictions-%s.json' % DATE)
    assert not os.path.exists(pred_path), 'predictions file exists: preregistration is write-once'
    os.makedirs(os.path.join(HERE, 'instances'), exist_ok=True)
    rows, order, numbers = [], [], []
    for c in CELLS:
        systems = []
        idx = 0
        for base, curve, target, count in c['groups']:
            for j in range(count):
                name = '%s_%s_%s_%s_%d' % (c['cell'], base[:4], curve, target[:3], j)
                spec = dict(name=name, n=c['n'], k=c['k'], base=base, curve=curve, target=target, seed=SEED, index=j if base == 'invariant' else 0)
                want = 'order-two' if target == 'order-two' else 'affine'
                t0 = time.time()
                inst, attempt = build(spec, want)
                path = os.path.join(HERE, 'instances', name + '.json')
                with open(path, 'w') as f:
                    json.dump(inst, f, separators=(',', ':'))
                st = inst['structure']
                k, e, Z = inst['k'], st['e'], inst['n_points']
                row = dict(cell=c['cell'], name=name, instance=name + '.json', n=c['n'], k=k, base=base, curve=curve, curve_family=inst['curve']['family'],
                           target=target, cls=inst['class'], frobenius_invariant=inst['frobenius_invariant'], trace_zero=inst['trace_zero'],
                           subfield=inst['scaled_subfield'], r_eff=st['r_eff'], m_eff=st['m_eff'], e=e, k_plus_e=k + e, Z=Z,
                           system_hash=inst['system_hash'], points_sha256=inst['points_sha256'], draw_attempt=attempt)
                if inst['class'] == 'affine':
                    row.update(pred_B=int(B_law(k + e)), pred_posthoc=int(posthoc(k, e)), pred_A=int(A_law(e)),
                               pred_A_Z=(5 if Z > 470 else 4))
                    pm = max(row['pred_B'], row['pred_posthoc'], row['pred_A'], row['pred_A_Z'])
                elif inst['scaled_subfield']:
                    row.update(pred_inverse_family=max(3, k // 2 + 1), pred_B_unscored=int(B_law(k + e)))
                    pm = row['pred_inverse_family']
                else:
                    row.update(pred_B_unscored=int(B_law(k + e)), pred_posthoc_unscored=int(posthoc(k, e)))
                    pm = max(row['pred_B_unscored'], row['pred_posthoc_unscored'])
                row['largest_predicted'] = int(pm)
                # certificate degree is LFD or (rarely) LFD + 1: the ladder cap
                Dcap = min(c['maxD'], pm + 1)
                row['max_degree'] = Dcap
                pl = EB.plan(2 * k, Dcap)
                est = RATIO * pl['columns'] ** 3 / RATE
                # time to the predicted degree plus the one above, geometric: dominated by the last degree
                row.update(columns_at_cap=pl['columns'], est_seconds_cap=round(est, 1), predicted_peak_gb=pl['predicted_peak_gb'],
                           basis_worst_gb=pl['basis_worst_gb'])
                rows.append(row)
                systems.append(dict(name=name, instance=name + '.json', max_degree=Dcap, est_seconds=round(est * 1.15 + 60),
                                    time_limit_s=int(max(1800, est * 6 + 600))))
                print('%s k=%d n=%d %s |Z|=%d e=%d  preds=%s  D<=%d N=%d est=%.0fs  (%.0fs build)' % (
                    name, k, c['n'], inst['class'], Z, e, {x: row[x] for x in row if x.startswith('pred')}, Dcap, pl['columns'], est, time.time() - t0), flush=True)
        order.append(dict(cell=c['cell'], purpose=c['purpose'], reserve=bool(c.get('reserve')), systems=systems, est_seconds=sum(s['est_seconds'] for s in systems)))
    # also record the per-(k, D) plan numbers
    for k in range(11, 17):
        for D in (5, 6, 7):
            numbers.append(EB.plan(2 * k, D))
    json.dump(numbers, open(os.path.join(HERE, 'plan-numbers.json'), 'w'), indent=1)
    json.dump(dict(order=order, rate_flops_per_s=RATE, flops_over_N3=RATIO), open(os.path.join(HERE, 'cells.json'), 'w'), indent=1)
    out = dict(written_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), seed=SEED,
               rule='written before any closure at k >= 11; instances built by ic_bounds_bin_big_run.build_instance; Z by field enumeration only',
               laws=dict(A='affine: 2/3/4/5 for e<=-1 / <=3 / <=8 / >=9', B='affine: 2+[k+e>=3]+[k+e>=9]+[k+e>=19]',
                         posthoc='affine: max(2+[k+e>=3]+[k+e>=9], A(e))  (post-hoc conjecture of the law README)',
                         A_Z='affine: 5 if |Z|>470 else 4 (k=10 grid rule)', inverse_family='scaled subfield F_(2^k) in F_(2^2k), order-two target: max(3, k/2+1)',
                         const='1 in the input span: LFD 2', order_two='recorded, unscored (B shown for reference)'),
               scoring='per distinct affine system: hit iff the certified exact LFD equals the prediction; unscored classes are recorded only',
               cases=rows)
    with open(pred_path, 'w') as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print('wrote', pred_path, len(rows), 'cases; sha256', hashlib.sha256(open(pred_path, 'rb').read()).hexdigest())


if __name__ == '__main__':
    main()
