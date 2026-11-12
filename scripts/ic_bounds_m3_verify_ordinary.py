#!/usr/bin/env python3
"""Ordinary-ring spot check: inputs built HERE from the verifier's own S_4 / S_3 (chains as explicit
inputs, full monomial space), LFD computed with the independent validator engine
experiments/ic-curve-validation/code/lfd_engine.py.  Split only where I_free = I_field (rational, radical).
Sample: seed 31337, stratified over k in {1,2} x encoding, all special rows (LFD 0/3/6) included.
  python3 scripts/ic_bounds_m3_verify_ordinary.py -> verify/ordinary.jsonl
"""
import json, random, sys, time, itertools, importlib.util
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ic_bounds_m3_verify_engine as V
ROOT = Path(__file__).resolve().parents[1]
EXP = ROOT / 'experiments' / 'ic-bounds-m3'
spec = importlib.util.spec_from_file_location('lfd_engine', ROOT / 'experiments/ic-curve-validation/code/lfd_engine.py')
L = importlib.util.module_from_spec(spec); spec.loader.exec_module(L)


def ordinary(sysd, p, k):
    m, zv = sysd['m'], sysd['zvar']
    r = m * k + (1 if zv else 0)
    inputs = []
    for i in range(m):
        for j in range(k):
            sq = [0] * r; sq[i * k + j] = 2
            f = {tuple(sq): 1}
            if j + 1 < k:
                nx = [0] * r; nx[i * k + j + 1] = 1; f[tuple(nx)] = p - 1
            else:
                f[tuple([0] * r)] = p - 1
            inputs.append(f)
    for g in sysd['inputs']:
        G = {}
        for mm, c in g.items():
            e = [0] * r
            for i in range(m):
                for j in range(k): e[i * k + j] = (mm[i] >> j) & 1
            if zv: e[-1] = mm[m]
            G[tuple(e)] = (G.get(tuple(e), 0) + c) % p
        inputs.append({e: c for e, c in G.items() if c})
    return inputs, r


def points(sysd, p, k, comps):
    m = sysd['m']; pts = []
    for x, h in comps:
        u = [pow(x[i], 1 << j, p) for i in range(m) for j in range(k)]
        if h is None: pts.append(tuple(u))
        else:
            roots = [z for z in range(p) if sum(c * pow(z, e, p) for e, c in enumerate(h)) % p == 0]
            assert len(roots) == len(h) - 1
            pts.extend(tuple(u + [z]) for z in roots)
    return pts


def main():
    cases = [json.loads(l) for l in open(EXP / 'cases.jsonl')]
    rows = {json.loads(l)['case_id']: json.loads(l) for l in open(EXP / 'rows.jsonl')}
    rng = random.Random(31337); terms = V.s4_terms()
    pick = []
    for k in (1, 2):
        ks = [c for c in cases if c['k'] == k]
        special = [c for c in ks if rows[c['case_id']]['direct_lfd'] in (0, 3, 6)]
        pick += [(c, 'direct') for c in special[:8]] + [(c, 'direct') for c in rng.sample(ks, 6)]
        okspl = [c for c in ks if rows[c['case_id']]['split_field_status'] != 'lower-bound-p (proof)']
        pick += [(c, 'split') for c in rng.sample(okspl, 6 if k == 1 else 4)]
    out = EXP / 'verify' / 'ordinary.jsonl'
    with open(out, 'w') as f:
        for c, enc in pick:
            p, k = c['prime'], c['k']
            sysd = V.build_system(enc, p, c['a'], c['b'], c['target'][0], k, terms)
            comps = V.quotient_algebra(sysd, p, k, V.subgroup(p, k))
            inputs, r = ordinary(sysd, p, k)
            pts = points(sysd, p, k, comps)
            t0 = time.time()
            res = L.exact_lfd(inputs, r, p, pts, Dmax=(9 if enc == 'direct' else 7), deadline=time.time() + 240)
            R = rows[c['case_id']]
            theirs = R['direct_lfd'] if enc == 'direct' else R['split_free_lfd']
            rec = dict(case_id=c['case_id'], k=k, enc=enc, p=p, validator_lfd=res.get('lfd'),
                       validator_largest_fall=res.get('lfd_largest_fall'), certified_at=res['certified_at'],
                       producer_lfd=theirs, agree=res.get('lfd') == theirs, seconds=round(time.time() - t0, 2))
            f.write(json.dumps(rec) + '\n'); f.flush(); print(rec, flush=True)


if __name__ == '__main__':
    main()
