#!/usr/bin/env python3
"""Independent re-computation of item-1 LFDs with the validator engine
experiments/ic-curve-validation/code/lfd_engine.py (ordinary polynomial ring, full monomial space,
no chain-quotient representation, no import of the closure engines used for the results).

Direct m = 3 systems: inputs = 3k chain equations u_{i,j}^2 - u_{i,j+1}, u_{i,k-1}^2 - 1 and the
normalized generator NF(S_4) written in the ordinary ring; points = chain lifts of Z (all
rational, ideal radical).  Split systems: 3k chain variables plus z; validated only where
dim R/I_free = dim R/I_field (all auxiliary roots rational and simple, so the ideal is radical with
rational zero set {(x, z)}), points enumerated as (chain lift of x, z).
Sample: every k = 1 case and every 4th k = 2 case (order of cases.jsonl).  Writes validation.json.
"""
import importlib.util
import json
import sys
import time
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from csd import ic_dreg as I  # noqa: E402
from csd import ic_bounds_m3 as M  # noqa: E402

spec = importlib.util.spec_from_file_location('lfd_engine', ROOT / 'experiments/ic-curve-validation/code/lfd_engine.py')
L = importlib.util.module_from_spec(spec)
spec.loader.exec_module(L)
OUT = ROOT / 'experiments' / 'ic-bounds-m3'


def chain_inputs(k, p, r):
    inputs = []
    for i in range(3):
        for j in range(k):
            sq = [0] * r
            sq[i * k + j] = 2
            f = {tuple(sq): 1}
            if j + 1 < k:
                nxt = [0] * r
                nxt[i * k + j + 1] = 1
                f[tuple(nxt)] = p - 1
            else:
                f[tuple([0] * r)] = p - 1
            inputs.append(f)
    return inputs


def to_exp(mask, e_z, r, k):
    e = [(mask >> b) & 1 for b in range(3 * k)]
    if r > 3 * k:
        e.append(e_z)
    return tuple(e)


def lift(E3, k, g0, p):
    out = []
    for e in E3:
        x = pow(g0, e, p)
        out += [pow(x, 1 << j, p) for j in range(k)]
    return out


def validate_case(c, rows):
    p, k = c['prime'], c['k']
    E = I.PrimeCurve(p, c['a'], c['b'])
    t = c['target'][0]
    g0, H = M.subgroup(p, k)
    N = 1 << k
    rec = dict(case_id=c['case_id'], k=k, prime=p, a=c['a'], b=c['b'], target=c['target'])
    # direct
    r = 3 * k
    g = M.direct_generator(E, t, k)
    G = {}
    for ms, cf in g.items():
        G[to_exp(ms, 0, r, k)] = cf
    inputs = chain_inputs(k, p, r) + ([G] if G else [])
    pts = []
    for E3 in product(range(N), repeat=3):
        if M.nf_value(p, k, g, E3, g0) == 0:
            pts.append(tuple(lift(E3, k, g0, p)))
    res = L.exact_lfd(inputs, r, p, pts, 3 * k + 1)
    rec['direct_validator'] = res['lfd']
    rec['direct_validator_largest_fall'] = res.get('lfd_largest_fall')
    rec['direct_ours'] = rows[c['case_id']]['direct_lfd']
    # split (only when I_free = I_field)
    row = rows[c['case_id']]
    if row.get('split_free_dim') is not None and row['split_free_dim'] == row.get('split_field_dim'):
        r = 3 * k + 1
        g1, g2, _ = M.split_generators(E, t, k)
        gens = []
        for gg in (g1, g2):
            d = {}
            for (ms, ez), cf in gg.items():
                d[to_exp(ms, ez, r, k)] = cf
            gens.append(d)
        loc = M.split_local_algebra(g1, g2, k, p, g0, True)
        spts = []
        for pt, h in loc.items():
            E3 = [(pt >> (i * k)) & (N - 1) for i in range(3)]
            roots = [z for z in range(p) if sum(cf * pow(z, i, p) for i, cf in enumerate(h)) % p == 0]
            assert len(roots) == len(h) - 1
            for z in roots:
                spts.append(tuple(lift(E3, k, g0, p) + [z]))
        res = L.exact_lfd(chain_inputs(k, p, r) + gens, r, p, spts, 3 * k + 4)
        rec['split_validator'] = res['lfd']
        rec['split_ours'] = row['split_free_lfd']
    return rec


def main():
    cases = [json.loads(l) for l in open(OUT / 'cases.jsonl')]
    rows = {json.loads(l)['case_id']: json.loads(l) for l in open(OUT / 'rows.jsonl')}
    sample = [c for c in cases if c['k'] == 1] + [c for c in cases if c['k'] == 2][::4]
    out = []
    t0 = time.time()
    for c in sample:
        rec = validate_case(c, rows)
        out.append(rec)
        print(rec, flush=True)
        if time.time() - t0 > 20 * 60:
            print('time budget reached')
            break
    agree_d = sum(r['direct_validator'] == r['direct_ours'] for r in out)
    sp = [r for r in out if 'split_validator' in r]
    agree_s = sum(r['split_validator'] == r['split_ours'] for r in sp)
    summary = dict(engine='experiments/ic-curve-validation/code/lfd_engine.py', cases=len(out),
                   direct_agree=agree_d, split_checked=len(sp), split_agree=agree_s, rows=out)
    (OUT / 'validation.json').write_text(json.dumps(summary, indent=1))
    print('direct agree', agree_d, '/', len(out), '; split agree', agree_s, '/', len(sp))


if __name__ == '__main__':
    main()
