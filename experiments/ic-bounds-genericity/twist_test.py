"""E1: c-dependence of the degree-4 closure in A_k^{(c)} (wrap factor c), k = 6, 7.
c = 1 is the chain algebra, c = 0 the no-wrap truncation, other c are the cosets mu^{-1} H.
Records dim, profile, pivot-set agreement and row-space agreement with c = 1."""
import json, sys, time, random
sys.path.insert(0, 'experiments/ic-bounds-genericity')
import numpy as np
import genericity_engine as G

instances = [
    dict(tag='p40961-s7', p=40961, a=21223, b=9887, t=4748),
    dict(tag='p12289-r3', p=12289, a=3899, b=9710, t=9895),
    dict(tag='p12289-j1728', p=12289, a=1, b=0, t=9362),
    dict(tag='p193-generic', p=193, a=22, b=23, t=1),
    dict(tag='p193-ordertwo', p=193, a=22, b=(-1 - 22) % 193, t=1),   # 1+A+B=0
    dict(tag='p193-tangent', p=193, a=190, b=23, t=1),                # A=-3
]
out = {}
for inst in instances:
    p = inst['p']
    for k in (6, 7):
        tg = G.TwistGrid(k, 4)
        g = G.generator(inst['a'], inst['b'], inst['t'], p, k)
        rng = random.Random(1)
        cs = [1, 0, 2, 3, p - 1] + [rng.randrange(2, p) for _ in range(4)]
        base = None
        rows = []
        for c in cs:
            t0 = time.time()
            r = G.closure(g, tg, p, c)
            if base is None:
                base = r
            rel = G.row_space_relation(base['basis'], r['basis'], p)
            rows.append(dict(c=c, dim=r['rank'], profile=r['profile'], closed=r['closed'], one=r['contains_one'],
                             rounds=r['rounds'], same_pivots=(r['pivots'] == base['pivots']),
                             same_space_as_c1=rel['equal'], a_in_b=rel['a_in_b'], b_in_a=rel['b_in_a'],
                             seconds=round(time.time() - t0, 2)))
            print(inst['tag'], k, rows[-1], flush=True)
        out[f"{inst['tag']}-k{k}"] = rows
json.dump(out, open('experiments/ic-bounds-genericity/twist_test.json', 'w'), indent=1)
