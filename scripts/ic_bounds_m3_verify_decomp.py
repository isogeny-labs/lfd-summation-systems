#!/usr/bin/env python3
"""Every rational decomposition P1+P2+P3 = R (target R = (t, y_R), x_i in H, P_i in E(F_p)) has
(x1,x2,x3) in Z = zeros of NF(S_4(.,.,.,t)) on H^3.  Independent group law; all 461 cases.
-> verify/decomp.json"""
import json, itertools, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ic_bounds_m3_verify_engine as V
from ic_bounds_m3_verify_symbolic import add, sqrt_mod
EXP = Path(__file__).resolve().parents[1] / 'experiments' / 'ic-bounds-m3'
terms = V.s4_terms()
rows = {json.loads(l)['case_id']: json.loads(l) for l in open(EXP / 'rows.jsonl')}
out = dict(cases=0, triples=0, triples_in_Z=0, count_matches_producer=0, problems=[])
for l in open(EXP / 'cases.jsonl'):
    c = json.loads(l); p, a, b, k = c['prime'], c['a'], c['b'], c['k']; R = tuple(c['target'])
    H = V.subgroup(p, k)
    g = V.nf(V.poly_s4_at(p, a, b, R[0], terms), 3, k, False, p)
    pts = {}
    for x in H:
        r = (x**3 + a * x + b) % p
        if r == 0: pts[x] = [(x, 0)]
        elif pow(r, (p - 1) // 2, p) == 1:
            y = sqrt_mod(r, p); pts[x] = [(x, y), (x, p - y)]
    trip = set()
    for x1, x2, x3 in itertools.product(pts, repeat=3):
        for P1 in pts[x1]:
            for P2 in pts[x2]:
                S = add(P1, P2, a, p)
                for P3 in pts[x3]:
                    if add(S, P3, a, p) == R: trip.add((x1, x2, x3))
    out['cases'] += 1; out['triples'] += len(trip)
    inZ = sum(1 for x in trip if sum(cc * V.mono_val(mm, x, p, 3) for mm, cc in g.items()) % p == 0)
    out['triples_in_Z'] += inZ
    if inZ != len(trip): out['problems'].append(c['case_id'])
    out['count_matches_producer'] += (len(trip) == rows[c['case_id']]['rational_ordered'])
(EXP / 'verify' / 'decomp.json').write_text(json.dumps(out, indent=1)); print(out)
