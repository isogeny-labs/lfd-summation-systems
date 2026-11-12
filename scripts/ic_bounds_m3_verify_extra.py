#!/usr/bin/env python3
"""Verifier extra checks (independent code):
 1. case data: p prime, 2^k | p-1, target on curve, planted P1+P2+P3 = R with x_i in H and (x1,x2,x3) in Z;
 2. S_4(x1,x2,x3,t) has 125 terms for random (a,b,t) mod a large prime;
 3. Corollary 2 hypotheses evaluated directly per case (k>=2): deg NF(g) = 6 <=> c333 != 0;
    deg NF(x1*g) <= 6; g not supported in one level set {x1 = c}; and x1*g not proportional to g;
 4. -> verify/extra.json
"""
import json, itertools, random
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ic_bounds_m3_verify_engine as V
from ic_bounds_m3_verify_symbolic import add

EXP = Path(__file__).resolve().parents[1] / 'experiments' / 'ic-bounds-m3'


def main():
    terms = V.s4_terms()
    cases = [json.loads(l) for l in open(EXP / 'cases.jsonl')]
    rows = {json.loads(l)['case_id']: json.loads(l) for l in open(EXP / 'rows.jsonl')}
    rep = dict(case_data_ok=0, planted=0, planted_in_Z=0, planted_sum_ok=0, cor=[], problems=[])
    for c in cases:
        p, a, b, k = c['prime'], c['a'], c['b'], c['k']
        t, y = c['target']
        ok = V.is_prime(p) and (p - 1) % (1 << k) == 0 and (y * y - (t**3 + a * t + b)) % p == 0 \
            and (4 * a**3 + 27 * b * b) % p != 0
        rep['case_data_ok'] += ok
        if not ok: rep['problems'].append(('case', c['case_id']))
        H = V.subgroup(p, k); Hs = set(H)
        g = V.nf(V.poly_s4_at(p, a, b, t, terms), 3, k, False, p)
        def gval(x): return sum(cc * V.mono_val(mm, x, p, 3) for mm, cc in g.items()) % p
        if c.get('planted'):
            rep['planted'] += 1
            P = c['planted']
            pts = [tuple(q) for q in P] if isinstance(P, list) else [tuple(q) for q in P['points']]
            xs = [q[0] for q in pts]
            if all(x in Hs for x in xs) and gval(tuple(xs)) == 0: rep['planted_in_Z'] += 1
            else: rep['problems'].append(('planted', c['case_id']))
            S = None
            for q in pts:
                assert (q[1]**2 - q[0]**3 - a * q[0] - b) % p == 0
                S = add(S, q, a, p)
            if S is not None and tuple(S) == (t % p, y % p): rep['planted_sum_ok'] += 1
            else: rep['problems'].append(('planted-sum', c['case_id']))
        if k < 2: continue
        deg = max(V.weight(mm, 3, False) for mm in g) if g else None
        c333 = (-8 * (5 * t**3 + 3 * a * t + 6 * b)) % p
        N = 1 << k
        ug = {}
        for mm, cc in g.items():
            key = ((mm[0] + 1) % N, mm[1], mm[2]); ug[key] = (ug.get(key, 0) + cc) % p
        ug = {kk: v for kk, v in ug.items() if v}
        degu = max(V.weight(mm, 3, False) for mm in ug) if ug else None
        vals = {x: gval(x) for x in itertools.product(H, repeat=3)}
        support_x1 = sorted({x[0] for x, v in vals.items() if v})
        prop = False
        nz = [x for x, v in vals.items() if v]
        if nz:
            x0 = nz[0]; lam = x0[0]
            prop = all((x[0] * v - lam * v) % p == 0 for x, v in vals.items())
        R = rows[c['case_id']]
        rec = dict(case_id=c['case_id'], k=k, p=p, deg=deg, c333_nonzero=c333 != 0, degNF_x1g=degu,
                   n_support_slices_x1=len(support_x1), x1g_proportional=prop,
                   producer_deg=R['generator_degree'], producer_c333=R['c333'] % p,
                   producer_applies=R['carry_lower_bound_7_applies'], producer_shift_prop=R['shift_proportional'])
        hyp = deg == 6 and degu is not None and degu <= 6 and len(support_x1) >= 2
        rec['hypotheses_hold'] = hyp
        if (deg == 6) != (c333 != 0): rep['problems'].append(('deg-c333', c['case_id']))
        if deg != R['generator_degree'] or c333 != R['c333'] % p: rep['problems'].append(('producer-mismatch', c['case_id']))
        if hyp != R['carry_lower_bound_7_applies'] and k >= 3: rep['problems'].append(('applies-mismatch', c['case_id']))
        if hyp and R['direct_lfd'] < 7: rep['problems'].append(('COUNTEREXAMPLE', c['case_id']))
        rep['cor'].append(rec)
    # 125 terms
    rng = random.Random(5); P = 1000003
    cnt = []
    for _ in range(5):
        a, b, t = rng.randrange(P), rng.randrange(P), rng.randrange(P)
        cnt.append(len(V.poly_s4_at(P, a, b, t, terms)))
    rep['terms_at_x4_t'] = cnt
    from collections import Counter
    rep['hyp_hold_by_k'] = dict(Counter((r['k'], r['hypotheses_hold']) for r in rep['cor']).__repr__() and
                                {f"k{k}-{h}": n for (k, h), n in Counter((r['k'], r['hypotheses_hold']) for r in rep['cor']).items()})
    rep['min_support_slices_when_deg6_k>=3'] = min(r['n_support_slices_x1'] for r in rep['cor'] if r['deg'] == 6 and r['k'] >= 3)
    rep['max_degNF_x1g'] = max(r['degNF_x1g'] or 0 for r in rep['cor'])
    rep['any_proportional'] = any(r['x1g_proportional'] for r in rep['cor'])
    summary = {kk: v for kk, v in rep.items() if kk != 'cor'}
    (EXP / 'verify' / 'extra.json').write_text(json.dumps(rep, indent=0))
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
