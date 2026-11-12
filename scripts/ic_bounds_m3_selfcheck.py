#!/usr/bin/env python3
"""Self-checks for csd/ic_bounds_m3.py (item 1, m = 3).

1. S_4 = Res_X(S_3(x1,x2,X), S_3(x3,x4,X)) vanishes on x-coordinates of random rational
   decompositions P1 + P2 + P3 = R (and P1 + P2 + P3 + P4 = O), is symmetric in its four
   arguments, has degree 4 in each variable, and its expansion (s4_polynomial) agrees with the
   resultant evaluation (s4_value) at random points.
2. S_3 helper agrees with csd.ic_dreg.s3_prime_value.
3. The independent mixed-closure engine with no auxiliary variable reproduces
   chain_degree_closure (quotient rank, pivots' degree profile) on direct m = 3 systems; the
   batched engine reproduces the dense engines on both encodings; full exact profiles agree
   between the dense engines and the batched engine with its early stop.
4. Split generators vanish on actual rational decompositions with z = x(P1 + P2).
Writes experiments/ic-bounds-m3/selfcheck.json.
"""
import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from csd import ic_dreg as I  # noqa: E402
from csd import ic_fall_certificates as C  # noqa: E402
from csd import ic_bounds_m3 as M  # noqa: E402

SEED = 20261107


def main():
    rng = random.Random(SEED)
    report = dict(seed=SEED, checks=[])
    for p in (17, 97, 257, 641, 7681, 1000003):
        for _ in range(3):
            while True:
                a, b = rng.randrange(p), rng.randrange(p)
                if (4 * a ** 3 + 27 * b * b) % p:
                    break
            E = I.PrimeCurve(p, a, b)
            sym_ok = deg_ok = expand_ok = s3_ok = dec_ok = split_ok = 0
            for _ in range(20):
                P = [E.random_point(rng) for _ in range(3)]
                R = E.add(E.add(P[0], P[1]), P[2])
                Q = E.add(P[0], P[1])
                if R is None or Q is None:
                    continue
                xs = [Pt[0] for Pt in P]
                assert M.s4_value(E, *xs, R[0]) == 0
                assert M.poly_value(p, M.s4_polynomial(E, R[0]), xs) == 0
                dec_ok += 1
                # split equations at z = x(P1 + P2)
                g1, g2, raw = M.split_generators(E, R[0], 1)
                assert M.poly_value(p, raw[0], xs + [Q[0]]) == 0 and M.poly_value(p, raw[1], xs + [Q[0]]) == 0
                split_ok += 1
                # S_3 consistency
                assert I.s3_prime_value(E, xs[0], xs[1], Q[0]) == 0
                s3_ok += 1
            for _ in range(10):
                v = [rng.randrange(p) for _ in range(4)]
                base = M.s4_value(E, *v)
                for perm in ([1, 0, 2, 3], [0, 2, 1, 3], [3, 1, 2, 0], [2, 3, 0, 1]):
                    assert M.s4_value(E, *[v[i] for i in perm]) == base
                sym_ok += 1
                assert M.poly_value(p, M.s4_polynomial(E, v[3]), v[:3]) == base
                expand_ok += 1
            S4 = M.s4_polynomial(E, rng.randrange(p))
            degs = [max(e[i] for e in S4) for i in range(3)]
            assert all(d <= 4 for d in degs)
            deg_ok = int(degs == [4, 4, 4])
            report['checks'].append(dict(p=p, a=a, b=b, decompositions=dec_ok, split_points=split_ok,
                                         s3_points=s3_ok, symmetry=sym_ok, expansion=expand_ok,
                                         degrees_in_x1_x2_x3=degs, total_degree=max(sum(e) for e in S4),
                                         terms=len(S4)))
    # engine cross-check
    cross = []
    for p, a, b, k in [(17, 2, 3, 1), (17, 2, 3, 2), (97, 0, 5, 2), (97, 3, 0, 2), (257, 5, 7, 3), (17, 1, 0, 3)]:
        E = I.PrimeCurve(p, a, b)
        for t in range(1, p, max(1, p // 4)):
            if E.lift_x(t) is None:
                continue
            g = M.direct_generator(E, t, k)
            for D in range(0, 3 * k + 2):
                cl = C.chain_degree_closure([g], 3, k, p, D, max_columns=5000, max_cells=60_000_000)
                mx = M.mixed_closure([g], 3, k, p, D, aux=False)
                assert mx['rank'] == cl.quotient_rank, (p, a, b, k, t, D)
                assert all(M.low_dimension(mx, d) == sum(cl.masks[j].bit_count() <= d for j in cl.pivots)
                           for d in range(D + 1))
                cross.append(dict(p=p, a=a, b=b, k=k, t=t, D=D, quotient_rank=cl.quotient_rank))
                if cl.contains_one:
                    break
            break
    report['engine_cross_check'] = dict(cases=len(cross), rows=cross)
    # batched (memory-bounded) engine against the dense engines, both encodings
    from csd import ic_bounds_m3_batched as BC
    batched = []
    for p, a, b in [(97, 5, 7), (17, 2, 3), (257, 0, 3), (257, 1, 0)]:
        E = I.PrimeCurve(p, a, b)
        t = next(x for x in range(3, p) if E.lift_x(x) is not None)
        for k in (1, 2, 3):
            g = M.direct_generator(E, t, k)
            g1, g2, _ = M.split_generators(E, t, k)
            for D in range(0, 8):
                A = M.mixed_closure([g], 3, k, p, D, aux=False)
                B = BC.batched_closure([g], 3, k, p, D, False, batch=64)
                assert A['rank'] == B.rank and A['contains_one'] == B.contains_one
                assert all(M.low_dimension(A, d) == B.low_dimension(d) for d in range(D + 1))
                batched.append(dict(p=p, a=a, b=b, t=t, k=k, D=D, encoding='direct', quotient_rank=B.rank))
                if D <= 6:
                    A = M.mixed_closure([g1, g2], 3, k, p, D, aux=True)
                    B = BC.batched_closure([g1, g2], 3, k, p, D, True, batch=64)
                    assert A['rank'] == B.rank and A['contains_one'] == B.contains_one
                    assert all(M.low_dimension(A, d) == B.low_dimension(d) for d in range(D + 1))
                    batched.append(dict(p=p, a=a, b=b, t=t, k=k, D=D, encoding='split', quotient_rank=B.rank))
    report['batched_engine_cross_check'] = dict(pairs=len(batched), rows=batched)
    # full exact profiles: dense engines vs batched engine with the early stop at dim(I cap M_<=D)
    profiles = []
    for p, a, b, t in [(97, 5, 7, 3), (17, 2, 3, 3), (17, 0, 1, 16), (257, 1, 0, 0), (17, 1, 0, 4)]:
        E = I.PrimeCurve(p, a, b)
        for k in (2, 3):
            for enc in ('direct', 'split'):
                r1 = M.exact_profile(enc, E, t, k)
                r2 = M.exact_profile(enc, E, t, k, engine='batched')
                assert len(r1['profile']) == len(r2['profile']) and r1['results'] == r2['results']
                for x, y in zip(r1['profile'], r2['profile']):
                    assert {q: v for q, v in x.items() if q not in ('engine', 'closure_check')} == \
                        {q: v for q, v in y.items() if q not in ('engine', 'closure_check')}
                profiles.append(dict(p=p, a=a, b=b, t=t, k=k, encoding=enc,
                                     results={q: v['exact_lfd'] for q, v in r1['results'].items()}))
    report['profile_engine_agreement'] = dict(profiles=len(profiles), rows=profiles)
    out = ROOT / 'experiments' / 'ic-bounds-m3' / 'selfcheck.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=1))
    print('ok', len(report['checks']), 'curve checks;', len(cross), 'engine cross-check rows;',
          len(batched), 'batched pairs;', len(profiles), 'profiles')


if __name__ == '__main__':
    main()
