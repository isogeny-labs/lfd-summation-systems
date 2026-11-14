"""Item (c): the univariate chain system {u_j^2 - u_{j+1}, u_{k-1}^2 - 1, q0 + q1 X + q2 X^2 reduced}
has LFD <= 2 for EVERY nonzero q of degree <= 2 (all degenerate cases included), exact, ordinary ring."""
from rt_lib import *
tot = 0; bad = []
for p, ks in ((5, (1, 2)), (13, (1, 2)), (17, (1, 2, 3, 4)), (41, (3,)), (97, (5,))):
    for k in ks:
        N = 1 << k
        H = mu(p, N)
        hist = {}; cases = {}
        qs = [(q0, q1, q2) for q0 in range(p) for q1 in range(p) for q2 in range(p) if (q0, q1, q2) != (0, 0, 0)]
        if p >= 41:
            # all monic / degenerate shapes: (X-b1)(X-b2) for all b1<=b2, irreducibles sampled, linear, constant
            qs = [((b1 * b2) % p, (-(b1 + b2)) % p, 1) for b1 in range(p) for b2 in range(b1, p)]
            qs += [((-b1) % p, 1, 0) for b1 in range(p)] + [(1, 0, 0)]
            nr = [x for x in range(2, p) if pow(x, (p - 1) // 2, p) == p - 1][0]
            qs += [((s * s - nr * d * d) % p, (-2 * s) % p, 1) for s in range(p) for d in range(1, p, 7)]   # irreducible
        for (q0, q1, q2) in qs:
            inputs = []
            for j in range(k):
                sq = [0] * k; sq[j] = 2
                f = {tuple(sq): 1}
                if j + 1 < k:
                    nx = [0] * k; nx[j + 1] = 1; f[tuple(nx)] = p - 1
                else:
                    f[tuple([0] * k)] = p - 1
                inputs.append(f)
            g = {}
            def addm(m, c):
                g[m] = (g.get(m, 0) + c) % p
            addm(tuple([0] * k), q0)
            m1 = [0] * k; m1[0] = 1; addm(tuple(m1), q1)
            if k >= 2:
                m2 = [0] * k; m2[1] = 1; addm(tuple(m2), q2)
            else:
                addm(tuple([0] * k), q2)
            g = {m: c for m, c in g.items() if c}
            inputs.append(g)
            roots = [x for x in H if (q0 + q1 * x + q2 * x * x) % p == 0]
            pts = [tuple(pow(x, 1 << j, p) for j in range(k)) for x in roots]
            res = E.exact_lfd(inputs, k, p, pts, 4)
            deg = 2 if q2 else (1 if q1 else 0)
            disc0 = deg == 2 and (q1 * q1 - 4 * q0 * q2) % p == 0
            ck = (deg, len(roots), 'double' if disc0 else '')
            cases[ck] = cases.get(ck, 0) + 1
            hist[res['lfd']] = hist.get(res['lfd'], 0) + 1
            tot += 1
            if res['lfd'] is None or res['lfd'] > 2:
                bad.append((p, k, q0, q1, q2, len(roots), res['lfd']))
        print(f'p={p} k2={k} polynomials={len(qs)} lfd_hist={hist} cases(deg,roots in H,double)={dict(sorted(cases.items()))}', flush=True)
print('total', tot, 'with LFD > 2:', bad[:20], len(bad))
