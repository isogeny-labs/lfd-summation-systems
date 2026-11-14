"""Negative control for the exact test: same chains, same pipeline, but a generator outside the
hypothesis of Theorem 1.  The generator is a random element of chain degree 2 (all squarefree chain
monomials of degree <= 2, i.e. x1^e1 x2^e2 with wt(e1)+wt(e2) <= 2; this includes x_i^3 = u_{i,0}u_{i,1},
so the fibre polynomials have degree 3 when the depth is >= 2).  If the pipeline can report an exact
LFD above min(k1,k2)+2 it does so here; for S_3 it never does (rt_run.py)."""
import itertools, random
from rt_lib import *

def control(p, k1, k2, rng):
    r = k1 + k2
    inputs, _, _ = build(p, 1, 1, 0, k1, k2)
    inputs = inputs[:-1]                                  # chains only
    mons = [tuple([0] * r)]
    for i in range(r):
        m = [0] * r; m[i] = 1; mons.append(tuple(m))
    for i, j in itertools.combinations(range(r), 2):
        m = [0] * r; m[i] = 1; m[j] = 1; mons.append(tuple(m))
    G = {m: rng.randrange(p) for m in mons}
    G = {m: c for m, c in G.items() if c}
    inputs.append(G)
    H1, H2 = mu(p, 1 << k1), mu(p, 1 << k2)
    pts = []
    for x1 in H1:
        for x2 in H2:
            u = [pow(x1, 1 << j, p) for j in range(k1)] + [pow(x2, 1 << j, p) for j in range(k2)]
            v = 0
            for m, c in G.items():
                t = c
                for j, e in enumerate(m):
                    if e: t = t * u[j] % p
                v = (v + t) % p
            if v == 0: pts.append(tuple(u))
    res = E.exact_lfd(inputs, r, p, pts, min(k1 + k2 + 1, 6))
    return res['lfd']

rng = random.Random(4242)
for (p, k1, k2, n) in ((17, 1, 2, 300), (17, 1, 3, 300), (17, 1, 4, 300), (13, 2, 2, 300), (17, 2, 2, 300), (17, 2, 3, 150), (17, 2, 4, 60)):
    hist = {}
    for _ in range(n):
        lfd = control(p, k1, k2, rng)
        hist[str(lfd)] = hist.get(str(lfd), 0) + 1
    b = min(k1, k2) + 2
    over = sum(c for l, c in hist.items() if l == 'None' or int(l) > b)
    print(f'p={p} (k1,k2)=({k1},{k2}) random degree-2 generators={n} exact LFD histogram={dict(sorted(hist.items()))} min+2={b} dimension bound={k1+k2+1} above min+2: {over}', flush=True)
