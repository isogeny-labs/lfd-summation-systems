"""Singular parameter values (4a^3 + 27b^2 = 0): outside T1, inside the kernel bound 2k+1 (which
holds for all a, b, t).  Is k+2 still an upper bound there?  Exact ordinary-ring LFD."""
import sys, json
from rt_lib import *
out = {}
worst = []
for p, k in ((5, 1), (13, 1), (5, 2), (13, 2), (17, 2), (29, 2), (17, 3), (41, 3), (17, 4)):
    hist = {}
    H = set(mu(p, 1 << k))
    for a in range(p):
        for b in range(p):
            if nonsingular(p, a, b):
                continue
            for t in range(p):
                if p == 41 and t not in H and t % 5:      # thin the p = 41 plane outside H
                    continue
                if k == 4 and t not in (0, 1, 2, 3) and (a, b) != (0, 0) and (t ** 3 + a * t + b) % p:
                    continue
                res, bound, z = exact(p, a, b, t, k, k, Dmax=2 * k + 1)
                # vanishing fibre: alpha = t in H with f(t) = f'(t) = 0
                vf = (t in H) and (t ** 3 + a * t + b) % p == 0 and (3 * t * t + a) % p == 0
                key = ('vanishing-fibre' if vf else 'no-vanishing-fibre', res['lfd'])
                hist[key] = hist.get(key, 0) + 1
                if res['lfd'] is None or res['lfd'] > k + 2:
                    worst.append((p, k, a, b, t, z, res['lfd'], vf))
    print(f'p={p} k={k} k+2={k+2} 2k+1={2*k+1}', {f'{a}:{b}': c for (a, b), c in sorted(hist.items(), key=str)}, flush=True)
print('instances with LFD > k+2 (p,k,a,b,t,|Z|,lfd,vanishing_fibre):')
for w in worst[:40]:
    print('  ', w)
print('count', len(worst), 'all of them have a vanishing fibre:', all(w[7] for w in worst))
