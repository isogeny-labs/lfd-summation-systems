"""E2: full (A,B)-invariance over F_p: dim V_4(mu^2 a, mu^3 b, mu t) for one representative mu of
every coset of H in F_p^*, k = 6 (192 cosets at p = 12289) and k = 7 (96 cosets).
Pivot sets are compared with mu = 1 (phi_mu is diagonal, so equal closures have equal pivots)."""
import json, sys, time
sys.path.insert(0, 'experiments/ic-bounds-genericity')
import genericity_engine as G

p = 12289
res = {}
for (a, b, t, tag) in [(3899, 9710, 9895, 'r3'), (1, 0, 9362, 'j1728')]:
    for k in (6, 7):
        tg = G.TwistGrid(k, 4)
        reps = G.coset_representatives(p, k)
        base = G.closure(G.generator(a, b, t, p, k), tg, p, 1)
        hist = {}
        bad = []
        t0 = time.time()
        for mu in reps:
            r = G.closure(G.generator(mu * mu * a % p, pow(mu, 3, p) * b % p, mu * t % p, p, k), tg, p, 1)
            key = (r['rank'], tuple(r['profile']), r['pivots'] == base['pivots'], r['closed'])
            hist[str(key)] = hist.get(str(key), 0) + 1
            if r['pivots'] != base['pivots'] or r['rank'] != base['rank']:
                bad.append(dict(mu=mu, dim=r['rank'], profile=r['profile']))
        res[f'{tag}-k{k}'] = dict(cosets=len(reps), base_dim=base['rank'], base_profile=base['profile'],
                                  histogram=hist, deviating=bad, seconds=round(time.time() - t0, 1))
        print(tag, k, res[f'{tag}-k{k}'], flush=True)
json.dump(res, open('experiments/ic-bounds-genericity/coset_test.json', 'w'), indent=1)
