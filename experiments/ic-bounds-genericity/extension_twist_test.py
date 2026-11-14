"""E4/E5/E10: (a) the free extension W_8 of V_4(7) is closed under the twisted shifts of A_8^{(c)}
for c in {0, 1, 2, random}, and equals V_4^{(c)}(8); (b) pivot sets of V_4(6), V_4(7) agree across
primes 193, 12289, 40961, 33554393 at random points; (c) density of the RREF bases."""
import json, sys, random, time
sys.path.insert(0, 'experiments/ic-bounds-genericity')
import numpy as np
import genericity_engine as G
from csd import elim
from csd.ic_bounds_fall import free_extension, make_grid

out = {}
# (a) twisted free-extension certificate at depth 8, primary instance
p, a, b, t = 40961, 21223, 9887, 4748
tg7, tg8 = G.TwistGrid(7, 4), G.TwistGrid(8, 4)
V7 = G.closure(G.generator(a, b, t, p, 7), tg7, p, 1)
res7 = dict(basis=V7['basis'], pivots=V7['pivots'])
W8, info = free_extension(res7, tg7.grid, tg8.grid, p)
rank8 = elim.fp_echelon(W8, p)[0]
rows = []
for c in (1, 0, 2, 777):
    t0 = time.time()
    outb = []
    r, piv, _ = elim.fp_echelon(W8, p, outb)
    B = outb[0]
    low = [i for i, col in enumerate(piv) if tg8.weights[col] < 4]
    closed = True
    for v in range(16):
        stacked = np.vstack([B, G.twisted_shift(B[low], tg8, v, c, p)])
        if elim.fp_echelon(stacked, p)[0] != r:
            closed = False
    V8c = G.closure(G.generator(a, b, t, p, 8), tg8, p, c)
    rel = G.row_space_relation(W8, V8c['basis'], p)
    rows.append(dict(c=c, dim_W8=int(r), W8_closed_twisted=closed, dim_V8c=V8c['rank'], profile_V8c=V8c['profile'],
                     V8c_equals_W8=rel['equal'], one_in_W8=(tg8.const_column in piv), seconds=round(time.time() - t0, 1)))
    print(rows[-1], flush=True)
out['depth8_twisted_free_extension'] = dict(predicted_dim=info['predicted_dim'], rows=rows)

# (b) pivot sets across primes
piv_ref = {}
rows = []
for k in (6, 7):
    tg = G.TwistGrid(k, 4)
    for q in (193, 12289, 40961, 33554393):
        rng = random.Random(q + k)
        for trial in range(2):
            A, B = rng.randrange(1, q), rng.randrange(1, q)
            r = G.closure(G.generator(A, B, 1, q, k), tg, q, 1)
            key = f'k{k}'
            if key not in piv_ref:
                piv_ref[key] = r['pivots']
            rows.append(dict(k=k, prime=q, A=A, B=B, dim=r['rank'], profile=r['profile'], pivots_equal_reference=r['pivots'] == piv_ref[key],
                             symmetric_difference=len(set(r['pivots']) ^ set(piv_ref[key]))))
            print(rows[-1], flush=True)
out['pivots_across_primes'] = rows

# (c) density
dens = []
for k, q in ((6, 40961), (7, 40961), (6, 33554393), (7, 33554393)):
    tg = G.TwistGrid(k, 4)
    rng = random.Random(k)
    r = G.closure(G.generator(rng.randrange(1, q), rng.randrange(1, q), 1, q, k), tg, q, 1)
    Bm = r['basis']
    nnz = int(np.count_nonzero(Bm))
    low = [i for i, col in enumerate(r['pivots']) if tg.weights[col] <= 3]
    dens.append(dict(k=k, prime=q, shape=list(Bm.shape), nonzero=nnz, density=round(nnz / Bm.size, 4),
                     nonzero_low_rows=int(np.count_nonzero(Bm[low])), low_rows=len(low),
                     weight_le3_columns=int(np.sum(tg.weights <= 3))))
    print(dens[-1], flush=True)
out['density'] = dens
json.dump(out, open('experiments/ic-bounds-genericity/extension_twist_test.json', 'w'), indent=1)
