"""E3: degree of the pivot denominator D along a random line in the normalized plane.

(A, B) = (A0 + s A1, B0 + s B1), t = 1, over a prime p < 2^25.  The degree-4 closure is
computed at s = 1..N; the pivot set must be the generic one at every s.  Selected RREF
entries are then reconstructed as rational functions of s (Newton interpolation + extended
Euclid), and D|_line is estimated as the lcm of the reconstructed denominators.
"""
import argparse, json, random, sys, time
sys.path.insert(0, 'experiments/ic-bounds-genericity')
import numpy as np
import genericity_engine as G

ap = argparse.ArgumentParser()
ap.add_argument('--k', type=int, default=6)
ap.add_argument('--N', type=int, default=1200)
ap.add_argument('--prime', type=int, default=33554393)
ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--entries', type=int, default=6, help='entries sampled per stored row')
ap.add_argument('--out', default='experiments/ic-bounds-genericity/line_degree_k6.json')
ap.add_argument('--store', default='/tmp/claude-0/-home-user/a6e79f2a-a5ec-5904-b377-343a14909fff/scratchpad/line_rows.npy')
args = ap.parse_args()
p, k, N = args.prime, args.k, args.N
rng = random.Random(args.seed)
A0, B0, A1, B1 = [rng.randrange(1, p) for _ in range(4)]
tg = G.TwistGrid(k, 4)
cols = len(tg.masks)

def point(s):
    return (A0 + s * A1) % p, (B0 + s * B1) % p

t0 = time.time()
base = G.closure(G.generator(*point(1), 1, p, k), tg, p, 1)
piv = base['pivots']
w = [tg.weights[c] for c in piv]
sel = [i for i, x in enumerate(w) if x == 2]
w3 = [i for i, x in enumerate(w) if x == 3]
w4 = [i for i, x in enumerate(w) if x == 4]
sel += [w3[i] for i in sorted(set(int(j) for j in np.linspace(0, len(w3) - 1, 28)))]
sel += [w4[i] for i in sorted(set(int(j) for j in np.linspace(0, len(w4) - 1, 30)))]
R = len(sel)
store = np.zeros((N, R, cols), dtype=np.int32)
store[0] = base['basis'][sel]
bad = []
for n in range(2, N + 1):
    r = G.closure(G.generator(*point(n), 1, p, k), tg, p, 1, verify=False)
    if r['pivots'] != piv:
        bad.append(n)
        continue
    store[n - 1] = r['basis'][sel]
closure_seconds = time.time() - t0
print(f'closures done: {closure_seconds:.0f}s, dim {base["rank"]}, profile {base["profile"]}, bad points {bad}', flush=True)
np.save(args.store, store)

# --- rational reconstruction -------------------------------------------------------------
xs = np.array([n for n in range(1, N + 1) if n not in bad], dtype=np.int64)
good = np.array([n - 1 for n in range(1, N + 1) if n not in bad])
M = len(xs)
# sample entries: per row, `entries` columns nonzero at s=1 (plus the pivot column is 1, skip it)
ent = []
rng2 = random.Random(args.seed + 7)
for ri, i in enumerate(sel):
    nz = [c for c in np.flatnonzero(store[0, ri]) if c != piv[i]]
    for c in rng2.sample(nz, min(args.entries, len(nz))):
        ent.append((ri, int(c)))
E = len(ent)
Y = np.stack([store[good, ri, c].astype(np.int64) for ri, c in ent], axis=1)  # M x E

def inv(a):
    return pow(int(a), -1, p)

# Newton divided differences, vectorised over entries
dd = Y.copy()
coef = np.zeros((M, E), dtype=np.int64)
coef[0] = dd[0]
for j in range(1, M):
    den = np.array([inv((xs[i] - xs[i - j]) % p) for i in range(j, M)], dtype=np.int64)
    dd[j:] = (dd[j:] - dd[j - 1:-1]) % p * den[:, None] % p
    coef[j] = dd[j]
# Newton -> monomial coefficients (ascending), vectorised
poly = np.zeros((M, E), dtype=np.int64)
poly[0] = coef[M - 1]
deg = 0
for j in range(M - 2, -1, -1):
    # poly = poly * (x - x_j) + a_j
    shifted = np.zeros_like(poly)
    shifted[1:deg + 2] = poly[0:deg + 1]
    poly = (shifted - (xs[j] * poly) % p) % p
    poly[0] = (poly[0] + coef[j]) % p
    deg += 1

def pdeg(a):
    nz = np.flatnonzero(a)
    return int(nz[-1]) if nz.size else -1

def ptrim(a):
    d = pdeg(a)
    return a[:d + 1].copy() if d >= 0 else np.zeros(1, dtype=np.int64)

def pmul(a, b):
    a, b = ptrim(a), ptrim(b)
    out = np.zeros(len(a) + len(b) - 1, dtype=object)
    for i, ai in enumerate(a):
        if ai:
            out[i:i + len(b)] += int(ai) * b.astype(object)
    return (out % p).astype(np.int64)

def pdivmod(a, b):
    a, b = ptrim(a) % p, ptrim(b) % p
    db = pdeg(b)
    ib = inv(b[db])
    q = np.zeros(max(len(a) - db, 1), dtype=np.int64)
    r = a.copy()
    while pdeg(r) >= db:
        dr = pdeg(r)
        c = r[dr] * ib % p
        q[dr - db] = c
        r[dr - db:dr + 1] = (r[dr - db:dr + 1] - c * b) % p
    return q, ptrim(r)

def pgcd(a, b):
    a, b = ptrim(a), ptrim(b)
    while pdeg(b) >= 0 and np.any(b):
        a, b = b, pdivmod(a, b)[1]
    return a

def pmonic(a):
    a = ptrim(a)
    return a * inv(a[-1]) % p

# modulus polynomial prod (x - x_i)
Mpoly = np.array([1], dtype=np.int64)
for x in xs:
    Mpoly = pmul(Mpoly, np.array([(-int(x)) % p, 1], dtype=np.int64))

def ratrec(F, bound_num):
    """Extended Euclid on (Mpoly, F): returns (num, den) with deg num <= bound_num, deg den <= M - 2 - bound_num."""
    r0, r1 = Mpoly.copy(), ptrim(F)
    t0, t1 = np.array([0], dtype=np.int64), np.array([1], dtype=np.int64)
    while pdeg(r1) > bound_num:
        q, r = pdivmod(r0, r1)
        r0, r1 = r1, r
        t0, t1 = t1, ptrim((np.pad(t0, (0, max(0, len(pmul(q, t1)) - len(t0)))) - pmul(q, t1)) % p)
    return r1, t1

results = []
den_lcm = np.array([1], dtype=np.int64)
t1 = time.time()
for e in range(E):
    F = poly[:, e]
    best = None
    # try the balanced bound, then verify by the degree budget: deg num + deg den < M - margin
    for bound in (M // 2, M // 3, 2 * M // 3):
        num, den = ratrec(F, bound)
        dn, dd_ = pdeg(num), pdeg(den)
        if dn + dd_ <= M - 12:
            best = (num, den)
            break
    ri, c = ent[e]
    if best is None:
        results.append(dict(row=int(sel[ri]), pivot_weight=int(w[sel[ri]]), col=c, reconstructed=False))
        continue
    num, den = best
    g = pgcd(num, den)
    if pdeg(g) > 0:
        num, den = pdivmod(num, g)[0], pdivmod(den, g)[0]
    den = pmonic(den)
    # lcm update
    gg = pgcd(den_lcm, den)
    den_lcm = pmul(den_lcm, pdivmod(den, gg)[0])
    results.append(dict(row=int(sel[ri]), pivot_weight=int(w[sel[ri]]), col=c, col_weight=int(tg.weights[c]),
                        reconstructed=True, deg_num=pdeg(num), deg_den=pdeg(den)))
rec_seconds = time.time() - t1
D = pmonic(den_lcm)
degD = pdeg(D)

def peval(a, x):
    v = 0
    for cf in reversed(a.tolist()):
        v = (v * x + int(cf)) % p
    return v

# intersections of the line with the known exceptional lines
s_ordertwo = (-(1 + A0 + B0)) * inv((A1 + B1) % p) % p if (A1 + B1) % p else None
s_tangent = (-(3 + A0)) * inv(A1) % p
vals = dict(D_at_ordertwo=peval(D, s_ordertwo) if s_ordertwo is not None else None, D_at_tangent=peval(D, s_tangent))
summary = dict(k=k, prime=p, N=N, M=M, line=dict(A0=A0, B0=B0, A1=A1, B1=B1), dim=base['rank'], profile=base['profile'],
               bad_points=bad, stored_rows=R, entries=E,
               reconstructed=sum(1 for r in results if r['reconstructed']),
               max_deg_den=max((r['deg_den'] for r in results if r['reconstructed']), default=None),
               max_deg_num=max((r['deg_num'] for r in results if r['reconstructed']), default=None),
               deg_lcm_den=degD, D_values_at_known_lines=vals, s_ordertwo=s_ordertwo, s_tangent=s_tangent,
               closure_seconds=round(closure_seconds), reconstruction_seconds=round(rec_seconds),
               entries_detail=results)
json.dump(summary, open(args.out, 'w'), indent=1)
print(json.dumps({kk: v for kk, v in summary.items() if kk != 'entries_detail'}, indent=1))
