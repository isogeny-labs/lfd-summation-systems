#!/usr/bin/env python3
"""Independent verification of item 5 ("Fall-space structure", experiments/ic-bounds-fall).

Written from the definitions only.  It does NOT import csd/ic_bounds_fall.py, the producer's
scripts, or any repository closure engine (ic_fall_certificates, ic_chain_shift_closure, ...).
Reused: the dense mod-p row-echelon class `Ech` of the independent validator engine
experiments/ic-curve-validation/code/lfd_engine.py (no closure logic from it).

Model.  Modulo the chains u_{i,j}^2 = u_{i,j+1}, u_{i,k-1}^2 = 1, a squarefree chain monomial is
an exponent pair (a1, a2) in (Z/2^k)^2, its degree is popcount(a1) + popcount(a2), and the chain
variable u_{i,j} acts as translation by 2^j in block i (mod 2^k).  The generator is
g = S_3(x1, x2, t) with x_i = u_{i,0}, built here from (a, b, t) by the textbook formula
S_3 = (x1 - x2)^2 t^2 - 2((x1 + x2)(x1 x2 + a) + 2b) t + (x1 x2 - a)^2 - 4b(x1 + x2).
V_D = least subspace of M_<=D containing g and closed under every translation applied to
its elements of weight < D (breadth-first: round r+1 adds all translates of V^(r) cap M_<=D-1).

Columns are sorted by weight DESCENDING (ties: (a1, a2) ascending), so in the RREF the rows
with pivot weight <= d span V cap M_<=d.

Subcommands (each writes one JSON into --out):
  closures   V_4(k) for the primary instance, dims, profiles, rank history, closedness, 1 notin
  free       free extension from base k0 (7 and 6), equality with V_4(k), exclusion certificate
  levels     level-set Moebius counts of the fall space, fall-space structure
  w          symbolic kernel (sympy) and the instance check of the weight-2 element w
  cross      cross-instance dims, pivot sets, free-extension checks
  mult       multiplier filtration on the unit-ideal instance
  d5         degree-5 closure and the local-to-global certificate
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'experiments/ic-curve-validation/code'))
from lfd_engine import Ech  # noqa: E402  (dense mod-p echelon only)

PC = np.array([bin(i).count('1') for i in range(1 << 13)], dtype=np.int64)

INSTANCES = {
    'p40961': dict(p=40961, a=21223, b=9887, target=(4748, 524)),
    'p12289r3': dict(p=12289, a=3899, b=9710, target=(9895, 11610)),
    'p12289j1728': dict(p=12289, a=1, b=0, target=(9362, 1126)),
    'p1073153': dict(p=1073153, a=281783, b=132345, target=(247294, 569236)),
}


# ------------------------------------------------------------------------------------------
# generator, points
# ------------------------------------------------------------------------------------------

def s3_poly(a, b, t, p):
    """S_3(x1, x2, t) as {(e1, e2): c mod p}, expanded by hand from the textbook formula."""
    c = {}

    def add(e, v):
        c[e] = (c.get(e, 0) + v) % p
    # t^2 (x1 - x2)^2
    add((2, 0), t * t); add((0, 2), t * t); add((1, 1), -2 * t * t)
    # -2 t ((x1 + x2)(x1 x2 + a) + 2b) = -2t (x1^2 x2 + x1 x2^2 + a x1 + a x2 + 2b)
    add((2, 1), -2 * t); add((1, 2), -2 * t); add((1, 0), -2 * a * t); add((0, 1), -2 * a * t)
    add((0, 0), -4 * b * t)
    # (x1 x2 - a)^2
    add((2, 2), 1); add((1, 1), -2 * a); add((0, 0), a * a)
    # -4b (x1 + x2)
    add((1, 0), -4 * b); add((0, 1), -4 * b)
    return {e: v for e, v in c.items() if v}


def on_curve(inst):
    p, a, b = inst['p'], inst['a'], inst['b']
    x, y = inst['target']
    return (y * y - (x ** 3 + a * x + b)) % p == 0


def subgroup(p, k):
    N = 1 << k
    assert (p - 1) % N == 0
    for c in range(2, p):
        z = pow(c, (p - 1) // N, p)
        if pow(z, N // 2, p) != 1:
            return z
    raise ValueError


def solutions(g, p, k):
    """ordered pairs (x1, x2) in H^2 with g(x1, x2) = 0, H = mu_{2^k}."""
    z = subgroup(p, k)
    N = 1 << k
    H = [pow(z, e, p) for e in range(N)]
    sols = []
    for e1, x1 in enumerate(H):
        for e2, x2 in enumerate(H):
            v = 0
            for (i, j), c in g.items():
                v = (v + c * pow(x1, i, p) * pow(x2, j, p)) % p
            if v == 0:
                sols.append((e1, e2))
    return sols


# ------------------------------------------------------------------------------------------
# group-algebra grid
# ------------------------------------------------------------------------------------------

class Grid:
    def __init__(self, k, D):
        self.k, self.D, self.N = k, D, 1 << k
        N = self.N
        a1, a2 = np.meshgrid(np.arange(N), np.arange(N), indexing='ij')
        a1, a2 = a1.ravel(), a2.ravel()
        w = PC[a1] + PC[a2]
        keep = w <= D
        a1, a2, w = a1[keep], a2[keep], w[keep]
        order = np.lexsort((a2, a1, -w))
        self.a1, self.a2, self.wt = a1[order], a2[order], w[order]
        self.n = len(self.wt)
        self.idx = np.full(N * N, -1, dtype=np.int64)
        self.idx[self.a1 * N + self.a2] = np.arange(self.n)
        self.low_start = int(np.argmax(self.wt <= D - 1))
        self.lev = (self.a1 | self.a2)  # level set as bitmask
        self.shifts = []
        for i in range(2):
            for j in range(k):
                b1 = self.a1[self.low_start:].copy()
                b2 = self.a2[self.low_start:].copy()
                if i == 0:
                    b1 = (b1 + (1 << j)) % N
                else:
                    b2 = (b2 + (1 << j)) % N
                t = self.idx[b1 * N + b2]
                assert np.all(t >= 0)
                self.shifts.append(t)

    def col(self, e1, e2):
        return int(self.idx[(e1 % self.N) * self.N + (e2 % self.N)])

    def translate(self, poly, h1, h2, p):
        row = np.zeros(self.n, dtype=np.float64)
        for (e1, e2), c in poly.items():
            j = self.col(e1 + h1, e2 + h2)
            assert j >= 0, 'translate leaves M_<=D'
            row[j] = (row[j] + c) % p
        return row

    def shift_rows(self, rows, s):
        out = np.zeros((rows.shape[0], self.n), dtype=rows.dtype)
        out[:, self.shifts[s]] = rows[:, self.low_start:]
        return out

    def one(self):
        r = np.zeros(self.n)
        r[self.col(0, 0)] = 1
        return r


def masks_of_weight(k, wmax):
    out = []
    for w in range(wmax + 1):
        for bits in itertools.combinations(range(2 * k), w):
            h1 = sum(1 << b for b in bits if b < k)
            h2 = sum(1 << (b - k) for b in bits if b >= k)
            out.append((h1, h2, w))
    return out


def ech_of(rows, n, p, cap=None):
    e = Ech(n, p, cap=cap or max(1, min(n, rows.shape[0] + 8)))
    if rows.shape[0]:
        e.add(np.mod(rows, p))
    return e


def rank_of(rows, p):
    if rows.shape[0] == 0 or rows.shape[1] == 0:
        return 0
    return ech_of(rows, rows.shape[1], p, cap=min(rows.shape) + 2).rk


def reduces_to_zero(E, rows, chunk=512):
    for a in range(0, rows.shape[0], chunk):
        if np.any(E.reduce(rows[a:a + chunk])):
            return False
    return True


def sha(arr):
    return hashlib.sha256(np.ascontiguousarray(np.asarray(arr, dtype=np.int64)).tobytes()).hexdigest()


# ------------------------------------------------------------------------------------------
# closure (BFS, own)
# ------------------------------------------------------------------------------------------

def closure(g, k, D, p, cap=None, done_cap=None, verify=True, log=print):
    G = Grid(k, D)
    t0 = time.time()
    V = Ech(G.n, p, cap=cap or G.n)
    gens = [G.translate(g, h1, h2, p) for h1, h2, _ in masks_of_weight(k, D - 2)]
    V.add(np.array(gens))
    low_count = lambda: int(np.sum(G.wt[V.piv] <= D - 1))
    hist, lowhist = [V.rk], [low_count()]
    done = Ech(G.n, p, cap=done_cap or (G.n - G.low_start))
    while True:
        low = V.E[G.wt[V.piv] <= D - 1]
        resid = done.reduce(low)
        resid = resid[np.any(resid != 0, axis=1)]
        if resid.shape[0] == 0:
            break
        new = done.add(resid)
        newrows = np.vstack(new)
        for s in range(2 * k):
            V.add(G.shift_rows(newrows, s))
        hist.append(V.rk)
        lowhist.append(low_count())
        log(f'  k={k} D={D} round {len(hist) - 1}: rank {V.rk}, low {lowhist[-1]}, {time.time() - t0:.1f}s')
    rounds = len(hist) - 1 - (1 if hist[-1] == hist[-2] and lowhist[-1] == lowhist[-2] else 0)
    closed = None
    if verify:
        low = V.E[G.wt[V.piv] <= D - 1]
        closed = all(reduces_to_zero(V, G.shift_rows(low, s)) for s in range(2 * k))
    piv = V.piv.copy()
    E = V.E.copy()
    pw = G.wt[piv]
    order = np.argsort(piv)
    return dict(G=G, V=V, E=E[order], piv=piv[order], pw=pw[order], dim=V.rk,
                profile=[int(np.sum(pw <= d)) for d in range(D + 1)],
                one_in=not np.any(V.reduce(G.one().reshape(1, -1))),
                closed=closed, rank_history=hist, low_history=lowhist,
                low_birth=[lowhist[0]] + [lowhist[i] - lowhist[i - 1] for i in range(1, len(lowhist))],
                n_generator_rows=len(gens), rounds=rounds, seconds=round(time.time() - t0, 1))


# ------------------------------------------------------------------------------------------
# free extension
# ------------------------------------------------------------------------------------------

def free_extension(base, k, k0, D=4):
    """rows tau_h(v), v a base RREF row of weight <= D - wt h, h on levels k0..k-1 (both blocks)."""
    G0 = base['G']
    G = Grid(k, D)
    rows = []
    bits = [(0, j) for j in range(k0, k)] + [(1, j) for j in range(k0, k)]
    for w in range(D + 1):
        for sub in itertools.combinations(bits, w):
            h1 = sum(1 << j for i, j in sub if i == 0)
            h2 = sum(1 << j for i, j in sub if i == 1)
            sel = np.flatnonzero(base['pw'] <= D - w)
            if sel.size == 0:
                continue
            cols = np.flatnonzero(G0.wt <= D - w)
            B = base['E'][sel][:, cols]
            assert not np.any(base['E'][sel][:, G0.wt > D - w])
            R = np.zeros((sel.size, G.n))
            tgt = G.idx[(G0.a1[cols] + h1) * G.N + (G0.a2[cols] + h2)]
            assert np.all(tgt >= 0)
            R[:, tgt] = B
            rows.append(R)
    return G, np.vstack(rows)


def embed(base, G):
    """read a basis of grid G0 (exponents < 2^k0) inside grid G."""
    G0 = base['G']
    R = np.zeros((base['E'].shape[0], G.n))
    tgt = G.idx[G0.a1 * G.N + G0.a2]
    R[:, tgt] = base['E']
    return R


def certificate(G, W, g, p, D=4):
    E = ech_of(W, G.n, p)
    gens = np.array([G.translate(g, h1, h2, p) for h1, h2, _ in masks_of_weight(G.k, D - 2)])
    contains_gen = reduces_to_zero(E, gens)
    pw = G.wt[E.piv]
    low = E.E[pw <= D - 1]
    closed = all(reduces_to_zero(E, G.shift_rows(low, s)) for s in range(2 * G.k))
    one = not np.any(E.reduce(G.one().reshape(1, -1)))
    return dict(dim=E.rk, contains_generator_rows=contains_gen, closed=closed, contains_one=one,
                profile=[int(np.sum(pw <= d)) for d in range(D + 1)])


def compare(A, B, p):
    ra, rb = rank_of(A, p), rank_of(B, p)
    rs = rank_of(np.vstack([A, B]), p)
    return dict(rank_a=ra, rank_b=rb, rank_sum=rs, equal=(ra == rb == rs), only_in_a=rs - rb, only_in_b=rs - ra)


def intersect_coord(E, inside, p):
    """basis of rowspace(E) cap span(columns 'inside'), via elimination with outside columns first."""
    perm = np.concatenate([np.flatnonzero(~inside), np.flatnonzero(inside)])
    nout = int(np.sum(~inside))
    M = E[:, perm]
    e = ech_of(M, M.shape[1], p)
    sel = e.piv >= nout
    out = np.zeros((int(sel.sum()), E.shape[1]))
    out[:, perm] = e.E[sel]
    return out


# ------------------------------------------------------------------------------------------
# subcommands
# ------------------------------------------------------------------------------------------

def inst_g(name):
    I = INSTANCES[name]
    return I, s3_poly(I['a'], I['b'], I['target'][0], I['p'])


def cmd_closures(args, out, log):
    I, g = inst_g('p40961')
    p = I['p']
    res = dict(instance=I, on_curve=on_curve(I), generator={f'{a},{b}': c for (a, b), c in sorted(g.items())})
    # cross-check the generator with the producer's instance.json (data comparison only)
    prod = json.loads((ROOT / 'experiments/ic-bounds-fall/p40961-s7/instance.json').read_text())
    res['generator_equals_producer'] = {(e1, e2): c for e1, e2, c in prod['generator_pairs']} == g
    per = {}
    for k in args.k:
        r = closure(g, k, 4, p, cap=4 * k * k + 172 * k + 200, done_cap=600, log=log)
        G = r['G']
        n_sol = len(solutions(g, p, k)) if k <= 10 else None
        per[k] = dict(monomials=G.n, dim_V4=r['dim'], profile=r['profile'], formula=4 * k * k + 172 * k - 571,
                      formula_low=4 * k + 87, one_in_V4=r['one_in'], closed=r['closed'],
                      rank_history=r['rank_history'], low_history=r['low_history'], low_birth=r['low_birth'],
                      generator_rows=r['n_generator_rows'], rounds=r['rounds'], seconds=r['seconds'],
                      solutions=n_sol, basis_sha256_own_order=sha(r['E']))
        log(f'k={k}: {per[k]}')
    res['per_k'] = per
    (out / 'closures.json').write_text(json.dumps(res, indent=1))


def cmd_free(args, out, log):
    I, g = inst_g(args.instance)
    p = I['p']
    res = {}
    cache = {}

    def V(k):
        if k not in cache:
            cache[k] = closure(g, k, 4, p, cap=4 * k * k + 172 * k + 200, done_cap=600, log=lambda m: None)
        return cache[k]
    for k0 in args.k0:
        base = V(k0)
        res[f'base{k0}'] = dict(dim=base['dim'], profile=base['profile'], one_in=base['one_in'])
        for k in args.k:
            if k <= k0:
                continue
            G, W = free_extension(base, k, k0)
            Vk = V(k)
            cert = certificate(G, W, g, p)
            cmp_ = compare(W, Vk['E'], p)
            inside = (G.a1 < (1 << k0)) & (G.a2 < (1 << k0))
            stable = intersect_coord(Vk['E'], inside, p)
            S_emb = embed(base, G)
            st = compare(stable, S_emb, p)
            # stable low part read back in the base grid: contained in the base?
            G0 = base['G']
            back = np.zeros((stable.shape[0], G0.n))
            src = G0.idx[G.a1[inside] * G0.N + G.a2[inside]]
            back[:, src] = stable[:, inside]
            in_base = reduces_to_zero(base['V'], back)
            rec = dict(rows=int(W.shape[0]), predicted=4 * k * k + 172 * k - 571, certificate=cert,
                       compare_W_vs_V4=cmp_, stable_low_dim=int(stable.shape[0]),
                       stable_vs_base=st, stable_low_in_base=in_base)
            res[f'k0={k0},k={k}'] = rec
            log(f'k0={k0} k={k}: {rec}')
    # top-level ranks of the fall spaces at k = 6, 7
    for k in (6, 7):
        if k in args.top:
            r = V(k)
            G = r['G']
            F = r['E'][r['pw'] <= 3]
            top = ((G.a1 >> (k - 1)) & 1) | ((G.a2 >> (k - 1)) & 1)
            res[f'fall_rank_on_top_level_k{k}'] = rank_of(F[:, top.astype(bool)], p)
    (out / f'free-{args.instance}.json').write_text(json.dumps(res, indent=1))
    log(json.dumps({k: v for k, v in res.items() if k.startswith('fall')}))


def level_analysis(F, G, p):
    k = G.k
    supp = np.flatnonzero(np.any(F != 0, axis=0))
    Fs = F[:, supp]
    lev = G.lev[supp]
    dimF = rank_of(Fs, p)
    dims, bases = {}, {}
    for L in range(1 << k):
        inside = (lev & ~L) == 0
        B = intersect_coord(Fs, inside, p) if inside.any() else np.zeros((0, Fs.shape[1]))
        dims[L], bases[L] = B.shape[0], B
    mob = {}
    for L in range(1 << k):
        s = 0
        sub = L
        while True:
            s += (-1) ** (PC[L] - PC[sub]) * dims[sub]
            if sub == 0:
                break
            sub = (sub - 1) & L
        mob[L] = int(s)
    graded = {}
    for L in range(1 << k):
        if dims[L] == 0:
            graded[L] = 0
            continue
        parts = [bases[L & ~(1 << i)] for i in range(k) if L >> i & 1]
        parts = [b for b in parts if b.shape[0]]
        lower = rank_of(np.vstack(parts), p) if parts else 0
        graded[L] = dims[L] - lower
    nz = {L: c for L, c in mob.items() if c}
    span = rank_of(np.vstack([bases[L] for L in nz if bases[L].shape[0]]), p)
    name = lambda L: ','.join(str(i) for i in range(k) if L >> i & 1)
    return dict(dim=dimF, moebius={name(L): c for L, c in sorted(nz.items())},
                graded={name(L): c for L, c in sorted(graded.items()) if c},
                sum_moebius=int(sum(nz.values())), sum_graded=int(sum(graded.values())),
                negative=int(sum(1 for c in nz.values() if c < 0)), span_of_pieces=span)


def cmd_levels(args, out, log):
    I, g = inst_g('p40961')
    p = I['p']
    res = {}
    cache = {k: closure(g, k, 4, p, cap=4 * k * k + 172 * k + 200, done_cap=600, log=lambda m: None)
             for k in sorted(set(args.k) | set(args.kstruct) | {7})}
    prod = json.loads((ROOT / 'experiments/ic-bounds-fall/p40961-s7/analysis.json').read_text())
    for k in args.k:
        r = cache[k]
        F = r['E'][r['pw'] <= 3]
        la = level_analysis(F, r['G'], p)
        pc = prod['level_decomposition'].get(str(k), {}).get('fall_space', {}).get('exact_level_counts')
        la['equals_producer_counts'] = (pc == la['moebius']) if pc is not None else None
        res[f'levels_k{k}'] = la
        log(f'levels k={k}: {la}')
    base = cache[7]
    G7 = base['G']
    F7 = base['E'][base['pw'] <= 3]
    two = base['E'][base['pw'] <= 2]   # span(g, w)
    for k in args.kstruct:
        r = cache[k]
        G = r['G']
        Fk = r['E'][r['pw'] <= 3]
        rows = [embed(dict(G=G7, E=F7), G)]
        for j in range(7, k):
            for i in range(2):
                h1, h2 = ((1 << j), 0) if i == 0 else (0, (1 << j))
                cols = np.flatnonzero(G7.wt <= 2)
                assert not np.any(two[:, G7.wt > 2])
                tgt = G.idx[(G7.a1[cols] + h1) * G.N + (G7.a2[cols] + h2)]
                assert np.all(tgt >= 0)
                R = np.zeros((two.shape[0], G.n))
                R[:, tgt] = two[:, cols]
                rows.append(R)
        res[f'fall_structure_k{k}'] = compare(np.vstack(rows), Fk, p)
        log(f'fall structure k={k}: {res[f"fall_structure_k{k}"]}')
    (out / 'levels.json').write_text(json.dumps(res, indent=1))


def cmd_w(args, out, log):
    import sympy as sp
    a, b, t, x1, x2 = sp.symbols('a b t x1 x2')
    S3 = sp.expand((x1 - x2) ** 2 * t ** 2 - 2 * ((x1 + x2) * (x1 * x2 + a) + 2 * b) * t + (x1 * x2 - a) ** 2 - 4 * b * (x1 + x2))
    qs = {(i, j): sp.Symbol(f'q{i}{j}') for i in range(3) for j in range(3)}
    Q = sum(qs[h] * x1 ** h[0] * x2 ** h[1] for h in qs)
    P = sp.Poly(sp.expand(Q * S3), x1, x2)
    killed = [(e1, e2) for e1 in range(5) for e2 in range(5) if PC[e1] + PC[e2] >= 3]
    eqs = [P.coeff_monomial(x1 ** e1 * x2 ** e2) for e1, e2 in killed]
    unk = [qs[h] for h in sorted(qs)]
    A = sp.Matrix([[sp.diff(e, u) for u in unk] for e in eqs])
    ns = A.nullspace()
    ns = [sp.simplify(v / v[-1]) if v[-1] != 0 else v for v in ns]
    Qw = {(0, 0): 0, (2, 2): 1, (1, 2): 2 * t, (2, 1): 2 * t, (1, 1): 2 * (a + 5 * t ** 2),
          (2, 0): -t * (11 * a * t + 8 * b + 17 * t ** 3) / (a + 3 * t ** 2),
          (0, 2): -t * (11 * a * t + 8 * b + 17 * t ** 3) / (a + 3 * t ** 2),
          (1, 0): 2 * (5 * a ** 2 * t + 2 * a * b + 15 * a * t ** 3 - 2 * b * t ** 2 + 16 * t ** 5) / (a + 3 * t ** 2),
          (0, 1): 2 * (5 * a ** 2 * t + 2 * a * b + 15 * a * t ** 3 - 2 * b * t ** 2 + 16 * t ** 5) / (a + 3 * t ** 2)}
    vec = sp.Matrix([Qw[h] for h in sorted(qs)])
    in_kernel = all(sp.simplify(x) == 0 for x in (A * vec))
    delta0 = sp.Matrix([1 if h == (0, 0) else 0 for h in sorted(qs)])
    rank_span = sp.Matrix.hstack(*(ns + [vec, delta0])).rank(simplify=True)
    res = dict(killed=killed, n_equations=len(eqs), kernel_dim=len(ns), formula_in_kernel=in_kernel,
               rank_of_kernel_plus_formula_plus_delta0=int(rank_span))
    # instance check
    I, g = inst_g('p40961')
    p = I['p']
    av, bv, tv = I['a'], I['b'], I['target'][0]
    qnum = {}
    for h, e in Qw.items():
        num, den = sp.fraction(sp.together(e))
        nv = int(num.subs({a: av, b: bv, t: tv})) % p
        dv = int(den.subs({a: av, b: bv, t: tv})) % p
        qnum[h] = nv * pow(dv, -1, p) % p
    res['a_plus_3t2_mod_p'] = (av + 3 * tv * tv) % p
    for k in args.k:
        r = closure(g, k, 4, p, cap=4 * k * k + 172 * k + 200, done_cap=600, log=lambda m: None)
        G = r['G']
        wrow = np.zeros(G.n)
        for h, c in qnum.items():
            if c:
                wrow = np.mod(wrow + c * G.translate(g, h[0], h[1], p), p)
        nz = np.flatnonzero(wrow)
        two = r['E'][r['pw'] <= 2]
        grow = G.translate(g, 0, 0, p)
        res[f'k{k}'] = dict(weight=int(G.wt[nz].max()), support=int(nz.size),
                            in_V4=not np.any(r['V'].reduce(wrow.reshape(1, -1))),
                            g_w_span_V4_cap_M2=compare(np.vstack([grow, wrow]), two, p))
    (out / 'w.json').write_text(json.dumps(res, indent=1, default=str))
    log(json.dumps(res, default=str))


def cmd_cross(args, out, log):
    res = {}
    bases = {}
    for name in args.instances:
        I, g = inst_g(name)
        p = I['p']
        rec = dict(instance=I, on_curve=on_curve(I), dims={}, profiles={}, solutions={})
        cache = {}
        for k in args.k:
            r = closure(g, k, 4, p, cap=4 * k * k + 172 * k + 200, done_cap=600, log=lambda m: None)
            cache[k] = r
            rec['dims'][k] = r['dim']
            rec['profiles'][k] = r['profile']
            rec['solutions'][k] = len(solutions(g, p, k)) if k <= 8 else None
            rec.setdefault('one_in', {})[k] = r['one_in']
            rec.setdefault('closed', {})[k] = r['closed']
        if 7 in cache:
            bases[name] = set(cache[7]['piv'].tolist())
            for k in args.kfree:
                if k in cache:
                    G, W = free_extension(cache[7], k, 7)
                    rec.setdefault('free', {})[k] = dict(cert=certificate(G, W, g, p), cmp=compare(W, cache[k]['E'], p))
        res[name] = rec
        log(f'{name}: {json.dumps(rec, default=str)[:600]}')
    names = list(bases)
    res['pivot_set_symdiff_vs_first'] = {n: len(bases[n] ^ bases[names[0]]) for n in names}
    (out / 'cross.json').write_text(json.dumps(res, indent=1, default=str))


def dft_matrix(z, N, p):
    e = np.outer(np.arange(N), np.arange(N)) % N
    pw = np.array([pow(z, i, p) for i in range(N)], dtype=object)
    return pw[e]


def matmul_mod(A, B, p):
    """exact (object-free) modular product via splitting into 16-bit halves."""
    A = np.asarray(A, dtype=np.int64) % p
    B = np.asarray(B, dtype=np.int64) % p
    lo, hi = B & 0xFFFF, B >> 16
    # A entries < 2^21, lo < 2^16: products < 2^37, sums over <= 2^12 terms < 2^49 (float-exact)
    r1 = np.rint(A.astype(np.float64) @ lo.astype(np.float64)).astype(np.int64) % p
    r2 = np.rint(A.astype(np.float64) @ hi.astype(np.float64)).astype(np.int64) % p
    return (r1 + (r2 * 65536) % p) % p


def cmd_mult(args, out, log):
    I, g = inst_g('p1073153')
    p = I['p']
    res = dict(instance=I, on_curve=on_curve(I))
    for k in args.k:
        N = 1 << k
        z = subgroup(p, k)
        Fm = dft_matrix(z, N, p).astype(np.int64)
        zi = pow(z, -1, p)
        Fi = dft_matrix(zi, N, p).astype(np.int64)
        Ninv = pow(N * N, -1, p)
        r = closure(g, k, 4, p, cap=4 * k * k + 172 * k + 200, done_cap=600, log=lambda m: None)
        G = r['G']

        def to_full(rows, Gx):
            X = np.zeros((rows.shape[0], N, N), dtype=np.int64)
            X[:, Gx.a1, Gx.a2] = rows.astype(np.int64)
            return X

        def fwd(X):   # X[r, e1, e2] -> hat[r, x1, x2] = sum X z^(e1 x1 + e2 x2)
            Y = matmul_mod(X.reshape(-1, N), Fm, p).reshape(X.shape)          # over e2
            Y = np.transpose(Y, (0, 2, 1))
            Y = matmul_mod(Y.reshape(-1, N), Fm, p).reshape(X.shape)          # over e1
            return np.transpose(Y, (0, 2, 1))

        def inv(Y):
            X = matmul_mod(Y.reshape(-1, N), Fi, p).reshape(Y.shape)
            X = np.transpose(X, (0, 2, 1))
            X = matmul_mod(X.reshape(-1, N), Fi, p).reshape(Y.shape)
            return np.transpose(X, (0, 2, 1)) * Ninv % p

        gfull = to_full(G.translate(g, 0, 0, p).reshape(1, -1), G)
        ghat = fwd(gfull)[0]
        nzero = int(np.sum(ghat == 0))
        ginv_hat = np.vectorize(lambda v: pow(int(v), -1, p))(ghat).astype(np.int64)
        wfull = (PC[np.arange(N)][:, None] + PC[np.arange(N)][None, :]).ravel()
        rec = dict(dim_V4=r['dim'], g_zero_evaluations=nzero, solutions=len(solutions(g, p, k)))

        def filt(rows):
            Yh = fwd(to_full(rows, G))
            Q = inv(Yh * ginv_hat[None] % p).reshape(rows.shape[0], -1)
            # sanity: Q * g == rows (check through the transform)
            dims = []
            for w in range(2 * k + 1):
                outside = wfull > w
                dims.append(rows.shape[0] - rank_of(Q[:, outside].astype(np.float64), p))
            return dims, Q
        Vrows = r['E']
        rec['dim_Q_cap_Bw'], Q = filt(Vrows)
        # exactness check: Q g recovers V on a sample of rows
        smp = Q[:5].reshape(5, N, N)
        back = inv(fwd(smp) * ghat[None] % p).reshape(5, -1)
        rows5 = to_full(Vrows[:5], G).reshape(5, -1)
        rec['multiplier_roundtrip_ok'] = bool(np.array_equal(back % p, rows5 % p))
        rec['max_multiplier_weight'] = int(next(w for w, d in enumerate(rec['dim_Q_cap_Bw']) if d == r['dim']))
        F = r['E'][r['pw'] <= 3]
        rec['dim_Qfall_cap_Bw'], _ = filt(F)
        # kappa_d(w) = dim(B_w cap M_<=d g^{-1}),  M_<=d g^{-1} = span of Dirac masses of weight <= d times g^{-1}
        for d in (4, 3):
            Gd = Grid(k, d)
            I_rows = np.zeros((Gd.n, N, N), dtype=np.int64)
            I_rows[np.arange(Gd.n), Gd.a1, Gd.a2] = 1
            ginv_full = inv(ginv_hat[None]).reshape(N, N)   # g^{-1} itself
            # (delta_m * g^{-1}) = translate of g^{-1} by m
            M = np.empty((Gd.n, N * N), dtype=np.int64)
            for i in range(Gd.n):
                M[i] = np.roll(np.roll(ginv_full, Gd.a1[i], axis=0), Gd.a2[i], axis=1).ravel()
            kap = []
            for w in range(min(2 * k, args.kappa_wmax[k == 7]) + 1):
                outside = wfull > w
                kap.append(Gd.n - rank_of(M[:, outside].astype(np.float64), p))
            rec[f'kappa{d}'] = kap
        res[k] = rec
        log(f'k={k}: {rec}')
    (out / 'mult.json').write_text(json.dumps(res, indent=1, default=str))


def cmd_d5(args, out, log):
    I, g = inst_g('p40961')
    p = I['p']
    res = {}
    for k in args.k:
        sols = solutions(g, p, k)
        r = closure(g, k, 5, p, cap=None, log=log)
        G = r['G']
        # evaluation of M_<=d on Z: monomial (e1, e2) at (x1, x2) = z^(s1 e1 + s2 e2)
        z = subgroup(p, k)
        ev = np.array([[pow(z, (s1 * e1 + s2 * e2) % (1 << k), p) for e1, e2 in zip(G.a1.tolist(), G.a2.tolist())] for s1, s2 in sols],
                      dtype=np.float64)
        rk = {d: rank_of(ev[:, G.wt <= d], p) for d in range(6)}
        nm = {d: int(np.sum(G.wt <= d)) for d in range(6)}
        res[k] = dict(solutions=len(sols), monomials=nm[5], dim_V5=r['dim'], target=nm[5] - len(sols),
                      eval_rank_below=rk[4], onto_below=(rk[4] == len(sols)),
                      ideal_profile=[nm[d] - rk[d] for d in range(6)], profile=r['profile'],
                      certificate=(r['dim'] == nm[5] - rk[5] and rk[4] == len(sols)),
                      closed=r['closed'], rank_history=r['rank_history'], rounds=r['rounds'], seconds=r['seconds'])
        log(f'k={k}: {res[k]}')
        (out / f'd5-k{k}.json').write_text(json.dumps(res[k], indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd')
    ap.add_argument('--k', type=int, nargs='*', default=[6, 7, 8])
    ap.add_argument('--k0', type=int, nargs='*', default=[7])
    ap.add_argument('--kstruct', type=int, nargs='*', default=[])
    ap.add_argument('--kfree', type=int, nargs='*', default=[8, 9])
    ap.add_argument('--top', type=int, nargs='*', default=[])
    ap.add_argument('--instance', default='p40961')
    ap.add_argument('--instances', nargs='*', default=['p40961', 'p12289r3', 'p12289j1728'])
    ap.add_argument('--kappa-wmax', type=int, nargs=2, default=[12, 6])
    ap.add_argument('--out', type=Path, default=ROOT / 'experiments/ic-bounds-fall/verify')
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    def log(m):
        print(m, flush=True)
    dict(closures=cmd_closures, free=cmd_free, levels=cmd_levels, w=cmd_w, cross=cmd_cross,
         mult=cmd_mult, d5=cmd_d5)[args.cmd](args, args.out, log)


if __name__ == '__main__':
    main()
