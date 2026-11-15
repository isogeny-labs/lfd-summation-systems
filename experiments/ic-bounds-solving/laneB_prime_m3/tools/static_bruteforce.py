#!/usr/bin/env python3
"""Ordinary-ring static Macaulay solving degree sd (and exact d_lf, d_ff via lfd_engine) at k=2 (6 variables).
sd: rows m*f (f in chains + {g}, deg(m f) <= D) processed degree by degree into a sparse echelon basis
(pivot = largest monomial in degree-then-lex order); 1 in Mac_D iff a pivot with leading monomial 1 exists.
No chain-quotient shortcut: all monomial multiples in the ordinary ring are generated and eliminated."""
import sys, os, time, json, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from build_m3 import *

def system(p, a, b, t, k):
    """Ordinary ring inputs, r = 3k vars u_{i,j} = var i*k+j. Returns (chains, g) as {exps: c}."""
    r = 3 * k
    chains = []
    for i in range(3):
        for j in range(k):
            sq = [0] * r; sq[i * k + j] = 2
            f = {tuple(sq): 1}
            if j + 1 < k:
                nx = [0] * r; nx[i * k + j + 1] = 1
                f[tuple(nx)] = p - 1
            else:
                f[tuple([0] * r)] = p - 1
            chains.append(f)
    g, delta = build_g(p, a, b, t, k)
    G = {}
    for idx in np.argwhere(g):
        e = [int(x) for x in idx]
        mask = e[0] | (e[1] << k) | (e[2] << 2 * k)
        ex = tuple((mask >> bit) & 1 for bit in range(r))
        G[ex] = int(g[tuple(e)])
    return chains, G, delta

def sd_sparse(chains, G, p, r, Dmax, verbose=False):
    gens = chains + [G]
    gdeg = [max(sum(e) for e in f) for f in gens]
    # monomial key: (deg, exps) tuple compared lexicographically
    def key(e): return (sum(e), e)
    pivots = {}
    def reduce_insert(row):
        while row:
            lead = max(row, key=key)
            c = row[lead]
            P = pivots.get(lead)
            if P is None:
                inv = pow(c, -1, p)
                pivots[lead] = {e: v * inv % p for e, v in row.items()}
                return True
            for e, v in P.items():
                nv = (row.get(e, 0) - c * v) % p
                if nv: row[e] = nv
                else: row.pop(e, None)
        return False
    def monos(r, d):
        for comb in itertools.combinations_with_replacement(range(r), d):
            e = [0] * r
            for i in comb: e[i] += 1
            yield tuple(e)
    zero = tuple([0] * r)
    prof = []
    for D in range(Dmax + 1):
        t0 = time.time()
        for f, d in zip(gens, gdeg):
            if D - d < 0: continue
            for m in monos(r, D - d):
                row = {}
                for e, c in f.items():
                    ee = tuple(x + y for x, y in zip(e, m))
                    row[ee] = (row.get(ee, 0) + c) % p
                row = {e: c for e, c in row.items() if c}
                reduce_insert(row)
        has1 = zero in pivots
        prof.append(dict(D=D, rank=len(pivots), one_in_Mac=has1, seconds=round(time.time() - t0, 2)))
        if verbose: print(prof[-1], flush=True)
        if has1:
            return D, prof
    return None, prof


def sd_sparse_int(chains, G, p, r, Dmax, verbose=False, log=None):
    """Same elimination with packed-integer monomial keys: key = deg << (4r) | 4-bit exponents (var 0 most significant).
    Needs every exponent <= 15 (true while Dmax <= 15). Order: degree, then lex."""
    S = 4 * r
    def pk(e):
        v = 0
        for i, x in enumerate(e):
            assert x <= 15
            v |= x << (4 * (r - 1 - i))
        return (sum(e) << S) | v
    gens = [{pk(e): c for e, c in f.items()} for f in chains + [G]]
    gdeg = [max(sum(e) for e in f) for f in chains + [G]]
    assert Dmax <= 15
    pivots = {}
    def reduce_insert(row):
        while row:
            lead = max(row)
            c = row[lead]
            P = pivots.get(lead)
            if P is None:
                inv = pow(c, -1, p)
                pivots[lead] = {e: v * inv % p for e, v in row.items()}
                return True
            for e, v in P.items():
                nv = (row.get(e, 0) - c * v) % p
                if nv: row[e] = nv
                else: row.pop(e, None)
        return False
    import itertools
    def monos(d):
        for comb in itertools.combinations_with_replacement(range(r), d):
            v = d << S
            for i in comb: v += 1 << (4 * (r - 1 - i))
            yield v
    prof = []
    for D in range(Dmax + 1):
        t0 = time.time()
        for f, d in zip(gens, gdeg):
            if D - d < 0: continue
            for m in monos(D - d):
                row = {}
                for e, c in f.items():
                    ee = e + m
                    row[ee] = (row.get(ee, 0) + c) % p
                reduce_insert({e: c for e, c in row.items() if c})
        has1 = 0 in pivots
        prof.append(dict(D=D, rank=len(pivots), one_in_Mac=has1, seconds=round(time.time() - t0, 2)))
        if verbose: print(prof[-1], flush=True)
        if log: print(prof[-1], file=log, flush=True)
        if has1:
            return D, prof
    return None, prof

def lfd_ff(chains, G, p, r, Dmax):
    sys.path.insert(0, '/home/user/ml-cryptanalysis/experiments/ic-curve-validation/code')
    import lfd_engine as L
    res = L.exact_lfd(chains + [G], r, p, [], Dmax)
    return dict(lfd=res['lfd'], falls=res['falls'], d_ff=(res['falls'][0] if res['falls'] else None),
                certified_at=res['certified_at'],
                per_degree=[{k: v for k, v in x.items()} for x in res['per_degree']])

if __name__ == '__main__':
    p, k, a, b, t = map(int, sys.argv[1:6])
    r = 3 * k
    chains, G, delta = system(p, a, b, t, k)
    t0 = time.time()
    sd, prof = (sd_sparse_int if len(sys.argv) > 6 else sd_sparse)(chains, G, p, r, 3 * k + delta, verbose=True)
    print('sd', sd, 'delta', delta, 'time', time.time() - t0)
