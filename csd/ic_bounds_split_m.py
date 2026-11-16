"""Split presentation with m summands and m - 2 free auxiliary abscissae: exact closures.

System (the generalization of the split encoding of csd.ic_bounds_m3 to m >= 3):
  variables  u_{i,j} (i < m, j < k) with the chain equations of the subgroup H of order 2^k,
             x_i = u_{i,0}, and a = m - 2 ordinary variables z_1 .. z_a (no chain, never reduced);
  inputs     g_1     = NF(S_3(x_1, x_2, z_1)),
             g_j     = NF(S_3(z_{j-1}, x_{j+1}, z_j))   (2 <= j <= m - 2),
             g_{m-1} = NF(S_3(z_{m-2}, x_m, t)).
For m = 3 this is exactly csd.ic_bounds_m3.split_generators (checked by scripts/run_ic_split_m.py).

V_D is computed modulo the chains by csd.ic_bounds_m3_batched.BatchedClosure, which is generic in
the monomial space; SplitSpace below supplies that space for several auxiliaries.  A closure that
does not contain 1 is accepted only after the engine's verification pass (every variable shift of
every basis row of degree < D reduces to zero, and every input lies in the span).

Exact value for a unit ideal: the least D with 1 in V_D, every smaller closure verified closed
without 1.
"""
from __future__ import annotations

import hashlib
import json
import math
import time
from itertools import combinations

import numpy as np

from .ic_bounds_m3 import chain_nf, padd, pmul, s3_quadratic
from .ic_bounds_m3_batched import TWO53, BatchedClosure


def compositions(total, parts):
    """Exponent tuples of length `parts` with sum exactly `total`."""
    if parts == 0:
        if total == 0:
            yield ()
        return
    for first in range(total + 1):
        for rest in compositions(total - first, parts - 1):
            yield (first,) + rest


class SplitSpace:
    """Monomials (mask, es) = u^mask z^es with weight(mask) + sum(es) <= D, degree descending,
    the constant monomial last.  Same interface as csd.ic_bounds_m3.MixedSpace."""

    def __init__(self, m, k, D, naux):
        self.m, self.k, self.D, self.naux = m, k, D, naux
        r = m * k
        masks_by_weight = [[sum(1 << j for j in c) for c in combinations(range(r), w)]
                           for w in range(min(r, D) + 1)]
        mons = []
        for d in range(D + 1):
            for e in range(d + 1):
                w = d - e
                if w > r:
                    continue
                for es in compositions(e, naux):
                    mons.extend((mask, es) for mask in masks_by_weight[w])
        mons.sort(key=lambda me: (me[0].bit_count() + sum(me[1]), me[1], me[0]), reverse=True)
        assert mons[-1] == (0, (0,) * naux)
        self.mons = mons
        self.index = {me: i for i, me in enumerate(mons)}
        self.deg = np.array([mask.bit_count() + sum(es) for mask, es in mons], dtype=np.int64)
        N = 1 << k
        cols = np.full((r + naux, len(mons)), -1, dtype=np.int64)
        for s, (mask, es) in enumerate(mons):
            if mask.bit_count() + sum(es) >= D:
                continue
            for j in range(r):
                blk, bit = divmod(j, k)
                sh = blk * k
                ex = ((mask >> sh) + (1 << bit)) % N
                cols[j, s] = self.index[((mask & ~((N - 1) << sh)) | (ex << sh), es)]
            for i in range(naux):
                cols[r + i, s] = self.index[(mask, es[:i] + (es[i] + 1,) + es[i + 1:])]
        self.cols = cols

    def __len__(self):
        return len(self.mons)

    def row(self, poly, p):
        v = np.zeros(len(self.mons), dtype=np.int64)
        for key, c in poly.items():
            i = self.index[key]
            v[i] = (v[i] + c) % p
        return v


def split_generators_m(E, t, k, m):
    """The m - 1 inputs as {(mask, (e_1, .., e_a)): c}; ordinary variables (x_1..x_m, z_1..z_a)."""
    assert m >= 3
    p, a = E.p, m - 2
    nv = m + a
    var = lambda i: {tuple(1 if j == i else 0 for j in range(nv)): 1}
    const = lambda c: {(0,) * nv: c % p}
    x = [var(i) for i in range(m)]
    z = [var(m + i) for i in range(a)]
    quad = lambda A, B, C, w: padd(p, (1, pmul(p, A, pmul(p, w, w))), (1, pmul(p, B, w)), (1, C))
    gens = [quad(*s3_quadratic(E, x[0], x[1], nv), z[0])]
    for j in range(1, a):
        gens.append(quad(*s3_quadratic(E, z[j - 1], x[j + 1], nv), z[j]))
    gens.append(quad(*s3_quadratic(E, z[a - 1], x[m - 1], nv), const(t)))
    out = []
    for g in gens:
        nf = chain_nf(p, g, k, m)
        out.append({(key[0], tuple(key[1:])): c for key, c in nf.items()})
    return out, gens


def generator_degree(g):
    return max(mask.bit_count() + sum(es) for mask, es in g)


def system_hash(p, m, k, gens):
    blob = json.dumps([p, m, k, [sorted((mask, list(es), c) for (mask, es), c in g.items()) for g in gens]],
                      separators=(',', ':'))
    return hashlib.sha256(blob.encode()).hexdigest()


def closure(gens, m, k, p, D, *, batch=768, log=None):
    """Exact V_D; returns a record.  Raises ValueError if float64 exactness cannot be guaranteed."""
    sp = SplitSpace(m, k, D, m - 2)
    n = len(sp)
    if n * (p - 1) ** 2 >= TWO53:
        raise ValueError(f'n (p-1)^2 = {n * (p - 1) ** 2} is not below 2^53')
    start = time.time()
    cl = BatchedClosure(sp, p, batch)
    rows = [sp.row(g, p) for g in gens if generator_degree(g) <= D]
    cl.add_inputs(rows)
    batches = cl.saturate(log)
    verified = None
    if not cl.contains_one:
        verified = bool(cl.verify_closed()) and all(cl.contains(r) for r in rows)
        if not verified:
            raise AssertionError('closure verification failed')
    return dict(D=D, monomials_le_D=n, inputs_used=len(rows), dim_V_D=int(cl.rank),
                low_profile=[cl.low_dimension(d) for d in range(D + 1)],
                contains_one=bool(cl.contains_one), verified_closed=verified,
                batches=int(batches), seconds=round(time.time() - start, 3))


def exact_unit_lfd(gens, m, k, p, *, d_min=3, d_max=9, max_columns=None, log=None):
    """Closures at D = d_min, d_min + 1, ... until 1 appears.  status: 'exact-unit' with the least
    D containing 1, or 'lower-bound' when the degree cap or the column cap stops the ladder."""
    profile, status, value = [], 'lower-bound', None
    for D in range(d_min, d_max + 1):
        n = sum(math.comb(m * k, w) * math.comb(D - w + m - 2, m - 2) for w in range(min(m * k, D) + 1))
        if max_columns is not None and n > max_columns:
            status = f'lower-bound (column cap at D = {D}: {n} > {max_columns})'
            break
        rec = closure(gens, m, k, p, D, log=log)
        assert rec['monomials_le_D'] == n
        profile.append(rec)
        if log:
            log(f'  D={D} cols={n} dim={rec["dim_V_D"]} one={rec["contains_one"]} {rec["seconds"]}s')
        if rec['contains_one']:
            status, value = 'exact-unit', D
            break
    lower = profile[-1]['D'] + (0 if profile and profile[-1]['contains_one'] else 1) if profile else d_min
    return dict(status=status, exact_lfd=value, lfd_lower_bound=lower, profile=profile)
