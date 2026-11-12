#!/usr/bin/env python3
"""Depth dependence of the direct three-summand prime chain system {3k chains, NF(S_4(x1,x2,x3,t))}.

Engine (``CompactClosure``): the same object as csd.ic_bounds_m3_batched.batched_closure (the
iterated degree-D closure V_D modulo the chains, columns = squarefree chain masks sorted by degree
descending), stored compactly.  The basis is kept in reduced row echelon form; since its pivot
columns form an identity block, only the non-pivot columns are stored (uint16 residues, p < 2^16),
and every reduction costs rank x (#non-pivot columns) instead of rank x |M_<=D|.  This is what makes
|M_<=7| = 16384 (k = 5) and 63004 (k = 6) fit.  Products are float64 with every accumulated sum
below 2^53 (asserted from n (p-1)^2 < 2^53).

Soundness (no early-stop shortcut is used unless stated):
  * every row added is an input or a variable shift x*h of a span element h of degree < D, so the
    span S is always contained in V_D (induction);
  * the run ends only after a full verification sweep in which EVERY variable shift of EVERY final
    basis row of degree < D reduces to zero (and the input reduces to zero): then S is closed under
    admissible products and contains the input, so V_D is contained in S.  Hence S = V_D exactly;
  * optional dimension sandwich: if dim(I cap M_<=D) is known exactly (|M_<=D| minus the rank of
    the evaluation of M_<=D on the solution set Z, R/I = F_p^Z because the chain ideal is radical
    with zero set H^3), the run may stop as soon as rank S reaches it (S subset V_D subset I cap
    M_<=D).  Rows record which check closed the run ("closed-sweep" or "rank=dim(I cap M_<=D)").
  * 1 in S (constant column a pivot) gives V_D = M_<=D.

Profiles: p_j = dim(V_D cap M_<=j) = number of pivots of degree <= j (columns sorted by degree
descending, RREF).

Subcommands:
  selfcheck  cross-check against csd.ic_bounds_m3.exact_profile (k = 2, 3, 4) on fixed instances
  run        closures V_D for given instances / depths / degrees -> JSONL (resumable)
  analyze    tables, recursion test, falling-behind depths -> summary.json, table.md
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import random
import resource
import sys
import time
from itertools import product
from math import comb

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from csd import elim  # noqa: E402
from csd import ic_dreg as I  # noqa: E402
from csd import ic_bounds_m3 as M  # noqa: E402

TWO53 = 2 ** 53
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments', 'ic-bounds-m3-depth')


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


# ---------------------------------------------------------------------------------------------
# exact mod-p helpers
# ---------------------------------------------------------------------------------------------

def rref_rows(E, piv, p):
    """E: int64 row echelon rows (unit pivots at columns piv, increasing). Returns RREF (int64)."""
    s = len(piv)
    if s <= 48:
        E = E.copy()
        for i in range(s - 1, 0, -1):
            c = E[:i, piv[i]].copy()
            if np.any(c):
                E[:i] = (E[:i] - np.outer(c, E[i])) % p
        return E
    h = s // 2
    top, bot = E[:h], rref_rows(E[h:], piv[h:], p)
    X = top[:, piv[h:]].astype(np.float64)
    top = np.mod(top.astype(np.float64) - X @ bot.astype(np.float64), p).astype(np.int64)
    top = rref_rows(top, piv[:h], p)
    return np.vstack([top, bot])


class CompactClosure:
    def __init__(self, space, p, batch=1024, chunk=4096, store=np.uint16):
        n = len(space)
        if n * (p - 1) ** 2 >= TWO53:
            raise ValueError('float64 exactness bound violated')
        assert p < 65536
        self.sp, self.p, self.n = space, p, n
        self.batch, self.chunk, self.store = batch, chunk, store
        self.deg = space.deg
        self.D = space.D
        self.nvar = len(space.cols)
        self.piv = []                       # pivot column of row i
        self.npc = np.arange(n)             # columns stored in N (non-pivot + not-yet-compacted pivots)
        self.inN = np.ones(n, dtype=bool)
        self.blocks = []                    # row blocks of N (capacity `chunk` rows each)
        self.nrows = 0
        self.dead = 0
        self.queue = []
        self.qpos = 0
        self.contains_one = False

    @property
    def rank(self):
        return self.n if self.contains_one else len(self.piv)

    def _block_rows(self, bi):
        return min(self.chunk, self.nrows - bi * self.chunk)

    def _matmul_N(self, X):
        """X (b x rank) float64 @ N (rank x |npc|), accumulated exactly; returns float64 (not reduced)."""
        acc = np.zeros((X.shape[0], len(self.npc)), dtype=np.float64)
        for bi, blk in enumerate(self.blocks):
            c = self._block_rows(bi)
            o = bi * self.chunk
            acc += X[:, o:o + c] @ blk[:c].astype(np.float64)
        return acc

    def row(self, i):
        return self.blocks[i // self.chunk][i % self.chunk]

    def reduce(self, P):
        """P: (b x n) float64 residues -> reduced rows in npc coordinates (float64 residues)."""
        Q = P[:, self.npc]
        if self.piv:
            X = P[:, np.asarray(self.piv)]
            Q = Q - self._matmul_N(X)
        return np.mod(Q, self.p)

    def _append(self, E):
        E = E.astype(self.store)
        pos = 0
        while pos < len(E):
            bi, off = divmod(self.nrows, self.chunk)
            if bi == len(self.blocks):
                self.blocks.append(np.zeros((self.chunk, len(self.npc)), dtype=self.store))
            take = min(self.chunk - off, len(E) - pos)
            self.blocks[bi][off:off + take] = E[pos:pos + take]
            pos += take
            self.nrows += take

    def absorb(self, Q):
        p = self.p
        nz = np.flatnonzero(np.any(Q != 0, axis=1))
        if not len(nz):
            return 0
        out = []
        s, lpiv, _ = elim.fp_echelon(Q[nz].astype(np.int64), p, out)
        if not s:
            return 0
        E = rref_rows(out[0], lpiv, p)
        lpiv = np.asarray(lpiv, dtype=np.int64)
        # eliminate the new pivot columns from the old rows (block by block, exact)
        Ef = E.astype(np.float64)
        for bi, blk in enumerate(self.blocks):
            cnt = self._block_rows(bi)
            c = blk[:cnt, lpiv].astype(np.float64)
            nzr = np.flatnonzero(np.any(c != 0, axis=1))
            if len(nzr):
                blk[nzr] = np.mod(blk[nzr].astype(np.float64) - c[nzr] @ Ef, p).astype(self.store)
        gpiv = self.npc[lpiv]
        self._append(E)
        for gi in gpiv.tolist():
            idx = len(self.piv)
            self.piv.append(gi)
            if gi == self.n - 1:
                self.contains_one = True
            if self.deg[gi] < self.D:
                self.queue.extend((idx, j) for j in range(self.nvar))
        self.dead += s
        if self.dead > max(512, len(self.npc) // 8):
            self.compact()
        return s

    def compact(self):
        pivset = np.zeros(self.n, dtype=bool)
        pivset[np.asarray(self.piv, dtype=np.int64)] = True
        keep = ~pivset[self.npc]
        self.inN[self.npc[~keep]] = False
        for bi in range(len(self.blocks)):
            self.blocks[bi] = np.ascontiguousarray(self.blocks[bi][:, keep])
        self.npc = self.npc[keep]
        self.dead = 0

    def stored_entries(self):
        return self.nrows * len(self.npc)

    def products(self, items):
        P = np.zeros((len(items), self.n), dtype=np.float64)
        cols = self.sp.cols
        for r, (i, j) in enumerate(items):
            row = self.row(i)
            sup = np.flatnonzero(row)
            dest = cols[j, self.npc[sup]]
            P[r, dest] = row[sup]
            g = self.piv[i]
            if not self.inN[g]:
                P[r, cols[j, g]] = 1.0
        return P

    def full_row(self, i):
        v = np.zeros(self.n, dtype=np.int64)
        v[self.npc] = self.row(i)
        v[self.piv[i]] = 1
        return v

    def add_rows(self, rows):
        if len(rows):
            self.absorb(self.reduce(np.asarray(rows, dtype=np.float64) % self.p))

    def saturate(self, target_rank=None, logevery=50):
        rounds = 0
        t0 = time.time()
        while self.qpos < len(self.queue) and not self.contains_one:
            if target_rank is not None and len(self.piv) >= target_rank:
                assert len(self.piv) == target_rank, (len(self.piv), target_rank)
                return 'rank=dim(I cap M_<=D)'
            items = self.queue[self.qpos:self.qpos + self.batch]
            self.qpos += len(items)
            self.absorb(self.reduce(self.products(items)))
            rounds += 1
            if rounds % logevery == 0:
                log(f'   D={self.D} batches={rounds} rank={len(self.piv)} queue={len(self.queue) - self.qpos} '
                    f'stored={self.nrows}x{len(self.npc)} {time.time() - t0:.0f}s rss={rss_gb():.2f}GB')
            if self.qpos > 200000:
                self.queue = self.queue[self.qpos:]
                self.qpos = 0
        if self.contains_one:
            return 'unit'
        if target_rank is not None and len(self.piv) == target_rank:
            return 'rank=dim(I cap M_<=D)'
        return None

    def sweep(self):
        """Verification sweep: all shifts of all low basis rows; returns number of new pivots."""
        low = [i for i, g in enumerate(self.piv) if self.deg[g] < self.D]
        items = [(i, j) for i in low for j in range(self.nvar)]
        added = 0
        for s in range(0, len(items), self.batch):
            Q = self.reduce(self.products(items[s:s + self.batch]))
            if np.any(Q):
                added += self.absorb(Q)
                if self.contains_one:
                    break
        return added

    def close(self, target_rank=None):
        sweeps = 0
        while True:
            res = self.saturate(target_rank)
            if res is not None:
                return res, sweeps
            added = self.sweep()
            sweeps += 1
            if self.contains_one:
                return 'unit', sweeps
            if added == 0:
                return 'closed-sweep', sweeps
            if target_rank is not None and len(self.piv) == target_rank:
                return 'rank=dim(I cap M_<=D)', sweeps

    def contains(self, vec):
        return not np.any(self.reduce(np.asarray(vec, dtype=np.float64)[None, :] % self.p))

    def low_profile(self):
        if self.contains_one:
            return [int(np.sum(self.deg <= d)) for d in range(self.D + 1)]
        pd = self.deg[np.asarray(self.piv, dtype=np.int64)] if self.piv else np.zeros(0, dtype=np.int64)
        return [int(np.sum(pd <= d)) for d in range(self.D + 1)]


def compact_closure(gens, k, p, D, *, target_rank=None, batch=1024, chunk=4096, check_inputs=True):
    sp = M.MixedSpace(3, k, D, aux=False)
    cl = CompactClosure(sp, p, batch=batch, chunk=chunk)
    rows = []
    for g in gens:
        if g and max(ms.bit_count() for ms in g) <= D:
            rows.append(sp.row(g, p))
    cl.add_rows(rows)
    how, sweeps = cl.close(target_rank)
    if check_inputs and not cl.contains_one:
        for r in rows:
            assert cl.contains(r)
    return cl, how, sweeps


# ---------------------------------------------------------------------------------------------
# instances
# ---------------------------------------------------------------------------------------------

def make_instance(p, seed, planted):
    """Seeded random nonsingular curve and target with 5t^3 + 3at + 6b != 0.
    planted: target = P1 + P2 + P3 with x_i in H_2 (order-4 subgroup), so it decomposes over the
    factor base at every depth k >= 2."""
    rng = random.Random(f'm3depth-{p}-{seed}-{planted}')
    while True:
        a, b = rng.randrange(1, p), rng.randrange(1, p)
        if (4 * a ** 3 + 27 * b ** 2) % p == 0:
            continue
        E = I.PrimeCurve(p, a, b)
        if planted:
            H2 = M.subgroup(p, 2)[1]
            pts = [E.lift_x(x) for x in H2]
            pts = [P for P in pts if P is not None]
            if not pts:
                continue
            P1, P2, P3 = (rng.choice(pts) for _ in range(3))
            P1 = P1 if rng.random() < .5 else E.neg(P1)
            R = E.add(E.add(P1, P2), P3)
            if R is None:
                continue
            t, y = R
            plant = [P1, P2, P3]
        else:
            t = rng.randrange(p)
            P = E.lift_x(t)
            if P is None:
                continue
            y = P[1]
            plant = None
        if (5 * t ** 3 + 3 * a * t + 6 * b) % p == 0:
            continue
        return dict(p=p, a=a, b=b, t=t, y=y, planted=plant, seed=seed,
                    target_class='planted-H2' if planted else 'random-rational')


INSTANCES = {
    # name: (p, seed, planted)
    'p97-r': (97, 1, False),
    'p193-r': (193, 1, False),
    'p257-r': (257, 1, False),
    'p7681-r': (7681, 1, False),
    'p97-pl': (97, 2, True),
    'p193-pl': (193, 2, True),
    'p257-pl': (257, 2, True),
    'p7681-pl': (7681, 2, True),
    # extra k <= 5 instances (32 | p - 1)
    'p353-r': (353, 3, False),
    'p449-pl': (449, 3, True),
    'p577-r': (577, 3, False),
    'p641-pl': (641, 3, True),
}


# ---------------------------------------------------------------------------------------------
# one (instance, k, D) closure
# ---------------------------------------------------------------------------------------------

def case_data(inst, k):
    p = inst['p']
    E = I.PrimeCurve(p, inst['a'], inst['b'])
    g = M.direct_generator(E, inst['t'], k)
    pts = M.direct_points(E, inst['t'], k, g)
    g0 = M.subgroup(p, k)[0]
    return E, g, pts, g0


def closure_row(name, inst, k, D, *, use_target=True, batch=1024, chunk=4096):
    p = inst['p']
    E, g, pts, g0 = case_data(inst, k)
    sp = M.MixedSpace(3, k, D, aux=False)
    masks = [ms for ms, _ in sp.mons]
    n = len(masks)
    t0 = time.time()
    ev_le = [M.evaluation_rank(pts, [ms for ms in masks if ms.bit_count() <= j], k, g0, p) if pts else 0
             for j in range(D + 1)]
    Mle = [sum(1 for ms in masks if ms.bit_count() <= j) for j in range(D + 1)]
    ideal_le = [Mle[j] - ev_le[j] for j in range(D + 1)]   # dim(I cap M_<=j)
    target = ideal_le[D] if use_target else None
    cl, how, sweeps = compact_closure([g], k, p, D, target_rank=target, batch=batch, chunk=chunk)
    prof = cl.low_profile()
    assert all(prof[j] <= ideal_le[j] for j in range(D + 1))
    row = dict(instance=name, p=p, a=inst['a'], b=inst['b'], t=inst['t'], k=k, D=D,
               monomials_le_j=Mle, profile=prof, ideal_dim_le_j=ideal_le, eval_rank_le_j=ev_le,
               solutions=len(pts), rank=cl.rank, contains_one=cl.contains_one,
               equals_ideal=(prof[D] == ideal_le[D]), closure_check=how, sweeps=sweeps,
               generator_degree=max(ms.bit_count() for ms in g), generator_terms=len(g),
               generator_hash=hashlib.sha256(json.dumps(sorted(g.items())).encode()).hexdigest()[:16],
               seconds=round(time.time() - t0, 1), max_rss_gb=round(rss_gb(), 2), engine='CompactClosure')
    return row


def load_rows(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return [json.loads(l) for l in f if l.strip()]


def cmd_run(args):
    path = args.out
    done = {(r['instance'], r['k'], r['D']) for r in load_rows(path)}
    for name in args.instances:
        inst = make_instance(*INSTANCES[name])
        for k in args.k:
            for D in args.D:
                if (name, k, D) in done:
                    continue
                n = len(M.MixedSpace(3, k, D, aux=False).mons)
                if n > args.max_columns:
                    log(f'skip {name} k={k} D={D}: n={n} > {args.max_columns}')
                    continue
                log(f'start {name} k={k} D={D} n={n}')
                row = closure_row(name, inst, k, D, use_target=not args.no_target, batch=args.batch, chunk=args.chunk)
                log(f'done {name} k={k} D={D} profile={row["profile"]} ideal={row["ideal_dim_le_j"][-1]} '
                    f'|Z|={row["solutions"]} check={row["closure_check"]} {row["seconds"]}s rss={row["max_rss_gb"]}GB')
                with open(path, 'a') as f:
                    f.write(json.dumps(row) + '\n')
    log('done')


# ---------------------------------------------------------------------------------------------
# selfcheck against the existing engines
# ---------------------------------------------------------------------------------------------

def cmd_selfcheck(args):
    out = []
    for name in args.instances:
        inst = make_instance(*INSTANCES[name])
        E = I.PrimeCurve(inst['p'], inst['a'], inst['b'])
        for k in args.k:
            ref = M.exact_profile('direct', E, inst['t'], k, max_columns=4000)
            for prow in ref['profile']:
                D = prow['D']
                if D < 1:
                    continue
                for use_target in (False, True):
                    row = closure_row(name, inst, k, D, use_target=use_target)
                    ok = row['rank'] == prow['quotient_rank'] and (
                        D == 0 or row['profile'][D - 1] == prow['quotient_low_rank'])
                    out.append(dict(instance=name, k=k, D=D, target=use_target, ours=row['rank'],
                                    ref=prow['quotient_rank'], ours_low=row['profile'][D - 1],
                                    ref_low=prow['quotient_low_rank'], ref_engine=prow['engine'],
                                    check=row['closure_check'], agree=ok))
                    log(name, k, D, use_target, row['rank'], prow['quotient_rank'], ok)
                    assert ok
            out.append(dict(instance=name, k=k, ref_result=ref['results']['direct']))
    with open(os.path.join(OUTDIR, 'selfcheck-k' + ''.join(map(str, args.k)) + '.json'), 'w') as f:
        json.dump(dict(rows=out, all_agree=all(r.get('agree', True) for r in out)), f, indent=1)
    log('selfcheck ok', len(out))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('selfcheck')
    s.add_argument('--instances', nargs='+', default=['p97-r', 'p257-pl'])
    s.add_argument('--k', nargs='+', type=int, default=[2, 3])
    r = sub.add_parser('run')
    r.add_argument('--instances', nargs='+', required=True)
    r.add_argument('--k', nargs='+', type=int, required=True)
    r.add_argument('--D', nargs='+', type=int, required=True)
    r.add_argument('--out', default=os.path.join(OUTDIR, 'closures.jsonl'))
    r.add_argument('--max-columns', type=int, default=70000)
    r.add_argument('--batch', type=int, default=1024)
    r.add_argument('--chunk', type=int, default=4096)
    r.add_argument('--no-target', action='store_true')
    a = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)
    {'selfcheck': cmd_selfcheck, 'run': cmd_run}[a.cmd](a)


if __name__ == '__main__':
    main()
