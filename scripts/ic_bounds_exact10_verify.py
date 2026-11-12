#!/usr/bin/env python3
"""Independent verification of item 4 ("exact value at depth 10") of the bounds program.

Own code: an exact degree-D chain closure over F_p written from the definition (V_D = least
subspace containing the inputs of degree <= D, closed under h -> u*h when deg h + 1 <= D), in the
chain quotient F_p[u_1,u_2]/(u_i^{2^k}-1) with squarefree chain masks as columns.  It does NOT
import csd/ic_bounds_exact10.py, csd/ic_chain_block_closure.py, csd/ic_fall_certificates.py,
csd/ic_chain_shift_closure.py or any producer script.  The only repository import is
csd/ic_dreg.py (PrimeCurve, s3_prime_coeffs: curve arithmetic and the S_3 coefficients).

Columns: masks e1 | e2 << k of popcount <= D, sorted by (weight descending, mask).  Rows are kept in
reduced row echelon form with leftmost pivots, so the rows whose pivot has weight < D span
V cap M_<D.  Multiplying by the chain variable u_{i,j} adds 2^j to block i modulo 2^k.

Two ways a degree is decided (both exact):
  'closed'  : the worklist of low rows is exhausted.  Every inserted row with pivot weight < D was
              shifted by all 2k variables and every shift was reduced and inserted, so the span is
              closed and contains the generator: it is V_D (every inserted row is in V_D).
  'reached' : rank = target = |M_<=D| - rank ev(M_<=D -> F_p^Z).  Every inserted row lies in V_D
              (shifts of low elements of V_D, and random combinations of such shifts), and
              V_D subset I cap M_<=D, whose dimension is the target (split semisimple quotient).
Random mixing (option --mix) only produces elements of V_D faster; a stall under mixing falls back
to the deterministic worklist over all current low rows, so no stall is ever reported as a bound.

Arithmetic: uint16 storage (p < 2^16), float64 matrix products in tiles (exact: r p^2 < 2^53).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import random
import resource
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from csd import ic_dreg as I  # noqa: E402  (curve arithmetic and S_3 coefficients only)

OUT = ROOT / 'experiments/ic-bounds-exact10/verify'


def maxrss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024


# ----------------------------------------------------------------------------- instances

def my_subgroup_generator(p, k):
    """First h = 2, 3, ... with g = h^((p-1)/2^k) of exact order 2^k (same convention as the lab)."""
    N = 1 << k
    assert (p - 1) % N == 0
    for h in range(2, p):
        g = pow(h, (p - 1) // N, p)
        if pow(g, N // 2, p) != 1:
            return g
    raise ValueError


def my_generator(p, a, b, x, k, c):
    """NF(S_3(c t1, c t2, x)) modulo t_i^{2^k} = 1 as {mask: coeff}; t^e has chain mask e (e <= 2)."""
    E = I.PrimeCurve(p, a, b)
    g = {}
    for (e1, e2), v in I.s3_prime_coeffs(E, x).items():
        m = e1 | (e2 << k)
        g[m] = (g.get(m, 0) + v * pow(c, e1 + e2, p)) % p
    return {m: v for m, v in g.items() if v}


def my_points(p, a, b, x, k, c, g0):
    """All (E1, E2) in [0,2^k)^2 with S_3(c g0^E1, c g0^E2, x) = 0, vectorised."""
    N = 1 << k
    pw = np.array([pow(g0, e, p) for e in range(N)], dtype=np.int64)
    x1 = (c * pw) % p
    X1 = x1[:, None]
    X2 = x1[None, :]
    A, B = a % p, b % p
    # S_3 = (x1-x2)^2 x^2 - 2((x1+x2)(x1 x2 + a) + 2b) x + (x1 x2 - a)^2 - 4b(x1+x2)
    d = (X1 - X2) % p
    s = (X1 + X2) % p
    pr = (X1 * X2) % p
    t1 = (d * d % p) * (x * x % p) % p
    t2 = (2 * x % p) * ((s * ((pr + A) % p) + 2 * B) % p) % p
    t3 = ((pr - A) % p) ** 2 % p
    t4 = (4 * B % p) * s % p
    val = (t1 - t2 + t3 - t4) % p
    E1, E2 = np.nonzero(val == 0)
    # cross-check with ic_dreg's scalar formula on a few points
    E = I.PrimeCurve(p, a, b)
    for i in range(min(5, len(E1))):
        assert I.s3_prime_value(E, int(x1[E1[i]]), int(x1[E2[i]]), x) == 0
    return [(int(u), int(v)) for u, v in zip(E1, E2)]


def ghash(g):
    return hashlib.sha256(json.dumps(sorted(g.items())).encode()).hexdigest()[:16]


# ----------------------------------------------------------------------------- linear algebra

def inv_mod(a, p):
    return pow(int(a), p - 2, p)


def rref_small(M, p):
    """RREF (leftmost pivots) of an int64 matrix with few rows. Returns (rows, pivots)."""
    M = np.asarray(M, dtype=np.int64) % p
    M = M[np.any(M != 0, axis=1)].copy()
    done = np.zeros(M.shape[0], dtype=bool)
    order = []
    n = M.shape[1]
    while True:
        cand = np.nonzero(~done)[0]
        if cand.size == 0:
            break
        sub = M[cand] != 0
        has = sub.any(axis=1)
        if not has.any():
            break
        lead = np.where(has, sub.argmax(axis=1), n)
        j = int(lead.argmin())
        i, c = int(cand[j]), int(lead[j])
        M[i] = (M[i] * inv_mod(M[i, c], p)) % p
        f = M[:, c].copy()
        f[i] = 0
        rows = np.nonzero(f)[0]
        if rows.size:
            M[rows] = (M[rows] - f[rows, None] * M[i][None, :]) % p
        done[i] = True
        order.append((i, c))
    return [M[i].copy() for i, _ in order], [c for _, c in order]


def rank_mod(M, p):
    M = np.array(M, dtype=np.int64) % p
    if M.shape[0] > M.shape[1]:
        M = M.T.copy()
    r = 0
    rows, cols = M.shape
    for c in range(cols):
        if r == rows:
            break
        nz = np.nonzero(M[r:, c])[0]
        if nz.size == 0:
            continue
        i = r + nz[0]
        if i != r:
            M[[r, i]] = M[[i, r]]
        M[r] = (M[r] * inv_mod(M[r, c], p)) % p
        f = M[:, c].copy()
        f[r] = 0
        M = (M - np.outer(f, M[r])) % p
        r += 1
    return r


class Space:
    def __init__(self, k, D):
        self.k, self.D = k, D
        n2 = 2 * k
        allm = np.arange(1 << n2, dtype=np.int64)
        w = np.array([bin(m).count('1') for m in range(1 << n2)], dtype=np.int64)
        keep = allm[w <= D]
        order = np.lexsort((keep, -w[keep]))
        self.masks = keep[order]
        self.w = w[self.masks]
        self.n = len(self.masks)
        self.col = np.full(1 << n2, -1, dtype=np.int64)
        self.col[self.masks] = np.arange(self.n)
        self.low_start = int(np.argmax(self.w < D))
        low = self.masks[self.low_start:]
        N = 1 << k
        e1, e2 = low & (N - 1), low >> k
        self.shift = []
        for i in range(2):
            for j in range(k):
                if i == 0:
                    m = ((e1 + (1 << j)) % N) | (e2 << k)
                else:
                    m = e1 | (((e2 + (1 << j)) % N) << k)
                c = self.col[m]
                assert (c >= 0).all()
                self.shift.append(c)
        self.one_col = int(self.col[0])


class Closure:
    def __init__(self, sp: Space, p, tile=2048, cap=None):
        assert p < (1 << 16)
        self.sp, self.p, self.tile = sp, p, tile
        cap = cap or sp.n
        self.R = np.zeros((cap, sp.n), dtype=np.uint16)
        self.piv = np.zeros(cap, dtype=np.int64)
        self.r = 0
        self.queue = []  # copies (uint16) of inserted rows with pivot weight < D, not yet shifted
        self.inserted = 0

    def reduce(self, B):
        """B float64 (b x n) -> B reduced against the RREF rows, entries in [0,p)."""
        p, r = self.p, self.r
        B %= p
        if r == 0:
            return B
        C = B[:, self.piv[:r]].copy()
        for s in range(0, self.sp.n, self.tile):
            T = self.R[:r, s:s + self.tile].astype(np.float64)
            B[:, s:s + self.tile] -= C @ T
            B[:, s:s + self.tile] %= p
        return B

    def insert_batch(self, B):
        """Reduce candidate rows, insert the new ones; return number of new rows."""
        p = self.p
        B = self.reduce(B.astype(np.float64))
        nzr = np.any(B != 0, axis=1)
        if not nzr.any():
            return 0
        M = B[nzr].astype(np.int64)
        rows, piv = rref_small(M, p)
        if not rows:
            return 0
        Nw = np.array(rows, dtype=np.int64)
        pv = np.array(piv, dtype=np.int64)
        r = self.r
        # back-reduce old rows by the new rows (keeps full RREF)
        if r:
            Nf = Nw.astype(np.float64)
            for s in range(0, r, self.tile):
                blk = self.R[s:min(r, s + self.tile)]
                coef = blk[:, pv].astype(np.float64)
                sel = np.nonzero(coef.any(axis=1))[0]
                if sel.size == 0:
                    continue
                upd = (blk[sel].astype(np.float64) - coef[sel] @ Nf) % p
                self.R[s + sel] = upd.astype(np.uint16)
        t = len(rows)
        assert r + t <= self.R.shape[0]
        self.R[r:r + t] = Nw.astype(np.uint16)
        self.piv[r:r + t] = pv
        self.r = r + t
        lowmask = self.sp.w[pv] < self.sp.D
        for i in np.nonzero(lowmask)[0]:
            self.queue.append(Nw[i].astype(np.uint16))
        self.inserted += t
        return t

    def shifts_of(self, rows):
        """All 2k shifts of rows (each supported in columns of weight < D)."""
        sp = self.sp
        L = np.asarray(rows)[:, sp.low_start:]
        out = np.zeros((len(rows) * len(sp.shift), sp.n), dtype=np.float64)
        for v, cmap in enumerate(sp.shift):
            out[v::len(sp.shift)][:, cmap] = L
        return out

    def low_rows(self):
        idx = np.nonzero(self.sp.w[self.piv[:self.r]] < self.sp.D)[0]
        return idx

    @property
    def low_dim(self):
        return int(len(self.low_rows()))

    @property
    def contains_one(self):
        return bool((self.piv[:self.r] == self.sp.one_col).any())

    def run_queue(self, target=None, batch_rows=None, deadline=None):
        """Deterministic worklist. Returns 'reached' if rank hits target, else 'closed'."""
        if batch_rows is None:
            batch_rows = max(1, 320 // len(self.sp.shift))
        while self.queue:
            if target is not None and self.r >= target:
                return 'reached'
            if deadline and time.time() > deadline:
                return 'timeout'
            take = self.queue[:batch_rows]
            del self.queue[:batch_rows]
            self.insert_batch(self.shifts_of(take))
        if target is not None and self.r >= target:
            return 'reached'
        return 'closed'

    def mix_round(self, b, terms, rng):
        sp = self.sp
        low = self.low_rows()
        if len(low) == 0:
            return 0
        B = np.zeros((b, sp.n), dtype=np.float64)
        nv = len(sp.shift)
        for _ in range(terms):
            ri = low[rng.integers(0, len(low), size=b)]
            vi = rng.integers(0, nv, size=b)
            cf = rng.integers(1, self.p, size=b).astype(np.float64)
            for v in range(nv):
                sel = np.nonzero(vi == v)[0]
                if sel.size == 0:
                    continue
                rows = self.R[ri[sel], sp.low_start:].astype(np.float64) * cf[sel, None]
                B[np.ix_(sel, sp.shift[v])] += rows
        return self.insert_batch(B)


def closure_degree(gen, k, D, p, prev=None, target=None, mix=0, terms=6, seed=1, deadline=None,
                   tile=2048, batch_rows=None):
    """Exact V_D. prev = (Space, Closure) of V_{D-1} to start from (all its rows are in V_D)."""
    sp = Space(k, D)
    cl = Closure(sp, p, tile=tile, cap=(target if target is not None and target > 0 else sp.n))
    if prev is None:
        row = np.zeros((1, sp.n))
        for m, v in gen.items():
            row[0, sp.col[m]] = v
        cl.insert_batch(row)
    else:
        psp, pcl = prev
        cmap = sp.col[psp.masks]
        for s in range(0, pcl.r, 2000):
            blk = pcl.R[s:min(pcl.r, s + 2000)]
            B = np.zeros((blk.shape[0], sp.n))
            B[:, cmap] = blk
            cl.insert_batch(B)
        # every row of V_{D-1} has weight <= D-1 < D: all are low in V_D; queue holds them already
    status = None
    if mix and target is not None:
        rng = np.random.default_rng(seed)
        stall = 0
        while cl.r < target and stall < 2:
            if deadline and time.time() > deadline:
                return sp, cl, 'timeout'
            got = cl.mix_round(mix, terms, rng)
            stall = stall + 1 if got == 0 else 0
        if cl.r >= target:
            status = 'reached'
        else:
            # deterministic fallback: re-shift every current low row
            cl.queue = [cl.R[i].copy() for i in cl.low_rows()]
    if status is None:
        status = cl.run_queue(target=target, batch_rows=batch_rows, deadline=deadline)
    if status == 'closed' and target is not None and cl.r == target:
        status = 'reached'
    return sp, cl, status


def eval_ranks(points, sp, k, p, g0):
    """rank ev(M_<=D) and rank ev(M_<D) on Z."""
    if not points:
        return 0, 0
    N = 1 << k
    pw = np.array([pow(g0, e, p) for e in range(N)], dtype=np.int64)
    m = sp.masks
    a1, a2 = m & (N - 1), m >> k
    E1 = np.array([u for u, _ in points])[:, None]
    E2 = np.array([v for _, v in points])[:, None]
    M = pw[(a1[None, :] * E1 + a2[None, :] * E2) % N]
    full = rank_mod(M, p)
    low = rank_mod(M[:, sp.low_start:], p)
    return full, low


def verify_vanishing(cl, points, k, p, g0, sample=64, seed=0):
    """Soundness check: random rows of V_D vanish on Z (proper ideals)."""
    if not points or cl.r == 0:
        return True
    N = 1 << k
    pw = np.array([pow(g0, e, p) for e in range(N)], dtype=np.int64)
    m = cl.sp.masks
    a1, a2 = m & (N - 1), m >> k
    E1 = np.array([u for u, _ in points])[:, None]
    E2 = np.array([v for _, v in points])[:, None]
    M = pw[(a1[None, :] * E1 + a2[None, :] * E2) % N]  # |Z| x n
    rng = np.random.default_rng(seed)
    idx = rng.choice(cl.r, size=min(sample, cl.r), replace=False)
    rows = cl.R[idx].astype(np.int64)
    val = (rows @ (M.T % p)) % p
    return bool((val == 0).all())


def profile(gen, k, p, points, g0, Dmax, mix=0, deadline=None, log=None, stop_when_reached=False):
    """dims per D=2..Dmax and exact LFD / falls."""
    recs, prev = [], None
    for D in range(2, Dmax + 1):
        sp0 = Space(k, D)
        full, low = eval_ranks(points, sp0, k, p, g0)
        target = sp0.n - full
        t0 = time.time()
        sp, cl, status = closure_degree(gen, k, D, p, prev=prev, target=target, mix=mix, seed=1000 + D,
                                        deadline=deadline)
        rec = dict(D=D, cols=sp.n, dimV=cl.r, lowV=cl.low_dim, one=cl.contains_one, target=target,
                   ev_rank=full, ev_low_rank=low, onto_lt=(low == len(points)), status=status,
                   vanish=verify_vanishing(cl, points, k, p, g0), sec=round(time.time() - t0, 1))
        rec['prevV'] = recs[-1]['dimV'] if recs else None
        rec['fall'] = bool(recs) and rec['lowV'] > recs[-1]['dimV']
        recs.append(rec)
        if log:
            log(rec)
        if status == 'timeout':
            break
        prev = (sp, cl)
        if status == 'reached' and rec['onto_lt'] and stop_when_reached:
            break
    out = dict(recs=recs)
    dstar = [r['D'] for r in recs if r['status'] == 'reached' and r['onto_lt']]
    if dstar:
        neq = [r['D'] for r in recs if r['dimV'] != r['target']]
        out['lfd'] = (max(neq) + 1) if neq else 2
        out['falls'] = [r['D'] for r in recs if r['fall'] and r['D'] >= 3]
    else:
        out['lfd'] = None
    return out, prev


# ----------------------------------------------------------------------------- commands

GRID = [
    dict(name='j1728-p7681-x936', p=7681, a=1, b=0, scale=1),
    dict(name='rand-p257-x58', p=257, a=2, b=3, scale=1),
    dict(name='rand-p7681-x2244', p=7681, a=5, b=7, scale=1),
    dict(name='rand-p12289-x11294', p=12289, a=3, b=11, scale=1),
    dict(name='j0-p7681-seed41', p=7681, a=0, b=7, scale=1),
    dict(name='coset5-p12289-seed42', p=12289, a=3, b=11, scale=5),
    dict(name='proper-p12289-seed43', p=12289, a=1, b=0, scale=1),
    dict(name='proper-coset7-p7681-seed44', p=7681, a=5, b=7, scale=7),
]


def producer_rows():
    rows = {}
    for f in sorted((ROOT / 'experiments/ic-bounds-exact10').glob('validation-*.jsonl')):
        if 'v2' in f.name:
            continue
        for line in open(f):
            r = json.loads(line)
            rows[(r['instance'], r['k'], r['degree'])] = r
    profs = {}
    for f in sorted((ROOT / 'experiments/ic-bounds-exact10').glob('profiles-*.jsonl')):
        if 'v2' in f.name:
            continue
        for line in open(f):
            r = json.loads(line)
            profs[(r['instance'], r['k'])] = r
    return rows, profs


def cmd_grid(args):
    rows, profs = producer_rows()
    OUT.mkdir(parents=True, exist_ok=True)
    out = open(OUT / f'grid{args.tag}.jsonl', 'a')
    deadline = time.time() + args.budget
    for inst in GRID:
        if args.instances and inst['name'] not in args.instances:
            continue
        for k in args.ks:
            ref = rows.get((inst['name'], k, 4))
            if ref is None:
                continue
            p, x, c = inst['p'], ref['target_x'], inst['scale']
            g0 = my_subgroup_generator(p, k)
            gen = my_generator(p, inst['a'], inst['b'], x, k, c)
            Z = my_points(p, inst['a'], inst['b'], x, k, c, g0)
            t0 = time.time()
            prof, _ = profile(gen, k, p, Z, g0, args.dmax, mix=args.mix, deadline=deadline)
            res = dict(instance=inst['name'], p=p, a=inst['a'], b=inst['b'], x=x, scale=c, k=k, g0=g0,
                       points=len(Z), gen_sha=ghash(gen), lfd=prof['lfd'], falls=prof.get('falls'),
                       recs=prof['recs'], sec=round(time.time() - t0, 1), maxrss_mb=maxrss_mb())
            cmp = []
            ok = len(Z) == ref['points']
            cmp.append(('points', len(Z), ref['points'], ok))
            for r in prof['recs']:
                pr = rows.get((inst['name'], k, r['D']))
                if pr is None:
                    continue
                ss = pr.get('split_scratch') or {}
                be = pr.get('block_engine') or {}
                for key, mine, theirs in [('dimV', r['dimV'], pr['ref_dim']),
                                          ('lowV', r['lowV'], pr['ref_low_dim']),
                                          ('one', r['one'], pr['ref_contains_one']),
                                          ('split_dim', r['dimV'], ss.get('exact_dimension')),
                                          ('split_low', r['lowV'], ss.get('low_dimension')),
                                          ('target', r['target'], ss.get('target')),
                                          ('block_rank_eq_when_not_one',
                                           True if r['one'] else (be.get('sym_rank', -1) + be.get('anti_rank', -1) == r['dimV']),
                                           True),
                                          ('block_one', r['one'], be.get('contains_one'))]:
                    cmp.append((f'D{r["D"]}.{key}', mine, theirs, mine == theirs))
                st_map = {'reached': 'reached', 'closed': 'closed-deficient'}
                mine_st = 'reached' if r['dimV'] == r['target'] else st_map.get(r['status'], r['status'])
                cmp.append((f'D{r["D"]}.status', mine_st, ss.get('status'), mine_st == ss.get('status')))
            pf = profs.get((inst['name'], k))
            if pf and prof['lfd'] is not None:
                cmp.append(('lfd', prof['lfd'], pf['split']['exact_lfd'], prof['lfd'] == pf['split']['exact_lfd']))
                cmp.append(('ref_lfd', prof['lfd'], pf['reference']['exact_lfd'], prof['lfd'] == pf['reference']['exact_lfd']))
                cmp.append(('falls', prof['falls'], pf['split']['strict_falls'], prof['falls'] == pf['split']['strict_falls']))
            res['compare'] = cmp
            res['all_agree'] = all(t[3] for t in cmp)
            out.write(json.dumps(res) + '\n')
            out.flush()
            print(json.dumps(dict(instance=inst['name'], k=k, points=len(Z), lfd=prof['lfd'], falls=prof.get('falls'),
                                  dims=[(r['D'], r['dimV'], r['lowV'], r['target'], r['status']) for r in prof['recs']],
                                  agree=res['all_agree'], bad=[t for t in cmp if not t[3]], sec=res['sec'],
                                  rss=maxrss_mb())), flush=True)


def load_system(depth_dir, idx):
    S = [json.loads(l) for l in open(ROOT / depth_dir / 'systems.jsonl')]
    C = [json.loads(l) for l in open(ROOT / depth_dir / 'cases.jsonl')]
    s = S[idx]
    case = [c for c in C if c['system_hash'] == s['system_hash']][0]
    return s, case


def cmd_midscale(args):
    OUT.mkdir(parents=True, exist_ok=True)
    out = open(OUT / f'midscale{args.tag}.jsonl', 'a')
    s, case = load_system(args.depth_dir, args.system)
    p, k = s['prime'], s['k']
    a, b, x, c = case['a'], case['b'], case['target'][0], case['base_scale']
    g0 = my_subgroup_generator(p, k)
    gen = my_generator(p, a, b, x, k, c)
    stored = {m: v % p for m, v in s['generators'][0]}
    # generators may be normalised by a scalar: compare up to scaling
    lam = None
    same = set(stored) == set(gen)
    if same:
        m0 = next(iter(gen))
        lam = stored[m0] * inv_mod(gen[m0], p) % p
        same = all(stored[m] == gen[m] * lam % p for m in gen)
    Z = my_points(p, a, b, x, k, c, g0)
    Zs = sorted(s['coordinate_solutions'])
    Zm = sorted(u | (v << k) for u, v in Z)
    head = dict(depth_dir=args.depth_dir, system=args.system, p=p, k=k, a=a, b=b, x=x, scale=c, g0=g0,
                stored_g0=s['generator_value'], generator_matches_stored=same, points=len(Z),
                points_match_stored=(Zs == Zm), gen_sha=ghash(gen))
    print(json.dumps(head), flush=True)
    deadline = time.time() + args.budget
    prev, recs = None, []
    for D in range(2, args.dmax + 1):
        sp0 = Space(k, D)
        full, low = eval_ranks(Z, sp0, k, p, g0)
        target = sp0.n - full
        t0 = time.time()
        use_mix = args.mix if D >= args.mix_from else 0
        sp, cl, status = closure_degree(gen, k, D, p, prev=prev, target=target, mix=use_mix, seed=2000 + D,
                                        deadline=deadline)
        rec = dict(D=D, cols=sp.n, dimV=cl.r, lowV=cl.low_dim, one=cl.contains_one, target=target,
                   onto_lt=(low == len(Z)), status=status, vanish=verify_vanishing(cl, Z, k, p, g0),
                   prevV=recs[-1]['dimV'] if recs else None, sec=round(time.time() - t0, 1), rss=maxrss_mb())
        recs.append(rec)
        print(json.dumps(rec), flush=True)
        prev = (sp, cl)
        if status == 'timeout':
            break
        if D == args.dmax and args.save:
            np.save(args.save, cl.R[:cl.r])
    head['recs'] = recs
    out.write(json.dumps(head) + '\n')


def cmd_dryrun(args):
    """Instance construction for depth 10, own code (target draw replicates the probe convention:
    random.Random(seed), x = randrange(1, p) until x^3+ax+b is a square; redraw for --want)."""
    res = []
    for p, a, b, seed, want in [(12289, 1, 0, 10, False), (13313, 1, 0, 11, True), (15361, 4021, 9917, 12, True)]:
        k = 10
        E = I.PrimeCurve(p, a, b)
        g0 = my_subgroup_generator(p, k)
        rng = random.Random(seed)
        attempts = 0
        while True:
            attempts += 1
            while True:
                x = rng.randrange(1, p)
                if E.lift_x(x) is not None:
                    break
            Z = my_points(p, a, b, x, k, 1, g0)
            if not want or Z:
                break
        res.append(dict(p=p, a=a, b=b, seed=seed, want=want, x=x, attempts=attempts, points=len(Z), g0=g0,
                        cols_D6=Space(k, 6).n if args.cols else None))
        print(json.dumps(res[-1]), flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'dryrun-k10.json').write_text(json.dumps(res, indent=1) + '\n')


def cmd_selftest(args):
    """Cross-check the chain-quotient engine against the independent ordinary-ring validator
    experiments/ic-curve-validation/code/lfd_engine.py (full polynomial ring with chain inputs)."""
    sys.path.insert(0, str(ROOT / 'experiments/ic-curve-validation/code'))
    import lfd_engine as LE
    OUT.mkdir(parents=True, exist_ok=True)
    res = []
    cases = [(257, 2, 3, 3), (7681, 1, 0, 3), (7681, 5, 7, 3), (257, 1, 0, 4), (7681, 0, 7, 4), (12289, 3, 11, 4)]
    for p, a, b, k in cases:
        E = I.PrimeCurve(p, a, b)
        rng = random.Random(p * 31 + k)
        for _ in range(args.per_case):
            while True:
                x = rng.randrange(1, p)
                if E.lift_x(x) is not None:
                    break
            g0 = my_subgroup_generator(p, k)
            gen = my_generator(p, a, b, x, k, 1)
            Z = my_points(p, a, b, x, k, 1, g0)
            mine, _ = profile(gen, k, p, Z, g0, args.dmax)
            inputs, G = LE.prime_inputs(sorted(gen.items()), k, p)
            pts = LE.prime_points(G, k, p)
            t0 = time.time()
            ref = LE.exact_lfd(inputs, 2 * k, p, pts, args.dmax)
            rec = dict(p=p, a=a, b=b, k=k, x=x, points=len(Z), ref_points=len(pts), my_lfd=mine['lfd'],
                       ref_lfd=ref['lfd'], my_falls=mine.get('falls'), ref_falls=[d for d in ref['falls'] if d >= 3],
                       my_one=[r['one'] for r in mine['recs']],
                       ref_one=[r['one_in_V'] for r in ref['per_degree'] if r['D'] >= 2],
                       sec=round(time.time() - t0, 1))
            rec['agree'] = (rec['my_lfd'] == rec['ref_lfd'] and rec['points'] == rec['ref_points']
                            and rec['my_falls'] == rec['ref_falls'])
            res.append(rec)
            print(json.dumps(rec), flush=True)
    (OUT / 'selftest.json').write_text(json.dumps(res, indent=1) + '\n')


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    g = sub.add_parser('grid')
    g.add_argument('--ks', type=int, nargs='+', default=[5, 6])
    g.add_argument('--instances', nargs='*', default=None)
    g.add_argument('--dmax', type=int, default=6)
    g.add_argument('--mix', type=int, default=0)
    g.add_argument('--budget', type=float, default=1400)
    g.add_argument('--tag', default='')
    m = sub.add_parser('midscale')
    m.add_argument('--depth-dir', default='experiments/ic-chain-depth8')
    m.add_argument('--system', type=int, default=0)
    m.add_argument('--dmax', type=int, default=5)
    m.add_argument('--mix', type=int, default=0)
    m.add_argument('--mix-from', type=int, default=5)
    m.add_argument('--budget', type=float, default=1400)
    m.add_argument('--save', default=None)
    m.add_argument('--tag', default='')
    t = sub.add_parser('selftest')
    t.add_argument('--dmax', type=int, default=5)
    t.add_argument('--per-case', type=int, default=2)
    d = sub.add_parser('dryrun')
    d.add_argument('--cols', action='store_true')
    args = ap.parse_args()
    dict(grid=cmd_grid, midscale=cmd_midscale, dryrun=cmd_dryrun, selftest=cmd_selftest)[args.cmd](args)


if __name__ == '__main__':
    main()
