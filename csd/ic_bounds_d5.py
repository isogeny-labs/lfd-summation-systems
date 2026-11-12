"""Symmetry-split, free-column closure engine for the degree-D chain closure V_D (item 3: d_5(k)).

Setting (as in csd/ic_chain_shift_closure.py).  Two summands, subgroup chains of depth k over F_p,
normalized generator g (chain degree <= 2).  Modulo the chains a squarefree monomial is a pair of
block exponents (a1, a2) in (Z/2^k)^2, its ordinary degree is the Hamming weight, and multiplication
by the chain variable u_{i,j} adds 2^j to a_i.  V_D is the least subspace of M_<=D containing the
rows NF(mu*g) (wt mu + deg g <= D) and closed under h -> u*h for h of weight < D.

Three ideas make depths 11-13 reachable where the depth-10 probe needed 8 GB:

1. Symmetry split.  The summand swap sigma fixes g, so V_D = V^+ (+) V^- (char p != 2) and both parts
   are kept separately on about |M_<=D|/2 orbit coordinates each.  For f in V^+- of weight < D the
   block-2 shift of f is +-sigma of the block-1 shift, so only the k block-1 variables are applied;
   each shifted row is projected to both parts (both projections lie in V_D because V_D is sigma-
   invariant).

2. Free-column storage.  Each part is a reduced row echelon basis (RREF) stored as its pivot columns
   and the coefficients on the non-pivot ("free") columns only, an r x (c - r) array instead of r x c.
   The full reduction of a row y is zero on the pivot columns, so the reduced row is
   y[free] - y[pivots] @ B_free: one exact float64 product of cost r * (c - r) per row.  At the
   depth-10 point (c = 21,700, d_5 = 18,569) this stores about 15 M entries instead of 400 M.

3. Exact tiled products.  Every product is a float64 GEMM whose inner dimension is tiled so that all
   partial sums stay below 2^53 (exact for residues mod p), reduced mod p after each tile; the
   storage is uint16 for p < 2^16.

Certificate (deterministic, not probabilistic).  Every row ever added lies in V_D: generator rows,
the variable shifts of low rows (weight < D) of the current space, and sigma-projections of those.
A closure phase processes every low pivot once; then a verification sweep re-reduces the k block-1
shifts of every final low row of both parts.  If the sweep adds nothing, the final space W contains
the generator rows and is closed under every admissible product (block-2 shifts by the sigma
symmetry), so V_D is contained in W; with W contained in V_D this gives W = V_D exactly.  If the
sweep adds a row, the closure phase resumes and a new sweep follows.  Reaching rank |M_<=D| - |Z|
earlier is the sound upper certificate of the shift engine (V_D = I cap M_<=D when M_<D evaluates
onto Z).

The per-weight pivot counts give the filtration profile dim(V_D cap M_<=w), w = 0..D, which is
recorded with every run (the leading-weight structure of the closure).
"""
from __future__ import annotations

import json
import math
import os
import time
from itertools import combinations
from pathlib import Path

import numpy as np

from . import elim
from . import ic_dreg as I
from . import ic_curve_sweep as S

TWO53 = 2 ** 53


# ---------------------------------------------------------------------------
# instances
# ---------------------------------------------------------------------------


def prime_points_fast(E: I.PrimeCurve, k: int, generator_value: int, xR: int, chunk: int = 256) -> list[int]:
    """Vectorised ic_curve_sweep.prime_points for the subgroup presentation (c1 = c2 = 1):
    every E1 | E2 << k with S_3(g^E1, g^E2, xR) = 0, in the same order."""
    p = E.p
    N = 1 << k
    H = np.asarray([pow(generator_value, e, p) for e in range(N)], dtype=np.int64)
    pw = {0: np.ones(N, dtype=np.int64), 1: H, 2: H * H % p}
    coeffs = I.s3_prime_coeffs(E, xR)
    out = []
    for s in range(0, N, chunk):
        e1 = np.arange(s, min(N, s + chunk))
        acc = np.zeros((len(e1), N), dtype=np.int64)
        for (a1, a2), c in coeffs.items():
            acc = (acc + (c * pw[a1][e1] % p)[:, None] * pw[a2][None, :]) % p
        r1, r2 = np.nonzero(acc == 0)
        out.extend((e1[r1] | (r2 << k)).tolist())
    return out


def evaluation_rank_fast(points: list[int], masks: list[int], k: int, generator_value: int, p: int) -> int:
    """Vectorised ic_curve_sweep.prime_evaluation_rank."""
    if not points or not masks:
        return 0
    N = 1 << k
    table = np.asarray([pow(generator_value, e, p) for e in range(N)], dtype=np.int64)
    pts = np.asarray(points, dtype=np.int64)
    E1, E2 = pts & (N - 1), pts >> k
    ms = np.asarray(masks, dtype=np.int64)
    a1, a2 = ms & (N - 1), ms >> k
    rank = 0
    mat = table[(E1[:, None] * a1[None, :] + E2[:, None] * a2[None, :]) % N]
    rank = elim.fp_echelon(mat, p)[0]
    return rank


def make_instance(p: int, a: int, b: int, seed: int, k: int, want_solutions: bool = False) -> dict:
    """The instance of scripts/run_ic_chain_depth10_probe.py, for any depth k: target drawn from
    random.Random(seed) by ic_curve_sweep.random_rational_target, subgroup presentation."""
    import random
    assert (p - 1) % (1 << k) == 0 and (4 * a ** 3 + 27 * b * b) % p
    rng = random.Random(seed)
    E = I.PrimeCurve(p, a, b)
    g0 = S.subgroup_generator_value(p, k)
    for attempt in range(200):
        R = S.random_rational_target(E, rng, p)
        Z = prime_points_fast(E, k, g0, R[0])
        if Z or not want_solutions:
            break
    else:
        raise SystemExit('No target with scalar solutions found.')
    g = I.pkm_generator(E, R[0], k)
    return dict(E=E, g0=g0, target=list(R), points=Z, generator=g, attempts=attempt + 1)


# ---------------------------------------------------------------------------
# monomials, shifts, symmetry
# ---------------------------------------------------------------------------


def monomials(k: int, D: int) -> list[int]:
    """Same order as ic_chain_shift_closure.monomials: weight descending, then mask descending."""
    r = 2 * k
    masks = [sum(1 << j for j in c) for d in range(min(r, D) + 1) for c in combinations(range(r), d)]
    return sorted(masks, key=lambda m: (m.bit_count(), m), reverse=True)


def column_count(k: int, D: int) -> int:
    return sum(math.comb(2 * k, d) for d in range(min(2 * k, D) + 1))


class Layout:
    """Columns of M_<=D, the k block-1 shift maps restricted to low columns, and the sigma split."""

    def __init__(self, k: int, D: int):
        self.k, self.D = k, D
        N = 1 << k
        self.masks = monomials(k, D)
        self.index = {m: i for i, m in enumerate(self.masks)}
        self.cols = len(self.masks)
        ms = np.asarray(self.masks, dtype=np.int64)
        self.weights = np.asarray([m.bit_count() for m in self.masks])
        self.low = np.flatnonzero(self.weights < D)            # columns of weight < D
        # block-1 shift v = j: (a1, a2) -> (a1 + 2^j mod N, a2); only low columns are ever shifted
        self.shift_low = np.zeros((k, len(self.low)), dtype=np.int64)
        lowm = ms[self.low]
        for j in range(k):
            new = (((lowm & (N - 1)) + (1 << j)) % N) | (lowm & ~np.int64(N - 1))
            self.shift_low[j] = [self.index[int(x)] for x in new]
        swap = np.asarray([self.index[(m >> k) | ((m & (N - 1)) << k)] for m in self.masks])
        self.swap = swap
        reps = np.asarray([i for i in range(self.cols) if i <= swap[i]])
        self.sym_reps = reps
        self.sym_partner = swap[reps]
        self.anti_reps = np.asarray([i for i in reps if swap[i] != i], dtype=np.int64)
        self.anti_partner = swap[self.anti_reps]
        self.part_cols = (len(self.sym_reps), len(self.anti_reps))
        self.part_weights = (self.weights[self.sym_reps], self.weights[self.anti_reps])
        one = self.cols - 1
        self.one_sym = int(np.flatnonzero(self.sym_reps == one)[0])

    def project(self, full: np.ndarray, p: int):
        """(sym coords, anti coords) of full rows: f_m + f_sigma(m) and f_m - f_sigma(m) on representatives."""
        s = (full[:, self.sym_reps] + full[:, self.sym_partner]) % p
        a = (full[:, self.anti_reps] - full[:, self.anti_partner]) % p
        return s, a

    def lift(self, part: int, coords: np.ndarray, p: int) -> np.ndarray:
        """A full vector whose projection to `part` is `coords` (and zero in the other part)."""
        inv2 = pow(2, -1, p)
        out = np.zeros((coords.shape[0], self.cols), dtype=np.int64)
        c = coords.astype(np.int64) * inv2 % p
        if part == 0:
            out[:, self.sym_reps] = c
            out[:, self.sym_partner] = c   # diagonal columns get c/2 once: projection 2*(c/2) = c
        else:
            out[:, self.anti_reps] = c
            out[:, self.anti_partner] = (-c) % p
        return out

    def shifted_block1(self, full_low: np.ndarray, j: int) -> np.ndarray:
        """u_{1,j} * f for rows f supported on low columns (given on self.low coordinates)."""
        out = np.zeros((full_low.shape[0], self.cols), dtype=np.int64)
        out[:, self.shift_low[j]] = full_low
        return out

    def generator_rows(self, g: dict[int, int], p: int, mus: list[int]) -> np.ndarray:
        N = 1 << self.k
        rows = np.zeros((len(mus), self.cols), dtype=np.int64)
        for r, mu in enumerate(mus):
            m1, m2 = mu & (N - 1), mu >> self.k
            for mask, c in g.items():
                col = self.index[((m1 + (mask & (N - 1))) % N) | (((m2 + (mask >> self.k)) % N) << self.k)]
                rows[r, col] = (rows[r, col] + c) % p
        return rows


# ---------------------------------------------------------------------------
# exact products and the free-column RREF
# ---------------------------------------------------------------------------


def mulmod(A: np.ndarray, B: np.ndarray, p: int, tile: int = 1024) -> np.ndarray:
    """(A @ B) mod p for residues, exact: the inner dimension is cut into tiles with
    tile * (p-1)^2 + p < 2^53, and the float64 accumulator is reduced after every tile."""
    n, m = A.shape
    step = max(1, min(tile, (TWO53 - p) // ((p - 1) ** 2)))
    out = np.zeros((n, B.shape[1]), dtype=np.float64)
    for s in range(0, m, step):
        e = min(m, s + step)
        out += A[:, s:e].astype(np.float64) @ B[s:e].astype(np.float64)
        np.fmod(out, p, out=out)
    return out.astype(np.int64)


class FreeRREF:
    """RREF basis of a subspace of F_p^c: row i is e_{piv[i]} + sum_j B[i, j] e_{free[j]}.

    Rows are in insertion order (not sorted by pivot).  Entries of B left of a row's pivot are zero
    (row echelon), which is preserved by every update.
    """

    def __init__(self, c: int, p: int):
        self.c, self.p = c, p
        self.dtype = np.uint16 if p < 65536 else np.uint32
        self.piv = np.zeros(0, dtype=np.int64)
        self.free = np.arange(c, dtype=np.int64)
        self.B = np.zeros((0, c), dtype=self.dtype)
        self.madds = 0  # multiply-adds of the exact products (machine-independent work counter)

    @property
    def rank(self) -> int:
        return len(self.piv)

    def residual(self, Y: np.ndarray, chunk: int = 1024) -> np.ndarray:
        """Reduced rows restricted to the free columns (the reduced rows vanish on the pivots)."""
        p = self.p
        Y = np.asarray(Y, dtype=np.int64) % p
        if not self.rank:
            return Y[:, self.free]
        out = np.empty((len(Y), len(self.free)), dtype=np.int64)
        for s in range(0, len(Y), chunk):
            y = Y[s:s + chunk]
            out[s:s + chunk] = (y[:, self.free] - mulmod(y[:, self.piv], self.B, p)) % p
        self.madds += len(Y) * self.rank * len(self.free)
        return out

    def add(self, Y: np.ndarray, chunk: int = 1024) -> int:
        """Insert rows; returns the rank increase."""
        added = 0
        for s in range(0, len(Y), chunk):
            added += self._add_reduced(self.residual(Y[s:s + chunk]))
        return added

    def _add_reduced(self, R: np.ndarray) -> int:
        p = self.p
        nz = np.flatnonzero(R.any(axis=1))
        if not len(nz):
            return 0
        out = []
        t, q, _ = elim.fp_echelon(R[nz], p, out)
        if not t:
            return 0
        E = out[0][:t] % p                     # row echelon, unit pivots at free positions q (ascending)
        q = np.asarray(q[:t], dtype=np.int64)
        E = rref_from_echelon(E, q, p)
        # clear the new pivot columns from the old rows
        if self.rank:
            for s0 in range(0, self.rank, 2048):
                blk = self.B[s0:s0 + 2048]
                coef = blk[:, q].astype(np.int64)
                rows = np.flatnonzero(coef.any(axis=1))
                if not len(rows):
                    continue
                upd = (blk[rows].astype(np.int64) - mulmod(coef[rows], E, p)) % p
                blk[rows] = upd.astype(self.dtype)
                self.madds += len(rows) * t * len(self.free)
        keep = np.ones(len(self.free), dtype=bool)
        keep[q] = False
        # one new array (old and new coexist once): rows r x (f - t) of the old basis, then the t new rows
        newB = np.empty((self.rank + t, int(keep.sum())), dtype=self.dtype)
        r = self.rank
        for s0 in range(0, r, 4096):
            e0 = min(r, s0 + 4096)
            newB[s0:e0] = self.B[s0:e0][:, keep]
        newB[self.rank:] = E[:, keep].astype(self.dtype)
        self.B = newB
        self.piv = np.concatenate([self.piv, self.free[q]])
        self.free = self.free[keep]
        return t

    def rows(self, which: np.ndarray) -> np.ndarray:
        """Full coordinates of the selected rows."""
        out = np.zeros((len(which), self.c), dtype=np.int64)
        out[np.arange(len(which)), self.piv[which]] = 1
        out[:, self.free] = self.B[which]
        return out

    def check_echelon(self) -> bool:
        """Structural check: distinct pivots, free/pivot partition, B zero left of each pivot."""
        if len(set(self.piv.tolist())) != self.rank:
            return False
        if sorted(self.piv.tolist() + self.free.tolist()) != list(range(self.c)):
            return False
        for s in range(0, self.rank, 4096):
            left = self.free[None, :] < self.piv[s:s + 4096, None]
            if (self.B[s:s + 4096].astype(np.int64) * left).any():
                return False
        return True


def rref_from_echelon(E: np.ndarray, q: np.ndarray, p: int, block: int = 128) -> np.ndarray:
    """Back-substitute a row echelon form with unit pivots at ascending columns q to RREF."""
    E = E.astype(np.int64) % p
    t = len(q)
    for e in range(t, 0, -block):
        s = max(0, e - block)
        if e < t:
            C = E[s:e][:, q[e:]]
            if C.any():
                E[s:e] = (E[s:e] - mulmod(C, E[e:], p)) % p
        for i in range(e - 2, s - 1, -1):
            c = E[i, q[i + 1:e]]
            if c.any():
                E[i] = (E[i] - mulmod(c[None, :], E[i + 1:e], p)[0]) % p
    return E


# ---------------------------------------------------------------------------
# the split closure
# ---------------------------------------------------------------------------


class SplitClosure:
    """Resumable split closure of one generator; see the module docstring."""

    def __init__(self, layout: Layout, g: dict[int, int], p: int, target: int, ckpt: Path | None = None,
                 log=None, save_every: float = 0.0):
        assert all(g.get((m >> layout.k) | ((m & ((1 << layout.k) - 1)) << layout.k), 0) == c
                   for m, c in g.items()), 'generator is not symmetric under the summand swap'
        assert p > 2
        self.L, self.g, self.p, self.target = layout, g, p, target
        self.ckpt = ckpt
        self.save_every = save_every  # seconds between checkpoints (pauses and phase ends always save)
        self.log = log or (lambda **kw: None)
        self.parts = [FreeRREF(layout.part_cols[0], p), FreeRREF(layout.part_cols[1], p)]
        deg = max(m.bit_count() for m in g)
        self.mus = [m for m in layout.masks if m.bit_count() + deg <= layout.D][::-1]
        self.state = dict(phase='generators', gen_done=0, rounds=0, sweeps=0, sweep_done=0, sweep_added=0,
                          processed=[[], []], history=[], seconds=0.0,
                          work=dict(generators=[0, 0.0], closure=[0, 0.0], verify=[0, 0.0]), peak_free_entries=0)
        self.started = time.time()

    def _madds(self) -> int:
        return self.parts[0].madds + self.parts[1].madds

    def _account(self, phase: str, t0: float, m0: int):
        w = self.state.setdefault('work', dict(generators=[0, 0.0], closure=[0, 0.0], verify=[0, 0.0]))
        w[phase][0] += self._madds() - m0
        w[phase][1] += time.time() - t0
        self.state['peak_free_entries'] = max(self.state.get('peak_free_entries', 0),
                                              int(sum(P.B.size for P in self.parts)))

    # -- bookkeeping --------------------------------------------------------
    @property
    def rank(self) -> int:
        return self.parts[0].rank + self.parts[1].rank

    def contains_one(self) -> bool:
        return self.L.one_sym in set(self.parts[0].piv.tolist())

    def low_rows(self, part: int) -> np.ndarray:
        P = self.parts[part]
        return np.flatnonzero(self.L.part_weights[part][P.piv] < self.L.D)

    def low_dimension(self) -> int:
        return len(self.low_rows(0)) + len(self.low_rows(1))

    def profile(self) -> list[int]:
        """dim(V cap M_<=w) for w = 0..D."""
        out = []
        for w in range(self.L.D + 1):
            out.append(int(sum((self.L.part_weights[i][self.parts[i].piv] <= w).sum() for i in (0, 1))))
        return out

    def part_profile(self, part: int) -> list[int]:
        pw = self.L.part_weights[part][self.parts[part].piv]
        return [int((pw == w).sum()) for w in range(self.L.D + 1)]

    def save(self, force: bool = False):
        if self.ckpt is None:
            return
        if not force and time.time() - getattr(self, '_last_save', 0.0) < self.save_every:
            return
        self._last_save = time.time()
        self.ckpt.mkdir(parents=True, exist_ok=True)
        arrays = {}
        for i, P in enumerate(self.parts):
            arrays[f'piv{i}'], arrays[f'free{i}'], arrays[f'B{i}'] = P.piv, P.free, P.B
        tmp = self.ckpt / 'parts.tmp.npz'
        with open(tmp, 'wb') as fh:
            np.savez(fh, **arrays)
        os.replace(tmp, self.ckpt / 'parts.npz')
        st = dict(self.state, seconds=self.state['seconds'] + time.time() - self.started)
        (self.ckpt / 'state.tmp.json').write_text(json.dumps(st))
        os.replace(self.ckpt / 'state.tmp.json', self.ckpt / 'state.json')

    def load(self) -> bool:
        if self.ckpt is None or not (self.ckpt / 'state.json').exists():
            return False
        self.state = json.loads((self.ckpt / 'state.json').read_text())
        with np.load(self.ckpt / 'parts.npz') as z:
            for i, P in enumerate(self.parts):
                P.piv, P.free, P.B = z[f'piv{i}'].astype(np.int64), z[f'free{i}'].astype(np.int64), z[f'B{i}']
        self.started = time.time()
        return True

    # -- row sources --------------------------------------------------------
    def _add_full(self, full: np.ndarray) -> int:
        s, a = self.L.project(full, self.p)
        return self.parts[0].add(s) + self.parts[1].add(a)

    def _shift_and_add(self, part: int, which: np.ndarray, chunk: int) -> int:
        """Add the k block-1 shifts of rows `which` of `part` (each projected to both parts)."""
        added = 0
        L = self.L
        for s0 in range(0, len(which), chunk):
            coords = self.parts[part].rows(which[s0:s0 + chunk])
            full = L.lift(part, coords, self.p)
            assert not full[:, L.weights >= L.D].any(), 'a low row has support of weight >= D'
            full_low = full[:, L.low]
            del full
            for j in range(L.k):
                added += self._add_full(L.shifted_block1(full_low, j))
        return added

    # -- driver ---------------------------------------------------------------
    def run(self, batch: int = 600, chunk: int = 256, deadline: float | None = None,
            early_stop: bool = True) -> str:
        """Advance the closure.  Returns 'reached', 'closed', or 'paused' (deadline; checkpoint saved)."""
        st = self.state
        L, p = self.L, self.p

        def out_of_time():
            return deadline is not None and time.time() > deadline

        if st['phase'] == 'generators':
            t0, m0 = time.time(), self._madds()
            gchunk = max(64, 20_000_000 // L.cols)
            while st['gen_done'] < len(self.mus):
                mus = self.mus[st['gen_done']:st['gen_done'] + gchunk]
                self._add_full(L.generator_rows(self.g, p, mus))
                st['gen_done'] += len(mus)
            st['history'].append(self.rank)
            st['phase'] = 'closure'
            self._account('generators', t0, m0)
            self.log(stage='generators', rows=len(self.mus), rank=self.rank)
            self.save()
        while True:
            if early_stop and self.rank >= self.target:
                st['phase'] = 'reached'
                self.save(force=True)
                return 'reached'
            if st['phase'] == 'closure':
                pending = []
                for part in (0, 1):
                    done = set(st['processed'][part])
                    low = self.low_rows(part)
                    piv = self.parts[part].piv[low]
                    pending += [(part, int(i), int(c)) for i, c in zip(low, piv) if int(c) not in done]
                if not pending:
                    st['phase'] = 'verify'
                    st['sweep_done'] = 0
                    st['sweep_added'] = 0
                    st['sweeps'] += 1
                    self.save()
                    continue
                if out_of_time():
                    self.save(force=True)
                    return 'paused'
                chosen = pending[:batch]
                t = time.time()
                m0 = self._madds()
                added = 0
                for part in (0, 1):
                    rows = np.asarray([i for (pp, i, c) in chosen if pp == part], dtype=np.int64)
                    if len(rows):
                        cols = self.parts[part].piv[rows].tolist()
                        added += self._shift_and_add(part, rows, chunk)
                        st['processed'][part] = sorted(set(st['processed'][part]) | set(int(c) for c in cols))
                st['rounds'] += 1
                st['history'].append(self.rank)
                self._account('closure', t, m0)
                self.log(stage='round', round=st['rounds'], processed=len(chosen), pending=len(pending) - len(chosen),
                         added=added, rank=self.rank, ranks=[self.parts[0].rank, self.parts[1].rank],
                         target=self.target, seconds=round(time.time() - t, 1))
                self.save()
                continue
            if st['phase'] == 'verify':
                # snapshot of the low rows: (part, row index) in a fixed order; rows only change if the
                # sweep adds something, in which case the sweep is restarted after a new closure phase
                order = [(part, int(i)) for part in (0, 1) for i in self.low_rows(part)]
                total = len(order)
                while st['sweep_done'] < total:
                    if out_of_time():
                        self.save(force=True)
                        return 'paused'
                    seg = order[st['sweep_done']:st['sweep_done'] + batch]
                    t = time.time()
                    m0 = self._madds()
                    before = self.rank
                    for part in (0, 1):
                        rows = np.asarray([i for (pp, i) in seg if pp == part], dtype=np.int64)
                        if len(rows):
                            st['sweep_added'] += self._shift_and_add(part, rows, chunk)
                    st['sweep_done'] += len(seg)
                    self._account('verify', t, m0)
                    self.log(stage='verify', sweep=st['sweeps'], done=st['sweep_done'], of=total,
                             added=self.rank - before, rank=self.rank, seconds=round(time.time() - t, 1))
                    if st['sweep_added']:
                        # not closed: the low rows changed; return to the closure phase
                        st['phase'] = 'closure'
                        st['history'].append(self.rank)
                        self.save()
                        break
                    self.save()
                else:
                    st['phase'] = 'closed'
                    st['closed_low_rows'] = total
                    self.save(force=True)
                    return 'closed'
                continue
            if st['phase'] in ('closed', 'reached'):
                return st['phase']
            raise ValueError(st['phase'])


def split_closure(g: dict[int, int], k: int, p: int, D: int, target: int, *, early_stop: bool = True,
                  batch: int = 600, ckpt: Path | None = None, log=None) -> dict:
    """One-call form for tests and small depths."""
    L = Layout(k, D)
    sc = SplitClosure(L, g, p, target, ckpt=ckpt, log=log)
    sc.load()
    status = sc.run(batch=batch, early_stop=early_stop)
    return dict(status=status, rank=sc.rank, ranks=[sc.parts[0].rank, sc.parts[1].rank],
                low_dimension=sc.low_dimension(), contains_one=sc.contains_one(), profile=sc.profile(),
                sym_profile=sc.part_profile(0), anti_profile=sc.part_profile(1),
                rounds=sc.state['rounds'], sweeps=sc.state['sweeps'], columns=L.cols,
                part_columns=list(L.part_cols), echelon_ok=all(P.check_echelon() for P in sc.parts))
