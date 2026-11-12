"""Blocked, symmetry-split upper certificate for the degree-D chain closure of large depth.

Decides whether 1 lies in the degree-D closure V_D of a two-summand chain system whose normalized
generator g is symmetric under swapping the two summands, starting from a known basis of V_{D-1}
(for example the verified V_5 of scripts/run_ic_chain_depth10_probe.py).

Representation.  Work modulo the chains in the squarefree monomial space M_<=D; multiplication by a
chain variable is the exponent shift of csd/ic_chain_shift_closure.py.  The swap sigma of the two
summands fixes g, so V_D is sigma-invariant and splits as V_D = V^+ (+) V^- into its symmetric and
antisymmetric parts; each part is stored as a reduced row echelon basis over its own coordinates
(about half of |M_<=D| each), in uint16 for p < 2^16.

Rows.  Every row added lies in V_D:
  * the rows of V_{D-1} (V_{D-1} is contained in V_D);
  * random linear combinations sum_v shift_v(R_v @ L) of the variable shifts of low rows L of the
    current space (weight < D), which is an admissible product of a derived lower-degree row.
Rows are projected to the two parts (both projections lie in V_D because V_D is sigma-invariant)
and reduced against the echelon bases with exact float64 products mod p (inner dimensions are
bounded so every accumulated sum stays below 2^53).

Soundness is one-sided: reaching full rank in the symmetric part proves 1 in V_D (then V_D is all of
M_<=D).  Random combinations only affect completeness: a round that adds less than the true span is
detected by the next round, and a failure to reach full rank is reported as undecided, never as a
lower bound.
"""
from __future__ import annotations

import json
import math
import time

import numpy as np

from . import elim
from .ic_chain_shift_closure import monomials, shift_permutation

TWO53 = 2 ** 53


def _modmat(a: np.ndarray, b: np.ndarray, p: int, max_inner: int, tile: int = 2048) -> np.ndarray:
    """Exact (a @ b) mod p for nonnegative residues, splitting the inner dimension into tiles small
    enough that every partial sum stays below 2^53 and every float64 tile fits in memory."""
    n = a.shape[1]
    step = max(1, min(max_inner, tile))
    out = np.zeros((a.shape[0], b.shape[1]), dtype=np.float64)
    for s in range(0, n, step):
        e = min(n, s + step)
        out += a[:, s:e].astype(np.float64) @ b[s:e].astype(np.float64)
        np.fmod(out, p, out=out)
    return out


class EchelonPart:
    """Reduced row echelon basis (pivots normalised to 1, zero above and below) over F_p."""

    def __init__(self, cols: int, p: int):
        self.cols = cols
        self.p = p
        self.basis = np.zeros((0, cols), dtype=np.uint16 if p < 65536 else np.uint32)
        self.pivots = np.zeros(0, dtype=np.int64)
        self.max_inner = max(1, int((TWO53 // (p * p)) // 2))

    @property
    def rank(self):
        return len(self.pivots)

    def reduce(self, rows: np.ndarray) -> np.ndarray:
        """rows - rows[:, pivots] @ basis (mod p); the result is zero on every pivot column."""
        p = self.p
        rows = np.asarray(rows, dtype=np.int64) % p
        if self.rank == 0 or len(rows) == 0:
            return rows
        coeff = rows[:, self.pivots]
        red = (rows.astype(np.float64) - _modmat(coeff, self.basis, p, self.max_inner)) % p
        return red.astype(np.int64)

    def add(self, rows: np.ndarray, chunk: int = 2000) -> int:
        """Insert rows; returns the rank increase."""
        p = self.p
        added = 0
        for s in range(0, len(rows), chunk):
            red = self.reduce(rows[s:s + chunk])
            nz = np.flatnonzero(red.any(axis=1))
            if not len(nz):
                continue
            out = []
            r, piv, _ = elim.fp_echelon(red[nz], p, out)
            if not r:
                continue
            new = out[0][:r] % p
            piv = np.asarray(piv[:r], dtype=np.int64)
            # reduced form of the new rows among themselves: fp_echelon returns a row echelon form
            # with unit pivots but nonzero entries above them; back-substitute bottom-up, one
            # vectorised step per row, keeping every intermediate sum exact.
            newf = new.astype(np.float64)
            for i in range(r - 2, -1, -1):
                c = newf[i, piv[i + 1:]]
                if c.any():
                    newf[i] = (newf[i] - _modmat(c[None, :], newf[i + 1:], p, self.max_inner)[0]) % p
            new = newf.astype(np.int64)
            # clear the new pivot columns from the old basis
            if self.rank:
                for s0 in range(0, self.rank, 2048):
                    blk = self.basis[s0:s0 + 2048]
                    c = blk[:, piv].astype(np.int64)
                    if not c.any():
                        continue
                    upd = (blk.astype(np.float64) - _modmat(c, new, p, self.max_inner)) % p
                    self.basis[s0:s0 + 2048] = upd.astype(self.basis.dtype)
            self.basis = np.vstack([self.basis, new.astype(self.basis.dtype)])
            self.pivots = np.concatenate([self.pivots, piv])
            added += r
        return added


class SymmetricSplit:
    """Coordinates of the sigma-symmetric and antisymmetric parts of M_<=D."""

    def __init__(self, k: int, masks: list[int], index: dict[int, int]):
        N = 1 << k
        swap = np.asarray([index[(m >> k) | ((m & (N - 1)) << k)] for m in masks])
        self.swap = swap
        reps = [i for i in range(len(masks)) if i <= swap[i]]
        self.sym_reps = np.asarray(reps)
        self.anti_reps = np.asarray([i for i in reps if swap[i] != i])
        self.partner_sym = swap[self.sym_reps]
        self.partner_anti = swap[self.anti_reps]

    def sym(self, rows: np.ndarray, p: int) -> np.ndarray:
        """Coordinates of (f + sigma f): f_m + f_{sigma m} on orbit representatives (diagonal: 2 f_m)."""
        return (rows[:, self.sym_reps] + rows[:, self.partner_sym]) % p

    def anti(self, rows: np.ndarray, p: int) -> np.ndarray:
        return (rows[:, self.anti_reps] - rows[:, self.partner_anti]) % p

    def full_from_sym(self, coords: np.ndarray, cols: int, p: int) -> np.ndarray:
        """A symmetric full vector with the given sym coordinates (scaled by 1/2 on pairs is unnecessary
        for spans; diagonal entries are halved so that sym(full) == coords)."""
        inv2 = pow(2, -1, p)
        out = np.zeros((coords.shape[0], cols), dtype=np.int64)
        c = coords.astype(np.int64) * inv2 % p
        out[:, self.sym_reps] = c
        out[:, self.partner_sym] = c
        return out

    def full_from_anti(self, coords: np.ndarray, cols: int, p: int) -> np.ndarray:
        inv2 = pow(2, -1, p)
        out = np.zeros((coords.shape[0], cols), dtype=np.int64)
        c = coords.astype(np.int64) * inv2 % p
        out[:, self.anti_reps] = c
        out[:, self.partner_anti] = (-c) % p
        return out


def block_closure_contains_one(k: int, p: int, D: int, prev_rows: np.ndarray, prev_masks: list[int], *,
                               seed: int = 0, combos_per_round: int = 4000, low_batch: int = 3000,
                               max_rounds: int = 30, log=None, state=None):
    """Decide 1 in V_D (upper certificate).  prev_rows: a basis of V_{D-1} over prev_masks.

    Returns a dict with contains_one, ranks per part and rounds.  `state` may carry a resumable
    dictionary (parts and round counter) between calls; the caller checkpoints it.
    """
    rng = np.random.default_rng(seed)
    masks = monomials(k, D)
    index = {m: i for i, m in enumerate(masks)}
    cols = len(masks)
    perms = [shift_permutation(k, masks, index, v) for v in range(2 * k)]
    weights = np.asarray([m.bit_count() for m in masks])
    split = SymmetricSplit(k, masks, index)
    one = cols - 1
    one_sym = int(np.flatnonzero(split.sym_reps == one)[0])
    if state is None:
        state = dict(round=0, sym=EchelonPart(len(split.sym_reps), p), anti=EchelonPart(len(split.anti_reps), p))
        pidx = np.asarray([index[m] for m in prev_masks])
        # embed V_{D-1} in batches: the dense int64 embedding of every row at once is
        # rows x |M_<=D| x 8 bytes (about 9 GB at k = 10, D = 6)
        for s0 in range(0, len(prev_rows), 2000):
            prev = np.zeros((len(prev_rows[s0:s0 + 2000]), cols), dtype=np.int64)
            prev[:, pidx] = prev_rows[s0:s0 + 2000]
            state['sym'].add(split.sym(prev, p))
            state['anti'].add(split.anti(prev, p))
            del prev
    sym, anti = state['sym'], state['anti']
    low_sym = split.sym_reps
    while state['round'] < max_rounds:
        state['round'] += 1
        t = time.time()
        if one_sym in set(sym.pivots.tolist()) or sym.rank == sym.cols:
            break
        lowcols = np.flatnonzero(weights < D)  # a low row is supported here, and x*m stays in M_<=D
        lowpos = -np.ones(cols, dtype=np.int64)
        lowpos[lowcols] = np.arange(len(lowcols))
        inv2 = pow(2, -1, p)

        def low_batches(part, reps, partners, sign):
            """Low rows of one part (pivot of weight < D), mapped to full low-column coordinates."""
            if not part.rank:
                return
            rep_w = weights[reps]
            low_idx = np.flatnonzero(rep_w[part.pivots] < D)
            low_reps_cols = np.flatnonzero(rep_w < D)      # part columns that are low
            tgt_a = lowpos[reps[low_reps_cols]]
            tgt_b = lowpos[partners[low_reps_cols]]
            for s0 in range(0, len(low_idx), low_batch):
                blk = part.basis[low_idx[s0:s0 + low_batch]][:, low_reps_cols].astype(np.int64) * inv2 % p
                out = np.zeros((blk.shape[0], len(lowcols)), dtype=np.int64)
                out[:, tgt_a] = blk
                out[:, tgt_b] = (sign * blk) % p  # diagonal entries: a = b and sign = +1 (sym only)
                yield out

        before = (sym.rank, anti.rank)
        remaining = (sym.cols - sym.rank) + (anti.cols - anti.rank)
        q = min(combos_per_round, remaining + 64)
        Y = np.zeros((q, cols), dtype=np.float64)
        n_low = 0
        for part, reps, partners, sign in ((sym, split.sym_reps, split.partner_sym, 1),
                                           (anti, split.anti_reps, split.partner_anti, -1)):
            for Lb in low_batches(part, reps, partners, sign):
                n_low += len(Lb)
                Lb = Lb.astype(np.float64)
                for perm in perms:
                    Rv = rng.integers(0, p, size=(q, len(Lb))).astype(np.float64)
                    comb = _modmat(Rv, Lb, p, sym.max_inner)
                    Y[:, perm[lowcols]] += comb
                    np.fmod(Y, p, out=Y)
        Y = Y.astype(np.int64)
        Ys = split.sym(Y, p)
        Ya = split.anti(Y, p)
        del Y
        da = sym.add(Ys)
        del Ys
        db = anti.add(Ya)
        del Ya
        if log:
            log(dict(round=state['round'], low_rows=n_low, combos=q, sym_rank=sym.rank, sym_cols=sym.cols,
                     anti_rank=anti.rank, anti_cols=anti.cols, added=[da, db], seconds=round(time.time() - t)))
        if (sym.rank, anti.rank) == before:
            break
    contains_one = one_sym in set(sym.pivots.tolist()) or sym.rank == sym.cols
    return dict(contains_one=bool(contains_one), sym_rank=sym.rank, sym_cols=sym.cols,
                anti_rank=anti.rank, anti_cols=anti.cols, rounds=state['round'], state=state)
