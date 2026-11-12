"""Symmetry-split degree-D closure for two-summand chain systems: unit AND proper ideals,
with a deterministic closure verification and a memory-lean input path.

Extends csd/ic_chain_block_closure.py (which decides only "1 in V_D" for a unit ideal) to the
certificate needed for proper ideals, and makes the closure exact when it stalls.

Setting.  Two summands, depth k, normalized generator g = NF(S_3(c x1, c x2, x_R)) modulo the
chains (squarefree masks over 2k chain variables, multiplication by a chain variable = exponent
shift mod 2^k, csd/ic_chain_shift_closure.py).  For D >= 2 the closure modulo the chains equals
the ordinary closure V_D (module docstring of csd/ic_fall_certificates.py).  The swap sigma of the
two summands fixes g whenever both summands use the same scale c, so V_D is sigma-invariant and
V_D = V^+ (+) V^- (char p != 2); each part is kept as a reduced row echelon basis (uint16 for
p < 2^16, uint32 above) over its own coordinates, via csd.ic_chain_block_closure.EchelonPart.

Rows added (each lies in V_D):
  * NF(mu * g) for squarefree mu with wt(mu) + deg g <= D (all inputs of degree <= D);
  * optionally the rows of a basis of V_{D-1} (subset of V_D), inserted in small chunks so the
    full-width int64 embedding is never materialised;
  * random combinations sum_v shift_v(R_v L) of the 2k shifts of low rows L (pivot weight < D).
Rows are projected to the two parts (both projections lie in V_D since V_D is sigma-invariant).

Certificates (dim W := rank(sym part) + rank(anti part) = dim(W' + sigma W') <= dim V_D):
  * reached:  dim W = target := |M_<=D| - rank(ev on M_<=D at Z) = dim(I cap M_<=D).  Since
              V_D subset I cap M_<=D, then V_D = I cap M_<=D.  With ev(M_<D) onto Z the
              local-to-global theorem (experiments/ic-curve-validation/code/lfd_engine.py
              docstring) gives LFD <= D.  Unit ideal: target = |M_<=D|, equivalently 1 in V_D.
  * closed:   every one of the 2k shifts of every low row of both parts re-reduces to zero
              against the final bases (one clean deterministic pass).  W contains all inputs
              and is closed under every admissible product, so W contains V_D; with W subset V_D,
              W = V_D exactly.  If moreover dim W < target, some h in I of degree <= D is not in
              V_D, so LFD >= D + 1 (unit ideal: 1 not in V_D).
Random combinations only affect speed: a stall is never reported as a bound unless the clean
verification pass has run.
"""
from __future__ import annotations

import math
import time

import numpy as np

from . import ic_chain_block_closure as B
from . import ic_chain_shift_closure as SC
from .ic_curve_sweep import prime_evaluation_rank


class LeanEchelonPart(B.EchelonPart):
    """csd.ic_chain_block_closure.EchelonPart with a configurable float64 tile height (the default
    4096-row float64 basis tile alone is 1 GB at depth 10, degree 6).  Same arithmetic."""

    def __init__(self, cols: int, p: int, tile: int = 2048):
        super().__init__(cols, p)
        self.tile = tile

    def reduce(self, rows: np.ndarray) -> np.ndarray:
        p = self.p
        rows = np.asarray(rows, dtype=np.int64) % p
        if self.rank == 0 or len(rows) == 0:
            return rows
        coeff = rows[:, self.pivots]
        red = (rows.astype(np.float64) - B._modmat(coeff, self.basis, p, self.max_inner, self.tile)) % p
        return red.astype(np.int64)

    def add(self, rows: np.ndarray, chunk: int = 1000) -> int:
        p = self.p
        added = 0
        for s in range(0, len(rows), chunk):
            red = self.reduce(rows[s:s + chunk])
            nz = np.flatnonzero(red.any(axis=1))
            if not len(nz):
                continue
            out = []
            r, piv, _ = B.elim.fp_echelon(red[nz], p, out)
            if not r:
                continue
            new = out[0][:r] % p
            piv = np.asarray(piv[:r], dtype=np.int64)
            newf = new.astype(np.float64)
            for i in range(r - 2, -1, -1):
                c = newf[i, piv[i + 1:]]
                if c.any():
                    newf[i] = (newf[i] - B._modmat(c[None, :], newf[i + 1:], p, self.max_inner, self.tile)[0]) % p
            new = newf.astype(np.int64)
            del newf
            if self.rank:
                for s0 in range(0, self.rank, self.tile):
                    blk = self.basis[s0:s0 + self.tile]
                    c = blk[:, piv].astype(np.int64)
                    if not c.any():
                        continue
                    upd = (blk.astype(np.float64) - B._modmat(c, new, p, self.max_inner, self.tile)) % p
                    self.basis[s0:s0 + self.tile] = upd.astype(self.basis.dtype)
            self.basis = np.vstack([self.basis, new.astype(self.basis.dtype)])
            self.pivots = np.concatenate([self.pivots, piv])
            added += r
        return added


def swap_mask(mask: int, k: int) -> int:
    N = 1 << k
    return (mask >> k) | ((mask & (N - 1)) << k)


def is_symmetric(generator: dict[int, int], k: int, p: int) -> bool:
    return all((generator.get(swap_mask(m, k), 0) - c) % p == 0 for m, c in generator.items())


def monomial_count(k: int, D: int) -> int:
    return sum(math.comb(2 * k, d) for d in range(min(2 * k, D) + 1))


def fixed_count(k: int, D: int) -> int:
    """Swap-fixed squarefree masks (a1 = a2) of weight <= D."""
    return sum(math.comb(k, w) for w in range(D // 2 + 1))


class SplitClosure:
    """Degree-D closure of {g} + chains, split into sigma-symmetric and antisymmetric parts."""

    def __init__(self, generator: dict[int, int], k: int, p: int, D: int, *, points=(), generator_value=None,
                 eval_max_cells: int = 400_000_000, tile: int = 2048):
        if D < 2:
            raise ValueError('The chain-quotient closure equals the ordinary closure only for D >= 2.')
        if p >= 1 << 25:
            raise ValueError('exact float64 tiles need p^2 * 2 < 2^53 per term; use p < 2^25.')
        if not is_symmetric(generator, k, p):
            raise ValueError('The generator is not swap-symmetric; the split does not apply.')
        self.generator = {int(m): int(c) % p for m, c in generator.items() if c % p}
        self.k, self.p, self.D = k, p, D
        self.tile = tile
        self.masks = SC.monomials(k, D)
        self.index = {m: i for i, m in enumerate(self.masks)}
        self.cols = len(self.masks)
        self.perms = [SC.shift_permutation(k, self.masks, self.index, v) for v in range(2 * k)]
        self.weights = np.asarray([m.bit_count() for m in self.masks])
        self.split = B.SymmetricSplit(k, self.masks, self.index)
        self.one = self.cols - 1
        assert self.masks[self.one] == 0
        self.one_sym = int(np.flatnonzero(self.split.sym_reps == self.one)[0])
        self.sym = LeanEchelonPart(len(self.split.sym_reps), p, tile)
        self.anti = LeanEchelonPart(len(self.split.anti_reps), p, tile)
        self.lowcols = np.flatnonzero(self.weights < D)
        # weight of every shifted low monomial is <= D, so every low column has an image
        for perm in self.perms:
            assert np.all(perm[self.lowcols] >= 0)
        self.points = list(points)
        self.generator_value = generator_value
        if self.points:
            assert generator_value is not None
            N = 1 << k
            for bits in self.points:
                E1, E2 = bits & (N - 1), bits >> k
                val = sum(c * pow(generator_value, ((m & (N - 1)) * E1 + (m >> k) * E2) % N, p)
                          for m, c in self.generator.items()) % p
                assert val == 0, 'a supplied point does not vanish on the generator'
            self.eval_rank_le = prime_evaluation_rank(self.points, self.masks, k, generator_value, p, eval_max_cells)
            self.eval_rank_lt = prime_evaluation_rank(self.points, [m for m in self.masks if m.bit_count() < D], k,
                                                      generator_value, p, eval_max_cells)
        else:
            self.eval_rank_le = self.eval_rank_lt = 0
        self.target = self.cols - self.eval_rank_le
        self.onto_lt = self.eval_rank_lt == len(self.points)
        self.rounds = 0
        self.passes = []
        self.history = []

    # ------------------------------------------------------------------ bookkeeping
    @property
    def rank(self) -> int:
        return self.sym.rank + self.anti.rank

    @property
    def contains_one(self) -> bool:
        return bool(self.one_sym in set(self.sym.pivots.tolist()))

    @property
    def low_dimension(self) -> int:
        """dim(W cap M_<D) (pivots of weight < D in both parts; columns sorted by weight descending)."""
        ws = self.weights[self.split.sym_reps][self.sym.pivots] if self.sym.rank else np.zeros(0)
        wa = self.weights[self.split.anti_reps][self.anti.pivots] if self.anti.rank else np.zeros(0)
        return int((ws < self.D).sum() + (wa < self.D).sum())

    @property
    def reached(self) -> bool:
        return self.rank >= self.target or self.contains_one

    def add_full(self, rows: np.ndarray, chunk: int = 1000) -> tuple[int, int]:
        """Insert full-coordinate rows (each must lie in V_D); returns rank increases (sym, anti)."""
        p = self.p
        da = db = 0
        for s in range(0, len(rows), chunk):
            blk = np.asarray(rows[s:s + chunk], dtype=np.int64) % p
            da += self.sym.add(self.split.sym(blk, p), chunk=chunk)
            db += self.anti.add(self.split.anti(blk, p), chunk=chunk)
        return da, db

    def add_generator_rows(self) -> tuple[int, int]:
        rows = SC.generator_rows(self.generator, self.k, self.D, self.masks, self.index, self.p)
        return self.add_full(rows)

    def add_previous(self, prev_rows: np.ndarray, prev_masks, chunk: int = 500) -> tuple[int, int]:
        """Insert a basis of V_{D-1} (any integer dtype) chunk by chunk without a full int64 copy."""
        pidx = np.asarray([self.index[m] for m in prev_masks])
        da = db = 0
        for s in range(0, len(prev_rows), chunk):
            blk = np.zeros((min(chunk, len(prev_rows) - s), self.cols), dtype=np.int64)
            blk[:, pidx] = np.asarray(prev_rows[s:s + chunk], dtype=np.int64)
            a, b = self.add_full(blk, chunk=chunk)
            da += a
            db += b
        return da, db

    def _low_rows(self, part, reps, partners, sign, batch):
        """Yield the low rows of one part as full-coordinate rows restricted to self.lowcols."""
        if not part.rank:
            return
        p = self.p
        inv2 = pow(2, -1, p)
        lowpos = -np.ones(self.cols, dtype=np.int64)
        lowpos[self.lowcols] = np.arange(len(self.lowcols))
        rep_w = self.weights[reps]
        low_idx = np.flatnonzero(rep_w[part.pivots] < self.D)
        low_part_cols = np.flatnonzero(rep_w < self.D)
        tgt_a = lowpos[reps[low_part_cols]]
        tgt_b = lowpos[partners[low_part_cols]]
        for s0 in range(0, len(low_idx), batch):
            blk = part.basis[low_idx[s0:s0 + batch]][:, low_part_cols].astype(np.int64) * inv2 % p
            out = np.zeros((blk.shape[0], len(self.lowcols)), dtype=np.int64)
            out[:, tgt_a] = blk
            out[:, tgt_b] = (sign * blk) % p  # diagonal (sym only): tgt_a == tgt_b and sign = +1
            yield out

    def _all_low_rows(self, batch):
        sp = self.split
        yield from self._low_rows(self.sym, sp.sym_reps, sp.partner_sym, 1, batch)
        yield from self._low_rows(self.anti, sp.anti_reps, sp.partner_anti, -1, batch)

    # ------------------------------------------------------------------ closure steps
    def random_round(self, rng, combos: int, low_batch: int = 1000, terms: int | None = None) -> tuple[int, int]:
        """Add q random elements of span{shift_v(L)} (L = low rows).  terms=None: dense random
        coefficients (cost q * n_low * |M_<D| * 2k); terms=t: each combination draws t random low
        rows per batch and per shift (cost q * t * |M_<D| * 2k * batches).  Either way every added row
        lies in V_D; the choice only changes how fast the rank grows."""
        p = self.p
        remaining = (self.sym.cols - self.sym.rank) + (self.anti.cols - self.anti.rank)
        q = int(min(combos, remaining + 64))
        if terms is None:
            Y = np.zeros((q, self.cols), dtype=np.float64)
            for Lb in self._all_low_rows(low_batch):
                Lb = Lb.astype(np.float64)
                for perm in self.perms:
                    R = rng.integers(0, p, size=(q, len(Lb))).astype(np.float64)
                    Y[:, perm[self.lowcols]] += B._modmat(R, Lb, p, self.sym.max_inner, self.tile)
                    np.fmod(Y, p, out=Y)
            return self.add_full(Y.astype(np.int64))
        Y = np.zeros((q, self.cols), dtype=np.int64)
        for Lb in self._all_low_rows(low_batch):
            nb = len(Lb)
            t = min(terms, nb)
            for perm in self.perms:
                idx = rng.integers(0, nb, size=(q, t))
                co = rng.integers(1, p, size=(q, t))
                acc = np.zeros((q, len(self.lowcols)), dtype=np.int64)
                for j in range(t):
                    acc += co[:, j:j + 1] * Lb[idx[:, j]]
                    acc %= p
                tgt = perm[self.lowcols]
                Y[:, tgt] = (Y[:, tgt] + acc) % p
        return self.add_full(Y)

    def basis_full_rows(self, batch: int = 1000):
        """Yield a basis of W = V^+ (+) V^- as full-coordinate int64 rows (sym rows, then anti rows)."""
        p, sp = self.p, self.split
        for part, conv in ((self.sym, sp.full_from_sym), (self.anti, sp.full_from_anti)):
            for s0 in range(0, part.rank, batch):
                yield conv(part.basis[s0:s0 + batch].astype(np.int64), self.cols, p)

    def verify_pass(self, low_batch: int = 500) -> int:
        """One deterministic pass: reduce every shift of every low row against both parts.

        Returns the number of nonzero residual rows; residuals are added (they lie in V_D) and the
        pass stops at the first batch that has any, so a return value of 0 means the bases were
        unchanged during the whole pass and the space is closed."""
        p = self.p
        sp = self.split
        for Lb in self._all_low_rows(low_batch):
            res = []
            for perm in self.perms:
                S = np.zeros((len(Lb), self.cols), dtype=np.int64)
                S[:, perm[self.lowcols]] = Lb
                rs = self.sym.reduce(sp.sym(S, p))
                ra = self.anti.reduce(sp.anti(S, p))
                bad = np.flatnonzero(rs.any(axis=1) | ra.any(axis=1))
                if len(bad):
                    res.append(S[bad])
            if res:
                R = np.vstack(res)
                self.add_full(R)
                return len(R)
        return 0

    def run(self, *, seed: int = 0, combos: int = 2000, max_rounds: int = 200, verify: bool = True,
            max_passes: int = 50, low_batch: int = 1000, verify_batch: int = 500, terms: int | None = None,
            log=None) -> dict:
        rng = np.random.default_rng(seed)
        started = time.time()
        closed = False
        while True:
            # random rounds until a round adds nothing (or the target is reached)
            stalled = False
            while not self.reached and self.rounds < max_rounds:
                self.rounds += 1
                t = time.time()
                before = self.rank
                da, db = self.random_round(rng, combos, low_batch, terms)
                self.history.append(self.rank)
                if log:
                    log(dict(stage='round', round=self.rounds, rank=self.rank, target=self.target, added=[da, db],
                             sym_rank=self.sym.rank, anti_rank=self.anti.rank, seconds=round(time.time() - t, 1)))
                if self.rank == before:
                    stalled = True
                    break
            if self.reached or not verify or len(self.passes) >= max_passes:
                break
            if not stalled and self.rounds >= max_rounds:
                break  # round budget exhausted while still growing: undecided, resumable
            t = time.time()
            missing = self.verify_pass(verify_batch)
            self.passes.append(missing)
            if log:
                log(dict(stage='verify', passes=len(self.passes), missing=missing, rank=self.rank,
                         seconds=round(time.time() - t, 1)))
            if missing == 0:
                closed = True
                break
        return self.result(closed=closed, seconds=time.time() - started)

    def result(self, closed: bool, seconds: float = 0.0) -> dict:
        reached = self.reached
        unit = not self.points
        out = dict(k=self.k, prime=self.p, degree=self.D, cols=self.cols, sym_cols=self.sym.cols,
                   anti_cols=self.anti.cols, points=len(self.points), eval_rank_le=self.eval_rank_le,
                   eval_rank_lt=self.eval_rank_lt, evaluation_onto_lt=self.onto_lt, target=self.target,
                   rank=self.rank if not self.contains_one else self.cols, sym_rank=self.sym.rank,
                   anti_rank=self.anti.rank, contains_one=self.contains_one,
                   low_dimension=self.low_dimension if not self.contains_one else monomial_count(self.k, self.D - 1),
                   reached=bool(reached), closed=bool(closed), rounds=self.rounds, verify_passes=list(self.passes),
                   rank_history=list(self.history), seconds=round(seconds, 2))
        exact_dim = reached or closed
        out['exact_dimension'] = (out['rank'] if exact_dim else None)
        if reached and (unit or self.onto_lt):
            out['status'] = 'reached'
            out['lfd_upper_bound'] = self.D
            out['certificate'] = ('1 in V_D' if unit else
                                  'dim V_D = |M_<=D| - |Z| and ev(M_<D) onto Z: V_D = I cap M_<=D, LFD <= D')
        elif reached:
            out['status'] = 'reached-not-onto'
            out['certificate'] = 'V_D = I cap M_<=D but ev(M_<D) is not onto Z: no LFD bound from this degree'
        elif closed:
            out['status'] = 'closed-deficient'
            out['lfd_lower_bound'] = self.D + 1
            out['certificate'] = ('closure verified closed and 1 not in V_D' if unit else
                                  'closure verified closed and dim V_D < dim(I cap M_<=D)')
        else:
            out['status'] = 'undecided'
            out['certificate'] = 'stalled without a clean verification pass: no bound'
        return out


def split_degree(generator, k, p, D, *, points=(), generator_value=None, prev=None, seed=0, combos=2000,
                 max_rounds=200, verify=True, log=None, tile=2048, terms=None, keep=False):
    """One degree: build the split closure (from the generator rows, plus prev=(rows, masks) if given).
    Returns the result dict, or (result, SplitClosure) when keep is set."""
    sc = SplitClosure(generator, k, p, D, points=points, generator_value=generator_value, tile=tile)
    sc.add_generator_rows()
    if prev is not None:
        sc.add_previous(*prev)
    res = sc.run(seed=seed, combos=combos, max_rounds=max_rounds, verify=verify, log=log, terms=terms)
    return (res, sc) if keep else res


def split_profile(generator, k, p, points, generator_value, *, max_degree=6, seed=0, combos=2000, terms=None,
                  log=None) -> dict:
    """Exact LFD from split closures at D = 2..max_degree (each degree from the generator rows).

    Falls are taken at D >= 3 only when both dim V_{D-1} and dim(V_D cap M_<D) are exact (each
    degree reached or verified closed).  Unit ideal: exact LFD = first D with 1 in V_D, given every
    earlier degree verified closed without 1.  Proper ideal: first D reaching the target with ev(M_<D)
    onto Z; exact LFD = last strict fall (all degrees before it exact)."""
    per = []
    falls = []
    prev_dim = None
    all_exact = True
    for D in range(2, max_degree + 1):
        r = split_degree(generator, k, p, D, points=points, generator_value=generator_value, seed=seed + D,
                         combos=combos, terms=terms, log=log)
        per.append(r)
        exact = r['exact_dimension'] is not None
        if D >= 3 and prev_dim is not None and exact and r['low_dimension'] > prev_dim:
            falls.append(D)
        if not exact:
            all_exact = False
        if r['status'] == 'reached':
            if not points:
                lfd = D if all(q['status'] == 'closed-deficient' for q in per[:-1]) else None
            else:
                lfd = max(falls, default=0) if all_exact else None
            return dict(status='exact' if lfd is not None else 'upper', exact_lfd=lfd, lfd_upper_bound=D,
                        strict_falls=falls, dims=[q['exact_dimension'] for q in per], per_degree=per)
        if r['status'] != 'closed-deficient' and r['status'] != 'reached-not-onto':
            return dict(status='undecided', exact_lfd=None, strict_falls=falls,
                        dims=[q['exact_dimension'] for q in per], per_degree=per)
        prev_dim = r['exact_dimension']
    lower = max_degree + 1 if per[-1]['status'] == 'closed-deficient' else None
    return dict(status='capped', exact_lfd=None, lfd_lower_bound=lower, strict_falls=falls,
                dims=[q['exact_dimension'] for q in per], per_degree=per)


# ---------------------------------------------------------------------- planning (item (c))
def memory_model(k: int, D: int, *, prime: int, prev_rank: int = 0, prev_cols: int = 0, combos: int = 2000,
                 low_batch: int = 1000, tile: int = 2048, chunk: int = 1000, full_rank: bool = True) -> dict:
    """Peak-memory estimate (bytes) of SplitClosure at (k, D) when both parts approach full rank.

    Terms: the two echelon parts (uint16 for p < 2^16 else uint32), the random-round buffers
    (Y and one combination, float64), one low batch (int64 + float64), the reduction temporaries
    of EchelonPart.add/reduce (rows int64 + float64 + output, and one float64 basis tile), and the
    previous basis held as uint16.  A model, not a measurement; calibrated at small k by
    scripts/ic_bounds_exact10_validate.py (maxrss column)."""
    cols = monomial_count(k, D)
    low = monomial_count(k, D - 1)
    fixed = fixed_count(k, D)
    sym = (cols + fixed) // 2
    anti = cols - sym
    w = 2 if prime < 65536 else 4
    parts = (sym * sym + anti * anti) * w if full_rank else 0
    rnd = combos * cols * 8 + combos * low * 8 + low_batch * low * 16 + combos * low_batch * 8
    red = chunk * sym * 8 * 3 + tile * sym * 8 * 2
    prev = prev_rank * prev_cols * 2
    total = parts + max(rnd, red) + red + prev
    return dict(k=k, D=D, cols=cols, low_cols=low, sym_cols=sym, anti_cols=anti, word_bytes=w,
                parts_bytes=parts, round_bytes=rnd, reduce_bytes=red, prev_bytes=prev, peak_bytes=total)
