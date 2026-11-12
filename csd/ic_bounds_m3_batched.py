"""Memory-bounded exact closure V_D for chain systems with an optional ordinary auxiliary variable.

Same object as csd.ic_bounds_m3.mixed_closure / csd.ic_fall_certificates.chain_degree_closure
(V_D of {chains} + inputs, represented modulo the chains), computed by a work-queue saturation
with a reduced row echelon basis and exact float64 matrix products, so that depth k = 4 with
m = 3 (|M_<=7| = 3302 chain monomials; 5282 mixed monomials at D = 6) fits in a 3 GB address
space.

Algorithm.  Columns are monomials sorted by degree descending.  Every vector ever added to the
span is reduced against the current basis first, so its leading column is a new pivot and the
added vectors have pairwise distinct leading columns.  Consequently every element of the span of
degree < D is a combination of added vectors of degree < D, and the span is closed under the
admissible products h -> x*h (deg h < D) as soon as x*v lies in the span for every added vector v
of degree < D and every variable x.  The queue holds exactly those products; the run ends when the
queue is empty (closed) or the constant monomial becomes a pivot (1 in V_D, then V_D = M_<=D).
An independent verification pass then re-reduces every variable shift of every final basis row
of degree < D and requires all of them to vanish.

Exactness of float64 products: entries are residues in [0, p); every accumulated product sum is
at most n * (p - 1)^2, required to be below 2^53.
"""
from __future__ import annotations

import numpy as np

from . import elim
from .ic_bounds_m3 import MixedSpace

TWO53 = 2 ** 53


class BatchedClosure:
    def __init__(self, space: MixedSpace, p: int, batch: int = 768):
        n = len(space)
        if n * (p - 1) ** 2 >= TWO53:
            raise ValueError('float64 exactness bound n*(p-1)^2 < 2^53 violated')
        self.sp, self.p, self.n, self.batch = space, p, n, batch
        self.B = np.zeros((0, n), dtype=np.float64)
        self.piv = np.zeros(0, dtype=np.int64)
        self.added = []            # int32 copies of added vectors of degree < D (product sources)
        self.queue = []            # (source index, variable)
        self.contains_one = False

    # -- reduction -----------------------------------------------------------------------------
    def _reduce(self, P):
        p = self.p
        if len(self.piv):
            P = P - P[:, self.piv] @ self.B
            P = np.mod(P, p)
        return P

    def _absorb(self, R):
        """R: float rows already reduced against the basis. Echelonize, insert, enqueue products."""
        p = self.p
        nz = np.flatnonzero(np.any(R != 0, axis=1))
        if not len(nz):
            return 0
        out = []
        s, newpiv, _ = elim.fp_echelon(R[nz].astype(np.int64), p, out)
        if not s:
            return 0
        Enew = out[0].astype(np.int64)
        # fully reduce the new rows among themselves (RREF)
        for i in range(s - 1, 0, -1):
            col = newpiv[i]
            c = Enew[:i, col].copy()
            if np.any(c):
                Enew[:i] = (Enew[:i] - np.outer(c, Enew[i])) % p
        Ef = Enew.astype(np.float64)
        if len(self.piv):
            self.B = np.mod(self.B - self.B[:, newpiv] @ Ef, p)
        self.B = np.vstack([self.B, Ef])
        self.piv = np.concatenate([self.piv, np.asarray(newpiv, dtype=np.int64)])
        D = self.sp.D
        nvar = len(self.sp.cols)
        for i in range(s):
            if newpiv[i] == self.n - 1:
                self.contains_one = True
            if self.sp.deg[newpiv[i]] < D:
                idx = len(self.added)
                self.added.append(Enew[i].astype(np.int32))
                self.queue.extend((idx, j) for j in range(nvar))
        return s

    def _product_rows(self, items):
        P = np.zeros((len(items), self.n), dtype=np.float64)
        for r, (idx, j) in enumerate(items):
            v = self.added[idx]
            sup = np.flatnonzero(v)
            dest = self.sp.cols[j, sup]
            P[r, dest] = v[sup]
        return P

    def add_inputs(self, rows):
        if len(rows):
            R = self._reduce(np.asarray(rows, dtype=np.float64) % self.p)
            self._absorb(R)

    def saturate(self, log=None, target_rank=None):
        """Run the queue.  If target_rank = dim(I cap M_<=D) is supplied, stop as soon as the rank
        reaches it: the span S satisfies S subset V_D subset I cap M_<=D, so equal dimensions give
        S = V_D exactly (no verification pass is needed for that closure)."""
        rounds = 0
        self.stopped_at_ideal_dimension = False
        while self.queue and not self.contains_one:
            if target_rank is not None and len(self.piv) >= target_rank:
                assert len(self.piv) == target_rank
                self.stopped_at_ideal_dimension = True
                break
            items = self.queue[:self.batch]
            del self.queue[:self.batch]
            R = self._reduce(self._product_rows(items))
            self._absorb(R)
            rounds += 1
            if log and rounds % 20 == 0:
                log(f'  batches={rounds} rank={len(self.piv)} queue={len(self.queue)}')
        return rounds

    def verify_closed(self):
        """Every variable shift of every basis row of degree < D reduces to zero."""
        if self.contains_one:
            return True
        D = self.sp.D
        low = [i for i in range(len(self.piv)) if self.sp.deg[self.piv[i]] < D]
        nvar = len(self.sp.cols)
        items = [(i, j) for i in low for j in range(nvar)]
        for start in range(0, len(items), self.batch):
            chunk = items[start:start + self.batch]
            P = np.zeros((len(chunk), self.n), dtype=np.float64)
            for r, (i, j) in enumerate(chunk):
                v = self.B[i]
                sup = np.flatnonzero(v)
                P[r, self.sp.cols[j, sup]] = v[sup]
            if np.any(self._reduce(P)):
                return False
        return True

    # -- queries -------------------------------------------------------------------------------
    @property
    def rank(self):
        # 1 in V_D generates every monomial of degree <= D (justified closure step)
        return self.n if self.contains_one else int(len(self.piv))

    def low_dimension(self, d):
        if self.contains_one:
            return int(np.sum(self.sp.deg <= d))
        return int(np.sum(self.sp.deg[self.piv] <= d)) if len(self.piv) else 0

    def contains(self, vec):
        return not np.any(self._reduce(np.asarray(vec, dtype=np.float64)[None, :] % self.p))


def batched_closure(generators, m, k, p, D, aux, *, batch=768, verify=True, log=None, target_rank=None):
    sp = MixedSpace(m, k, D, aux)
    cl = BatchedClosure(sp, p, batch)
    rows = []
    for g in generators:
        if not g:
            continue
        deg = max((key[0].bit_count() + key[1]) if isinstance(key, tuple) else key.bit_count() for key in g)
        if deg <= D:
            rows.append(sp.row(g, p))
    cl.add_inputs(rows)
    cl.saturate(log, target_rank)
    if target_rank is not None and not cl.contains_one and cl.rank == target_rank:
        cl.stopped_at_ideal_dimension = True
    if cl.stopped_at_ideal_dimension:
        verify = False  # exact by the dimension sandwich
    closed = cl.verify_closed() if verify else None
    if verify:
        assert closed, 'closure verification failed'
        for r in rows:
            assert cl.contains(r)
    cl.verified_closed = closed
    return cl
