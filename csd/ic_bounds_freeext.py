"""Free high extension W_{k0+1}(S) of a two-summand chain closure and its certificates.

Setting and notation as in docs/proofs/N5-reduction.md (Sections 4-6) and
formal/ICFallDegrees/DepthReduction.lean.  At depth k an exponent pair a = (a1, a2), a_i in
[0, 2^k), is written a_i = l_i + 2^{k0} q_i with the low part l_i < 2^{k0} and the high pattern
q_i.  For a subspace S of the depth-k0 algebra inside M_{<=D}

    W_k(S) = (+)_q tau_q (S cap M_{<= D - wt q}),       tau_q : m_l -> m_{l + 2^{k0} q},

the sum over the high patterns q (here only k = k0 + 1, so q in {0,1}^2).  The kernel theorem
DepthReduction.W_closed / ChainModelMap.not_one_of_base needs (a) g in S, (b) W_{k0+1}(S)
closed at depth k0+1 under every u_{i,j} (j < k0+1) applied to its elements of weight < D,
(c) 1 notin S; it then gives 1 notin V_D(k) for every k >= k0+1, and V_D(k) subset W_k(S).

This module builds W_{k0+1}(S) for S a closure computed by the independent engine
csd/ic_bounds_d5_verify.py (its Layout / Part / Closure classes; reduced row echelon, columns by
weight descending, so the rows of S whose pivot has weight <= d span S cap M_{<=d}), and
certifies (b) by the deterministic full sweep of that engine run on W (every admissible shift of
every basis row of weight < D reduces to zero against W).  It imports nothing else.

The swap split is used when the generator is swap-symmetric (it always is for S_3): S and W are
then sigma-invariant, W = W^+ (+) W^-, and closedness under the block-1 shifts u_{0,j} implies
closedness under u_{1,j} = sigma u_{0,j} sigma.  split=False runs the same checks in full
coordinates with all 2k shifts (used as a cross-check on the small cases).
"""
from __future__ import annotations

import numpy as np

from csd import ic_bounds_d5_verify as V

F64 = np.float64


def popcount(x) -> np.ndarray:
    x = np.asarray(x, dtype=np.int64)
    c = np.zeros_like(x)
    y = x.copy()
    while np.any(y):
        c += y & 1
        y >>= 1
    return c


class Space:
    """A subspace held as engine Parts on a Layout, with the Closure helpers (sweeps, profile,
    contains_one) available through a Closure object whose __init__ is bypassed."""

    def __init__(self, k: int, p: int, D: int, split: bool):
        C = object.__new__(V.Closure)
        C.k, C.p, C.D = k, p, D
        C.L = V.Layout(k, D, split)
        C.parts = [V.Part(P['w'], p) for P in C.L.parts]
        C.processed = [set() for _ in C.parts]
        C.log = lambda **kw: None
        C.rounds = 0
        self.C = C

    def add_full(self, H: np.ndarray):
        for part, Y in zip(self.C.parts, self.C.L.project(H, self.C.p)):
            part.add(Y)

    def residual_rows(self, H: np.ndarray) -> int:
        """number of full-coordinate rows of H that do NOT lie in the space."""
        bad = 0
        for part, Y in zip(self.C.parts, self.C.L.project(H, self.C.p)):
            bad += int(np.count_nonzero(np.any(part.residual(Y) != 0, axis=1)))
        return bad


def as_space(C: V.Closure) -> Space:
    S = object.__new__(Space)
    S.C = C
    return S


def generator_matrix(g: dict, L: V.Layout, D: int, k: int, p: int, which: str = 'all') -> np.ndarray:
    """rows NF(m g) for squarefree m of weight <= D - deg g (which='all') or g alone ('g')."""
    N = 1 << k
    gdeg = max(bin(e1).count('1') + bin(e2).count('1') for e1, e2 in g)
    if which == 'all':
        ms = [(int(b1), int(b2)) for b1, b2, w in zip(L.a1, L.a2, L.w) if w <= D - gdeg]
    else:
        ms = [(0, 0)]
    terms = list(g.items())
    H = np.zeros((len(ms), L.c), dtype=F64)
    for r, (b1, b2) in enumerate(ms):
        cols = L.index(np.array([(b1 + e1) % N for (e1, e2), _ in terms]),
                       np.array([(b2 + e2) % N for (e1, e2), _ in terms]))
        for col, (_, v) in zip(cols, terms):
            H[r, col] = (H[r, col] + v) % p
    return H


def free_extension(S: V.Closure, chunk: int = 1024) -> Space:
    """W_{k0+1}(S) for S at depth k0 (same p, D, split)."""
    k0, p, D = S.k, S.p, S.D
    L0 = S.L
    W = Space(k0 + 1, p, D, L0.split)
    L1 = W.C.L
    for q1 in (0, 1):
        for q2 in (0, 1):
            wq = q1 + q2
            if wq > D:
                continue
            ok = L0.w <= D - wq
            colmap = -np.ones(L0.c, dtype=np.int64)
            colmap[ok] = L1.index(L0.a1[ok] + (q1 << k0), L0.a2[ok] + (q2 << k0))
            for pi, P in enumerate(S.parts):
                idx = np.flatnonzero(P.w[P.piv] <= D - wq)
                for a in range(0, len(idx), chunk):
                    F = L0.lift(pi, P.rows(idx[a:a + chunk]))
                    assert not np.any(F[:, ~ok]), 'row of S cap M_<=d with support of weight > d'
                    H = np.zeros((F.shape[0], L1.c), dtype=F64)
                    H[:, colmap[ok]] = F[:, ok]
                    W.add_full(np.mod(H, p))
    return W


def profile_recursion(prof: list) -> list:
    """p_j(k+1) = p_j(k) + 2 p_{j-1}(k) + p_{j-2}(k): the profile of W_{k+1}(S) from that of S."""
    def at(j):
        return prof[j] if j >= 0 else 0
    return [at(j) + 2 * at(j - 1) + at(j - 2) for j in range(len(prof))]


def same_space(A: Space, B: Space) -> bool:
    """equality of two spaces on the same Layout: canonical RREFs identical (pivots, free
    columns and coefficients), checked part by part."""
    for P, Q in zip(A.C.parts, B.C.parts):
        if P.rank != Q.rank:
            return False
        oa, ob = np.argsort(P.piv), np.argsort(Q.piv)
        if not np.array_equal(P.piv[oa], Q.piv[ob]) or not np.array_equal(P.free, Q.free):
            return False
        if not np.array_equal(P.B[oa], Q.B[ob]):
            return False
    return True


def shifted_generator_witnesses(S: V.Closure, g: dict, D: int, limit: int = 4000, want: int = 1,
                                order: str = 'weight') -> dict:
    """Search low exponent pairs e (e_i + max exponent of g < 2^{k0}: no carry out of the low
    window, so tau_e(g) is the same element at every depth >= k0) with every monomial of
    tau_e(g) of weight <= D, and test tau_e(g) notin S.  Candidates in order of decreasing
    wt(e) (the ones not reachable by generator rows), ties by (e1, e2)."""
    k0, p = S.k, S.p
    N = 1 << k0
    terms = list(g.items())
    emax = max(max(e1, e2) for (e1, e2), _ in terms)
    E = np.arange(N - emax)
    pc = {s: popcount(E + s) for s in range(emax + 1)}
    E1, E2 = np.meshgrid(E, E, indexing='ij')
    E1, E2 = E1.ravel(), E2.ravel()
    mw = np.zeros(E1.shape, dtype=np.int64)
    for (s1, s2), _ in terms:
        mw = np.maximum(mw, pc[s1][E1] + pc[s2][E2])
    cand = np.flatnonzero(mw <= D)
    we = popcount(E1[cand]) + popcount(E2[cand])
    o = np.lexsort((E2[cand], E1[cand], -we))
    cand = cand[o]
    sp = as_space(S)
    found, tested = [], 0
    gdeg = max(bin(a).count('1') + bin(b).count('1') for a, b in g)
    for a in range(0, min(len(cand), limit), 256):
        blk = cand[a:min(a + 256, limit)]
        H = np.zeros((len(blk), S.L.c), dtype=F64)
        for r, c in enumerate(blk):
            e1, e2 = int(E1[c]), int(E2[c])
            cols = S.L.index(np.array([e1 + s1 for (s1, s2), _ in terms]),
                             np.array([e2 + s2 for (s1, s2), _ in terms]))
            for col, (_, v) in zip(cols, terms):
                H[r, col] = (H[r, col] + v) % p
        res = [part.residual(Y) for part, Y in zip(S.parts, S.L.project(H, p))]
        out = np.zeros(len(blk), dtype=bool)
        for R in res:
            out |= np.any(R != 0, axis=1)
        tested += len(blk)
        for r in np.flatnonzero(out):
            c = blk[r]
            e1, e2 = int(E1[c]), int(E2[c])
            found.append(dict(e=[e1, e2], wt_e=int(popcount(e1) + popcount(e2)),
                              max_weight=int(mw[c]),
                              support=sorted([[e1 + s1, e2 + s2] for (s1, s2), _ in terms]),
                              coefficients={f'{s1},{s2}': int(v) for (s1, s2), v in sorted(terms)}))
            if len(found) >= want:
                break
        if len(found) >= want:
            break
    return dict(candidates_low_weight=int(len(cand)), tested=tested, found=found,
                generator_degree=gdeg, rule='e_i + %d < 2^%d, all weights of tau_e(g) <= %d' % (emax, k0, D))


def free_extension_direct(S: V.Closure) -> Space:
    """W_{k0+1}(S) assembled without elimination over the whole space (split layout only).

    The four components tau_q(S cap M_{<=D-wt q}) have disjoint column supports.  In the swap
    split of the depth-(k0+1) layout:
      * q = (0,0) and q = (1,1): the orbit representatives (a1 <= a2) of the component are the
        images of the depth-k0 representatives, and the column order (weight descending, then
        key) is preserved, so the rows of S^{+/-} with pivot weight <= D - wt q map to reduced
        echelon rows unchanged;
      * q in {(1,0),(0,1)} (exchanged by sigma): every orbit of the pair of components has its
        representative in tau_{(0,1)}, and in part coordinates the component of W^{+/-} is
        {sigma v : v in S cap M_{<=D-1}} = S cap M_{<=D-1} (S is sigma-invariant), placed by
        (m1, m2) -> (m1, m2 + 2^{k0}) in the depth-k0 full order; its RREF is computed here by
        eliminating the full-coordinate rows of S cap M_{<=D-1} (small).
    The union of the component RREFs is the RREF of W (disjoint supports, order preserved within
    each component).  Cross-checked against free_extension() on small cases (same_space)."""
    k0, p, D = S.k, S.p, S.D
    L0 = S.L
    assert L0.split
    W = Space(k0 + 1, p, D, True)
    L1 = W.C.L
    posW = []
    for P in L1.parts:
        pos = -np.ones(L1.c, dtype=np.int64)
        pos[P['reps']] = np.arange(len(P['reps']))
        posW.append(pos)
    # RREF of S cap M_{<=D-1} in full depth-k0 coordinates restricted to weight <= D-1 columns
    lowcols = np.flatnonzero(L0.w <= D - 1)
    M = V.Part(L0.w[lowcols], p)
    for pi, P in enumerate(S.parts):
        idx = np.flatnonzero(P.w[P.piv] <= D - 1)
        for a in range(0, len(idx), 1024):
            F = L0.lift(pi, P.rows(idx[a:a + 1024]))
            assert not np.any(np.delete(F, lowcols, axis=1))
            M.add(np.mod(F[:, lowcols], p))
    mixed_map_full = L1.index(L0.a1[lowcols], L0.a2[lowcols] + (1 << k0))
    for wi, (WP, PW) in enumerate(zip(W.C.parts, L1.parts)):
        pos = posW[wi]
        SP, SPd = S.parts[wi], L0.parts[wi]
        blocks = []   # (pivot target cols, Part, row idx, colmap from that Part's columns)
        for q, dmax in (((0, 0), D), ((1, 1), D - 2)):
            if dmax < 0:
                continue
            full = SPd['reps']
            ok = L0.w[full] <= dmax
            cmap = -np.ones(len(full), dtype=np.int64)
            cmap[ok] = pos[L1.index(L0.a1[full[ok]] + (q[0] << k0), L0.a2[full[ok]] + (q[1] << k0))]
            assert np.all(cmap[ok] >= 0)
            idx = np.flatnonzero(SP.w[SP.piv] <= dmax)
            blocks.append((SP, idx, cmap, ok))
        cm = pos[mixed_map_full]
        assert np.all(cm >= 0)
        blocks.append((M, np.arange(M.rank), cm, np.ones(len(cm), dtype=bool)))
        piv = np.concatenate([cmap[P.piv[idx]] for P, idx, cmap, ok in blocks])
        assert len(np.unique(piv)) == len(piv) and np.all(piv >= 0)
        free = np.setdiff1d(np.arange(PW['c']), piv)
        B = np.zeros((len(piv), len(free)), dtype=np.float32)
        where = -np.ones(PW['c'], dtype=np.int64)
        where[free] = np.arange(len(free))
        r0 = 0
        for P, idx, cmap, ok in blocks:
            fc = cmap[P.free]                     # target columns of the block's free columns
            good = fc >= 0
            for a in range(0, len(idx), 2048):
                rows = P.B[idx[a:a + 2048]]
                assert not np.any(rows[:, ~good]), 'row with support outside its component'
                B[r0 + a:r0 + a + len(rows)][:, where[fc[good]]] = rows[:, good]
            r0 += len(idx)
        assert np.all(where[np.concatenate([cmap[P.free][cmap[P.free] >= 0] for P, _, cmap, _ in blocks])] >= 0)
        WP.piv, WP.free, WP.B = piv.astype(np.int64), free, B
        WP.peak_entries = B.size
    return W
