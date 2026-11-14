"""N8-growth: the splitting property FEI_k(w) alone, for larger cases (blocked elimination with float64 BLAS).

FEI_k(w): for every polynomial f of weight <= w with exponents < 2^k vanishing on C, every coordinate i and
level j, the part of f on the monomials with a_i >= 2^k - 2^j vanishes on C.  By Proposition 2.1 of
docs/proofs/N8-growth.md, FEI_k(D-1) <=> the no-wrap ideal space S_D(k) is closed under the degree-D rule,
which implies 1 not in V_D(k).

Method: S_w(k) = kernel of the evaluation matrix of the monomials of weight <= w on points of C.  The kernel
is computed from cols + 64 points and then every basis vector is checked on 4(2^k - 1) + 65 further points
(more than 4(2^k - 1), so a vector passing the check lies in the ideal; a polynomial of the ideal is in the
computed kernel).  A top-run part lies in the ideal iff it reduces to zero against the echelon basis of S.
usage: ic_bounds_growth_fei.py p a b t w k [out.json]
"""
import json, sys, time
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from ic_bounds_growth_lib import Layout, curve_points, rref

def mm(A, B, p, step=192):
    """(A @ B) mod p exactly in float64: entries < p <= 2^16, partial sums of `step` products < 2^53."""
    out = np.zeros((A.shape[0], B.shape[1]))
    for s in range(0, A.shape[1], step):
        out += A[:, s:s + step] @ B[s:s + step]
        np.mod(out, p, out=out)
    return out

def rref_blocked(M, p, bs=192):
    """Reduced row echelon form mod p of a float64 integer matrix, in place. Returns pivot list (col -> row order)."""
    rows, cols = M.shape
    r = 0; piv = []
    for c0 in range(0, cols, bs):
        if r == rows: break
        c1 = min(cols, c0 + bs)
        panel = M[:, c0:c1].astype(np.int64)
        # choose pivots in the panel by elimination on rows r.. (small int64 work)
        work = panel[r:].copy(); loc = []; prow = []
        used = np.zeros(work.shape[0], dtype=bool)
        for c in range(c1 - c0):
            cand = np.nonzero((work[:, c] != 0) & ~used)[0]
            if cand.size == 0: continue
            s = cand[0]; used[s] = True
            work[s] = work[s] * pow(int(work[s, c]), p - 2, p) % p
            col = work[:, c].copy(); col[s] = 0
            nz = np.nonzero(col)[0]
            if nz.size: work[nz] = (work[nz] - np.outer(col[nz], work[s])) % p
            loc.append(c); prow.append(r + s)
        if not loc: continue
        n = len(loc)
        App = panel[prow][:, loc]
        aug = np.concatenate([App, np.eye(n, dtype=np.int64)], axis=1)
        R, pv = rref(aug, p); assert pv == list(range(n))
        inv = R[:, n:].astype(np.float64)
        top = mm(inv, M[prow, c0:], p)                    # new pivot rows on columns >= c0
        coef = panel[:, loc].astype(np.float64)
        coef[prow] = 0
        M[:, c0:] -= mm(coef, top, p)
        np.mod(M[:, c0:], p, out=M[:, c0:])
        M[prow, c0:] = top
        # move pivot rows to positions r.. in order
        for t_ in range(n):
            pr = prow[t_]
            if pr != r + t_:
                M[[r + t_, pr]] = M[[pr, r + t_]]
                prow = [(pr if q == r + t_ else q) for q in prow]
                prow[t_] = r + t_
        piv += [c0 + c for c in loc]; r += n
    return piv

def main():
    p, a, b, t, w, k = map(int, sys.argv[1:7])
    t00 = time.time()
    L = Layout(k, w); N = L.N
    pts = curve_points(p, a, b, t)
    need = 4 * (N - 1) + 1
    rec = dict(p=p, a=a, b=b, t=t, w=w, k=k, columns=L.n, points=len(pts), points_needed=need, exact=len(pts) >= need)
    rng = np.random.default_rng(11)
    sel = rng.permutation(len(pts))
    x1 = np.array([q[0] for q in pts], dtype=np.int64); x2 = np.array([q[1] for q in pts], dtype=np.int64)
    def evalm(idx):
        P1 = np.ones((len(idx), N), dtype=np.int64); P2 = np.ones((len(idx), N), dtype=np.int64)
        for e in range(1, N):
            P1[:, e] = P1[:, e - 1] * x1[idx] % p; P2[:, e] = P2[:, e - 1] * x2[idx] % p
        return (P1[:, L.a1] * P2[:, L.a2] % p).astype(np.float64)
    nrow = min(len(pts), L.n + 64)
    E = evalm(sel[:nrow])
    piv = rref_blocked(E, p)
    rank = len(piv); rec['rank_eval_sample'] = rank
    free = np.array([c for c in range(L.n) if c not in set(piv)], dtype=np.int64)
    # kernel basis: row for free column f: e_f - sum_r E[r, f] e_{piv[r]}
    Kb = np.zeros((len(free), L.n))
    Kb[np.arange(len(free)), free] = 1
    Kb[:, piv] = (-E[:rank][:, free].T) % p
    del E
    # check on all points
    bad = 0
    ncheck = min(len(pts), need + 64)          # more than 4(2^k - 1) points: exact
    rec['points_checked'] = ncheck
    for s in range(0, ncheck, 1024):
        Ec = evalm(sel[::-1][s:min(ncheck, s + 1024)])
        bad += int(np.count_nonzero(mm(Kb, Ec.T.copy(), p)))
    rec['kernel_vectors_failing_on_all_points'] = bad
    rec['dim_S'] = int(Kb.shape[0])
    # echelon form of the kernel basis (pivot = free column, already reduced: Kb[:, free] = I)
    def in_S(Y):
        """rows Y reduce to zero against Kb (which is reduced w.r.t. the free columns)."""
        res = Y - mm(Y[:, free].copy(), Kb, p)
        np.mod(res, p, out=res)
        return ~np.any(res != 0, axis=1)
    fail = {}; checked = 0
    for i in range(2):
        ai = L.a1 if i == 0 else L.a2
        for j in range(k):
            mask = ai >= N - (1 << j)
            if not mask.any(): continue
            rows = np.nonzero(np.any(Kb[:, mask] != 0, axis=1))[0]
            if rows.size == 0: continue
            Y = Kb[rows] * mask
            ok = in_S(Y); checked += rows.size
            if not ok.all(): fail[f'{i},{j}'] = int((~ok).sum())
    rec['FEI'] = (not fail) and bad == 0; rec['FEI_fail'] = fail; rec['top_run_parts_checked'] = checked
    rec['seconds'] = round(time.time() - t00, 1)
    print(rec, flush=True)
    if len(sys.argv) > 7:
        json.dump(rec, open(sys.argv[7], 'w'), indent=1)

main()
