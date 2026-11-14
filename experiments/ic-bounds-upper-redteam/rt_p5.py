"""Independent recomputation of the Proposition 5 example (k = 8, p = 257) in the group algebra
F_p[(Z/2^k)^2] with Hamming-weight degree, and of the arithmetic in E2 / Section 6.
Own code for monomials, evaluation, kernels and shifts; only rank_mod is taken from lfd_engine."""
import sys, math
import numpy as np
from rt_lib import s3_eval, mu, E

def wt(e): return bin(e).count('1')

def mons(k, D):
    N = 1 << k
    return [(e1, e2) for e1 in range(N) for e2 in range(N) if wt(e1) + wt(e2) <= D]

def evalmat(Z, ms, p):
    M = np.zeros((len(Z), len(ms)), dtype=np.int64)
    for i, (x1, x2) in enumerate(Z):
        for j, (e1, e2) in enumerate(ms):
            M[i, j] = pow(x1, e1, p) * pow(x2, e2, p) % p
    return M

def kernel(M, p):
    """Basis of {v : M v = 0} over F_p by column-style Gauss-Jordan on M (rows = points)."""
    A = M.copy() % p
    rows, cols = A.shape
    piv = []
    r = 0
    for c in range(cols):
        if r == rows: break
        nz = np.flatnonzero(A[r:, c])
        if nz.size == 0: continue
        i = r + int(nz[0])
        if i != r: A[[r, i]] = A[[i, r]]
        A[r] = A[r] * pow(int(A[r, c]), -1, p) % p
        for i2 in np.flatnonzero(A[:, c]):
            if i2 != r:
                A[i2] = (A[i2] - A[i2, c] * A[r]) % p
        piv.append(c); r += 1
    free = [c for c in range(cols) if c not in set(piv)]
    K = np.zeros((len(free), cols), dtype=np.int64)
    for n, fcol in enumerate(free):
        K[n, fcol] = 1
        for ri, pc in enumerate(piv):
            K[n, pc] = (-A[ri, fcol]) % p
    assert not np.any((M @ K.T) % p)
    return K, len(piv)

def shift_embed(K, ms_from, ms_to, k, var):
    """(variable var) * rows, re-expressed on ms_to; var = None is the identity embedding."""
    N = 1 << k
    idx = {m: i for i, m in enumerate(ms_to)}
    out = np.zeros((K.shape[0], len(ms_to)), dtype=np.int64)
    for j, (e1, e2) in enumerate(ms_from):
        if var is None: m = (e1, e2)
        elif var[0] == 0: m = ((e1 + (1 << var[1])) % N, e2)
        else: m = (e1, (e2 + (1 << var[1])) % N)
        out[:, idx[m]] += K[:, j]          # a shift permutes monomials, so no collisions
    return out

def instance(p, a, b, t, k, Dtop):
    H = mu(p, 1 << k)
    Z = [(x1, x2) for x1 in H for x2 in H if s3_eval(p, a, b, t, x1, x2) == 0]
    print(f'p={p} (a,b,t)=({a},{b},{t}) k={k} |Z|={len(Z)}', flush=True)
    prev = None
    for D in range(2, Dtop + 1):
        ms = mons(k, D)
        assert len(ms) == sum(math.comb(2 * k, j) for j in range(D + 1))
        K, rank = kernel(evalmat(Z, ms, p), p)
        line = f'  D={D} |M|={len(ms)} eval_rank={rank} dim I_D={K.shape[0]}'
        if prev is not None:
            pms, pK = prev
            gens = [shift_embed(pK, pms, ms, k, None)]
            for i in range(2):
                for j in range(k):
                    gens.append(shift_embed(pK, pms, ms, k, (i, j)))
            span = E.rank_mod(np.vstack(gens) % p, p)
            line += f'  dim(I_{D-1} + sum v I_{D-1})={span}  (2k+1)*dim I_{D-1}={(2*k+1)*pK.shape[0]}  generation={"holds" if span == K.shape[0] else "FAILS"}'
        print(line, flush=True)
        prev = (ms, K)

if __name__ == '__main__':
    instance(257, 1, 3, 190, 8, 4)
    instance(257, 5, 7, 11, 8, 4)
    # thresholds D*(k) and E2 arithmetic
    Mk = lambda k, D: sum(math.comb(2 * k, j) for j in range(D + 1))
    print('D*(k):', [(k, next(D for D in range(1, 2 * k + 1) if Mk(k, D) >= 2 ** k), (k + 14) // 4, k + 2) for k in (8, 16, 24, 32, 48, 64, 96, 128)])
    f4 = lambda k: 4 * k * k + 172 * k - 571
    f5 = lambda k: (8 * k ** 3 + 576 * k * k + 4648 * k - 56373)
    print('E2 D=4:', [(k, f4(k), Mk(k, 4), f4(k) > Mk(k, 4)) for k in range(2, 12)])
    print('E2 D=5:', [(k, f5(k) // 3 if f5(k) % 3 == 0 else f5(k) / 3, Mk(k, 5), f5(k) > 3 * Mk(k, 5)) for k in range(4, 16)])
    Q5 = [0, 0, 2, 98, 823, -18791]
    for k, prof in {10: (0, 0, 2, 138, 3163, 18569), 11: (0, 0, 2, 142, 3441, 25033), 12: (0, 0, 2, 146, 3727, 32057), 13: (0, 0, 2, 150, 4021, 39657)}.items():
        got = tuple(sum(Q5[i] * math.comb(2 * k, j - i) for i in range(j + 1)) for j in range(6))
        assert got == prof, (k, got, prof)
    print('Q5 reproduces the four D=5 profiles; f5 polynomial matches Q5:', all(sum(Q5[i] * math.comb(2 * k, 5 - i) for i in range(6)) * 3 == f5(k) for k in range(3, 40)))
    print('k=16 numbers: f4(16)=', f4(16), ' p4 of V5 =', 2 * math.comb(32, 2) + 98 * 32 + 823, ' threshold (M5-2^16)/33 =', (Mk(16, 5) - 2 ** 16) / 33, ' M4(32)=', Mk(16, 4))
    import bisect
    H = lambda x: -x * math.log2(x) - (1 - x) * math.log2(1 - x)
    lo, hi = 1e-9, 0.5
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if H(mid) < 0.5 else (lo, mid)
    print('2*Hinv(1/2)=', 2 * lo, ' 2^(8*H(1/8))=', 2 ** (8 * H(1 / 8)), ' M6(26)=', Mk(13, 6), ' M6(28)=', Mk(14, 6))
