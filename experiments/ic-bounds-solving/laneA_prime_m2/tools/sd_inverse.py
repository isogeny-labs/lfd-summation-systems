"""g^{-1} in A = F_p[(Z/N)^2], N = 2^k, by evaluation on mu_N x mu_N, pointwise inversion and inverse DFT
(two 1-D radix-2 modular transforms, exact int64 arithmetic, p < 2^31).  Reports maxwt(supp g^{-1}),
zero counts per weight class, and the (N-1,N-1) coefficient.  Prediction P1: sd = 2 + maxwt.
Usage: python3 sd_inverse.py k p a b t   (single instance, prints JSON)    [pcc:unverified]
"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from common import *


def bitrev(N):
    k = N.bit_length() - 1
    r = np.zeros(N, dtype=np.int64)
    for i in range(N):
        r[i] = int(format(i, '0%db' % k)[::-1], 2) if k else 0
    return r


def ntt_rows(a, w, p, rev):
    """In place DFT along axis 0 of a (shape (N, M), values in [0,p)): out[e] = sum_s a[s] w^{e s}."""
    N = a.shape[0]
    a[:] = a[rev]
    ln = 2
    while ln <= N:
        half = ln // 2
        wl = pow(w, N // ln, p)
        tw = np.empty(half, dtype=np.int64); tw[0] = 1
        for j in range(1, half): tw[j] = tw[j - 1] * wl % p
        v = a.reshape(N // ln, ln, -1)
        u0 = v[:, :half].copy()
        t = v[:, half:] * tw[None, :, None] % p
        v[:, :half] = (u0 + t) % p
        v[:, half:] = (u0 - t) % p
        ln *= 2


def batch_inv(x, p):
    """elementwise inverse mod p by Fermat exponentiation (all entries must be nonzero)."""
    e = p - 2
    result = np.ones_like(x); base = x.copy()
    while e:
        if e & 1: result = result * base % p
        e >>= 1
        if e: base = base * base % p
    return result


def inverse_coeffs(a, b, t, p, k, chunk_cols=1 << 20):
    """returns array C (N x N int64) with g^{-1} = sum C[e1,e2] x1^e1 x2^e2, or None if g has a zero on H x H."""
    N = 1 << k
    w = root_of_unity(p, N); wi = pow(w, -1, p)
    H = np.empty(N, dtype=np.int64); H[0] = 1
    for i in range(1, N): H[i] = H[i - 1] * w % p
    Psi = np.empty((N, N), dtype=np.int64)
    for s in range(N):
        Psi[s] = g_eval(a, b, t, p, H[s], H)  # g at (w^s, w^s2) using normalised coefficients
    if (Psi == 0).any():
        return None
    for s in range(0, N, max(1, (1 << 22) // N)):
        Psi[s:s + max(1, (1 << 22) // N)] = batch_inv(Psi[s:s + max(1, (1 << 22) // N)], p)
    rev = bitrev(N)
    ntt_rows(Psi, wi, p, rev)                      # axis 0
    Psi = np.ascontiguousarray(Psi.T)
    ntt_rows(Psi, wi, p, rev)
    Psi = Psi.T
    scale = pow(N * N % p, -1, p)
    return np.ascontiguousarray(Psi) * scale % p


def check_inverse_exact(C, a, b, t, p, k):
    """exact cyclic convolution g * C == 1 (9-term g), O(9 N^2)."""
    N = 1 << k
    out = np.zeros((N, N), dtype=np.int64)
    for (e1, e2), c in gen_coeffs(a, b, t, p).items():
        out = (out + c * np.roll(np.roll(C, e1, axis=0), e2, axis=1)) % p
    ok = out[0, 0] == 1 and np.count_nonzero(out) == 1
    return bool(ok)


def analyse(C, k):
    N = 1 << k
    pc = popcount_arr(N)
    W = pc[:, None] + pc[None, :]
    nz = C != 0
    maxwt = int(W[nz].max()) if nz.any() else -1
    zeros_per_w = [int(np.count_nonzero(~nz & (W == w))) for w in range(2 * k + 1)]
    return dict(maxwt=maxwt, zeros_per_weight=zeros_per_w, zero_top_coeff=bool(C[N - 1, N - 1] == 0),
                nzeros=int(np.count_nonzero(~nz)))


def sd_inverse(a, b, t, p, k, verify=True):
    C = inverse_coeffs(a, b, t, p, k)
    if C is None: return None
    res = analyse(C, k)
    if verify: res['inverse_verified'] = check_inverse_exact(C, a, b, t, p, k)
    res['sd_inverse'] = 2 + res['maxwt']
    return res


if __name__ == '__main__':
    k, p, a, b, t = map(int, sys.argv[1:6])
    t0 = time.time(); r = sd_inverse(a, b, t, p, k); r['sec'] = round(time.time() - t0, 2)
    print(json.dumps(r))
