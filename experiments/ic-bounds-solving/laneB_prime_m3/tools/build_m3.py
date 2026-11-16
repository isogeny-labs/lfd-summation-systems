#!/usr/bin/env python3
"""Lane B: build g = NF(S_4(x1,x2,x3,t)) in A = F_p[(Z/2^k)^3] (array (N,N,N)), plus NTT helpers.
Independent of the repository code except in the explicit cross-check (--selftest)."""
import sys, os, time, json, random
import numpy as np

# ---------- small number theory ----------
def is_prime(n):
    if n < 2: return False
    for q in (2, 3, 5, 7, 11, 13):
        if n % q == 0: return n == q
    i = 17
    while i * i <= n:
        if n % i == 0: return False
        i += 2
    return True

def primes_1mod(N, lo, hi=None, count=1, descending=False):
    out = []
    if descending:
        m = (hi - 1) // N * N + 1
        while m > lo and len(out) < count:
            if is_prime(m): out.append(m)
            m -= N
    else:
        m = (lo // N + 1) * N + 1
        while len(out) < count:
            if is_prime(m): out.append(m)
            m += N
    return out

def root_N(p, N):
    """Same convention as the repository: first h with g=h^((p-1)/N) of exact order N."""
    for h in range(2, p):
        g = pow(h, (p - 1) // N, p)
        if pow(g, N // 2, p) != 1:
            return g
    raise ValueError

# ---------- polynomials as {(e1,e2,e3): c} ----------
def padd(p, *terms):
    out = {}
    for s, f in terms:
        for e, c in f.items():
            out[e] = (out.get(e, 0) + s * c) % p
    return {e: c for e, c in out.items() if c}

def pmul(p, f, g):
    out = {}
    for e1, c1 in f.items():
        for e2, c2 in g.items():
            e = (e1[0] + e2[0], e1[1] + e2[1], e1[2] + e2[2])
            out[e] = (out.get(e, 0) + c1 * c2) % p
    return {e: c for e, c in out.items() if c}

def s3_abc(p, a, b, u, v):
    one = {(0, 0, 0): 1}
    uv = pmul(p, u, v)
    d = padd(p, (1, u), (-1, v))
    A = pmul(p, d, d)
    s = padd(p, (1, u), (1, v))
    B = padd(p, (-2, pmul(p, s, padd(p, (1, uv), (a, one)))), (-4 * b, one))
    w = padd(p, (1, uv), (-a, one))
    C = padd(p, (1, pmul(p, w, w)), (-4 * b, s))
    return A, B, C

def s4_poly(p, a, b, t):
    """S_4(x1,x2,x3,t) = Res_X(S_3(x1,x2,X), S_3(x3,t,X)) = (A C' - A' C)^2 - (A B' - A' B)(B C' - B' C)."""
    x1, x2, x3 = {(1, 0, 0): 1}, {(0, 1, 0): 1}, {(0, 0, 1): 1}
    tc = {(0, 0, 0): t % p}
    A, B, C = s3_abc(p, a, b, x1, x2)
    A2, B2, C2 = s3_abc(p, a, b, x3, tc)
    m = lambda f, g: pmul(p, f, g)
    u = padd(p, (1, m(A, C2)), (-1, m(A2, C)))
    v = padd(p, (1, m(A, B2)), (-1, m(A2, B)))
    w = padd(p, (1, m(B, C2)), (-1, m(B2, C)))
    return padd(p, (1, m(u, u)), (-1, m(v, w)))

def popcount_table(N):
    return np.array([bin(i).count('1') for i in range(N)], dtype=np.int16)

def build_g(p, a, b, t, k):
    """Return (g array (N,N,N) int64 in [0,p), delta). Exponents reduced mod N."""
    N = 1 << k
    S = s4_poly(p, a, b, t)
    g = np.zeros((N, N, N), dtype=np.int64)
    for (e1, e2, e3), c in S.items():
        g[e1 % N, e2 % N, e3 % N] = (g[e1 % N, e2 % N, e3 % N] + c) % p
    pc = popcount_table(N)
    W = pc[:, None, None] + pc[None, :, None] + pc[None, None, :]
    nz = g != 0
    delta = int(W[nz].max()) if nz.any() else -1
    return g, delta

# ---------- NTT over axes (p < 2^31, int64) ----------
def _ntt_last(x, w, p):
    N = x.shape[-1]
    if N == 1:
        return x
    w2 = w * w % p
    even = _ntt_last(x[..., 0::2], w2, p)
    odd = _ntt_last(x[..., 1::2], w2, p)
    tw = np.empty(N // 2, dtype=np.int64)
    c = 1
    for j in range(N // 2):
        tw[j] = c
        c = c * w % p
    t = odd * tw % p
    return np.concatenate([(even + t) % p, (even - t) % p], axis=-1)

def dft3(a, w, p):
    """v[E] = sum_e a[e] w^{e.E}  (w of order N)."""
    for ax in range(3):
        a = np.ascontiguousarray(np.moveaxis(_ntt_last(np.ascontiguousarray(np.moveaxis(a, ax, -1)), w, p), -1, ax))
    return a

def idft3(v, w, p):
    N = v.shape[0]
    wi = pow(w, -1, p)
    a = dft3(v, wi, p)
    Ninv3 = pow(N ** 3, -1, p)
    return a * Ninv3 % p

def modinv_arr(x, p):
    r = np.ones_like(x)
    base = x.copy()
    e = p - 2
    while e:
        if e & 1: r = r * base % p
        base = base * base % p
        e >>= 1
    return r

def random_instance(rng, p):
    """Random (a,b,t) nonsingular curve, t the x-coordinate of an F_p-point (rhs nonzero square)."""
    while True:
        a, b = rng.randrange(p), rng.randrange(p)
        if (4 * a ** 3 + 27 * b * b) % p == 0: continue
        t = rng.randrange(p)
        r = (t ** 3 + a * t + b) % p
        if r and pow(r, (p - 1) // 2, p) == 1:
            return a, b, t

def locus_instance(rng, p):
    """5t^3 + 3at + 6b = 0 locus (t root), t an x-coordinate of a point."""
    while True:
        a, t = rng.randrange(p), rng.randrange(1, p)
        # b = -(5t^3+3at)/6
        b = (-(5 * t ** 3 + 3 * a * t)) * pow(6, -1, p) % p
        if (4 * a ** 3 + 27 * b * b) % p == 0: continue
        r = (t ** 3 + a * t + b) % p
        if r and pow(r, (p - 1) // 2, p) == 1:
            return a, b, t

# ---------- self test ----------
def selftest():
    sys.path.insert(0, '/home/user/ml-cryptanalysis')
    from csd import ic_dreg as I, ic_bounds_m3 as M
    rng = random.Random(1)
    rec = {}
    for (p, k) in [(97, 2), (97, 3), (257, 4), (193, 3)]:
        a, b, t = random_instance(rng, p)
        g, delta = build_g(p, a, b, t, k)
        E = I.PrimeCurve(p, a, b)
        ref = M.direct_generator(E, t, k)  # {mask: c}, mask = e1 | e2<<k | e3<<2k
        mine = {}
        N = 1 << k
        for idx in np.argwhere(g):
            e1, e2, e3 = map(int, idx)
            mine[e1 | (e2 << k) | (e3 << 2 * k)] = int(g[e1, e2, e3])
        same = mine == {int(m): int(c) % p for m, c in ref.items()}
        refdelta = max(bin(m).count('1') for m in ref)
        # vanishing set vs brute force with the repository's numeric s4_value
        w = root_N(p, N)
        v = dft3(g, w, p)
        H = [pow(w, e, p) for e in range(N)]
        bad = 0
        for E1 in range(N):
            for E2 in range(N):
                for E3 in range(N):
                    zr = M.s4_value(E, H[E1], H[E2], H[E3], t) == 0
                    zm = v[E1, E2, E3] == 0
                    bad += (zr != zm)
        # inverse check by direct cyclic convolution (small k)
        unit = not (v == 0).any()
        invok = None
        if unit:
            gi = idft3(modinv_arr(v, p), w, p)
            # convolution in A: sum over e: g[e] gi[f-e]
            prod = np.zeros_like(g)
            for idx in np.argwhere(g):
                e1, e2, e3 = map(int, idx)
                prod = (prod + g[e1, e2, e3] * np.roll(gi, (e1, e2, e3), axis=(0, 1, 2))) % p
            one = np.zeros_like(g); one[0, 0, 0] = 1
            invok = bool((prod == one).all())
        rec[f'p{p}_k{k}'] = dict(a=a, b=b, t=t, coeffs_equal_repo=same, delta=delta, repo_delta=refdelta,
                                 vanishing_mismatches=int(bad), n_zero=int((v == 0).sum()), unit=unit, inverse_ok=invok)
        print(p, k, rec[f'p{p}_k{k}'], flush=True)
    # a unit instance with Z empty, searched at p=193,k=2 (rate small) -> larger prime
    p, k = 1009 if False else primes_1mod(4, 4000)[0], 2
    found = None
    for _ in range(2000):
        a, b, t = random_instance(rng, p)
        g, delta = build_g(p, a, b, t, k)
        v = dft3(g, root_N(p, 4), p)
        if not (v == 0).any():
            found = (a, b, t); break
    E = I.PrimeCurve(p, a, b)
    w = root_N(p, 4); H = [pow(w, e, p) for e in range(4)]
    zs = sum(M.s4_value(E, H[i], H[j], H[l], t) == 0 for i in range(4) for j in range(4) for l in range(4))
    rec['unit_instance'] = dict(p=p, k=2, abt=found, brute_zero_count=int(zs))
    print(rec['unit_instance'])
    out = os.path.join(os.path.dirname(__file__), '..', 'data', 'selftest_build_m3.json')
    json.dump(rec, open(out, 'w'), indent=1, default=int)

if __name__ == '__main__':
    if '--selftest' in sys.argv:
        selftest()
