"""Shared helpers for lane A (prime two-summand chain system). [pcc:unverified]"""
import os, sys, random
import numpy as np
sys.path.insert(0, '/home/user/ml-cryptanalysis/experiments/ic-curve-validation/code')
import lfd_engine as LE  # reference engine (independent, exact)


def is_prime(n):
    if n < 2: return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0: return n == q
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True


def prime_above(lo, k, hi=2**31 - 1):
    """smallest prime p >= lo with 2^k | p-1 (p < 2^31)"""
    N = 1 << k
    p = ((lo - 1) // N + 1) * N + 1
    while not is_prime(p):
        p += N
    assert p < hi
    return p


def prime_below(hi, k):
    """largest prime p <= hi with 2^k | p-1"""
    N = 1 << k
    p = ((hi - 1) // N) * N + 1
    while not is_prime(p):
        p -= N
    return p


def root_of_unity(p, N):
    """primitive N-th root of unity (N power of 2) in F_p"""
    assert (p - 1) % N == 0
    e = (p - 1) // N
    for g in range(2, p):
        w = pow(g, e, p)
        if pow(w, N // 2, p) == p - 1:
            return w


def S3(x1, x2, t, a, b, p):
    """Semaev S_3 (direct formula); works on ints or numpy int64 arrays (p < 2^31)."""
    d = (x1 - x2) % p
    s = (x1 + x2) % p
    pr = (x1 * x2) % p
    t2 = t * t % p
    A = d * d % p * t2 % p
    B = 2 * t % p * ((s * ((pr + a) % p) + 2 * b) % p) % p
    C = (pr - a) % p
    C = C * C % p
    return (A - B + C - 4 * b % p * s) % p


def gen_coeffs(a, b, t, p):
    """normalised generator g in F_p[(Z/N)^2]: dict (e1,e2)->coeff (README section 1)."""
    return {(2, 2): 1,
            (2, 1): (-2 * t) % p, (1, 2): (-2 * t) % p,
            (2, 0): t * t % p, (0, 2): t * t % p,
            (1, 1): (-2 * (t * t + a)) % p,
            (1, 0): (-(2 * a * t + 4 * b)) % p, (0, 1): (-(2 * a * t + 4 * b)) % p,
            (0, 0): (a * a - 4 * b * t) % p}


def g_eval(a, b, t, p, al, be):
    """evaluate normalised g at (al, be) via its coefficients (numpy arrays ok)."""
    tot = 0
    for (e1, e2), c in gen_coeffs(a, b, t, p).items():
        tot = (tot + c * pow_arr(al, e1, p) % p * pow_arr(be, e2, p)) % p
    return tot


def pow_arr(x, e, p):
    r = np.ones_like(x) if isinstance(x, np.ndarray) else 1
    for _ in range(e):
        r = r * x % p
    return r


def random_curve(p, rng):
    while True:
        a, b = rng.randrange(1, p), rng.randrange(1, p)
        if (4 * a ** 3 + 27 * b * b) % p:
            return a, b


def zeros_count(a, b, t, p, k, chunk=1 << 22):
    """number of (alpha,beta) in mu_N x mu_N with S_3(alpha,beta,t)=0 (0 => unit instance)."""
    N = 1 << k
    w = root_of_unity(p, N)
    H = np.empty(N, dtype=np.int64); H[0] = 1
    for i in range(1, N): H[i] = H[i - 1] * w % p
    rows = max(1, chunk // N)
    z = 0
    for s in range(0, N, rows):
        al = H[s:s + rows, None]
        v = S3(al, H[None, :], t, a, b, p)
        z += int(np.count_nonzero(v == 0))
    return z


def unit_instances(p, k, n, rng, max_try=2000):
    out, tries = [], 0
    while len(out) < n and tries < max_try:
        tries += 1
        a, b = random_curve(p, rng)
        t = rng.randrange(p)
        if zeros_count(a, b, t, p, k) == 0:
            out.append((a, b, t))
    return out, tries


def popcount_arr(n):
    x = np.arange(n, dtype=np.int64)
    c = np.zeros(n, dtype=np.int64)
    while x.any():
        c += x & 1; x >>= 1
    return c
