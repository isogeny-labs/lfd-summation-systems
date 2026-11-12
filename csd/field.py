"""Prime-field helpers used by the toy-scale curve pipeline.

Everything here is exact integer arithmetic. Nothing is performance
critical beyond the point-count sweep, which lives in curve.py and uses
numpy.
"""
from __future__ import annotations

import math
import random
from functools import lru_cache

_SMALL_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)


def is_prime(n: int) -> bool:
    """Deterministic Miller-Rabin for n < 3.3e24 (bases up to 37)."""
    if n < 2:
        return False
    for p in _SMALL_PRIMES:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in _SMALL_PRIMES:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def factor(n: int) -> dict[int, int]:
    """Trial-division factorisation; fine for the n < 2**40 used here."""
    n = abs(n)
    out: dict[int, int] = {}
    if n < 2:
        return out
    for p in (2, 3):
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
    f = 5
    while f * f <= n:
        for q in (f, f + 2):
            while n % q == 0:
                out[q] = out.get(q, 0) + 1
                n //= q
        f += 6
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def largest_prime_factor(n: int) -> int:
    return max(factor(n)) if abs(n) > 1 else 1


def legendre(a: int, p: int) -> int:
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def sqrt_mod(a: int, p: int) -> int | None:
    """Tonelli-Shanks. Returns a root or None if `a` is a non-residue."""
    a %= p
    if a == 0:
        return 0
    if legendre(a, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while legendre(z, p) != -1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r


def primitive_root_of_unity(n: int, p: int) -> int | None:
    """A primitive n-th root of unity in F_p, or None when n does not divide p-1."""
    if (p - 1) % n:
        return None
    e = (p - 1) // n
    for g in range(2, p):
        z = pow(g, e, p)
        if z != 1 and all(pow(z, n // q, p) != 1 for q in factor(n)):
            return z
    return None


def multiplicative_order(a: int, n: int, limit: int) -> int | None:
    """Smallest k <= limit with a**k == 1 mod n, else None."""
    x = 1
    for k in range(1, limit + 1):
        x = x * a % n
        if x == 1:
            return k
    return None


@lru_cache(maxsize=None)
def primes_in_range(lo: int, hi: int) -> tuple[int, ...]:
    sieve = bytearray([1]) * (hi + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, math.isqrt(hi) + 1):
        if sieve[i]:
            sieve[i * i :: i] = bytearray(len(range(i * i, hi + 1, i)))
    return tuple(i for i in range(max(lo, 2), hi + 1) if sieve[i])


def random_prime(bits: int, rng: random.Random) -> int:
    while True:
        c = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if is_prime(c):
            return c


def fundamental_discriminant(D: int) -> tuple[int, int]:
    """Write a negative discriminant D = f**2 * D_K with D_K fundamental.

    Returns (D_K, f). D must be negative and congruent to 0 or 1 mod 4.
    """
    assert D < 0 and D % 4 in (0, 1), D
    n = -D
    s, m = 1, 1
    for q, e in factor(n).items():
        s *= q ** (e // 2)
        if e % 2:
            m *= q
    # n = s**2 * m with m squarefree.
    if (-m) % 4 == 1:
        D_K, f = -m, s
    else:
        assert s % 2 == 0, (D, s, m)
        D_K, f = -4 * m, s // 2
    assert D_K * f * f == D, (D, D_K, f)
    return D_K, f


def class_number(D: int) -> int:
    """h(D) for a negative discriminant D by counting reduced forms.

    Vectorised over (a, b); fine for |D| up to a few million.
    """
    import numpy as np

    assert D < 0 and D % 4 in (0, 1)
    A = math.isqrt(-D // 3) + 1
    a = np.arange(1, A + 1, dtype=np.int64)
    b = np.arange(-A, A + 1, dtype=np.int64)
    aa, bb = np.meshgrid(a, b, indexing="ij")
    num = bb * bb - D
    ok = (num % (4 * aa) == 0) & (np.abs(bb) <= aa)
    cc = np.where(ok, num // (4 * aa), 0)
    ok &= cc >= aa
    # reduced: |b| <= a <= c, and if |b| == a or a == c then b >= 0
    ok &= ~(((np.abs(bb) == aa) | (aa == cc)) & (bb < 0))
    # primitive forms only
    if ok.any():
        g = np.gcd(np.gcd(aa, np.abs(bb)), cc)
        ok &= g == 1
    return int(ok.sum())
