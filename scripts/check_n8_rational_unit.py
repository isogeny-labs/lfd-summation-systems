#!/usr/bin/env python3
"""Exact finite-field checks for the rational N8 unit-ideal family.

Uses only the standard library. The proof (in docs/proofs/N8-rational-unit-family.md)
is symbolic and applies to every k. This script merely checks examples.
"""
from math import isqrt


def prime(n):
    if n < 2:
        return False
    return all(n % d for d in range(2, isqrt(n) + 1))


def prime_above_2n(n):
    p = 2 * n + 1
    while not (p % n == 1 and prime(p)):
        p += n
    return p


def mul(u, v, p, d):
    a, b = u
    c, e = v
    return ((a * c + d * b * e) % p, (a * e + b * c) % p)


def pow2(u, exponent, p, d):
    result = (1, 0)
    while exponent:
        if exponent & 1:
            result = mul(result, u, p, d)
        u = mul(u, u, p, d)
        exponent >>= 1
    return result


def test(k):
    n = 1 << k
    p = prime_above_2n(n)
    d = next(a for a in range(2, p) if pow(a, (p - 1) // 2, p) == p - 1)
    # K=F_p[alpha]/(alpha^2-d), a=alpha, t=1, b=-1-alpha.
    roots = [x for x in range(1, p) if pow(x, n, p) == 1]
    assert len(roots) == n and (p - 1) % n == 0 and p > 2 * n
    alpha = (0, 1)
    # disc(E_a) = (a+3)^2 (4a+3) != 0 in K.
    disc = mul(mul((3, 1), (3, 1), p, d), (3, 4), p, d)
    assert disc != (0, 0)
    # a=alpha, a^p=-alpha, c=a-a^p=2alpha and c^-1=alpha/(2d).
    c_inverse = (0, pow(2 * d, -1, p))
    for x in roots:
        for z in roots:
            # u=q-alpha, q=(x-1)z-x-2 in Fp; S3=u^2.
            u = (((x - 1) * z - x - 2) % p, p - 1)
            square = mul(u, u, p, d)
            assert square != (0, 0)
            # Explicit ideal certificate: u*(u^(p-1)-1)/c = 1.
            vraw = pow2(u, p - 1, p, d)
            v = mul(((vraw[0] - 1) % p, vraw[1]), c_inverse, p, d)
            assert mul(u, v, p, d) == (1, 0)
            assert mul(square, mul(v, v, p, d), p, d) == (1, 0)
    return k, n, p, n * n


def main():
    for k in range(2, 9):
        k, n, p, checked = test(k)
        print(f"k={k} N={n} p={p} checked={checked} nonsingular=1 no-roots=1")
    print("PASS: seven finite-field instances, 87,376 tested root pairs")


if __name__ == "__main__":
    main()
