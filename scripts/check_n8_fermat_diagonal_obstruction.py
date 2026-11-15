#!/usr/bin/env python3
"""Exact regression for an all-depth mathematical E6 bad-prime lemma.

The proof for every k uses Fermat factors and order of 2; this checks k=2..5,
including composite Fermat F_5, without probability or numerical tolerances.
"""
from math import isqrt

def factors(n: int) -> list[int]:
    result = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            result.append(d)
            while n % d == 0:
                n //= d
        d = 3 if d == 2 else d + 2
    if n > 1:
        result.append(n)
    return result

def prime_trial(p: int) -> bool:
    if p < 2:
        return False
    for d in range(2, isqrt(p) + 1):
        if p % d == 0:
            return False
    return True

def main() -> None:
    for k in range(2, 6):
        n = 1 << k
        fermat = (1 << n) + 1
        primes = factors(fermat)
        assert primes
        assert all(prime_trial(p) for p in primes)
        for p in primes:
            assert fermat % p == 0
            assert (p - 1) % (2 * n) == 0, (k, p, "not admissible")
            assert pow(2, n, p) == p - 1
            assert pow(4, n, p) == 1
            assert pow(4, n // 2, p) != 1
            assert (4 * 4 - 4 - 4 - 8) % p == 0
            assert 2187 % p != 0  # E6 y^2=x^3+6x-7 nonsingular.
            print(f"PASS k={k} p={p} ord_p(4)={n}, h(4,4)=0")
    print("PASS all-depth Fermat-factor obstruction: exact k=2..5 regressions")

if __name__ == "__main__":
    main()
