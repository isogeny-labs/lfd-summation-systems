#!/usr/bin/env python3
"""Independent finite-field check of N8's exceptional 2-torsion model.

Validates the exact Semaev identity, the Möbius involution and the
nonvanishing denominator at its second application. These finite tests
DO NOT establish the function-field valuation/spread theorem or certify Lean.
"""
from random import Random


def s3(p: int, a: int, b: int, x: int, z: int, t: int) -> int:
    return (
        (x - z) ** 2 * t**2
        - 2 * ((x + z) * (x * z + a) + 2 * b) * t
        + (x * z - a) ** 2
        - 4 * b * (x + z)
    ) % p


def check() -> int:
    rng = Random(0x3208)
    total = 0
    for p in (3, 5, 7, 11, 13, 17, 31, 101, 131, 257):
        for _ in range(3500):
            t, a = rng.randrange(p), rng.randrange(p)
            b = (-t**3 - a * t) % p
            if (4 * a**3 + 27 * b**2) % p == 0:
                continue
            derivative = (3 * t * t + a) % p
            assert derivative != 0
            x = rng.randrange(p)
            if x == t:
                continue
            z = ((t * x + 2 * t * t + a) * pow((x - t) % p, -1, p)) % p
            assert z != t
            x_recovered = ((t * z + 2 * t * t + a) * pow((z - t) % p, -1, p)) % p
            assert x_recovered == x
            assert s3(p, a, b, x, z, t) == 0
            total += 1
    return total


if __name__ == "__main__":
    n = check()
    print(f"N8 two-torsion Mobius regression: {n} / {n} exact finite-field checks")
