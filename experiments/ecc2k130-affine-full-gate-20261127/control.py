"""Finite controls for the full affine degree-131 cover exclusion."""

from math import gcd

p = 131
divisors = tuple(d for d in range(2, p) if (p - 1) % d == 0)
assert divisors == (2, 5, 10, 13, 26, 65, 130)

# E: y^2 + x*y = x^3 + 1 over F_2.
base = [(x, y) for x in (0, 1) for y in (0, 1)
        if (y * y + x * y - x * x * x - 1) % 2 == 0]
assert base == [(0, 1), (1, 0), (1, 1)]
assert len(base) + 1 == 4

for ell, expected_order, minimum_codimension in ((5, 4, 104), (13, 12, 120)):
    order = next(k for k in range(1, ell) if pow(2, k, ell) == 1)
    assert order == expected_order
    assert (p - 1) * (ell - 1) // ell == minimum_codimension
    assert order * minimum_codimension > 2 * p - 2

checked = 0
for d in divisors:
    H = {pow(2, (p - 1) // d * k, p) for k in range(d)}
    assert len(H) == d
    for u in range(1, p):
        if u in H:
            continue
        assert sum((u * x) % p == x for x in range(p)) == 1
        for e in range(1, d + 1):
            if d % e:
                continue
            I = {pow(2, (p - 1) // e * k, p) for k in range(e)}
            assert I <= H
            orbits = {tuple(sorted((a * x) % p for a in I)) for x in range(p)}
            assert len(orbits) == 1 + (p - 1) // e
            assert sum(tuple(sorted((u * x) % p for x in orbit)) == orbit
                       for orbit in orbits) == 1
            checked += 1

# The only nontrivial divisor of 130 without a 5- or 13-factor is 2.
assert tuple(d for d in divisors if gcd(d, 5 * 13) == 1) == (2,)
print(f"four base points; {checked} affine multiplier/inertia cases; "
      "order-5 and order-13 conductor minima 416 and 1440 exceed 260")
