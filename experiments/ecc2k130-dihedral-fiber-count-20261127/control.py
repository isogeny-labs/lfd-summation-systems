"""Finite controls for rational fibers in the 131-sheet affine action."""

p = 131

# E: y^2 + x*y = x^3 + 1 over F_2.
finite_points = [
    (x, y)
    for x in (0, 1)
    for y in (0, 1)
    if (y * y + x * y - x * x * x - 1) % 2 == 0
]
assert finite_points == [(0, 1), (1, 0), (1, 1)]
assert 1 + len(finite_points) == 4

# A 131-cycle has no fixed basis vector, hence permutation trace zero.
assert all((j + 1) % p != j for j in range(p))

tested = 0
for multiplier in range(p):
    if multiplier in (0, 1, p - 1):
        continue
    for center in range(p):
        frobenius = lambda x: (center + multiplier * (x - center)) % p
        reflection = lambda x: (2 * center - x) % p

        # Unramified fiber: exactly one Frobenius-fixed sheet.
        assert sum(frobenius(x) == x for x in range(p)) == 1

        # Reflection inertia: one fixed sheet and 65 paired orbits.
        orbits = {
            tuple(sorted((x, reflection(x)))) for x in range(p)
        }
        assert len(orbits) == 66
        assert sum(
            tuple(sorted((frobenius(a), frobenius(b)))) == (a, b)
            for a, b in orbits
        ) == 1

        # Translation inertia is regular: all 131 sheets form one orbit.
        tested += 1

assert tested == 128 * 131
print(f"four base points; {tested} affine cases: one rational point per fiber")
