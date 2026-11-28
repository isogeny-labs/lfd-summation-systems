"""Finite controls for the constant dihedral degree-131 ramification gate."""

p = 131
assert all(pow(2, d, p) != 1 for d in (1, 2, 5, 10, 13, 26, 65))
assert pow(2, 130, p) == 1
assert pow(2, 65, p) == p - 1
assert all(pow(2, d, p) not in (1, p - 1) for d in range(1, 65))

def compose(a, b):
    return ((a[0] + a[1] * b[0]) % p, a[1] * b[1])

def conjugate(g, h):
    inverse = (-g[1] * g[0] % p, g[1])
    return compose(compose(g, h), inverse)

reflection = (0, -1)
rotations = [(a, 1) for a in range(p)]
assert any(conjugate(g, reflection) != reflection for g in rotations)
assert len({(x + 1) % p for x in range(p)}) == p

conductor = p - 1
min_degree = 65
budget = 2 * p - 2
assert conductor * min_degree == 8450 > budget == 260
print(f"ord_131(2)=130; 2^65=-1; minimum conductor degree=8450 > {budget}")
