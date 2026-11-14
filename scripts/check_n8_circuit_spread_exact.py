#!/usr/bin/env python3
"""Independent exact finite-circuit tests for the N8 rational section.

For c=2,8,-9 and exponents (u,v) in {0,1,2}^2, clear the common
(X-1)^V denominator in Q(X), enumerate ALL dependent minimal subsets
of size 2..6, and independently verify:
  * the sole dependency has nonzero coefficients at every column;
  * deleting any one column makes the family independent;
  * every pair satisfies the sharper 4(m-1)(m-2) exponent bound.

This finite regression is not a proof of GrowthSpread.SpreadBound.
"""
from functools import lru_cache
from itertools import combinations, product

import sympy as sp

X = sp.Symbol("X")
EXPONENTS = tuple(product(range(3), repeat=2))


def check_parameter(c):
    @lru_cache(None)
    def matrix(indices):
        vectors = [EXPONENTS[j] for j in indices]
        v_max = max((v for _, v in vectors), default=0)
        polynomials = [
            sp.Poly(
                sp.expand(X ** u * (X + c) ** v * (X - 1) ** (v_max - v)),
                X,
                domain=sp.QQ,
            )
            for u, v in vectors
        ]
        degree = max((p.degree() for p in polynomials), default=-1)
        return sp.Matrix([
            [p.nth(i) for p in polynomials] for i in range(degree + 1)
        ])

    @lru_cache(None)
    def rank(indices):
        return matrix(indices).rank()

    by_size = {}
    for m in range(2, 7):
        count = 0
        for indices in combinations(range(len(EXPONENTS)), m):
            if rank(indices) == m:
                continue
            if not all(rank(tuple(j for j in indices if j != omitted)) == m-1
                       for omitted in indices):
                continue
            M = matrix(indices)
            assert M.rank() == m - 1
            dependencies = M.nullspace()
            assert len(dependencies) == 1
            assert all(value != 0 for value in dependencies[0])
            for i, j in combinations(indices, 2):
                u, v = EXPONENTS[i]
                U, V = EXPONENTS[j]
                assert 2 * (abs(u - U) + abs(v - V)) <= 4 * (m - 1) * (m - 2)
            count += 1
        by_size[m] = count
    assert by_size == {2: 0, 3: 0, 4: 4, 5: 9, 6: 20}, (c, by_size)
    return sum(by_size.values())


if __name__ == "__main__":
    total = sum(check_parameter(c) for c in (2, 8, -9))
    assert total == 99
    print(f"N8 independent exact minimal circuits: {total} passed")
