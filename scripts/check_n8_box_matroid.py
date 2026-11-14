#!/usr/bin/env python3
"""Exact bounded exponent-box matroid regression for the E6 two-torsion section.

Finite tests are not a proof of specialization for all but finitely many p.
"""
from math import comb
from random import Random
from sympy import Matrix, ZZ
from sympy.polys.domains import GF
from sympy.polys.matrices import DomainMatrix


def integer_matroid(N, a=6):
    labels = [(u, v) for u in range(N) for v in range(N)]
    columns = []
    for u, v in labels:
        col = [0] * (2 * N - 1)
        for i in range(v + 1):
            for j in range(N - v):
                col[u + i + j] += (comb(v, i) * (a + 2) ** (v - i)
                                  * comb(N - 1 - v, j) * (-1) ** (N - 1 - v - j))
        columns.append(col)
    return labels, [list(row) for row in zip(*columns)]


def rank(A, cols, p=None):
    if not cols:
        return 0
    selected = [[row[j] for j in cols] for row in A]
    return DomainMatrix.from_list(selected, ZZ if p is None else GF(p)).rank()


def test():
    for N in (2, 4):
        labels, A = integer_matroid(N)
        all_cols = tuple(range(N * N))
        assert rank(A, all_cols) == 2 * N - 1
        assert rank(A, all_cols, 257) == 2 * N - 1
        rng = Random(N + 20261113)
        for i in range(200):
            size = rng.randint(0, N * N)
            subset = tuple(sorted(rng.sample(all_cols, size)))
            assert rank(A, subset) == rank(A, subset, 257), (N, i, subset)
    labels, A = integer_matroid(4)
    witness = [(0, 1), (0, 3), (1, 0), (1, 1), (2, 3), (3, 2), (3, 3)]
    idx = [labels.index(e) for e in witness]
    det = int(Matrix([[row[j] for j in idx] for row in A]).det())
    assert det == 578207808 == 2 ** 6 * 3 ** 12 * 17
    assert rank(A, idx) == 7 and rank(A, idx, 17) == 6
    assert rank(A, idx, 257) == 7
    print("PASS bounded-box matroid: k=1,2; 400 matched p=257 subsets")
    print("PASS p=17 rank-drop witness: det=578207808=2^6*3^12*17, rank 7 -> 6")


if __name__ == "__main__":
    test()
