#!/usr/bin/env python3
"""Exact integer-matrix certificates for prime-specialization of bounded closure.

An integer matrix M represents a rational span of bounded-degree rows.
For b outside its Q-column span, select an invertible rank-r minor of M
and a nonzero rank-(r+1) minor of [M|b]. Every prime not dividing their
product preserves nonmembership. For b inside the span, return an
integer identity M*z=d*b; every prime not dividing d preserves membership.

This tool does NOT construct N8's entire closure or prove an LFD theorem.
No floating-point arithmetic, probabilistic rank or modular inference.
"""
from dataclasses import dataclass
from itertools import combinations
from math import gcd
from functools import reduce

from sympy import Matrix, ZZ, ilcm, factorint


def minor_witness(A, size):
    """Return first square nonsingular minor (row indices, column indices, det)."""
    for rows in combinations(range(A.rows), size):
        for cols in combinations(range(A.cols), size):
            delta = int(A.extract(rows, cols).det(method="domain-ge"))
            if delta:
                return rows, cols, delta
    raise ValueError(f"no minor of size {size}")


@dataclass(frozen=True)
class Outside:
    rank: int
    basis_rows: tuple
    basis_cols: tuple
    delta: int
    aug_rows: tuple
    aug_cols: tuple
    gamma: int

    @property
    def excluded_divisor(self):
        return abs(self.delta * self.gamma)


@dataclass(frozen=True)
class Inside:
    integer_vector: tuple
    denominator: int

    @property
    def excluded_divisor(self):
        return self.denominator


def certify(M, b):
    M = Matrix(M)
    b = Matrix(b)
    assert b.cols == 1 and b.rows == M.rows
    assert all(x.q == 1 for x in M) and all(x.q == 1 for x in b)
    rank = M.rank()
    aug = M.row_join(b)
    if aug.rank() > rank:
        if rank:
            rows, cols, delta = minor_witness(M, rank)
        else:
            rows, cols, delta = (), (), 1
        arows, acols, gamma = minor_witness(aug, rank + 1)
        cert = Outside(rank, rows, cols, delta, arows, acols, gamma)
    else:
        sol, parameters = M.gauss_jordan_solve(b)
        if parameters.rows:
            sol = sol.subs({p: 0 for p in parameters})
        d = reduce(ilcm, [int(x.q) for x in sol], 1)
        z = tuple(int(d*x) for x in sol)
        cert = Inside(z, int(d))
    verify(M, b, cert)
    return cert


def verify(M, b, cert):
    M, b = Matrix(M), Matrix(b)
    if isinstance(cert, Inside):
        assert cert.denominator > 0
        assert M * Matrix(cert.integer_vector) == cert.denominator * b
    else:
        assert cert.rank == M.rank()
        assert int(M.extract(cert.basis_rows, cert.basis_cols).det()) == cert.delta
        assert cert.delta != 0
        aug = M.row_join(b)
        assert int(aug.extract(cert.aug_rows, cert.aug_cols).det()) == cert.gamma
        assert cert.gamma != 0
        assert aug.rank() == cert.rank + 1


def modular_membership(M, b, p):
    from sympy.polys.domains import GF
    from sympy.polys.matrices import DomainMatrix
    M, b = Matrix(M), Matrix(b)
    rows = [[int(M[i,j]) % p for j in range(M.cols)] for i in range(M.rows)]
    augmented = [row + [int(b[i,0]) % p] for i,row in enumerate(rows)]
    def rank(x):
        # DomainMatrix over an exact finite field, never symbolic mod-p simplification.
        if not x or not x[0]:
            return 0
        return DomainMatrix.from_list(x, GF(p)).rank()
    return rank(rows) == rank(augmented)


def selftest():
    cases = [
        (Matrix([[2,0],[0,3]]), Matrix([1,1]), True),
        (Matrix([[2],[4],[6]]), Matrix([1,2,3]), True),
        (Matrix([[2],[4],[6]]), Matrix([1,2,4]), False),
        (Matrix([[0,0],[0,0]]), Matrix([0,1]), False),
        (Matrix([[0,0],[0,0]]), Matrix([0,0]), True),
        (Matrix([[2,4],[3,6],[5,10]]), Matrix([1,0,0]), False),
    ]
    primes = [2,3,5,7,11,13,17,19,23,31]
    count = 0
    for M,b,expect in cases:
        cert = certify(M,b)
        assert isinstance(cert, Inside) == expect
        for p in primes:
            if cert.excluded_divisor % p == 0:
                continue
            assert modular_membership(M,b,p) == expect, (M,b,p,cert)
            count += 1
    # A prime outside the certificate's exceptional set preserves
    # the true rational membership decision, never merely a guessed rank.
    print("PASS integer rank/minor specialization:", len(cases),
          "certificates,", count, "good-prime checks")
    print("Representative obstruction:", certify(Matrix([[2],[4],[6]]),Matrix([1,2,4])))


if __name__ == "__main__":
    selftest()
