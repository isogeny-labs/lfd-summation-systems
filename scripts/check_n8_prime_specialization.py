#!/usr/bin/env python3
"""Exact depth-k=2 prime-specialization regression for the N8 fixed E6 curve.

The two-torsion S3 relation is h(x,z)^2 with
h=(x-1)z-x-8. Work in Z[x,z]/(x^4-1,z^4-1).
The multiplication matrix of h is integral and has a nonzero determinant.
A prime dividing the determinant may destroy the unit-ideal certificate.
No numerical approximation is used.
"""
from sympy import Matrix, factorint

N = 4
TERMS = {(1, 1): 1, (0, 1): -1, (1, 0): -1, (0, 0): -8}

def multiplication_matrix(n=N):
    mat = [[0 for _ in range(n*n)] for _ in range(n*n)]
    for i in range(n):
        for j in range(n):
            column = i*n + j
            for (a,b), coeff in TERMS.items():
                row = ((i+a)%n)*n + ((j+b)%n)
                mat[row][column] += coeff
    return Matrix(mat)

def root_pairs(p, n=N):
    roots = [a for a in range(1,p) if pow(a,n,p)==1]
    return [(x,z) for x in roots for z in roots
            if ((x-1)*z-x-8)%p == 0]

def main():
    det = int(multiplication_matrix().det())
    assert det == 279791370354825, det
    factors = factorint(abs(det))
    assert factors == {3:14, 5:2, 7:2, 17:1, 53:2}, factors
    assert root_pairs(17) == [(4,4)]
    assert root_pairs(53) == [(23,52),(52,23)]
    for p in (97,113,193,257):
        assert not root_pairs(p), (p,root_pairs(p))
    assert (4**4-1)%17 == 0
    assert (4*6**3+27*(-7)**2)%17 != 0
    assert (17-1)%8 == 0
    print("PASS N8 E6: k=2 exact grid determinant", det)
    print("factorization:", factors)
    print("admissible counterexample: p=17, x=z=4, 8|(p-1), nonsingular")
    print("good sample primes: 97, 113, 193, 257")

if __name__=="__main__":
    main()
