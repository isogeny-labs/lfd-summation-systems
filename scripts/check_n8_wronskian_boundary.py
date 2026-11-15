"""Exact characteristic-p Wronskian tests at maximal circuit size.

These independent computations exercise p+1-element circuits over F_p(x^p)
and nonconstant p-th-power coefficient multipliers. They are NOT proofs
of the N8 function-field height bound or of any Lean declaration.
"""
from __future__ import annotations
from math import factorial
import json
import random


def determinant_mod(A, p):
    A = [[z % p for z in row] for row in A]
    n = len(A)
    determinant = 1
    for j in range(n):
        i = next((i for i in range(j,n) if A[i][j]), None)
        if i is None:
            return 0
        if i != j:
            A[j], A[i] = A[i], A[j]
            determinant = -determinant
        z = A[j][j]
        determinant = determinant*z % p
        inv = pow(z,-1,p)
        for i in range(j+1,n):
            mult = A[i][j]*inv % p
            for k in range(j,n):
                A[i][k] = (A[i][k]-mult*A[j][k]) % p
    return determinant % p


def derivative_at_one(coefficients, order, p):
    return sum(c * factorial(i)//factorial(i-order)
               for i,c in enumerate(coefficients) if i>=order) % p


def twist(coefficients, e, p):
    """Multiply by 1 + e*x^p, a scalar from F_p(x^p)."""
    out = [0]*(len(coefficients)+p)
    for i,c in enumerate(coefficients):
        out[i] = c
        out[i+p] = e*c
    return out


def run(p, repeats=60):
    # p+1 members: 1,x,...,x^(p-1), and sum(x^i). Every p-subset
    # independent over F_p(x^p); all p+1 are dependent.
    columns = [[int(j==i) for j in range(p)] for i in range(p)]
    columns.append([1]*p)
    for omit in range(p+1):
        selected = [c for j,c in enumerate(columns) if j!=omit]
        assert determinant_mod([[selected[j][i] for j in range(p)]
                                 for i in range(p)],p) != 0
    rng = random.Random(81307+p)
    nonconstant=0
    for _ in range(repeats):
        cols=[[rng.randrange(p) for _ in range(p)] for _ in range(p)]
        es=[rng.randrange(p-1) for _ in range(p)]
        # Ordinary derivatives of f(x^p) are zero; the Wronskian
        # of the twisted family is the untwisted Wronskian
        # times product_i (1+e_i x^p).
        polys=[twist(c,e,p) for c,e in zip(cols,es)]
        matrix=[[derivative_at_one(f,i,p) for f in polys] for i in range(p)]
        actual=determinant_mod(matrix,p)
        C=[[cols[j][i] for j in range(p)] for i in range(p)]
        expected=determinant_mod(C,p)
        for j in range(p):
            expected=expected*factorial(j)*(1+es[j]) % p
        assert actual==expected,(p,actual,expected)
        nonconstant+=int(any(es))
    return dict(characteristic=p,circuit_size=p+1,all_p_minors=p+1,
                twisted_wronskians=repeats,nonconstant_trials=nonconstant)


if __name__=="__main__":
    print(json.dumps({str(p):run(p) for p in (3,5,7,11)},indent=2))
