#!/usr/bin/env python3
"""Exact Q[X] N8 four-term Wronskian cofactor regression at five parameters.

Finite exact tests only: they do not establish nonvanishing for arbitrary
circuits, SpreadBound, CurveModel or unconditional LFD growth.
"""
from fractions import Fraction as F
from check_n8_explicit_wronskian import poly, add, sub, scale, mul, power, diff, det


def check(c):
    assert c not in (F(-1), F(0))
    x, xm1, xc = poly(0, 1), poly(-1, 1), poly(c, 1)
    D = mul(mul(x, xm1), xc)
    exponents = [(1, 1), (1, 0), (0, 1), (0, 0)]
    coefficients = [F(1), F(-1), F(-1), -c]
    r = len(exponents) - 1
    P = []
    for u, v in exponents:
        Q = sub(scale(mul(xc, xm1), u), scale(x, v * (c + 1)))
        terms = [poly(1)]
        for j in range(r - 1):
            terms.append(add(mul(D, diff(terms[-1])),
                             mul(sub(Q, scale(diff(D), j)), terms[-1])))
        assert all(len(terms[j]) - 1 <= 3 * j for j in range(r))
        P.append(terms)
    cofactors = []
    for omitted in range(len(exponents)):
        A = det([[P[i][j] for j in range(r)]
                 for i in range(4) if i != omitted])
        assert A != poly(0) and len(A) - 1 <= 3 * r * (r - 1) // 2
        cofactors.append(A)
    for j, (u, v) in enumerate(exponents):
        u0, v0 = exponents[0]
        fjnum = mul(power(x, u), power(xc, v))
        fjden = power(xm1, v)
        f0num = mul(power(x, u0), power(xc, v0))
        f0den = power(xm1, v0)
        lhs = scale(mul(mul(f0num, fjden), cofactors[j]), coefficients[0])
        rhs = scale(mul(mul(fjnum, f0den), cofactors[0]),
                    (-1) ** j * coefficients[j])
        assert lhs == rhs, (c, j)
    return [len(q) - 1 for q in cofactors]


if __name__ == "__main__":
    for c in [F(5, 2), F(6), F(8), F(-9), F(17, 3)]:
        degrees = check(c)
        assert degrees == [5, 5, 4, 4], (c, degrees)
        print("PASS", c, degrees, "4-term circuit and exact cofactor identities")
