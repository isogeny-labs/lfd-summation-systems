#!/usr/bin/env python3
"""Independent exact Q(X) test of arbitrary-row Wronskian determinant normalization.

The genuine derivative rows are built by repeated quotient-rule differentiation
of rational fractions; the polynomial cofactor rows are computed *separately*
by the polynomial logarithmic recurrence. This finite experiment checks the
actual determinant factorization for five rational parameters and four
small indexed families, not arbitrary-size Wronskian nonvanishing or SpreadBound.
"""
from fractions import Fraction as F
from itertools import permutations

from check_n8_explicit_wronskian import poly, add, sub, neg, mul, scale, diff, power

ZERO = (poly(0), poly(1))
ONE = (poly(1), poly(1))


def qadd(a, b):
    return add(mul(a[0], b[1]), mul(b[0], a[1])), mul(a[1], b[1])


def qmul(a, b):
    return mul(a[0], b[0]), mul(a[1], b[1])


def qder(a):
    return sub(mul(diff(a[0]), a[1]), mul(a[0], diff(a[1]))), power(a[1], 2)


def qpow(a, n):
    result = ONE
    for _ in range(n):
        result = qmul(result, a)
    return result


def qeq(a, b):
    return mul(a[0], b[1]) == mul(a[1], b[0])


def qdet(rows):
    result = ZERO
    for sigma in permutations(range(len(rows))):
        term = ONE
        for i, j in enumerate(sigma):
            term = qmul(term, rows[i][j])
        inversions = sum(sigma[i] > sigma[j] for i in range(len(sigma))
                         for j in range(i + 1, len(sigma)))
        if inversions % 2:
            term = neg(term[0]), term[1]
        result = qadd(result, term)
    return result


def check(c, exponents):
    x, one = poly(0, 1), poly(1)
    xm1, xc = poly(-1, 1), poly(c, 1)
    D = mul(mul(x, xc), xm1)
    r = len(exponents)
    functions = []
    numerators = []
    for u, v in exponents:
        functions.append((mul(power(x, u), power(xc, v)), power(xm1, v)))
        Q = sub(scale(mul(xc, xm1), F(u)), scale(x, F(v) * (c + 1)))
        P = [one]
        for i in range(r - 1):
            P.append(add(mul(D, diff(P[-1])),
                         mul(sub(Q, scale(diff(D), i)), P[-1])))
        numerators.append(P)

    # Independent direct quotient-rule derivatives.
    W = []
    derivatives = list(functions)
    for i in range(r):
        W.append(list(derivatives))
        derivatives = [qder(f) for f in derivatives]

    # Check every scaled row/column equality, then the full determinant.
    d = (D, one)
    for i in range(r):
        for j in range(r):
            assert qeq(qmul(W[i][j], qpow(d, i)),
                       qmul(functions[j], (numerators[j][i], one))), (c, exponents, i, j)

    lhs = qmul((power(D, r * (r - 1) // 2), one), qdet(W))
    cofactor = qdet([[(numerators[j][i], one) for j in range(r)] for i in range(r)])
    factors = ONE
    for f in functions:
        factors = qmul(factors, f)
    assert qeq(lhs, qmul(cofactor, factors)), (c, exponents)


if __name__ == "__main__":
    parameters = [F(5, 2), F(6), F(8), F(-9), F(17, 3)]
    families = [
        [(0, 0)],
        [(1, 0), (0, 1)],
        [(0, 0), (1, 0), (0, 1)],
        [(1, 1), (1, 0), (0, 1)],
    ]
    for c in parameters:
        for e in families:
            check(c, e)
    print("PASS: 20 exact independent rational Wronskian determinant normalizations (r=1,2,3)")
