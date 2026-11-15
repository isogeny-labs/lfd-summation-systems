#!/usr/bin/env python3
"""Exact independent checks of N8 rational derivatives and the tight 2s degree bound.

This tests polynomial identities over Q[X], not a general theorem about
characteristic-zero Wronskians or GrowthSpread.SpreadBound. No floating point.
"""
from fractions import Fraction as F

from check_n8_explicit_wronskian import (
    add, sub, mul, scale, poly, power, diff,
)


def degree(f):
    return len(f) - 1


def check(c, u, v, max_order=4):
    x = poly(0, 1)
    xc = poly(c, 1)
    xm1 = poly(-1, 1)
    D = mul(mul(x, xc), xm1)
    Q = sub(scale(mul(xc, xm1), u), scale(x, v * (c + 1)))
    N = mul(power(x, u), power(xc, v))
    R = power(xm1, v)

    # Exactly D * (f'/f numerator) = Q * N * R for f=N/R.
    assert mul(D, sub(mul(diff(N), R), mul(N, diff(R)))) == mul(Q, mul(N, R))

    P = [poly(1)]
    for s in range(max_order + 1):
        assert degree(P[s]) <= 2 * s, (c, u, v, s, degree(P[s]))
        Qminus = sub(Q, scale(diff(D), s))
        Pnext = add(mul(D, diff(P[s])), mul(Qminus, P[s]))

        # Differentiate N P_s/(R D^s) by the quotient rule,
        # and compare exactly with N P_(s+1)/(R D^(s+1)).
        numer = mul(N, P[s])
        denom = mul(R, power(D, s))
        numerator_derivative = sub(mul(diff(numer), denom),
                                   mul(numer, diff(denom)))
        lhs = mul(numerator_derivative, mul(R, power(D, s + 1)))
        rhs = mul(mul(N, Pnext), power(denom, 2))
        assert lhs == rhs, (c, u, v, s)

        P.append(Pnext)


if __name__ == "__main__":
    cases = [F(5, 2), F(6), F(8), F(-9), F(17, 3)]
    for c in cases:
        for u in range(5):
            for v in range(5):
                check(c, u, v)
    print("PASS: 125 exact rational logarithmic derivative identities")
    print("PASS: 625 exact higher derivative identities and 625 tight degree ceilings")
