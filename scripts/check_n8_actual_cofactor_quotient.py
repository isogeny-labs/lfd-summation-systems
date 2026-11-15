#!/usr/bin/env python3
"""Independent exact N8 genuine polynomial-cofactor quotient regression.

The relation x*z = x + z + c, z=(x+c)/(x-1), gives
z = -c*1 - x + x*z. For each omitted column check both polynomial
cofactors, the nonzero denominator, and the derived rational quotient.
"""
from fractions import Fraction
from itertools import product

import sympy as sp

x = sp.symbols("x")
base = ((0, 0), (1, 0), (1, 1))
replacement = (0, 1)
parameters = (sp.Integer(2), sp.Integer(6), sp.Integer(8),
              sp.Rational(5, 2), sp.Integer(-9), sp.Rational(17, 3))
checks = 0

for c in parameters:
    D = x * (x + c) * (x - 1)

    def monomial(e):
        u, v = e
        return x ** u * ((x + c) / (x - 1)) ** v

    def numerator_rows(e, r):
        u, v = e
        Q = u * (x + c) * (x - 1) - v * (c + 1) * x
        P = sp.Integer(1)
        result = []
        for s in range(r):
            result.append(P)
            P = sp.expand(D * sp.diff(P, x) + (Q - s * sp.diff(D, x)) * P)
        return result

    def cofactor(exponents):
        r = len(exponents)
        columns = [numerator_rows(e, r) for e in exponents]
        return sp.factor(sp.Matrix(r, r, lambda s, j: columns[j][s]).det())

    A = cofactor(base)
    assert A != 0
    assert sp.simplify(
        monomial(replacement) - sum(
            coeff * monomial(e) for coeff, e in zip((-c, -1, 1), base)
        )) == 0

    for i, gamma in enumerate((-c, -1, 1)):
        updated = list(base)
        updated[i] = replacement
        B = cofactor(updated)
        assert B != 0
        assert gamma != 0
        assert sp.cancel(B * monomial(replacement) -
                         gamma * A * monomial(base[i])) == 0
        assert sp.cancel(monomial(replacement) / monomial(base[i])
                         - gamma * A / B) == 0
        checks += 2

assert checks == 36
print(f"N8 actual polynomial cofactor cross and quotient: {checks} exact identities passed")
