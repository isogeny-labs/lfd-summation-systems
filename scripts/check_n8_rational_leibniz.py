#!/usr/bin/env python3
"""Independent exact SymPy checks for genuine all-order rational Leibniz rule.

This does not replace the pinned Lean proof. It tests exact nontrivial
rational numerators/denominators, including two-torsion monomials.
"""
from math import comb

import sympy as sp

X = sp.symbols("X")
rational = [
    sp.Integer(0),
    sp.Integer(1),
    X,
    X**2 + 3*X + 2,
    (X + 8) / (X - 1),
    X**2 * (X + 8) / (X - 1),
    (X - 2) / (X**2 + 5),
    (X**3 + 3*X + 1) / ((X - 1) * (X + 8)),
]
count = 0
for f in rational:
    for g in rational:
        for n in range(4):
            lhs = sp.diff(f * g, X, n)
            rhs = sum(
                comb(n, k) * sp.diff(f, X, n - k) * sp.diff(g, X, k)
                for k in range(n + 1)
            )
            d = sp.cancel(lhs - rhs)
            assert d == 0, (f, g, n, d)
            count += 1

print(f"N8 higher rational Leibniz: {count} exact identities passed")
