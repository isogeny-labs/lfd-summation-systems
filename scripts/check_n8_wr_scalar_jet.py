#!/usr/bin/env python3
"""Exact independent tests for the common-scalar Wronskian identity.

Lean on this PR proves only the triangular jet determinant. This script
separately tests the full prospective determinant identity over Q(X).
"""
import sympy as sp

X = sp.symbols("X")
families = [
    [sp.Integer(1), X, X**2],
    [X, (X+8)/(X-1), X**2 + 3],
]
scalars = [X+1, (X+8)/(X-1)]
checked = 0
for ff in families:
    for h in scalars:
        for r in range(1,4):
            f = ff[:r]
            wr = sp.det(sp.Matrix([[sp.diff(z, X, i) for z in f] for i in range(r)]))
            wrh = sp.det(sp.Matrix([[sp.diff(h*z, X, i) for z in f] for i in range(r)]))
            assert sp.cancel(wrh - h**r * wr) == 0, (f,h,r)
            checked += 1
print(f"N8 scalar Wronskian: {checked} exact determinant identities passed")
