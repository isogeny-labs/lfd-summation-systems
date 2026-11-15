#!/usr/bin/env python3
"""Independent SymPy checks of the all-size monomial Wronskian determinant."""
import math
import sympy as sp

X=sp.symbols("X")
cases=0
for r in range(0,11):
    polynomials=[X**j for j in range(r)]
    if r:
        W=sp.Matrix([[sp.diff(p, X, i) for p in polynomials] for i in range(r)])
        determinant=sp.factor(W.det(method="berkowitz"))
    else:
        determinant=sp.Integer(1)
    expected=math.prod(math.factorial(j) for j in range(r))
    assert determinant == expected, (r,determinant,expected)
    cases+=1
print("N8 exact monomial Wronskian determinants passed:",cases,"(r=0..10)")
