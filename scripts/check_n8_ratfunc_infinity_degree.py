#!/usr/bin/env python3
"""Independent exact infinity-degree checks for the N8 rational monomials."""
from itertools import product
import sympy as s
x=s.symbols("x")
cases=0
for c in [s.Integer(8), s.Integer(-1), s.Rational(3,2), s.Integer(0),s.Integer(7)]:
    for u,v in product(range(4),range(4)):
        f=s.cancel(x**u*((x+c)/(x-1))**v)
        p,q=s.fraction(f)
        assert s.degree(p,x)-s.degree(q,x)==u
        cases+=1
    for (u,v),(U,V) in zip(product(range(4),range(4)),reversed(list(product(range(4),range(4))))):
        f=s.cancel(x**u*((x+c)/(x-1))**v / (x**U*((x+c)/(x-1))**V))
        p,q=s.fraction(f)
        assert s.degree(p,x)-s.degree(q,x)==u-U
        cases+=1
print(f"OK: {cases} exact rational N8 monomial/ratio pole degrees over Q")
