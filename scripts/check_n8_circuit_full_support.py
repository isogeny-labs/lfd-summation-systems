#!/usr/bin/env python3
"""Independent exact minimal-circuit regression on the fixed E_6 rational model.

For z=(x+8)/(x-1), the four distinct monomials 1,x,z,xz obey
xz-z-x-8=0, and each 3-element subfamily is independent over Q.
"""
from itertools import combinations
import sympy as s

x=s.symbols("x")
z=(x+8)/(x-1)
f=[s.Integer(1),x,z,x*z]
co=[-8,-1,-1,1]
assert s.cancel(sum(co[i]*f[i] for i in range(4)))==0
assert all(c!=0 for c in co)
p=[s.expand((x-1)*q) for q in f]
def vec(q):
    P=s.Poly(q,x)
    return [P.nth(i) for i in range(3)]
M=s.Matrix.hstack(*(s.Matrix(vec(q)) for q in p))
assert M.rank()==3
for ids in combinations(range(4),3):
    assert M[:,list(ids)].det()!=0, ids
print("OK: fixed E6 four-term circuit has full nonzero support; all four 3-subsets independent")
