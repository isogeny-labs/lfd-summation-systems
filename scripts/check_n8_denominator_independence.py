#!/usr/bin/env python3
"""Independent exact rank regressions for rational-to-polynomial basis transfer."""
import sympy as sp
x=sp.Symbol("x")
families=[
 [1/(x-1), x/(x+2), x*x/(x+3)],
 [(x+8)/(x-1), x*(x+8)/(x-1), x*x*(x+8)/(x-1)],
 [1/(x-1), x/(x-1), (x+1)/(x-1)],
 [x/(x*x+1), x*x/(x+1), 1/(x-2)],
]
count=0
for fs in families:
 nd=[sp.fraction(sp.cancel(v)) for v in fs]
 H=sp.prod(d for _,d in nd)
 p=[sp.cancel(H*v) for v in fs]
 assert all(sp.denom(y)==1 for y in p)
 # Independence over Q can be tested on the coefficient matrix of
 # the corresponding polynomials: multiplying every member by the
 # same nonzero rational H cannot change the constant-linear rank.
 def rank(f):
  num=[sp.together(q).as_numer_denom()[0] for q in f]
  den=[sp.together(q).as_numer_denom()[1] for q in f]
  common=sp.prod(den)
  polys=[sp.Poly(sp.cancel(common*q),x) for q in f]
  deg=max(z.degree() for z in polys if not z.is_zero)
  M=sp.Matrix([[z.nth(i) for z in polys] for i in range(deg+1)])
  return M.rank()
 assert rank(fs)==rank(p)
 count+=1
print(f"N8 independent numerator transfer: {count} exact rational/poly ranks matched")
