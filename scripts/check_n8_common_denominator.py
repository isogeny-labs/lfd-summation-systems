#!/usr/bin/env python3
"""Exact independent regression for the common-denominator Wronskian bridge."""
import sympy as sp
from itertools import combinations
x = sp.Symbol("x")
families = [
 [x/(x-1), (x+8)/(x-1), 1/(x+2)],
 [1/(x*x+1), x/(x+1), (x*x+2)/(x*x-3)],
 [x*x/(x-1), (x+8)**2/(x-1)**2, (x+5)/(x+3)],
 [sp.Integer(0), x/(x+1), sp.Integer(1)],
]
tests=0
for fs in families:
 for r in range(1,4):
  for indices in combinations(range(3),r):
   chosen = [sp.cancel(fs[j]) for j in indices]
   num_den = [sp.fraction(q) for q in chosen]
   H = sp.prod(d for _,d in num_den)
   assert H != 0
   poly = [sp.cancel(H*q) for q in chosen]
   assert all(sp.denom(v)==1 for v in poly)
   Wrat = sp.det(sp.Matrix([[sp.diff(q,x,i) for q in chosen] for i in range(r)]))
   Wpoly = sp.det(sp.Matrix([[sp.diff(q,x,i) for q in poly] for i in range(r)]))
   assert sp.cancel(H**r*Wrat-Wpoly)==0
   tests += 1
print(f"N8 common polynomial denominator: {tests} exact Wronskian identities passed")
