#!/usr/bin/env python3
"""Exact independent N8 rational-monomial cofactor Wronskian tests."""
import random
import sympy as s

x = s.symbols("x")
D = x*(x+8)*(x-1)
rng = random.Random(140611)
checks=0
for r in range(1,5):
    trials=0
    while trials < 7:
        es = rng.sample([(u,v) for u in range(5) for v in range(5)],r)
        fs=[s.cancel(x**u*((x+8)/(x-1))**v) for u,v in es]
        H=(x-1)**max(v for _,v in es)
        polys=[s.Poly(s.cancel(H*f),x) for f in fs]
        k=max(p.degree() for p in polys)
        M=s.Matrix([[p.nth(i) for p in polys] for i in range(k+1)])
        if M.rank()!=r:
            continue
        Ps=[]
        for u,v in es:
            Q=u*(x+8)*(x-1)-9*v*x
            p=s.Integer(1)
            row=[]
            for n in range(r):
                row.append(s.expand(p))
                p=s.expand(D*s.diff(p,x)+(Q-n*s.diff(D))*p)
            Ps.append(row)
        A=s.Matrix([[Ps[j][i] for j in range(r)] for i in range(r)])
        a=s.Poly(s.expand(A.det(method="domain-ge")),x)
        assert not a.is_zero
        assert a.degree() <= r*(r-1), (es,a.degree(),r)
        W=s.det(s.Matrix([[s.diff(f,x,i) for f in fs] for i in range(r)]))
        assert s.cancel(W * D**(r*(r-1)//2) - s.prod(fs)*a.as_expr())==0
        trials+=1
        checks+=1
print(f"OK: {checks} independent rational-monomial cofactors nonzero with exact degree and Wronskian identity")
