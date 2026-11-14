#!/usr/bin/env python3
"""Independent exact rank-preserving rational denominator and Wronskian tests."""
import random
import sympy as s
x = s.symbols("x")
rng = random.Random(20261113)
cases=0
for r in range(1,5):
    for t in range(10):
        while True:
            U=s.Matrix([[rng.randint(-3,3) for _ in range(r)] for __ in range(r)])
            if U.det()!=0: break
        degrees=list(range(r))
        ps=[sum(U[i,j]*x**degrees[i] for i in range(r)) for j in range(r)]
        H=(x-1)**(t%3+1)*(x+2)**(1+t%2)
        fs=[s.cancel(p/H) for p in ps]
        Wp=s.det(s.Matrix([[s.diff(p,x,i) for p in ps] for i in range(r)]))
        Wf=s.det(s.Matrix([[s.diff(f,x,i) for f in fs] for i in range(r)]))
        assert s.simplify(Wp) != 0
        assert s.cancel(Wf*H**r-Wp)==0
        assert s.simplify(Wf)!=0
        cases+=1
print(f"OK: {cases} exact rational Wronskian independence/denominator tests")
