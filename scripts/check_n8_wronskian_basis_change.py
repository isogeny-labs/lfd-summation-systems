#!/usr/bin/env python3
"""Independent exact Wronskian basis change for polynomial subspaces over QQ."""
import random
import sympy as s
x=s.symbols('x')
rng=random.Random(102819)
cases=0
for r in range(1,6):
    ds=[2*j+1 for j in range(r)]
    base=[x**d + (j+1)*x**(d-1) for j,d in enumerate(ds)]
    W=s.det(s.Matrix([[s.diff(p,x,i) for p in base] for i in range(r)]))
    assert s.expand(W) != 0
    for t in range(8):
        U=s.Matrix([[rng.randint(-4,4) for _ in range(r)] for __ in range(r)])
        if U.det()==0: continue
        new=[s.expand(sum(U[i,j]*base[i] for i in range(r))) for j in range(r)]
        V=s.det(s.Matrix([[s.diff(p,x,i) for p in new] for i in range(r)]))
        assert s.expand(V-U.det()*W)==0
        assert s.expand(V)!=0
        cases+=1
print(f"OK: {cases} exact polynomial Wronskian basis changes across dimensions 1..5")
