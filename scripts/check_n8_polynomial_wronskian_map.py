#!/usr/bin/env python3
"""Independent exact polynomial/rational Wronskian embedding regression."""
import sympy as sp
X = sp.symbols("X")
families = [
  [sp.Integer(1), X, X**2, X**3],
  [X**4 + X + 1, X**2 - 3, 1 + X**3, X**4 + 2 * X],
  [X + 8, X - 1, X**2 + X + 7, X**3 - X + 1],
]
scalars = [sp.Integer(1), X + 1, (X + 8)/(X - 1), 1/(X**2 + 1)]
count = 0
for p in families:
    for r in range(1,5):
        a = p[:r]
        w = sp.expand(sp.det(sp.Matrix([[sp.diff(f,X,i) for f in a] for i in range(r)])))
        # Exact equality after embedding into Q(X), and after common rational scaling.
        for h in scalars:
            wh = sp.det(sp.Matrix([[sp.diff(h*f,X,i) for f in a] for i in range(r)]))
            assert sp.cancel(wh - h**r * w) == 0, (a,r,h)
            count += 1
print(f"N8 polynomial-to-rational Wronskian: {count} exact determinant identities passed")
