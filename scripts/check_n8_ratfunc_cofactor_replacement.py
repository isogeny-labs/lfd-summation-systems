#!/usr/bin/env python3
"""Exact independent rational Wronskian column-replacement regression.

For 1<=r<=4, construct K-linearly independent rational functions f_j
over Q(X), choose a constant linear combination g=sum c_j f_j, and
verify det(W with column i replaced by derivatives of g)=c_i det(W).
"""
import random
import sympy as s

x = s.symbols("x")
rng = random.Random(20261113)
count = 0
for r in range(1, 5):
    for trial in range(8):
        H = (x - 1) ** (1 + trial % 2)
        ps = [x ** j + (j + 1) * (x + 2) ** (j - 1) if j > 0 else s.Integer(1)
              for j in range(r)]
        f = [s.cancel(p / H) for p in ps]
        c = [rng.choice([-3, -2, -1, 1, 2, 3]) for _ in range(r)]
        g = s.cancel(sum(c[j] * f[j] for j in range(r)))
        A = s.Matrix([[s.diff(f[j], x, i) for j in range(r)] for i in range(r)])
        W = s.cancel(A.det())
        assert W != 0
        for j in range(r):
            B = A.copy()
            for i in range(r):
                B[i, j] = s.diff(g, x, i)
            assert s.cancel(B.det() - c[j] * W) == 0
            count += 1
print(f"OK: {count} exact rational Wronskian cofactor column replacements")
