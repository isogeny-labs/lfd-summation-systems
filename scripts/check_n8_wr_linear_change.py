#!/usr/bin/env python3
"""Independent exact checks of the genuine rational Wronskian basis-change law.

This is a regression, not a substitute for the Lean proof. Uses SymPy's
rational differentiation and matrix determinants, with no floating point.
"""
from sympy import Matrix, Rational, Symbol, cancel

x = Symbol("x")
z = (x + 8) / (x - 1)
checks = 0
for n in (1, 2, 3):
    for seed in range(9):
        f = [x ** ((i + seed) % 4) * z ** ((2 * i + seed // 2) % 3)
             for i in range(n)]
        U = Matrix(n, n, lambda i, j: Rational(
            ((i + 1) * (j + 2) + seed * (i - j)) % 7 - 3))
        if seed % 3 == 0 and n > 1:
            for j in range(n):
                U[1, j] = U[0, j]
        g = [sum(U[i, j] * f[i] for i in range(n)) for j in range(n)]
        Wf = Matrix(n, n, lambda i, j: f[j].diff(x, i))
        Wg = Matrix(n, n, lambda i, j: g[j].diff(x, i))
        error = cancel(Wg.det(method="berkowitz") -
                       U.det(method="berkowitz") * Wf.det(method="berkowitz"))
        assert error == 0, (n, seed, error)
        checks += 1
print(f"PASS: {checks} exact genuine rational Wronskian basis changes (r=1..3)")
