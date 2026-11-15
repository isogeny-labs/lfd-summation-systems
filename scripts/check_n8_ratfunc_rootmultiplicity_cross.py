#!/usr/bin/env python3
"""Independent exact Q(X) regression for N8 finite-pole fraction identities."""
from itertools import product
from sympy import Poly, Rational, cancel, symbols

x = symbols("x")
X1 = Poly(x - 1, x, domain="QQ")
checks = 0

def order_at_one(expr):
    p = Poly(expr, x, domain="QQ")
    assert not p.is_zero
    n = 0
    while p.eval(1) == 0:
        p = p.exquo(X1)
        n += 1
    return n

for u, v, k in product(range(6), repeat=3):
    A = (x - 1)**u * (x + 2)**(1 + u % 2)
    B = (x - 1)**v * (x + 3)**(1 + v % 2)
    R = (x - 1)**k * (x + 5)**(1 + k % 2)
    gamma = Rational(k + 1, k + 2)
    C = A * R
    D = B * R
    assert cancel(A/B - C/D) == 0
    assert order_at_one(A) + order_at_one(D) == order_at_one(C) + order_at_one(B)
    # Here A/B = gamma * C/(gamma*D), so the integer order is unchanged
    assert cancel(A/B - gamma*C/(gamma*D)) == 0
    assert order_at_one(A)-order_at_one(B) == order_at_one(gamma*C)-order_at_one(gamma*D)
    checks += 1

print(f"PASS: {checks} exact rational fraction identities, cross-products and X=1 root orders")
