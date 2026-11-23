#!/usr/bin/env python3
"""Depth thresholds for the fixed-m growth statement.
Lambda(w): largest l with 2^(l-1) < 7776 * binom(2l, <=w)^2 (paper, Section 6);
k_D = D + Lambda(D-1). The m-summand transfer excludes 1 from V_D once
k >= k_{max(D, C_m)}, C_m = 4(m-2) (or k >= 14 max(D,C_m) + 10).
It exceeds the constant carry floor m(m-1)+1 (Theorem 'Carry bound for every
number of summands') only at D >= m(m-1)+1."""
from math import comb
def binle(n, w): return sum(comb(n, i) for i in range(w + 1))
def Lam(w):
    l, best = 1, None
    for l in range(1, 4000):
        if 2 ** (l - 1) < 7776 * binle(2 * l, w) ** 2: best = l
    return best
def kD(D): return D + Lam(D - 1)
assert [kD(D) for D in range(2, 7)] == [27, 39, 52, 64, 77]
print("m  C_m  floor  first k (k_D form) excluding 1 at D=C_m | first k where bound > floor: k_D form, linear form")
for m in range(3, 9):
    C = 4 * (m - 2); F = m * (m - 1) + 1
    print(m, C, F, kD(max(C, 2)), "|", kD(F), 14 * F + 10)
