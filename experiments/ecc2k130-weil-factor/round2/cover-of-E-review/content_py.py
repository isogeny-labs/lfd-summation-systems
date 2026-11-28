#!/usr/bin/env python3
"""Independent cross-check (no polynomial division): P_A(alpha) = P_W'(alpha)/(2 alpha + 1),
since P_W = (x^2+x+2) P_A and alpha is a simple root of x^2+x+2. Exact arithmetic in Z[alpha],
alpha^2 = -alpha - 2. Then content = gcd of coordinates; compared with l*|U_l| (Lucas)."""
from math import gcd
def mul(u, v):  # (a + b al)(c + d al), al^2 = -al - 2
    a, b = u; c, d = v
    return (a*c - 2*b*d, a*d + b*c - b*d)
def pw(u, n):
    r = (1, 0)
    while n:
        if n & 1: r = mul(r, u)
        u = mul(u, u); n >>= 1
    return r
def s(l):  # alpha^l + beta^l, beta = -1 - alpha
    a, b = pw((0, 1), l); return 2*a - b
def U(n):
    a, b = 0, 1
    for _ in range(n-1): a, b = b, -b - 2*a
    return b if n else 0
def divide_by_sqrtm7(v):  # 2 al + 1 = sqrt(-7); (x)/(2al+1) = x*(2al+1)/(-7)
    a, b = mul(v, (1, 2)); assert a % 7 == 0 and b % 7 == 0
    return (-a//7, -b//7)
for l in [3, 5, 11, 13, 17, 19, 23, 29, 31, 37, 131]:
    al = (0, 1)
    t1 = pw(al, 2*l - 1); t1 = (2*l*t1[0], 2*l*t1[1])
    t2 = pw(al, l - 1);   t2 = (l*s(l)*t2[0], l*s(l)*t2[1])
    dPW = (t1[0] - t2[0], t1[1] - t2[1])
    PAal = divide_by_sqrtm7(dPW)
    c = gcd(*PAal)
    print(f"l={l:3d}  content(P_A(alpha)) = {c}   l*|U_l| = {l*abs(U(l))}   match={c == l*abs(U(l))}")
print("U_131 =", U(131), "; U_131/263 =", U(131)//263 if U(131) % 263 == 0 else None)
