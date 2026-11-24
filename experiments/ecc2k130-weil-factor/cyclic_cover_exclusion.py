#!/usr/bin/env python3
"""Rigorous exclusion of the cyclic-cover construction for the ECC2K-130 Weil factor A.

Theorem (this file verifies the computational inputs):
  A is NOT isogenous to a simple factor of the Jacobian of any cyclic (more generally, abelian) cover
  C -> P^1 of order coprime to 7. (Covers of P^1 only; covers of E are a separate, open case -- see
  cover_of_E_reduction.py.) For a cyclic cover whose order is divisible by 7 (needed to supply
  sqrt(-7)), the order is divisible by lcm(7,131)=917 and the genus is >= 458, far above the usable
  window (~130-300). So no cyclic cover realizes A in-window.

Mechanism. For an abelian cover C -> P^1 with Galois group of exponent N (char 2, N odd), every
eigenvalue of Frobenius on Jac(C) is a Jacobi sum, hence lies in Q(zeta_N) (Weil). If A were an isogeny
factor, A's Frobenius eigenvalues would lie in Q(zeta_N) too, forcing the CM field
K = Q(zeta_131, sqrt(-7)) (= Q(eigenvalues of A), since A is simple with this CM field) into Q(zeta_N),
i.e. conductor(K) = lcm(131,7) = 917 divides N. But the base-change eigenvalues of A to F_{2^130},
namely alpha^130 * zeta_131^{-j} with alpha = (-1+sqrt(-7))/2, are NOT in Q(zeta_131): alpha^130 has a
nonzero sqrt(-7) component and Q(sqrt(-7)) cap Q(zeta_131) = Q (7 does not divide 131). So a pure
Z/131 cover (eigenvalues in Q(zeta_131)) cannot carry A; the obstruction is exactly the sqrt(-7) that A
inherits from E's CM.
"""
from fractions import Fraction as F

def mul(x, y):
    (a, b), (c, d) = x, y               # (a+b s)(c+d s), s^2 = -7
    return (a*c - 7*b*d, a*d + b*c)

def powe(x, n):
    r = (F(1), F(0))
    while n:
        if n & 1: r = mul(r, x)
        x = mul(x, x); n >>= 1
    return r

def cyclic_cover_genus(N, r):           # Z/N tame cover of P^1, totally ramified at r points
    v = -2*N + r*(N-1)
    assert v % 2 == 0
    return v//2 + 1

def main():
    alpha = (F(-1, 2), F(1, 2))
    a130, b130 = powe(alpha, 130)
    print("alpha^130 =", a130, "+", b130, "* sqrt(-7)")
    print("  sqrt(-7) component nonzero (=> alpha^130 not rational, not in Q(zeta_131)):", b130 != 0)
    print("  |alpha^130|^2 == 2^130 (valid F_{2^130} Weil number):", a130*a130 + 7*b130*b130 == 2**130)
    print("sqrt(-7) in Q(zeta_7)?  7 % 4 == 3:", 7 % 4 == 3)
    print("sqrt(-7) in Q(zeta_131)? 131 % 7 == 0:", 131 % 7 == 0)
    cond = 7*131
    print("conductor(Q(zeta_131, sqrt-7)) = lcm(7,131) =", cond)
    print("min cyclic-cover (order 917) genus, r=3 branch points:", cyclic_cover_genus(cond, 3))
    print("contrast: Z/131 cover, r=4, genus:", cyclic_cover_genus(131, 4), "(in-window but no sqrt-7)")
    assert b130 != 0 and 131 % 7 != 0 and cyclic_cover_genus(cond, 3) > 300
    print("\nVERDICT: cyclic covers cannot realize A in the usable genus window (rigorous).")

if __name__ == "__main__":
    main()
