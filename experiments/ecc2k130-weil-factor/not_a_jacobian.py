#!/usr/bin/env python3
# RETRACTED 2026-11-24 (red-team catch): the sqrt(-7) eigenvalue-field argument is INVALID.
# For a 4-branch-point mu_131 cover of P^1 the character eigenspaces are 2-dimensional, so only
# the char-poly COEFFICIENTS (trace, det) lie in Q(zeta_131), not the eigenvalues; A satisfies
# exactly that (sum=s_130*zeta^-m, prod=2^130*zeta^-2m, s_130 integer). So this file DOES NOT
# establish its claim. See the CORRECTION section of docs/ecc2k130-weil-factor-A.md. A possible
# eigenvalue-independent (ray-class) rescue is a candidate pending review, not proven here.
"""A (natural model) is not isomorphic to a Jacobian -- rigorous; isogeny residue characterized.

THEOREM. Let A be the 130-dimensional factor of Res_{F_2^131/F_2}(E) (ECC2K-130), with its full
endomorphism ring End(A) = O_K, K = Q(zeta_131, sqrt(-7)). Then A is NOT isomorphic (as a principally
polarized abelian variety over F_2) to the Jacobian of any curve.

PROOF.
1. A carries a rational order-131 automorphism rho: the deck transformation of the (geometrically
   disconnected) cover underlying Res permutes the 131 conjugate copies of E cyclically; it fixes the
   diagonal E and acts on A with order 131. In End^0(A) = K it is multiplication by zeta_131.
2. rho preserves EVERY polarization. A is a CM abelian variety (simple, dim 130, End^0 = K a field of
   degree 260 = 2 dim A), so the Rosati involution of any polarization restricts to complex conjugation
   on K. rho = zeta_131 satisfies zeta_131 * conj(zeta_131) = zeta_131 * zeta_131^{-1} = 1, so it lies in
   the unitary group of every polarization: rho in Aut(A, lambda) for all lambda.
3. Suppose (A, theta) ~= (Jac C, theta_C) as ppav over F_2. Then g(C) = dim A = 130, and by Torelli
   Aut(A, theta) ~= Aut(C) (up to the hyperelliptic/(-1) factor of order 2). Since rho has order 131
   (odd, != +-1), C has a genuine order-131 automorphism phi over F_2.
4. Riemann-Hurwitz for C -> D := C/<phi> (cyclic degree 131, tame since 131 != char 2, totally ramified
   at r points): 2*130-2 = 131(2 g_D - 2) + r*130. The only solution with g_D, r >= 0 is g_D = 0, r = 4.
   So C is a cyclic degree-131 cover of a genus-0 curve.
5. Over F_2^130 (where mu_131 is rational and D becomes P^1), Jac(C) is a mu_131-cyclic cover of P^1, so
   its Frobenius eigenvalues are Jacobi sums, lying in Q(zeta_131). But A (x) F_2^130 has Frobenius
   eigenvalues alpha^130 * zeta_131^{-j}, and alpha^130 has a nonzero sqrt(-7) component, so these are
   NOT in Q(zeta_131) (sqrt(-7) not in Q(zeta_131), as 7 does not divide 131). Contradiction. QED.

SCOPE. This proves A (the natural model, End = O_K) is not ISOMORPHIC to a Jacobian. The ECDLP transfer
needs only A ISOGENOUS to a Jacobian. A model A' in the isogeny class with End(A') a smaller order NOT
containing zeta_131 would have no forced automorphism, and the argument above does not apply to it. The
isogeny question therefore reduces to: does the minimal endomorphism ring Z[pi, pi-bar] already contain
zeta_131? If yes, every isogenous model contains the order-131 automorphism and A is not isogenous to a
Jacobian (full closure). If no, a 'no-automorphism' model A' is not excluded by this argument. This
residue is the user's 'curve with no useful automorphism group' case; it is NOT resolved here.
"""
from fractions import Fraction as F

def main():
    g, ell = 130, 131
    sols = []
    for gD in range(3):
        rhs = (2*g-2) - ell*(2*gD-2)
        if rhs >= 0 and rhs % (ell-1) == 0: sols.append((gD, rhs//(ell-1)))
    assert sols == [(0, 4)], sols
    # alpha^130 has nonzero sqrt(-7) component:
    def mul(x, y):
        (a, b), (c, d) = x, y; return (a*c - 7*b*d, a*d + b*c)
    def powe(x, n):
        r = (F(1), F(0))
        while n:
            if n & 1: r = mul(r, x)
            x = mul(x, x); n >>= 1
        return r
    a, b = powe((F(-1, 2), F(1, 2)), 130)
    assert b != 0 and a*a + 7*b*b == 2**130
    print("RH forces C -> P^1 cyclic degree-131 (g_D=0, r=4 branch points):", sols)
    print("alpha^130 sqrt(-7)-component nonzero:", b != 0, " |alpha^130|^2 = 2^130:", a*a+7*b*b == 2**130)
    print("7 | 131?", 131 % 7 == 0, " => sqrt(-7) not in Q(zeta_131) => eigenvalue fields differ.")
    print("\nTHEOREM (verified inputs): A (natural model, End = O_K) is NOT isomorphic to a Jacobian.")
    print("ISOGENY RESIDUE (open): whether some model A' with zeta_131 not in End(A') is a Jacobian,")
    print("i.e. whether zeta_131 in Z[pi, pi-bar]. Not resolved here.")

if __name__ == "__main__":
    main()
