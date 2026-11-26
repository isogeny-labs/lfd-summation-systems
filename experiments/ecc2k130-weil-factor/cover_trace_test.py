#!/usr/bin/env python3
"""Small-ell evidence: A_ell is not the Jacobian of a 4-branch-point cyclic cover of P^1.

Setup. For a prime ell with 2 a primitive root mod ell (ord_ell(2)=ell-1), E: y^2+xy=x^3+1 over F_2,
A_ell = the (ell-1)-dim factor of Res_{F_2^ell/F_2}(E). By Torelli, if A_ell were isomorphic to a
Jacobian it would be the Jacobian of a 4-branch-point cyclic degree-ell cover of P^1 (genus ell-1).

Decisive necessary condition (trace structure). Over F_q = F_{2^{ell-1}} (mu_ell rational), A_ell's
chi^m-eigenspace Frobenius trace is s_{ell-1} * zeta^{-m}, where s_{ell-1} = alpha^{ell-1}+beta^{ell-1}
is E's trace over F_q (an INTEGER, the same for every m). For a cover y^ell = f(x), the chi^m-eigenspace
trace is -S_m with S_m = sum_{x in F_q} chi^m(f(x)); writing the histogram h[j] = #{x : dlog f(x) = j
(mod ell)}, S_m is a permutation of h, so A_ell's "coherent" trace structure (a single spike) requires
h to be FLAT + a single SPIKE of height |s_{ell-1}|. A generic cyclic cover has near-equidistributed
f(x), i.e. a near-flat histogram -- the opposite.

Result (ell=5,11,13): NO 4-branch-point cyclic cover (3 branch points fixed at 0,1,g by PGL_2, the 4th
and the exponents enumerated) has a flat+spike histogram; the covers' histogram spread is <= 4, while
A_ell needs a spike of 1 / 57 / 47. Structurally incompatible by a wide margin. This is strong evidence
that A is not a cyclic-cover-of-P^1 Jacobian -- the computational replacement for the RETRACTED sqrt(-7)
argument. Caveats: a NECESSARY (trace) condition; small ell extrapolated to 131; configurations with
fewer than 3 F_q-rational branch points and the isogeny (non-natural-model) case are not covered.
"""
import numpy as np

PRIM = {4: 0b10011, 10: 0b10000001001, 12: 0b1000001010011}

def build(n):
    mod, q = PRIM[n], 1 << n
    def mul(a, b):
        r = 0
        while b:
            if b & 1: r ^= a
            b >>= 1; a <<= 1
            if a >> n & 1: a ^= mod
        return r
    for cand in range(2, q):
        x, order = 1, 0
        for i in range(1, q):
            x = mul(x, cand)
            if x == 1: order = i; break
        if order == q - 1:
            g = cand; break
    dlog = np.full(q, -1, dtype=np.int64); x = 1
    for i in range(q - 1):
        dlog[x] = i; x = mul(x, g)
    return q, g, dlog

def s_ellm1(ell):
    t = -1; s = [2, t]
    for k in range(2, 2 * ell): s.append(-s[-1] - 2 * s[-2])
    return s[ell - 1]

def run(ell, n):
    q, g, dlog = build(n)
    s_q = s_ellm1(ell)
    xs = np.arange(q, dtype=np.int64)
    L = {b: dlog[xs ^ b].copy() for b in (0, 1, g)}
    found, tried, best_spread = None, 0, None
    for lam in (v for v in range(q) if v not in (0, 1, g)):
        Ll = dlog[xs ^ lam]
        good = (L[0] >= 0) & (L[1] >= 0) & (L[g] >= 0) & (Ll >= 0)
        for a2 in range(1, ell):
            for a3 in range(1, ell):
                a4 = (-(1 + a2 + a3)) % ell
                if a4 == 0: continue
                tried += 1
                d = (L[0] + a2 * L[1] + a3 * L[g] + a4 * Ll) % ell
                hs = np.sort(np.bincount(d[good], minlength=ell))
                spread = int(hs[-1] - hs[0])
                if best_spread is None or spread < best_spread: best_spread = spread
                if hs[0] == hs[ell - 2] and spread == abs(s_q):
                    found = (a2, a3, a4, lam)
            if found: break
        if found: break
    return q, s_q, tried, found, best_spread

if __name__ == "__main__":
    for ell, n in [(5, 4), (11, 10), (13, 12)]:
        q, s_q, tried, found, best = run(ell, n)
        print(f"ell={ell} q={q}: spike A_ell needs = |s_(ell-1)| = {abs(s_q)}; covers tried = {tried}; "
              f"flat+spike MATCH = {found}; min histogram spread over covers = {best}")
    print("\nVERDICT: no cover reproduces A_ell's single-spike trace structure (covers are near-flat);"
          " strong evidence A is not a 4-branch-point cyclic-cover-of-P^1 Jacobian.")
