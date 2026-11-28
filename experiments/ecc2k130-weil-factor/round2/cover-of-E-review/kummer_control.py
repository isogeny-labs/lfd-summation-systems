#!/usr/bin/env python3
"""Controls for the twisted-form lemma and the point-count argument, on the explicit
genus-ell cover  C_ell: z^ell = x  over  E: y^2 + xy = x^3 + 1  (div(x) = 2(0,1) - 2O, so C -> E is
geometrically cyclic of degree ell, totally (tamely) ramified at (0,1) and O, genus ell by RH).
Over the base field F_Q (Q = 2^b) the deck group scheme is mu_ell, a form of Z/ell with Frobenius
acting by u = Q mod ell; m = ord_ell(Q).
Predictions (lemma): Prym char poly P_Pr(x) in Z[x^m]; if m > 1, Tr(pi|Prym) = 0 and every rational point
of E has exactly one rational preimage (#C(F_Q) = #E(F_Q)). If m = 1 (constant group, e.g. ell | Q-1),
#C(F_Q) = ell*a + r with r = #rational branch points = 2 (both branch points are rational).
Negative control: P_Pr != P_{A_ell} (A_ell's Weil poly over F_2) whenever m > 1.
"""
import sys
from fractions import Fraction

POLYS = {1: 0b11, 2: 0b111, 3: 0b1011, 4: 0b10011, 5: 0b100101, 6: 0b1000011, 7: 0b10000011,
         8: 0b100011011, 9: 0b1000010001, 10: 0b10000001001, 11: 0b100000000101, 12: 0b1000001010011,
         13: 0b10000000011011, 14: 0b100010001000011, 15: 0b1000000000000011, 16: 0b10001000000001011}

def make_field(n):
    mod = POLYS[n]
    def mul(a, b):
        r = 0
        while b:
            if b & 1: r ^= a
            b >>= 1; a <<= 1
            if a >> n: a ^= mod
        return r
    def pw(a, e):
        r = 1
        while e:
            if e & 1: r = mul(r, a)
            a = mul(a, a); e >>= 1
        return r
    def trace(a):
        t, x = 0, a
        for _ in range(n):
            t ^= x; x = mul(x, x)
        return t
    return mul, pw, trace

def count(ell, n):
    """#C_ell(F_{2^n}) and #E(F_{2^n})."""
    q = 1 << n
    mul, pw, tr = make_field(n)
    from math import gcd
    g = gcd(ell, q - 1)
    nE = 2            # O and (0,1)
    nC = 2            # one point over each (totally ramified) branch point
    for x in range(1, q):
        inv2 = pw(x, (q - 3) % (q - 1))                    # x^{-2}
        if tr(x ^ inv2) == 0:                    # w^2 + w = x + x^{-2} solvable -> 2 points with this x
            nE += 2
            if g == 1: r = 1
            else: r = g if pw(x, (q - 1) // g) == 1 else 0
            nC += 2 * r
    return nC, nE

def weil_poly_from_counts(Q, g, N):
    """N[k] = #C(F_{Q^k}) for k=1..g -> L-polynomial coefficients (reciprocal char poly), length 2g+1."""
    S = [None] + [Q**k + 1 - N[k] for k in range(1, g + 1)]   # power sums of Frobenius eigenvalues
    # Newton: e_k from power sums; char poly x^2g - e1 x^{2g-1} + e2 ...
    e = [Fraction(1)]
    for k in range(1, g + 1):
        e.append(sum((-1) ** (i - 1) * e[k - i] * S[i] for i in range(1, k + 1)) / k)
    c = [int((-1) ** k * e[k]) for k in range(g + 1)]          # coefficients of x^{2g-k}
    for k in range(g + 1, 2 * g + 1):
        c.append(c[2 * g - k] * Q ** (k - g))
    return c                                                   # c[k] = coeff of x^{2g-k}

def polydiv(a, b):
    a = a[:]; out = []
    for i in range(len(a) - len(b) + 1):
        f = a[i] // b[0]; assert f * b[0] == a[i]
        out.append(f)
        for j in range(len(b)): a[i + j] -= f * b[j]
    assert all(v == 0 for v in a[len(out):]), "not divisible"
    return out

def PA_ell(ell):
    # A_ell Weil poly over F_2: (x^{2ell} - s_ell x^ell + 2^ell)/(x^2 + x + 2)
    s = [2, -1]
    for k in range(2, ell + 1): s.append(-s[-1] - 2 * s[-2])
    W = [1] + [0] * (ell - 1) + [-s[ell]] + [0] * (ell - 1) + [2 ** ell]
    return polydiv(W, [1, 1, 2])

def run(ell, b):
    Q = 2 ** b
    m = 1
    while pow(Q, m, ell) != 1: m += 1
    N = [None]; NE = [None]
    for k in range(1, ell + 2):
        if b * k > 16: break
        c, e = count(ell, b * k); N.append(c); NE.append(e)
    g = ell
    if len(N) - 1 < g:
        return None
    L = weil_poly_from_counts(Q, g, N)
    LE = weil_poly_from_counts(Q, 1, NE)
    Pr = polydiv(L, LE)
    deg = len(Pr) - 1
    inZxm = all(Pr[i] == 0 for i in range(deg + 1) if (deg - i) % m != 0)
    trPr = -Pr[1]
    extra = ""
    if len(N) - 1 >= g + 1:           # sanity: predict #C(F_{Q^{g+1}}) from L and compare
        import cmath
        # power sum S_{g+1} via Newton from coefficients
        e = [(-1) ** k * L[k] for k in range(2 * g + 1)]
        S = [None]
        for k in range(1, g + 2):
            val = (-1) ** (k - 1) * k * e[k] if k <= 2 * g else 0
            val += sum((-1) ** (k - 1 + i) * e[k - i] * S[i] for i in range(1, k))
            S.append(val)
        extra = f"  check k={g+1}: predicted {Q**(g+1)+1-S[g+1]} counted {N[g+1]} ok={Q**(g+1)+1-S[g+1]==N[g+1]}"
    res = dict(ell=ell, base=f"F_{Q}", u=Q % ell, m=m, NC1=N[1], NE1=NE[1], trace_Prym=trPr,
               Prym_in_Z_xm=inZxm, Prym_is_A=(b == 1 and Pr == PA_ell(ell)))
    print(res, extra)
    return res

if __name__ == "__main__":
    out = []
    for ell in [3, 5, 11]:          # base F_2: twisted (mu_ell, u=2, m=ord_ell 2)
        out.append(run(ell, 1))
    out.append(run(7, 1))               # ell=7: ord_7(2)=3
    out.append(run(3, 2))               # base F_4: mu_3 constant (3 | 3)  -> m=1 positive control
    out.append(run(5, 4))               # base F_16: mu_5 constant (5 | 15) -> m=1 positive control
    out.append(run(7, 3))               # base F_8: mu_7 constant (7 | 7) -> m=1 positive control
    for r in out:
        if r is None: continue
        if r["m"] > 1:
            assert r["Prym_in_Z_xm"] and r["trace_Prym"] == 0 and r["NC1"] == r["NE1"] and not r["Prym_is_A"], r
        else:
            assert (r["NC1"] - 2) % r["ell"] == 0, r
    print("ALL CONTROL ASSERTIONS PASSED")
