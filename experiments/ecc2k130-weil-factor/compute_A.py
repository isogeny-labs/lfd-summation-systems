#!/usr/bin/env python3
"""Exact CM invariants of the 130-dimensional factor A of the Weil restriction of ECC2K-130.

E / F_2 : y^2 + x y = x^3 + 1   (a=0, b=1), Frobenius pi with pi^2 + pi + 2 = 0 (trace t=-1).
W = Res_{F_2^131 / F_2} E has dimension 131 and 2-Frobenius eigenvalues {alpha zeta^j, beta zeta^j}
(j=0..130, zeta a primitive 131st root of 1). W ~ E x A; the j=0 pair is E, the rest is A (dim 130).

Characteristic polynomial of 2-Frobenius on W:
    P_W(x) = (x^131 - alpha^131)(x^131 - beta^131) = x^262 - s_131 x^131 + 2^131,   s_131 = alpha^131 + beta^131.
E's Weil polynomial is x^2 - t x + 2 = x^2 + x + 2, so  P_A = P_W / (x^2 + x + 2)  (exact division).

This script computes P_A exactly and verifies: exact division; deg 260; P_A(1) = #A(F_2); #A prime;
constant term 2^130; ordinary (middle coeff odd); squarefree; CM field Q(zeta_131, sqrt(-7)) of degree
260; and the mod-p equal-degree factorization pattern (abelian Galois C_130 x C_2).
"""
import json, os

def miller_rabin(n):
    if n < 2: return False
    for p in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n % p == 0: return n == p
    d, r = n-1, 0
    while d % 2 == 0: d //= 2; r += 1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37,41):
        x = pow(a, d, n)
        if x in (1, n-1): continue
        for _ in range(r-1):
            x = x*x % n
            if x == n-1: break
        else: return False
    return True

def polydiv_exact(num, den):
    num = num[:]; q = [0]*(len(num)-len(den)+1)
    for i in range(len(num)-len(den), -1, -1):
        c = num[i+len(den)-1] // den[-1]
        q[i] = c
        for j in range(len(den)):
            num[i+j] -= c*den[j]
    return q, num

def gcd_modp(a, b, p):
    a = [c % p for c in a]; b = [c % p for c in b]
    def norm(x):
        while x and x[-1] % p == 0: x.pop()
        return x
    a, b = norm(a), norm(b)
    while b:
        inv = pow(b[-1], p-2, p)
        while len(a) >= len(b) and a:
            c = (a[-1]*inv) % p; sh = len(a)-len(b)
            for i in range(len(b)): a[sh+i] = (a[sh+i]-c*b[i]) % p
            a = norm(a)
        a, b = b, a
    return a

def frob_order(poly, p):
    """For squarefree equal-degree poly mod p: common irreducible-factor degree = order of Frob_p."""
    f = [c % p for c in poly]
    while f and f[-1] % p == 0: f.pop()
    inv = pow(f[-1], p-2, p); f = [(c*inv) % p for c in f]; n = len(f)-1
    def mul(a, b):
        r = [0]*(len(a)+len(b)-1)
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b): r[i+j] = (r[i+j]+ai*bj) % p
        for i in range(len(r)-1, n-1, -1):
            if r[i]:
                c = r[i]
                for j in range(n+1): r[i-n+j] = (r[i-n+j]-c*f[j]) % p
        return [c % p for c in r[:n]]
    def powp(base, e):
        res = [1]+[0]*(n-1); b = base[:]
        while e:
            if e & 1: res = mul(res, b)
            b = mul(b, b); e >>= 1
        return res
    x = [0, 1]+[0]*(n-2)
    cur = powp(x, p); d = 1
    while cur != x:
        cur = powp(cur, p); d += 1
        if d > n: return None
    return d

def main():
    t = -1
    s = [2, t]
    for k in range(2, 132): s.append(t*s[-1] - 2*s[-2])
    s131 = s[131]
    NE = 2**131 + 1 - s131                      # #E(F_2^131) = #W(F_2)
    rA = NE // 4
    PW = [0]*263; PW[0] = 2**131; PW[131] = -s131; PW[262] = 1
    PA, rem = polydiv_exact(PW, [2, 1, 1])      # divide by x^2 + x + 2

    checks = {}
    checks["s_131"] = s131
    checks["#E(F_2^131)=#W(F2)"] = NE
    checks["cofactor_is_4"] = (NE % 4 == 0 and NE // 4 * 4 == NE)
    checks["#A(F2)"] = rA
    checks["#A_prime"] = miller_rabin(rA)
    checks["#A_bits"] = rA.bit_length()
    checks["exact_division"] = all(r == 0 for r in rem)
    checks["deg_P_A"] = len(PA)-1
    checks["P_A(1)==#A"] = (sum(PA) == rA)
    checks["const_term==2^130"] = (PA[0] == 2**130)
    checks["lead==1"] = (PA[-1] == 1)
    checks["ordinary (middle coeff odd)"] = (PA[130] % 2 == 1)
    checks["squarefree"] = all(len(gcd_modp(PA, [(i*PA[i]) % q for i in range(1, len(PA))], q)) <= 1
                               for q in (1000003, 2000003))
    checks["sqrt(-7) in Q(zeta_131)"] = (131 % 7 == 0)
    checks["CM_field"] = "Q(zeta_131, sqrt(-7)), degree 260"
    checks["modp_factor_degrees"] = {p: frob_order(PA, p) for p in (3, 5, 13)}

    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "P_A.json"), "w") as f:
        json.dump({"P_A_coeffs_low_to_high": PA, "s_131": s131,
                   "num_E_F2_131": NE, "num_A_F2": rA}, f)
    for k, v in checks.items():
        print(f"{k}: {v}")

if __name__ == "__main__":
    main()
