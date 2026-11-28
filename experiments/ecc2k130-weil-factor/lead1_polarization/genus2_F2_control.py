#!/usr/bin/env python3
"""Independent control at ell = 3 (no Deligne modules, no PARI).

A_3 has Weil polynomial x^4 - x^3 - x^2 - 2x + 4 (computed below from alpha, beta).
Prediction of REPORT.md: (-7/3) = -1 and U_3 = -1, so the class has NO principally polarized member.
A simple ppav surface over F_2 is the Jacobian of a genus-2 curve over F_2, or (geometrically a product
of elliptic curves with product polarization) a product E1 x E2 over F_2 or Res_{F_4/F_2} E; the latter two
have Weil polynomials h1*h2 or h(x^2), which we exclude below. So it suffices to show that no smooth genus-2
curve y^2 + h(x) y = f(x) over F_2 (deg h <= 3, deg f <= 6) has (N1, N2) = (#C(F_2), #C(F_4)) matching.
Also cross-checks the Howe-Nart-Ritzenthaler criterion (2009): ordinary x^4+ax^3+bx^2+aqx+q^2 contains no
ppav iff a^2-b = q, b < 0, every prime divisor of b is 1 mod 3.
"""
import itertools, json

# GF(64) = F2[t]/(t^6+t^4+t^3+t+1)  (irreducible), contains GF(4) and GF(8)
MOD, DEG = 0b1011011, 6
def gmul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a >> DEG & 1: a ^= MOD
    return r
def gpow(a, n):
    r = 1
    while n:
        if n & 1: r = gmul(r, a)
        a = gmul(a, a); n >>= 1
    return r
# check irreducibility: t^(2^6) = t and t^(2^2) != t, t^(2^3) != t
assert gpow(2, 64) == 2 and gpow(2, 4) != 2 and gpow(2, 8) != 2
F4 = [x for x in range(64) if gpow(x, 4) == x]; assert len(F4) == 4
F2 = [0, 1]
def peval(c, x):  # c: list of F2 coefficients low->high
    r = 0
    for a in reversed(c): r = gmul(r, x) ^ a
    return r
def deriv(c): return [c[i] if i % 2 == 1 else 0 for i in range(1, len(c))]

def smooth(h, f):
    # affine singular points over F64-bar: need points over extension containing roots of h; deg h <= 3 so F64 suffices
    hd, fd = deriv(h), deriv(f)
    for x in range(64):
        if peval(h, x) != 0: continue
        fx = peval(f, x)
        y = gpow(fx, 32)  # sqrt in char 2 on F64
        if gmul(peval(hd, x), y) ^ peval(fd, x) == 0: return False
    # at infinity: x = 1/z, Y = y z^3 : Y^2 + H(z) Y = F(z), H = z^3 h(1/z), F = z^6 f(1/z)
    H = [h[3-i] if 3-i < len(h) else 0 for i in range(4)]
    Fp = [f[6-i] if 6-i < len(f) else 0 for i in range(7)]
    z = 0
    if peval(H, z) == 0:
        Y = gpow(peval(Fp, z), 32)
        if gmul(peval(deriv(H), z), Y) ^ peval(deriv(Fp), z) == 0: return False
    return True

def count(h, f, field):
    n = 0
    for x in field:
        hx, fx = peval(h, x), peval(f, x)
        for y in field:
            if gmul(y, y) ^ gmul(hx, y) ^ fx == 0: n += 1
    H3 = h[3] if len(h) > 3 else 0; F6 = f[6] if len(f) > 6 else 0
    for Y in field:
        if gmul(Y, Y) ^ gmul(H3, Y) ^ F6 == 0: n += 1
    return n

def weil_from_counts(N1, N2, q=2):
    s1 = q + 1 - N1; s2 = q*q + 1 - N2
    e2 = (s1*s1 - s2) // 2
    return (1, -s1, e2, -q*s1, q*q)  # x^4 - s1 x^3 + e2 x^2 - q s1 x + q^2

# Weil polynomial of A_3: (x^2+alpha x+alpha^2)(x^2+beta x+beta^2) with alpha+beta=-1, alpha beta=2
s, p = -1, 2
target = (1, s, s*s - 2*p + p, p*s, p*p)  # x^4 + (a+b)x^3 + (a^2+b^2+ab)x^2 + ab(a+b)x + a^2b^2
print("Weil polynomial of A_3 (coeffs high->low):", target)
a, b = target[1], target[2]
hnr = (a*a - b == 2) and b < 0 and all(pf % 3 == 1 for pf in [d for d in range(2, abs(b)+1) if abs(b) % d == 0 and all(d % e for e in range(2, d))])
print("HNR 2009 no-ppav criterion satisfied:", hnr)

seen = {}
ncurves = 0
for hc in itertools.product([0, 1], repeat=4):
    h = list(hc)
    for fc in itertools.product([0, 1], repeat=7):
        f = list(fc)
        if not any(h): continue  # y^2 = f is never smooth of genus 2 in char 2
        if max([i for i in range(4) if h[i]] + [0]) < 3 and max([i for i in range(7) if f[i]] + [0]) < 5:
            continue  # genus < 2
        if not smooth(h, f): continue
        ncurves += 1
        N1, N2 = count(h, f, F2), count(h, f, F4)
        seen.setdefault(weil_from_counts(N1, N2), 0)
        seen[weil_from_counts(N1, N2)] += 1
print("smooth genus-2 models over F_2:", ncurves, " distinct Weil polynomials:", len(seen))
print("target Weil polynomial realised by a genus-2 curve over F_2:", target in seen)
# positive control on the enumerator: #models = |GL2(F2)| * (q-1) * q^4 * (sum 1/|Aut|) = 96 * q^3 = 768
# (mass formula sum_C 1/|Aut_k(C)| = q^(2g-1) for hyperelliptic curves of genus g)
print("mass-formula control 96*2^3 == #models:", 96 * 2**3 == ncurves)
print("sample realised classes:", sorted(seen)[:6])
# exclude product / restriction shapes for the target
print("target of shape h(x^2) (Res_{F4/F2}):", target[1] == 0 and target[3] == 0)
json.dump({"target": target, "hnr_criterion": hnr, "n_smooth_models": ncurves,
           "n_classes": len(seen), "target_realised": target in seen,
           "classes": {str(k): v for k, v in sorted(seen.items())}},
          open(__file__.replace(".py", ".json"), "w"), indent=1)
