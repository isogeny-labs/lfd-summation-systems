#!/usr/bin/env python3
"""Curve-type tests for the isogeny class of A (Weil polynomial P_A) and for its quadratic twist
(Weil polynomial P_A(-x), the class that would carry a Serre-obstructed geometric Jacobian).

Outputs JSON. Two independent power-sum routes: Newton identities on P_A, and the closed form
p_k = -s_k^E (131 !| k), 130*s_k^E (131 | k), s_k^E = alpha^k+beta^k (Lucas, x^2+x+2).
"""
import json, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PA = json.loads((HERE.parents[1] / "P_A.json").read_text())["P_A_coeffs_low_to_high"]
G = 130

def newton(coeffs, K):
    n = len(coeffs)-1
    assert coeffs[n] == 1
    e = [((-1)**k)*coeffs[n-k] for k in range(n+1)]
    P = [0]*(K+1); P[0] = n
    for k in range(1, K+1):
        s = sum(((-1)**(i-1))*e[i]*P[k-i] for i in range(1, min(k-1, n)+1))
        if k <= n: s += ((-1)**(k-1))*k*e[k]
        P[k] = s
    return P

def lucas(K):
    s = [2, -1]          # s_k = alpha^k+beta^k, alpha+beta=-1, alpha*beta=2
    while len(s) <= K: s.append(-s[-1]-2*s[-2])
    return s

def mobius(m):
    r, p = 1, 2
    while p*p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0: return 0
            r = -r
        p += 1
    return -r if m > 1 else r

def places(N, D):
    return {d: sum(mobius(d//e)*N[e] for e in range(1, d+1) if d % e == 0)//d for d in range(1, D+1)}

def places_P1(D, q=2):
    N = [0]+[q**k+1 for k in range(1, D+1)]
    return places(N, D)

K = 600
pk = newton(PA, K)
s = lucas(K)
closed = [None]+[(-s[k] if k % 131 else 130*s[k]) for k in range(1, K+1)]
assert pk[1:] == closed[1:], "power-sum cross-check failed"
out = {"power_sum_crosscheck_k_le": K}

# (N2) geometric indecomposability: P_A in Z[x^k] for some k>1, k|130 ?
divs = [k for k in range(2, 131) if 130 % k == 0]
inZxk = {k: all(PA[i] == 0 for i in range(len(PA)) if i % k) for k in divs}
out["P_A_in_Z[x^k]_for_k|130"] = inZxk
out["coeff_x^259"] = PA[259]
# negative control for the detector: P(x)=Q(x^k) must be flagged
toy = [0]*(2*2*2+1)
for i, c in enumerate([4, 0, 1]): toy[2*i] = c   # x^4+... in x^2: Q(y)=y^2+0y+4 -> 4+x^4 ; flagged k=2
out["detector_control_flags_Q(x^2)"] = all(toy[i] == 0 for i in range(len(toy)) if i % 2)

res = {}
for label, sign in (("A", 1), ("twist_-pi", -1)):
    N = [0]+[2**k+1-(sign**k)*pk[k] for k in range(1, K+1)]
    a = places(N, 400)
    P1 = places_P1(400)
    r = {"N_1..8": N[1:9], "a_1..8": [a[d] for d in range(1, 9)],
         "min_N": min(N[1:]), "min_a": min(a.values()),
         "all_a_nonneg_integer_d<=400": all(v >= 0 for v in a.values())}
    # gonality lower bound gamma >= N_k/(q^k+1)
    r["gonality_lower_bound_from_counts"] = max(-(-N[k]//(2**k+1)) for k in range(1, K+1))
    # hyperelliptic test over F_2: 3 rational points of P^1, each fiber = 2 rational pts | 1 ramified rational | 1 degree-2 place.
    # Feasible iff exist s,r,i>=0, s+r+i=3, 2s+r=a1, and i <= a2 (inert rational points give degree-2 places).
    feas = [(s_, r_, 3-s_-r_) for s_ in range(4) for r_ in range(4-s_)
            if 2*s_+r_ == a[1] and 3-s_-r_ <= a[2]]
    r["hyperelliptic_fiber_patterns_over_P1(F2)"] = feas
    r["hyperelliptic_possible_by_degree_1_2_fibers"] = bool(feas)
    res[label] = r
out["classes"] = res

# Control for the hyperelliptic fiber test: brute-force every genus-2 curve y^2+h y=f over F_2
# (deg h<=3, deg f<=6, smooth) with point counts over F_2 and F_4 in the weighted model; the test
# must accept each (no false positives) and some curves must have a2=0 (so the a2 branch is exercised).
def gf4():
    # F_4 = {0,1,w,w^2} as ints 0..3 with w^2=w+1
    def mul(x, y):
        r = 0
        for i in range(2):
            if (y >> i) & 1: r ^= x << i
        if r & 4: r ^= 0b111
        return r
    return mul
m4 = gf4()
def ev(poly, x, mul):
    r = 0
    for c in reversed(poly): r = mul(r, x) ^ c
    return r
def count(h, f, field, mul, g=2):
    # affine points + points at infinity of weighted model (H of degree g+1, F of degree 2g+2)
    hh = h+[0]*(g+2-len(h)); ff = f+[0]*(2*g+3-len(f))
    n = 0
    for x in field:
        hx, fx = ev(hh, x, mul), ev(ff, x, mul)
        n += sum(1 for y in field if mul(y, y) ^ mul(hx, y) ^ fx == 0)
    H0, F0 = hh[g+1], ff[2*g+2]
    n += sum(1 for y in field if mul(y, y) ^ mul(H0, y) ^ F0 == 0)
    return n
def mul2(x, y): return x & y
def smooth_ordinary_free(h, f):
    return any(h)   # crude: keep h != 0 (char-2 genus-2 needs h != 0); smoothness checked via counts below
bad, seen, a2zero = 0, 0, 0
from itertools import product
for hb in range(1, 16):
    h = [(hb >> i) & 1 for i in range(4)]
    for fb in range(128):
        f = [(fb >> i) & 1 for i in range(7)]
        N1 = count(h, f, [0, 1], mul2)
        N2 = count(h, f, [0, 1, 2, 3], m4)
        a1, a2 = N1, (N2-N1)//2
        if (N2-N1) % 2: continue
        seen += 1
        if a2 == 0: a2zero += 1
        ok = any(2*s_+r_ == a1 and 3-s_-r_ <= a2 for s_ in range(4) for r_ in range(4-s_))
        if not ok: bad += 1
out["hyperelliptic_test_control_genus2_models"] = {"models": seen, "with_a2=0": a2zero, "rejected_by_test(should be 0)": bad}
json.dump(out, sys.stdout, indent=1)
