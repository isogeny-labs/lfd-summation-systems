#!/usr/bin/env python3
"""Computable Jacobian necessary-condition tests for the ECC2K-130 Weil factor A.

If A were isogenous to the Jacobian of a (genus-130) curve C / F_2, then C's zeta numerator is the
reciprocal of A's Frobenius characteristic polynomial P_A, and the implied point counts
    #C(F_{2^k}) = 2^k + 1 - p_k,   p_k = sum of k-th powers of P_A's roots,
and the implied closed-point counts
    a_d = (1/d) sum_{e | d} mu(d/e) #C(F_{2^e})
must ALL be non-negative integers (a_d counts degree-d places). A single negative value would prove no
curve has this zeta function, so A would not be a Jacobian. This script computes p_k by Newton's
identities from P_A (cross-checked against the exact structural formula p_k = -s_k^E for 131 nmid k,
130 s_k^E for 131 | k, with s^E the trace sequence of E), then checks #C and a_d.
"""
import json, os

def power_sums(coeffs, upto):
    n = len(coeffs)-1; c = coeffs
    e = [((-1)**k)*c[n-k] for k in range(n+1)]       # elementary symmetric of the roots
    P = [0]*(upto+1); P[0] = n
    for k in range(1, upto+1):
        s = 0
        for i in range(1, min(k-1, n)+1):            # Newton: i = 1..k-1
            s += ((-1)**(i-1))*e[i]*P[k-i]
        if k <= n:
            s += ((-1)**(k-1))*k*e[k]
        P[k] = s
    return P

def mobius(m):
    if m == 1: return 1
    res, mm, p = 1, m, 2
    while p*p <= mm:
        if mm % p == 0:
            mm //= p
            if mm % p == 0: return 0
            res = -res
        p += 1
    if mm > 1: res = -res
    return res

def main(B=400):
    here = os.path.dirname(os.path.abspath(__file__))
    PA = json.load(open(os.path.join(here, "P_A.json")))["P_A_coeffs_low_to_high"]
    P = power_sums(PA, 4*B)
    # structural cross-check
    sE = [2, -1]
    for k in range(2, 4*B+1): sE.append(-sE[-1] - 2*sE[-2])
    pred = lambda k: 130*sE[k] if k % 131 == 0 else -sE[k]
    assert all(P[k] == pred(k) for k in range(1, 4*B)), "Newton vs structural power-sum mismatch"

    C = lambda k: 2**k + 1 - P[k]
    div = lambda d: [i for i in range(1, d+1) if d % i == 0]
    a = {}
    for d in range(1, B+1):
        val = sum(mobius(d//e)*C(e) for e in div(d))
        assert val % d == 0, f"non-integer place count at d={d}"
        a[d] = val // d
    neg_C = [k for k in range(1, 4*B) if C(k) < 0]
    neg_a = [(d, a[d]) for d in a if a[d] < 0]
    print(f"power sums cross-checked against structural formula: OK (k=1..{4*B-1})")
    print(f"#C(F_2)={C(1)}  #C(F_4)={C(2)}  #C(F_8)={C(3)}  #C(F_16)={C(4)}")
    print(f"min #C(F_2^k), k=1..{4*B-1}: {min(C(k) for k in range(1,4*B))}")
    print(f"negative #C(F_2^k): {neg_C if neg_C else 'NONE'}")
    print(f"closed-place counts a_d, d=1..12: {[(d,a[d]) for d in range(1,13)]}")
    print(f"min a_d, d=1..{B}: {min(a.values())}")
    print(f"negative a_d: {neg_a if neg_a else 'NONE'}")
    print("VERDICT:", "A passes all computable Jacobian necessary conditions (not ruled out)"
          if not neg_C and not neg_a else "A is NOT a Jacobian (a necessary condition fails)")

if __name__ == "__main__":
    main()
