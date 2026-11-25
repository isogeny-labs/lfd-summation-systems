#!/usr/bin/env python3
"""Index-calculus cost model for the HYPOTHETICAL genus-130 Jacobian carrying A.

Conditional on A ~ Jac(C) for a genus-130 curve C / F_2 (the open window question; see
docs/ecc2k130-weil-factor-A.md), C's place counts a_d and effective-divisor counts B_m are fixed by
P_A. This script computes the standard smooth-divisor index-calculus cost from those EXACT counts and
compares it to the matched Pollard rho reference (~2^61). It makes the stakes precise; it does not
assert A is a Jacobian.

Model (Enge-Gaudry style): factor base = places of degree <= t, N(t) = sum_{d<=t} a_d. A relation is a
random reduced divisor of degree ~ g that is t-smooth; P_smooth = (t-smooth deg-g divisors)/B_g.
Relation collection = N(t)/P_smooth trials; sparse linear algebra ~ N(t)*g. Cost reported in trials and,
honestly, also multiplied by a poly(g) per-trial smoothness-test factor.
"""
import json, math, os

def load_PA():
    here = os.path.dirname(os.path.abspath(__file__))
    return json.load(open(os.path.join(here, "P_A.json")))["P_A_coeffs_low_to_high"]

def main(tmax=26, g=130):
    PA = load_PA(); n = 260; c = PA
    e = [((-1)**k)*c[n-k] for k in range(n+1)]
    l = [((-1)**k)*e[k] for k in range(n+1)]          # L(T)=prod(1-lambda_i T) coeffs
    # effective divisor counts B_m
    B = [sum(l[k]*(2**(m-k+1)-1) for k in range(min(m, n)+1)) for m in range(g+1)]
    assert all(b >= 0 for b in B), "negative effective-divisor count"
    # power sums -> point counts -> place counts a_d
    P = [0]*(g+1); P[0] = n
    for k in range(1, g+1):
        s = sum(((-1)**(i-1))*e[i]*P[k-i] for i in range(1, min(k-1, n)+1))
        if k <= n: s += ((-1)**(k-1))*k*e[k]
        P[k] = s
    def Cc(k): return 2**k + 1 - P[k]
    def mu(m):
        if m == 1: return 1
        r, mm, p = 1, m, 2
        while p*p <= mm:
            if mm % p == 0:
                mm //= p
                if mm % p == 0: return 0
                r = -r
            p += 1
        if mm > 1: r = -r
        return r
    a = {d: sum(mu(d//ee)*Cc(ee) for ee in range(1, d+1) if d % ee == 0)//d for d in range(1, tmax+1)}

    from math import comb
    def smooth_count(t):
        s = [0]*(g+1); s[0] = 1
        for d in range(1, t+1):
            ad = a[d]
            if ad == 0: continue
            factor = [0]*(g+1)                         # (1-T^d)^{-a_d} = sum_j C(ad+j-1,j) T^{dj}
            j = 0
            while d*j <= g:
                factor[d*j] = comb(ad+j-1, j)
                j += 1
            ns = [0]*(g+1)
            for i in range(g+1):
                if s[i]:
                    si = s[i]
                    for k in range(0, g+1-i, 1):
                        if factor[k]:
                            ns[i+k] += si*factor[k]
            s = ns
        return s[g]

    Bg = B[g]
    smooth_test_ops = g*g        # poly(g) smoothness/factoring cost per trial (~g^2), a conservative unit
    print(f"{'t':>3} {'N(t)=|FB|':>12} {'P_smooth':>11} {'rel trials':>11} {'+smooth op':>11} {'sparse LA':>10} {'max(op)':>8}")
    best = None
    for t in range(1, tmax+1):
        N = sum(a[d] for d in range(1, t+1))
        if N <= 0: continue
        sm = smooth_count(t)
        if sm <= 0: continue
        Ps = sm/Bg
        rel_trials = math.log2(N/Ps)
        rel_ops = rel_trials + math.log2(smooth_test_ops)
        la = math.log2(N) + math.log2(g)
        mx = max(rel_ops, la)
        if best is None or mx < best[-1]: best = (t, N, Ps, rel_trials, rel_ops, la, mx)
        if t % 2 == 0 or t <= 12:
            print(f"{t:>3} {N:>12} {Ps:>11.2e} {rel_trials:>11.2f} {rel_ops:>11.2f} {la:>10.2f} {mx:>8.2f}")
    t, N, Ps, rt, ro, la, mx = best
    print()
    print(f"OPTIMAL t={t}: |FB|={N}, relation {rt:.1f} trials (log2), with ~g^2 smooth test {ro:.1f} ops,"
          f" sparse LA {la:.1f}")
    print(f"=> IC cost ~ 2^{rt:.0f} trials  (~2^{mx:.0f} operations with the smoothness-test factor)"
          f"  vs Pollard rho ~ 2^61")
    print(f"   gap below rho: {61-rt:.0f} bits (trials) / {61-mx:.0f} bits (operations) -- IF A is a Jacobian")

if __name__ == "__main__":
    main()
