#!/usr/bin/env python3
"""Small-ell test: is A_ell the Jacobian of a 4-branch-point cyclic cover of P^1? (corrected 2026-11-26)

For a cover y^ell = f(x), f = x(x-1)^{a2}(x-g)^{a3}(x-lam)^{a4}, sum a_i = 0 mod ell, over F_q, q=2^{ell-1}:
#C(F_q) = q+1+sum_m (S_m + eps), S_m = sum_{f(x)!=0} chi^m(f(x)), eps = 1 (infinity unramified, leading coeff 1).
So the chi^m-eigenspace trace is -(S_m+eps). S_m is the Galois conjugate of S_1, so the cover matches A_ell's
traces {s*zeta^-j} (s = E's trace over F_q) iff -(S_1+eps) = s*zeta^k, i.e. the histogram h' = h + eps*e_0
(h[j] = #{x : dlog f(x) = j mod ell}) is constant except one entry offset by exactly -s.

Earlier version (PR #70 first draft) reported the MINIMUM spread and called it the largest, looked only for a
spike up, and omitted eps. This version is sign-aware, includes infinity, and reports the actual trace integers
of every cover with the coherent 'integer x root of unity' trace structure.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collections import Counter
import numpy as np
from ct_lib import build, s_of

def run(ell, n):
    q, g, dlog, mul, pw = build(n)
    s = s_of(ell)
    xs = np.arange(q, dtype=np.int64)
    L = {b: dlog[xs ^ b] for b in (0, 1, g)}
    tried = 0; coherent = Counter(); smin, smax = 10**9, 0
    for lam in range(q):
        if lam in (0, 1, g): continue
        Ll = dlog[xs ^ lam]
        good = (L[0] >= 0) & (L[1] >= 0) & (L[g] >= 0) & (Ll >= 0)
        for a2 in range(1, ell):
            for a3 in range(1, ell):
                a4 = (-(1 + a2 + a3)) % ell
                if a4 == 0: continue
                tried += 1
                d = (L[0] + a2 * L[1] + a3 * L[g] + a4 * Ll) % ell
                h = np.bincount(d[good], minlength=ell).astype(np.int64); h[0] += 1
                srt = np.sort(h)
                smin = min(smin, int(srt[-1] - srt[0])); smax = max(smax, int(srt[-1] - srt[0]))
                if srt[1] == srt[-1]: delta = int(srt[0] - srt[1])
                elif srt[0] == srt[-2]: delta = int(srt[-1] - srt[0])
                else: continue
                coherent[-delta] += 1          # trace integer
    return q, s, tried, coherent, smin, smax

if __name__ == "__main__":
    for ell, n in [(5, 4), (11, 10), (13, 12)]:
        q, s, tried, coh, smin, smax = run(ell, n)
        print(f"ell={ell} q={q}: A needs trace integer {s}; covers={tried}; spread min/max={smin}/{smax}; "
              f"coherent-trace covers (trace integer: count)={dict(sorted(coh.items()))}; MATCH={coh.get(s, 0)}")
