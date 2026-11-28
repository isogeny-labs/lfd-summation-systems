#!/usr/bin/env python3
"""2-torsion tests on the explicit Deligne module (construct_ppav.py) and on P_A.

(i)  P_A mod 2 == x^130 * Phi_131(x) mod 2  (etale part of A[2] has Frobenius char poly Phi_131).
(ii) On the explicit lattice: T contains O_K with odd index 263*p2, so T/2T = O_K/2O_K as F-modules.
     Using the Frobenius matrix of construct_ppav (basis zeta^j, alpha*zeta^j), compute the stable rank of
     F^k mod 2 (= dim of the etale part) and check F - 1 is invertible on the etale part:
     A'[2](F_2) = 0 and H^1(F_2, A'[2](F2bar)) = coker(F-1) = 0  ->  a rational symmetric theta divisor is
     not obstructed (torsor under A'[2] is trivial).
Negative control: the same detector on E x E-type toy Frobenius with an F_2-rational 2-torsion point
(F = identity block) must report H^1 != 0.
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PA = json.loads((HERE.parents[1] / "P_A.json").read_text())["P_A_coeffs_low_to_high"]

def polmod2(c): return [x & 1 for x in c]
def pmul2(a, b):
    r = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y: r[i+j] ^= 1
    return r
phi = [1]*131
target = pmul2([0]*130+[1], phi)
out = {"P_A_mod2_equals_x^130*Phi131": polmod2(PA) == target}

# Frobenius on O_K basis (copied from construct_ppav.py): F = pi = alpha*zeta acting on (c, d) with x = c + alpha d
n = 260
F = [[0]*n for _ in range(n)]
for j in range(130):
    shift = [int(i == j+1) if j < 129 else -1 for i in range(130)]
    for i in range(130):
        F[130+i][j] = shift[i]
        F[i][130+j] = -2*shift[i]
        F[130+i][130+j] = -shift[i]
def tobits(M): return [sum(((M[i][j] & 1) << j) for j in range(len(M[0]))) for i in range(len(M))]
def mulbits(A, B, m):
    # rows as bitmasks; (A*B)[i] = XOR of B[k] for k in A[i]
    out = []
    for r in A:
        acc, k = 0, 0
        while r:
            if r & 1: acc ^= B[k]
            r >>= 1; k += 1
        out.append(acc)
    return out
def rank(rows, m):
    rows = list(rows); r = 0
    for col in range(m):
        piv = next((i for i in range(r, len(rows)) if (rows[i] >> col) & 1), None)
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(len(rows)):
            if i != r and (rows[i] >> col) & 1: rows[i] ^= rows[r]
        r += 1
    return r
def analyse(Fb, m):
    P = Fb; ranks = []
    for k in range(1, m+1):
        ranks.append(rank(P, m)); P = mulbits(P, Fb, m)
    et = ranks[-1]
    # etale part = image of F^m ; restrict F - 1 to it: rank of (F-1) F^m equals et iff invertible there
    Fm = Fb
    for _ in range(m-1): Fm = mulbits(Fm, Fb, m)
    I = [1 << i for i in range(m)]
    FmI = [a ^ b for a, b in zip(Fb, I)]
    return et, rank(mulbits(Fm, FmI, m), m)
Fb = tobits(F)
et, r = analyse(Fb, n)
out["explicit_lattice_mod2"] = {"etale_dim": et, "rank_(F-1)_on_etale": r,
                                "A'[2](F2)=0_and_H1=0": et == 130 and r == et}
# negative control: toy 4-dim with F = diag(1,1,0,0) (two rational 2-torsion points): H^1 != 0 expected
toy = tobits([[1,0,0,0],[0,1,0,0],[0,0,0,0],[0,0,0,0]])
et2, r2 = analyse(toy, 4)
out["negative_control_toy"] = {"etale_dim": et2, "rank_(F-1)_on_etale": r2, "H1_nonzero_detected": r2 < et2}
print(json.dumps(out, indent=1))
