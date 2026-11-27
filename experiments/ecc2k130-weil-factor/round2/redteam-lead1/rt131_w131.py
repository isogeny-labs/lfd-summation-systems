#!/usr/bin/env python3
"""Red-team (redteam-lead1): w131-adic valuations (the unique place of K above 131) of
   m = mnum/131 (PR #57 record) and of lambda_nat, by expanding in the uniformizer w = 1 - zeta:
   O_w = Z_131[alpha][w], Eisenstein (Phi_131(1-w)), so v(sum a_j w^j) = min_j (130 v_131(a_j) + j).
   Also the 2-adic CM-type statement at ell = 131 (pure arithmetic): v_2(sigma(pi)) for all sigma.
   Writes nothing except rt131_w131.json."""
import json
from pathlib import Path
import flint
HERE = Path(__file__).resolve().parent
REC = json.loads((HERE.parent.parent / "explicit_ppav.json").read_text())
ELL, N = 131, 130
phi = flint.fmpz_poly([1] * ELL)
onemx = flint.fmpz_poly([1, -1])
E = phi(onemx)                      # Phi_131(1 - w), Eisenstein at 131
assert int(E.coeffs()[-1]) == 1 and all(int(c) % 131 == 0 for c in E.coeffs()[:-1]) and int(E.coeffs()[0]) % (131 * 131) != 0
def vp(n, p):
    n = abs(int(n))
    if n == 0: return 10 ** 9
    k = 0
    while n % p == 0: n //= p; k += 1
    return k
def wval(c, d):
    C = flint.fmpz_poly(c)(onemx) % E
    D = flint.fmpz_poly(d)(onemx) % E
    cc = [int(t) for t in C.coeffs()] + [0] * N
    dd = [int(t) for t in D.coeffs()] + [0] * N
    return min(130 * min(vp(cc[j], 131), vp(dd[j], 131)) + j for j in range(N))
mnum = REC["polarization"]["m_numerator_coefficients_1_alpha_over_zeta"]
unit = REC["positive_unit"]["unit_coefficients_1_alpha_over_zeta"]
v_mnum = wval(mnum[0], mnum[1]); v_m = v_mnum - 130
v_unit = wval(unit[0], unit[1])
# lambda_nat * 917 = -(2 - z - 1/z)(2 alpha + 1) ; z^-1 = z^130 = -(1 + ... + z^129)
two_m = [2 + 1, -1] + [0] * (N - 2)          # 2 - z - z^-1 with z^-1 = -(1+z+...+z^129): 2 - z + 1 + z + ... = adjust below
zinv = [-1] * N
two_m = [2 - zinv[0]] + [(-1 if j == 1 else 0) - zinv[j] for j in range(1, N)]
c_nat = [-t for t in two_m]; d_nat = [-2 * t for t in two_m]
v_nat = wval(c_nat, d_nat) - 130             # minus v_w(917) = v_w(131) = 130 (7 is a unit at w)
# control: v_w(1 - zeta) must be 1, v_w(131) = 130
v_ctrl1 = wval([1, -1] + [0] * (N - 2), [0] * N)
v_ctrl2 = wval([131] + [0] * (N - 1), [0] * N)
print(f"controls: v_w(1-zeta) = {v_ctrl1} (expect 1), v_w(131) = {v_ctrl2} (expect 130)")
print(f"PR #57: v_w131(mnum) = {v_mnum}, so v_w131(m) = {v_m}; v_w131(unit) = {v_unit} (expect 0)")
print(f"natural polarization: v_w131(lambda_nat) = {v_nat} (expect -128, even)")
print(f"=> PR #57 lambda = lambda0 * conj(u m): v_w131(lambda) = {v_nat + v_m + v_unit} (odd => (r,-7)_w131 = -1 = eps)")
# CM type at ell = 131: sigma = (a, s) acts by zeta -> zeta^a, alpha -> alpha (s=+1) or beta (s=-1).
# v_2(sigma(pi)) = v_2(sigma(alpha)) since zeta is a unit; in Z_2, x^2+x+2 has roots of valuation 1 and 0.
def two_adic_roots(prec=40):
    m = 2 ** prec
    # alpha_2: root with v=1 ; beta_2 = -1 - alpha_2 (unit). Newton/Hensel from x = 2 (f(2)=8, f'(2)=5 odd)
    x = 2
    for _ in range(10):
        x = (x - (x * x + x + 2) * pow(2 * x + 1, -1, m)) % m
    return x, (-1 - x) % m
a2, b2 = two_adic_roots()
assert (a2 * a2 + a2 + 2) % 2 ** 40 == 0
ty = {s: vp(a2 if s == 1 else b2, 2) for s in (1, -1)}
print(f"CM type at ell=131: v_2(sigma(pi)) = {ty[1]} for all sigma=(a,+1), {ty[-1]} for all sigma=(a,-1): Phi = {{s = +1}} (induced from Q(sqrt-7))")
json.dump(dict(v_w131_mnum=v_mnum, v_w131_m=v_m, v_w131_unit=v_unit, v_w131_lambda_nat=v_nat,
               v_w131_lambda_pr57=v_nat + v_m + v_unit, controls=[v_ctrl1, v_ctrl2], cm_type_v2=ty),
          open(HERE / "rt131_w131.json", "w"), indent=1)
