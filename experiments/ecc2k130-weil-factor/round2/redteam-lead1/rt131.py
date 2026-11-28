#!/usr/bin/env python3
"""Red-team (redteam-lead1): independent re-verification at ell = 131 of
  (1) the natural restricted-product polarization lambda_nat = (2 - z - 1/z)/(131 sqrt-7):
      its Gram on the basis (z^j, alpha z^j) is +-[[0,G],[-G,0]] (G = Cartan A_130), det = 131^2;
  (2) the explicit ppav of PR #57 (../../explicit_ppav.json, read-only):
      rebuild lambda = lambda0 * conj(u*mnum)/131 from the recorded unit u and numerator mnum,
      build the Gram of Tr(lambda x ybar) on I = O_K + Z v/D OURSELVES (exact integers),
      check integral / alternating / det = 1 (python-flint exact det), pi- and pibar-stability, positivity;
  (3) the local bookkeeping of lead 1's Theorem A on this certificate:
      p2-adic valuations of m at all 130 places above p2 (self-glued v0 must carry valuation 1, odd),
      263-adic valuations at the 260 places above 263, and where v/D sits (End non-maximal exactly at v0).
Writes rt131.json next to this file.  No other file is written.
"""
import json, sys, time
from fractions import Fraction
from pathlib import Path

import flint
import mpmath

HERE = Path(__file__).resolve().parent
REC = json.loads((HERE.parent.parent / "explicit_ppav.json").read_text())
ELL = 131
N = ELL - 1          # 130
P2 = 146505763881528721

# ---------- arithmetic in O_K = Z[z]/(Phi_131) [alpha]/(alpha^2+alpha+2); element = (c, d) = c + d*alpha ----------
def zmul(a, b):
    """product of two length-130 integer lists modulo Phi_131 (z^130 = -(1+...+z^129))"""
    fa = flint.fmpz_poly(a); fb = flint.fmpz_poly(b)
    pr = fa * fb
    coeffs = [int(c) for c in pr.coeffs()] + [0] * (2 * N)
    out = [0] * (ELL)
    for k, c in enumerate(coeffs[:2 * N + 1]):
        out[k % ELL] += c
    top = out[N]
    return [out[j] - top for j in range(N)]

def zadd(a, b): return [x + y for x, y in zip(a, b)]
def zsub(a, b): return [x - y for x, y in zip(a, b)]
def zscal(s, a): return [s * x for x in a]
def zmono(k):
    k %= ELL
    if k == N: return [-1] * N
    v = [0] * N; v[k] = 1; return v
ONE = zmono(0); ZERO = [0] * N

def kmul(x, y):
    c1, d1 = x; c2, d2 = y
    dd = zmul(d1, d2)
    return (zsub(zmul(c1, c2), zscal(2, dd)), zsub(zadd(zmul(c1, d2), zmul(d1, c2)), dd))

def zconj(a):
    """z -> z^-1 on a length-130 list"""
    out = [0] * ELL
    for j, c in enumerate(a):
        out[(-j) % ELL] += c
    top = out[N]
    return [out[j] - top for j in range(N)]

def kconj(x):
    """z -> 1/z, alpha -> beta = -1 - alpha"""
    c, d = x
    cc, dc = zconj(c), zconj(d)
    return (zsub(cc, dc), zscal(-1, dc))

def coords(x): return list(x[0]) + list(x[1])

# Tr_{Q(z)/Q}(z^k) = 130 if k = 0 mod 131 else -1 ; Tr_{Q(alpha)/Q}: 1 -> 2, alpha -> -1
def ztrace(a): return N * a[0] - sum(a[1:])
def ktrace(x):
    c, d = x
    return 2 * ztrace(c) - ztrace(d)

# ---------- (1) natural polarization ----------
t0 = time.time()
z, zi = zmono(1), zmono(-1)
two_minus = zsub(zsub(zscal(2, ONE), z), zi)              # 2 - z - 1/z
sqrtm7 = (ONE, zscal(2, ONE))                              # 2 alpha + 1 = sqrt-7
# lambda_nat * 917 = (2 - z - 1/z) * (-(2 alpha + 1))  since 1/sqrt-7 = -(2alpha+1)/7 ; lambda_nat = that / (7*131)
lnat917 = kmul((two_minus, ZERO), (zscal(-1, ONE), zscal(-2, ONE)))
basis = [(zmono(j), ZERO) for j in range(N)] + [(ZERO, zmono(j)) for j in range(N)]
cbasis = [kconj(b) for b in basis]
def gram_on_OK(lam_num, den):
    """Gram[i][j] = Tr(lam e_i conj(e_j)), lam = lam_num/den.  With x = lam e_i = c + d alpha:
       Tr(x z^-b) = 131 (2c-d)_b - sum(2c-d) ;  Tr(x beta z^-b) = 131 (4d-c)_b - sum(4d-c)."""
    rows = []
    for i in range(2 * N):
        c, d = kmul(lam_num, basis[i])
        u1 = [2 * a - b for a, b in zip(c, d)]; s1 = sum(u1)
        u2 = [4 * b - a for a, b in zip(c, d)]; s2 = sum(u2)
        row = [Fraction(ELL * u1[b] - s1, den) for b in range(N)] + [Fraction(ELL * u2[b] - s2, den) for b in range(N)]
        rows.append(row)
    return rows
def gram_check_slow(lam_num, den, pairs):
    return all(gram_on_OK(lam_num, den)[i][j] == Fraction(ktrace(kmul(kmul(lam_num, basis[i]), cbasis[j])), den) for i, j in pairs)
G_nat = gram_on_OK(lnat917, 917)
print('    fast trace formula == direct trace on sample entries:', all(G_nat[i][j] == Fraction(ktrace(kmul(kmul(lnat917, basis[i]), cbasis[j])), 917) for i, j in [(0,0),(0,1),(3,130),(130,3),(200,77),(259,259),(5,135),(140,9)]))
cartan_ok = True
for a in range(N):
    for b in range(N):
        g = 2 if a == b else (-1 if abs(a - b) == 1 else 0)
        if G_nat[a][b] != 0 or G_nat[N + a][N + b] != 0: cartan_ok = False
        if G_nat[a][N + b] != -g or G_nat[N + a][b] != g: cartan_ok = False
Gn = flint.fmpz_mat([[int(x) for x in row] for row in G_nat]) if all(x.denominator == 1 for row in G_nat for x in row) else None
det_nat = int(Gn.det()) if Gn is not None else None
print(f"(1) natural polarization lambda_nat=(2-z-1/z)/(131 sqrt-7): Gram == [[0,-G],[G,0]] (G=Cartan A_130): {cartan_ok}; "
      f"integral {Gn is not None}; det = {det_nat} (131^2 = {131**2}) [{time.time()-t0:.1f}s]")

# ---------- (2) PR #57 certificate rebuilt from its recorded unit and m ----------
pol = REC["polarization"]; lat = REC["lattice"]
mnum = tuple(pol["m_numerator_coefficients_1_alpha_over_zeta"])
unit = tuple(REC["positive_unit"]["unit_coefficients_1_alpha_over_zeta"])
assert pol["m_denominator"] == 131
r = kmul(unit, mnum)                                    # r/131 = u*m
lam0_917 = kmul((zscal(-1, ONE), ZERO), lnat917)        # lambda0 = -lambda_nat  (Gram [[0,G],[-G,0]] as in PR #57)
rbar = kconj(r)
lam_num = kmul(lam0_917, rbar)                          # lambda = lam_num / (917*131)
den = 917 * 131
print("    r real (conj r == r):", kconj(r) == r, " ; lambda antisymmetric:", kconj(lam_num) == tuple(zscal(-1, c) for c in lam_num))
t0 = time.time()
Om = gram_on_OK(lam_num, den)
om_int = all(x.denominator == 1 for row in Om for x in row)
D = lat["degree"]; piv = lat["pivot"]; v = lat["generator_numerator"]
assert D == 263 * P2 and v[piv] == 1
OmI = [[int(x) for x in row] for row in Om] if om_int else None
# Gram on I: basis = std basis with column piv replaced by v/D
Omv = [sum(OmI[i][k] * v[k] for k in range(2 * N)) for i in range(2 * N)]   # Omega * v
vOm = [sum(v[k] * OmI[k][j] for k in range(2 * N)) for j in range(2 * N)]   # v^T Omega
GI = [row[:] for row in OmI]
integral_I = True
for i in range(2 * N):
    if Omv[i] % D or vOm[i] % D: integral_I = False
    GI[i][piv] = Fraction(Omv[i], D); GI[piv][i] = Fraction(vOm[i], D)
vOmv = sum(v[k] * Omv[k] for k in range(2 * N))
GI[piv][piv] = Fraction(vOmv, D * D)
if integral_I and vOmv % (D * D) == 0:
    GI = [[int(x) for x in row] for row in GI]
    M = flint.fmpz_mat(GI)
    detI = int(M.det())
    alt = all(GI[i][j] == -GI[j][i] for i in range(2 * N) for j in range(2 * N))
else:
    detI, alt = None, None
detO = int(flint.fmpz_mat(OmI).det()) if om_int else None
print(f"(2) PR #57 lattice I = O_K + Z v/D, D = 263*p2: Gram on O_K integral {om_int}, det on O_K = {detO} (D^2 = {D*D}: {detO == D*D}); "
      f"Gram on I integral {integral_I}, alternating {alt}, det = {detI} [{time.time()-t0:.1f}s]")

# pi- and pibar-stability of I: pi*v = s*v mod D for an integer s (v[piv] = 1)
pi = (ZERO, zmono(1)); pib = kmul((zscal(-1, ONE), zscal(-1, ONE)), (zmono(-1), ZERO))   # beta * z^-1
assert kmul(pi, pib) == (zscal(2, ONE), ZERO)
vv = (v[:N], v[N:])
def stab(x):
    w = coords(kmul(x, vv))
    s = w[piv] % D
    return all((w[k] - s * v[k]) % D == 0 for k in range(2 * N)), s
sp, s_pi = stab(pi); sb, s_pib = stab(pib); sz, _ = stab((zmono(1), ZERO))
print(f"    I stable under pi: {sp} (scalar {s_pi}, s^2 = 2 mod D: {s_pi*s_pi % D == 2}), under pibar: {sb}, under zeta: {sz}")

# positivity: Im phi_k(lambda) at phi_k: z -> exp(2 pi i k/131), sqrt-7 -> +i sqrt7 (Phi_+), k = 1..130
mpmath.mp.dps = 700   # coefficients of lam_num reach ~400 bits: need >> 120 digits against cancellation
alpha_c = (-1 + mpmath.mpc(0, mpmath.sqrt(7))) / 2
def emb(x, k):
    zz = mpmath.exp(2j * mpmath.pi * k / ELL)
    c, d = x
    cv = mpmath.polyval(list(reversed(c)), zz); dv = mpmath.polyval(list(reversed(d)), zz)
    return cv + dv * alpha_c
ims = [mpmath.im(emb(lam_num, k)) / den for k in range(1, ELL)]
npos = sum(1 for t in ims if t > 0); nneg = sum(1 for t in ims if t < 0)
print(f"    positivity on Phi_+: #Im>0 = {npos}, #Im<0 = {nneg}, min |Im| = {mpmath.nstr(min(abs(t) for t in ims), 5)}")

# ---------- (3) local bookkeeping of Theorem A on this certificate ----------
# Z_q = Z_{p2}[alpha] (alpha^2+alpha+2 irreducible mod p2), work mod p2^2.  gamma = beta/alpha = (alpha-1)/2
def qmul(x, y, m):
    a, b = x; c, d = y
    return ((a * c - 2 * b * d) % m, (a * d + b * c - b * d) % m)
def qpow(x, e, m):
    r_ = (1, 0)
    while e:
        if e & 1: r_ = qmul(r_, x, m)
        x = qmul(x, x, m); e >>= 1
    return r_
def zeval(cs, w, m):
    acc = (0, 0)
    for c in reversed(cs):
        acc = qmul(acc, w, m); acc = ((acc[0] + c) % m, acc[1])
    return acc
def keval(x, w, m):
    c, d = x
    cv = zeval(c, w, m); dv = zeval(d, w, m)
    return ((cv[0] + qmul(dv, (0, 1), m)[0]) % m, (cv[1] + qmul(dv, (0, 1), m)[1]) % m)
def qval(x, p, m):
    if x[0] % m == 0 and x[1] % m == 0: return ">=2"
    if x[0] % p == 0 and x[1] % p == 0: return 1
    return 0
m2 = P2 * P2
inv2 = pow(2, -1, m2)
gamma = ((-1 * inv2) % m2, inv2)                 # (alpha - 1)/2
x0 = qpow(gamma, 66, m2)
Nexp = P2 * pow(P2, -1, ELL)                     # = 1 mod 131, = 0 mod p2
omega = qpow(x0, Nexp, m2)
assert qpow(omega, ELL, m2) == (1, 0) and omega != (1, 0)
# check: v0 <-> k with omega^(2k) = gamma mod p2 ; with this omega, k = 1
selfk = [k for k in range(1, ELL) if (lambda t: (t[0] - gamma[0]) % P2 == 0 and (t[1] - gamma[1]) % P2 == 0)(qpow(omega, 2 * k, m2))]
m_val = {}; v_res = {}
for k in range(1, ELL):
    w = qpow(omega, k, m2)
    m_val[k] = qval(keval(mnum, w, m2), P2, m2)
    vv_k = keval(vv, w, m2)
    v_res[k] = (vv_k[0] % P2, vv_k[1] % P2) != (0, 0)
odd_places = [k for k in range(1, ELL) if m_val[k] == 1]
nonzero_v = [k for k in range(1, ELL) if v_res[k]]
print(f"(3) above p2: self-glued place v0 <-> k = {selfk}; places where v_p2(m) = 1: {odd_places}; "
      f"valuation multiset {sorted(set(map(str, m_val.values())))}; places where v/D is NOT integral (v != 0 mod place): {nonzero_v}")
# 263: split completely.  roots of x^2+x+2 mod 263 and a 131st root of unity mod 263^2
p3 = 263; m3 = p3 * p3
def sqrt_mod(a, p):
    for t in range(p):
        if (t * t - a) % p == 0: return t
s7 = sqrt_mod(-7 % p3, p3)
roots_a = [((-1 + s7) * pow(2, -1, p3)) % p3, ((-1 - s7) * pow(2, -1, p3)) % p3]
def hensel_alpha(a0):
    a = a0
    for _ in range(3):
        f = (a * a + a + 2) % m3; fp = (2 * a + 1) % m3
        a = (a - f * pow(fp, -1, m3)) % m3
    return a
g3 = next(g for g in range(2, p3) if pow(g, (p3 - 1) // ELL, p3) != 1)
rho0 = pow(g3, (p3 - 1) // ELL, p3)
rho = pow(rho0, p3 * pow(p3, -1, ELL), m3)
assert pow(rho, ELL, m3) == 1 and rho % p3 == rho0
cnt263 = {}; v263 = []
for ai, a0 in enumerate(roots_a):
    a = hensel_alpha(a0)
    for k in range(1, ELL):
        zz = pow(rho, k, m3)
        def ev(x):
            c, d = x
            cv = sum(cc * pow(zz, j, m3) for j, cc in enumerate(c)) % m3
            dv = sum(cc * pow(zz, j, m3) for j, cc in enumerate(d)) % m3
            return (cv + dv * a) % m3
        val = ev(mnum)
        vv_ = 0 if val % m3 == 0 else (1 if val % p3 == 0 else 0)
        if val % m3 == 0: vv_ = ">=2"
        cnt263[(ai, k)] = vv_
        if ev(vv) % p3: v263.append((ai, k))
ones263 = [key for key, t in cnt263.items() if t == 1]
print(f"    above 263 (split): places with v(m) = 1: {ones263} ; places where v/D not integral: {v263}")
res = dict(cartan_natural=cartan_ok, det_natural=det_nat, gram_OK_integral=om_int, det_on_OK=detO, det_on_I=detI,
           alternating_I=alt, pi_stable=sp, pibar_stable=sb, zeta_stable=sz, positivity=[npos, nneg],
           v0_index=selfk, p2_places_with_v_m_eq_1=odd_places, p2_places_v_nonintegral=nonzero_v,
           places263_v_m_eq_1=[list(t) for t in ones263], places263_v_nonintegral=[list(t) for t in v263])
(HERE / "rt131.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
print("written rt131.json")
