#!/usr/bin/env python3
"""ell = 131, PR #57 module I = O_K + Z v/DEG (../../explicit_ppav.json, read only).
The inclusion O_K -> I is an F_2-isogeny from the O_K-variety to a principally polarized one, with
kernel I/O_K = Z/DEG cyclic (DEG = 263 * p2).  Recompute independently (own arithmetic in
Z[zeta_131][alpha]) the scalar c with pi*v = c*v mod DEG*O_K and pibar*v = c'*v, compare with the JSON's
recorded value and with lambda from residues.gp / residues_check.py, and get the field of definition of the
kernel points: F_{2^m}, m = lcm(ord_263(c), ord_p2(c)).  Also [Z[alpha]I : I] (needed for the cover bound)."""
import json, pathlib, math
here = pathlib.Path(__file__).resolve().parent
d = json.load(open(here.parents[1] / "explicit_ppav.json"))
ell, n = 131, 130
DEG = int(d["lattice"]["degree"]); v = [int(x) for x in d["lattice"]["generator_numerator"]]
p1, p2 = 263, 146505763881528721
assert p1 * p2 == DEG
# element = (a, b) with a, b lists of length 130 (coeffs of zeta^j), meaning a + b*alpha
def zmul(a, b, mod):          # product in Z[zeta]/(Phi_131), coefficients mod 'mod'
    c = [0] * (2 * n)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y: c[i + j] += x * y
    # reduce: zeta^131 = 1, then zeta^130 = -(1+...+zeta^129)
    r = [0] * ell
    for k, x in enumerate(c): r[k % ell] += x
    top = r[n]
    return [(r[j] - top) % mod for j in range(n)]
def emul(x, y, mod):           # (a+b al)(c+d al), al^2 = -al - 2
    a, b = x; c, e = y
    ac, ae, bc, be = zmul(a, c, mod), zmul(a, e, mod), zmul(b, c, mod), zmul(b, e, mod)
    return ([(ac[j] - 2 * be[j]) % mod for j in range(n)], [(ae[j] + bc[j] - be[j]) % mod for j in range(n)])
zeta = [0] * n; zeta[1] = 1; one = [0] * n; one[0] = 1; zero = [0] * n
PI = (zero, zeta)                                            # alpha*zeta
# pibar = 2/pi = beta*zeta^-1, beta = -1 - alpha ; zeta^-1 = zeta^130 = -(1+...+zeta^129)
zinv = [-1] * n
PIB = ([-x for x in zinv], [-x for x in zinv])
ALPHA = (zero, one)
V = (v[:n], v[n:])
def prop(w, u, p):            # find c with w = c*u mod p (u nonzero mod p) ; None if not proportional
    uu = u[0] + u[1]; ww = w[0] + w[1]
    k = next(i for i, x in enumerate(uu) if x % p)
    c = ww[k] * pow(uu[k], -1, p) % p
    return c if all((ww[i] - c * uu[i]) % p == 0 for i in range(len(uu))) else None
out = {}
for p in (p1, p2):
    Vp = ([x % p for x in V[0]], [x % p for x in V[1]])
    cp = prop(emul(PI, Vp, p), Vp, p); cbp = prop(emul(PIB, Vp, p), Vp, p); ap = prop(emul(ALPHA, Vp, p), Vp, p)
    zp = prop(emul(([0]+[0]*0 if False else zeta, zero), Vp, p), Vp, p)
    out[p] = (cp, cbp, ap, zp)
def order(x, p, fac):
    o = p - 1
    for q in fac:
        while o % q == 0 and pow(x, o // q, p) == 1: o //= q
    return o
f1 = [2, 131]; f2 = [2, 3, 5, 17, 31, 4001, 289510489]
assert math.prod([2,131]) == p1 - 1 and 16*3*5*17*31*4001*289510489 == p2 - 1
c1, c2 = out[p1][0], out[p2][0]
js = int(d["lattice"]["F_and_V_on_cyclic_quotient"])
print("pi on I/O_K:  mod 263: c=%s  pibar: %s  alpha-stable: %s  zeta-stable: %s" % (out[p1][0], out[p1][1], out[p1][2] is not None, out[p1][3] is not None))
print("              mod p2 : c=%s  pibar: %s  alpha-stable: %s  zeta-stable: %s" % (out[p2][0], out[p2][1], out[p2][2] is not None, out[p2][3] is not None))
print("JSON F_and_V_on_cyclic_quotient = %d ; mod 263 = %d, mod p2 = %d ; agrees with recomputed: %s" % (js, js % p1, js % p2, (js % p1, js % p2) == (c1, c2)))
lam = 14975904886999252
print("c mod p2 == lambda (residues.gp / residues_check.py): %s ; c^2 == 2 mod p2: %s ; c^2 mod 263 == 2: %s" % (c2 == lam, c2 * c2 % p2 == 2, c1 * c1 % p1 == 2))
o1, o2 = order(c1, p1, f1), order(c2, p2, f2)
m = o1 * o2 // math.gcd(o1, o2)
print("ord_263(c) = %d ; ord_p2(c) = %d = (p2-1)/%d ; kernel points of O_K -> I rational over F_{2^m} iff %d | m" % (o1, o2, (p2 - 1) // o2, m))
print("negative control: the 263-part alone needs only F_{2^%d}" % o1)
print("[Z[alpha] I : I] = %s" % ("1" if out[p2][2] is not None else "p2 (alpha*v not proportional to v mod p2; proportional mod 263: %s)" % (out[p1][2] is not None)))
