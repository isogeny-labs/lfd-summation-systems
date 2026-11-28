#!/usr/bin/env python3
"""Independent (no PARI) certificate that zeta_ell is NOT in Z[pi, pibar] for A_ell.

O_K = Z[z]/Phi_ell (x) Z[w]/(w^2-w+2); alpha=w-1, beta=-w; pi=alpha z, pibar=beta z^-1.
Z[pi,pibar] = Z-span{pi^i (0<=i<2d), pibar^j (1<=j<2d)}, d=ell-1 (pi*pibar=2).
Certificate: for a prime p, if zeta mod p is not in the F_p-span of the generators mod p, then
zeta is not in Z[pi,pibar] + p O_K, hence not in Z[pi,pibar]. Pure integer arithmetic.
"""
import sys

def mul(a, b, ell, p):
    # elements: (c0, c1) polynomials in z (lists len ell, mod z^ell-1), value c0 + c1*w
    def pm(x, y):
        r = [0]*ell
        for i, xi in enumerate(x):
            if xi:
                for j, yj in enumerate(y):
                    if yj: r[(i+j) % ell] = (r[(i+j) % ell] + xi*yj) % p
        return r
    a0, a1 = a; b0, b1 = b
    t0, t1, t2 = pm(a0, b0), [(u+v) % p for u, v in zip(pm(a0, b1), pm(a1, b0))], pm(a1, b1)
    # w^2 = w - 2
    return ([(u - 2*v) % p for u, v in zip(t0, t2)], [(u+v) % p for u, v in zip(t1, t2)])

def coords(a, ell, p):
    # reduce mod Phi_ell: z^(ell-1) = -(1+z+...+z^(ell-2))
    out = []
    for c in a:
        c = c[:]; top = c[ell-1]
        r = [(c[i] - top) % p for i in range(ell-1)]
        out += r
    return out

def rank(rows, p):
    rows = [r[:] for r in rows]; rk = 0; ncol = len(rows[0])
    for col in range(ncol):
        piv = next((i for i in range(rk, len(rows)) if rows[i][col]), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        inv = pow(rows[rk][col], p-2, p)
        rows[rk] = [x*inv % p for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][col]:
                f = rows[i][col]; rows[i] = [(x - f*y) % p for x, y in zip(rows[i], rows[rk])]
        rk += 1
    return rk

def check(ell, p):
    d = ell-1
    z = ([0,1]+[0]*(ell-2), [0]*ell)
    one = ([1]+[0]*(ell-1), [0]*ell)
    alpha = ([p-1]+[0]*(ell-1), [1]+[0]*(ell-1))      # w - 1
    beta = ([0]*ell, [p-1]+[0]*(ell-1))                 # -w
    zinv = ([0]*(ell-1)+[1], [0]*ell)
    pi, pb = mul(alpha, z, ell, p), mul(beta, zinv, ell, p)
    gens, x = [], one
    for _ in range(2*d): gens.append(coords(x, ell, p)); x = mul(x, pi, ell, p)
    x = pb
    for _ in range(1, 2*d): gens.append(coords(x, ell, p)); x = mul(x, pb, ell, p)
    r0 = rank(gens, p); r1 = rank(gens + [coords(z, ell, p)], p)
    return r0, r1

if __name__ == "__main__":
    for ell, p in [(11, 23), (19, 457), (131, 263), (5, 3), (13, 3)]:
        r0, r1 = check(ell, p)
        print(f"ell={ell} p={p}: rank mod p of Z[pi,pibar] gens = {r0} / {2*(ell-1)};"
              f" with zeta = {r1}  => zeta NOT in Z[pi,pibar]: {r1 > r0}")
