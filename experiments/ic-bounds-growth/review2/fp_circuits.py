#!/usr/bin/env python3
"""Independent probe of Theorem 3.3 of docs/proofs/N8-growth.md in characteristic p, with
dependence over F^p (review 2, 2026-11-12).  Imports nothing from the repository.

Model: the plane curve C: S_3(x1, x2, t) = 0 over F_p (t != 0, R not 2-torsion), function field
F = Fbar_p(C), xi = (x1, x2).  d/dx1 is a derivation of F with kernel F^p (x1 is separating).
At a point P = (alpha, beta) of C with dS/dx2 != 0, s = x1 - alpha is a uniformizer and x2 = y(s)
is a power series (Newton iteration).  For r <= p functions f_j, the Wronskian with ordinary
derivatives at P is prod_{i<r} i! * det([s^i] f_j) (i < r < p, so i! is a unit); it is nonzero at
some point iff the f_j are independent over F^p (Wronskian criterion, ordinary derivatives,
kernel F^p).  A nonzero value at one point is a rigorous independence certificate; vanishing at
all tested points is only a candidate for dependence and is re-tested over F_{p^2}.

Reduction to residues: multiplying all monomials by xi^(-e0) and single ones by p-th powers
preserves F^p-(in)dependence, so F^p-circuits correspond to sets of exponent residues mod p.

Checks:
  m=2  every residue e != 0: xi^e not in F^p (Theorem 3.3 with m=2: bound 0).
  m=3  all pairs {e1, e2} of nonzero residues: is {1, xi^e1, xi^e2} F^p-dependent?  Every
       dependent triple must satisfy ||2 da||_p + ||2 db||_p <= 12 for its three pairs.
  m=4  random bases {0, e1, e2} with all pairwise spreads > 36, every residue e3 scanned:
       a dependent {0, e1, e2, e3} would contain a circuit violating Theorem 3.3.
  control: the nine monomials of S_3 are F^p-dependent, every eight of them independent.

Usage: python3 fp_circuits.py p a b t [NPOINTS] [M4_BASES] [SEED] [m3|nom3]
"""
import json
import random
import sys

import numpy as np


def S3(x1, x2, a, b, t, p):
    return ((x1 - x2) ** 2 * t * t - 2 * ((x1 + x2) * (x1 * x2 + a) + 2 * b) * t
            + (x1 * x2 - a) ** 2 - 4 * b * (x1 + x2)) % p


def dS3_dx2(x1, x2, a, b, t, p):
    return (-2 * (x1 - x2) * t * t - 2 * ((x1 * x2 + a) + (x1 + x2) * x1) * t
            + 2 * (x1 * x2 - a) * x1 - 4 * b) % p


# ---------- truncated power series over a ring given by (add, mul) on python ints mod p ----------
def ser_mul(A, B, n, p):
    C = [0] * n
    for i, ai in enumerate(A[:n]):
        if ai:
            for j in range(n - i):
                C[i + j] = (C[i + j] + ai * B[j]) % p
    return C


def ser_inv(A, n, p):
    inv0 = pow(A[0], p - 2, p)
    B = [0] * n
    B[0] = inv0
    for k in range(1, n):
        s = sum(A[i] * B[k - i] for i in range(1, k + 1)) % p
        B[k] = (-s * inv0) % p
    return B


def local_series(alpha, beta, a, b, t, p, n):
    """x2 = y(s) with S3(alpha + s, y(s)) = 0 mod s^n, y(0) = beta (Newton iteration)."""
    X1 = [alpha % p, 1] + [0] * (n - 2)
    Y = [beta % p] + [0] * (n - 1)
    tt = [t % p] + [0] * (n - 1)
    aa = [a % p] + [0] * (n - 1)
    bb = [b % p] + [0] * (n - 1)

    def add(*xs):
        return [sum(v) % p for v in zip(*xs)]

    def neg(x):
        return [(-v) % p for v in x]

    def sc(c, x):
        return [(c * v) % p for v in x]

    for _ in range(n + 1):
        d = add(X1, neg(Y))
        s_ = add(X1, Y)
        prod = ser_mul(X1, Y, n, p)
        val = add(ser_mul(ser_mul(d, d, n, p), ser_mul(tt, tt, n, p), n, p),
                  neg(sc(2, ser_mul(add(ser_mul(s_, add(prod, aa), n, p), sc(2, bb)), tt, n, p))),
                  ser_mul(add(prod, neg(aa)), add(prod, neg(aa)), n, p),
                  neg(sc(4, ser_mul(bb, s_, n, p))))
        der = add(sc(-2, ser_mul(d, ser_mul(tt, tt, n, p), n, p)),
                  sc(-2, ser_mul(add(add(prod, aa), ser_mul(s_, X1, n, p)), tt, n, p)),
                  sc(2, ser_mul(add(prod, neg(aa)), X1, n, p)),
                  neg(sc(4, bb)))
        Y = add(Y, neg(ser_mul(val, ser_inv(der, n, p), n, p)))
    return X1, Y


def residue_table(X1, Y, p, n):
    """F[a, b, i] = [s^i] (x1^a x2^b), a, b in [0, p)."""
    P1 = np.zeros((p, n), dtype=np.int64)
    P2 = np.zeros((p, n), dtype=np.int64)
    cur1 = [1] + [0] * (n - 1)
    cur2 = [1] + [0] * (n - 1)
    for e in range(p):
        P1[e] = cur1
        P2[e] = cur2
        cur1 = ser_mul(cur1, X1, n, p)
        cur2 = ser_mul(cur2, Y, n, p)
    F = np.zeros((p, p, n), dtype=np.int64)
    for i in range(n):
        acc = np.zeros((p, p), dtype=np.int64)
        for u in range(i + 1):
            acc = (acc + np.outer(P1[:, u], P2[:, i - u])) % p
        F[:, :, i] = acc
    return F


def dist(n, p):
    r = n % p
    return min(r, p - r)


def spread(e, f, p):
    return dist(2 * (e[0] - f[0]), p) + dist(2 * (e[1] - f[1]), p)


def det_mod(M, p):
    M = [list(map(int, row)) for row in M]
    n = len(M)
    det = 1
    for c in range(n):
        piv = next((r for r in range(c, n) if M[r][c] % p), None)
        if piv is None:
            return 0
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            det = -det
        det = det * M[c][c] % p
        inv = pow(M[c][c], p - 2, p)
        for r in range(c + 1, n):
            f = M[r][c] * inv % p
            if f:
                M[r] = [(x - f * y) % p for x, y in zip(M[r], M[c])]
    return det % p


# ---------- F_{p^2} slow path for candidate dependences ----------
class Fp2:
    def __init__(self, p):
        self.p = p
        self.nr = next(c for c in range(2, p) if pow(c, (p - 1) // 2, p) == p - 1)

    def mul(self, x, y):
        p = self.p
        return ((x[0] * y[0] + self.nr * x[1] * y[1]) % p, (x[0] * y[1] + x[1] * y[0]) % p)

    def add(self, x, y):
        return ((x[0] + y[0]) % self.p, (x[1] + y[1]) % self.p)

    def inv(self, x):
        p = self.p
        nrm = (x[0] * x[0] - self.nr * x[1] * x[1]) % p
        ni = pow(nrm, p - 2, p)
        return (x[0] * ni % p, (-x[1]) * ni % p)


def fp2_wronskian_nonzero(p, a, b, t, exps, rng, tries=40):
    """True if some point of C(F_{p^2}) gives a nonzero Wronskian of {xi^e : e in exps}."""
    K = Fp2(p)
    r = len(exps)
    n = r
    elems = [(u, v) for u in range(p) for v in range(p)]
    zero = (0, 0)

    def poly_eval_S(x1, x2):
        # evaluate S3 in F_{p^2}
        def sc(c, x):
            return ((c * x[0]) % p, (c * x[1]) % p)
        mul, add = K.mul, K.add
        d = add(x1, sc(-1, x2))
        s_ = add(x1, x2)
        pr = mul(x1, x2)
        v = mul(mul(d, d), ((t * t) % p, 0))
        v = add(v, sc(-2 * t, add(mul(s_, add(pr, (a % p, 0))), ((2 * b) % p, 0))))
        q = add(pr, ((-a) % p, 0))
        v = add(v, mul(q, q))
        v = add(v, sc(-4 * b, s_))
        return v

    def poly_eval_Sy(x1, x2):
        def sc(c, x):
            return ((c * x[0]) % p, (c * x[1]) % p)
        mul, add = K.mul, K.add
        d = add(x1, sc(-1, x2))
        pr = mul(x1, x2)
        v = sc(-2 * t * t, d)
        v = add(v, sc(-2 * t, add(add(pr, (a % p, 0)), mul(add(x1, x2), x1))))
        v = add(v, sc(2, mul(add(pr, ((-a) % p, 0)), x1)))
        v = add(v, ((-4 * b) % p, 0))
        return v

    def smul(A, B):
        C = [zero] * n
        for i in range(n):
            if A[i] != zero:
                for j in range(n - i):
                    C[i + j] = K.add(C[i + j], K.mul(A[i], B[j]))
        return C

    def sadd(*xs):
        out = [zero] * n
        for x in xs:
            out = [K.add(u, v) for u, v in zip(out, x)]
        return out

    def ssc(c, x):
        return [((c * u[0]) % p, (c * u[1]) % p) for u in x]

    def sinv(A):
        i0 = K.inv(A[0])
        B = [zero] * n
        B[0] = i0
        for k in range(1, n):
            s = zero
            for i in range(1, k + 1):
                s = K.add(s, K.mul(A[i], B[k - i]))
            B[k] = K.mul(((-s[0]) % p, (-s[1]) % p), i0)
        return B

    for _ in range(tries):
        alpha = rng.choice(elems)
        roots = [y for y in elems if poly_eval_S(alpha, y) == zero and poly_eval_Sy(alpha, y) != zero]
        if not roots:
            continue
        beta = rng.choice(roots)
        X1 = [alpha, (1, 0)] + [zero] * (n - 2)
        Y = [beta] + [zero] * (n - 1)
        tt = [(t % p, 0)] + [zero] * (n - 1)
        aa = [(a % p, 0)] + [zero] * (n - 1)
        bb = [(b % p, 0)] + [zero] * (n - 1)
        for _ in range(n + 1):
            d = sadd(X1, ssc(-1, Y))
            s_ = sadd(X1, Y)
            pr = smul(X1, Y)
            val = sadd(smul(smul(d, d), smul(tt, tt)), ssc(-2, smul(sadd(smul(s_, sadd(pr, aa)), ssc(2, bb)), tt)),
                       smul(sadd(pr, ssc(-1, aa)), sadd(pr, ssc(-1, aa))), ssc(-4, smul(bb, s_)))
            der = sadd(ssc(-2, smul(d, smul(tt, tt))), ssc(-2, smul(sadd(sadd(pr, aa), smul(s_, X1)), tt)),
                       ssc(2, smul(sadd(pr, ssc(-1, aa)), X1)), ssc(-4, bb))
            Y = sadd(Y, ssc(-1, smul(val, sinv(der))))
        cols = []
        for (ea, eb) in exps:
            f = [(1, 0)] + [zero] * (n - 1)
            for _ in range(ea):
                f = smul(f, X1)
            for _ in range(eb):
                f = smul(f, Y)
            cols.append(f)
        # determinant over F_{p^2} of the r x r matrix [s^i] f_j by Gaussian elimination
        M = [[cols[j][i] for j in range(r)] for i in range(r)]
        nonzero = True
        for c in range(r):
            piv = next((rr for rr in range(c, r) if M[rr][c] != zero), None)
            if piv is None:
                nonzero = False
                break
            M[c], M[piv] = M[piv], M[c]
            iv = K.inv(M[c][c])
            for rr in range(c + 1, r):
                f = K.mul(M[rr][c], iv)
                if f != zero:
                    M[rr] = [K.add(x, K.mul(((-f[0]) % p, (-f[1]) % p), y)) for x, y in zip(M[rr], M[c])]
        if nonzero:
            return True
    return False


def main():
    p, a, b, t = (int(x) for x in sys.argv[1:5])
    npts = int(sys.argv[5]) if len(sys.argv) > 5 else 12
    m4_bases = int(sys.argv[6]) if len(sys.argv) > 6 else 30
    seed = int(sys.argv[7]) if len(sys.argv) > 7 else 1
    do_m3 = (sys.argv[8] != "nom3") if len(sys.argv) > 8 else True
    rng = random.Random(seed)
    assert t % p and (4 * a ** 3 + 27 * b ** 2) % p and (t ** 3 + a * t + b) % p
    n = 9  # series length: derivatives up to order 8 (control with nine terms)
    pts = [(al, be) for al in range(p) for be in range(p)
           if S3(al, be, a, b, t, p) == 0 and dS3_dx2(al, be, a, b, t, p) != 0]
    rng.shuffle(pts)
    pts = pts[:npts]
    tables = []
    for (al, be) in pts:
        X1, Y = local_series(al, be, a, b, t, p, n)
        assert all(v == 0 for v in [S3(al, be, a, b, t, p)])
        tables.append(residue_table(X1, Y, p, n))
    rep = {"p": p, "a": a, "b": b, "t": t, "points": len(pts)}

    # control: the nine monomials of S_3 (exponents {0,1,2}^2) are dependent; every eight independent
    E9 = [(u, v) for u in range(3) for v in range(3)]

    def wr(exps, T):
        M = [[int(T[e][i]) for e in exps] for i in range(len(exps))]
        return det_mod(M, p)
    rep["control_S3_nine_dependent"] = all(wr(E9, T) == 0 for T in tables)
    rep["control_every_eight_independent"] = all(any(wr([e for e in E9 if e != x], T) for T in tables) for x in E9)

    # m = 2
    nonzero_res = [(u, v) for u in range(p) for v in range(p) if (u, v) != (0, 0)]
    m2_fail = [e for e in nonzero_res if all(T[e][1] == 0 for T in tables)]
    rep["m2_residues_possibly_in_Fp"] = m2_fail

    # m = 3: {1, f, g} dependent iff f' g'' - f'' g' vanishes identically (row chunks)
    R = p * p
    cand = []
    if do_m3:
        C1 = [T[:, :, 1].reshape(R) for T in tables]
        C2 = [T[:, :, 2].reshape(R) for T in tables]
        chunk = max(1, 4000000 // R)
        for lo in range(1, R, chunk):
            hi = min(R, lo + chunk)
            mask = np.ones((hi - lo, R), dtype=bool)
            rows = np.arange(lo, hi)
            mask &= (np.arange(R)[None, :] > rows[:, None])
            for c1, c2 in zip(C1, C2):
                M = (np.outer(c1[lo:hi], c2) - np.outer(c2[lo:hi], c1)) % p
                mask &= (M == 0)
                if not mask.any():
                    break
            for (i, j) in np.argwhere(mask):
                cand.append((lo + int(i), int(j)))
    rep["m3_exhaustive"] = do_m3
    rep["m3_candidates_after_Fp_points"] = int(len(cand))
    confirmed, refuted, viol = [], 0, []
    for (i, j) in cand[:200]:
        e1 = (int(i) // p, int(i) % p)
        e2 = (int(j) // p, int(j) % p)
        if fp2_wronskian_nonzero(p, a, b, t, [(0, 0), e1, e2], rng):
            refuted += 1
        else:
            confirmed.append((e1, e2))
            sp = max(spread((0, 0), e1, p), spread((0, 0), e2, p), spread(e1, e2, p))
            if sp > 12:
                viol.append((e1, e2, sp))
    rep["m3_candidates_refuted_over_Fp2"] = refuted
    rep["m3_dependent_triples"] = confirmed
    rep["m3_violations"] = viol

    # m = 4 extension scan
    far = 36
    bases, tries = [], 0
    while len(bases) < m4_bases and tries < 100000:
        tries += 1
        e1 = rng.choice(nonzero_res)
        e2 = rng.choice(nonzero_res)
        if e1 == e2:
            continue
        if min(spread((0, 0), e1, p), spread((0, 0), e2, p), spread(e1, e2, p)) <= far:
            continue
        bases.append((e1, e2))
    scanned, dep_ext = 0, []
    for (e1, e2) in bases:
        alive = np.ones(R, dtype=bool)
        alive[0] = False
        alive[e1[0] * p + e1[1]] = False
        alive[e2[0] * p + e2[1]] = False
        for T in tables:
            C = T.reshape(R, n)
            u, v = C[e1[0] * p + e1[1]], C[e2[0] * p + e2[1]]
            # det [[u1,v1,w1],[u2,v2,w2],[u3,v3,w3]] expanded along the w column
            m1 = (u[2] * v[3] - u[3] * v[2]) % p
            m2 = (u[1] * v[3] - u[3] * v[1]) % p
            m3 = (u[1] * v[2] - u[2] * v[1]) % p
            d = (C[:, 1] * m1 - C[:, 2] * m2 + C[:, 3] * m3) % p
            alive &= (d == 0)
        scanned += R - 3
        for idx in np.flatnonzero(alive):
            e3 = (int(idx) // p, int(idx) % p)
            if not fp2_wronskian_nonzero(p, a, b, t, [(0, 0), e1, e2, e3], rng):
                dep_ext.append((e1, e2, e3))
    rep["m4_bases"] = len(bases)
    rep["m4_extensions_scanned"] = scanned
    rep["m4_dependent_extensions"] = dep_ext
    print(json.dumps(rep))
    ok = (rep["control_S3_nine_dependent"] and rep["control_every_eight_independent"] and not m2_fail
          and not viol and not dep_ext)
    print("NO VIOLATION FOUND" if ok else "CHECK FAILED OR VIOLATION FOUND")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
