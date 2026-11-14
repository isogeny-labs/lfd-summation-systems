#!/usr/bin/env python3
"""Independent divisor computations on E: y^2 = x^3 + a x + b over F_{p^2} (review 2, 2026-11-12).
Imports nothing from the repository.

Engine: at every F_{p^2}-rational place Q of E (affine points and O) a uniformizer s is chosen
(x - x0 if y0 != 0; y if y0 = 0; z = -x/y at O) and x, y are expanded as Laurent series in s
(Newton iteration; at O the series w(z) = z^3 + a z w^2 + b w^3 of Silverman IV.1).  The functions
xi1 = x(P), xi2 = x(R - P) = lambda^2 - x - t, lambda = (-y - yR)/(x - t), are expanded and their
valuations read off.  Precision is tracked absolutely; a valuation that is not determined within
the precision raises an error.

  mode "lemma32"  (addendum (a)): div(xi1), div(xi2) at every rational place versus the predicted
                  (P0) + (-P0) - 2(O) and (R - P0) + (R + P0) - 2(R), for the configurations
                  generic / b = 0 / R 2-torsion / R = +-2P0 / R 2-torsion and b = 0 / t = 0 / t = b = 0,
                  each on two curves; reports ord_O(xi2), ord_R(xi1) and the orders of monomials.
  mode "twisted"  (addendum (b)): a relation sum_e phi_e xi^(m_e) = 0 with coefficients
                  phi_e = c_e h^p xi^(-p u_e) that are non-constant p-th powers (h has zeros and
                  poles outside S), m_e = e + p u_e, c_e the coefficients of S_3.  Checks: the relation
                  vanishes; {xi^(m_e)} is F^p-dependent with every 8 independent (Wronskian with
                  ordinary derivatives at points); rank 9 over the constants; ord(phi_e) = 0 mod p at
                  every place; H(v) <= n0 |S|; the pairwise bound and its residues mod p; and the local
                  Wronskian inequality of step (d) of Theorem 3.1, computed from Laurent series, at
                  every place of S, at the zeros and poles of h and at ramification points of x.

Usage: python3 divisor_check.py lemma32 [OUT.json]
       python3 divisor_check.py twisted p a b t [OUT.json]
"""
import json
import random
import sys


class Fq:
    """F_{p^2} = F_p[w]/(w^2 - nr); element u + v w encoded as u + p v."""

    def __init__(self, p):
        self.p = p
        q = p * p
        self.q = q
        nr = next(c for c in range(2, p) if pow(c, (p - 1) // 2, p) == p - 1)
        self.MUL = [0] * (q * q)
        self.ADD = [0] * (q * q)
        for z1 in range(q):
            u1, v1 = z1 % p, z1 // p
            base = z1 * q
            for z2 in range(q):
                u2, v2 = z2 % p, z2 // p
                self.MUL[base + z2] = (u1 * u2 + nr * v1 * v2) % p + p * ((u1 * v2 + u2 * v1) % p)
                self.ADD[base + z2] = (u1 + u2) % p + p * ((v1 + v2) % p)
        self.NEG = [(-(z % p)) % p + p * ((-(z // p)) % p) for z in range(q)]
        self.INV = [0] * q
        for z in range(1, q):
            u, v = z % p, z // p
            nrm = (u * u - nr * v * v) % p
            ni = pow(nrm, p - 2, p)
            self.INV[z] = (u * ni) % p + p * ((-v * ni) % p)
        self.SQRT = {}
        for z in range(q):
            self.SQRT.setdefault(self.MUL[z * q + z], []).append(z)

    def k(self, n):
        return n % self.p

    def m(self, x, y):
        return self.MUL[x * self.q + y]

    def a(self, x, y):
        return self.ADD[x * self.q + y]

    def s(self, x, y):
        return self.ADD[x * self.q + self.NEG[y]]

    def i(self, x):
        assert x != 0
        return self.INV[x]


class LS:
    """Laurent series t^v (c0 + c1 t + ...), coefficients known for exponents < v + len(c)."""
    __slots__ = ("v", "c")

    def __init__(self, v, c):
        self.v, self.c = v, c

    def prec(self):
        return self.v + len(self.c)


def norm(F, A):
    c = A.c
    i = 0
    while i < len(c) and c[i] == 0:
        i += 1
    return LS(A.v + i, c[i:])


def val(F, A):
    B = norm(F, A)
    if not B.c:
        raise ArithmeticError("valuation not determined within precision (>= %d)" % B.v)
    return B.v


def const(F, x, L):
    return LS(0, [x] + [0] * (L - 1))


def add(F, A, B):
    v = min(A.v, B.v)
    N = min(A.prec(), B.prec())
    out = []
    for e in range(v, N):
        x = A.c[e - A.v] if 0 <= e - A.v < len(A.c) else 0
        y = B.c[e - B.v] if 0 <= e - B.v < len(B.c) else 0
        out.append(F.a(x, y))
    return LS(v, out)


def neg(F, A):
    return LS(A.v, [F.NEG[x] for x in A.c])


def sub(F, A, B):
    return add(F, A, neg(F, B))


def scal(F, k, A):
    return LS(A.v, [F.m(k, x) for x in A.c])


def mul(F, A, B):
    A, B = norm(F, A), norm(F, B)
    n = min(len(A.c), len(B.c))
    out = [0] * n
    MUL, ADD, q = F.MUL, F.ADD, F.q
    for i in range(n):
        ai = A.c[i]
        if ai:
            base = ai * q
            for j in range(n - i):
                bj = B.c[j]
                if bj:
                    out[i + j] = ADD[out[i + j] * q + MUL[base + bj]]
    return LS(A.v + B.v, out)


def inv(F, A):
    A = norm(F, A)
    if not A.c:
        raise ArithmeticError("inverse of a series that is zero to precision")
    n = len(A.c)
    i0 = F.i(A.c[0])
    out = [0] * n
    out[0] = i0
    for k in range(1, n):
        s = 0
        for i in range(1, k + 1):
            if A.c[i] and out[k - i]:
                s = F.a(s, F.m(A.c[i], out[k - i]))
        out[k] = F.m(F.NEG[s], i0)
    return LS(-A.v, out)


def div(F, A, B):
    return mul(F, A, inv(F, B))


def pw(F, A, n, L):
    if n < 0:
        return pw(F, inv(F, A), -n, L)
    R = const(F, 1, L)
    base = A
    while n:
        if n & 1:
            R = mul(F, R, base)
        n >>= 1
        if n:
            base = mul(F, base, base)
    return R


def deriv(F, A):
    p = F.p
    return LS(A.v - 1, [F.m(F.k(A.v + i), x) for i, x in enumerate(A.c)])


# ---------------- the curve ----------------
class Curve:
    def __init__(self, F, a, b):
        self.F, self.a, self.b = F, F.k(a), F.k(b)
        assert (4 * a ** 3 + 27 * b ** 2) % F.p

    def rhs(self, x):
        F = self.F
        return F.a(F.a(F.m(F.m(x, x), x), F.m(self.a, x)), self.b)

    def points(self):
        F = self.F
        pts = []
        for x in range(F.q):
            for y in F.SQRT.get(self.rhs(x), []):
                pts.append((x, y))
        return pts

    def add_pts(self, P, Q):
        F = self.F
        if P is None:
            return Q
        if Q is None:
            return P
        (x1, y1), (x2, y2) = P, Q
        if x1 == x2:
            if F.a(y1, y2) == 0:
                return None
            lam = F.m(F.a(F.m(F.k(3), F.m(x1, x1)), self.a), F.i(F.m(F.k(2), y1)))
        else:
            lam = F.m(F.s(y2, y1), F.i(F.s(x2, x1)))
        x3 = F.s(F.s(F.m(lam, lam), x1), x2)
        y3 = F.s(F.m(lam, F.s(x1, x3)), y1)
        return (x3, y3)

    def negp(self, P):
        return None if P is None else (P[0], self.F.NEG[P[1]])

    def local(self, Q, L):
        """(X, Y) Laurent series of x, y in a uniformizer at Q (None = O); also returns dX/ds."""
        F = self.F
        two_inv = F.i(F.k(2))
        if Q is None:
            # w = z^3 + a z w^2 + b w^3 (power series in z), x = z / w, y = -1 / w
            z = LS(1, [1] + [0] * (L + 6))
            w = pw(F, z, 3, L + 8)
            for _ in range(L + 4):
                w = add(F, add(F, pw(F, z, 3, L + 8), scal(F, self.a, mul(F, z, mul(F, w, w)))),
                        scal(F, self.b, pw(F, w, 3, L + 8)))
                w = LS(w.v, w.c[:L + 6])
            X = div(F, z, w)
            Y = neg(F, inv(F, w))
            X, Y = LS(X.v, X.c[:L]), LS(Y.v, Y.c[:L])
            return X, Y
        x0, y0 = Q
        if y0 != 0:
            X = LS(0, [x0, 1] + [0] * (L - 2))
            Y = const(F, y0, L)
            R = self.rhs_series(X, L)
            for _ in range(10):
                Y = scal(F, two_inv, add(F, Y, div(F, R, Y)))
                Y = LS(Y.v, Y.c[:L])
            return X, Y
        Y = LS(0, [0, 1] + [0] * (L - 2))
        X = const(F, x0, L)
        for _ in range(10):
            G = sub(F, self.rhs_series(X, L), mul(F, Y, Y))
            dG = add(F, scal(F, F.k(3), mul(F, X, X)), const(F, self.a, L))
            X = sub(F, X, div(F, G, dG))
            X = LS(X.v, X.c[:L])
        return X, Y

    def rhs_series(self, X, L):
        F = self.F
        return add(F, add(F, mul(F, mul(F, X, X), X), scal(F, self.a, X)), const(F, self.b, L))


def xi2_series(F, X, Y, t, yR, L):
    lam = div(F, sub(F, neg(F, Y), const(F, yR, L)), sub(F, X, const(F, t, L)))
    return sub(F, sub(F, mul(F, lam, lam), X), const(F, t, L))


def divisor(E, func, L):
    """{place: order} over all rational places; func(X, Y) -> Laurent series."""
    out = {}
    for Q in [None] + E.points():
        X, Y = E.local(Q, L)
        o = val(E.F, func(X, Y))
        if o:
            out[Q] = o
    return out


def predicted(E, R, P0):
    pred1, pred2 = {}, {}
    for Q, o in ((P0, 1), (E.negp(P0), 1), (None, -2)):
        pred1[Q] = pred1.get(Q, 0) + o
    RmP0 = E.add_pts(R, E.negp(P0))
    RpP0 = E.add_pts(R, P0)
    for Q, o in ((RmP0, 1), (RpP0, 1), (R, -2)):
        pred2[Q] = pred2.get(Q, 0) + o
    return {k: v for k, v in pred1.items() if v}, {k: v for k, v in pred2.items() if v}


def fmt(F, Q):
    if Q is None:
        return "O"
    return "(%s,%s)" % tuple(("%d" % (z % F.p)) if z < F.p else ("%d+%dw" % (z % F.p, z // F.p)) for z in Q)


def run_lemma32(out_path):
    configs = []
    # find parameters over small primes for each configuration, two curves each
    for p in (13, 17, 19, 23):
        found = {}
        for a in range(p):
            for b in range(p):
                if (4 * a ** 3 + 27 * b ** 2) % p == 0:
                    continue
                for t in range(p):
                    cub = (t ** 3 + a * t + b) % p
                    kinds = []
                    if t and b and cub and a and (4 * b * t - a * a) % p:
                        kinds.append("generic")
                    if t and b == 0 and cub:
                        kinds.append("b=0")
                    if t and cub == 0 and b:
                        kinds.append("R 2-torsion")
                    if t and b == 0 and cub == 0:
                        kinds.append("R 2-torsion and b=0")
                    if t and b and a and (4 * b * t - a * a) % p == 0:
                        kinds.append("R=+-2P0")
                    if t == 0 and b:
                        kinds.append("t=0")
                    if t == 0 and b == 0:
                        kinds.append("t=0,b=0")
                    for kd in kinds:
                        found.setdefault(kd, [])
                        if len(found[kd]) < 1 and (a, b) not in [(c[1], c[2]) for c in configs if c[3] == kd]:
                            found[kd].append((p, a, b, t, kd))
        for kd, lst in found.items():
            for c in lst:
                if sum(1 for d in configs if d[4] == kd) < 2:
                    configs.append(c)
    report = []
    L = 14
    for (p, a, b, t, kind) in configs:
        F = Fq(p)
        E = Curve(F, a, b)
        y0 = F.SQRT[F.k(b)][0]
        P0 = (0, y0)
        yR = F.SQRT[E.rhs(F.k(t))][0]
        R = (F.k(t), yR)
        d1 = divisor(E, lambda X, Y: X, L)
        d2 = divisor(E, lambda X, Y: xi2_series(F, X, Y, F.k(t), yR, L), L)
        pr1, pr2 = predicted(E, R, P0)
        # orders of monomials at O and R
        mono = {}
        for (ea, eb) in ((1, 0), (0, 1), (3, 2), (5, 7)):
            res = {}
            for name, Q in (("O", None), ("R", R)):
                X, Y = E.local(Q, L + 20)
                f = mul(F, pw(F, X, ea, L + 20), pw(F, xi2_series(F, X, Y, F.k(t), yR, L + 20), eb, L + 20))
                res[name] = val(F, f)
            mono["%d,%d" % (ea, eb)] = res
        row = {"p": p, "a": a, "b": b, "t": t, "kind": kind,
               "div_xi1": {fmt(F, Q): o for Q, o in d1.items()},
               "div_xi2": {fmt(F, Q): o for Q, o in d2.items()},
               "div_xi1_matches": d1 == pr1, "div_xi2_matches": d2 == pr2,
               "deg0": sum(d1.values()) == 0 and sum(d2.values()) == 0,
               "zeros_total": [sum(o for o in d1.values() if o > 0), sum(o for o in d2.values() if o > 0)],
               "ord_O_xi2": d2.get(None, 0), "ord_R_xi1": d1.get(R, 0),
               "monomial_orders_at_O_R": mono,
               "lemma32_conclusion": all(v["O"] == -2 * int(k.split(",")[0]) and v["R"] == -2 * int(k.split(",")[1])
                                         for k, v in mono.items()),
               "S_size": len({None, R, P0, E.negp(P0), E.add_pts(R, P0), E.add_pts(R, E.negp(P0))})}
        report.append(row)
        print(json.dumps({k: row[k] for k in ("p", "a", "b", "t", "kind", "div_xi1", "div_xi2", "div_xi1_matches",
                                               "div_xi2_matches", "ord_O_xi2", "ord_R_xi1", "lemma32_conclusion",
                                               "S_size")}))
    ok = all(r["div_xi1_matches"] and r["div_xi2_matches"] and r["deg0"] and r["zeros_total"] == [2, 2]
             for r in report)
    ok = ok and all(r["lemma32_conclusion"] == (r["t"] != 0) for r in report)
    print("LEMMA 3.2 CHECKS PASS (conclusion holds exactly when t != 0)" if ok else "LEMMA 3.2 CHECK FAILED")
    if out_path:
        with open(out_path, "w") as fh:
            json.dump(report, fh, indent=1)
    return 0 if ok else 1


# ---------------- addendum (b): non-constant F^p coefficients ----------------
def det_valuation(F, M):
    """valuation of det M for a square matrix of Laurent series (DVR elimination)."""
    n = len(M)
    M = [[norm(F, x) for x in row] for row in M]
    total = 0
    for col in range(n):
        best, bv = None, None
        for r in range(col, n):
            x = M[r][col]
            if x.c:
                if bv is None or x.v < bv:
                    best, bv = r, x.v
        if best is None:
            raise ArithmeticError("singular to precision")
        M[col], M[best] = M[best], M[col]
        total += bv
        piv_inv = inv(F, M[col][col])
        for r in range(col + 1, n):
            if M[r][col].c:
                f = mul(F, M[r][col], piv_inv)
                M[r] = [norm(F, sub(F, M[r][j], mul(F, f, M[col][j]))) if j > col else M[r][j] for j in range(n)]
    return total


def rank_mod(F, rows):
    rows = [list(r) for r in rows]
    rk, ncol = 0, len(rows[0])
    for c in range(ncol):
        piv = next((r for r in range(rk, len(rows)) if rows[r][c]), None)
        if piv is None:
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        iv = F.i(rows[rk][c])
        for r in range(len(rows)):
            if r != rk and rows[r][c]:
                f = F.m(rows[r][c], iv)
                rows[r] = [F.s(x, F.m(f, y)) for x, y in zip(rows[r], rows[rk])]
        rk += 1
    return rk


def run_twisted(p, a, b, t, out_path):
    F = Fq(p)
    E = Curve(F, a, b)
    assert t % p and (t ** 3 + a * t + b) % p, "need t != 0 and R not 2-torsion"
    rng = random.Random(7)
    k = F.k
    tt = k(t)
    yR = F.SQRT[E.rhs(tt)][0]
    R = (tt, yR)
    y0 = F.SQRT[k(b)][0]
    P0 = (0, y0)
    Sset = {None, R, P0, E.negp(P0), E.add_pts(R, P0), E.add_pts(R, E.negp(P0))}
    xS = {Q[0] for Q in Sset if Q is not None}
    # S_3 coefficients c_e, e in {0,1,2}^2
    A, B, T = a % p, b % p, t % p
    coeff = {(2, 0): T * T, (0, 2): T * T, (1, 1): -2 * T * T - 2 * A, (2, 1): -2 * T, (1, 2): -2 * T,
             (1, 0): -2 * A * T - 4 * B, (0, 1): -2 * A * T - 4 * B, (0, 0): A * A - 4 * B * T, (2, 2): 1}
    E9 = sorted(coeff)
    cE = {e: k(coeff[e]) for e in E9}
    twists = [(0, 0), (3, 0), (0, 5), (7, 2), (1, 1), (0, 9), (6, 6), (2, 8), (9, 0)]
    u = {e: tw for e, tw in zip(E9, twists)}
    m = {e: (e[0] + p * u[e][0], e[1] + p * u[e][1]) for e in E9}
    # h = (x1 - alpha) / (x2 - beta): zeros and poles outside S
    alpha = next(c for c in range(1, p) if k(c) not in xS)
    xi2_at_S = {tt, 0}
    for Q in (P0, E.negp(P0)):
        RQ = E.add_pts(R, E.negp(Q))
        if RQ is not None:
            xi2_at_S.add(RQ[0])
    beta = next(c for c in range(1, p) if k(c) not in xi2_at_S and c != alpha)
    Lbig = 40

    def funcs(X, Y, L):
        x2 = xi2_series(F, X, Y, tt, yR, L)
        h = div(F, sub(F, X, const(F, k(alpha), L)), sub(F, x2, const(F, k(beta), L)))
        hp = pw(F, h, p, L)
        phi = {e: scal(F, cE[e], mul(F, hp, mul(F, pw(F, X, -p * u[e][0], L), pw(F, x2, -p * u[e][1], L))))
               for e in E9}
        mon = {e: mul(F, pw(F, X, m[e][0], L), pw(F, x2, m[e][1], L)) for e in E9}
        v = {e: mul(F, phi[e], mon[e]) for e in E9}
        return x2, h, phi, mon, v

    places = [None] + E.points()
    rep = {"p": p, "a": a, "b": b, "t": t, "alpha": alpha, "beta": beta,
           "exponents_m_e": {str(e): m[e] for e in E9}, "S_size": len(Sset)}
    ord_phi_bad, ordv = [], {}
    rel_ok = True
    hzp = {}
    for Q in places:
        X, Y = E.local(Q, Lbig)
        x2, h, phi, mon, v = funcs(X, Y, Lbig)
        for e in E9:
            o = val(F, phi[e])
            if o % p:
                ord_phi_bad.append((fmt(F, Q), str(e), o))
        oh = val(F, h)
        if oh and Q not in Sset:
            hzp[fmt(F, Q)] = oh
        ordv[Q] = {e: val(F, v[e]) for e in E9}
        s = v[E9[0]]
        for e in E9[1:]:
            s = add(F, s, v[e])
        # exact arithmetic: every known coefficient of sum_e v_e must vanish, and the known
        # precision must exceed the smallest valuation by a margin
        if norm(F, s).c or s.prec() < min(ordv[Q].values()) + 12:
            rel_ok = False
    rep["relation_vanishes_at_every_place"] = rel_ok
    rep["ord_phi_not_divisible_by_p"] = ord_phi_bad
    rep["zeros_poles_of_h_outside_S"] = hzp
    rep["deg0_each_v_e"] = all(sum(ordv[Q][e] for Q in places) == 0 for e in E9)
    H = -sum(min(ordv[Q].values()) for Q in places)
    n0 = 28
    rep["H_v"] = H
    rep["n0_times_S"] = n0 * len(Sset)
    rep["H_le_bound"] = H <= n0 * len(Sset)
    pair_ok, cong_ok = True, True
    for e in E9:
        for f in E9:
            dO = ordv[None][e] - ordv[None][f]
            dR = ordv[R][e] - ordv[R][f]
            if abs(dO) + abs(dR) > 2 * H:
                pair_ok = False
            if (dO + 2 * (m[e][0] - m[f][0])) % p or (dR + 2 * (m[e][1] - m[f][1])) % p:
                cong_ok = False
    rep["pairwise_orders_le_2H"] = pair_ok
    rep["orders_at_O_R_congruent_to_minus_2_delta_mod_p"] = cong_ok

    def dist(n):
        r = n % p
        return min(r, p - r)
    spreads_p = [dist(2 * (m[e][0] - m[f][0])) + dist(2 * (m[e][1] - m[f][1])) for e in E9 for f in E9]
    spreads_int = [abs(2 * (m[e][0] - m[f][0])) + abs(2 * (m[e][1] - m[f][1])) for e in E9 for f in E9]
    rep["thm33_mod_p_max_spread"] = max(spreads_p)
    rep["integer_max_spread"] = max(spreads_int)
    rep["thm33_bound_6_8_7"] = 6 * 8 * 7

    # F^p-dependence of {xi^(m_e)} and minimality; rank over the constants
    gen = [Q for Q in E.points() if Q[1] != 0 and Q not in Sset and Q[0] != k(alpha)]
    rng.shuffle(gen)
    test_pts = gen[:10]
    nine_zero, eight_nonzero = True, {str(e): False for e in E9}
    evals = []
    for Q in test_pts:
        X, Y = E.local(Q, 12)
        x2 = xi2_series(F, X, Y, tt, yR, 12)
        mon = {e: norm(F, mul(F, pw(F, X, m[e][0], 12), pw(F, x2, m[e][1], 12))) for e in E9}
        assert all(mon[e].v == 0 for e in E9)
        evals.append([mon[e].c[0] for e in E9])
        Mx = [[mon[e].c[i] for e in E9] for i in range(9)]
        if rank_mod(F, Mx) == 9:
            nine_zero = False
        for x in E9:
            sub8 = [ee for ee in E9 if ee != x]
            if rank_mod(F, [[mon[ee].c[i] for ee in sub8] for i in range(8)]) == 8:
                eight_nonzero[str(x)] = True
    more = gen[10:40]
    for Q in more:
        X, Y = E.local(Q, 4)
        x2 = xi2_series(F, X, Y, tt, yR, 4)
        evals.append([norm(F, mul(F, pw(F, X, m[e][0], 4), pw(F, x2, m[e][1], 4))).c[0] for e in E9])
    rep["Fp_dependent_nine_Wronskian_zero_at_10_points"] = nine_zero
    rep["every_eight_independent"] = all(eight_nonzero.values())
    rep["rank_over_constants_at_%d_points" % len(evals)] = rank_mod(F, evals)

    # local Wronskian inequality of step (d): ord_P(W) >= -n0 ord_P(dx) + (A_P - mu_P) - n0 [P in S]
    check_places = list(Sset) + [Q for Q in E.points() if fmt(F, Q) in hzp] + \
        [Q for Q in E.points() if Q[1] == 0][:2] + test_pts[:2]
    local_rows = []
    for Q in check_places:
        Lw = 70
        X, Y = E.local(Q, Lw)
        x2, h, phi, mon, v = funcs(X, Y, Lw)
        ords = {e: val(F, v[e]) for e in E9}
        mu = min(ords.values())
        iP = next(e for e in E9 if ords[e] == mu)
        cols = [v[e] for e in E9 if e != iP]
        dX = deriv(F, X)
        dX_inv = inv(F, dX)
        rows = []
        cur = cols
        for i in range(8):
            rows.append(cur)
            cur = [mul(F, deriv(F, c), dX_inv) for c in cur]
        Wv = det_valuation(F, rows)
        AP = sum(ords.values())
        lower = -n0 * val(F, dX) + (AP - mu) - n0 * (1 if Q in Sset else 0)
        local_rows.append({"place": fmt(F, Q), "in_S": Q in Sset, "ord_W": Wv, "lower_bound": lower,
                           "ord_dx": val(F, dX), "holds": Wv >= lower})
    rep["local_wronskian_inequality"] = local_rows
    ok = (rel_ok and not ord_phi_bad and rep["deg0_each_v_e"] and rep["H_le_bound"] and pair_ok and cong_ok
          and nine_zero and rep["every_eight_independent"] and rep["rank_over_constants_at_%d_points" % len(evals)] == 9
          and all(r["holds"] for r in local_rows) and hzp)
    print(json.dumps(rep, indent=1, default=str))
    print("TWISTED-RELATION CHECKS PASS" if ok else "TWISTED-RELATION CHECK FAILED")
    if out_path:
        with open(out_path, "w") as fh:
            json.dump(rep, fh, indent=1, default=str)
    return 0 if ok else 1


if __name__ == "__main__":
    if sys.argv[1] == "lemma32":
        sys.exit(run_lemma32(sys.argv[2] if len(sys.argv) > 2 else None))
    if sys.argv[1] == "twisted":
        p, a, b, t = (int(x) for x in sys.argv[2:6])
        sys.exit(run_twisted(p, a, b, t, sys.argv[6] if len(sys.argv) > 6 else None))
