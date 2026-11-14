"""Literal, budget-enforced execution of the proof of Theorem 1 (LFD <= min(k1,k2) + 2).

For an instance (p, a, b, t, k1, k2) with k1 <= k2 every row is built only with the three
operations of the quotient closure V_D, D = k1 + 2:

  gen(m)        tau_m g for a monomial m with wt(m) + 2 <= D,
  lin(...)      F_p-linear combinations,
  mulvar(r, v)  product with a chain variable, allowed only if wt(r) <= budget - 1.

The class Row refuses any other construction.  The script then checks, for EVERY monomial m of
the group algebra, that the row built by the induction of the proof equals m - pi(m), where
pi = s o q is computed independently from the solution set (brute-force roots of the fibre
polynomials) and is verified to be the interpolation projector on Z.  This is the statement
m - pi(m) in V_{max(D, wt m)}, which implies I ∩ M_<=D' ⊆ V_D' for every D' >= D.

Usage: python3 selector_bound_check.py [--out results-selector.jsonl]
"""
from __future__ import annotations

import argparse
import json
import random
import sys
import time

import numpy as np

from chain_algebra import s3_coeffs, s3_value, subgroup, wt


class Ctx:
    def __init__(self, p, a, b, t, k1, k2):
        self.p, self.a, self.b, self.t, self.k1, self.k2 = p, a, b, t, k1, k2
        self.N1, self.N2 = 1 << k1, 1 << k2
        self.D = k1 + 2
        self.W = np.array([[wt(e1) + wt(e2) for e2 in range(self.N2)] for e1 in range(self.N1)])
        g = np.zeros((self.N1, self.N2), dtype=np.int64)
        for (e1, e2), c in s3_coeffs(p, a, b, t).items():
            g[e1 % self.N1, e2 % self.N2] = (g[e1 % self.N1, e2 % self.N2] + c) % p
        self.g = g
        self.gdeg = int(self.W[g != 0].max()) if g.any() else 0
        self.mulvar_count = 0

    def weight(self, arr):
        nz = arr != 0
        return int(self.W[nz].max()) if nz.any() else -1


class Row:
    """An element of A together with the budget at which it is certified to lie in V_budget."""
    __slots__ = ('arr', 'budget')

    def __init__(self, arr, budget, _token=None):
        assert _token is Row._TOKEN, 'rows are built only by gen / lin / mulvar'
        self.arr, self.budget = arr, budget

    _TOKEN = object()

    @staticmethod
    def gen(ctx, m, budget):
        e1, e2 = m
        assert wt(e1) + wt(e2) + max(ctx.gdeg, 2) <= budget, 'generator multiple over budget'
        return Row(np.roll(np.roll(ctx.g, e1, axis=0), e2, axis=1), budget, Row._TOKEN)

    @staticmethod
    def zero(ctx, budget):
        return Row(np.zeros((ctx.N1, ctx.N2), dtype=np.int64), budget, Row._TOKEN)

    @staticmethod
    def lin(ctx, terms, budget):
        """sum c_i r_i; every r_i must be certified at a budget <= the target budget (V is monotone)."""
        acc = np.zeros((ctx.N1, ctx.N2), dtype=np.int64)
        for c, r in terms:
            assert r.budget <= budget
            acc = (acc + (c % ctx.p) * r.arr) % ctx.p
        return Row(acc, budget, Row._TOKEN)

    @staticmethod
    def lin_stack(ctx, combined, rows, budget):
        """A linear combination of the certified rows `rows`, precomputed as one matrix product
        (`combined` = coefficients @ stacked rows); only the budgets are checked here."""
        assert all(r.budget <= budget for r in rows)
        return Row(combined.reshape(ctx.N1, ctx.N2) % ctx.p, budget, Row._TOKEN)

    def mulvar(self, ctx, v, budget):
        assert self.budget <= budget
        assert ctx.weight(self.arr) <= budget - 1, 'variable product over budget'
        ctx.mulvar_count += 1
        i, j = v
        return Row(np.roll(self.arr, 1 << j, axis=i), budget, Row._TOKEN)


def run_instance(p, a, b, t, k1, k2, check_all_degrees=True):
    assert k1 <= k2
    assert (4 * a ** 3 + 27 * b * b) % p, 'singular curve'
    ctx = Ctx(p, a, b, t, k1, k2)
    D, N1, N2 = ctx.D, ctx.N1, ctx.N2
    inv = lambda x: pow(int(x) % p, p - 2, p)
    _, H1 = subgroup(p, k1)
    _, H2 = subgroup(p, k2)
    X2 = (1, 0)                      # the variable x_2 = u_{2,0}
    co = s3_coeffs(p, a, b, t)
    info = dict(p=p, a=a, b=b, t=t, k1=k1, k2=k2, D=D, fibres={0: 0, 1: 0, 2: 0}, cases={})

    def e_alpha(al):
        arr = np.zeros((N1, N2), dtype=np.int64)
        ia, iN = inv(al), inv(N1)
        for e in range(N1):
            arr[e, 0] = pow(ia, e, p) * iN % p
        return arr

    def elem(poly2, al):
        """(sum_e c_e x_2^e) * e_alpha as an array (for comparisons only)."""
        arr = np.zeros((N1, N2), dtype=np.int64)
        ea = e_alpha(al)[:, 0]
        for e, c in poly2.items():
            arr[:, e % N2] = (arr[:, e % N2] + c * ea) % p
        return arr

    # ------------------------------------------------------------------ per-fibre derivations
    basis = []            # (alpha index, r) with r in {0: 1, 1: x_2}
    corr = {}             # (alpha index, r, j) -> Row for (u_{2,j} * r - pi_alpha(u_{2,j} r)) e_alpha
    pi_coef = {}          # alpha index -> function e2 -> coefficients of pi_alpha(x_2^e2) on R_alpha
    empty_rows = {}       # alpha index -> Row for e_alpha (empty fibre)
    for ai, al in enumerate(H1):
        ea = e_alpha(al)
        assert ctx.weight(ea) <= k1
        # step 2: e_alpha * g, a combination of tau_m g with wt(m) <= k1
        R = Row.lin(ctx, [(int(ea[e, 0]), Row.gen(ctx, (e, 0), D)) for e in range(N1) if ea[e, 0]], D)
        q = [0, 0, 0]
        for (e1, e2), c in co.items():
            q[e2] = (q[e2] + c * pow(al, e1, p)) % p
        assert any(q), 'identically vanishing fibre polynomial'
        assert np.array_equal(R.arr, elem({0: q[0], 1: q[1], 2: q[2]}, al)), 'e_alpha g != q_alpha e_alpha'
        assert ctx.weight(R.arr) <= k1 + 1
        roots = [E for E, x in enumerate(H2) if (q[0] + q[1] * x + q[2] * x * x) % p == 0]
        assert [E for E, x in enumerate(H2) if s3_value(p, a, b, t, al, x) == 0] == roots
        info['fibres'][len(roots)] += 1
        deg = 2 if q[2] else (1 if q[1] else 0)
        if deg == 0:
            case = 'const'
            empty_rows[ai] = Row.lin(ctx, [(inv(q[0]), R)], D)
            assert not roots
        elif deg == 1:
            beta = -q[0] * inv(q[1]) % p
            w = Row.lin(ctx, [(inv(q[1]), R)], D)          # (x - beta) e_alpha
            ws = [w]
            bj = beta
            for j in range(k2):
                nxt = Row.lin(ctx, [(1, ws[j].mulvar(ctx, (1, j), D)), (bj, ws[j])], D)
                bj = bj * bj % p
                ws.append(nxt)                             # (rho(u_j^2) - beta^(2^(j+1))) e_alpha
            T = ws[k2]                                     # (1 - beta^N) e_alpha
            if (1 - bj) % p:
                case = 'linear-empty'
                empty_rows[ai] = Row.lin(ctx, [(inv(1 - bj), T)], D)
                assert not roots
            else:
                case = 'linear-point'
                assert len(roots) == 1 and H2[roots[0]] == beta
                basis.append((ai, 0))
                for j in range(k2):
                    corr[(ai, 0, j)] = ws[j]
                pi_coef[ai] = (lambda beta: (lambda e: [pow(beta, e, p)]))(beta)
        else:
            s, n = -q[1] * inv(q[2]) % p, q[0] * inv(q[2]) % p     # x^2 = s x - n modulo q
            Q1 = Row.lin(ctx, [(inv(q[2]), R)], D)                 # (rho(x^2) - s x + n) e_alpha

            def mulx(ab):                                          # x (a + b x) mod q
                return (-ab[1] * n % p, (ab[0] + ab[1] * s) % p)

            ab = (0, 1)                                            # x^(2^j) mod q
            w = Row.zero(ctx, D)                                   # (u_j - a_j - b_j x) e_alpha
            ws, wxs, abs_ = [], [], []
            for j in range(k2 + 1):
                wx = Row.lin(ctx, [(1, w.mulvar(ctx, X2, D)), (ab[1], Q1)], D)   # (x u_j - x(a+bx) mod q) e_alpha
                ws.append(w); wxs.append(wx); abs_.append(ab)
                if j == k2:
                    break
                w = Row.lin(ctx, [(1, w.mulvar(ctx, (1, j), D)), (ab[0], w), (ab[1], wx)], D)
                xab = mulx(ab)
                ab = ((ab[0] * ab[0] + ab[1] * xab[0]) % p, (ab[0] * ab[1] + ab[1] * xab[1]) % p)
            # ws[k2] = (1 - A - B x) e_alpha, wxs[k2] = (x - x(A + B x) mod q) e_alpha
            A, B = abs_[k2]
            v1 = ((1 - A) % p, (-B) % p)
            xAB = mulx((A, B))
            v2 = ((-xAB[0]) % p, (1 - xAB[1]) % p)
            assert np.array_equal(ws[k2].arr, elem({0: v1[0], 1: v1[1]}, al))
            assert np.array_equal(wxs[k2].arr, elem({0: v2[0], 1: v2[1]}, al))
            det = (v1[0] * v2[1] - v1[1] * v2[0]) % p
            if det:
                case = 'quadratic-empty'
                assert not roots
                di = inv(det)
                empty_rows[ai] = Row.lin(ctx, [(v2[1] * di, ws[k2]), (-v1[1] * di, wxs[k2])], D)
            elif any(v1) or any(v2):
                case = 'quadratic-point'
                v, src = (v1, ws[k2]) if any(v1) else (v2, wxs[k2])
                assert v[1], 'gcd of degree one must be monic linear'
                L = Row.lin(ctx, [(inv(v[1]), src)], D)            # (x - beta) e_alpha
                beta = -v[0] * inv(v[1]) % p
                assert len(roots) == 1 and H2[roots[0]] == beta
                basis.append((ai, 0))
                for j in range(k2):
                    corr[(ai, 0, j)] = Row.lin(ctx, [(1, ws[j]), (abs_[j][1], L)], D)
                pi_coef[ai] = (lambda beta: (lambda e: [pow(beta, e, p)]))(beta)
            else:
                case = 'quadratic-pair'
                assert len(roots) == 2
                basis.append((ai, 0)); basis.append((ai, 1))
                for j in range(k2):
                    corr[(ai, 0, j)] = ws[j]
                    corr[(ai, 1, j)] = wxs[j]

                def red(e, s=s, n=n):
                    r, base = (1, 0), (0, 1)
                    while e:
                        if e & 1:
                            xb = (-base[1] * n % p, (base[0] + base[1] * s) % p)
                            r = ((r[0] * base[0] + r[1] * xb[0]) % p, (r[0] * base[1] + r[1] * xb[1]) % p)
                        xb = (-base[1] * n % p, (base[0] + base[1] * s) % p)
                        base = ((base[0] * base[0] + base[1] * xb[0]) % p, (base[0] * base[1] + base[1] * xb[1]) % p)
                        e >>= 1
                    return [r[0], r[1]]
                pi_coef[ai] = red
        info['cases'][case] = info['cases'].get(case, 0) + 1
        if ai in empty_rows:
            assert np.array_equal(empty_rows[ai].arr, ea)

    # ------------------------------------------------------------------ the projector pi = s o q
    Z = [(E1, E2) for E1, x1 in enumerate(H1) for E2, x2 in enumerate(H2) if s3_value(p, a, b, t, x1, x2) == 0]
    assert len(basis) == len(Z)
    bidx = {bk: i for i, bk in enumerate(basis)}
    barr = np.array([elem({r: 1}, H1[ai]).ravel() for ai, r in basis], dtype=np.int64).reshape((len(basis), N1 * N2))
    assert all(ctx.weight(x.reshape(N1, N2)) <= k1 + 1 for x in barr)

    def evaluate(arr, pt):
        E1, E2 = pt
        v1 = np.array([H1[(E1 * e) % N1] for e in range(N1)], dtype=np.int64)
        v2 = np.array([H2[(E2 * e) % N2] for e in range(N2)], dtype=np.int64)
        return int((v1 @ (arr % p) % p) @ v2 % p)

    def pi_vec(m):
        """Coefficients of pi(x1^e1 x2^e2) on the basis {r e_alpha}."""
        e1, e2 = m
        c = np.zeros(len(basis), dtype=np.int64)
        for ai, f in pi_coef.items():
            sc = pow(H1[ai], e1, p)
            for r, val in enumerate(f(e2)):
                c[bidx[(ai, r)]] = sc * val % p
        return c

    # pi(m) agrees with m on Z (sampled monomials x all points), so pi is the interpolation projector
    rng = random.Random(1)
    for _ in range(min(40, N1 * N2)):
        m = (rng.randrange(N1), rng.randrange(N2))
        arr = (pi_vec(m) @ barr % p).reshape(N1, N2) if len(basis) else np.zeros((N1, N2), dtype=np.int64)
        for pt in Z:
            assert evaluate(arr, pt) == H1[(pt[0] * m[0]) % N1] * H2[(pt[1] * m[1]) % N2] % p

    # correction rows C(v, b) = v b - pi(v b) for block-2 variables (block-1 corrections are 0)
    if basis:
        C2 = [np.array([corr[(ai, r, j)].arr.ravel() for ai, r in basis], dtype=np.int64) for j in range(k2)]
        for j in range(k2):
            for (ai, r), row in zip(basis, C2[j]):
                assert corr[(ai, r, j)].budget == D
                expect = np.roll(barr[bidx[(ai, r)]].reshape(N1, N2), 1 << j, axis=1).ravel() % p
                # v b - pi(v b): pi(v b) = sum over the same fibre
                f = pi_coef[ai]
                vb = f((1 << j) + r)            # pi_alpha(x^(2^j + r))
                sub = np.zeros(N1 * N2, dtype=np.int64)
                for rr, val in enumerate(vb):
                    sub = (sub + val * barr[bidx[(ai, rr)]]) % p
                assert np.array_equal(row % p, (expect - sub) % p), 'correction row mismatch'

    # ------------------------------------------------------------------ the induction (star)
    # rows[m] = m - pi(m), certified at budget max(D, wt m)
    mons = sorted(((e1, e2) for e1 in range(N1) for e2 in range(N2)), key=lambda m: (wt(m[0]) + wt(m[1]), m))
    rows = {}
    one = np.zeros((N1, N2), dtype=np.int64); one[0, 0] = 1
    rows[(0, 0)] = Row.lin(ctx, [(1, r) for r in empty_rows.values()], D)
    maxdeg = (k1 + k2) if check_all_degrees else D
    checked = 0
    for m in mons:
        r = wt(m[0]) + wt(m[1])
        if r > maxdeg:
            break
        budget = max(D, r)
        if r == 0:
            row = rows[m]
        else:
            if m[0]:
                j = (m[0] & -m[0]).bit_length() - 1
                v, mp = (0, j), (m[0] ^ (1 << j), m[1])
            else:
                j = (m[1] & -m[1]).bit_length() - 1
                v, mp = (1, j), (m[0], m[1] ^ (1 << j))
            prev = rows[mp]
            assert ctx.weight(prev.arr) <= max(r - 1, k1 + 1)
            terms = [(1, prev.mulvar(ctx, v, budget))]
            if v[0] == 1 and basis:
                c = pi_vec(mp)
                acc = (c @ C2[v[1]]) % p
                # a combination of certified correction rows (budget D <= budget)
                terms.append((1, Row.lin_stack(ctx, acc, [corr[(ai, rr, v[1])] for ai, rr in basis], D)))
            row = Row.lin(ctx, terms, budget)
            rows[m] = row
        expect = np.zeros((N1, N2), dtype=np.int64); expect[m] = 1
        if basis:
            expect = (expect - (pi_vec(m) @ barr % p).reshape(N1, N2)) % p
        assert np.array_equal(row.arr % p, expect), ('induction row mismatch', m)
        assert ctx.weight(row.arr) <= max(r, k1 + 1)
        checked += 1
    info.update(points=len(Z), monomials_checked=checked, max_degree_checked=maxdeg,
                variable_products=ctx.mulvar_count, ok=True)
    return info


def instances():
    out = []
    rng = random.Random(20261110)
    # exhaustive small planes: every nonsingular (a, b) and every t
    for p, k in [(5, 1), (5, 2), (7, 1), (13, 2), (17, 2), (17, 3), (17, 4)]:
        for a in range(p):
            for b in range(p):
                if (4 * a ** 3 + 27 * b * b) % p == 0:
                    continue
                for t in range(p):
                    out.append((p, a, b, t, k, k))
    # random instances, deeper
    for p, k, count in [(41, 3, 150), (97, 4, 60), (97, 5, 40), (193, 5, 30), (193, 6, 6), (257, 6, 4),
                        (12289, 5, 12), (40961, 4, 12)]:
        got = 0
        while got < count:
            a, b, t = rng.randrange(p), rng.randrange(p), rng.randrange(p)
            if (4 * a ** 3 + 27 * b * b) % p == 0:
                continue
            if got % 3 == 0:                 # force t into the subgroup (a fibre polynomial drops degree)
                t = subgroup(p, k)[1][rng.randrange(1 << k)]
            out.append((p, a, b, t, k, k))
            got += 1
    # unequal depths k1 < k2
    for p, k1, k2, count in [(17, 1, 4, 40), (17, 2, 4, 40), (97, 2, 5, 20), (97, 3, 5, 20), (193, 3, 6, 10), (193, 4, 6, 6)]:
        got = 0
        while got < count:
            a, b, t = rng.randrange(p), rng.randrange(p), rng.randrange(p)
            if (4 * a ** 3 + 27 * b * b) % p == 0:
                continue
            out.append((p, a, b, t, k1, k2))
            got += 1
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='results-selector.jsonl')
    ap.add_argument('--limit', type=int, default=0)
    args = ap.parse_args()
    todo = instances()
    if args.limit:
        todo = todo[:args.limit]
    t0 = time.time()
    summary = {}
    with open(args.out, 'w') as fh:
        for inst in todo:
            info = run_instance(*inst)
            key = f"p={info['p']},k1={info['k1']},k2={info['k2']}"
            s = summary.setdefault(key, dict(instances=0, unit=0, proper=0, cases={}, fibres={0: 0, 1: 0, 2: 0}))
            s['instances'] += 1
            s['unit' if info['points'] == 0 else 'proper'] += 1
            for c, v in info['cases'].items():
                s['cases'][c] = s['cases'].get(c, 0) + v
            for c, v in info['fibres'].items():
                s['fibres'][c] += v
            fh.write(json.dumps(info, sort_keys=True) + '\n')
    json.dump(dict(total=len(todo), seconds=round(time.time() - t0, 1), all_ok=True, by_family=summary),
              open(args.out.replace('.jsonl', '-summary.json'), 'w'), indent=1, sort_keys=True)
    print(json.dumps(dict(total=len(todo), seconds=round(time.time() - t0, 1), all_ok=True)))
    for k, v in summary.items():
        print(k, v)


if __name__ == '__main__':
    sys.exit(main())
