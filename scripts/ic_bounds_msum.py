#!/usr/bin/env python3
"""Chain degree of NF(S_{m+1}) for m = 2..8 summands, and the last fall degree it forces.

For the direct presentation G_{m,k} u {g}, g = NF(S_{m+1}(x_1, ..., x_m, t)), over F_p with 2^k | p - 1:

  * upper bound  LFD <= mk + 1          (kernel-checked: ChainQuotient.reduced_inputs_lfd_bound);
  * lower bound  LFD >= deg g + 1       when h = NF(x_1 g) has degree <= deg g and is not a multiple of g
                                        (kernel-checked general statement: CarryLower.chain_carry_lower).

So LFD = mk + 1 exactly whenever deg g = mk, and LFD >= deg g + 1 in general.  This script computes
g exactly, without expanding S_{m+1}: the chain quotient is the algebra of functions on H^m (H the
subgroup of order N = 2^k), so g is the inverse discrete Fourier transform of the values of S_{m+1} on
the grid H^m.  The values come from Semaev's recursion
    S_{n}(x_1, ..., x_{n-1}, X) = Res_Y( S_{n-1}(x_1, ..., x_{n-2}, Y), S_3(x_{n-1}, X, Y) ),
carried out with numbers: the univariate polynomial S_n(x_1..x_{n-1}, X) is interpolated from its values
at 2^{n-2} + 1 nodes, and a resultant with a quadratic is a norm from F_p[Y]/(quadratic).

  python3 scripts/ic_bounds_msum.py selftest
  python3 scripts/ic_bounds_msum.py run [--max-log2-points 21] [--out experiments/ic-bounds-msum]
  python3 scripts/ic_bounds_msum.py one --m 5 --k 2 --prime 97 --a 1 --b 0 --t 3

numpy only.  Exact: all arithmetic is modulo p <= 65537 in int64 (products stay below 2^63).
"""
import argparse
import hashlib
import json
import random
import sys
import time
from itertools import product
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


# ----------------------------------------------------------------------------------------------------
# arithmetic
# ----------------------------------------------------------------------------------------------------
def s3_abc(a, b, p, u, v):
    """S_3(u, v, Y) = A Y^2 + B Y + C for y^2 = x^3 + a x + b (u, v numbers or int64 arrays)."""
    A = (u - v) % p
    A = A * A % p
    uv = u * v % p
    B = (-2 * (((u + v) % p) * ((uv + a) % p) % p + 2 * b)) % p
    w = (uv - a) % p
    C = (w * w - 4 * b * ((u + v) % p)) % p
    return A, B, C


def vpow(x, e, p):
    """x ** e mod p for an int64 array x and a non-negative integer e."""
    r = np.ones_like(x)
    base = x % p
    while e:
        if e & 1:
            r = r * base % p
        base = base * base % p
        e >>= 1
    return r


def res_quad_batch(P, A, B, C, p):
    """Res_Y(P_T(Y), A_K Y^2 + B_K Y + C_K) with formal degrees (d, 2), for every row T of P and every K.

    P: (T, d + 1) coefficients, low degree first.  A, B, C: (K,).  Returns (T, K).
    A != 0:  A^d * P(r1) P(r2), computed as a norm in F_p[Y]/(Y^2 - sY - q), s = -B/A, q = -C/A.
    A == 0, B != 0:  p_d * B^d * P(-C/B).     A == B == 0:  p_d^2 * C^d.   (All are the Sylvester value.)
    """
    T, d1 = P.shape
    d = d1 - 1
    K = A.shape[0]
    out = np.zeros((T, K), dtype=np.int64)
    inv = lambda x: pow(int(x), p - 2, p)
    mA = A != 0
    if mA.any():
        iA = np.array([inv(x) for x in A[mA]], dtype=np.int64)
        s = (-B[mA] * iA) % p
        q = (-C[mA] * iA) % p
        al = np.zeros((T, s.shape[0]), dtype=np.int64)
        be = np.zeros_like(al)
        for j in range(d, -1, -1):          # acc = acc * Y + c_j in the quotient ring
            al, be = (al * s + be) % p, (al * q + P[:, j:j + 1]) % p
        val = (be * be + (al * be % p) * s - (al * al % p) * q) % p
        out[:, mA] = val * vpow(A[mA], d, p) % p
    mB = (~mA) & (B != 0)
    if mB.any():
        r = (-C[mB] * np.array([inv(x) for x in B[mB]], dtype=np.int64)) % p
        acc = np.zeros((T, r.shape[0]), dtype=np.int64)
        for j in range(d, -1, -1):
            acc = (acc * r + P[:, j:j + 1]) % p
        out[:, mB] = (P[:, d:d + 1] * vpow(B[mB], d, p) % p) * acc % p
    mC = (~mA) & (B == 0)
    if mC.any():
        out[:, mC] = (P[:, d:d + 1] * P[:, d:d + 1] % p) * vpow(C[mC], d, p) % p
    return out


def inverse_vandermonde(n, p):
    """Inverse modulo p of V[r, j] = r^j for the nodes r = 0..n-1 (coefficients = V^-1 . values)."""
    M = [[pow(r, j, p) for j in range(n)] + [int(i == r) for i in range(n)] for r in range(n)]
    for c in range(n):
        piv = next(r for r in range(c, n) if M[r][c])
        M[c], M[piv] = M[piv], M[c]
        iv = pow(M[c][c], p - 2, p)
        M[c] = [x * iv % p for x in M[c]]
        for r in range(n):
            if r != c and M[r][c]:
                f = M[r][c]
                M[r] = [(x - f * y) % p for x, y in zip(M[r], M[c])]
    return np.array([row[n:] for row in M], dtype=np.int64)


def subgroup(p, k):
    N = 1 << k
    assert (p - 1) % N == 0, 'need 2^k | p - 1'
    z = next(z for z in range(2, p) if pow(z, (p - 1) // 2, p) == p - 1)
    g0 = pow(z, (p - 1) // N, p)
    H = [pow(g0, e, p) for e in range(N)]
    assert len(set(H)) == N
    return g0, H


def grid_values(p, a, b, t, m, k, chunk=4096):
    """S_{m+1}(x_1, ..., x_m, t) on H^m, as an int64 array of shape (N,) * m (axis i <-> exponent of x_i)."""
    g0, H = subgroup(p, k)
    N = 1 << k
    Hn = np.array(H, dtype=np.int64)
    if m == 1:
        raise ValueError('m >= 2')
    # P[j]: coefficients in X of S_{j+1}(x_1, ..., x_j, X), shape (N^j, 2^(j-1) + 1)
    u = np.repeat(Hn, N)
    v = np.tile(Hn, N)
    A, B, C = s3_abc(a, b, p, u, v)
    P = np.stack([C, B, A], axis=1)                      # j = 2
    for j in range(2, m - 1 if m > 2 else 2):            # build j + 1 while j + 1 <= m - 1
        d = 1 << j                                        # degree in X of S_{j+2}(x_1..x_{j+1}, X)
        assert p > d + 1, 'prime too small for the interpolation nodes'
        nodes = np.arange(d + 1, dtype=np.int64)
        xs = np.repeat(Hn, d + 1)
        Xs = np.tile(nodes, N)
        qa, qb, qc = s3_abc(a, b, p, xs, Xs)             # S_3(x_{j+1}, X_r, Y)
        Vinv = inverse_vandermonde(d + 1, p)
        rows = []
        for lo in range(0, P.shape[0], chunk):
            vals = res_quad_batch(P[lo:lo + chunk], qa, qb, qc, p)        # (T, N * (d + 1))
            vals = vals.reshape(-1, d + 1)                                 # (T * N, d + 1)
            rows.append(vals @ Vinv.T % p)
        P = np.concatenate(rows, axis=0)
    if m == 2:
        tt = int(t) % p
        vals = (P[:, 2] * tt % p * tt + P[:, 1] * tt + P[:, 0]) % p
        return vals.reshape((N,) * 2)
    qa, qb, qc = s3_abc(a, b, p, Hn, np.full(N, int(t) % p, dtype=np.int64))   # S_3(x_m, t, Y)
    out = []
    for lo in range(0, P.shape[0], max(1, chunk * 8)):
        out.append(res_quad_batch(P[lo:lo + chunk * 8], qa, qb, qc, p))
    return np.concatenate(out, axis=0).reshape((N,) * m)


def normal_form(vals, p, k):
    """Coefficients of the chain normal form: c[f] = N^-m sum_e vals[e] g0^(-e.f) (inverse DFT on each axis)."""
    g0, _ = subgroup(p, k)
    N = 1 << k
    gi = pow(g0, p - 2, p)
    iN = pow(N, p - 2, p)
    D = np.array([[pow(gi, e * f, p) * iN % p for e in range(N)] for f in range(N)], dtype=np.int64)
    c = vals
    for ax in range(vals.ndim):
        c = np.moveaxis(np.tensordot(D, c, axes=([1], [ax])) % p, 0, ax)
    return c


def weight_grid(m, k):
    N = 1 << k
    pc = np.array([bin(e).count('1') for e in range(N)], dtype=np.int64)
    w = np.zeros((N,) * m, dtype=np.int64)
    for ax in range(m):
        shape = [1] * m
        shape[ax] = N
        w = w + pc.reshape(shape)
    return w


def analyse(p, a, b, t, m, k):
    t0 = time.time()
    vals = grid_values(p, a, b, t, m, k)
    N = 1 << k
    c = normal_form(vals, p, k)
    w = weight_grid(m, k)
    nz = c != 0
    assert nz.any(), 'S_{m+1} vanishes on the whole grid'
    deg = int(w[nz].max())
    top = int(c[(N - 1,) * m])
    h = np.roll(c, 1, axis=0)                           # NF(x_1 g): exponent of x_1 shifted by one, modulo N
    deg_h = int(w[h != 0].max())
    idx = tuple(int(i[0]) for i in np.nonzero(nz))
    ratio = int(h[idx]) * pow(int(c[idx]), p - 2, p) % p
    proportional = bool(np.array_equal(h, c * ratio % p))
    carry = deg_h <= deg and not proportional
    lower = deg + 1 if carry else None
    upper = m * k + 1
    hist = {int(d): int(((w == d) & nz).sum()) for d in sorted(set(w[nz].tolist()))}
    rec = dict(p=p, a=a, b=b, t=t, m=m, k=k, summation_polynomial=m + 1, variables=m * k, points=N ** m,
               zeros=int((vals == 0).sum()), nf_terms=int(nz.sum()), nf_degree=deg, nf_degree_bound=m * min(k, m - 1),
               nf_terms_by_degree=hist, top_coefficient=top, shift_degree=deg_h, shift_proportional=proportional,
               carry_hypotheses=bool(carry), lfd_lower=lower, lfd_upper=upper,
               lfd_exact=(upper if lower == upper else None),
               nf_sha256=hashlib.sha256(np.ascontiguousarray(c).tobytes()).hexdigest(), seconds=round(time.time() - t0, 2))
    return rec, vals, c


# ----------------------------------------------------------------------------------------------------
# instances
# ----------------------------------------------------------------------------------------------------
def rational_abscissa(p, a, b, rng):
    while True:
        x = rng.randrange(1, p)
        r = (x * x * x + a * x + b) % p
        if r and pow(r, (p - 1) // 2, p) == 1:
            return x


def instances(k, count, seed=20261112):
    rng = random.Random(seed * 131 + k)
    primes = [q for q in (97, 193, 257, 7681, 12289, 65537) if (q - 1) % (1 << k) == 0]
    out = []
    for p in primes:
        curves = [(1, 0), (0, 1)]
        while len(curves) < 4:
            a, b = rng.randrange(1, p), rng.randrange(1, p)
            if (4 * a ** 3 + 27 * b * b) % p:
                curves.append((a, b))
        for a, b in curves:
            for _ in range(2):
                out.append((p, a, b, rational_abscissa(p, a, b, rng)))
    rng.shuffle(out)
    return out[:count]


GRID = [(m, k) for m in range(2, 9) for k in range(1, 7)]


def cmd_run(a):
    out = ROOT / a.out
    out.mkdir(parents=True, exist_ok=True)
    path = out / 'results.jsonl'
    done = set()
    if path.exists():
        for line in path.read_text().splitlines():
            r = json.loads(line)
            done.add((r['m'], r['k'], r['p'], r['a'], r['b'], r['t']))
    with path.open('a') as f:
        for m, k in GRID:
            lg = m * k
            if lg > a.max_log2_points or (k > m + 1 and m > 2):
                continue
            count = 48 if lg <= 16 else 8 if lg <= 21 else 2
            for p, ca, cb, t in instances(k, count):
                if p <= (1 << (m - 2)) + 1 or (m, k, p, ca, cb, t) in done:
                    continue
                rec, _, _ = analyse(p, ca, cb, t, m, k)
                f.write(json.dumps(rec, sort_keys=True) + '\n')
                f.flush()
                print(json.dumps({x: rec[x] for x in ('m', 'k', 'p', 'a', 'b', 't', 'nf_degree', 'lfd_lower', 'lfd_upper', 'seconds')}), flush=True)
    summarize(out)


def summarize(out):
    rows = [json.loads(l) for l in (out / 'results.jsonl').read_text().splitlines()]
    cells = {}
    for r in rows:
        c = cells.setdefault((r['m'], r['k']), dict(m=r['m'], k=r['k'], summation_polynomial=r['m'] + 1, variables=r['m'] * r['k'],
                                                   nf_degree_bound=r['nf_degree_bound'], lfd_upper=r['lfd_upper'], instances=0,
                                                   nf_degree_hist={}, full_degree=0, lfd_exact=0, lfd_lower_hist={}, carry_failed=0))
        c['instances'] += 1
        c['nf_degree_hist'][str(r['nf_degree'])] = c['nf_degree_hist'].get(str(r['nf_degree']), 0) + 1
        c['full_degree'] += int(r['nf_degree'] == r['nf_degree_bound'])
        c['lfd_exact'] += int(r['lfd_exact'] is not None)
        c['carry_failed'] += int(not r['carry_hypotheses'])
        if r['lfd_lower'] is not None:
            c['lfd_lower_hist'][str(r['lfd_lower'])] = c['lfd_lower_hist'].get(str(r['lfd_lower']), 0) + 1
    table = [cells[key] for key in sorted(cells)]
    (out / 'summary.json').write_text(json.dumps(dict(cells=table, records=len(rows),
                                                      script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()), indent=1) + '\n')
    lines = ['| m | S | k | variables | instances | deg NF (bound) | full degree | LFD lower (count) | LFD upper | exact |', '|---|---|---|---|---|---|---|---|---|---|']
    for c in table:
        lines.append('| %d | S_%d | %d | %d | %d | %s (%d) | %d | %s | %d | %d |' % (
            c['m'], c['summation_polynomial'], c['k'], c['variables'], c['instances'],
            ', '.join('%s: %d' % kv for kv in sorted(c['nf_degree_hist'].items(), key=lambda kv: int(kv[0]))), c['nf_degree_bound'],
            c['full_degree'], ', '.join('%s: %d' % kv for kv in sorted(c['lfd_lower_hist'].items(), key=lambda kv: int(kv[0]))),
            c['lfd_upper'], c['lfd_exact']))
    (out / 'table.md').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


# ----------------------------------------------------------------------------------------------------
# self-test
# ----------------------------------------------------------------------------------------------------
def curve_add(p, a, P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    if P[0] == Q[0] and (P[1] + Q[1]) % p == 0:
        return None
    lam = ((3 * P[0] * P[0] + a) * pow(2 * P[1], p - 2, p) if P == Q else (Q[1] - P[1]) * pow(Q[0] - P[0], p - 2, p)) % p
    x = (lam * lam - P[0] - Q[0]) % p
    return x, (lam * (P[0] - x) - P[1]) % p


def scalar_value(p, a, b, xs, t):
    """S_{m+1}(xs, t) for arbitrary numbers, by the same recursion with plain integers (slow reference)."""
    one = lambda x: np.array([x], dtype=np.int64)
    A, B, C = s3_abc(a, b, p, xs[0], xs[1])
    P = np.array([[C, B, A]], dtype=np.int64)
    for j in range(2, len(xs)):
        d = 1 << j                           # degree in X of S_{j+2}(x_1, ..., x_{j+1}, X)
        nodes = np.arange(d + 1, dtype=np.int64)
        qa, qb, qc = s3_abc(a, b, p, np.full(d + 1, xs[j], dtype=np.int64), nodes)
        P = res_quad_batch(P, qa, qb, qc, p) @ inverse_vandermonde(d + 1, p).T % p
    d = P.shape[1] - 1
    return int(sum(int(P[0, j]) * pow(int(t), j, p) for j in range(d + 1)) % p)


def cmd_selftest(_a):
    ok = True
    rng = random.Random(7)
    sys.path.insert(0, str(ROOT))
    from csd import ic_dreg as I
    from csd import ic_bounds_m3 as M3
    # 1. S_3 and S_4 values against the repository's formulas, on the grid
    for p, a, b, k in [(97, 1, 0, 3), (193, 5, 7, 2), (257, 0, 1, 3)]:
        E = I.PrimeCurve(p, a, b)
        g0, H = subgroup(p, k)
        t = rational_abscissa(p, a, b, rng)
        v2 = grid_values(p, a, b, t, 2, k)
        v3 = grid_values(p, a, b, t, 3, k)
        e2 = all(int(v2[i, j]) == I.s3_prime_value(E, H[i], H[j], t) for i in range(1 << k) for j in range(1 << k))
        e3 = all(int(v3[i, j, l]) == M3.s4_value(E, H[i], H[j], H[l], t) for i, j, l in product(range(1 << k), repeat=3))
        # 2. normal form against the repository's expansion of S_4
        c3 = normal_form(v3, p, k)
        nf = M3.direct_generator(E, t, k)
        N = 1 << k
        mine = {int(i) | (int(j) << k) | (int(l) << (2 * k)): int(c3[i, j, l]) for i, j, l in zip(*np.nonzero(c3))}
        e4 = mine == nf
        print(json.dumps(dict(check='S3,S4 values and NF(S4)', p=p, k=k, s3=e2, s4=e3, nf=e4)))
        ok &= e2 and e3 and e4
    # 3. top coefficient of NF(S_4) at k = 2 is -8 (5 t^3 + 3 a t + 6 b)
    for p, a, b in [(97, 1, 0), (193, 5, 7), (7681, 11, 13)]:
        t = rational_abscissa(p, a, b, rng)
        rec, _, _ = analyse(p, a, b, t, 3, 2)
        e = rec['top_coefficient'] == (-8 * (5 * t ** 3 + 3 * a * t + 6 * b)) % p
        print(json.dumps(dict(check='top coefficient of NF(S_4), k = 2', p=p, ok=e)))
        ok &= e
    # 4. planted sums: S_{m+1}(x(P_1), ..., x(P_m), x(P_1 + ... + P_m)) = 0, and symmetry in the x_i
    for m in range(2, 9):
        p, a, b = 7681, 11, 13
        pts = []
        while len(pts) < m:
            x = rational_abscissa(p, a, b, rng)
            r = (x ** 3 + a * x + b) % p
            y = next(y for y in range(1, p) if y * y % p == r)
            pts.append((x, y))
        S = None
        for P_ in pts:
            S = curve_add(p, a, S, P_)
        if S is None:
            continue
        xs = [P_[0] for P_ in pts]
        z = scalar_value(p, a, b, xs, S[0])
        other = scalar_value(p, a, b, xs, (S[0] + 1) % p)
        perm = xs[:]
        rng.shuffle(perm)
        sym = scalar_value(p, a, b, perm, (S[0] + 1) % p) == other
        print(json.dumps(dict(check='planted sum and symmetry', m=m, zero_at_sum=z == 0, value_off_sum=other, symmetric=sym)))
        ok &= z == 0 and sym
    # 5. grid values against the scalar reference for m = 4, 5
    for m, k in [(4, 2), (5, 2), (4, 3)]:
        p, a, b = 257, 3, 5
        t = rational_abscissa(p, a, b, rng)
        v = grid_values(p, a, b, t, m, k)
        g0, H = subgroup(p, k)
        e = True
        for _ in range(40):
            ix = tuple(rng.randrange(1 << k) for _ in range(m))
            e &= int(v[ix]) == scalar_value(p, a, b, [H[i] for i in ix], t)
        print(json.dumps(dict(check='grid against scalar reference', m=m, k=k, ok=e)))
        ok &= e
    # 6. three summands against the certified records of experiments/ic-bounds-m3 (direct encoding)
    cases = {}
    for line in (ROOT / 'experiments/ic-bounds-m3/cases.jsonl').read_text().splitlines():
        r = json.loads(line)
        cases[r['case_id']] = r
    n = bad = exact_hits = 0
    for line in (ROOT / 'experiments/ic-bounds-m3/results.jsonl').read_text().splitlines():
        r = json.loads(line)
        d = (r.get('results') or {}).get('direct') or {}
        if r.get('encoding') != 'direct' or d.get('exact_lfd') is None:
            continue
        cse = cases[r['case_id']]
        if cse['k'] > 3:
            continue
        rec, _, _ = analyse(cse['prime'], cse['a'], cse['b'], cse['target'][0], 3, cse['k'])
        n += 1
        lo, up, ex = rec['lfd_lower'], rec['lfd_upper'], d['exact_lfd']
        good = ((lo is None or lo <= ex) and ex <= up and (rec['lfd_exact'] is None or rec['lfd_exact'] == ex)
                and [rec['nf_degree']] == r['generator_degrees'] and rec['zeros'] == r['geometry']['scalar_solutions_H3'])
        exact_hits += int(rec['lfd_exact'] is not None)
        bad += int(not good)
    print(json.dumps(dict(check='three-summand records', compared=n, inconsistent=bad, exact_from_degree=exact_hits)))
    ok &= n > 0 and bad == 0
    print('SELFTEST', 'PASSED' if ok else 'FAILED')
    sys.exit(0 if ok else 1)


def cmd_one(a):
    rec, _, _ = analyse(a.prime, a.a, a.b, a.t, a.m, a.k)
    print(json.dumps(rec, indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['selftest', 'run', 'one', 'summarize'])
    ap.add_argument('--max-log2-points', type=int, default=21)
    ap.add_argument('--out', default='experiments/ic-bounds-msum')
    for name in ('m', 'k', 'prime', 'a', 'b', 't'):
        ap.add_argument('--' + name, type=int)
    a = ap.parse_args()
    if a.cmd == 'summarize':
        summarize(ROOT / a.out)
    else:
        {'selftest': cmd_selftest, 'run': cmd_run, 'one': cmd_one}[a.cmd](a)


if __name__ == '__main__':
    main()
