#!/usr/bin/env python3
"""Independent check of the universal specialization identities behind the
fixed-m growth argument (docs/proofs/N8-higher-summand-growth.md, PR #45).

Written from scratch for this check: it does not import the PR #45 resultant
engine. Pure Python, modular arithmetic, formal-degree Sylvester resultants.

S_2(x1,x2) = x1 - x2,
S_3 = (x1-x2)^2 x3^2 - 2((x1+x2)(x1 x2 + a) + 2b) x3 + (x1 x2 - a)^2 - 4b(x1+x2),
S_{n+1}(x_1..x_{n+1}) = Res_X(S_n(x_1..x_{n-1},X), S_3(x_n,x_{n+1},X)),
formal degrees 2^(n-2) and 2 in X (the paper's definition, Section 5).

S_j(x_1..x_{j-1}, w) is computed as its coefficient vector in w (formal degree
2^(j-2)) by evaluating the formal-degree determinant at 2^(j-2)+1 nodes and
interpolating. Specializing coefficients commutes with a determinant of fixed
formal size, so values are exact for every specialization, including those
where a leading coefficient vanishes.

Checks (each at random points over several primes, Schwartz-Zippel style,
plus targeted degenerate points):
  SYM  S_n is symmetric (n = 4..7), a check that the recursion is the right one.
  A    S_n(U,1,1) = c^(2^(n-3)) S_{n-2}(U)^4 on E_a: y^2 = x^3 + a x - 1 - a, c=(a+3)^2.
  B    S_4(x,y,1,t) = (t-1)^4 S_3(x,y,u)^2, u = (t+a+2)/(t-1).
  EVEN m=2r+2:  S_{m+1}(x,y,1^(2r),t) = c^(r 4^r) S_3(x,y,t)^(4^r).
  ODD  m=2r+3:  S_{m+1}(x,y,1^(2r+1),t) = c^(2r 4^r) (t-1)^(4^(r+1)) S_3(x,y,u)^(2 4^r).
  A_mutant_detected (control) c^(2^(n-3)+1) in place of c^(2^(n-3)) must fail.
  TRAP (positive control) a specialize-first resultant with actual degrees
       must DISAGREE with A at n=4 (it loses the factor (x-y)^4).
  UNIT the existence step: on E_0 over a small prime with p > 4((2N)^m + 7),
       a rational target outside the m-fold sum set of points with abscissa in
       H has S_{m+1}(h, x(R)) != 0 on all of H^m, and a target inside it has a zero.
"""
import random, sys, json, itertools
from functools import lru_cache

def inv(x, p): return pow(x % p, p - 2, p)

def det(M, p):
    M = [r[:] for r in M]; n = len(M); d = 1
    for c in range(n):
        piv = next((r for r in range(c, n) if M[r][c] % p), None)
        if piv is None: return 0
        if piv != c: M[c], M[piv] = M[piv], M[c]; d = -d
        d = d * M[c][c] % p; ic = inv(M[c][c], p)
        for r in range(c + 1, n):
            f = M[r][c] * ic % p
            if f:
                for j in range(c, n): M[r][j] = (M[r][j] - f * M[c][j]) % p
    return d % p

def sylvester_res(f, g, p):
    """Res of f (coeff list low->high, formal degree len(f)-1) and g likewise."""
    d, e = len(f) - 1, len(g) - 1; n = d + e
    fh, gh = f[::-1], g[::-1]
    M = []
    for i in range(e): M.append([0] * i + fh + [0] * (n - d - 1 - i))
    for i in range(d): M.append([0] * i + gh + [0] * (n - e - 1 - i))
    return det(M, p)

def s3_coeffs_in_x3(x1, x2, a, b, p):
    return [((x1 * x2 - a) ** 2 - 4 * b * (x1 + x2)) % p,
            (-2 * ((x1 + x2) * (x1 * x2 + a) + 2 * b)) % p,
            ((x1 - x2) ** 2) % p]

def s3_coeffs_in_x2(x1, x3, a, b, p):  # S_3 symmetric: same formula with roles swapped
    return s3_coeffs_in_x3(x1, x3, a, b, p)

def interpolate(nodes, vals, p):
    n = len(nodes); coeffs = [0] * n
    for i in range(n):
        num = [1]; den = 1
        for j in range(n):
            if j == i: continue
            num = [(c1 - nodes[j] * c0) % p for c0, c1 in zip(num + [0], [0] + num)]
            den = den * (nodes[i] - nodes[j]) % p
        s = vals[i] * inv(den, p) % p
        for t in range(n): coeffs[t] = (coeffs[t] + s * num[t]) % p
    return coeffs

def poly_last(prefix, a, b, p):
    """Coefficients in w of S_j(prefix, w), j = len(prefix)+1, formal degree 2^(j-2)."""
    j = len(prefix) + 1
    if j == 2: return [prefix[0] % p, p - 1]
    if j == 3: return s3_coeffs_in_x3(prefix[0], prefix[1], a, b, p)
    f = poly_last(prefix[:-1], a, b, p)          # S_{j-1}(x_1..x_{j-2}, z), deg 2^(j-3)
    D = 2 ** (j - 2); nodes = list(range(D + 1)); vals = []
    for w in nodes:
        g = s3_coeffs_in_x3(prefix[-1], w, a, b, p)  # S_3(x_{j-1}, w, z) in z
        vals.append(sylvester_res(f, g, p))
    return interpolate(nodes, vals, p)

def S(xs, a, b, p):
    c = poly_last(list(xs[:-1]), a, b, p); x = xs[-1] % p
    return sum(ci * pow(x, i, p) for i, ci in enumerate(c)) % p

def S_naive(xs, a, b, p):
    """Specialize first, then resultant with ACTUAL degrees (the invalid shortcut)."""
    def strip(c):
        c = c[:]
        while len(c) > 1 and c[-1] % p == 0: c.pop()
        return c
    n = len(xs)
    if n <= 3: return S(xs, a, b, p)
    f = strip(poly_last(list(xs[:n - 2]), a, b, p))
    g = strip(s3_coeffs_in_x3(xs[n - 2], xs[n - 1], a, b, p))
    if len(f) == 1 and len(g) == 1: return 1
    return sylvester_res(f, g, p)

PRIMES = [257, 65537, 2**31 - 1, 2**61 - 1]
results = {}
def record(name, ok, n):
    r = results.setdefault(name, [0, 0]); r[0] += ok; r[1] += n

def main(trials=6, mmax=8, seed=20261123):
    rng = random.Random(seed)
    for p in PRIMES:
        for _ in range(trials):
            a = rng.randrange(p); b = (-1 - a) % p
            while (a + 3) % p == 0 or (4 * a + 3) % p == 0: a = rng.randrange(p); b = (-1 - a) % p
            c = (a + 3) ** 2 % p
            # SYM on a general curve (b free)
            bb = rng.randrange(p)
            for n in range(4, 8):
                xs = [rng.randrange(p) for _ in range(n)]
                v = S(xs, a, bb, p); ys = xs[:]; rng.shuffle(ys)
                record("SYM", int(S(ys, a, bb, p) == v), 1)
            # A
            for n in range(4, mmax + 2):
                U = [rng.randrange(p) for _ in range(n - 2)]
                for UU in (U, [U[0]] * (n - 2), [1] + U[1:]):  # generic, diagonal, one 1
                    lhs = S(UU + [1, 1], a, b, p)
                    rhs = pow(c, 2 ** (n - 3), p) * pow(S(UU, a, b, p), 4, p) % p
                    record("A", int(lhs == rhs), 1)
                    # mutation control: a wrong scalar exponent must be detected
                    mut = pow(c, 2 ** (n - 3) + 1, p) * pow(S(UU, a, b, p), 4, p) % p
                    if rhs != 0: record("A_mutant_detected", int(lhs != mut), 1)
            # B
            for _k in range(3):
                x, y = rng.randrange(p), rng.randrange(p); t = rng.randrange(2, p)
                if (t - 1) % p == 0: continue
                u = (t + a + 2) * inv(t - 1, p) % p
                lhs = S([x, y, 1, t], a, b, p)
                rhs = pow(t - 1, 4, p) * pow(S([x, y, u], a, b, p), 2, p) % p
                record("B", int(lhs == rhs), 1)
            # parity formulas, m = 2..mmax, with generic and x=y points
            for m in range(2, mmax + 1):
                for (x, y) in ((rng.randrange(p), rng.randrange(p)), (7, 7), (1, 5)):
                    t = rng.randrange(2, p)
                    if (t - 1) % p == 0: continue
                    xs = [x, y] + [1] * (m - 2) + [t]
                    lhs = S(xs, a, b, p)
                    if m % 2 == 0:
                        r = (m - 2) // 2
                        rhs = pow(c, r * 4 ** r, p) * pow(S([x, y, t], a, b, p), 4 ** r, p) % p
                        record("EVEN", int(lhs == rhs), 1)
                    else:
                        r = (m - 3) // 2; u = (t + a + 2) * inv(t - 1, p) % p
                        rhs = (pow(c, 2 * r * 4 ** r, p) * pow(t - 1, 4 ** (r + 1), p)
                               * pow(S([x, y, u], a, b, p), 2 * 4 ** r, p)) % p
                        record("ODD", int(lhs == rhs), 1)
            # TRAP positive control: specialize-first must lose (x-y)^4 at n=4
            x, y = rng.randrange(p), rng.randrange(p)
            if x != y:
                formal = S([x, y, 1, 1], a, b, p); naive = S_naive([x, y, 1, 1], a, b, p)
                record("TRAP_detects", int(formal != naive and naive == pow(c, 2, p)), 1)
    unit_check()
    ok = all(v[0] == v[1] for v in results.values())
    print(json.dumps({k: {"passed": v[0], "total": v[1]} for k, v in results.items()}, indent=1))
    print("ALL PASS" if ok else "FAILURE"); return 0 if ok else 1

def ec_add(P, Q, a, p):
    if P is None: return Q
    if Q is None: return P
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2 and (y1 + y2) % p == 0: return None
    lam = ((3 * x1 * x1 + a) * inv(2 * y1, p) if P == Q else (y2 - y1) * inv(x2 - x1, p)) % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)

def unit_check():
    """E_0: y^2 = x^3 - 1, m = 3, k = 2 (N = 4). Needs p = 1 mod 4 and p > 4((2N)^3 + 7) = 2076.
    Points with abscissa in H may lie over F_{p^2}; we use p = 1 mod 12 ... and simply test
    the two directions of the defining property on rational data, plus the full H^3 scan."""
    m, k = 3, 2; N = 2 ** k; a, b = 0, -1
    p = next(q for q in range(4 * ((2 * N) ** m + 7) + 1, 10 ** 6)
             if all(q % d for d in range(2, int(q ** .5) + 1)) and q % N == 1)
    g = next(h for h in range(2, p) if pow(h, (p - 1) // 2, p) != 1)
    H = sorted({pow(g, (p - 1) // N * i, p) for i in range(N)})
    pts = []
    for x in range(p):
        r = (x ** 3 + a * x + b) % p
        if r == 0: pts.append((x, 0)); continue
        if pow(r, (p - 1) // 2, p) == 1:
            y = next(y for y in range(1, p) if y * y % p == r)  # small p, fine
            pts += [(x, y), (x, p - y)]
    # rational points with abscissa in H (the algebraic-closure ones over F_{p^2}
    # are excluded only through S != 0 on H^3, which is what we test)
    hp = [P for P in pts if P[0] in H]
    sums = set()
    for P1 in hp:
        for P2 in hp:
            for P3 in hp:
                s = ec_add(ec_add(P1, P2, a, p), P3, a, p)
                if s is not None: sums.add(s[0])
    bad = set(sums) | {0, 1, (-2) % p}
    outside = [P for P in pts if P[0] not in bad]
    inside = [P for P in pts if P[0] in sums]
    def all_nonzero(t):
        return all(S([h1, h2, h3, t], a, b, p) != 0 for h1 in H for h2 in H for h3 in H)
    # an outside target: S must be nonzero on H^3 unless a zero comes from F_{p^2}-points;
    # we scan several outside targets and require at least one unit instance (the
    # existence claim), and we report how many outside rational-sum targets are units.
    n_unit = sum(all_nonzero(P[0]) for P in outside[:40])
    record("UNIT_exists", int(n_unit >= 1), 1)
    # inside targets must have a zero on H^3 (defining property, forward direction)
    n_zero = sum(not all_nonzero(P[0]) for P in inside[:20])
    record("UNIT_inside_has_zero", int(n_zero == min(20, len(inside))), 1)
    print(f"UNIT: p={p}, |H-points rational|={len(hp)}, outside sampled={min(40,len(outside))},"
          f" units among them={n_unit}; inside sampled={min(20,len(inside))}, with zero={n_zero}",
          file=sys.stderr)

if __name__ == "__main__":
    sys.exit(main(*(int(x) for x in sys.argv[1:])))
