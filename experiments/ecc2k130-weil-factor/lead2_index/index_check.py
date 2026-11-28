#!/usr/bin/env python3
"""Independent (no PARI) cross-check of [O_F (x) Z[zeta_ell] : Z[pi,pibar]] for primes ell <= LMAX.

Pure-Python big-integer arithmetic. Element of O_F[z]/Phi_ell stored as two length-(ell-1) integer
coefficient vectors (c0 + c1*w), w^2 = w - 2.  alpha = w - 1, beta = -w, pi = alpha z, pibar = beta z^{-1}.

Method A (exact): integer row-echelon (Hermite) basis of the Z-span of the generators
  pi^i (0<=i<2d), pibar^j (1<=j<2d); index = |product of pivots|.
Method B (local structure): for each prime p | U_ell, F_p-rank of the generators (expected ell), and
  exact membership tests f*e_k in R for every standard basis vector e_k (exponent of O/R divides f), and
  zeta in R.  Order f^(ell-2) + exponent | f + p-rank ell  ==>  O/R = (Z/f)^(ell-2).
U_ell is factored by trial division (independent of PARI).
"""
import sys

def U(n):
    a, b = 0, 1
    if n == 0: return 0
    for _ in range(n - 1): a, b = b, -b - 2 * a
    return b

def factor(n):
    n = abs(n); out = {}; p = 2
    while p * p <= n:
        while n % p == 0: out[p] = out.get(p, 0) + 1; n //= p
        p += 1
    if n > 1: out[n] = out.get(n, 0) + 1
    return out

def isprime(n): return n > 1 and factor(n) == {n: 1}

def polmul(x, y, ell):
    """x, y: lists of length ell (coeffs of z^0..z^{ell-1} in Z[z]/(z^ell-1))."""
    r = [0] * ell
    for i, xi in enumerate(x):
        if xi:
            for j, yj in enumerate(y):
                if yj: r[(i + j) % ell] += xi * yj
    return r

def mul(a, b, ell):
    a0, a1 = a; b0, b1 = b
    t0 = polmul(a0, b0, ell); t2 = polmul(a1, b1, ell)
    t1 = [u + v for u, v in zip(polmul(a0, b1, ell), polmul(a1, b0, ell))]
    return ([u - 2 * v for u, v in zip(t0, t2)], [u + v for u, v in zip(t1, t2)])

def coords(a, ell):
    out = []
    for c in a:
        top = c[ell - 1]
        out += [c[i] - top for i in range(ell - 1)]   # z^{ell-1} = -(1+...+z^{ell-2})
    return out

def gens(ell):
    d = ell - 1
    zero = [0] * ell
    one = ([1] + [0] * (ell - 1), zero[:])
    z = ([0, 1] + [0] * (ell - 2), zero[:])
    zinv = ([0] * (ell - 1) + [1], zero[:])
    alpha = ([-1] + [0] * (ell - 1), [1] + [0] * (ell - 1))
    beta = (zero[:], [-1] + [0] * (ell - 1))
    pi, pb = mul(alpha, z, ell), mul(beta, zinv, ell)
    G, x = [], one
    for _ in range(2 * d): G.append(coords(x, ell)); x = mul(x, pi, ell)
    x = pb
    for _ in range(1, 2 * d): G.append(coords(x, ell)); x = mul(x, pb, ell)
    return G, coords(z, ell)

def egcd(a, b):
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b:
        q, r = divmod(a, b); a, b = b, r
        x0, x1 = x1, x0 - q * x1; y0, y1 = y1, y0 - q * y1
    return a, x0, y0

def bareiss_det(M):
    """Exact determinant of a square integer matrix (fraction-free Bareiss)."""
    A = [r[:] for r in M]; n = len(A); sign = 1; prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if A[i][k]), None)
            if sw is None: return 0
            A[k], A[sw] = A[sw], A[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]

def hnf(G, n):
    """Upper-triangular basis {col: row} of the Z-span L of the rows of G.
    D := |det(first n rows)| (exact Bareiss) is a nonzero multiple of [Z^n : L], so D*e_j in L for all j;
    we seed the pivots with D*e_j and reduce every non-pivot entry mod D (adds elements of D*Z^n, so the
    lattice is unchanged). This is the classical modular HNF; it does not use the conjectured index."""
    D = abs(bareiss_det(G[:n]))
    assert D != 0, "first n generators not independent"
    piv = {j: [D if k == j else 0 for k in range(n)] for j in range(n)}
    for v in G:
        v = [x % D for x in v]
        for c in range(n):
            if v[c] == 0: continue
            r = piv[c]
            g, s_, t_ = egcd(r[c], v[c])
            a, b = r[c] // g, v[c] // g
            newr = [s_ * x + t_ * y for x, y in zip(r, v)]
            v = [(a * y - b * x) % D for x, y in zip(r, v)]
            newr = [newr[k] if k == c else newr[k] % D for k in range(n)]
            newr[c] = abs(newr[c])
            piv[c] = newr
    cols = sorted(piv)
    for i in reversed(range(len(cols))):
        ci = cols[i]
        for j in range(i):
            cj = cols[j]; r = piv[cj]; q = r[ci] // piv[ci][ci]
            if q: piv[cj] = [x - q * y for x, y in zip(r, piv[ci])]
    return piv

def member(piv, v, n):
    v = v[:]
    for c in range(n):
        if v[c] == 0: continue
        if c not in piv: return False
        q, r = divmod(v[c], piv[c][c])
        if r: return False
        v = [x - q * y for x, y in zip(v, piv[c])]
    return True

def rank_mod(G, p):
    rows = [[x % p for x in r] for r in G]; rk = 0; ncol = len(rows[0])
    for col in range(ncol):
        pv = next((i for i in range(rk, len(rows)) if rows[i][col]), None)
        if pv is None: continue
        rows[rk], rows[pv] = rows[pv], rows[rk]
        inv = pow(rows[rk][col], p - 2, p)
        rows[rk] = [x * inv % p for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][col]:
                f = rows[i][col]; rows[i] = [(x - f * y) % p for x, y in zip(rows[i], rows[rk])]
        rk += 1
    return rk

def main(lmax):
    print("ell|U_ell|factor(|U|)|index(pyHNF)|index==|U|^(ell-2)|p-ranks|f*e_k in R all k|zeta in R|structure (Z/f)^(ell-2) proved-from-data")
    allok = True
    for ell in range(3, lmax + 1):
        if not isprime(ell): continue
        n = 2 * (ell - 1)
        G, zc = gens(ell)
        piv = hnf(G, n)
        assert len(piv) == n, "rank deficient"
        idx = 1
        for c in piv: idx *= piv[c][c]
        u = U(ell); f = abs(u); fac = factor(f)
        ok1 = idx == f ** (ell - 2)
        ranks = {p: rank_mod(G, p) for p in fac}
        e = all(member(piv, [f if k == j else 0 for k in range(n)], n) for j in range(n))
        zin = member(piv, zc, n)
        struct = ok1 and e and all(r == ell for r in ranks.values())
        ok = ok1 and (zin == (f == 1)) and (struct or f == 1)
        allok &= ok
        print(f"{ell}|{u}|{fac if fac else 1}|{idx}|{ok1}|{ranks}|{e}|{zin}|{struct or f == 1}")
        sys.stdout.flush()
    print("ALL OK" if allok else "MISMATCH")

if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 60)
