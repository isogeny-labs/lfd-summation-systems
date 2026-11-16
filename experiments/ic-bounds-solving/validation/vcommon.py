"""Validator's own helpers for the solving-degree study (written independently of the lanes).

Everything here is written from the definitions in PROTOCOL.md / SD1 only:
  * prime chain system: R = F_p[u_{i,j}], i < m, j < k, chains u_{i,j}^2 - u_{i,j+1}, u_{i,k-1}^2 - 1;
  * A = F_p[(Z/N)^m], N = 2^k, NF(u_{i,j}) = x_i^{2^j}; weight of x^e = sum of popcounts;
  * the normalized generator g = NF(S_3(x1, x2, t)) (m = 2) or NF(Res_X(S_3(x1,x2,X), S_3(x3,t,X))) (m = 3),
    obtained here as the unique element of A with the right values on H^m (A ~ F_p^{H^m} by the DFT);
  * g^{-1} by pointwise inversion and my own inverse DFT (matmul with 16-bit splitting, exact in int64 for p < 2^31).
No code from laneA/B/C/D or from the repository engines is imported.
"""
import numpy as np

# ------------------------------------------------------------------ number theory
def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def root_of_unity(p, N):
    """an element of exact order N (N a power of two dividing p-1); the smallest base is used."""
    assert (p - 1) % N == 0, (p, N)
    for h in range(2, p):
        w = pow(h, (p - 1) // N, p)
        if N == 1 or pow(w, N // 2, p) == p - 1:
            return w
    raise ValueError


def popcount(x):
    return bin(x).count('1')


# ------------------------------------------------------------------ S_3 and S_4 (numeric, vectorised)
def S3(x1, x2, x3, a, b, p):
    """Semaev S_3 for y^2 = x^3 + a x + b (Semaev 2004):
    (x1-x2)^2 x3^2 - 2((x1+x2)(x1 x2 + a) + 2b) x3 + (x1 x2 - a)^2 - 4b(x1+x2)."""
    x1 = np.asarray(x1, dtype=np.int64) % p
    x2 = np.asarray(x2, dtype=np.int64) % p
    x3 = np.asarray(x3, dtype=np.int64) % p
    d = (x1 - x2) % p
    s = (x1 + x2) % p
    pr = x1 * x2 % p
    A = d * d % p
    B = (-2 * ((s * ((pr + a) % p) % p + 2 * b) % p)) % p
    C = ((pr - a) % p) ** 2 % p
    C = (C - 4 * b % p * s) % p
    return (A * (x3 * x3 % p) % p + B * x3 % p + C) % p


def S3_quadratic_coeffs(x1, x2, a, b, p):
    """coefficients (A, B, C) of S_3(x1, x2, X) as A X^2 + B X + C."""
    x1 = np.asarray(x1, dtype=np.int64) % p
    x2 = np.asarray(x2, dtype=np.int64) % p
    d = (x1 - x2) % p
    s = (x1 + x2) % p
    pr = x1 * x2 % p
    A = d * d % p
    B = (-2 * ((s * ((pr + a) % p) % p + 2 * b) % p)) % p
    C = ((pr - a) % p) ** 2 % p
    C = (C - 4 * b % p * s) % p
    return A, B, C


def S4(x1, x2, x3, t, a, b, p):
    """S_4(x1,x2,x3,t) = Res_X(S_3(x1,x2,X), S_3(x3,t,X)) via the 2x2 resultant
    (A C' - A' C)^2 - (A B' - A' B)(B C' - B' C)  (polynomial identity, used as the definition of ghat for m = 3)."""
    A, B, C = S3_quadratic_coeffs(x1, x2, a, b, p)
    A2, B2, C2 = S3_quadratic_coeffs(x3, np.asarray(t), a, b, p)
    u = (A * C2 - A2 * C) % p
    v = (A * B2 - A2 * B) % p
    w = (B * C2 - B2 * C) % p
    return (u * u % p - v * w % p) % p


# ------------------------------------------------------------------ exact DFT on (Z/N)^m, p < 2^31
def dft_matrix(w, N, p):
    """W[E, e] = w^{E e} mod p as int64."""
    pw = np.empty(N, dtype=np.int64)
    pw[0] = 1
    for i in range(1, N):
        pw[i] = pw[i - 1] * w % p
    idx = (np.arange(N)[:, None] * np.arange(N)[None, :]) % N
    return pw[idx]


def _matmul_mod(W, X, p):
    """exact (W @ X) mod p for 0 <= W < 2^31, 0 <= X < 2^31, N <= 2^16 rows summed: split X into 16-bit halves."""
    N = W.shape[1]
    assert N <= (1 << 16)
    lo = X & 0xFFFF
    hi = X >> 16
    Ylo = (W @ lo) % p          # terms < 2^47, N <= 2^16 terms -> < 2^63
    Yhi = (W @ hi) % p          # terms < 2^46
    return (Ylo + ((Yhi << 16) % p)) % p


def _ntt_axis0(X, w, p):
    """iterative radix-2 Cooley-Tukey along axis 0 of X (shape (N, M)), N a power of two, w of exact order N.
    out[E] = sum_e X[e] w^{E e} mod p.  Exact in int64 for p < 2^31 (every product is reduced at once)."""
    N, M = X.shape
    k = N.bit_length() - 1
    rev = np.zeros(N, dtype=np.int64)
    for i in range(N):
        r = 0
        x = i
        for _ in range(k):
            r = (r << 1) | (x & 1)
            x >>= 1
        rev[i] = r
    Y = X[rev].copy()
    ln = 2
    while ln <= N:
        half = ln // 2
        wl = pow(w, N // ln, p)
        tw = np.empty(half, dtype=np.int64)
        tw[0] = 1
        for j in range(1, half):
            tw[j] = tw[j - 1] * wl % p
        V = Y.reshape(N // ln, ln, M)
        u = V[:, :half, :]
        t = V[:, half:, :] * tw[None, :, None] % p
        V[:, :half, :], V[:, half:, :] = (u + t) % p, (u - t) % p
        ln *= 2
    return Y


def dft_axes(X, w, p, method='auto'):
    """v[E] = sum_e X[e] w^{E.e} over every axis of X (shape (N,)*m); 'matmul' (16-bit split) or 'ntt'."""
    m = X.ndim
    N = X.shape[0]
    use_ntt = (method == 'ntt') or (method == 'auto' and N > 64)
    W = None if use_ntt else dft_matrix(w, N, p)
    Y = X.astype(np.int64) % p
    for ax in range(m):
        Y = np.moveaxis(Y, ax, 0)
        sh = Y.shape
        flat = np.ascontiguousarray(Y.reshape(N, -1))
        Y = (_ntt_axis0(flat, w, p) if use_ntt else _matmul_mod(W, flat, p)).reshape(sh)
        Y = np.moveaxis(Y, 0, ax)
    return Y


def idft_axes(V, w, p):
    N = V.shape[0]
    m = V.ndim
    wi = pow(w, -1, p)
    a = dft_axes(V, wi, p)
    return a * pow(pow(N, m, p), -1, p) % p


def modinv_array(x, p):
    """elementwise inverse by Fermat (all entries nonzero)."""
    r = np.ones_like(x)
    base = x % p
    e = p - 2
    while e:
        if e & 1:
            r = r * base % p
        base = base * base % p
        e >>= 1
    return r


def H_points(p, k):
    N = 1 << k
    w = root_of_unity(p, N)
    H = np.empty(N, dtype=np.int64)
    H[0] = 1
    for i in range(1, N):
        H[i] = H[i - 1] * w % p
    return w, H


def ghat_values(m, k, p, a, b, t):
    """values of the generator on H^m (array (N,)*m): entry [E1,...,Em] = ghat(x_i = w^{E_i})."""
    w, H = H_points(p, k)
    if m == 2:
        V = S3(H[:, None], H[None, :], t, a, b, p)
    elif m == 3:
        V = S4(H[:, None, None], H[None, :, None], H[None, None, :], t, a, b, p)
    else:
        raise ValueError(m)
    return w, V


def weight_array(k, m):
    N = 1 << k
    pc = np.array([popcount(i) for i in range(N)], dtype=np.int64)
    W = pc
    for _ in range(m - 1):
        W = W[..., None] + pc
    return W  # shape (N,)*m


def generator_and_inverse(m, k, p, a, b, t):
    """returns dict with g (coefficient array), delta = maxwt supp g, unit flag, ginv array, maxwt(supp ginv),
    c_top, c_top from the character sum, nzeros of ginv, convolution check."""
    N = 1 << k
    w, V = ghat_values(m, k, p, a, b, t)
    g = idft_axes(V, w, p)
    Wt = weight_array(k, m)
    nzg = g != 0
    delta = int(Wt[nzg].max())
    out = dict(m=m, k=k, p=p, a=a, b=b, t=t, w=int(w), delta=delta, n_terms_g=int(nzg.sum()),
               n_zeros_on_H=int((V == 0).sum()))
    out['unit'] = out['n_zeros_on_H'] == 0
    if not out['unit']:
        return out, g, None
    inv = modinv_array(V, p)
    ginv = idft_axes(inv, w, p)
    nz = ginv != 0
    out['maxwt'] = int(Wt[nz].max())
    out['sd_predicted'] = delta + out['maxwt']
    top = tuple([N - 1] * m)
    out['c_top'] = int(ginv[top])
    out['c_top_zero'] = out['c_top'] == 0
    out['nzeros_ginv'] = int((~nz).sum())
    # character sum for c_top: N^{-m} sum_alpha (prod alpha_i) / g(alpha)
    _, H = H_points(p, k)
    prod = H
    for _ in range(m - 1):
        prod = prod[..., None] * H % p
    cs = int((prod * inv % p).sum() % p) * pow(pow(N, m, p), -1, p) % p
    out['c_top_charsum'] = cs
    out['c_top_charsum_agrees'] = (cs == out['c_top'])
    # exact cyclic convolution g * ginv == 1 (independent of the DFT)
    conv = np.zeros_like(ginv)
    for idx in np.argwhere(nzg):
        c = int(g[tuple(idx)])
        conv = (conv + c * np.roll(ginv, tuple(int(i) for i in idx), axis=tuple(range(m)))) % p
    one = np.zeros_like(ginv)
    one[(0,) * m] = 1
    out['convolution_ok'] = bool((conv == one).all())
    return out, g, ginv


# ------------------------------------------------------------------ ordinary ring: monomials, chains, lifts
class MonoTable:
    """all monomials in n variables of total degree <= Dmax, indexed in (degree, lex) ascending order:
    index 0 is the constant monomial.  Packed as 4 bits per variable (exponents <= 15)."""

    def __init__(self, n, Dmax):
        assert Dmax <= 15
        self.n, self.Dmax = n, Dmax
        monos = []

        def rec(prefix, remaining, pos):
            if pos == n:
                monos.append(tuple(prefix))
                return
            for e in range(remaining + 1):
                prefix.append(e)
                rec(prefix, remaining - e, pos + 1)
                prefix.pop()
        rec([], Dmax, 0)
        monos.sort(key=lambda t: (sum(t), t))
        self.tuples = monos
        self.packed = [self.pack(t) for t in monos]
        self.deg = [sum(t) for t in monos]
        self.index = {pk: i for i, pk in enumerate(self.packed)}
        self.count_upto = [0] * (Dmax + 1)
        for d in self.deg:
            self.count_upto[d] += 1
        for d in range(1, Dmax + 1):
            self.count_upto[d] += self.count_upto[d - 1]

    @staticmethod
    def pack(t):
        v = 0
        for i, e in enumerate(t):
            v |= e << (4 * i)
        return v

    def var(self, i):
        return 1 << (4 * i)

    def is_squarefree(self, idx):
        return all(e <= 1 for e in self.tuples[idx])


def chain_polys(m, k, table):
    """the chains as sparse rows {monomial index: coeff} (coeff in Z, reduce mod p later)."""
    rows = []
    for i in range(m):
        for j in range(k):
            v = i * k + j
            sq = table.index[2 * table.var(v)]
            if j < k - 1:
                nxt = table.index[table.var(i * k + j + 1)]
                rows.append({sq: 1, nxt: -1})
            else:
                rows.append({sq: 1, 0: -1})
    return rows


def lift_generator(g, m, k, table, p):
    """reduced lift iota(g): x^e -> prod u_{i,j}^{bit_j(e_i)}; returns sparse row and its ordinary degree."""
    row = {}
    deg = 0
    for idx in np.argwhere(g):
        c = int(g[tuple(idx)]) % p
        if c == 0:
            continue
        pk = 0
        d = 0
        for i in range(m):
            e = int(idx[i])
            for j in range(k):
                if (e >> j) & 1:
                    pk |= 1 << (4 * (i * k + j))
                    d += 1
        row[table.index[pk]] = (row.get(table.index[pk], 0) + c) % p
        deg = max(deg, d)
    return {a: c for a, c in row.items() if c}, deg


def mul_monomial(row, mu_packed, table, p):
    """mu * row (row sparse by index); requires the products to be within Dmax."""
    out = {}
    for idx, c in row.items():
        pk = table.packed[idx] + mu_packed
        j = table.index[pk]
        out[j] = (out.get(j, 0) + c) % p
    return {a: c for a, c in out.items() if c}


class SparsePivots:
    """incremental echelon form over F_p with pivots on the largest monomial index (degree-compatible order)."""

    def __init__(self, p):
        self.p = p
        self.piv = {}

    def reduce_lead(self, row):
        """reduce the leading term repeatedly; returns the reduced row (possibly empty) whose lead is not a pivot."""
        p = self.p
        piv = self.piv
        while row:
            lead = max(row)
            pr = piv.get(lead)
            if pr is None:
                return row
            c = row[lead]
            for mnm, v in pr.items():
                nv = (row.get(mnm, 0) - c * v) % p
                if nv:
                    row[mnm] = nv
                else:
                    row.pop(mnm, None)
        return row

    def insert(self, row):
        """returns True if the row increased the span."""
        row = self.reduce_lead(dict(row))
        if not row:
            return False
        lead = max(row)
        inv = pow(row[lead], -1, self.p)
        self.piv[lead] = {a: c * inv % self.p for a, c in row.items()}
        return True

    def reduce_full_by(self, row, leads_subset):
        """subtract pivots for EVERY monomial in the row that is a pivot lead (not only the leading one).
        Terminates because pivot tails are strictly smaller than their lead."""
        p = self.p
        piv = self.piv
        row = dict(row)
        while True:
            hit = [mnm for mnm in row if mnm in piv]
            if not hit:
                return row
            mnm = max(hit)
            c = row[mnm]
            for a, v in piv[mnm].items():
                nv = (row.get(a, 0) - c * v) % p
                if nv:
                    row[a] = nv
                else:
                    row.pop(a, None)


def dense_rank_mod_p(M, p):
    """rank of an int64 matrix mod p (p small enough that p^2 * rows fits int64: p < 2^20 is plenty)."""
    M = M.copy() % p
    rows, cols = M.shape
    r = 0
    for c in range(cols):
        if r >= rows:
            break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0:
            continue
        piv = r + int(nz[0])
        if piv != r:
            M[[r, piv]] = M[[piv, r]]
        inv = pow(int(M[r, c]), -1, p)
        M[r] = M[r] * inv % p
        others = np.nonzero(M[:, c])[0]
        others = others[others != r]
        if len(others):
            f = M[others, c][:, None]
            M[others] = (M[others] - f * M[r][None, :]) % p
        r += 1
    return r
