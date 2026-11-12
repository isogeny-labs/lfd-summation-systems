#!/usr/bin/env python3
"""Independent exact last-fall-degree (LFD) engine, written from the definition in
formal/ICFallDegrees/Closure.lean.  It does NOT import any closure engine of the
repository (ic_fall_certificates, ic_chain_shift_closure, ic_boolean_closure,
ic_curve_sweep profiles).  Arithmetic is exact mod p (p prime) or GF(2) (p = 2),
in the ORDINARY polynomial ring K[z_1..z_r] (full monomial space of degree <= D).

Definitions implemented
-----------------------
V_{F,D}  : least K-subspace of P_{<=D} that contains every input of total degree <= D and
           is closed under h -> x*h for every variable x whenever deg h + 1 <= D.
           (variable-wise rule; closure_naive() below implements the literal Lean rule
           "a*f for every a, f in V with deg a + deg f <= D" to cross-check it.)
I_D      : I cap P_{<=D}, the ideal elements of degree <= D.  For a zero-dimensional ideal
           whose zero set Z is all K-rational (true for chain/Boolean systems) and which
           is radical, I_D = kernel of evaluation P_{<=D} -> K^Z, so
           dim I_D = |M_{<=D}| - rank(eval on M_{<=D}).
LFD      : least d with  I_D = V_D  for all D >= d  (this is HasLFDBound: every h in I lies
           in V_{max(d, deg h)}; V_D is contained in I_D by soundness).
Certificate of  LFD <= D0 : dim V_{D0} = dim I_{D0}  and evaluation of the monomials of
           degree < D0 is onto K^Z.  Proof: for h in I_{D0+1} write h_top = sum x_i g_i,
           choose q_i with deg q_i < D0 and g_i = q_i on Z; then g_i - q_i in I_{D0} = V_{D0},
           x_i(g_i - q_i) in V_{D0+1}, and h - sum x_i (g_i - q_i) in I_{D0} = V_{D0} <= V_{D0+1};
           induct.  Unit ideal (Z empty): LFD = least D with 1 in V_D (then all monomials of
           degree <= D are in V_D).
Fall at D (D >= 1): V_D cap P_{<=D-1} strictly larger than V_{D-1}.  LFD = largest fall.
"""
import itertools
import math
import numpy as np


# --------------------------------------------------------------------------------------
# monomial spaces
# --------------------------------------------------------------------------------------

def monomials(r, D):
    out = []

    def rec(i, left, cur):
        if i == r:
            out.append(tuple(cur))
            return
        for e in range(left + 1):
            cur.append(e)
            rec(i + 1, left - e, cur)
            cur.pop()

    rec(0, D, [])
    out.sort(key=lambda m: (sum(m), m), reverse=True)  # degree DESCENDING: leading column = top degree
    return out


class Space:
    """Columns = monomials of degree <= D, sorted by degree descending."""

    def __init__(self, r, D):
        self.r, self.D = r, D
        self.cols = monomials(r, D)
        self.n = len(self.cols)
        self.idx = {m: i for i, m in enumerate(self.cols)}
        self.deg = np.array([sum(m) for m in self.cols], dtype=np.int64)
        below = np.flatnonzero(self.deg <= D - 1)
        self.low_start = int(below[0]) if below.size else self.n  # first column of degree <= D-1
        self.mul = []
        for i in range(r):
            arr = np.empty(self.n - self.low_start, dtype=np.int64)
            for t in range(self.low_start, self.n):
                m = list(self.cols[t])
                m[i] += 1
                arr[t - self.low_start] = self.idx[tuple(m)]
            self.mul.append(arr)

    def shift(self, rows, i):
        """x_i * rows, rows have support only in columns of degree <= D-1."""
        out = np.zeros((rows.shape[0], self.n), dtype=rows.dtype)
        out[:, self.mul[i]] = rows[:, self.low_start:]
        return out

    def row_of(self, poly, p):
        """poly: dict exps->coeff.  Returns a row (float) or raises if degree > D."""
        row = np.zeros(self.n, dtype=np.float64)
        for m, c in poly.items():
            c %= p
            if c:
                row[self.idx[m]] = c
        return row


# --------------------------------------------------------------------------------------
# row echelon over F_p (dense, BLAS-assisted, exact)
# --------------------------------------------------------------------------------------

class Ech:
    """Reduced row echelon basis; pivot = first nonzero column (= leading monomial in the
    degree-descending column order)."""

    def __init__(self, n, p, cap=None):
        self.n, self.p = n, p
        self.dt = np.float32 if p == 2 else np.float64
        self.cap = cap or n
        self.Ebuf = np.zeros((self.cap, n), dtype=self.dt)
        self.rk = 0
        self.pivbuf = np.zeros(self.cap, dtype=np.int64)
        assert self.cap * (p - 1) ** 2 < (2 ** 24 if self.dt == np.float32 else 2 ** 53)

    @property
    def E(self):
        return self.Ebuf[:self.rk]

    @property
    def piv(self):
        return self.pivbuf[:self.rk]

    def reduce(self, N):
        N = N.astype(self.dt, copy=False)
        if self.rk == 0:
            return N
        return np.mod(N - N[:, self.piv] @ self.E, self.p)

    def add(self, N, chunk=192):
        """Insert rows N (values in [0,p)).  Returns list of the newly created RREF rows."""
        p = self.p
        created = []
        for a in range(0, N.shape[0], chunk):
            C = self.reduce(N[a:a + chunk])
            C = C[np.any(C != 0, axis=1)]
            if C.shape[0] == 0:
                continue
            L = np.zeros((C.shape[0], self.n), dtype=self.dt)
            pl = []
            cnt = 0
            for row in C:
                v = row.copy()
                if cnt:
                    v = np.mod(v - v[pl] @ L[:cnt], p)
                nz = np.flatnonzero(v)
                if nz.size == 0:
                    continue
                lead = int(nz[0])
                inv = pow(int(v[lead]), -1, p)
                v = np.mod(v * inv, p)
                if cnt:
                    L[:cnt] = np.mod(L[:cnt] - np.outer(L[:cnt, lead], v), p)
                L[cnt] = v
                pl.append(lead)
                cnt += 1
            L = L[:cnt]
            pl = np.array(pl, dtype=np.int64)
            for a2 in range(0, self.rk, 1024):   # blockwise to bound temporaries
                blk = self.Ebuf[a2:min(a2 + 1024, self.rk)]
                blk[...] = np.mod(blk - blk[:, pl] @ L, p)
            self.Ebuf[self.rk:self.rk + cnt] = L
            self.pivbuf[self.rk:self.rk + cnt] = pl
            self.rk += cnt
            created.append(L.copy())
        return created

    def contains(self, vec):
        r = self.reduce(vec.reshape(1, -1))
        return not np.any(r)


def rank_mod(M, p):
    M = np.asarray(M)
    if M.shape[0] == 0 or M.shape[1] == 0:
        return 0
    e = Ech(M.shape[1], p, cap=min(M.shape[0], M.shape[1]))
    e.add(np.mod(M, p))
    return e.rk


# --------------------------------------------------------------------------------------
# V_{F,D}
# --------------------------------------------------------------------------------------

def closure(inputs, r, D, p, space=None, max_rounds=10 ** 6):
    """inputs: list of dict exps->coeff.  Returns dict with dim, dim_low (= dim V_D cap P_{<=D-1}),
    contains_one, rounds, and the Ech object."""
    sp = space or Space(r, D)
    V = Ech(sp.n, p)
    done = Ech(sp.n, p, cap=max(1, sp.n - sp.low_start))
    rows = []
    for f in inputs:
        d = max((sum(m) for m, c in f.items() if c % p), default=-1)
        if d <= D:
            rows.append(sp.row_of(f, p))
    if rows:
        V.add(np.array(rows))
    rounds = 0
    while rounds < max_rounds:
        lowmask = V.piv >= sp.low_start  # pivot of degree <= D-1  <=>  row has degree <= D-1
        low = V.E[lowmask]
        if low.shape[0] == 0:
            break
        resid = done.reduce(low)
        resid = resid[np.any(resid != 0, axis=1)]
        if resid.shape[0] == 0:
            break
        new = done.add(resid)
        rounds += 1
        for blk in new:
            for i in range(r):
                V.add(sp.shift(blk, i))
    lowcount = int(np.sum(V.piv >= sp.low_start))
    one = bool(sp.n - 1 in set(V.piv.tolist()))
    return dict(dim=V.rk, dim_low=lowcount, contains_one=one, rounds=rounds, V=V, space=sp)


def closure_naive(inputs, r, D, p):
    """Literal Lean rule on a small instance: close under a*f for every MONOMIAL a and every
    element f of V with deg a + deg f <= D.  Pure Python, exact integer arithmetic, own echelon."""
    cols = monomials(r, D)
    idx = {m: i for i, m in enumerate(cols)}
    n = len(cols)
    deg = [sum(m) for m in cols]

    piv = {}  # lead col -> row (list of ints mod p), normalised lead = 1

    def insert(row):
        row = [x % p for x in row]
        while True:
            lead = next((j for j, x in enumerate(row) if x), None)
            if lead is None:
                return False
            if lead in piv:
                c = row[lead]
                pr = piv[lead]
                row = [(a - c * b) % p for a, b in zip(row, pr)]
            else:
                inv = pow(row[lead], -1, p)
                piv[lead] = [(x * inv) % p for x in row]
                return True

    def to_row(poly):
        row = [0] * n
        for m, c in poly.items():
            row[idx[m]] = (row[idx[m]] + c) % p
        return row

    for f in inputs:
        d = max((sum(m) for m, c in f.items() if c % p), default=-1)
        if d <= D:
            insert(to_row(f))
    monos = {e: [m for m in cols if sum(m) <= e] for e in range(D + 1)}
    changed = True
    while changed:
        changed = False
        for lead, row in list(piv.items()):
            # row degree = degree of its leading column (graded order)
            dr = deg[lead]
            for a in monos[D - dr] if D - dr >= 0 else []:
                if sum(a) == 0:
                    continue
                new = [0] * n
                for j, x in enumerate(row):
                    if x:
                        m = tuple(u + v for u, v in zip(cols[j], a))
                        new[idx[m]] = x
                if insert(new):
                    changed = True
        # also rows' linear combinations are covered because every pivot row (a basis) is multiplied
        # but elements of V of low degree that are NOT basis rows are combinations of basis rows
        # with pivot degree <= their degree, so a*f is a combination of a*basis rows with budget
        # deg a + deg(pivot) <= deg a + deg f <= D.
    low = sum(1 for lead in piv if deg[lead] <= D - 1)
    one = (n - 1) in piv
    return dict(dim=len(piv), dim_low=low, contains_one=one)


# --------------------------------------------------------------------------------------
# evaluation and exact LFD
# --------------------------------------------------------------------------------------

def eval_matrix(points, sp, p):
    """rows = points (tuples of r field elements), columns = monomials of sp."""
    r = sp.r
    ex = np.array(sp.cols, dtype=np.int64)  # n x r
    M = np.ones((len(points), sp.n), dtype=np.int64)
    for a, pt in enumerate(points):
        for j in range(r):
            pw = np.array([pow(int(pt[j]), e, p) if (pt[j] or e) else 1 for e in range(sp.D + 1)], dtype=np.int64)
            M[a] = (M[a] * pw[ex[:, j]]) % p
    return M


def exact_lfd(inputs, r, p, points, Dmax, verbose=False, deadline=None, closure_fn=None):
    """Scan D = 0..Dmax.  Returns dict with per-degree data and the verdict."""
    import time
    Z = len(points)
    per = []
    sp_max = Space(r, Dmax)
    Mfull = eval_matrix(points, sp_max, p) if Z else np.zeros((0, sp_max.n), dtype=np.int64)
    rank_eval = []
    for d in range(Dmax + 1):
        cols = sp_max.deg <= d
        rank_eval.append(rank_mod(Mfull[:, cols], p) if Z else 0)
    nmono = [int(np.sum(sp_max.deg <= d)) for d in range(Dmax + 1)]
    cert = None
    for D in range(Dmax + 1):
        t0 = time.time()
        res = (closure_fn or closure)(inputs, r, D, p)
        dimI = nmono[D] - rank_eval[D]
        rec = dict(D=D, monomials=nmono[D], dimV=res['dim'], dimV_low=res['dim_low'], dimI=dimI,
                   rank_eval=rank_eval[D], one_in_V=res['contains_one'], rounds=res.get('rounds'),
                   seconds=round(time.time() - t0, 2))
        assert res['dim'] <= dimI, ('V_D not inside I_D: soundness violated', rec)
        per.append(rec)
        if verbose:
            print(rec, flush=True)
        onto_below = (D >= 1 and rank_eval[D - 1] == Z) or (Z == 0)
        if rec['dimV'] == dimI and onto_below:
            cert = D
            break
        if deadline and time.time() > deadline:
            break
    falls = [D for D in range(1, len(per)) if per[D]['dimV_low'] > per[D - 1]['dimV']]
    out = dict(points=Z, per_degree=per, falls=falls, certified_at=cert)
    if cert is not None:
        neq = [rec['D'] for rec in per if rec['dimV'] != rec['dimI']]
        out['lfd'] = (max(neq) + 1) if neq else 0
        out['lfd_largest_fall'] = max(falls) if falls else 0
    else:
        out['lfd'] = None
    return out


# --------------------------------------------------------------------------------------
# system builders (own code)
# --------------------------------------------------------------------------------------

def prime_inputs(gen_items, k, p):
    """Ordinary-ring inputs in r = 2k variables u_{i,j} = variable i*k+j:
    chain u_{i,j}^2 - u_{i,j+1}, u_{i,k-1}^2 - 1, and the generator G (squarefree chain monomials)."""
    r = 2 * k
    inputs = []
    for i in range(2):
        for j in range(k):
            sq = [0] * r
            sq[i * k + j] = 2
            f = {tuple(sq): 1}
            if j + 1 < k:
                nxt = [0] * r
                nxt[i * k + j + 1] = 1
                f[tuple(nxt)] = (p - 1)
            else:
                f[tuple([0] * r)] = p - 1
            inputs.append(f)
    G = {}
    for mask, c in gen_items:
        e = tuple((mask >> b) & 1 for b in range(r))
        G[e] = (G.get(e, 0) + c) % p
    G = {e: c for e, c in G.items() if c}
    inputs.append(G)
    return inputs, G


def mu_N(p, N):
    return [x for x in range(1, p) if pow(x, N, p) == 1]


def prime_points(G, k, p):
    """Zero set of chain + G in F_p^{2k}: x_i in mu_{2^k}, u_{i,j} = x_i^{2^j}."""
    N = 1 << k
    mus = mu_N(p, N)
    assert len(mus) == N
    pts = []
    for x1 in mus:
        u1 = [pow(x1, 1 << j, p) for j in range(k)]
        for x2 in mus:
            u2 = [pow(x2, 1 << j, p) for j in range(k)]
            u = u1 + u2
            v = 0
            for e, c in G.items():
                t = c
                for j, ej in enumerate(e):
                    if ej:
                        t = t * pow(u[j], ej, p) % p
                v = (v + t) % p
            if v == 0:
                pts.append(tuple(u))
    return pts


def binary_inputs(equations, r):
    inputs = []
    for j in range(r):
        sq = [0] * r
        sq[j] = 2
        lin = [0] * r
        lin[j] = 1
        inputs.append({tuple(sq): 1, tuple(lin): 1})  # z^2 - z = z^2 + z over GF(2)
    for eq in equations:
        f = {}
        for mask in eq:
            e = tuple((mask >> b) & 1 for b in range(r))
            f[e] = (f.get(e, 0) + 1) % 2
        f = {e: c for e, c in f.items() if c}
        inputs.append(f)
    return inputs


def binary_points(equations, r):
    pts = []
    for bits in range(1 << r):
        ok = True
        for eq in equations:
            v = 0
            for mask in eq:
                v ^= int((bits & mask) == mask)
            if v:
                ok = False
                break
        if ok:
            pts.append(tuple((bits >> b) & 1 for b in range(r)))
    return pts
