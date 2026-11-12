#!/usr/bin/env python3
"""Verifier engine for ITEM 'Three or more summands' (written independently of
csd/ic_bounds_m3.py, csd/ic_bounds_m3_batched.py and every repository closure engine).

Closure V_{F,D} computed modulo the chains (squarefree masks; x_i^e with e mod 2^k, degree
= popcount; multiplication by u_{i,j} adds 2^j to e_i mod 2^k), plus an optional free variable z
(split encoding, exponent unrestricted).  Exact arithmetic mod p in float64 (all intermediate
integers < 2^53).  Incremental reduced row echelon form, columns ordered by degree DESCENDING, so
V cap M_<=e is spanned by the rows whose pivot has degree <= e.

Exact LFD (local-to-global, R/I computed by CRT over the chain variety H^m):
  direct/m2 : R/I = F_p^Z, Z = zeros of NF(g) on H^m   (R/C = F_p^{H^m}, so (g) is radical)
  split     : R/I = (+)_{x in H^3} F_p[z]/gcd(g1(x,z), g2(x,z))
For each D: dim(I cap M_<=D) = |M_<=D| - rank(ev: M_<=D -> R/I); the closure saturation stops
early only when its rank reaches that number (then V_D = I cap M_<=D exactly), otherwise it runs
to a fixed point.  Certified at the first D with ev(M_<D) onto R/I and V_D = I cap M_<=D;
LFD = last strict fall <= D.  Unit ideal: least D with 1 in V_D.
"""
import json, math, sys, time, hashlib, argparse, random
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
VDIR = ROOT / 'experiments' / 'ic-bounds-m3' / 'verify'

# ------------------------------------------------------------------ field / polynomials

def s4_terms():
    d = json.loads((VDIR / 's4_coeffs.json').read_text())
    return d['terms']


def poly_s4_at(p, a, b, t, terms):
    """S_4(x1,x2,x3,t) mod p as dict (e1,e2,e3)->c."""
    out = {}
    for (ea, eb, e1, e2, e3, e4), c in terms:
        v = c * pow(a, ea, p) * pow(b, eb, p) * pow(t, e4, p) % p
        if v:
            out[(e1, e2, e3)] = (out.get((e1, e2, e3), 0) + v) % p
    return {m: c for m, c in out.items() if c}


def s3_dict(p, a, b):
    """S_3(u,v,w) as dict (eu,ev,ew)->c, expanded by hand."""
    out = {}
    def add(m, c):
        out[m] = (out.get(m, 0) + c) % p
    # (u-v)^2 w^2
    add((2, 0, 2), 1); add((1, 1, 2), -2); add((0, 2, 2), 1)
    # -2((u+v)(uv+a)+2b) w = -2(u^2 v + u v^2 + a u + a v + 2b) w
    add((2, 1, 1), -2); add((1, 2, 1), -2); add((1, 0, 1), -2 * a); add((0, 1, 1), -2 * a); add((0, 0, 1), -4 * b)
    # (uv-a)^2 - 4b(u+v)
    add((2, 2, 0), 1); add((1, 1, 0), -2 * a); add((0, 0, 0), a * a); add((1, 0, 0), -4 * b); add((0, 1, 0), -4 * b)
    return {m: c % p for m, c in out.items() if c % p}


def subgroup(p, k):
    assert (p - 1) % (1 << k) == 0
    for h in range(2, p):
        g = pow(h, (p - 1) >> k, p)
        if k == 0 or pow(g, 1 << (k - 1), p) != 1:
            H = [pow(g, e, p) for e in range(1 << k)]
            assert len(set(H)) == 1 << k
            return H
    raise ValueError


def is_prime(n):
    return n > 1 and all(n % q for q in range(2, int(n ** .5) + 1))

# polynomial (list low->high) helpers mod p
def ptrim(f):
    while f and f[-1] == 0: f.pop()
    return f

def pmod(f, g, p):
    f = f[:]; inv = pow(g[-1], p - 2, p); dg = len(g) - 1
    while len(f) - 1 >= dg and f:
        c = f[-1] * inv % p; s = len(f) - 1 - dg
        for i in range(dg + 1): f[s + i] = (f[s + i] - c * g[i]) % p
        ptrim(f)
    return f

def pgcd(f, g, p):
    f, g = ptrim(f[:]), ptrim(g[:])
    while g:
        f, g = g, pmod(f, g, p)
    if f:
        inv = pow(f[-1], p - 2, p); f = [c * inv % p for c in f]
    return f

def pmulmod(f, g, h, p):
    r = [0] * (len(f) + len(g) - 1) if f and g else []
    for i, x in enumerate(f):
        if x:
            for j, y in enumerate(g): r[i + j] = (r[i + j] + x * y) % p
    return pmod(ptrim(r), h, p)

def zpow_mod(e, h, p):
    """z^e mod h."""
    res, base = [1], pmod([0, 1], h, p)
    res = pmod(res, h, p)
    while e:
        if e & 1: res = pmulmod(res, base, h, p)
        base = pmulmod(base, base, h, p); e >>= 1
    return res

# ------------------------------------------------------------------ linear algebra mod p

def rank_mod(M, p):
    M = np.array(M, dtype=np.float64) % p
    r = 0; rows, cols = M.shape
    for c in range(cols):
        if r == rows: break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0: continue
        i = r + nz[0]
        if i != r: M[[r, i]] = M[[i, r]]
        M[r] = (M[r] * pow(int(M[r, c]), p - 2, p)) % p
        col = M[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr):
            M[nzr] = np.mod(M[nzr] - np.outer(col[nzr], M[r]), p)
        r += 1
    return r


class Echelon:
    """Reduced row echelon basis mod p, pivot = leftmost nonzero (columns degree-descending)."""

    def __init__(self, n, p, cap):
        self.n, self.p = n, p
        self.B = np.zeros((cap, n))
        self.piv = []
        self.r = 0

    def grow(self, need):
        if need > self.B.shape[0]:
            nb = np.zeros((max(need, 2 * self.B.shape[0]) if self.B.shape[0] else need, self.n))
            nb = nb[:min(nb.shape[0], self.n)] if need <= self.n else nb
            nb[:self.r] = self.B[:self.r]; self.B = nb

    def add(self, N):
        """Insert rows of N; return the inserted reduced rows (RREF, leftmost pivots)."""
        p = self.p
        if N.shape[0] == 0: return np.zeros((0, self.n)), []
        if self.r:
            pv = np.array(self.piv)
            N = np.mod(N - N[:, pv] @ self.B[:self.r], p)
        # Gauss-Jordan on N (small), leftmost pivot each step (in place)
        N = N[np.any(N != 0, axis=1)].copy()
        pivs = []; r = 0; n = self.n
        while r < N.shape[0]:
            sub = N[r:] != 0
            nzr = np.any(sub, axis=1)
            if not nzr.any(): break
            lead = np.argmax(sub, axis=1); lead[~nzr] = n + 1
            i = r + int(np.argmin(lead)); c = int(lead[i - r])
            if i != r: N[[r, i]] = N[[i, r]]
            N[r] = np.mod(N[r] * pow(int(N[r, c]), p - 2, p), p)
            col = N[:, c].copy(); col[r] = 0
            nz = np.nonzero(col)[0]
            if len(nz): N[nz] = np.mod(N[nz] - np.outer(col[nz], N[r]), p)
            pivs.append(c); r += 1
        rows = list(N[:r])
        if not rows: return np.zeros((0, self.n)), []
        Nr = np.array(rows)
        if self.r:
            coef = self.B[:self.r][:, pivs]
            if np.any(coef):
                self.B[:self.r] = np.mod(self.B[:self.r] - coef @ Nr, p)
        self.grow(self.r + len(rows))
        self.B[self.r:self.r + len(rows)] = Nr
        self.r += len(rows); self.piv.extend(pivs)
        return Nr, pivs

# ------------------------------------------------------------------ monomial spaces

def popc(e):
    return bin(e).count('1')


class Space:
    def __init__(self, m, k, zvar, D):
        N = 1 << k
        mons = []
        def rec(i, cur, w):
            if i == m:
                if zvar:
                    for ez in range(D - w + 1): mons.append(tuple(cur) + (ez,))
                else:
                    mons.append(tuple(cur))
                return
            for e in range(N):
                pw = popc(e)
                if w + pw <= D:
                    cur.append(e); rec(i + 1, cur, w + pw); cur.pop()
        rec(0, [], 0)
        def wt(mm): return sum(popc(e) for e in mm[:m]) + (mm[m] if zvar else 0)
        mons.sort(key=lambda mm: (-wt(mm), mm))
        self.m, self.k, self.zvar, self.D = m, k, zvar, D
        self.mons = mons
        self.idx = {mm: i for i, mm in enumerate(mons)}
        self.w = np.array([wt(mm) for mm in mons])
        self.n = len(mons)
        # variable maps for columns of weight <= D-1
        src = np.nonzero(self.w <= D - 1)[0]
        self.vmaps = []
        for i in range(m):
            for j in range(k):
                tg = []
                for c in src:
                    mm = list(mons[c]); mm[i] = (mm[i] + (1 << j)) % N; tg.append(self.idx[tuple(mm)])
                self.vmaps.append((src, np.array(tg, dtype=np.int64)))
        if zvar:
            tg = []
            for c in src:
                mm = list(mons[c]); mm[m] += 1; tg.append(self.idx[tuple(mm)])
            self.vmaps.append((src, np.array(tg, dtype=np.int64)))

    def vec(self, poly, p):
        v = np.zeros(self.n)
        for mm, c in poly.items(): v[self.idx[mm]] = c % p
        return v


def nf(poly, m, k, zvar, p):
    """Chain normal form: exponents of the m block variables mod 2^k (z untouched)."""
    N = 1 << k; out = {}
    for mm, c in poly.items():
        key = tuple(e % N for e in mm[:m]) + (tuple(mm[m:]) if zvar else ())
        out[key] = (out.get(key, 0) + c) % p
    return {kk: c for kk, c in out.items() if c}


def weight(mm, m, zvar):
    return sum(popc(e) for e in mm[:m]) + (mm[m] if zvar else 0)

# ------------------------------------------------------------------ systems

def build_system(enc, p, a, b, t, k, terms):
    if enc == 'direct':
        g = nf(poly_s4_at(p, a, b, t, terms), 3, k, False, p)
        return dict(m=3, zvar=False, inputs=[g])
    s3 = s3_dict(p, a, b)
    if enc == 'm2':
        g = {}
        for (eu, ev, ew), c in s3.items():
            key = (eu, ev); g[key] = (g.get(key, 0) + c * pow(t, ew, p)) % p
        g = nf({kk: c for kk, c in g.items() if c}, 2, k, False, p)
        return dict(m=2, zvar=False, inputs=[g])
    if enc == 'split':
        g1 = {(eu, ev, 0, ew): c for (eu, ev, ew), c in s3.items()}            # S_3(x1,x2,z)
        g2 = {}
        for (eu, ev, ew), c in s3.items():                                    # S_3(z,x3,t)
            key = (0, 0, ev, eu); g2[key] = (g2.get(key, 0) + c * pow(t, ew, p)) % p
        g2 = {kk: c for kk, c in g2.items() if c}
        return dict(m=3, zvar=True, inputs=[nf(g1, 3, k, True, p), nf(g2, 3, k, True, p)])
    raise ValueError(enc)


def mono_val(mm, x, p, m):
    v = 1
    for i in range(m): v = v * pow(x[i], mm[i], p) % p
    return v


def quotient_algebra(sysd, p, k, H, field=False):
    """Return list of components (x, h_x) with h_x the local polynomial in z (split) or None."""
    import itertools
    m = sysd['m']
    pts = list(itertools.product(H, repeat=m))
    if not sysd['zvar']:
        g = sysd['inputs'][0]
        Z = [x for x in pts if sum(c * mono_val(mm, x, p, m) for mm, c in g.items()) % p == 0]
        return [(x, None) for x in Z]
    comps = []
    for x in pts:
        locs = []
        for g in sysd['inputs']:
            f = {}
            for mm, c in g.items():
                f[mm[m]] = (f.get(mm[m], 0) + c * mono_val(mm, x, p, m)) % p
            d = max(f) if f else 0
            locs.append(ptrim([f.get(i, 0) for i in range(d + 1)]))
        h = pgcd(locs[0], locs[1], p)
        if not h:
            raise RuntimeError('infinite quotient at %r' % (x,))
        if field and len(h) > 1:
            zp = zpow_mod(p, h, p)
            zp = zp + [0] * max(0, 2 - len(zp)); zp[1] = (zp[1] - 1) % p
            h = pgcd(h, ptrim(zp), p)
        if len(h) > 1: comps.append((x, h))
    return comps


def eval_rows(sp, comps, p, m):
    """Matrix |M_<=D| x dim A of images of monomials in R/I."""
    if not comps: return np.zeros((sp.n, 0))
    if comps[0][1] is None:
        M = np.zeros((sp.n, len(comps)))
        for j, (x, _) in enumerate(comps):
            for i, mm in enumerate(sp.mons): M[i, j] = mono_val(mm, x, p, m)
        return M
    blocks = []
    for x, h in comps:
        d = len(h) - 1
        maxz = max(mm[m] for mm in sp.mons)
        zp = [zpow_mod(e, h, p) for e in range(maxz + 1)]
        Bk = np.zeros((sp.n, d))
        for i, mm in enumerate(sp.mons):
            v = mono_val(mm, x, p, m)
            if v:
                r = zp[mm[m]]
                for j, c in enumerate(r): Bk[i, j] = v * c % p
        blocks.append(Bk)
    return np.hstack(blocks)


def closure(sp, inputs, seeds, p, target=None, deadline=None, verify_closed=False):
    """V_D from seeds (rows already in V_D) and inputs; returns Echelon."""
    D = sp.D
    E = Echelon(sp.n, p, cap=min(sp.n, 64))
    init = [sp.vec(g, p) for g in inputs if max(weight(mm, sp.m, sp.zvar) for mm in g) <= D]
    if seeds is not None and seeds.shape[0]:
        init = list(seeds) + init
    queue = []
    if init:
        Nr, pv = E.add(np.array(init))
        low = Nr[sp.w[np.array(pv)] <= D - 1] if len(pv) else Nr
        if low.shape[0]: queue.append(low)
    nvars = len(sp.vmaps)
    while queue:
        if target is not None and E.r >= target: break
        if deadline and time.time() > deadline: raise TimeoutError
        L = queue.pop()
        if L.shape[0] > 24:
            queue.append(L[24:]); L = L[:24]
        P = np.zeros((L.shape[0] * nvars, sp.n))
        for vi, (src, tg) in enumerate(sp.vmaps):
            P[vi * L.shape[0]:(vi + 1) * L.shape[0]][:, tg] = L[:, src]
        Nr, pv = E.add(P)
        if len(pv):
            low = Nr[sp.w[np.array(pv)] <= D - 1]
            if low.shape[0]: queue.append(low)
    if verify_closed and (target is None or E.r < target):
        lowB = E.B[:E.r][sp.w[np.array(E.piv)] <= D - 1] if E.r else np.zeros((0, sp.n))
        for s in range(0, lowB.shape[0], 64):
            L = lowB[s:s + 64]
            for src, tg in sp.vmaps:
                P = np.zeros((L.shape[0], sp.n)); P[:, tg] = L[:, src]
                red = np.mod(P - P[:, np.array(E.piv)] @ E.B[:E.r], p)
                assert not np.any(red), 'closure not closed'
    return E


def exact_profile(enc, p, a, b, t, k, terms, Dmax=14, deadline=None, verify_closed=False, field=False):
    sysd = build_system(enc, p, a, b, t, k, terms)
    m, zvar = sysd['m'], sysd['zvar']
    H = subgroup(p, k)
    comps = quotient_algebra(sysd, p, k, H)
    dimA = sum(1 if h is None else len(h) - 1 for _, h in comps)
    out = dict(enc=enc, p=p, a=a, b=b, t=t, k=k, dimA=dimA,
               input_degrees=[max(weight(mm, m, zvar) for mm in g) if g else None for g in sysd['inputs']],
               profile=[])
    if zvar:
        out['dimA_field'] = sum(len(h) - 1 for _, h in quotient_algebra(sysd, p, k, H, field=True))
    prev_rank, prev_basis, prev_sp = 0, None, None
    falls = []
    for D in range(0, Dmax + 1):
        sp = Space(m, k, zvar, D)
        Ev = eval_rows(sp, comps, p, m)
        rk_le = rank_mod(Ev, p) if dimA else 0
        rk_lt = rank_mod(Ev[sp.w <= D - 1], p) if dimA and D else 0
        target = sp.n - rk_le
        seeds = None
        if prev_basis is not None and prev_basis.shape[0]:
            seeds = np.zeros((prev_basis.shape[0], sp.n))
            cols = np.array([sp.idx[mm] for mm in prev_sp.mons])
            seeds[:, cols] = prev_basis
        E = closure(sp, sysd['inputs'], seeds, p, target=target, deadline=deadline, verify_closed=verify_closed)
        rank = E.r
        assert rank <= target, 'closure exceeds I cap M_<=D: unsound'
        low = int(np.sum(sp.w[np.array(E.piv, dtype=int)] <= D - 1)) if rank else 0
        exact_V_eq_I = rank == target
        if D >= 1 and low > prev_rank: falls.append(D)
        unit = dimA == 0 and rank == sp.n
        out['profile'].append(dict(D=D, n=sp.n, rank=rank, low=low, target=target, rank_ev_lt=rk_lt))
        if dimA == 0 and rank == sp.n:
            out.update(status='unit', lfd=D, falls=falls, cert_D=D); return out
        if dimA and rk_lt == dimA and exact_V_eq_I:
            out.update(status='proper', lfd=(falls[-1] if falls else 0), falls=falls, cert_D=D); return out
        prev_rank, prev_basis, prev_sp = rank, E.B[:E.r].copy(), sp
    out.update(status='capped', lower=(falls[-1] if falls else None), falls=falls)
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--enc'); ap.add_argument('--p', type=int); ap.add_argument('--a', type=int)
    ap.add_argument('--b', type=int); ap.add_argument('--t', type=int); ap.add_argument('--k', type=int)
    ap.add_argument('--verify-closed', action='store_true')
    A = ap.parse_args()
    t0 = time.time()
    r = exact_profile(A.enc, A.p, A.a, A.b, A.t, A.k, s4_terms(), verify_closed=A.verify_closed)
    r['seconds'] = round(time.time() - t0, 2)
    print(json.dumps(r))
