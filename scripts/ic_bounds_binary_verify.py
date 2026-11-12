#!/usr/bin/env python3
"""Independent verifier for item 6 (experiments/ic-bounds-binary).

Written for the verification pass; it does NOT import csd.ic_bounds_binary, the
producer's scripts, csd.ic_boolean_closure, csd.ic_local_certificate or any other
closure engine of the repository.  Everything below is re-implemented:

* GF(2^n) arithmetic from the stored modulus (vectorised shift-and-add, numpy);
* the systems, from the stored specification (modulus, b, basis, target):
  the value table of S_3 (m2, m3-split) or S_4 = Res_X(S_3(x1,x2,X), S_3(x3,xR,X))
  (m3-direct) on every coordinate assignment, then the algebraic normal form by
  the Moebius transform; equation j = bit j (polynomial basis) of the ANF
  coefficients.  The ANF of a function is unique, so this is the Weil descent
  reduced modulo x_i^2 = x_i.  The system hash is recomputed in the stored format;
* Z = zeros of the value table (no Boolean equations used);
* the Boolean closure V_D (least subspace of M_<=D containing the inputs of degree
  <= D and closed under h -> x_i h (normal form) for deg h <= D - 1), by a small C
  kernel compiled at run time: semi-echelon basis over squarefree monomials ordered
  by degree descending, every basis row whose leading monomial has degree <= D-1 is
  shifted by every variable exactly once (rows are never modified after insertion,
  so on termination every low basis row has all shifts in the span: V_D is closed);
* evaluation ranks of M_<D on Z, the certificate (codim V_D = |Z| and M_<D onto Z),
  falls (dim V_D cap M_<D > dim V_{D-1}) and the exact LFD (last fall <= D*), or the
  least D with 1 in V_D for unit ideals.

A pure-Python naive closure (repeated full saturation, no semi-echelon shortcut) is
used to cross-check the C kernel on small systems (--naive-max-r).

Usage:
  ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  python3 scripts/ic_bounds_binary_verify.py rows --family m3-direct --case-ids ID ... --out F.jsonl
  python3 scripts/ic_bounds_binary_verify.py sample --plan PLAN --out F.jsonl
"""
from __future__ import annotations

import argparse
import ctypes
import hashlib
import json
import os
import subprocess
import sys
import time
import resource

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
EXP = os.path.join(ROOT, 'experiments', 'ic-bounds-binary')

C_SRC = r"""
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
typedef uint64_t u64;

static int W_, C_;
static u64 *R_; static int *LEAD_; static int *PIV_; static int N_, CAP_;

static int reduce_(u64 *v) {
  int w = 0;
  for (;;) {
    while (w < W_ && v[w] == 0) w++;
    if (w >= W_) return -1;
    int c = w * 64 + __builtin_ctzll(v[w]);
    int p = PIV_[c];
    if (p < 0) return c;
    const u64 *row = R_ + (size_t)p * W_;
    for (int j = w; j < W_; j++) v[j] ^= row[j];
  }
}

static int insert_(u64 *v) {
  int c = reduce_(v);
  if (c < 0) return -1;
  if (N_ >= CAP_) return -2;
  memcpy(R_ + (size_t)N_ * W_, v, (size_t)W_ * 8);
  LEAD_[N_] = c; PIV_[c] = N_; N_++;
  return c;
}

/* closure at degree D. coldeg[C], tgt[r*C] (-1 where degree would exceed D),
   inputs[nin*W]. out_rows[cap*W], out_lead[cap]. stats[0]=shifts, [1]=processed rows.
   returns dim or -2 on capacity overflow. */
int closure(int C, int r, int D, const int *coldeg, const int *tgt, int nin, const u64 *inputs,
            u64 *out_rows, int *out_lead, int cap, long long *stats) {
  W_ = (C + 63) / 64; C_ = C; R_ = out_rows; LEAD_ = out_lead; CAP_ = cap; N_ = 0;
  PIV_ = (int *)malloc(sizeof(int) * C);
  for (int i = 0; i < C; i++) PIV_[i] = -1;
  u64 *tmp = (u64 *)malloc((size_t)W_ * 8);
  long long shifts = 0, processed = 0;
  for (int e = 0; e < nin; e++) {
    memcpy(tmp, inputs + (size_t)e * W_, (size_t)W_ * 8);
    if (insert_(tmp) == -2) { free(tmp); free(PIV_); return -2; }
  }
  for (int idx = 0; idx < N_; idx++) {
    if (coldeg[LEAD_[idx]] > D - 1) continue;
    processed++;
    for (int i = 0; i < r; i++) {
      const u64 *src = R_ + (size_t)idx * W_;
      memset(tmp, 0, (size_t)W_ * 8);
      for (int w = 0; w < W_; w++) {
        u64 x = src[w];
        while (x) {
          int b = __builtin_ctzll(x);
          int t = tgt[(size_t)i * C + w * 64 + b];
          if (t < 0) { free(tmp); free(PIV_); return -3; }
          tmp[t >> 6] ^= (1ULL << (t & 63));
          x &= x - 1;
        }
      }
      shifts++;
      if (insert_(tmp) == -2) { free(tmp); free(PIV_); return -2; }
    }
  }
  stats[0] = shifts; stats[1] = processed;
  free(tmp); free(PIV_);
  return N_;
}

/* GF(2) rank of nrows x C bit matrix (rows packed, W words) ; destroys input */
int rank(int C, int nrows, u64 *rows, u64 *work, int cap) {
  W_ = (C + 63) / 64; C_ = C; R_ = work; CAP_ = cap; N_ = 0;
  LEAD_ = (int *)malloc(sizeof(int) * cap);
  PIV_ = (int *)malloc(sizeof(int) * C);
  for (int i = 0; i < C; i++) PIV_[i] = -1;
  for (int e = 0; e < nrows; e++) insert_(rows + (size_t)e * W_);
  free(LEAD_); free(PIV_);
  return N_;
}
"""

_LIB = None


def lib():
    global _LIB
    if _LIB is None:
        d = os.environ.get('VERIFY_BUILD', '/tmp/ic_bounds_binary_verify_build')
        os.makedirs(d, exist_ok=True)
        h = hashlib.sha256(C_SRC.encode()).hexdigest()[:12]
        so = os.path.join(d, f'vclosure_{h}.so')
        if not os.path.exists(so):
            src = os.path.join(d, f'vclosure_{h}.c')
            open(src, 'w').write(C_SRC)
            subprocess.check_call(['gcc', '-O2', '-march=native', '-shared', '-fPIC', '-o', so, src])
        L = ctypes.CDLL(so)
        P = np.ctypeslib.ndpointer
        L.closure.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, P(np.int32), P(np.int32), ctypes.c_int,
                              P(np.uint64), P(np.uint64), P(np.int32), ctypes.c_int, P(np.int64)]
        L.closure.restype = ctypes.c_int
        L.rank.argtypes = [ctypes.c_int, ctypes.c_int, P(np.uint64), P(np.uint64), ctypes.c_int]
        L.rank.restype = ctypes.c_int
        _LIB = L
    return _LIB


# ---------------------------------------------------------------------------
# GF(2^n)
# ---------------------------------------------------------------------------

class GF:
    def __init__(self, modulus: int):
        self.mod = modulus
        self.n = modulus.bit_length() - 1

    def mul(self, a, b):
        """vectorised (numpy int64 arrays or ints)"""
        a = np.asarray(a, dtype=np.int64).copy()
        b = np.asarray(b, dtype=np.int64)
        a, b = np.broadcast_arrays(a, b)
        a = a.copy()
        out = np.zeros(a.shape, dtype=np.int64)
        top = 1 << self.n
        for j in range(self.n):
            out ^= np.where((b >> j) & 1, a, 0)
            a = a << 1
            a = np.where(a & top, a ^ self.mod, a)
        return out

    def smul(self, a: int, b: int) -> int:
        out = 0
        while b:
            if b & 1:
                out ^= a
            b >>= 1
            a <<= 1
            if a >> self.n:
                a ^= self.mod
        return out

    def sqrt(self, a: int) -> int:
        for _ in range(self.n - 1):
            a = self.smul(a, a)
        return a


def span_list(basis):
    out = [0]
    for v in basis:
        out += [u ^ v for u in out]
    return np.array(out, dtype=np.int64)


def s3(F, b, x1, x2, x3):
    t = F.mul(x1, x2) ^ F.mul(x1, x3) ^ F.mul(x2, x3)
    return F.mul(t, t) ^ F.mul(F.mul(x1, x2), x3) ^ b


def s4(F, b, x1, x2, x3, x4):
    """Res_X of A X^2 + B X + C, A=(u+v)^2, B=uv, C=(uv)^2+b, pairs (x1,x2), (x3,x4)."""
    A1 = F.mul(x1 ^ x2, x1 ^ x2); B1 = F.mul(x1, x2); C1 = F.mul(B1, B1) ^ b
    A2 = F.mul(x3 ^ x4, x3 ^ x4); B2 = F.mul(x3, x4); C2 = F.mul(B2, B2) ^ b
    u = F.mul(A1, C2) ^ F.mul(A2, C1)
    return F.mul(u, u) ^ F.mul(F.mul(A1, B2) ^ F.mul(A2, B1), F.mul(B1, C2) ^ F.mul(B2, C1))


def anf(table: np.ndarray, nv: int) -> np.ndarray:
    a = table.astype(np.int64).copy()
    for i in range(nv):
        a = a.reshape(-1, 2, 1 << i)
        a[:, 1, :] ^= a[:, 0, :]
        a = a.reshape(-1)
    return a


def equations_from_anf(coef: np.ndarray, n: int, remap=None):
    eqs = []
    nz = np.nonzero(coef)[0]
    vals = coef[nz]
    masks = nz if remap is None else remap[nz]
    for j in range(n):
        sel = masks[(vals >> j) & 1 == 1]
        if len(sel):
            eqs.append(tuple(sorted(int(x) for x in sel)))
    return eqs


def hash_system(kind, variables, eqs):
    system = dict(kind=kind, variables=variables, equations=sorted(set(tuple(sorted(e)) for e in eqs if e)))
    return hashlib.sha256(json.dumps(system, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def build(row):
    """Return (r, equations, Z sorted list, my system hash)."""
    F = GF(row['modulus'])
    n = F.n
    assert n == row['n']
    b = row['b']
    xR = row['target'][0]
    basis = row['basis']
    k = len(basis)
    V = span_list(basis)
    fam = row['family']
    if fam == 'm2':
        nv = 2 * k
        p = np.arange(1 << nv, dtype=np.int64)
        val = s3(F, b, V[p & ((1 << k) - 1)], V[p >> k], np.int64(xR))
        eqs = equations_from_anf(anf(val, nv), n)
        Z = np.nonzero(val == 0)[0].tolist()
        return nv, eqs, Z, hash_system('binary', nv, eqs)
    if fam == 'm3-direct':
        nv = 3 * k
        p = np.arange(1 << nv, dtype=np.int64)
        m = (1 << k) - 1
        val = s4(F, b, V[p & m], V[(p >> k) & m], V[p >> (2 * k)], np.int64(xR))
        eqs = equations_from_anf(anf(val, nv), n)
        Z = np.nonzero(val == 0)[0].tolist()
        return nv, eqs, Z, hash_system('binary-m3-direct', nv, eqs)
    if fam == 'm3-split':
        r = 3 * k + n
        m = (1 << k) - 1
        # first equation block: variables x1 (0..k-1), x2 (k..2k-1), z (3k..3k+n-1)
        nv1 = 2 * k + n
        p = np.arange(1 << nv1, dtype=np.int64)
        val1 = s3(F, b, V[p & m], V[(p >> k) & m], p >> (2 * k))
        remap1 = (p & ((1 << (2 * k)) - 1)) | ((p >> (2 * k)) << (3 * k))
        eqs = equations_from_anf(anf(val1, nv1), n, remap1)
        # second block: z, x3 (2k..3k-1), xR
        nv2 = k + n
        q = np.arange(1 << nv2, dtype=np.int64)
        val2 = s3(F, b, q >> k, V[q & m], np.int64(xR))
        remap2 = ((q & m) << (2 * k)) | ((q >> k) << (3 * k))
        eqs += equations_from_anf(anf(val2, nv2), n, remap2)
        # Z: for each z, pairs (c1,c2) with val1 = 0 and c3 with val2 = 0
        Z = []
        z1 = val1.reshape(1 << n, 1 << (2 * k))  # index p = c12 | z << 2k
        z2 = val2.reshape(1 << n, 1 << k)
        for z in range(1 << n):
            a = np.nonzero(z1[z] == 0)[0]
            c = np.nonzero(z2[z] == 0)[0]
            if len(a) and len(c):
                for c3 in c.tolist():
                    Z.extend(int(x) | (c3 << (2 * k)) | (z << (3 * k)) for x in a.tolist())
        Z.sort()
        return r, eqs, Z, hash_system('binary-m3-split', r, eqs)
    raise ValueError(fam)


# ---------------------------------------------------------------------------
# closure
# ---------------------------------------------------------------------------

def popcount_arr(x):
    x = x.astype(np.int64)
    c = np.zeros_like(x)
    while np.any(x):
        c += x & 1
        x >>= 1
    return c


class Space:
    def __init__(self, r, D):
        self.r, self.D = r, D
        allm = np.arange(1 << r, dtype=np.int64)
        pc = popcount_arr(allm)
        sel = allm[pc <= D]
        deg = pc[pc <= D]
        order = np.lexsort((sel, -deg))  # degree descending, then mask
        self.cols = sel[order]
        self.deg = deg[order].astype(np.int32)
        self.C = len(self.cols)
        self.W = (self.C + 63) // 64
        idx = np.full(1 << r, -1, dtype=np.int32)
        idx[self.cols] = np.arange(self.C, dtype=np.int32)
        self.idx = idx
        t = np.empty((r, self.C), dtype=np.int32)
        for i in range(r):
            tm = self.cols | (1 << i)
            ok = popcount_arr(tm) <= D
            t[i] = np.where(ok, idx[np.where(ok, tm, 0)], -1)
        self.tgt = np.ascontiguousarray(t.reshape(-1))

    def pack(self, polys):
        out = np.zeros((max(1, len(polys)), self.W), dtype=np.uint64)
        for e, poly in enumerate(polys):
            for mask in poly:
                c = int(self.idx[mask])
                out[e, c >> 6] ^= np.uint64(1 << (c & 63))
        return out


def closure_dims(r, eqs, D):
    """dim V_D, dim V_D cap M_<D, 1 in V_D, stats."""
    S = Space(r, D)
    inputs = [e for e in eqs if max(bin(m).count('1') for m in e) <= D]
    packed = S.pack(inputs)
    rows = np.zeros((S.C, S.W), dtype=np.uint64)
    lead = np.zeros(S.C, dtype=np.int32)
    stats = np.zeros(2, dtype=np.int64)
    n = lib().closure(S.C, r, D, S.deg, S.tgt, len(inputs), packed, rows, lead, S.C, stats)
    if n < 0:
        raise RuntimeError(f'closure error {n}')
    ld = S.deg[lead[:n]]
    low = int(np.sum(ld <= D - 1))
    one = bool(np.any(lead[:n] == S.C - 1)) if S.deg[S.C - 1] == 0 else False
    return dict(D=D, columns=S.C, dim=n, low=low, one=one, shifts=int(stats[0]), processed=int(stats[1]))


def eval_rank(r, Z, Dlt):
    """rank of evaluation of squarefree monomials of degree < Dlt on Z."""
    if not Z or Dlt <= 0:
        return 0
    Zarr = np.array(Z, dtype=np.int64)
    allm = np.arange(1 << r, dtype=np.int64) if r <= 22 else None
    pc = popcount_arr(allm)
    mons = allm[pc < Dlt]
    nz = len(Z)
    W = (nz + 63) // 64
    M = np.zeros((len(mons), W), dtype=np.uint64)
    # bit j of row m: m subset of Z[j]
    for j0 in range(0, nz, 64):
        chunk = Zarr[j0:j0 + 64]
        bits = ((mons[:, None] & chunk[None, :]) == mons[:, None])
        weights = (np.uint64(1) << np.arange(len(chunk), dtype=np.uint64))
        M[:, j0 // 64] = np.bitwise_or.reduce(bits.astype(np.uint64) * weights[None, :], axis=1)
    work = np.zeros((min(len(mons), nz) + 1, W), dtype=np.uint64)
    M = np.ascontiguousarray(M)
    return lib().rank(nz, len(mons), M, work, min(len(mons), nz) + 1)


def certify(r, eqs, Z, Dmax, extra=0):
    """Exact LFD by the definitions; scans D = 1..(certificate degree + extra) up to Dmax."""
    from math import comb
    out = dict(dims=[0], lows=[0], evalrank_lt=[0], falls=[], certificate_degree=None, exact_lfd=None,
               unit=None, closure=[])
    # V_0: constant inputs
    const = [e for e in eqs if e == (0,)]
    out['dims'][0] = 1 if const else 0
    stop = None
    for D in range(1, Dmax + 1):
        t0 = time.time()
        c = closure_dims(r, eqs, D)
        c['seconds'] = round(time.time() - t0, 2)
        out['closure'].append(c)
        out['dims'].append(c['dim'])
        out['lows'].append(c['low'])
        if c['low'] > out['dims'][D - 1]:
            out['falls'].append(D)
        er = eval_rank(r, Z, D)
        out['evalrank_lt'].append(er)
        size = sum(comb(r, i) for i in range(D + 1))
        codim = size - c['dim']
        er_le = eval_rank(r, Z, D + 1)
        assert codim >= er_le, ('soundness violated', D, codim, er_le)
        if not Z:
            if c['one']:
                assert c['dim'] == size
                if out['exact_lfd'] is None:
                    out['exact_lfd'] = D
                    out['certificate_degree'] = D
                    out['unit'] = True
                    stop = D + extra
        else:
            if out['certificate_degree'] is None and er == len(Z) and codim == len(Z):
                out['certificate_degree'] = D
                out['exact_lfd'] = max(out['falls']) if out['falls'] else 0
                out['unit'] = False
                stop = D + extra
        if stop is not None and D >= stop:
            break
    if out['certificate_degree'] is not None and out['falls'] and max(out['falls']) > out['certificate_degree']:
        out['fall_after_certificate'] = True  # would contradict the local-to-global lemma
    return out


# ---------------------------------------------------------------------------
# naive cross-check (pure Python, small r)
# ---------------------------------------------------------------------------

def naive_dims(r, eqs, D):
    from itertools import combinations
    cols = [m for d in range(D, -1, -1) for m in
            sorted(sum(1 << i for i in c) for c in combinations(range(r), d))]
    idx = {m: i for i, m in enumerate(cols)}
    C = len(cols)
    deg = [bin(m).count('1') for m in cols]

    def to_int(poly):
        v = 0
        for m in poly:
            v ^= 1 << (C - 1 - idx[m])  # high bit = first column (highest degree)
        return v

    def rref(vs):
        piv = {}
        for v in vs:
            for p in sorted(piv, reverse=True):
                if v >> p & 1:
                    v ^= piv[p]
            if v:
                p = v.bit_length() - 1
                for q in list(piv):
                    if piv[q] >> p & 1:
                        piv[q] ^= v
                piv[p] = v
        return piv

    def mul(v, i):
        out = 0
        while v:
            p = v.bit_length() - 1
            v ^= 1 << p
            m = cols[C - 1 - p] | (1 << i)
            out ^= 1 << (C - 1 - idx[m])
        return out

    basis = rref([to_int(e) for e in eqs if max(bin(m).count('1') for m in e) <= D])
    while True:
        low = [v for p, v in basis.items() if deg[C - 1 - p] <= D - 1]
        new = list(basis.values()) + [mul(v, i) for v in low for i in range(r)]
        nb = rref(new)
        if len(nb) == len(basis):
            break
        basis = nb
    low = sum(1 for p in basis if deg[C - 1 - p] <= D - 1)
    return len(basis), low


# ---------------------------------------------------------------------------
# drivers
# ---------------------------------------------------------------------------

def load(family):
    return [json.loads(l) for l in open(os.path.join(EXP, f'results-{family}.jsonl'))]


def verify_row(row, naive_max_r=0, extra=0):
    t0 = time.time()
    r, eqs, Z, h = build(row)
    tb = time.time() - t0
    zsha = hashlib.sha256(json.dumps(Z).encode()).hexdigest()[:16]
    Dmax = (row.get('certificate_degree') or row.get('scanned_through') or 8) + extra
    cert = certify(r, eqs, Z, Dmax, extra=extra)
    res = dict(case_id=row['case_id'], family=row['family'], n=row['n'], k=row['k'], a=row['a'], b=row['b'],
               base_kind=row['base_kind'], target_class=row['target_class'], target=row['target'], r=r,
               hash_match=(h == row['system_hash']), my_system_hash=h, equations=len(set(eqs)),
               solutions=len(Z), stored_solutions=row['solutions'], points_sha256_mine=zsha,
               points_sha256_stored=row.get('points_sha256'),
               lfd=cert['exact_lfd'], stored_lfd=row['exact_lfd'],
               dstar=cert['certificate_degree'], stored_dstar=row.get('certificate_degree'),
               dims=cert['dims'], stored_dims=row.get('quotient_ranks'),
               falls=cert['falls'], stored_falls=row.get('strict_falls'), evalrank_lt=cert['evalrank_lt'],
               unit=cert['unit'], stored_status=row['status'], build_seconds=round(tb, 2),
               closure=cert['closure'], fall_after_certificate=cert.get('fall_after_certificate', False))
    sd = row.get('quotient_ranks') or []
    res['dims_match'] = all(sd[i] == cert['dims'][i] for i in range(1, min(len(sd), len(cert['dims']))))
    res['lfd_match'] = res['lfd'] == res['stored_lfd']
    res['dstar_match'] = res['dstar'] == res['stored_dstar']
    if naive_max_r and r <= naive_max_r:
        nv = []
        for D in range(1, (cert['certificate_degree'] or Dmax) + 1):
            d, lo = naive_dims(r, eqs, D)
            nv.append([d, lo])
        res['naive'] = nv
        res['naive_match'] = all(nv[i][0] == cert['dims'][i + 1] and nv[i][1] == cert['lows'][i + 1]
                                 for i in range(len(nv)))
    res['total_seconds'] = round(time.time() - t0, 2)
    res['peak_rss_mb'] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['rows', 'all'])
    ap.add_argument('--family', required=True)
    ap.add_argument('--case-ids', nargs='+', default=None)
    ap.add_argument('--k', type=int, nargs='+', default=None)
    ap.add_argument('--n', type=int, nargs='+', default=None)
    ap.add_argument('--out', required=True)
    ap.add_argument('--naive-max-r', type=int, default=0)
    ap.add_argument('--extra', type=int, default=0, help='also scan this many degrees past D*')
    ap.add_argument('--extra-max-r', type=int, default=99, help='apply --extra only when r <= this')
    a = ap.parse_args()
    allrows = load(a.family)
    rows = {r['case_id']: r for r in allrows}
    if a.mode == 'rows':
        ids = a.case_ids
    else:
        ids = [r['case_id'] for r in allrows if (a.k is None or r['k'] in a.k) and (a.n is None or r['n'] in a.n)]
    done, by_hash = set(), {}
    if os.path.exists(a.out):
        for l in open(a.out):
            x = json.loads(l)
            done.add(x['case_id'])
            if x.get('closure_from') is None:
                by_hash[x['my_system_hash']] = x
    for cid in ids:
        if cid in done:
            continue
        row = rows[cid]
        if row['system_hash'] in by_hash:
            # same system already certified in this run: rebuild + hash-check only, reuse the closure
            r, eqs, Z, h = build(row)
            prev = by_hash[row['system_hash']]
            res = dict(prev)
            res.update(case_id=cid, n=row['n'], k=row['k'], a=row['a'], b=row['b'], base_kind=row['base_kind'],
                       target_class=row['target_class'], target=row['target'], hash_match=(h == row['system_hash']),
                       my_system_hash=h, solutions=len(Z), stored_solutions=row['solutions'],
                       points_sha256_mine=hashlib.sha256(json.dumps(Z).encode()).hexdigest()[:16],
                       points_sha256_stored=row.get('points_sha256'), stored_lfd=row['exact_lfd'],
                       stored_dstar=row.get('certificate_degree'), stored_dims=row.get('quotient_ranks'),
                       stored_falls=row.get('strict_falls'), stored_status=row['status'], closure_from=prev['case_id'])
            res['lfd_match'] = res['lfd'] == res['stored_lfd']
            res['dstar_match'] = res['dstar'] == res['stored_dstar']
            sd = row.get('quotient_ranks') or []
            res['dims_match'] = all(sd[i] == res['dims'][i] for i in range(1, min(len(sd), len(res['dims']))))
            res['solutions_match_prev'] = len(Z) == prev['solutions']
        else:
            r0 = row['variables']
            res = verify_row(row, a.naive_max_r, a.extra if r0 <= a.extra_max_r else 0)
            res['closure_from'] = None
            if res['hash_match']:
                by_hash[res['my_system_hash']] = res
        res['points_match'] = res['points_sha256_mine'] == res['points_sha256_stored']
        with open(a.out, 'a') as f:
            f.write(json.dumps(res, sort_keys=True) + '\n')
        print(cid, res['family'], 'n', res['n'], 'k', res['k'], '|Z|', res['solutions'], '/', res['stored_solutions'],
              'hash', res['hash_match'], 'pts', res['points_match'], 'LFD', res['lfd'], '/', res['stored_lfd'],
              'D*', res['dstar'], '/', res['stored_dstar'], 'dims', res['dims_match'], 'naive', res.get('naive_match'),
              'reuse' if res.get('closure_from') else f"{res['total_seconds']}s {res['peak_rss_mb']}MB", flush=True)


if __name__ == '__main__':
    main()
