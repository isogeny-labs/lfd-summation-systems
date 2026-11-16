"""Static Macaulay solving degree sd = least D with 1 in Mac_D (unit instances), computed in the
ordinary ring.  Two exact modes (both mod p, both independent of g^{-1}):

 mode 'dense' : literal Mac_D = span{m*f : f in chains+{g}, deg(m f) <= D} in the full monomial
                space P_{<=D} of K[u_{i,j}] (r = 2k vars), reduced row echelon via lfd_engine.Ech.
                Feasible for k = 3 (and k = 4 only for D <= 7).
 mode 'nf'    : the same row space, with the chain rows eliminated by the rewriting system
                u_{i,j}^2 -> u_{i,j+1}, u_{i,k-1}^2 -> 1 (confluent: coprime leading terms; rewriting never
                raises degree, so ker(NF) cap P_{<=D} = span of the chain rows of degree <= D).
                Rows m*g for EVERY ordinary monomial m with deg m <= D-2 are NF-reduced (monomial
                products, then rewriting) and the resulting vectors on the 2^{2k} squarefree monomials are
                echelonised.  Cross-checked against 'dense' at k = 3 (all D) in self_test().
Usage: python3 static_bruteforce.py k p [n_instances] [seed]   (writes data/bf_k{k}_p{p}.json)
[pcc:unverified]
"""
import sys, os, json, time, random
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from common import *


def g_terms_mask(a, b, t, p, k):
    return [(e1 | (e2 << k), c) for (e1, e2), c in gen_coeffs(a, b, t, p).items() if c % p]


# ---------------------------------------------------------------- dense mode
def mac_dense_contains_one(a, b, t, p, k, D, want_dim=False):
    r = 2 * k
    sp = LE.Space(r, D)
    inputs, G = LE.prime_inputs(g_terms_mask(a, b, t, p, k), k, p)
    ech = LE.Ech(sp.n, p)
    for f in inputs:
        df = max(sum(m) for m in f)
        if df > D: continue
        mons = [m for m in LE.monomials(r, D - df)]
        rows = np.zeros((len(mons), sp.n))
        for i, m in enumerate(mons):
            for e, c in f.items():
                rows[i, sp.idx[tuple(x + y for x, y in zip(m, e))]] = c % p
        ech.add(rows)
    one = np.zeros(sp.n); one[sp.n - 1] = 1  # constant monomial is last column
    res = ech.contains(one)
    return (res, ech.rk) if want_dim else res


# ---------------------------------------------------------------- nf mode
def nf_exps(exps, k):
    """rewrite ordinary exponent vector (len 2k) to squarefree mask via u_j^2->u_{j+1}, u_{k-1}^2->1."""
    mask = 0
    for i in range(2):
        carry = 0
        for j in range(k):
            e = exps[i * k + j] + carry
            if e & 1: mask |= 1 << (i * k + j)
            carry = e >> 1  # at level k-1 the carry becomes a power of 1
    return mask


def mulmod_row(m_exps, gt, k, p, ncols):
    row = np.zeros(ncols, dtype=np.int64)
    for gm, c in gt:
        ge = [(gm >> i) & 1 for i in range(2 * k)]
        row[nf_exps([x + y for x, y in zip(m_exps, ge)], k)] += c
    return row % p


class Echelon:
    def __init__(self, ncols, p):
        self.n, self.p, self.piv = ncols, p, {}

    def reduce(self, v):
        v = v.copy() % self.p
        while True:
            nz = np.flatnonzero(v)
            if nz.size == 0: return None
            lead = int(nz[0])
            if lead in self.piv:
                v = (v - v[lead] * self.piv[lead]) % self.p
            else:
                return lead, v

    def add(self, v):
        r = self.reduce(v)
        if r is None: return False
        lead, v = r
        self.piv[lead] = v * pow(int(v[lead]), -1, self.p) % self.p
        return True

    def contains(self, v):
        return self.reduce(v) is None


_NF_CACHE = {}


def nf_targets(k, Dmax):
    """For every ordinary monomial m with deg m <= Dmax-2 (r = 2k vars): the 9-tuple of squarefree masks
    NF(m * x^gmask) over the nine monomials of g.  Identical tuples give identical rows for every instance
    (rows are c_gm-weighted scatters), so keep the least degree per tuple.  Returns list of
    (mindeg, tuple) sorted by mindeg."""
    key = (k, Dmax)
    if key in _NF_CACHE: return _NF_CACHE[key]
    r = 2 * k
    gms = [e1 | (e2 << k) for (e1, e2) in sorted(gen_coeffs(1, 1, 1, 10 ** 9 + 7))]
    gexp = [[(gm >> i) & 1 for i in range(r)] for gm in gms]
    best = {}
    for m in LE.monomials(r, Dmax - 2):
        d = sum(m)
        tg = tuple(nf_exps([x + y for x, y in zip(m, ge)], k) for ge in gexp)
        if tg not in best or best[tg] > d: best[tg] = d
    out = (gms, sorted((d, tg) for tg, d in best.items()))
    _NF_CACHE[key] = out
    return out


def sd_nf(a, b, t, p, k, Dmax=None):
    """returns (sd, dims) ; sd=None if 1 not in Mac_D up to Dmax (default 2k+2)"""
    Dmax = Dmax or 2 * k + 2
    ncols = 1 << (2 * k)
    gms, rows_all = nf_targets(k, Dmax)
    coef = gen_coeffs(a, b, t, p)
    cv = [coef[(gm & ((1 << k) - 1), gm >> k)] % p for gm in gms]
    ech = Echelon(ncols, p)
    one = np.zeros(ncols, dtype=np.int64); one[0] = 1
    dims = {}; pos = 0
    for D in range(2, Dmax + 1):
        while pos < len(rows_all) and rows_all[pos][0] <= D - 2:
            row = np.zeros(ncols, dtype=np.int64)
            for c, tgt in zip(cv, rows_all[pos][1]): row[tgt] += c
            ech.add(row % p); pos += 1
        dims[D] = len(ech.piv)
        if ech.contains(one): return D, dims
    return None, dims


def lfd_per_degree(a, b, t, p, k, Dmax, deadline=None):
    """exact d_lf, d_ff, per-degree closure dims with lfd_engine (ordinary ring) -- unit instances only."""
    inputs, G = LE.prime_inputs(g_terms_mask(a, b, t, p, k), k, p)
    res = LE.exact_lfd(inputs, 2 * k, p, [], Dmax, deadline=deadline)
    falls = res['falls']
    return dict(d_lf=res['lfd'], d_ff=(min(falls) if falls else None), falls=falls,
                dimV=[x['dimV'] for x in res['per_degree']], one_in_V=[x['one_in_V'] for x in res['per_degree']],
                seconds=[x['seconds'] for x in res['per_degree']])


def self_test():
    """dense == nf on k=3 for every D: dim Mac_D(dense) = (#monos<=D - #squarefree<=D) + dim NF-rows, same 1-membership."""
    import math
    rng = random.Random(5)
    for p in (97, 193, 113 if False else 257):
        k = 3
        for (a, b, t) in unit_instances(p, k, 3, rng)[0]:
            sd, dims = sd_nf(a, b, t, p, k, Dmax=10)
            for D in range(2, 9):
                one, rk = mac_dense_contains_one(a, b, t, p, k, D, want_dim=True)
                nmono = math.comb(2 * k + D, D)
                nsf = sum(math.comb(2 * k, w) for w in range(0, min(D, 2 * k) + 1))
                assert one == (sd is not None and D >= sd), (p, a, b, t, D, one, sd)
                assert rk == (nmono - nsf) + dims.get(D, 0), (D, rk, nmono - nsf, dims.get(D))
            print("self_test k=3 p=%d (a,b,t)=%s sd=%s: dense==nf for D=2..8 ok" % (p, (a, b, t), sd), flush=True)


if __name__ == '__main__':
    if sys.argv[1] == 'selftest':
        self_test(); sys.exit()
    k, p = int(sys.argv[1]), int(sys.argv[2])
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    rng = random.Random(seed * 1000003 + p * 31 + k)
    inst, tries = unit_instances(p, k, n, rng)
    out = dict(k=k, p=p, seed=seed, unit_tries=tries, rows=[])
    for (a, b, t) in inst:
        t0 = time.time()
        sd, dims = sd_nf(a, b, t, p, k)
        out['rows'].append(dict(p=p, a=a, b=b, t=t, k=k, sd_bruteforce=sd, mac_rank=dims, sec=round(time.time() - t0, 2)))
        print(k, p, (a, b, t), 'sd', sd, round(time.time() - t0, 1), flush=True)
    json.dump(out, open(os.path.join(os.path.dirname(__file__), '..', 'data', 'bf_k%d_p%d.json' % (k, p)), 'w'), indent=1)
