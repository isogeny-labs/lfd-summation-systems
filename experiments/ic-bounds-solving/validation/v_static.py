#!/usr/bin/env python3
"""Validator check 2: own static-Macaulay brute force in the ORDINARY ring F_p[u_{i,j}] versus own g^{-1}.

Mac_D = span{ mu * f : f in chains u {ghat}, deg(mu f) <= D }, literally: every chain multiple and every
generator multiple of degree <= D is a row over ALL monomials of degree <= D.  Elimination is exact Gaussian
elimination over F_p with the pivot on the largest monomial in a degree-compatible order, performed as
   (1) insert every chain multiple of degree <= D into an echelon table (pivots have distinct leads);
   (2) fully reduce every generator multiple by those pivots;
   (3) 1 is in Mac_D  <=>  e_1 is in the span of the reduced generator rows
       (a nonzero combination of the chain pivots has a lead that is a pivot lead, while the reduced
        generator rows contain no pivot lead, so the chain part of any representation of e_1 is zero).
Nothing about NF(Mac_D) = g M_{<=D-delta} or about the kernel is assumed; the observations
   chain pivot count == #monomials<=D - #squarefree<=D and all chain-pivot leads non-squarefree
are RECORDED per degree as a test of Lemma 2.1(a), not used.

--dense-check additionally runs one fully dense elimination of the whole Macaulay matrix (feasible only for
tiny sizes, e.g. m = 2, k <= 2) and compares the membership answer.
"""
import argparse, json, os, random, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vcommon import (MonoTable, SparsePivots, chain_polys, lift_generator, mul_monomial, generator_and_inverse,
                     dense_rank_mod_p, is_prime)

_TABLES = {}


def table_for(n, Dmax):
    key = (n, Dmax)
    if key not in _TABLES:
        _TABLES[key] = MonoTable(n, Dmax)
    return _TABLES[key]


def static_sd(m, k, p, g, delta, Dmax=None, dense_check=False):
    n = m * k
    if Dmax is None:
        Dmax = delta + m * k
    T = table_for(n, Dmax)
    chains = chain_polys(m, k, T)
    ghat, deg_ghat = lift_generator(g, m, k, T, p)
    assert deg_ghat == delta, (deg_ghat, delta)
    sqfree_upto = [0] * (Dmax + 1)
    for idx in range(len(T.tuples)):
        if T.is_squarefree(idx):
            sqfree_upto[T.deg[idx]] += 1
    for d in range(1, Dmax + 1):
        sqfree_upto[d] += sqfree_upto[d - 1]
    by_deg = {}
    for idx, d in enumerate(T.deg):
        by_deg.setdefault(d, []).append(T.packed[idx])

    cp = SparsePivots(p)
    profile = []
    sd = None
    for D in range(0, Dmax + 1):
        t0 = time.time()
        if D >= 2:
            for mu in by_deg[D - 2]:
                for c in chains:
                    cp.insert(mul_monomial(c, mu, T, p))
        n_chain = len(cp.piv)
        leads_sqfree = sum(1 for L in cp.piv if T.is_squarefree(L))
        lemma_a_count_ok = (n_chain == T.count_upto[D] - sqfree_upto[D])
        # generator rows: every mu with deg mu <= D - delta
        rows = {}
        support_sqfree = True
        if D >= delta:
            for d in range(0, D - delta + 1):
                for mu in by_deg[d]:
                    r = cp.reduce_full_by(mul_monomial(ghat, mu, T, p), None)
                    if not r:
                        continue
                    if support_sqfree and any(not T.is_squarefree(i) for i in r):
                        support_sqfree = False
                    rows[tuple(sorted(r.items()))] = r
        one_in = False
        rankG = 0
        if rows:
            cols = sorted({i for r in rows.values() for i in r})
            cidx = {c: j for j, c in enumerate(cols)}
            G = np.zeros((len(rows), len(cols)), dtype=np.int64)
            for ri, r in enumerate(rows.values()):
                for i, c in r.items():
                    G[ri, cidx[i]] = c
            rankG = dense_rank_mod_p(G, p)
            if 0 in cidx:
                e1 = np.zeros((1, len(cols)), dtype=np.int64)
                e1[0, cidx[0]] = 1
                one_in = dense_rank_mod_p(np.vstack([G, e1]), p) == rankG
        rec = dict(D=D, chain_pivots=n_chain, chain_leads_squarefree=leads_sqfree, lemma_a_count_ok=lemma_a_count_ok,
                   gen_rows_distinct=len(rows), gen_support_squarefree=support_sqfree, rank_gen=rankG,
                   rank_Mac=n_chain + rankG, monomials=T.count_upto[D], one_in_Mac=one_in, sec=round(time.time() - t0, 2))
        if dense_check:
            rec['dense_one_in_Mac'] = dense_membership(m, k, p, T, chains, ghat, delta, D)
        profile.append(rec)
        if one_in:
            sd = D
            break
    return sd, profile


def dense_membership(m, k, p, T, chains, ghat, delta, D):
    """fully dense: all rows mu*f with deg(mu f) <= D over all monomials of degree <= D; is e_1 in the row space?"""
    ncol = T.count_upto[D]
    rows = []
    for idx in range(ncol):
        d = T.deg[idx]
        mu = T.packed[idx]
        if d + 2 <= D:
            for c in chains:
                rows.append(mul_monomial(c, mu, T, p))
        if d + delta <= D:
            rows.append(mul_monomial(ghat, mu, T, p))
    M = np.zeros((len(rows), ncol), dtype=np.int64)
    for ri, r in enumerate(rows):
        for i, c in r.items():
            M[ri, i] = c % p
    r1 = dense_rank_mod_p(M, p)
    e1 = np.zeros((1, ncol), dtype=np.int64)
    e1[0, 0] = 1
    r2 = dense_rank_mod_p(np.vstack([M, e1]), p)
    return bool(r1 == r2)


def random_unit_instance(m, k, p, rng, max_tries=5000):
    for _ in range(max_tries):
        a = rng.randrange(1, p)
        b = rng.randrange(1, p)
        if (4 * a ** 3 + 27 * b * b) % p == 0:
            continue
        t = rng.randrange(p)
        info, g, ginv = generator_and_inverse(m, k, p, a, b, t)
        if info['unit']:
            return info, g, ginv
    raise RuntimeError('no unit instance')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--m', type=int, required=True)
    ap.add_argument('--k', type=int, required=True)
    ap.add_argument('--primes', type=str, required=True)
    ap.add_argument('--per-prime', type=int, default=10)
    ap.add_argument('--seed', type=int, required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--dense-check', action='store_true')
    ap.add_argument('--instances', type=str, default=None,
                    help='optional JSON list of {p,a,b,t} to run instead of random draws')
    a = ap.parse_args()
    rng = random.Random(a.seed)
    primes = [int(x) for x in a.primes.split(',')]
    N = 1 << a.k
    todo = []
    if a.instances:
        for inst in json.load(open(a.instances)):
            todo.append((inst['p'], inst['a'], inst['b'], inst['t']))
    else:
        for p in primes:
            assert is_prime(p) and (p - 1) % N == 0, p
            for _ in range(a.per_prime):
                todo.append((p, None, None, None))
    n_mismatch = 0
    with open(a.out, 'a') as fo:
        for (p, aa, bb, tt) in todo:
            t0 = time.time()
            if aa is None:
                info, g, ginv = random_unit_instance(a.m, a.k, p, rng)
            else:
                info, g, ginv = generator_and_inverse(a.m, a.k, p, aa, bb, tt)
                if not info['unit']:
                    info['sd_static'] = None
                    fo.write(json.dumps(info) + '\n')
                    continue
            sd, profile = static_sd(a.m, a.k, p, g, info['delta'], dense_check=a.dense_check)
            info['sd_static'] = sd
            info['match'] = (sd == info['sd_predicted'])
            info['lemma_a_all_degrees'] = all(r['lemma_a_count_ok'] and r['chain_leads_squarefree'] == 0 for r in profile)
            info['gen_support_squarefree_all'] = all(r['gen_support_squarefree'] for r in profile)
            info['profile'] = profile
            if a.dense_check:
                info['dense_agrees'] = all(r['dense_one_in_Mac'] == r['one_in_Mac'] for r in profile)
            info['seconds'] = round(time.time() - t0, 2)
            n_mismatch += (not info['match'])
            fo.write(json.dumps(info) + '\n')
            fo.flush()
            print(f"m={a.m} k={a.k} p={p} a={info['a']} b={info['b']} t={info['t']} delta={info['delta']} maxwt={info['maxwt']} "
                  f"pred={info['sd_predicted']} sd_static={sd} match={info['match']} c_top_zero={info['c_top_zero']} "
                  f"lemmaA={info['lemma_a_all_degrees']} {info['seconds']}s", flush=True)
    print('mismatches:', n_mismatch)


if __name__ == '__main__':
    main()
