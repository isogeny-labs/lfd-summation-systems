#!/usr/bin/env python3
"""Random symmetric affine-bilinear Boolean systems with the shape of the summation systems.

The two-summand descent of S_3(x1, x2, x_R), x1 = sum a_i v_i, x2 = sum b_j v_j in V,
is (structure.jsonl, every one of the 2277 systems): Boolean degree 2 with every
quadratic monomial a_i b_j, symmetric under a <-> b.  For x_R != 0 the quadratic
parts phi(x1 x2), phi(u) = u^2 + x_R u, have rank <= n - 1 (ker phi = {0, x_R}),
so the span contains one affine-linear equation.  The model draws, for given
(k, m, lin, homogeneous):

  m equations  f = a^T M b + c.(a + b) + e,  M symmetric uniform, c, e uniform
                (homogeneous=True: c = 0, the order-two shape x1 x2 = const)
  lin equations  random affine-linear forms symmetric in a <-> b (l.(a + b) + e)

and certifies the exact LFD with csd.ic_local_certificate.local_profile (Z by
brute force over all 2^(2k) assignments).  Variables: a_i = bit i, b_j = bit k + j.

Usage: python3 bilinear_model.py --k 7 --m 11 --lin 1 --seeds 0 1 2 --out FILE
"""
import argparse, json, os, random, sys, time, resource, hashlib

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..')))
from csd.ic_local_certificate import local_profile  # noqa: E402


def system(k, m, lin, homogeneous, seed):
    rng = random.Random(int.from_bytes(hashlib.sha256(
        ('bilinear-model:%d:%d:%d:%d:%d' % (k, m, lin, homogeneous, seed)).encode()).digest()[:8], 'big'))
    eqs = []
    for _ in range(m):
        row = set()
        def tog(mask):
            row.symmetric_difference_update({mask})
        for i in range(k):
            for j in range(i, k):
                if rng.getrandbits(1):
                    tog((1 << i) | (1 << (k + j)))
                    if i != j:
                        tog((1 << j) | (1 << (k + i)))
        if not homogeneous:
            for i in range(k):
                if rng.getrandbits(1):
                    tog(1 << i); tog(1 << (k + i))
        if rng.getrandbits(1):
            tog(0)
        if row:
            eqs.append(sorted(row))
    for _ in range(lin):
        row = set()
        for i in range(k):
            if rng.getrandbits(1):
                row ^= {1 << i, 1 << (k + i)}
        if rng.getrandbits(1):
            row ^= {0}
        if row:
            eqs.append(sorted(row))
    return eqs


def points(eqs, r):
    out = []
    for p in range(1 << r):
        if all(sum((p & mm) == mm for mm in e) % 2 == 0 for e in eqs):
            out.append(p)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--k', type=int, required=True)
    ap.add_argument('--m', type=int, required=True)
    ap.add_argument('--lin', type=int, default=1)
    ap.add_argument('--homogeneous', type=int, default=0)
    ap.add_argument('--seeds', type=int, nargs='+', required=True)
    ap.add_argument('--max-degree', type=int, default=6)
    ap.add_argument('--max-columns', type=int, default=24000)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    r = 2 * a.k
    with open(a.out, 'a') as f:
        for s in a.seeds:
            t = time.time()
            eqs = system(a.k, a.m, a.lin, a.homogeneous, s)
            Z = points(eqs, r)
            prof = local_profile(eqs, r, Z, max_degree=a.max_degree, max_columns=a.max_columns,
                                 max_bytes=900_000_000)
            row = dict(model='symmetric-affine-bilinear', k=a.k, m=a.m, lin=a.lin, homogeneous=a.homogeneous,
                       seed=s, Z=len(Z), status=prof['status'], lfd=prof['exact_lfd'],
                       falls=prof['strict_falls'], dims=prof['quotient_ranks'],
                       seconds=round(time.time() - t, 2),
                       peak_rss_mb=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1))
            f.write(json.dumps(row) + '\n'); f.flush()
            print(json.dumps(row))


if __name__ == '__main__':
    main()
