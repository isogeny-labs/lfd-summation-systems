#!/usr/bin/env sage -python
"""Cold single-target screen with one-addition queries and lazy point lifts."""

import json
from random import Random
from time import perf_counter

from sage.all import GF, matrix, vector

from replay import E, F, K, b6, bits, cold_bsgs, cold_rho, generator, prime, roots_at, s3


def lazy_oracle(Q, cache, as_roots):
    """Return actual signed relations, lifting only encountered x-values."""
    if Q.is_zero():
        return []
    t, tq = Q[0], Q[0]**256
    assert t != tq
    A = F(t+tq)
    B = F(b6+b6**256)
    N = F(t*tq)
    C = B/A
    D = F(C*C+t*C+b6)
    alpha, beta = N/(A*A), D/(A*A)

    def lift(x):
        if x not in cache:
            ys = roots_at(x)
            cache[x] = E(x, ys[0]) if ys and x != 0 else None
        return cache[x]

    out = {}
    for v in as_roots.get(beta/(alpha*alpha), ()):
        u = alpha*v
        s, p = u.sqrt(), A*u+C
        xs = ([K(p.sqrt())]*2 if s == 0 else
              sorted((K(s*w) for w in as_roots.get(p/(s*s), ())), key=bits))
        if len(xs) != 2:
            continue
        x1, x2 = xs
        P1, P2 = lift(x1), lift(x2)
        if P1 is None or P2 is None:
            continue
        for P, sign1 in ((P1, 1), (-P1, -1)):
            R = Q-P
            if R == P2:
                sign2 = 1
            elif R == -P2:
                sign2 = -1
            else:
                continue
            assert P+R == Q and s3(x1, x2, t) == 0
            out[(bits(x1), sign1, bits(x2), sign2)] = (x1, sign1, x2, sign2)
    assert len(out) <= 2
    return list(out.values())


def run_walk(G, T, target_scalar, seed, max_queries=1024):
    rng = Random(seed)
    t0 = perf_counter()
    as_roots = {}
    for v in F:
        as_roots.setdefault(v*v+v, (v, v+F.one()))
    assert len(as_roots) == 128
    steps = []
    for _ in range(8):
        da, db = rng.randrange(prime), rng.randrange(prime)
        steps.append((da*G+db*T, da, db))
    a, b = rng.randrange(prime), rng.randrange(prime)
    Q = a*G+b*T
    setup_seconds = perf_counter()-t0

    cache, rows, rhs, row_end, elapsed = {}, [], [], [], []
    hit_queries = 0
    query_start = perf_counter()
    for _ in range(max_queries):
        j = 0 if Q.is_zero() else bits(Q[0]) % 8
        point, da, db = steps[j]
        Q += point
        a, b = (a+da) % prime, (b+db) % prime
        hits = lazy_oracle(Q, cache, as_roots)
        hit_queries += bool(hits)
        for x1, s1, x2, s2 in hits:
            rows.append((x1, s1, x2, s2, -2*b))
            rhs.append(2*a)
        row_end.append(len(rows))
        elapsed.append(perf_counter()-query_start)

    Fr = GF(prime)

    def prefix_matrix(q):
        prefix = rows[:row_end[q-1]]
        xs = sorted({x for x1, _, x2, _, _ in prefix for x in (x1, x2)}, key=bits)
        col = {x: i for i, x in enumerate(xs)}
        A = matrix(Fr, len(prefix), len(xs)+1)
        for i, (x1, s1, x2, s2, target_coeff) in enumerate(prefix):
            A[i, col[x1]] += s1
            A[i, col[x2]] += s2
            A[i, len(xs)] = target_coeff
        return A, len(xs)

    def determined(q):
        A, n = prefix_matrix(q)
        return A.rank() > A.matrix_from_columns(range(n)).rank()

    assert determined(max_queries)
    lo, hi = 1, max_queries
    while lo < hi:
        mid = (lo+hi)//2
        if determined(mid):
            hi = mid
        else:
            lo = mid+1
    t1 = perf_counter()
    A, used_x = prefix_matrix(lo)
    answer = A.solve_right(vector(Fr, rhs[:row_end[lo-1]]))
    assert int(answer[-1]) == target_scalar
    assert A*answer == vector(Fr, rhs[:row_end[lo-1]])
    solve_seconds = perf_counter()-t1
    return {
        "seed": seed, "queries_examined": max_queries,
        "queries_with_relation": hit_queries,
        "first_target_queries": lo, "first_target_rows": row_end[lo-1],
        "used_abscissas_at_target": used_x,
        "target_recovered": int(answer[-1]),
        "walk_setup_seconds": setup_seconds,
        "query_seconds_at_target": elapsed[lo-1],
        "solve_seconds": solve_seconds,
        "optimistic_cold_seconds": setup_seconds+elapsed[lo-1]+solve_seconds,
        "lazy_lifts_across_full_stream": len(cache),
    }


def main():
    G = generator()
    runs = {}
    for i, d in enumerate((3425, 17123, 29113)):
        T = d*G
        runs[str(d)] = {
            "target": [bits(T[0]), bits(T[1])],
            "cold_bsgs": cold_bsgs(G, T),
            "cold_rho": cold_rho(G, T, 20261127+i),
            "lazy_walk": run_walk(G, T, d, 20261127+i),
        }
    print(json.dumps({"field_modulus": str(K.modulus()), "curve_order": E.cardinality(),
                      "prime_order": prime, "generator": [bits(G[0]), bits(G[1])],
                      "runs": runs}, indent=2, sort_keys=True, default=int))


if __name__ == "__main__":
    main()
