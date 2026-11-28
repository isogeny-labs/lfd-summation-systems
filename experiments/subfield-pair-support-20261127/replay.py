#!/usr/bin/env sage -python
"""Exact two-summand relation and cold baseline replay over F_(2^16)."""

import json
from collections import Counter, defaultdict
from math import isqrt
from random import Random
from time import perf_counter

from sage.all import EllipticCurve, GF, PolynomialRing, matrix, vector


R = PolynomialRing(GF(2), "w")
w = R.gen()
K = GF(2**16, "z", modulus=w**16 + w**5 + w**3 + w**2 + 1)
z = K.gen()
F = K.subfield(8, "u")
a2, b6 = z**11, z
E = EllipticCurve(K, [1, a2, 0, 0, b6])
order, prime = 65986, 32993
assert E.cardinality() == order == 2 * prime
assert (1 / b6) ** 256 != 1 / b6  # The j-invariant forbids F_(2^8) descent.
Yring = PolynomialRing(K, "Y")
Y = Yring.gen()


def bits(x):
    return sum(int(c) << i for i, c in enumerate(x.polynomial().list()))


def field_element(n):
    return sum(((n >> i) & 1) * z**i for i in range(16))


def roots_at(x):
    return sorted((y for y, mult in
                   (Y**2 + x*Y + x**3 + a2*x*x + b6).roots()), key=bits)


def canonical_base(xs):
    out = []
    for x in sorted(xs, key=bits):
        assert x.trace() == 0  # Both candidate spaces share this parity control.
        if x == 0:
            continue  # The point at x=0 has order two and no prime-part log.
        ys = roots_at(x)
        if ys:
            assert len(ys) == 2
            P = E(x, ys[0])
            assert -P == E(x, ys[1])
            assert prime * P != E(0)
            out.append(P)
    return out


def generator():
    for n in range(1, 2**16):
        x = field_element(n)
        ys = roots_at(x)
        if ys:
            G = 2 * E(x, ys[0])
            if G != E(0):
                assert prime * G == E(0)
                return G
    raise AssertionError("missing generator")


def binary_rank(vals):
    pivots = {}
    for x in vals:
        v = bits(x)
        while v:
            p = v.bit_length()-1
            if p not in pivots:
                pivots[p] = v
                break
            v ^= pivots[p]
    return len(pivots)


def product_rank(basis):
    return binary_rank([basis[i]*basis[j] for i in range(len(basis))
                        for j in range(i, len(basis))])


def pair_table(base):
    signed = [(P, i, 1) for i, P in enumerate(base)]
    signed += [(-P, i, -1) for i, P in enumerate(base)]
    table = defaultdict(list)
    for i, (P, pi, si) in enumerate(signed):
        for Q, qi, ti in signed[i:]:
            S = P + Q
            assert prime * S == E(0)  # Both summands share the same 2-coset.
            table[S].append((pi, si, qi, ti, P[0], Q[0]))
    return table


def s3(x, y, t):
    return (x*y + x*t + y*t)**2 + x*y*t + b6


def relation_run(base, G, target_scalar, seed, base_seconds, max_queries=1024):
    t0 = perf_counter()
    table = pair_table(base)
    table_seconds = perf_counter() - t0
    T = target_scalar * G
    rng = Random(seed)
    rows, rhs, row_end, query_cumulative = [], [], [], []
    found_queries = 0
    query_start = perf_counter()
    for q in range(max_queries):
        a = rng.randrange(prime)
        b = rng.randrange(1, prime)
        Q = a*G + b*T
        hits = table.get(Q, ())
        found_queries += bool(hits)
        for i, si, j, sj, x, y in hits:
            if Q != E(0):
                assert s3(x, y, Q[0]) == 0
            row = [0] * (len(base) + 1)
            row[i] += si
            row[j] += sj
            row[-1] = -2*b
            rows.append(row)
            rhs.append(2*a)
        row_end.append(len(rows))
        query_cumulative.append(perf_counter() - query_start)
    Fr = GF(prime)
    if matrix(Fr, rows).rank() != len(base) + 1:
        raise AssertionError("full log system did not reach rank")
    lo, hi = 1, len(row_end)
    while lo < hi:
        mid = (lo + hi) // 2
        if matrix(Fr, rows[:row_end[mid-1]]).rank() == len(base) + 1:
            hi = mid
        else:
            lo = mid + 1
    needed = lo
    t1 = perf_counter()
    M = matrix(Fr, rows[:row_end[needed-1]])
    rank = M.rank()
    answer = M.solve_right(vector(Fr, rhs[:row_end[needed-1]]))
    assert int(answer[-1]) == target_scalar
    assert M * answer == vector(Fr, rhs[:row_end[needed-1]])
    la_seconds = perf_counter() - t1
    query_seconds = query_cumulative[needed-1]
    return {
        "base_logs": len(base), "signed_points": 2*len(base),
        "unordered_pair_additions": len(base)*(2*len(base)+1),
        "support": len(table), "pair_multiplicities": dict(sorted(Counter(map(len, table.values())).items())),
        "queries_examined": max_queries, "queries_with_relation": found_queries,
        "full_rank_queries": needed, "full_rank_rows": row_end[needed-1],
        "rank": int(rank), "recovered_scalar": int(answer[-1]),
        "base_seconds": base_seconds, "table_seconds": table_seconds,
        "query_seconds_at_full_rank": query_seconds,
        "linear_solve_seconds": la_seconds,
        "cold_total_seconds": base_seconds + table_seconds + query_seconds + la_seconds,
    }


def cold_bsgs(G, T):
    m = isqrt(prime) + 1
    t0 = perf_counter()
    babies = {}
    Q = E(0)
    for j in range(m):
        babies.setdefault(Q, j)
        Q += G
    step = -(m*G)
    Q = T
    for i in range(m + 1):
        if Q in babies:
            d = (i*m + babies[Q]) % prime
            assert d*G == T
            return {"recovered_scalar": d, "baby_steps": m,
                    "giant_steps": i, "seconds": perf_counter()-t0}
        Q += step
    raise AssertionError("BSGS failed")


def cold_rho(G, T, seed):
    rng = Random(seed)
    t0 = perf_counter()

    def step(state):
        X, a, b = state
        bucket = 0 if X.is_zero() else bits(X[0]) % 3
        if bucket == 0:
            return X + G, (a + 1) % prime, b
        if bucket == 1:
            return X + T, a, (b + 1) % prime
        return 2*X, (2*a) % prime, (2*b) % prime

    for restart in range(4):
        a, b = rng.randrange(prime), rng.randrange(prime)
        start = (a*G + b*T, a, b)
        slow, fast = step(start), step(step(start))
        for i in range(1, 10001):
            if slow[0] == fast[0]:
                den = (slow[2] - fast[2]) % prime
                if den:
                    d = ((fast[1] - slow[1]) * pow(den, -1, prime)) % prime
                    assert d*G == T
                    return {"recovered_scalar": d, "restart": restart,
                            "floyd_iterations": i,
                            "walk_steps": 3*i,
                            "seconds": perf_counter() - t0}
                break
            slow = step(slow)
            fast = step(step(fast))
    raise AssertionError("rho failed after four deterministic restarts")


def main():
    t0 = perf_counter()
    sub = canonical_base([K(x) for x in F])
    sub_seconds = perf_counter() - t0
    t0 = perf_counter()
    mono = canonical_base([sum(((i >> j) & 1)*z**j for j in range(8))
                           for i in range(256)])
    mono_seconds = perf_counter() - t0
    assert len(sub) == 128 and len(mono) == 129
    mono = mono[:128]  # Drop one complete signed x-class; equalize unknown count.
    dims = {"subfield": product_rank([K(F.gen()**i) for i in range(8)]),
            "monomial_control": product_rank([z**i for i in range(8)])}
    assert dims == {"subfield": 8, "monomial_control": 15}
    G = generator()
    scalars = [20261127 % prime, 17123, 29113]
    runs = {}
    for i, target_scalar in enumerate(scalars):
        T = target_scalar * G
        runs[str(target_scalar)] = {
            "target": [bits(T[0]), bits(T[1])],
            "cold_bsgs": cold_bsgs(G, T),
            "cold_rho": cold_rho(G, T, 20261127+i),
            "subfield": relation_run(sub, G, target_scalar, 20261127+i, sub_seconds),
            "monomial_control": relation_run(mono, G, target_scalar, 20261127+i, mono_seconds),
        }
    result = {
        "field_modulus": str(K.modulus()), "curve": {"a1": 1, "a2": bits(a2), "a6": bits(b6)},
        "order": order, "prime_subgroup_order": prime,
        "generator": [bits(G[0]), bits(G[1])],
        "product_space_dimensions": dims,
        "target_scalars_for_verification": scalars,
        "runs": runs,
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
