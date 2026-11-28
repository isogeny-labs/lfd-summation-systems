#!/usr/bin/env sage -python
"""Online single-target relation solve with sparse graph-style elimination."""

import json
from random import Random
from time import perf_counter

from lazy_walk import lazy_oracle
from replay import E, F, K, bits, cold_bsgs, cold_rho, generator, prime


class OnlineTargetElimination:
    """Eliminate base-log columns while retaining the one target column."""

    def __init__(self):
        self.pivots = {}
        self.rows = 0

    def add(self, x1, sign1, x2, sign2, target_coeff, rhs):
        self.rows += 1
        coeff = {}
        for x, sign in ((x1, sign1), (x2, sign2)):
            coeff[x] = (coeff.get(x, 0)+sign) % prime
            if not coeff[x]:
                del coeff[x]
        target_coeff %= prime
        rhs %= prime
        while coeff:
            pivot = min(coeff, key=bits)
            if pivot not in self.pivots:
                inv = pow(coeff[pivot], -1, prime)
                self.pivots[pivot] = (
                    {x: (c*inv) % prime for x, c in coeff.items()},
                    target_coeff*inv % prime, rhs*inv % prime)
                return None
            pivot_coeff, pivot_target, pivot_rhs = self.pivots[pivot]
            factor = coeff[pivot]
            for x, c in pivot_coeff.items():
                coeff[x] = (coeff.get(x, 0)-factor*c) % prime
                if not coeff[x]:
                    del coeff[x]
            target_coeff = (target_coeff-factor*pivot_target) % prime
            rhs = (rhs-factor*pivot_rhs) % prime
        if target_coeff:
            return rhs*pow(target_coeff, -1, prime) % prime
        assert rhs == 0  # A dependent, consistent relation.
        return None


def run_online(G, T, target_scalar, seed, max_queries=1024):
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

    cache = {}
    elim = OnlineTargetElimination()
    hit_queries = 0
    t1 = perf_counter()
    for query in range(1, max_queries+1):
        slot = 0 if Q.is_zero() else bits(Q[0]) % 8
        point, da, db = steps[slot]
        Q += point
        a, b = (a+da) % prime, (b+db) % prime
        hits = lazy_oracle(Q, cache, as_roots)
        hit_queries += bool(hits)
        for x1, s1, x2, s2 in hits:
            answer = elim.add(x1, s1, x2, s2, -2*b, 2*a)
            if answer is not None:
                assert answer == target_scalar and answer*G == T
                query_seconds = perf_counter()-t1
                return {
                    "seed": seed, "queries": query,
                    "queries_with_relation": hit_queries,
                    "failed_queries": query-hit_queries,
                    "relation_rows": elim.rows,
                    "base_pivots": len(elim.pivots),
                    "lazy_lifted_abscissas": len(cache),
                    "recovered_scalar": answer,
                    "increment_setup_seconds": setup_seconds,
                    "online_query_and_solve_seconds": query_seconds,
                    "cold_total_seconds": setup_seconds+query_seconds,
                }
    raise AssertionError("target not determined within frozen query cap")


def main():
    G = generator()
    runs = {}
    for i, scalar in enumerate((3425, 17123, 29113)):
        T = scalar*G
        runs[str(scalar)] = {
            "target": [bits(T[0]), bits(T[1])],
            "cold_bsgs": cold_bsgs(G, T),
            "cold_rho": cold_rho(G, T, 20261127+i),
            "online_walk": run_online(G, T, scalar, 20261127+i),
        }
    print(json.dumps({"field_modulus": str(K.modulus()),
                      "curve_order": E.cardinality(),
                      "prime_order": prime,
                      "generator": [bits(G[0]), bits(G[1])],
                      "runs": runs}, indent=2, sort_keys=True, default=int))


if __name__ == "__main__":
    main()
