#!/usr/bin/env sage -python
"""Exact online quadratic-S3 cost panel across three quadratic extensions."""

import json
import platform
from math import isqrt
from random import Random
from time import perf_counter, process_time

from sage.all import EllipticCurve, GF, PolynomialRing, is_prime
from sage.version import version as sage_version


SEED = 20261127
CASES = (
    dict(n=16, k=8, modulus=(16, 5, 3, 2, 0), a2_exp=11,
         b6_bits=2, prime=32993),
    dict(n=20, k=10, modulus=(20, 10, 9, 7, 6, 5, 4, 1, 0), a2_exp=11,
         b6_bits=878565, prime=524189),
    dict(n=24, k=12, modulus=(24, 16, 15, 14, 13, 10, 9, 7, 5, 3, 0),
         a2_exp=9, b6_bits=4389102, prime=8387251),
)


class OnlineElimination:
    """Sparse elimination of base logs with one retained target column."""

    def __init__(self, prime, bits):
        self.prime, self.bits = prime, bits
        self.pivots = {}
        self.rows = 0

    def add(self, x1, sign1, x2, sign2, target_coeff, rhs):
        r = self.prime
        self.rows += 1
        coeff = {}
        for x, sign in ((x1, sign1), (x2, sign2)):
            coeff[x] = (coeff.get(x, 0)+sign) % r
            if not coeff[x]:
                del coeff[x]
        target_coeff, rhs = target_coeff % r, rhs % r
        while coeff:
            pivot = min(coeff, key=self.bits)
            if pivot not in self.pivots:
                inv = pow(coeff[pivot], -1, r)
                self.pivots[pivot] = (
                    {x: c*inv % r for x, c in coeff.items()},
                    target_coeff*inv % r, rhs*inv % r)
                return None
            pivot_coeff, pivot_target, pivot_rhs = self.pivots[pivot]
            factor = coeff[pivot]
            for x, c in pivot_coeff.items():
                coeff[x] = (coeff.get(x, 0)-factor*c) % r
                if not coeff[x]:
                    del coeff[x]
            target_coeff = (target_coeff-factor*pivot_target) % r
            rhs = (rhs-factor*pivot_rhs) % r
        if target_coeff:
            return rhs*pow(target_coeff, -1, r) % r
        assert rhs == 0
        return None


class Case:
    def __init__(self, spec):
        self.spec = spec
        self.n, self.k = spec["n"], spec["k"]
        self.q, self.r = 1 << self.k, spec["prime"]
        ring = PolynomialRing(GF(2), "w")
        w = ring.gen()
        modulus = sum(w**e for e in spec["modulus"])
        self.K = GF(1 << self.n, "z", modulus=modulus)
        self.z = self.K.gen()
        self.F = self.K.subfield(self.k, "u")
        self.a2 = self.z**spec["a2_exp"]
        self.b6 = self.element(spec["b6_bits"])
        self.E = EllipticCurve(self.K, [1, self.a2, 0, 0, self.b6])
        y_ring = PolynomialRing(self.K, "Y")
        self.Y = y_ring.gen()
        assert is_prime(self.r) and self.E.cardinality() == 2*self.r
        assert (1/self.b6)**self.q != 1/self.b6

    def element(self, n):
        return sum(((n >> i) & 1)*self.z**i for i in range(self.n))

    @staticmethod
    def bits(x):
        return sum(int(c) << i for i, c in enumerate(x.polynomial().list()))

    def roots_at(self, x):
        poly = self.Y*self.Y+x*self.Y+x**3+self.a2*x*x+self.b6
        return sorted((y for y, mult in poly.roots()), key=self.bits)

    def generator(self):
        for n in range(1, 1 << self.n):
            x = self.element(n)
            ys = self.roots_at(x)
            if ys:
                G = 2*self.E(x, ys[0])
                if not G.is_zero():
                    assert self.r*G == self.E(0)
                    return G
        raise AssertionError("missing subgroup generator")

    def verify_base_coset(self):
        count = 0
        for f in self.F:
            x = self.K(f)
            if x == 0:
                continue
            ys = self.roots_at(x)
            if ys:
                P = self.E(x, ys[0])
                assert self.r*P != self.E(0)
                count += 1
        assert count == self.q//2
        return count

    def oracle(self, Q, cache, as_roots):
        if Q.is_zero():
            return []
        t, tq = Q[0], Q[0]**self.q
        assert t != tq
        A = self.F(t+tq)
        B = self.F(self.b6+self.b6**self.q)
        N = self.F(t*tq)
        C = B/A
        D = self.F(C*C+t*C+self.b6)
        alpha, beta = N/(A*A), D/(A*A)

        def lift(x):
            if x not in cache:
                ys = self.roots_at(x)
                cache[x] = self.E(x, ys[0]) if ys and x != 0 else None
            return cache[x]

        out = {}
        for v in as_roots.get(beta/(alpha*alpha), ()):
            u = alpha*v
            s, p = u.sqrt(), A*u+C
            xs = ([self.K(p.sqrt())]*2 if s == 0 else
                  sorted((self.K(s*w) for w in as_roots.get(p/(s*s), ())),
                         key=self.bits))
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
                assert P+R == Q
                assert (x1*x2+x1*t+x2*t)**2+x1*x2*t+self.b6 == 0
                out[(self.bits(x1), sign1, self.bits(x2), sign2)] = (
                    x1, sign1, x2, sign2)
        assert len(out) <= 2
        return list(out.values())

    def online(self, G, T, scalar, seed, cap=8192):
        rng = Random(seed)
        wall0, cpu0 = perf_counter(), process_time()
        as_roots = {}
        for v in self.F:
            as_roots.setdefault(v*v+v, (v, v+self.F.one()))
        assert len(as_roots) == self.q//2
        steps = []
        for _ in range(8):
            da, db = rng.randrange(self.r), rng.randrange(self.r)
            steps.append((da*G+db*T, da, db))
        a, b = rng.randrange(self.r), rng.randrange(self.r)
        Q = a*G+b*T
        cache, elim, hits = {}, OnlineElimination(self.r, self.bits), 0
        for query in range(1, cap+1):
            slot = 0 if Q.is_zero() else self.bits(Q[0]) % 8
            point, da, db = steps[slot]
            Q += point
            a, b = (a+da) % self.r, (b+db) % self.r
            relations = self.oracle(Q, cache, as_roots)
            hits += bool(relations)
            for x1, s1, x2, s2 in relations:
                answer = elim.add(x1, s1, x2, s2, -2*b, 2*a)
                if answer is not None:
                    assert answer == scalar and answer*G == T
                    return {
                        "seed": seed, "queries": query,
                        "failed_queries": query-hits,
                        "relation_rows": elim.rows,
                        "base_pivots": len(elim.pivots),
                        "lifted_abscissas": len(cache),
                        "recovered_scalar": answer,
                        "wall_seconds": perf_counter()-wall0,
                        "cpu_seconds": process_time()-cpu0,
                    }
        raise AssertionError("online target unresolved at query cap")

    def bsgs(self, G, T):
        m = isqrt(self.r)+1
        wall0, cpu0 = perf_counter(), process_time()
        table, Q = {}, self.E(0)
        for j in range(m):
            table.setdefault(Q, j)
            Q += G
        step, Q = -(m*G), T
        for i in range(m+1):
            if Q in table:
                scalar = (i*m+table[Q]) % self.r
                assert scalar*G == T
                return {"recovered_scalar": scalar, "baby_steps": m,
                        "giant_steps": i,
                        "wall_seconds": perf_counter()-wall0,
                        "cpu_seconds": process_time()-cpu0}
            Q += step
        raise AssertionError("BSGS failed")

    def rho(self, G, T, seed):
        rng = Random(seed)
        wall0, cpu0 = perf_counter(), process_time()

        def step(state):
            X, a, b = state
            slot = 0 if X.is_zero() else self.bits(X[0]) % 3
            if slot == 0:
                return X+G, (a+1) % self.r, b
            if slot == 1:
                return X+T, a, (b+1) % self.r
            return 2*X, 2*a % self.r, 2*b % self.r

        for restart in range(4):
            a, b = rng.randrange(self.r), rng.randrange(self.r)
            start = (a*G+b*T, a, b)
            slow, fast = step(start), step(step(start))
            for i in range(1, 100001):
                if slow[0] == fast[0]:
                    den = (slow[2]-fast[2]) % self.r
                    if den:
                        scalar = (fast[1]-slow[1])*pow(den, -1, self.r) % self.r
                        assert scalar*G == T
                        return {"recovered_scalar": scalar,
                                "restart": restart, "walk_steps": 3*i,
                                "wall_seconds": perf_counter()-wall0,
                                "cpu_seconds": process_time()-cpu0}
                    break
                slow, fast = step(slow), step(step(fast))
        raise AssertionError("rho failed at frozen cap")


def main():
    cases = []
    for spec in CASES:
        case = Case(spec)
        base_x = case.verify_base_coset()
        G = case.generator()
        rng = Random(SEED+case.r)
        scalars = [rng.randrange(1, case.r) for _ in range(3)]
        runs = []
        for i, scalar in enumerate(scalars):
            T = scalar*G
            runs.append({"scalar_for_verification": scalar,
                         "target": [case.bits(T[0]), case.bits(T[1])],
                         "bsgs": case.bsgs(G, T),
                         "rho": case.rho(G, T, SEED+i),
                         "online": case.online(G, T, scalar, SEED+i)})
        cases.append({"n": case.n, "k": case.k, "q": case.q,
                      "modulus": str(case.K.modulus()),
                      "a2_bits": case.bits(case.a2),
                      "b6_bits": spec["b6_bits"],
                      "curve_order": int(case.E.cardinality()),
                      "prime_order": case.r,
                      "base_abscissas": base_x,
                      "generator": [case.bits(G[0]), case.bits(G[1])],
                      "runs": runs})
    print(json.dumps({"seed": SEED, "sage_version": sage_version,
                      "python_platform": platform.platform(), "cases": cases},
                     indent=2, sort_keys=True, default=int))


if __name__ == "__main__":
    main()
