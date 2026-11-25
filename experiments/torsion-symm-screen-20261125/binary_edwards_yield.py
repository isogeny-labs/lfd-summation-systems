#!/usr/bin/env python3
"""Exact small-field relation-yield screen for the published binary Edwards t map."""
from __future__ import annotations

from collections import defaultdict
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from csd.ic_dreg import binary_field, binary_modulus  # noqa: E402
from csd.field import factor, is_prime  # noqa: E402

HERE = Path(__file__).resolve().parent
CELLS = ((5, 2), (7, 3), (9, 4))
O = (0, 0)
T4 = (1, 0)


class Edwards:
    def __init__(self, n):
        self.n = n
        self.F = binary_field(n)
        self.q = 1 << n
        self.mul = self.F.mul
        self.sq = self.F.sq

    def on_curve(self, P):
        x, y = P
        m, s = self.mul, self.sq
        lhs = x ^ y ^ s(x) ^ s(y)
        rhs = m(x, y) ^ m(m(x, y), x ^ y) ^ m(s(x), s(y))
        return lhs == rhs

    def neg(self, P):
        return P[1], P[0]

    def add(self, P, Q):
        x1, y1 = P
        x2, y2 = Q
        m, s = self.mul, self.sq
        u, v = x1 ^ s(x1), y1 ^ s(y1)
        xy = m(x1 ^ y1, x2 ^ y2)
        nx = x1 ^ x2 ^ xy ^ m(u, m(x2, y1 ^ y2 ^ 1) ^ m(y1, y2))
        ny = y1 ^ y2 ^ xy ^ m(v, m(y2, x1 ^ x2 ^ 1) ^ m(x1, x2))
        dx = 1 ^ m(u, x2 ^ y2)
        dy = 1 ^ m(v, x2 ^ y2)
        if dx == 0 or dy == 0:
            raise ArithmeticError("complete Edwards addition denominator vanished")
        return m(nx, self.F.inv(dx)), m(ny, self.F.inv(dy))

    def scalar(self, k, P):
        R = O
        while k:
            if k & 1:
                R = self.add(R, P)
            P = self.add(P, P)
            k >>= 1
        return R

    def s3(self, t1, t2, tr):
        m, s = self.mul, self.sq
        a = m(t1, t2)
        e = t1 ^ t2
        inv = m(s(t1) ^ t1, s(t2) ^ t2)
        direct = m(s(a) ^ m(a, e) ^ a ^ 1, s(tr)) \
            ^ m(s(a) ^ m(a, e) ^ a, tr) ^ s(e)
        quotient = s(e) ^ m(inv, s(tr) ^ tr) ^ s(tr)
        if direct != quotient:
            raise AssertionError("direct/invariant f3 mismatch")
        return direct

    def points(self):
        return [(x, y) for x in range(self.q) for y in range(self.q)
                if self.on_curve((x, y))]


def subspace(j):
    return tuple(range(1 << j))


def rank_mod(rows, p):
    basis = {}
    for vector in rows:
        row = list(vector)
        for pivot in sorted(basis, reverse=True):
            if row[pivot]:
                coefficient = row[pivot]
                row = [(x - coefficient*y) % p for x, y in zip(row, basis[pivot])]
        pivot = next((i for i in range(len(row)-1, -1, -1) if row[i]), None)
        if pivot is not None:
            inv = pow(row[pivot], -1, p)
            basis[pivot] = [(x*inv) % p for x in row]
    return len(basis)


def cell(n, l):
    started = time.perf_counter()
    E = Edwards(n)
    points = E.points()
    point_set = set(points)
    if O not in point_set or T4 not in point_set:
        raise AssertionError("rational torsion points absent")
    if E.scalar(2, T4) != (1, 1) or E.scalar(4, T4) != O:
        raise AssertionError("T4 does not have order four")
    # Check the published complete law, inverse and t-translation on all points.
    for P in points:
        if E.add(P, O) != P or E.add(P, E.neg(P)) != O:
            raise AssertionError("identity/inverse failure")
        Q = E.add(P, T4)
        if Q not in point_set or (Q[0] ^ Q[1]) != (P[0] ^ P[1] ^ 1):
            raise AssertionError("t(P+T4)=t(P)+1 failure")
    for P in points:
        t1 = P[0] ^ P[1]
        for Q in points:
            t2 = Q[0] ^ Q[1]
            S = E.add(P, Q)
            if E.s3(t1, t2, S[0] ^ S[1]) != 0:
                raise AssertionError("published f3 fails on a checked point sum")
    order = len(points)
    odd_primes = [p for p in factor(order) if p & 1 and is_prime(p)]
    if not odd_primes:
        raise AssertionError("no odd prime subgroup")
    r = max(odd_primes)
    cofactor = order // r
    generator = next((E.scalar(cofactor, P) for P in points
                      if E.scalar(cofactor, P) != O), None)
    if generator is None or E.scalar(r, generator) != O:
        raise AssertionError("prime-order generator selection failed")
    target_count = min(24, r - 1)
    targets = [(i, E.scalar(i, generator)) for i in range(1, target_count+1)]
    # All t orbits are kept intact; exclude the torsion-only orbit {0,1}.
    W = set(subspace(l+1))
    A = set(subspace(l))
    B = {x ^ (1 << l) for x in A}
    if A & B or A | B != W:
        raise AssertionError("colored bases do not partition matched union")
    fibers = defaultdict(list)
    for P in points:
        fibers[P[0] ^ P[1]].append(P)
    columns = []
    for t in sorted(W):
        if t > (t ^ 1) or t in (0, 1):
            continue
        orbit_points = fibers[t] + fibers[t ^ 1]
        if not orbit_points:
            continue
        rep = min(orbit_points)
        if E.scalar(cofactor, rep) == O:
            continue
        columns.append({"t_pair": [t, t ^ 1], "point": rep,
                        "color": "A" if t in A else "B"})
    if not columns:
        raise AssertionError("empty matched base")
    # Independent S3 check on every factor-base coordinate pair and target.
    for c1 in columns:
        for c2 in columns:
            for _, R in targets:
                E.s3(c1["t_pair"][0], c2["t_pair"][0], R[0] ^ R[1])
    torsion = [E.scalar(j, T4) for j in range(4)]
    signed = [(P, E.neg(P)) for P in (c["point"] for c in columns)]
    rows = {"baseline": [], "colored": []}
    admitted = {"baseline": 0, "colored": 0}
    per_target = []
    for scalar, R in targets:
        found = {"baseline": set(), "colored": set()}
        for i, c1 in enumerate(columns):
            for j in range(i, len(columns)):
                c2 = columns[j]
                color_cross = c1["color"] != c2["color"]
                for si, P in enumerate(signed[i]):
                    for sj, Q in enumerate(signed[j]):
                        S = E.add(P, Q)
                        for z, T in enumerate(torsion):
                            if E.add(S, T) != R:
                                continue
                            t1, t2, tr = (c1["t_pair"][0], c2["t_pair"][0], R[0] ^ R[1])
                            if not any(E.s3(t1 ^ a, t2 ^ b, tr) == 0
                                       for a in (0, 1) for b in (0, 1)):
                                raise AssertionError("checked relation lacks a summation-polynomial lift")
                            vec = [0] * (len(columns)+1)
                            vec[i] = (vec[i] + (1 if si == 0 else -1)) % r
                            vec[j] = (vec[j] + (1 if sj == 0 else -1)) % r
                            vec[-1] = (-scalar) % r
                            witness = tuple(vec)
                            found["baseline"].add(witness)
                            if color_cross:
                                found["colored"].add(witness)
        if not found["colored"] <= found["baseline"]:
            raise AssertionError("colored rows not contained in matched union")
        for arm in rows:
            if found[arm]:
                admitted[arm] += 1
            rows[arm].extend(sorted(found[arm]))
        per_target.append({"scalar": scalar, "target": R,
                           "baseline_rows": len(found["baseline"]),
                           "colored_rows": len(found["colored"]),
                           "baseline_rank": rank_mod(rows["baseline"], r),
                           "colored_rank": rank_mod(rows["colored"], r)})
    return {"n": n, "l": l, "modulus": binary_modulus(n), "curve": "E_(1,1)",
            "group_order": order, "r": r, "cofactor": cofactor,
            "generator": generator, "torsion_point": T4,
            "total_columns": len(columns), "columns": columns,
            "target_count": target_count, "admitted_targets": admitted,
            "distinct_rows": {arm: len(set(v)) for arm, v in rows.items()},
            "relation_rows": {arm: sorted(set(v)) for arm, v in rows.items()},
            "final_rank": {arm: rank_mod(v, r) for arm, v in rows.items()},
            "per_target": per_target, "elapsed_s": time.perf_counter()-started}


def main():
    source = Path(__file__).read_bytes()
    output = {"schema": "torsion-symm-binary-yield/v1",
              "protocol": "PROTOCOL_BINARY_20261125.md",
              "source_sha256": hashlib.sha256(source).hexdigest(),
              "cells": [cell(n, l) for n, l in CELLS]}
    path = HERE / "binary-edwards-yield.json"
    path.write_text(json.dumps(output, sort_keys=True, indent=2)+"\n")
    for item in output["cells"]:
        print(item["n"], item["group_order"], item["r"], item["total_columns"],
              item["admitted_targets"], item["final_rank"], round(item["elapsed_s"], 3))


if __name__ == "__main__":
    main()
