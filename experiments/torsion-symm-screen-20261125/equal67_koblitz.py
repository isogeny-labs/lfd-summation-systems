#!/usr/bin/env python3
"""Frozen equal-size K1/F2^15 x/u-frame relation and Boolean-system screen."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from csd.gf2_ghs import GF2n
from csd import ic_dreg as I

HERE = Path(__file__).resolve().parent
BASIS = [1, 892, 1694, 3678, 5682, 8344, 22848]
COFACTOR, PRIME, Q = 154, 211, (15115, 1293)


def sha_json(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def clock_start():
    return time.perf_counter(), time.process_time()


def elapsed(start):
    return {"wall_s": time.perf_counter() - start[0], "cpu_s": time.process_time() - start[1]}


def eval_equations(equations, mask):
    return all(sum((mask & m) == m for m in eq) % 2 == 0 for eq in equations)


def bp_to_equations(poly, n):
    return sorted({tuple(sorted(mask for mask, c in poly.items() if c >> k & 1))
                   for k in range(n)} - {()})


def sym_equations(F, target_x):
    """Native m=2 formula, with the same 13-bit layout as koblitz_symmetrised.rs."""
    uR = F.inv(target_x ^ 1)
    assert uR
    wR = F.sq(uR) ^ uR
    b = BASIS[1:]
    w_basis = [F.sq(v) ^ v for v in b]
    w1, w2 = I._bp_block(w_basis, 0), I._bp_block(w_basis, 6)
    s = I._bp_add(I._bp_add(I._bp_block(b, 0), I._bp_block(b, 6)),
                  I._bp_add(I._bp_const(uR), {1 << 12: 1}))
    product = I._bp_mul(F, I._bp_mul(F, w1, w2), I._bp_const(wR))
    poly = I._bp_add(I._bp_add(product, w1), I._bp_add(w2, I._bp_add(I._bp_const(wR), s)))
    return sym_eqs_from_poly(poly, F.N), uR, wR


def sym_eqs_from_poly(poly, n):
    return [list(eq) for eq in bp_to_equations(poly, n)]


def sym_value(F, c1, c2, epsilon, uR, wR):
    v1 = v2 = 0
    for j, b in enumerate(BASIS[1:]):
        if c1 >> j & 1:
            v1 ^= b
        if c2 >> j & 1:
            v2 ^= b
    w1, w2 = F.sq(v1) ^ v1, F.sq(v2) ^ v2
    return F.mul(F.mul(w1, w2), wR) ^ w1 ^ w2 ^ wR ^ v1 ^ v2 ^ epsilon ^ uR


def rank_mod(rows, p):
    pivots = {}
    for row in rows:
        v = [x % p for x in row]
        for j in sorted(pivots):
            if v[j]:
                a = v[j]
                v = [(x - a*y) % p for x, y in zip(v, pivots[j])]
        j = next((j for j, x in enumerate(v) if x), None)
        if j is not None:
            inv = pow(v[j], -1, p)
            pivots[j] = [(x*inv) % p for x in v]
    return len(pivots)


def preflight():
    start = clock_start()
    F, E = GF2n(15, 0x8003), None
    E = I.BinaryCurve(F, 1, 1)
    assert E.on_curve(Q) and E.mul(PRIME, Q) is None
    assert E.mul(COFACTOR, (2, 23020)) == Q
    V = I.span(BASIS)
    assert len(V) == 128 and len(set(V)) == 128
    assert all(F.sq(v) in V for v in V)
    assert 0 in V and 1 in V
    x_sets = {"direct": sorted(x for x in V if E.lift_x(x) is not None),
              "symmetrized": sorted(F.inv(u) ^ 1 for u in V if u and E.lift_x(F.inv(u) ^ 1) is not None)}
    assert all(len(s) == 67 for s in x_sets.values())
    assert len(set(x_sets["direct"]) & set(x_sets["symmetrized"])) == 7
    targets = [E.mul(i, Q) for i in range(1, 25)]
    assert len(set(targets)) == 24 and all(E.on_curve(t) for t in targets)
    receipt = {"basis": BASIS, "basis_sha256": sha_json(BASIS),
               "curve": {"n": 15, "modulus": 0x8003, "a": 1, "b": 1,
                         "order": 32494, "cofactor": COFACTOR, "prime_order": PRIME},
               "generator": list(Q), "targets": [list(t) for t in targets],
               "bases": {k: {"abscissae": v, "sha256": sha_json(v), "size": len(v)}
                         for k, v in x_sets.items()},
               "intersection_size": 7, "timing": elapsed(start)}
    return F, E, x_sets, targets, receipt


def arm_relations(E, xs, targets):
    start = clock_start()
    reps = {x: min(p, E.neg(p)) for x in xs if (p := E.lift_x(x)) is not None}
    columns = [x for x in xs if E.mul(COFACTOR, reps[x]) is not None]
    assert len(columns) == 60
    col_index = {x: j for j, x in enumerate(columns)}
    points = []
    for x in xs:
        p = reps[x]
        points.append((p, x, 1))
        if E.neg(p) != p:
            points.append((E.neg(p), x, -1))
    assert len(points) == 133
    target_index = {t: j for j, t in enumerate(targets)}
    per_target = [set() for _ in targets]
    checked_pairs = [0 for _ in targets]
    attempts = 0
    for p1, x1, s1 in points:
        for p2, x2, s2 in points:
            attempts += 1
            j = target_index.get(E.add(p1, p2))
            if j is None:
                continue
            checked_pairs[j] += 1
            row = [0] * (len(columns) + len(targets))
            if x1 in col_index:
                row[col_index[x1]] += s1
            if x2 in col_index:
                row[col_index[x2]] += s2
            row[len(columns) + j] -= 1
            per_target[j].add(tuple(x % PRIME for x in row))
    cumulative = []
    summary = []
    for i, rows in enumerate(per_target):
        cumulative.extend(sorted(rows))
        summary.append({"scalar": i+1, "checked_ordered_pairs": checked_pairs[i],
                        "distinct_rows": len(rows), "rank_through_target": rank_mod(cumulative, PRIME),
                        "rows": [list(r) for r in sorted(rows)]})
    return {"columns": columns, "columns_sha256": sha_json(columns),
            "point_count": len(points), "ordered_pair_attempts": attempts,
            "targets": summary, "admitted_targets": sum(bool(r) for r in per_target),
            "distinct_rows_total": len(set(cumulative)), "final_rank": rank_mod(cumulative, PRIME),
            "timing": elapsed(start)}


def systems(F, E, targets):
    out = []
    for i, R in enumerate(targets[:4], 1):
        start = clock_start()
        direct = sorted({tuple(eq) for eq in I.binary_descent(E, [BASIS, BASIS], R[0]) if eq})
        direct = [list(eq) for eq in direct]
        sym, uR, wR = sym_equations(F, R[0])
        direct_roots = []
        sym_roots = []
        xs = I.span(BASIS)
        for bits in range(1 << 14):
            x1, x2 = xs[bits & 127], xs[bits >> 7]
            zero = I.s3_binary_value(F, 1, x1, x2, R[0]) == 0
            assert eval_equations(direct, bits) == zero
            if zero:
                direct_roots.append(bits)
        for bits in range(1 << 13):
            zero = sym_value(F, bits & 63, (bits >> 6) & 63, (bits >> 12) & 1, uR, wR) == 0
            assert eval_equations(sym, bits) == zero
            if zero:
                sym_roots.append(bits)
        out.append({"scalar": i, "target": list(R), "presentations": {
            "direct": {"variables": 14, "equations": direct, "roots": direct_roots,
                       "input_degree": max((m.bit_count() for eq in direct for m in eq), default=0),
                       "equations_sha256": sha_json(direct), "roots_sha256": sha_json(direct_roots)},
            "symmetrized": {"variables": 13, "equations": sym, "roots": sym_roots,
                            "input_degree": max((m.bit_count() for eq in sym for m in eq), default=0),
                            "equations_sha256": sha_json(sym), "roots_sha256": sha_json(sym_roots)}},
            "timing": elapsed(start)})
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=HERE / "equal67-koblitz.json")
    args = parser.parse_args()
    F, E, bases, targets, receipt = preflight()
    result = {"schema": "torsion-symm-equal67/v1",
              "protocol": "PROTOCOL_EQUAL67_20261125.md",
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "preflight": receipt, "relations": {}, "systems": systems(F, E, targets)}
    for name in ("direct", "symmetrized"):
        result["relations"][name] = arm_relations(E, bases[name], targets)
        print(name, "admitted", result["relations"][name]["admitted_targets"],
              "rank", result["relations"][name]["final_rank"], flush=True)
    args.out.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    print(args.out, flush=True)


if __name__ == "__main__":
    main()
