#!/usr/bin/env python3
"""Frozen equal-11 three-summand Koblitz x/u comparison."""
from __future__ import annotations

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
BASIS = [1, 2, 4, 56]
Q, T = (135, 123), (0, 1)
COFACTOR, PRIME = 14, 37
TERMS = [
    (2,2,2,0,0), (2,2,0,2,0), (2,0,2,2,0), (0,2,2,2,0),
    (2,1,1,1,0), (1,2,1,1,0), (1,1,2,1,0), (1,1,1,2,0),
    (1,1,1,1,1), (2,0,0,0,0), (0,2,0,0,0), (0,0,2,0,0),
    (0,0,0,2,0), (1,0,0,0,0), (0,1,0,0,0), (0,0,1,0,0),
    (0,0,0,1,0), (0,0,0,0,1),
]


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def start():
    return time.perf_counter(), time.process_time()


def elapsed(t):
    return {"wall_s": time.perf_counter()-t[0], "cpu_s": time.process_time()-t[1]}


def rank(rows):
    pivots = {}
    for row in rows:
        v = [x % PRIME for x in row]
        for j in sorted(pivots):
            if v[j]:
                f = v[j]
                v = [(x-f*y) % PRIME for x, y in zip(v, pivots[j])]
        j = next((j for j, x in enumerate(v) if x), None)
        if j is not None:
            inv = pow(v[j], -1, PRIME)
            pivots[j] = [(x*inv) % PRIME for x in v]
    return len(pivots)


def poly_equations(poly, n):
    return [list(eq) for eq in sorted({tuple(sorted(mask for mask, c in poly.items() if c >> k & 1))
                                       for k in range(n)} - {()})]


def sym_equations(F, xR):
    uR = F.inv(xR ^ 1)
    assert uR
    wR = F.sq(uR) ^ uR
    w_basis = [F.sq(b) ^ b for b in BASIS[1:]]
    w = [I._bp_block(w_basis, 3*i) for i in range(3)]
    s = I._bp_const(uR)
    for i in range(3):
        s = I._bp_add(s, I._bp_block(BASIS[1:], 3*i))
    s = I._bp_add(s, {1 << 9: 1})
    factors = w + [I._bp_const(wR), s]
    result = {}
    for exponents in TERMS:
        term = I._bp_const(1)
        for factor, exponent in zip(factors, exponents):
            if exponent:
                power = factor if exponent == 1 else I._bp_sq(F, factor)
                term = I._bp_mul(F, term, power)
        result = I._bp_add(result, term)
    return poly_equations(result, F.N), uR, wR


def sym_value(F, bits, uR, wR):
    vals = []
    for i in range(3):
        c = (bits >> (3*i)) & 7
        vals.append((BASIS[1] if c & 1 else 0) ^
                    (BASIS[2] if c & 2 else 0) ^
                    (BASIS[3] if c & 4 else 0))
    w = [F.sq(v) ^ v for v in vals] + [wR]
    s = vals[0] ^ vals[1] ^ vals[2] ^ ((bits >> 9) & 1) ^ uR
    value = 0
    for exponents in TERMS:
        term = 1
        for factor, exponent in zip(w+[s], exponents):
            for _ in range(exponent):
                term = F.mul(term, factor)
        value ^= term
    return value


def eval_eq(equations, bits):
    return all(sum((bits & m) == m for m in eq) % 2 == 0 for eq in equations)


def preflight():
    t = start()
    F = GF2n(9, 515)
    E = I.BinaryCurve(F, 1, 1)
    assert E.mul(14, (2, 47)) == Q and E.mul(PRIME, Q) is None
    V = I.span(BASIS)
    assert len(set(V)) == 16
    xs = {"direct": sorted(x for x in V if E.lift_x(x) is not None),
          "symmetrized": sorted(F.inv(u) ^ 1 for u in V if u and E.lift_x(F.inv(u) ^ 1) is not None)}
    assert all(len(x) == 11 for x in xs.values())
    assert len(set(xs["direct"]) & set(xs["symmetrized"])) == 1
    targets = [E.mul(i, Q) for i in range(1, 25)]
    return F, E, xs, targets, {"curve": {"n": 9, "modulus": 515, "a": 1, "b": 1,
                                         "order": 518, "cofactor": 14, "prime_order": 37},
                               "basis": BASIS, "generator": list(Q),
                               "targets": [list(x) for x in targets],
                               "bases": {name: {"abscissae": v, "sha256": digest(v), "size": len(v)}
                                         for name, v in xs.items()},
                               "timing": elapsed(t)}


def systems(F, E, targets):
    out = []
    values = I.span(BASIS)
    for i, R in enumerate(targets[:4], 1):
        t = start()
        direct = sorted({tuple(eq) for eq in I.binary_descent(E, [BASIS]*3, R[0]) if eq})
        direct = [list(eq) for eq in direct]
        sym, uR, wR = sym_equations(F, R[0])
        dr, sr = [], []
        for bits in range(1 << 12):
            x1, x2, x3 = values[bits & 15], values[(bits >> 4) & 15], values[bits >> 8]
            zero = I.s4_binary_value(F, 1, x1, x2, x3, R[0]) == 0
            assert eval_eq(direct, bits) == zero
            if zero:
                dr.append(bits)
        for bits in range(1 << 10):
            zero = sym_value(F, bits, uR, wR) == 0
            assert eval_eq(sym, bits) == zero
            if zero:
                sr.append(bits)
        out.append({"scalar": i, "target": list(R), "presentations": {
            "direct": {"variables": 12, "equations": direct, "roots": dr,
                       "input_degree": max((m.bit_count() for eq in direct for m in eq), default=0),
                       "equations_sha256": digest(direct), "roots_sha256": digest(dr)},
            "symmetrized": {"variables": 10, "equations": sym, "roots": sr,
                            "input_degree": max((m.bit_count() for eq in sym for m in eq), default=0),
                            "equations_sha256": digest(sym), "roots_sha256": digest(sr)}},
                    "timing": elapsed(t)})
    return out


def relations(E, xs, targets):
    t = start()
    reps = {x: min(p, E.neg(p)) for x in xs if (p := E.lift_x(x)) is not None}
    cols = [x for x in xs if E.mul(COFACTOR, reps[x]) is not None]
    assert len(cols) == 10
    col = {x: j for j, x in enumerate(cols)}
    points = []
    for x, p in reps.items():
        points.append((p, x, 1))
        if E.neg(p) != p:
            points.append((E.neg(p), x, -1))
    assert len(points) == 21
    by_target = {R: j for j, R in enumerate(targets)}
    rows = [set() for _ in targets]
    matches = [0] * len(targets)
    for p1, x1, s1 in points:
        for p2, x2, s2 in points:
            p12 = E.add(p1, p2)
            for p3, x3, s3 in points:
                j = by_target.get(E.add(p12, p3))
                if j is None:
                    continue
                matches[j] += 1
                row = [0] * len(cols)
                for x, sign in ((x1, s1), (x2, s2), (x3, s3)):
                    if x in col:
                        row[col[x]] += sign
                rows[j].add(tuple(x % PRIME for x in row))
    all_rows = []
    summary = []
    for j, r in enumerate(rows):
        all_rows.extend(sorted(r))
        summary.append({"scalar": j+1, "checked_ordered_triples": matches[j],
                        "distinct_rows": len(r), "base_rank_through_target": rank(all_rows),
                        "rows": [list(x) for x in sorted(r)]})
    return {"columns": cols, "columns_sha256": digest(cols), "signed_points": len(points),
            "ordered_triple_attempts": len(points)**3,
            "targets": summary, "admitted_targets": sum(bool(r) for r in rows),
            "distinct_rows_total": len(set(all_rows)), "base_rank": rank(all_rows),
            "timing": elapsed(t)}


def main():
    F, E, xs, targets, pre = preflight()
    out = {"schema": "torsion-symm-equal11-m3/v1",
           "protocol": "PROTOCOL_EQUAL11_M3_20261125.md",
           "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           "native_sym_terms": TERMS, "native_sym_terms_sha256": digest(TERMS),
           "preflight": pre, "systems": systems(F, E, targets),
           "relations": {a: relations(E, xs[a], targets) for a in ("direct", "symmetrized")}}
    path = HERE / "equal11-m3.json"
    path.write_text(json.dumps(out, sort_keys=True, indent=2) + "\n")
    print(json.dumps({a: {"admitted": r["admitted_targets"], "rank": r["base_rank"]}
                      for a, r in out["relations"].items()}, sort_keys=True))
    print(path)


if __name__ == "__main__":
    main()
