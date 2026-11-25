#!/usr/bin/env python3
"""Matched frozen binary S3 and boundary-patched two-torsion invariant screen."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from csd.gf2_ghs import GF2n  # noqa: E402
from csd import ic_dreg as I  # noqa: E402
from csd.ic_local_certificate import local_profile  # noqa: E402

HERE = Path(__file__).resolve().parent
CASES = ROOT / "experiments/ic-curve-sweep/cases.jsonl"
IDS = ("1bac99e1e73e08b00999", "a5eb99758876c7f2fee6", "7c65f50a88d2c5e980b2")
R = 37  # #E(F_512)=518=14*37, checked independently below.


def rank_mod(rows, p):
    pivots = {}
    for row in rows:
        v = [int(x) % p for x in row]
        for j in sorted(pivots):
            if v[j]:
                a = v[j]
                v = [(u - a * q) % p for u, q in zip(v, pivots[j])]
        j = next((j for j, x in enumerate(v) if x), None)
        if j is not None:
            inverse = pow(v[j], -1, p)
            pivots[j] = [(x * inverse) % p for x in v]
    return len(pivots)


def anf(values, variables):
    coeff = list(values)
    for bit in range(variables):
        for mask in range(1 << variables):
            if mask & (1 << bit):
                coeff[mask] ^= coeff[mask ^ (1 << bit)]
    return [mask for mask, x in enumerate(coeff) if x]


def invariant_equations(F, base, z):
    values = []
    boundary = []
    for bits in range(64):
        x1, x2 = base[bits & 7], base[bits >> 3]
        u = x1 ^ x2
        v = F.mul(x1, x2)
        direct = I.s3_binary_value(F, 1, x1, x2, z)
        if v:
            inv = F.inv(v)
            s = v ^ inv
            w = F.mul(F.sq(u), inv)
            quotient = s ^ F.mul(F.sq(z), w) ^ z
            assert direct == F.mul(v, quotient)
            value = quotient
        else:
            boundary.append((bits, x1, x2, direct))
            value = direct
        assert (value == 0) == (direct == 0)
        values.append(value)
    equations = sorted(set(tuple(anf([(v >> j) & 1 for v in values], 6))
                           for j in range(F.N)) - {()})
    for bits, value in enumerate(values):
        assert all(sum((bits & monom) == monom for monom in eq) % 2 == 0 for eq in equations) == (value == 0)
    return equations, [bits for bits, value in enumerate(values) if value == 0], boundary


def relation_rows(E, base, targets):
    points = {}
    for x in base:
        point = E.lift_x(x)
        if point is not None:
            points[x] = min(point, E.neg(point))
    cols = [x for x in base if x and x in points]
    all_rows = []
    records = []
    for target_index, target in enumerate(targets):
        checked = []
        rows = []
        for x1 in base:
            if x1 not in points:
                continue
            for x2 in base:
                if x2 not in points:
                    continue
                signs1 = (1,) if points[x1] == E.neg(points[x1]) else (1, -1)
                signs2 = (1,) if points[x2] == E.neg(points[x2]) else (1, -1)
                for sign1 in signs1:
                    p1 = points[x1] if sign1 == 1 else E.neg(points[x1])
                    for sign2 in signs2:
                        p2 = points[x2] if sign2 == 1 else E.neg(points[x2])
                        if E.add(p1, p2) != target:
                            continue
                        row = [0] * (len(cols) + len(targets))
                        for x, sign in ((x1, sign1), (x2, sign2)):
                            if x:
                                row[cols.index(x)] += sign
                        row[len(cols) + target_index] -= 1
                        row = [x % R for x in row]
                        checked.append((x1, x2, sign1, sign2))
                        rows.append(row)
        unique = sorted({tuple(row) for row in rows})
        all_rows.extend(unique)
        records.append({"target": list(target), "checked_ordered_signed": checked,
                        "distinct_rows": [list(row) for row in unique],
                        "rank_through_target": rank_mod(all_rows, R)})
    return {"columns": cols, "point_representatives": {str(x): list(p) for x, p in points.items()},
            "targets": records, "total_distinct_rows": len(set(map(tuple, all_rows))),
            "final_rank": rank_mod(all_rows, R)}


def main():
    rows = {row["case_id"]: row for line in CASES.open() if (row := json.loads(line)).get("case_id") in IDS}
    assert set(rows) == set(IDS)
    first = rows[IDS[0]]
    assert all((row["n"], row["modulus"], row["a"], row["b"], row["basis"], row["base_kind"])
               == (9, 515, 1, 1, [302, 252, 1], "frobenius-invariant") for row in rows.values())
    F = GF2n(9, 515)
    E = I.BinaryCurve(F, 1, 1)
    base = I.span(first["basis"])
    assert len(base) == 8 and set(base) == {0, 1, 252, 253, 302, 303, 466, 467}
    assert all(F.inv(x) in base for x in base if x)
    group_points = {None}
    for x in range(512):
        point = E.lift_x(x)
        if point is not None:
            group_points.update((point, E.neg(point)))
    assert len(group_points) == 518

    output = {"schema": "torsion-symm-frozen-boolean/v1",
              "protocol": "PROTOCOL_FROZEN_BOOLEAN_20261125.md",
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "cases_sha256": hashlib.sha256(CASES.read_bytes()).hexdigest(),
              "base": base, "curve_order": 518, "prime_order": R, "cases": []}
    for case_id in IDS:
        case = rows[case_id]
        target = tuple(case["target"])
        assert E.on_curve(target)
        start = time.perf_counter(); cpu_start = time.process_time()
        direct = sorted(set(tuple(eq) for eq in I.binary_descent(E, [first["basis"]] * 2, target[0]) if eq))
        invariant, inv_points, boundary = invariant_equations(F, base, target[0])
        points = [bits for bits in range(64) if I.s3_binary_value(F, 1, base[bits & 7], base[bits >> 3], target[0]) == 0]
        assert points == inv_points
        assert all(all(sum((bits & m) == m for m in eq) % 2 == 0 for eq in direct) == (bits in points)
                   for bits in range(64))
        setup = {"wall_s": time.perf_counter()-start, "cpu_s": time.process_time()-cpu_start}
        presentations = {}
        for name, equations in (("direct", direct), ("patched_invariant", invariant)):
            start = time.perf_counter(); cpu_start = time.process_time()
            result = local_profile(equations, 6, points, max_degree=8, max_columns=64)
            presentations[name] = {"equations": [list(eq) for eq in equations],
                                   "generator_sha256": hashlib.sha256(json.dumps(equations, separators=(",", ":")).encode()).hexdigest(),
                                   "input_degree": max((m.bit_count() for eq in equations for m in eq), default=0),
                                   "closure": result, "closure_wall_s": time.perf_counter()-start,
                                   "closure_cpu_s": time.process_time()-cpu_start}
        output["cases"].append({"case_id": case_id, "target_class": case["target_class"],
                                "target": list(target), "solution_masks": points,
                                "boundary_pairs": boundary, "setup": setup,
                                "presentations": presentations})
        print(case_id, "solutions", len(points), "LFD", presentations["direct"]["closure"]["exact_lfd"],
              presentations["patched_invariant"]["closure"]["exact_lfd"], flush=True)
    start = time.perf_counter(); cpu_start = time.process_time()
    output["relations"] = relation_rows(E, base, [tuple(rows[x]["target"]) for x in IDS])
    output["relation_wall_s"] = time.perf_counter()-start
    output["relation_cpu_s"] = time.process_time()-cpu_start
    (HERE / "frozen-boolean-torsion.json").write_text(json.dumps(output, sort_keys=True, indent=2)+"\n")
    print("relation rank", output["relations"]["final_rank"])


if __name__ == "__main__":
    main()
