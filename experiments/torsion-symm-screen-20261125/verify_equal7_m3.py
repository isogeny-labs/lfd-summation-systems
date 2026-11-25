#!/usr/bin/env python3
"""Independent group-law and rank audit of the equal-7 three-summand panel."""
from __future__ import annotations

from collections import Counter
import hashlib
from itertools import combinations_with_replacement, product
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from csd.gf2_ghs import GF2n
from csd import ic_dreg as I
from verify_equal67 import rank

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "equal7-m3.json"
OUT = HERE / "equal7-m3-verification.json"
BASIS = [1, 2, 4]


def roots_with_lifts(F, E, system, target):
    V = I.span(BASIS)
    direct_valid = 0
    for bits in system["presentations"]["direct"]["roots"]:
        pts = [E.lift_x(V[(bits >> (3*i)) & 7]) for i in range(3)]
        if any(p is None for p in pts):
            continue
        if any(E.add(E.add(a, b), c) in (target, E.neg(target))
               for a in (pts[0], E.neg(pts[0]))
               for b in (pts[1], E.neg(pts[1]))
               for c in (pts[2], E.neg(pts[2]))):
            direct_valid += 1
    sym_valid = 0
    for bits in system["presentations"]["symmetrized"]["roots"]:
        vals = []
        for i in range(3):
            code = (bits >> (2*i)) & 3
            vals.append((BASIS[1] if code & 1 else 0) ^
                        (BASIS[2] if code & 2 else 0))
        parity = (bits >> 6) & 1
        valid = False
        for e1, e2 in product((0, 1), repeat=2):
            us = [vals[0]^e1, vals[1]^e2, vals[2]^parity^e1^e2]
            if not all(us):
                continue
            pts = [E.lift_x(F.inv(u)^1) for u in us]
            if any(p is None for p in pts):
                continue
            valid |= any(E.add(E.add(a, b), c) in (target, E.add(target, (0, 1)))
                         for a in (pts[0], E.neg(pts[0]))
                         for b in (pts[1], E.neg(pts[1]))
                         for c in (pts[2], E.neg(pts[2])))
        if valid:
            sym_valid += 1
    return {"scalar": system["scalar"],
            "direct_roots": len(system["presentations"]["direct"]["roots"]),
            "direct_group_lifts": direct_valid,
            "symmetrized_roots": len(system["presentations"]["symmetrized"]["roots"]),
            "symmetrized_group_lifts": sym_valid}


def arm_rows(E, xs, saved, targets, logs):
    reps = {x: min(p, E.neg(p)) for x in xs if (p := E.lift_x(x)) is not None}
    signed = [(p, x, 1) for x, p in reps.items()]
    signed += [(E.neg(p), x, -1) for x, p in reps.items() if E.neg(p) != p]
    assert len(signed) == 13
    cols = saved["columns"]
    assert len(cols) == 6
    col = {x: k for k, x in enumerate(cols)}
    target_index = {p: k for k, p in enumerate(targets)}
    rows = [set() for _ in targets]
    counts = [0] * len(targets)
    for triple in combinations_with_replacement(range(len(signed)), 3):
        pp = [signed[k] for k in triple]
        j = target_index.get(E.add(E.add(pp[0][0], pp[1][0]), pp[2][0]))
        if j is None:
            continue
        multiplicities = Counter(triple)
        counts[j] += 6 // (2 if 2 in multiplicities.values() else 6 if 3 in multiplicities.values() else 1)
        row = [0] * len(cols)
        for _, x, sign in pp:
            if x in col:
                row[col[x]] += sign
        rows[j].add(tuple(v % 37 for v in row))
    running = []
    log_cols = [logs[E.mul(14, reps[x])] for x in cols]
    for j, record in enumerate(saved["targets"]):
        assert counts[j] == record["checked_ordered_triples"]
        assert rows[j] == set(map(tuple, record["rows"]))
        running.extend(rows[j])
        assert rank(running, 37) == record["base_rank_through_target"]
        for row in rows[j]:
            assert sum(a*b for a, b in zip(row, log_cols)) % 37 == 14*(j+1) % 37
    assert rank(running, 37) == saved["base_rank"]
    return {"admitted_targets": sum(bool(r) for r in rows),
            "distinct_rows": len(set(running)), "base_rank": rank(running, 37),
            "checked_ordered_triples": sum(counts), "subgroup_log_rows_checked": len(running)}


def main():
    data = json.loads(SOURCE.read_text())
    F = GF2n(9, 515)
    E = I.BinaryCurve(F, 1, 1)
    Q = (135, 123)
    assert E.mul(14, (2, 47)) == Q and E.mul(37, Q) is None
    order = 1
    for x in range(512):
        p = E.lift_x(x)
        if p is not None:
            order += 1 if p == E.neg(p) else 2
    assert order == 518
    selected = None
    for a in range(2, 512):
        for b in range(a+1, 512):
            candidate = [1, a, b]
            span = I.span(candidate)
            if len(set(span)) != 8:
                continue
            direct = {x for x in span if E.lift_x(x) is not None}
            if len(direct) < 4:
                continue
            sym = {F.inv(u)^1 for u in span if u and E.lift_x(F.inv(u)^1) is not None}
            if len(direct) != len(sym) or direct == sym:
                continue
            if min(sum(E.mul(14, E.lift_x(x)) is not None for x in direct),
                   sum(E.mul(14, E.lift_x(x)) is not None for x in sym)) < 3:
                continue
            selected = candidate
            break
        if selected is not None:
            break
    assert selected == BASIS
    V = I.span(BASIS)
    xs = {"direct": data["preflight"]["bases"]["direct"]["abscissae"],
          "symmetrized": data["preflight"]["bases"]["symmetrized"]["abscissae"]}
    assert xs["direct"] == sorted(x for x in V if E.lift_x(x) is not None)
    assert xs["symmetrized"] == sorted(F.inv(u)^1 for u in V if u and E.lift_x(F.inv(u)^1) is not None)
    for x in xs["symmetrized"]:
        p = E.lift_x(x)
        if x:
            q = E.add(p, (0, 1))
            assert q[0] == F.inv(x)
            assert F.inv(q[0]^1) == F.inv(x^1)^1
    targets = [E.mul(k, Q) for k in range(1, 25)]
    assert [list(p) for p in targets] == data["preflight"]["targets"]
    logs = {E.mul(k, Q): k for k in range(37)}
    assert len(logs) == 37
    arms = {a: arm_rows(E, xs[a], data["relations"][a], targets, logs)
            for a in ("direct", "symmetrized")}
    root_checks = [roots_with_lifts(F, E, s, tuple(s["target"])) for s in data["systems"]]
    result = {"schema": "torsion-symm-equal7-m3-verification/v1",
              "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "curve_order": order, "base_selection_replayed": selected,
              "arms": arms, "first_four_root_lifts": root_checks}
    OUT.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"arms": arms, "roots": root_checks}, sort_keys=True))


if __name__ == "__main__":
    main()
