#!/usr/bin/env python3
"""Independent point/group and subgroup-log checks for the equal-67 panel."""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from csd.gf2_ghs import GF2n
from csd import ic_dreg as I

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "equal67-koblitz.json"
OUT = HERE / "equal67-verification.json"
BASIS = [1, 892, 1694, 3678, 5682, 8344, 22848]


def rank(rows, p):
    piv = {}
    for row in rows:
        v = [x % p for x in row]
        while True:
            j = next((i for i, x in enumerate(v) if x), None)
            if j is None:
                break
            if j not in piv:
                inverse = pow(v[j], -1, p)
                piv[j] = [(x*inverse) % p for x in v]
                break
            factor = v[j]
            v = [(x-factor*y) % p for x, y in zip(v, piv[j])]
    return len(piv)


def arm_check(E, xs, saved, targets, dlogs):
    reps = {x: min(p, E.neg(p)) for x in xs if (p := E.lift_x(x)) is not None}
    pts = [(p, x, 1) for x, p in reps.items()]
    pts += [(E.neg(p), x, -1) for x, p in reps.items() if E.neg(p) != p]
    cols = saved["columns"]
    assert len(pts) == 133 and len(cols) == 60
    col = {x: j for j, x in enumerate(cols)}
    projected_logs = {x: dlogs[E.mul(154, reps[x])] for x in cols}
    target_index = {p: j for j, p in enumerate(targets)}
    found = [set() for _ in targets]
    count = [0] * len(targets)
    for (p1, x1, s1), (p2, x2, s2) in itertools.combinations_with_replacement(pts, 2):
        j = target_index.get(E.add(p1, p2))
        if j is None:
            continue
        count[j] += 1 if p1 == p2 else 2
        row = [0] * (len(cols) + len(targets))
        if x1 in col:
            row[col[x1]] += s1
        if x2 in col:
            row[col[x2]] += s2
        row[len(cols)+j] -= 1
        found[j].add(tuple(x % 211 for x in row))
    all_rows = []
    for j, entry in enumerate(saved["targets"]):
        assert entry["scalar"] == j+1
        assert entry["checked_ordered_pairs"] == count[j]
        assert set(map(tuple, entry["rows"])) == found[j]
        for row in found[j]:
            value = sum(row[k]*projected_logs[x] for k, x in enumerate(cols))
            value += row[len(cols)+j] * ((154 * (j+1)) % 211)
            assert value % 211 == 0
        all_rows.extend(found[j])
        assert rank(all_rows, 211) == entry["rank_through_target"]
    assert rank(all_rows, 211) == saved["final_rank"]
    return {"admitted_targets": sum(bool(s) for s in found),
            "distinct_rows": len(set(all_rows)), "rank": rank(all_rows, 211),
            "checked_ordered_pairs": sum(count), "subgroup_log_rows_checked": len(all_rows)}


def sym_root_check(F, E, targets):
    T = (0, 1)
    out = []
    for j, target in enumerate(targets):
        uR = F.inv(target[0] ^ 1)
        wR = F.sq(uR) ^ uR
        hits = boundary = nonrational = 0
        for bits in range(1 << 13):
            c1, c2, eps = bits & 63, (bits >> 6) & 63, (bits >> 12) & 1
            v = []
            for c in (c1, c2):
                value = 0
                for k, b in enumerate(BASIS[1:]):
                    if c >> k & 1:
                        value ^= b
                v.append(value)
            w1, w2 = F.sq(v[0]) ^ v[0], F.sq(v[1]) ^ v[1]
            polynomial = F.mul(F.mul(w1, w2), wR) ^ w1 ^ w2 ^ wR ^ v[0] ^ v[1] ^ eps ^ uR
            if polynomial:
                continue
            hits += 1
            has_lift = False
            has_exception = False
            for e1 in (0, 1):
                u1, u2 = v[0] ^ e1, v[1] ^ e1 ^ eps
                if not u1 or not u2:
                    boundary += 1
                    has_exception = True
                    continue
                p1, p2 = E.lift_x(F.inv(u1) ^ 1), E.lift_x(F.inv(u2) ^ 1)
                if p1 is None or p2 is None:
                    nonrational += 1
                    has_exception = True
                    continue
                for a in (p1, E.neg(p1)):
                    for b in (p2, E.neg(p2)):
                        if E.add(a, b) in (target, E.add(target, T)):
                            has_lift = True
            assert has_lift or has_exception
        out.append({"scalar": j+1, "polynomial_roots": hits,
                    "boundary_choices": boundary, "nonrational_choices": nonrational})
    return out


def main():
    data = json.loads(SOURCE.read_text())
    F, E = GF2n(15, 0x8003), None
    E = I.BinaryCurve(F, 1, 1)
    Q, T = (15115, 1293), (0, 1)
    assert E.on_curve(Q) and E.on_curve(T) and E.mul(211, Q) is None
    assert E.add(T, T) is None
    # Count every rational abscissa, including the unique order-two lift.
    order = 1
    for x in range(1 << 15):
        p = E.lift_x(x)
        if p is not None:
            order += 1 if p == E.neg(p) else 2
    assert order == 32494
    V = I.span(BASIS)
    assert len(V) == 128 and all(F.sq(v) in V for v in V)
    xs = {"direct": data["preflight"]["bases"]["direct"]["abscissae"],
          "symmetrized": data["preflight"]["bases"]["symmetrized"]["abscissae"]}
    assert xs["direct"] == sorted(x for x in V if E.lift_x(x) is not None)
    assert xs["symmetrized"] == sorted(F.inv(u) ^ 1 for u in V if u and E.lift_x(F.inv(u) ^ 1) is not None)
    for x in xs["symmetrized"]:
        p = E.lift_x(x)
        translated = E.add(p, T)
        if x == 0:
            assert translated is None
        else:
            assert translated[0] == F.inv(x)
            assert F.inv(translated[0] ^ 1) == F.inv(x ^ 1) ^ 1
    dlogs = {E.mul(k, Q): k for k in range(211)}
    assert len(dlogs) == 211
    targets = [E.mul(k, Q) for k in range(1, 25)]
    assert data["preflight"]["targets"] == [list(t) for t in targets]
    arms = {name: arm_check(E, xs[name], data["relations"][name], targets, dlogs)
            for name in ("direct", "symmetrized")}
    roots = sym_root_check(F, E, targets)
    out = {"schema": "torsion-symm-equal67-verification/v1",
           "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
           "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           "curve_order": order, "torsion_translation_checked": len(xs["symmetrized"]),
           "arms": arms, "symmetrized_roots": roots}
    OUT.write_text(json.dumps(out, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"arms": arms, "root_counts": [x["polynomial_roots"] for x in roots]}, sort_keys=True))


if __name__ == "__main__":
    main()
