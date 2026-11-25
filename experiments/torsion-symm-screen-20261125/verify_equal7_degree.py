#!/usr/bin/env python3
"""Independent bitset closure at the decisive equal-7 target's LFD boundary."""
from __future__ import annotations

from itertools import combinations
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "equal7-m3.json"
OUT = HERE / "equal7-degree-verification.json"


def masks_through(variables, degree):
    return sorted((sum(1 << j for j in combo)
                   for d in range(min(variables, degree)+1)
                   for combo in combinations(range(variables), d)),
                  key=lambda x: (x.bit_count(), x), reverse=True)


def echelon(rows, columns):
    pivots = {}
    for original in rows:
        row = original
        while row:
            k = (row & -row).bit_length()-1
            if k not in pivots:
                pivots[k] = row
                break
            row ^= pivots[k]
    return pivots


def closure(equations, variables, degree):
    masks = masks_through(variables, degree)
    index = {m: j for j, m in enumerate(masks)}
    rows = []
    for equation in equations:
        if equation and max(m.bit_count() for m in equation) <= degree:
            row = 0
            for m in equation:
                row ^= 1 << index[m]
            rows.append(row)
    history = []
    for _ in range(20):
        pivots = echelon(rows, len(masks))
        history.append(len(pivots))
        if len(masks)-1 in pivots:
            return {"contains_one": True, "rank_history": history,
                    "normal_form_columns": len(masks)}
        shifted = []
        for row in pivots.values():
            support = [masks[j] for j in range(len(masks)) if row >> j & 1]
            for variable in range(variables):
                products = [mask | (1 << variable) for mask in support]
                if any(mask.bit_count() > degree for mask in products):
                    continue
                product = 0
                for mask in products:
                    product ^= 1 << index[mask]
                shifted.append(product)
        next_pivots = echelon(list(pivots.values()) + shifted, len(masks))
        if len(next_pivots) == len(pivots):
            return {"contains_one": False, "rank_history": history,
                    "normal_form_columns": len(masks)}
        rows = list(next_pivots.values())
    raise RuntimeError("closure did not stabilize in 20 rounds")


def main():
    data = json.loads(SOURCE.read_text())
    target = data["systems"][2]
    assert target["scalar"] == 3
    checks = {}
    for arm, degrees in (("direct", (6, 7)), ("symmetrized", (3, 4))):
        system = target["presentations"][arm]
        assert not system["roots"]
        checks[arm] = {str(D): closure(system["equations"], system["variables"], D)
                       for D in degrees}
    assert not checks["direct"]["6"]["contains_one"]
    assert checks["direct"]["7"]["contains_one"]
    assert not checks["symmetrized"]["3"]["contains_one"]
    assert checks["symmetrized"]["4"]["contains_one"]
    result = {"schema": "torsion-symm-equal7-independent-degree/v1",
              "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "target_scalar": 3, "checks": checks}
    OUT.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    print(json.dumps(checks, sort_keys=True))


if __name__ == "__main__":
    main()
