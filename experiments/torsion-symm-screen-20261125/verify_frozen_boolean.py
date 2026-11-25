#!/usr/bin/env python3
"""Independent dense 64-column Boolean closure replay of frozen torsion data."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "frozen-boolean-torsion.json"
MONOMIALS = sorted(range(64), key=lambda mask: (mask.bit_count(), mask))
INDEX = {mask: index for index, mask in enumerate(MONOMIALS)}


def reduce_row(row, pivots):
    while row:
        pivot = row.bit_length() - 1
        if pivot not in pivots:
            return row
        row ^= pivots[pivot]
    return 0


def add_row(row, pivots):
    row = reduce_row(row, pivots)
    if not row:
        return False
    pivots[row.bit_length() - 1] = row
    return True


def product_by_variable(row, variable):
    out = 0
    while row:
        term = row & -row
        mask = MONOMIALS[term.bit_length() - 1]
        out ^= 1 << INDEX[mask | (1 << variable)]
        row ^= term
    return out


def closure(equations, degree):
    pivots = {}
    for eq in equations:
        if all(mask.bit_count() <= degree for mask in eq):
            add_row(sum(1 << INDEX[mask] for mask in eq), pivots)
    while True:
        if reduce_row(1, pivots) == 0:
            # The unit generates the complete normal-form space at this budget.
            return {INDEX[mask]: 1 << INDEX[mask] for mask in MONOMIALS if mask.bit_count() <= degree}
        changed = False
        for row in list(pivots.values()):
            # Charge ordinary degree before reducing x_i^2 to x_i.
            if MONOMIALS[row.bit_length() - 1].bit_count() >= degree:
                continue
            for variable in range(6):
                if add_row(product_by_variable(row, variable), pivots):
                    changed = True
        if not changed:
            break
    return pivots


def evaluation_rank(points, degree):
    pivots = {}
    for mask in range(64):
        if mask.bit_count() >= degree:
            continue
        row = sum(1 << j for j, bits in enumerate(points) if bits & mask == mask)
        add_row(row, pivots)
    return len(pivots)


def main():
    source = json.loads(SOURCE.read_text())
    records = []
    for case in source["cases"]:
        points = case["solution_masks"]
        for form, item in case["presentations"].items():
            equations = item["equations"]
            expected = item["closure"]
            ranks = []
            unit = []
            for degree in range(expected["scanned_through"] + 1):
                pivots = closure(equations, degree)
                ranks.append(len(pivots))
                unit.append(reduce_row(1, pivots) == 0)
            assert ranks == expected["quotient_ranks"], (case["case_id"], form, ranks, expected["quotient_ranks"])
            if points:
                D = expected["local_upper_bound"]
                assert evaluation_rank(points, D) == len(points)
                assert sum(m.bit_count() <= D for m in range(64)) - ranks[D] == len(points)
                assert not any(unit)
            else:
                D = expected["exact_lfd"]
                assert not any(unit[:D]) and unit[D]
            records.append({"case_id": case["case_id"], "form": form,
                            "dense_closure_ranks": ranks, "unit_by_degree": unit,
                            "certificate_degree": D,
                            "evaluation_rank_below_certificate": evaluation_rank(points, D) if points else None,
                            "verified_exact_lfd": expected["exact_lfd"]})
    output = {"schema": "torsion-symm-frozen-boolean-independent-closure/v1",
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "input_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              "records": records}
    (HERE / "frozen-boolean-verification.json").write_text(json.dumps(output, sort_keys=True, indent=2)+"\n")
    print("independent dense closure checks", len(records))


if __name__ == "__main__":
    main()
