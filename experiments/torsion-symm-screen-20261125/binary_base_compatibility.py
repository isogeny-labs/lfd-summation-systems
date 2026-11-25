#!/usr/bin/env python3
"""Screen frozen binary x-factor bases for rational two-torsion closure."""
from __future__ import annotations

from collections import defaultdict
import hashlib
import json
from pathlib import Path
import statistics
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from csd.gf2_ghs import GF2n  # noqa: E402
from csd.ic_dreg import BinaryCurve  # noqa: E402

HERE = Path(__file__).resolve().parent
CASES = ROOT / "experiments/ic-curve-sweep/cases.jsonl"


def span(basis):
    values = {0}
    for b in basis:
        values |= {x ^ b for x in list(values)}
    return values


def main():
    rows = [json.loads(line) for line in CASES.read_text().splitlines()]
    unique = {}
    for row in rows:
        if row.get("kind") not in ("binary", "koblitz") or "basis" not in row:
            continue
        key = (row["kind"], row["n"], row["modulus"], row["a"], row["b"],
               tuple(row["basis"]), row["base_kind"])
        unique.setdefault(key, row)
    records = []
    groups = defaultdict(list)
    checked_curves = set()
    for key, row in sorted(unique.items(), key=lambda item: str(item[0])):
        F = GF2n(row["n"], row["modulus"])
        curve_key = (row["n"], row["modulus"], row["a"], row["b"])
        if curve_key not in checked_curves:
            E = BinaryCurve(F, row["a"], row["b"])
            T2 = (0, F.sqrt(row["b"]))
            if not E.on_curve(T2) or E.add(T2, T2) is not None:
                raise AssertionError("rational order-two point check failed")
            tested = 0
            for x in range(1, 1 << row["n"]):
                P = E.lift_x(x)
                if P is None:
                    continue
                Q = E.add(P, T2)
                if Q is None or Q[0] != F.mul(T2[1], F.inv(x)):
                    raise AssertionError("two-torsion x-translation identity failed")
                tested += 1
                if tested == 16:
                    break
            checked_curves.add(curve_key)
        V = span(row["basis"])
        nonzero = V - {0}
        sqrt_b = F.sqrt(row["b"])
        inside = sum(F.mul(sqrt_b, F.inv(x)) in nonzero for x in nonzero)
        item = {"case_id": row["case_id"], "kind": row["kind"], "n": row["n"],
                "k": row["k"], "a": row["a"], "b": row["b"],
                "basis": row["basis"], "base_kind": row["base_kind"],
                "sqrt_b": sqrt_b,
                "nonzero_size": len(nonzero), "inverse_inside": inside,
                "closed": inside == len(nonzero)}
        records.append(item)
        groups[(row["kind"], row["n"], row["k"])].append(inside/len(nonzero))
    output = {"schema": "torsion-symm-frozen-binary-base-compatibility/v1",
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "cases_sha256": hashlib.sha256(CASES.read_bytes()).hexdigest(),
              "identity": "translation by (0,sqrt(b)) sends nonzero x to sqrt(b)/x",
              "checked_curve_count": len(checked_curves),
              "records": records,
              "closed_count": sum(x["closed"] for x in records),
              "summary": [{"kind": kind, "n": n, "k": k, "bases": len(values),
                           "median_inverse_coverage": statistics.median(values)}
                          for (kind, n, k), values in sorted(groups.items())]}
    path = HERE / "binary-base-compatibility.json"
    path.write_text(json.dumps(output, sort_keys=True, indent=2)+"\n")
    print("unique bases", len(records), "closed", output["closed_count"])
    for item in output["summary"]:
        if item["kind"] == "koblitz" and item["n"] in (5, 7, 9):
            print(item)


if __name__ == "__main__":
    main()
