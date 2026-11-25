#!/usr/bin/env python3
"""Rank the equal-67 checked equations on factor-base columns only."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from verify_equal67 import rank

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "equal67-koblitz.json"
OUT = HERE / "equal67-rank-audit.json"


def main():
    data = json.loads(SOURCE.read_text())
    output = {"schema": "torsion-symm-equal67-rank-audit/v1",
              "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "field_order": 211, "arms": {}}
    for arm in ("direct", "symmetrized"):
        saved = data["relations"][arm]
        n = len(saved["columns"])
        equations = []
        through = []
        for j, target in enumerate(saved["targets"]):
            assert target["scalar"] == j+1
            for row in target["rows"]:
                assert len(row) == n+24
                assert row[n+j] == 210 and all(row[n+k] == 0 for k in range(24) if k != j)
                # [154]P1+[154]P2=[154]R; target R=(j+1)Q.
                rhs = (154*(j+1)) % 211
                equations.append((row[:n], rhs))
            base_rank = rank([a for a, _ in equations], 211)
            augmented_rank = rank([a+[b] for a, b in equations], 211)
            assert base_rank == augmented_rank  # independently checked actual logs exist
            through.append({"scalar": j+1, "base_rank": base_rank,
                            "rhs_augmented_rank": augmented_rank,
                            "new_rows": len(target["rows"])})
        output["arms"][arm] = {"base_columns": n, "equations": len(equations),
                                "admitted_targets": saved["admitted_targets"],
                                "base_rank": through[-1]["base_rank"],
                                "rank_through_target": through}
    OUT.write_text(json.dumps(output, sort_keys=True, indent=2) + "\n")
    print(json.dumps({k: {"base_rank": v["base_rank"],
                          "admitted_targets": v["admitted_targets"]}
                      for k, v in output["arms"].items()}, sort_keys=True))


if __name__ == "__main__":
    main()
