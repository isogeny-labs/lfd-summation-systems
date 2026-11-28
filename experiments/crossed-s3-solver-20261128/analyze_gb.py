#!/usr/bin/env python3
"""Verify custody and summarize the frozen GB-only replay."""

import hashlib
import json
import re
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROUND = re.compile(r"^\s*(\d+)\s+\d+\s+\d+\s+(\d+)\s+x\s+(\d+)\s+")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse(row):
    for key in ("stdout", "stderr", "output"):
        path = HERE / row[key]
        assert path.exists()
        assert digest(path) == row[key + "_sha256"]
    assert digest(HERE / row["input"]) == row["sha256"]
    assert row["status"] == "completed" and row["returncode"] == 0
    log = (HERE / row["stderr"]).read_text()
    output = (HERE / row["output"]).read_text()
    rounds = [tuple(map(int, match.groups())) for line in log.splitlines()
              if (match := ROUND.match(line))]
    basis = re.search(r"#length of basis:\s+(\d+)", output)
    assert basis is not None
    assert "#invalid equations               0" in log
    return {"max_f4_degree": max((degree for degree, _, _ in rounds),
                                  default=None),
            "f4_rounds": len(rounds),
            "cumulative_matrix_cells": sum(r * c for _, r, c in rounds),
            "largest_matrix_cells": max((r * c for _, r, c in rounds),
                                        default=0),
            "basis_elements": int(basis.group(1)),
            "rows_reduced": int(re.search(r"#rows reduced\s+(\d+)", log).group(1)),
            "zero_reductions": int(re.search(r"#zero reductions\s+(\d+)", log).group(1))}


def main():
    runs = json.loads((HERE / "runs-gb.json").read_text())
    assert len(runs) == 81
    grouped = {}
    for row in runs:
        key = (row["n"], row["a"], row["b"], row["target"])
        arm = row["arm"]
        grouped.setdefault(key, {}).setdefault(arm, []).append({**row, **parse(row)})
    assert len(grouped) == 9
    cases = []
    for key in sorted(grouped):
        arms = grouped[key]
        assert set(arms) == {"tree", "crossed", "direct"}
        summary = {}
        for arm, repetitions in arms.items():
            assert sorted(x["repetition"] for x in repetitions) == [0, 1, 2]
            geometry = ("max_f4_degree", "f4_rounds",
                        "cumulative_matrix_cells", "largest_matrix_cells",
                        "basis_elements", "rows_reduced", "zero_reductions")
            assert all(tuple(x[k] for k in geometry) ==
                       tuple(repetitions[0][k] for k in geometry)
                       for x in repetitions)
            summary[arm] = {k: repetitions[0][k] for k in geometry}
            summary[arm].update({
                "median_wall_s": statistics.median(x["wall_s"] for x in repetitions),
                "median_child_cpu_s": statistics.median(x["child_cpu_s"]
                                                        for x in repetitions),
                "walls_s": [x["wall_s"] for x in repetitions],
                "child_cpus_s": [x["child_cpu_s"] for x in repetitions]})
        tree = summary["tree"]
        crossed = summary["crossed"]
        direct = summary["direct"]
        cases.append({"n": key[0], "a": key[1], "b": key[2],
                      "target": key[3], "relations": arms["tree"][0]["relations"],
                      "arms": summary,
                      "crossed_over_tree_wall": crossed["median_wall_s"] /
                      tree["median_wall_s"],
                      "crossed_over_tree_matrix_cells":
                      crossed["cumulative_matrix_cells"] /
                      tree["cumulative_matrix_cells"],
                      "crossed_minus_tree_f4_degree":
                      crossed["max_f4_degree"] - tree["max_f4_degree"],
                      "direct_over_tree_wall": direct["median_wall_s"] /
                      tree["median_wall_s"]})
    result = {"solver": "msolve 0.11.0 -g 2 -c 0 -v 2 -t 1",
              "scope": "GB stage only; no parametrization or extraction",
              "completed_runs": len(runs), "cases": cases,
              "crossed_faster_cases": sum(c["crossed_over_tree_wall"] < 1
                                          for c in cases),
              "crossed_lower_f4_degree_cases": sum(
                  c["crossed_minus_tree_f4_degree"] < 0 for c in cases),
              "median_crossed_over_tree_wall": statistics.median(
                  c["crossed_over_tree_wall"] for c in cases),
              "median_crossed_over_tree_matrix_cells": statistics.median(
                  c["crossed_over_tree_matrix_cells"] for c in cases)}
    (HERE / "analysis-gb.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: result[k] for k in (
        "completed_runs", "crossed_faster_cases",
        "crossed_lower_f4_degree_cases",
        "median_crossed_over_tree_wall",
        "median_crossed_over_tree_matrix_cells")}, indent=2))


if __name__ == "__main__":
    main()
