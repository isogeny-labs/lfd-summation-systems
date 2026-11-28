"""Independent full Boolean-grid audit of the frozen lifted equations."""

import json
import time
from pathlib import Path

from run import affine_data, exact_pairs, source_check
from csd import ic_dreg as I
from csd import ic_rootfind as rf

HERE = Path(__file__).resolve().parent


def main():
    source_check()
    saved = json.loads((HERE / "results.json").read_text())
    total_tested = 0
    start = time.perf_counter()
    for cell in saved["cells"]:
        n, k, delta = cell["n"], cell["k"], cell["delta"]
        F = I.GF2n(n, (1 << n) | (1 << 5) | 1)
        E = I.BinaryCurve(F, 0, cell["b"])
        S = rf.Solver(E, [1 << i for i in range(k)])
        for category in ("empty", "positive"):
            for case in cell["selected"][category]:
                T = tuple(case["T"])
                data = affine_data(S, T)
                assert data is not None and data["kernel_dim"] == delta
                direct = []
                for gmask in range(1 << delta):
                    e1, e2, _ = data["evaluation"](gmask)
                    for zmask in range(1 << k):
                        if F.sq(zmask) ^ F.mul(e1, zmask) ^ e2 == 0:
                            direct.append([gmask, zmask])
                assert direct == case["pairs"]
                pairs, rows = exact_pairs(S, T, data)
                assert pairs == direct and [list(map(list, row)) for row in rows] == case["rows"]
                assert rows == S.d2_reference(T)
                assert rows == S.lrf(2, T)[0]
                total_tested += (1 << delta) * (1 << k)
                print(f"n={n} k={k} {category} attempt={case['attempt']}: "
                      f"{len(direct)} exact Boolean solutions", flush=True)
    print(f"audited {total_tested} Boolean assignments in "
          f"{time.perf_counter()-start:.2f}s")


if __name__ == "__main__":
    main()
