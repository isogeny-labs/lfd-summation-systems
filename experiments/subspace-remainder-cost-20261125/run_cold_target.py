"""Cold same-target image query versus conservative generic BSGS."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import resource
import subprocess
import sys
from pathlib import Path
from time import perf_counter, process_time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from csd.gf2_ghs import GF2n  # noqa: E402
from csd.ic_dreg import BinaryCurve  # noqa: E402
from csd.ic_rootfind import Solver, _lib  # noqa: E402
from make_cold_manifest import make  # noqa: E402

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "cold_target_manifest.json"
OUT = HERE / "cold_target_results.json"
CAP_SECONDS = 45


def sha(data):
    return hashlib.sha256(data).hexdigest()


def clock():
    return perf_counter(), process_time()


def elapsed(t):
    return perf_counter() - t[0], process_time() - t[1]


def bsgs(E, G, Q, q):
    width = math.isqrt(q) + 1
    table = {}
    P = None
    for j in range(width):
        table.setdefault(P, j)
        P = E.add(P, G)
    step = E.neg(E.mul(width, G))
    X = Q
    for i in range(width + 1):
        if X in table:
            result = (i*width + table[X]) % q
            assert E.mul(result, G) == Q
            return result, len(table), i+1
        X = E.add(X, step)
    raise AssertionError("BSGS exhausted subgroup order")


def worker(n, k, target_index, arm):
    cell = next(x for x in json.loads(MANIFEST.read_text())["cells"]
                if (x["n"], x["k"]) == (n, k))
    target = cell["targets"][target_index]
    Q, G = tuple(target["Q"]), tuple(cell["G"])
    t0 = clock()
    E = BinaryCurve(GF2n(n, cell["modulus"]), 0, cell["b"])
    if arm == "bsgs":
        scalar, table_size, giant_steps = bsgs(E, G, Q, cell["q"])
        assert scalar == target["secret_check_only"]
        outcome = {"scalar": scalar, "table_size": table_size,
                   "giant_steps": giant_steps}
    elif arm == "image":
        S = Solver(E, [1 << i for i in range(k)])
        fb = sorted((int(x), int(y)) for x, y in zip(S.fx, S.fy))
        assert len(fb) == cell["factor_base_points"]
        assert sha(json.dumps(fb, separators=(",", ":")).encode()) == cell["factor_base_sha256"]
        fbset = set(fb)
        triples, stats, _ = S.lrf_image(3, Q)
        for triple in triples:
            assert all(p in fbset for p in triple)
            assert E.add(E.add(triple[0], triple[1]), triple[2]) == Q
        outcome = {"relations": [[list(p) for p in triple] for triple in triples],
                   "stats": stats,
                   "image_table_bytes": S.image_table_bytes}
    else:
        raise ValueError(arm)
    wall, cpu = elapsed(t0)
    print(json.dumps({"status": "complete", "arm": arm, "n": n, "k": k,
                      "target_index": target_index, "wall_seconds": wall,
                      "cpu_seconds": cpu, "peak_rss_raw": resource.getrusage(
                          resource.RUSAGE_SELF).ru_maxrss,
                      "peak_rss_unit": "bytes" if sys.platform == "darwin" else "KiB",
                      "outcome": outcome}, separators=(",", ":")), flush=True)


def child(n, k, i, arm):
    cmd = [sys.executable, str(Path(__file__).resolve()), "--worker",
           str(n), str(k), str(i), arm]
    t0 = perf_counter()
    try:
        process = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True,
                                 timeout=CAP_SECONDS)
    except subprocess.TimeoutExpired as exc:
        return {"status": "timeout_censored", "wall_lower_bound_seconds": perf_counter()-t0,
                "stdout_prefix": str(exc.stdout or "")[-1000:],
                "stderr_tail": str(exc.stderr or "")[-1000:]}
    if process.returncode != 0:
        return {"status": "worker_error", "wall_lower_bound_seconds": perf_counter()-t0,
                "stderr_tail": process.stderr[-1500:]}
    return json.loads(process.stdout)


def run():
    t0 = clock()
    frozen = json.loads(MANIFEST.read_text())
    assert frozen == make()
    common_wall, common_cpu = elapsed(t0)
    _lib()  # Common native build before timing either arm.
    source = {name: sha((ROOT / name).read_bytes()) for name in
              ("csd/_rootfind.c", "csd/ic_rootfind.py",
               "csd/gf2_ghs.py", "csd/ic_dreg.py")}
    source["experiments/subspace-remainder-cost-20261125/run_cold_target.py"] = sha(
        Path(__file__).read_bytes())
    if OUT.exists():
        data = json.loads(OUT.read_text())
        assert data["source_sha256"] == source
        assert data["manifest_sha256"] == sha(MANIFEST.read_bytes())
    else:
        data = {"manifest_sha256": sha(MANIFEST.read_bytes()),
                "source_sha256": source,
                "common_manifest_check_wall_seconds": common_wall,
                "common_manifest_check_cpu_seconds": common_cpu,
                "host": {"platform": platform.platform(), "machine": platform.machine(),
                         "python": sys.version,
                         "compiler": os.popen("cc --version").read().splitlines()[:1]},
                "trials": []}
    done = {(x["n"], x["k"], x["target_index"], x["repetition"]) for x in data["trials"]}
    for cell in frozen["cells"]:
        n, k = cell["n"], cell["k"]
        for target in cell["targets"]:
            i = target["index"]
            for rep in range(3):
                if (n, k, i, rep) in done:
                    continue
                order = ("bsgs", "image") if rep % 2 == 0 else ("image", "bsgs")
                arms = {name: child(n, k, i, name) for name in order}
                b, a = arms["bsgs"], arms["image"]
                bsgs_faster = (b["status"] == a["status"] == "complete" and
                               b["wall_seconds"] < a["wall_seconds"] and
                               b["cpu_seconds"] < a["cpu_seconds"])
                data["trials"].append({"n": n, "k": k, "target_index": i,
                                       "repetition": rep, "arms": arms,
                                       "bsgs_faster": bsgs_faster})
                OUT.write_text(json.dumps(data, sort_keys=True, indent=2) + "\n")
                print(n, k, i, rep, {name: v["status"] for name, v in arms.items()},
                      "bsgs_faster", bsgs_faster, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker", nargs=4)
    args = ap.parse_args()
    if args.worker:
        n, k, i, arm = args.worker
        worker(int(n), int(k), int(i), arm)
    else:
        run()


if __name__ == "__main__":
    main()
