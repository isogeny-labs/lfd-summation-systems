"""Frozen three-summand matched cost pilot with timeout receipts."""
from __future__ import annotations

import argparse
import hashlib
import json
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
from make_m3_manifest import make  # noqa: E402
from run import rank_mod  # noqa: E402

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "m3_manifest.json"
OUT = HERE / "m3_pilot_results.json"
LIMIT_SECONDS = 45


def sha(data):
    return hashlib.sha256(data).hexdigest()


def stamp():
    return perf_counter(), process_time()


def elapsed(t):
    return perf_counter() - t[0], process_time() - t[1]


def emit(row):
    print(json.dumps(row, separators=(",", ":")), flush=True)


def worker(n, k, arm):
    fixture = next(x for x in json.loads(MANIFEST.read_text())["cells"]
                   if (x["n"], x["k"]) == (n, k))
    t0 = stamp()
    E = BinaryCurve(GF2n(n, fixture["modulus"]), 0, fixture["b"])
    S = Solver(E, [1 << i for i in range(k)])
    fb = sorted((int(x), int(y)) for x, y in zip(S.fx, S.fy))
    assert len(fb) == fixture["factor_base_points"]
    assert sha(json.dumps(fb, separators=(",", ":")).encode()) == fixture["factor_base_sha256"]
    fbset = set(fb)
    sw, sc = elapsed(t0)
    emit({"kind": "setup", "wall_seconds": sw, "cpu_seconds": sc,
          "factor_base_points": len(fb)})
    for i, target in enumerate(fixture["targets"]):
        R = tuple(target["target"])
        q0 = stamp()
        if arm == "lrf":
            rows, stats, _ = S.lrf(3, R)
        elif arm == "image":
            rows, stats, _ = S.lrf_image(3, R)
        elif arm == "brute":
            rows, _ = S.brute(3, R)
            stats = {}
        else:
            raise ValueError(arm)
        qw, qc = elapsed(q0)
        v0 = stamp()
        for triple in rows:
            assert all(p in fbset for p in triple)
            assert E.add(E.add(triple[0], triple[1]), triple[2]) == R
        relations = [[list(p) for p in triple] for triple in rows]
        if target["kind"] == "planted":
            assert target["witness"] in relations
        vw, vc = elapsed(v0)
        emit({"kind": "target", "index": i, "query_wall_seconds": qw,
              "query_cpu_seconds": qc, "verification_wall_seconds": vw,
              "verification_cpu_seconds": vc, "relations": relations, "stats": stats})
    tw, tc = elapsed(t0)
    emit({"kind": "finish", "wall_seconds": tw, "cpu_seconds": tc,
          "image_table_bytes": getattr(S, "image_table_bytes", 0),
          "peak_rss_raw": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
          "peak_rss_unit": "bytes" if sys.platform == "darwin" else "KiB"})


def child(n, k, arm):
    command = [sys.executable, str(Path(__file__).resolve()), "--worker",
               str(n), str(k), arm]
    t0 = perf_counter()
    proc = subprocess.Popen(command, cwd=ROOT, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True)
    try:
        stdout, stderr = proc.communicate(timeout=LIMIT_SECONDS)
        status = "complete" if proc.returncode == 0 else "worker_error"
    except subprocess.TimeoutExpired:
        proc.kill()
        stdout, stderr = proc.communicate()
        status = "timeout_censored"
    receipts = [json.loads(line) for line in stdout.splitlines() if line.strip()]
    if status == "complete":
        assert receipts[-1]["kind"] == "finish"
        wall, cpu = receipts[-1]["wall_seconds"], receipts[-1]["cpu_seconds"]
    else:
        wall = perf_counter() - t0
        cpu = sum(r.get("cpu_seconds", 0) + r.get("query_cpu_seconds", 0)
                  + r.get("verification_cpu_seconds", 0) for r in receipts)
    return {"status": status, "wall_seconds": wall, "cpu_seconds": cpu,
            "receipts": receipts, "stderr_tail": stderr[-1000:]}


def elimination(fixture, reference):
    t0 = stamp()
    E = BinaryCurve(GF2n(fixture["n"], fixture["modulus"]), 0, fixture["b"])
    h, q = fixture["cofactor"], fixture["q"]
    vectors, rhs_values, columns, projected = [], [], {}, {}
    receipts = [r for r in reference["receipts"] if r["kind"] == "target"]
    for target, receipt in zip(fixture["targets"], receipts):
        if target["kind"] != "random":
            break
        for triple in receipt["relations"]:
            v = {}
            for raw in triple:
                p = tuple(raw)
                if p not in projected:
                    projected[p] = E.mul(h, p)
                if projected[p] is None:
                    continue
                rep = min(p, E.neg(p))
                j = columns.setdefault(rep, len(columns))
                v[j] = v.get(j, 0) + (1 if p == rep else -1)
            vectors.append(v)
            rhs_values.append(-h * (target["u"] + 113 * target["v"]))
    rank = rank_mod(vectors, q)
    augmented = [v | {len(columns): rhs} for v, rhs in zip(vectors, rhs_values)]
    augmented_rank = rank_mod(augmented, q)
    assert rank == augmented_rank
    wall, cpu = elapsed(t0)
    return {"rows": len(vectors), "columns": len(columns), "rank": rank,
            "augmented_rank": augmented_rank, "wall_seconds": wall,
            "cpu_seconds": cpu}


def run():
    frozen = json.loads(MANIFEST.read_text())
    assert frozen == make()
    _lib()
    source = {name: sha((ROOT / name).read_bytes()) for name in
              ("csd/_rootfind.c", "csd/ic_rootfind.py")}
    source["experiments/subspace-remainder-cost-20261125/run_m3_pilot.py"] = sha(
        Path(__file__).read_bytes())
    if OUT.exists():
        result = json.loads(OUT.read_text())
        assert result["source_sha256"] == source
        assert result["manifest_sha256"] == sha(MANIFEST.read_bytes())
    else:
        result = {"manifest_sha256": sha(MANIFEST.read_bytes()),
                  "source_sha256": source,
                  "host": {"platform": platform.platform(),
                           "machine": platform.machine(), "python": sys.version},
                  "cells": []}
    done = {(x["n"], x["k"], x["repetition"]) for x in result["cells"]}
    for fixture in frozen["cells"]:
        n, k = fixture["n"], fixture["k"]
        for rep in range(3):
            if (n, k, rep) in done:
                continue
            order = ("lrf", "brute") if rep % 2 == 0 else ("brute", "lrf")
            arms = {name: child(n, k, name) for name in order}
            arms["image"] = child(n, k, "image")
            lists = [[r["relations"] for r in arm["receipts"] if r["kind"] == "target"]
                     for arm in arms.values()]
            prefix = min(map(len, lists))
            assert all(x[:prefix] == lists[0][:prefix] for x in lists[1:])
            lin = elimination(fixture, arms["lrf"]) if arms["lrf"]["status"] == "complete" else None
            crossover = all(a["status"] == "complete" for a in arms.values()) and all(
                arms["image"][metric] < arms[name][metric]
                for name in order for metric in ("wall_seconds", "cpu_seconds"))
            result["cells"].append({"n": n, "k": k, "repetition": rep,
                                    "arms": arms, "shared_completed_prefix": prefix,
                                    "linear_algebra": lin, "crossover": crossover})
            OUT.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
            print(n, k, rep, {a: v["status"] for a, v in arms.items()},
                  "rank", lin["rank"] if lin else None,
                  "crossover", crossover, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker", nargs=3)
    args = ap.parse_args()
    if args.worker:
        worker(int(args.worker[0]), int(args.worker[1]), args.worker[2])
    else:
        run()


if __name__ == "__main__":
    main()
