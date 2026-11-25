#!/usr/bin/env python3
"""Bounded exact iterated LFD for the frozen equal-7 m3 Boolean systems."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from csd.ic_local_certificate import local_profile

HERE = Path(__file__).resolve().parent
SYSTEMS = HERE / "equal7-m3.json"
STORE = HERE / "equal7-m3-degree.jsonl"
ARMS = ("direct", "symmetrized")
LIMIT_BYTES = 512 * 1024 * 1024


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def append(path, value):
    with path.open("a") as out:
        out.write(json.dumps(value, sort_keys=True) + "\n")
        out.flush()
        os.fsync(out.fileno())
    print(value["id"], value.get("status"), flush=True)


def load(path):
    records = {}
    if path.exists():
        for line in path.read_text().splitlines():
            value = json.loads(line)
            if value["id"] in records:
                raise ValueError("duplicate record: " + value["id"])
            records[value["id"]] = value
    return records


def worker(scalar, arm):
    # local_profile applies this cap to every closure and evaluation matrix.
    # macOS rejects RLIMIT_AS even when setting a lower soft limit.
    data = json.loads(SYSTEMS.read_text())
    system = data["systems"][scalar - 1]["presentations"][arm]
    result = local_profile(system["equations"], system["variables"], system["roots"],
                           max_degree=8, max_columns=512, max_bytes=LIMIT_BYTES)
    print(json.dumps(result, sort_keys=True), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker", action="store_true")
    ap.add_argument("--scalar", type=int, choices=(1, 2, 3, 4))
    ap.add_argument("--arm", choices=ARMS)
    ap.add_argument("--store", type=Path, default=STORE)
    args = ap.parse_args()
    if args.worker:
        if args.scalar is None or args.arm is None:
            ap.error("worker requires scalar and arm")
        worker(args.scalar, args.arm)
        return
    data = json.loads(SYSTEMS.read_text())
    assert len(data["systems"]) == 4
    records = load(args.store)
    if not records:
        header = {"id": "header", "schema": "torsion-symm-equal7-m3-degree/v1",
                  "protocol": "PROTOCOL_EQUAL7_M3_20261125.md", "systems_sha256": digest(SYSTEMS),
                  "source_sha256": digest(Path(__file__)), "degree_cap": 8,
                  "normal_form_column_cap": 512, "memory_cap_bytes": LIMIT_BYTES,
                  "wall_cap_s_per_cell": 60}
        append(args.store, header)
        records["header"] = header
    if records["header"]["systems_sha256"] != digest(SYSTEMS):
        raise ValueError("systems changed after degree store was frozen")
    for scalar in (1, 2, 3, 4):
        for arm in ARMS:
            key = f"cell:{scalar}:{arm}"
            if key in records:
                continue
            system = data["systems"][scalar-1]["presentations"][arm]
            command = [sys.executable, str(Path(__file__).resolve()), "--worker",
                       "--scalar", str(scalar), "--arm", arm]
            before = resource.getrusage(resource.RUSAGE_CHILDREN)
            start = time.perf_counter()
            proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            timed_out = False
            try:
                stdout, stderr = proc.communicate(timeout=60)
            except subprocess.TimeoutExpired:
                timed_out = True
                proc.terminate()
                try:
                    stdout, stderr = proc.communicate(timeout=5)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    stdout, stderr = proc.communicate()
            after = resource.getrusage(resource.RUSAGE_CHILDREN)
            try:
                profile = json.loads(stdout)
            except json.JSONDecodeError:
                profile = None
            value = {"id": key, "scalar": scalar, "arm": arm, "command": command,
                     "equations_sha256": system["equations_sha256"],
                     "roots_sha256": system["roots_sha256"],
                     "root_count": len(system["roots"]), "input_degree": system["input_degree"],
                     "status": "timeout" if timed_out else ("completed" if proc.returncode == 0 else "error"),
                     "exit_code": proc.returncode, "wall_s": time.perf_counter()-start,
                     "child_cpu_s": (after.ru_utime+after.ru_stime)-(before.ru_utime+before.ru_stime),
                     "profile": profile, "stderr_tail": stderr[-4000:],
                     "stdout_tail": stdout[-1000:] if profile is None else None}
            append(args.store, value)
            records[key] = value


if __name__ == "__main__":
    main()
