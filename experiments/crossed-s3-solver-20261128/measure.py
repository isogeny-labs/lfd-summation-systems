#!/usr/bin/env python3
"""Run frozen matched msolve panel with raw records and censoring."""

import hashlib
import json
import resource
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROWS = json.loads((HERE / "manifest.json").read_text())
LIMIT = 2 * 1024**3
ARMS = ("tree", "crossed", "direct")


def limit_memory():
    resource.setrlimit(resource.RLIMIT_AS, (LIMIT, LIMIT))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    raw = HERE / "raw"
    raw.mkdir(exist_ok=True)
    grouped = {}
    for row in ROWS:
        key = (row["n"], row["a"], row["b"], row["target"])
        grouped.setdefault(key, {})[row["arm"]] = row
    assert all(set(group) == set(ARMS) for group in grouped.values())
    version = subprocess.run(["msolve", "-V"], capture_output=True, text=True)
    assert version.returncode == 0
    runs = []
    for key in sorted(grouped):
        for rep in range(3):
            for arm in ARMS[rep:] + ARMS[:rep]:
                row = grouped[key][arm]
                source = HERE / row["input"]
                assert sha(source) == row["sha256"]
                tag = f"n{key[0]}_a{key[1]}_b{key[2]}_t{key[3]}_{arm}_r{rep}"
                stdout_path = raw / f"{tag}.stdout"
                stderr_path = raw / f"{tag}.stderr"
                output_path = raw / f"{tag}.out"
                cmd = ["msolve", "-f", str(source), "-o", str(output_path),
                       "-v", "2", "-t", "1", "--random-seed", "12345"]
                ru0 = resource.getrusage(resource.RUSAGE_CHILDREN)
                start = time.monotonic()
                try:
                    proc = subprocess.run(cmd, capture_output=True, timeout=20,
                                          preexec_fn=limit_memory)
                    status = "completed" if proc.returncode == 0 else "error"
                    rc = proc.returncode
                    stdout, stderr = proc.stdout, proc.stderr
                except subprocess.TimeoutExpired as exc:
                    status, rc = "timeout_20s", None
                    stdout, stderr = exc.stdout or b"", exc.stderr or b""
                elapsed = time.monotonic() - start
                ru1 = resource.getrusage(resource.RUSAGE_CHILDREN)
                stdout_path.write_bytes(stdout)
                stderr_path.write_bytes(stderr)
                record = {**row, "repetition": rep, "status": status,
                          "returncode": rc, "wall_s": elapsed,
                          "child_cpu_s": (ru1.ru_utime + ru1.ru_stime -
                                          ru0.ru_utime - ru0.ru_stime),
                          "child_maxrss_raw": ru1.ru_maxrss,
                          "stdout": str(stdout_path.relative_to(HERE)),
                          "stderr": str(stderr_path.relative_to(HERE)),
                          "output": str(output_path.relative_to(HERE))
                          if output_path.exists() else None,
                          "stdout_sha256": sha(stdout_path),
                          "stderr_sha256": sha(stderr_path),
                          "output_sha256": sha(output_path)
                          if output_path.exists() else None,
                          "msolve_version": version.stdout.strip() or
                          version.stderr.strip(),
                          "command": cmd, "timeout_s": 20,
                          "memory_limit_bytes": LIMIT}
                runs.append(record)
                (HERE / "runs.json").write_text(json.dumps(runs, indent=2) + "\n")
                print(json.dumps({k: record[k] for k in
                                  ("n", "a", "target", "arm", "repetition",
                                   "status", "wall_s")}), flush=True)


if __name__ == "__main__":
    main()
