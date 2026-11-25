"""Independent replay of the frozen same-target cost receipts."""
from __future__ import annotations

import gzip
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from csd.gf2_ghs import GF2n  # noqa: E402
from csd.ic_dreg import BinaryCurve  # noqa: E402
from make_cold_manifest import make  # noqa: E402

HERE = Path(__file__).resolve().parent
RAW = HERE / "cold_target_results.json.gz"
OUT = HERE / "cold_target_verification.json"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def verify():
    manifest_raw = (HERE / "cold_target_manifest.json").read_bytes()
    frozen = json.loads(manifest_raw)
    assert frozen == make()
    with gzip.open(RAW, "rt") as stream:
        results = json.load(stream)
    assert results["manifest_sha256"] == sha(manifest_raw)
    for path, digest in results["source_sha256"].items():
        assert sha((ROOT / path).read_bytes()) == digest, path
    assert len(results["trials"]) == 27
    seen = set()
    summary = []
    for cell in frozen["cells"]:
        n, k = cell["n"], cell["k"]
        E = BinaryCurve(GF2n(n, cell["modulus"]), 0, cell["b"])
        G = tuple(cell["G"])
        fb = []
        for x in range(1 << k):
            p = E.lift_x(x)
            if p is not None:
                fb.append(p)
                if E.neg(p) != p:
                    fb.append(E.neg(p))
        fb.sort()
        assert len(fb) == cell["factor_base_points"]
        assert sha(json.dumps(fb, separators=(",", ":")).encode()) == cell["factor_base_sha256"]
        fbset = set(fb)
        for target in cell["targets"]:
            Q, s = tuple(target["Q"]), target["secret_check_only"]
            assert E.mul(s, G) == Q
            sets = []
            for rep in range(3):
                key = n, k, target["index"], rep
                assert key not in seen
                seen.add(key)
                row = next(r for r in results["trials"]
                           if (r["n"], r["k"], r["target_index"], r["repetition"]) == key)
                b, a = row["arms"]["bsgs"], row["arms"]["image"]
                assert b["status"] == a["status"] == "complete"
                assert b["outcome"]["scalar"] == s
                triples = a["outcome"]["relations"]
                assert triples == sorted(triples)
                for triple in triples:
                    assert len(triple) == 3
                    points = [tuple(p) for p in triple]
                    assert points == sorted(points)
                    assert all(p in fbset for p in points)
                    assert E.add(E.add(points[0], points[1]), points[2]) == Q
                sets.append(triples)
                assert row["bsgs_faster"] == (
                    b["wall_seconds"] < a["wall_seconds"] and
                    b["cpu_seconds"] < a["cpu_seconds"])
            assert sets[0] == sets[1] == sets[2]
        rows = [r for r in results["trials"] if (r["n"], r["k"]) == (n, k)]
        summary.append({"n": n, "k": k, "pairs": len(rows),
                        "bsgs_faster": sum(r["bsgs_faster"] for r in rows),
                        "relation_counts": [len(next(r for r in rows if r["target_index"] == i)
                                                ["arms"]["image"]["outcome"]["relations"])
                                            for i in range(3)]})
    assert len(seen) == len(results["trials"])
    return {"raw_sha256": sha(RAW.read_bytes()), "cells": summary}


if __name__ == "__main__":
    data = verify()
    OUT.write_text(json.dumps(data, sort_keys=True, indent=2) + "\n")
    print(json.dumps(data, sort_keys=True))
