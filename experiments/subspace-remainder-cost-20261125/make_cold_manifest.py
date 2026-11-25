"""Freeze matched cold target scalars and public points before timing."""
from __future__ import annotations

import hashlib
import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from csd.gf2_ghs import GF2n  # noqa: E402
from csd.ic_dreg import BinaryCurve  # noqa: E402

HERE = Path(__file__).resolve().parent
PARENT = HERE / "manifest.json"
PROTOCOL = HERE / "COLD_TARGET_PROTOCOL.md"
OUT = HERE / "cold_target_manifest.json"
CASES = ((23, 11), (29, 13), (29, 14))


def sha(data):
    return hashlib.sha256(data).hexdigest()


def make():
    parent = json.loads(PARENT.read_text())
    cells = []
    for n, k in CASES:
        old = next(x for x in parent["cells"] if (x["n"], x["k"]) == (n, k))
        E = BinaryCurve(GF2n(n, old["modulus"]), 0, old["b"])
        G = tuple(old["G"])
        assert E.mul(old["q"], G) is None
        targets = []
        for i in range(3):
            secret = random.Random(20264125 + 100*n + i).randrange(old["q"]//2, old["q"])
            Q = E.mul(secret, G)
            assert Q is not None
            targets.append({"index": i, "secret_check_only": secret, "Q": list(Q)})
        cells.append({key: old[key] for key in
                      ("n", "k", "modulus", "b", "order", "q", "cofactor", "G",
                       "factor_base_points", "factor_base_sha256")} | {"targets": targets})
    return {"parent_manifest_sha256": sha(PARENT.read_bytes()),
            "protocol_sha256": sha(PROTOCOL.read_bytes()),
            "source_sha256": sha(Path(__file__).read_bytes()), "cells": cells}


def main():
    data = make()
    if OUT.exists():
        assert json.loads(OUT.read_text()) == data, "cold target manifest changed"
    else:
        OUT.write_text(json.dumps(data, sort_keys=True, indent=2) + "\n")
    for c in data["cells"]:
        print(c["n"], c["k"], c["factor_base_points"], len(c["targets"]))


if __name__ == "__main__":
    main()
