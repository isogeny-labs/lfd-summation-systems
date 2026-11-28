#!/usr/bin/env python3
"""Compare a fresh Sage replay with the frozen exact records."""

import json
import sys
from pathlib import Path


def strip_times(obj):
    if isinstance(obj, dict):
        return {k: strip_times(v) for k, v in obj.items()
                if "seconds" not in k}
    if isinstance(obj, list):
        return [strip_times(v) for v in obj]
    return obj


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: verify.py REPLAY_JSON")
    frozen = json.loads((Path(__file__).parent / "results.json").read_text())
    fresh = json.loads(Path(sys.argv[1]).read_text())
    assert strip_times(frozen) == strip_times(fresh)
    print("exact subfield pair-support replay verified")


if __name__ == "__main__":
    main()
