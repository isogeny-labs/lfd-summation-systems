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
    if len(sys.argv) not in (2, 3):
        raise SystemExit("usage: verify.py REPLAY_JSON [FROZEN_JSON]")
    frozen_name = sys.argv[2] if len(sys.argv) == 3 else "results.json"
    frozen = json.loads((Path(__file__).parent / frozen_name).read_text())
    fresh = json.loads(Path(sys.argv[1]).read_text())
    assert strip_times(frozen) == strip_times(fresh)
    print("exact subfield pair-support replay verified")


if __name__ == "__main__":
    main()
