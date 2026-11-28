#!/usr/bin/env sage -python
"""Replay the frozen panel and compare every deterministic field and count."""

import json
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent


def without_timing(value):
    if isinstance(value, dict):
        return {key: without_timing(item) for key, item in value.items()
                if not key.endswith("_seconds") and key != "python_platform"}
    if isinstance(value, list):
        return [without_timing(item) for item in value]
    return value


def main():
    frozen = json.loads((HERE / "results.json").read_text())
    fresh = json.loads(subprocess.check_output(
        [sys.executable, str(HERE / "panel.py")], text=True))
    assert without_timing(fresh) == without_timing(frozen), (
        "deterministic panel output differs from the frozen record")
    print("Verified all 9 targets, relation counts, baseline counts, and fixtures")


if __name__ == "__main__":
    main()
