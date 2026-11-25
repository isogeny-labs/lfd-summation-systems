"""Bind the frozen depth-14/depth-18 prime inputs to the exact m=2 FFD theorem."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = Path(__file__).resolve().parents[2]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def is_prime(p: int) -> bool:
    if p < 2:
        return False
    for divisor in range(2, int(p**0.5) + 1):
        if p % divisor == 0:
            return False
    return True


def s3_coefficients(p: int, a: int, b: int, t: int) -> dict[tuple[int, int], int]:
    raw = {
        (2, 0): t * t, (0, 2): t * t, (1, 1): -2 * t * t - 2 * a,
        (2, 1): -2 * t, (1, 2): -2 * t,
        (1, 0): -2 * a * t - 4 * b, (0, 1): -2 * a * t - 4 * b,
        (0, 0): a * a - 4 * b * t, (2, 2): 1,
    }
    return {exponents: coefficient % p for exponents, coefficient in raw.items()
            if coefficient % p}


def parse_input(path: Path):
    lines = path.read_text().splitlines()
    header = list(map(int, lines[0].split()))
    p, m, k, aux, D = header[:5]
    assert m == 2 and aux == 0 and k >= 2 and len(lines) == 3 and lines[1] == "1"
    values = list(map(int, lines[2].split()))
    assert len(values) == 1 + 2 * values[0]
    encoded = dict(zip(values[1::2], values[2::2]))
    assert len(encoded) == values[0]
    assert all(0 < coefficient < p for coefficient in encoded.values())
    decoded = {}
    for mask, coefficient in encoded.items():
        e1, e2 = mask & ((1 << k) - 1), mask >> k
        assert e1 <= 2 and e2 <= 2 and mask == e1 + (e2 << k)
        decoded[(e1, e2)] = coefficient
    return p, k, D, decoded


def main() -> None:
    prior = json.loads((HERE / "audit-prime-ffd.json").read_text())
    theorem_sources = prior["theorem_source_sha256"]
    assert prior["theorem"] == "ICFallDegrees.PrimeFieldConclusions.ordinary_ffd"
    assert all(sha(REPO / name) == digest for name, digest in theorem_sources.items())

    cases = []
    for p in (65537, 163841):
        identity_path = HERE / f"publication-depth14-p{p}-identity.json"
        input_path = HERE / f"publication-depth14-p{p}-D6.txt"
        identity = json.loads(identity_path.read_text())
        assert identity["p"] == p and identity["k"] == 14 and identity["D"] == 6
        assert identity["input_sha256"] == sha(input_path)
        assert identity["target_on_curve"] and identity["target_nonsquare"]
        cases.append(("depth14", input_path, identity_path,
                      p, 14, identity["a"], identity["b"], identity["target"][0]))

    instance_path = HERE / "publication-depth18-instance.json"
    d7 = HERE / "publication-depth18-p7340033-D7.txt"
    d8 = HERE / "publication-depth18-p7340033-D8.txt"
    instance = json.loads(instance_path.read_text())
    assert (instance["prime"], instance["k"], instance["a"], instance["b"]) == (7340033, 18, 1, 0)
    assert instance["unit_nonresidue"]
    assert d7.read_text().splitlines()[1:] == d8.read_text().splitlines()[1:]
    p7, k7, D7, generator7 = parse_input(d7)
    assert (p7, k7, D7) == (7340033, 18, 7)
    assert generator7 == s3_coefficients(7340033, instance["a"], instance["b"], instance["target"][0])
    instance_generator = {tuple(map(int, key.split(","))): value for key, value in instance["generator"].items()}
    assert instance_generator == generator7
    cases.append(("depth18", d8, instance_path, 7340033, 18,
                  instance["a"], instance["b"], instance["target"][0]))

    cells = []
    for name, path, identity_path, p, k, a, b, t in cases:
        assert is_prime(p) and p >= 5
        parsed_p, parsed_k, D, generator = parse_input(path)
        assert (parsed_p, parsed_k) == (p, k)
        assert generator == s3_coefficients(p, a, b, t)
        condition = (3 * t * t + a) % p
        degree = 3 if condition == 0 else 4
        cells.append({
            "name": name, "field": f"GF({p})", "n": k, "m": 2,
            "D_of_frozen_input": D, "input_sha256": sha(path),
            "identity_sha256": sha(identity_path),
            "theorem_condition_residue": condition,
            "ffd_lower_bound": degree, "ffd_upper_bound": degree,
        })
    assert [cell["ffd_lower_bound"] for cell in cells] == [4, 4, 4]
    out = {
        "status": "passed", "scope": "three frozen publication prime m=2 systems",
        "theorem": prior["theorem"], "theorem_source_sha256": theorem_sources,
        "theorem_condition": prior["theorem_condition"],
        "cells": cells, "audit_source_sha256": sha(Path(__file__)),
    }
    (HERE / "audit-publication-prime-ffd.json").write_text(json.dumps(out, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"status": out["status"], "cells": len(cells),
                      "ffd": [cell["ffd_lower_bound"] for cell in cells]}))


if __name__ == "__main__":
    main()
