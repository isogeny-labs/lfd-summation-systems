#!/usr/bin/env sage -python
"""Independent Sage point-addition replay of the selected relation sets."""

import hashlib
import itertools
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing

import generate

HERE = Path(__file__).resolve().parent


def digest(rows):
    return hashlib.sha256(json.dumps(sorted(rows), separators=(",", ":"))
                          .encode()).hexdigest()


def field(n, modulus):
    ring = PolynomialRing(GF(2), "z")
    z = ring.gen()
    defining = sum(z**i for i in range(n + 1) if modulus & (1 << i))
    ext = GF(2**n, "z", modulus=defining)
    gen = ext.gen()

    def element(value):
        return sum((gen**i for i in range(n) if value & (1 << i)), ext(0))

    lookup = {element(value): value for value in range(1 << n)}
    assert len(lookup) == 1 << n
    return ext, element, lookup


def main():
    results = []
    for n, modulus, k in generate.GRID:
        ext, element, to_int = field(n, modulus)
        for a in (0, 1):
            curve = EllipticCurve(ext, [ext(1), ext(a), ext(0), ext(0), ext(1)])
            lifts = {}
            for point in curve.points():
                if not point.is_zero():
                    lifts.setdefault(to_int[point[0]], []).append(point)
            base, _, _, selected = generate.point_data(n, modulus, a, 1, k)
            assert all(x in lifts for x in base)
            for target, expected in selected:
                actual = set()
                for quad in itertools.product(base, repeat=4):
                    if any(not (sum(pts, curve(0))).is_zero() and
                           to_int[sum(pts, curve(0))[0]] == target
                           for pts in itertools.product(*(lifts[x] for x in quad))):
                        actual.add(quad)
                assert actual == expected, (n, a, target,
                                            sorted(actual ^ expected)[:3])
                results.append({"n": n, "modulus": modulus, "k": k,
                                "a": a, "b": 1, "target": target,
                                "curve_order": int(curve.order()),
                                "base": base,
                                "relations": len(actual),
                                "tuple_sha256": digest(actual)})
    (HERE / "sage-replay.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps({"checked_queries": len(results), "mismatches": 0}))


if __name__ == "__main__":
    main()
