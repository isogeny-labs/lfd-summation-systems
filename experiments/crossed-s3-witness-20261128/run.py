#!/usr/bin/env python3
"""Exhaustive O-aware crossed-S3 witness screen on frozen prime-field cells."""

import hashlib
import itertools
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
CELLS = ((7, 1, 1), (7, 1, 0), (11, 1, 1), (11, 1, 0),
         (13, 1, 1), (13, 1, 0))


def add(p, a, left, right):
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if x1 == x2:
        slope = (3 * x1 * x1 + a) * pow(2 * y1, -1, p) % p
    else:
        slope = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (slope * slope - x1 - x2) % p
    y3 = (slope * (x1 - x3) - y1) % p
    return x3, y3


def s3(p, a, b, x, y, z):
    return (((x - y) ** 2 * z ** 2
             - 2 * ((x + y) * (x * y + a) + 2 * b) * z
             + (x * y - a) ** 2 - 4 * b * (x + y)) % p)


def label(point):
    return None if point is None else point[0]


def root_match(p, a, b, u, v, target):
    """x-only S3 root, extended to the identity and target O."""
    if target is None:
        return u == v
    if u is None:
        return v == target
    if v is None:
        return u == target
    return s3(p, a, b, u, v, target) == 0


def main():
    out = ROOT / "result-v3.json"
    if out.exists():
        raise FileExistsError(out)
    records = []
    for p, a, b in CELLS:
        assert p > 3 and (4 * a ** 3 + 27 * b ** 2) % p != 0
        points = [(x, y) for x in range(p) for y in range(p)
                  if (y * y - x ** 3 - a * x - b) % p == 0]
        assert all((y * y - x ** 3 - a * x - b) % p == 0 for x, y in points)
        xs = sorted({x for x, _ in points})
        lifts = {x: [pt for pt in points if pt[0] == x] for x in xs}
        pair_points = {}
        pair_labels = {}
        pair_mismatches = []
        for x, y in itertools.product(xs, repeat=2):
            sums = {add(p, a, left, right)
                    for left, right in itertools.product(lifts[x], lifts[y])}
            pair_points[x, y] = sums
            direct_labels = {label(pt) for pt in sums}
            polynomial_labels = {z for z in range(p) if s3(p, a, b, x, y, z) == 0}
            if x == y:
                polynomial_labels.add(None)
            if direct_labels != polynomial_labels:
                pair_mismatches.append({"x": x, "y": y,
                                        "direct": sorted(direct_labels,
                                                         key=lambda z: -1 if z is None else z),
                                        "polynomial": sorted(polynomial_labels,
                                                             key=lambda z: -1 if z is None else z)})
            pair_labels[x, y] = polynomial_labels
        target_labels = [None] + xs
        checked = true_relations = cancelling = naive_missed = 0
        failures = []
        first_naive_miss = None
        for x1, x2, x3, x4 in itertools.product(xs, repeat=4):
            direct_labels = {label(add(p, a, left, right))
                             for left, right in itertools.product(
                                 pair_points[x1, x2], pair_points[x3, x4])}
            families = (((x1, x2), (x3, x4)),
                        ((x1, x3), (x2, x4)))
            for target in target_labels:
                expected = target in direct_labels
                witnesses = []
                naive = []
                for left, right in families:
                    witnesses.append(any(root_match(p, a, b, u, v, target)
                                         for u, v in itertools.product(
                                             pair_labels[left], pair_labels[right])))
                    naive.append(any(u is not None and v is not None and
                                     target is not None and
                                     s3(p, a, b, u, v, target) == 0
                                     for u, v in itertools.product(
                                         pair_labels[left], pair_labels[right])))
                crossed = all(witnesses)
                checked += 1
                true_relations += expected
                if expected and (None in pair_labels[x1, x2] or
                                 None in pair_labels[x3, x4] or
                                 None in pair_labels[x1, x3] or
                                 None in pair_labels[x2, x4]):
                    cancelling += 1
                if expected and not all(naive):
                    naive_missed += 1
                    if first_naive_miss is None:
                        first_naive_miss = [x1, x2, x3, x4, target]
                if witnesses != [expected, expected] or crossed != expected:
                    failures.append([x1, x2, x3, x4, target,
                                     expected, witnesses, crossed])
        records.append({"p": p, "a": a, "b": b, "group_order": len(points) + 1,
                        "x_count": len(xs), "checked": checked,
                        "true_relations": true_relations,
                        "true_relations_with_cancelling_pair": cancelling,
                        "naive_affine_only_missed": naive_missed,
                        "first_naive_miss": first_naive_miss,
                        "pair_mismatches": pair_mismatches,
                        "mismatches": failures})
        print(f"p={p} a={a} b={b}: {checked} tuple-target cells, "
              f"{len(pair_mismatches)} pair-root mismatches, "
              f"{len(failures)} O-aware mismatches, {naive_missed} affine-only misses")
    result = {"schema": "crossed-s3-witness/v1", "cells": records,
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "protocol_sha256": hashlib.sha256((ROOT / "PROTOCOL.md").read_bytes()).hexdigest(),
              "amendment_sha256": hashlib.sha256((ROOT / "AMENDMENT.md").read_bytes()).hexdigest(),
              "amendment_a2_sha256": hashlib.sha256((ROOT / "AMENDMENT-A2.md").read_bytes()).hexdigest()}
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    assert all(not row["pair_mismatches"] and not row["mismatches"]
               for row in records)


if __name__ == "__main__":
    main()
