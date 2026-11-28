#!/usr/bin/env python3
"""Exhaustive O-aware crossed-S3 witness screen over frozen binary curves."""

import hashlib
import itertools
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
FIELDS = ((3, 0b1011), (4, 0b10011))
CURVES = ((0, 1), (1, 1), (0, 2), (1, 2))


def mul(x, y, modulus):
    n = modulus.bit_length() - 1
    out = 0
    while y:
        if y & 1:
            out ^= x
        y >>= 1
        x <<= 1
        if x & (1 << n):
            x ^= modulus
    return out


def power(x, exponent, modulus):
    out = 1
    while exponent:
        if exponent & 1:
            out = mul(out, x, modulus)
        x = mul(x, x, modulus)
        exponent >>= 1
    return out


def inverse(x, modulus):
    assert x
    return power(x, (1 << (modulus.bit_length() - 1)) - 2, modulus)


def curve_rhs(x, a, b, modulus):
    x2 = mul(x, x, modulus)
    return mul(x2, x, modulus) ^ mul(a, x2, modulus) ^ b


def on_curve(point, a, b, modulus):
    if point is None:
        return True
    x, y = point
    return mul(y, y, modulus) ^ mul(x, y, modulus) == curve_rhs(x, a, b, modulus)


def add(left, right, a, modulus):
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2 and y1 ^ y2 == x1:
        return None
    if x1 == x2:
        slope = x1 ^ mul(y1, inverse(x1, modulus), modulus)
        x3 = mul(slope, slope, modulus) ^ slope ^ a
        y3 = mul(x1, x1, modulus) ^ mul(slope ^ 1, x3, modulus)
    else:
        slope = mul(y1 ^ y2, inverse(x1 ^ x2, modulus), modulus)
        x3 = mul(slope, slope, modulus) ^ slope ^ x1 ^ x2 ^ a
        y3 = mul(slope, x1 ^ x3, modulus) ^ x3 ^ y1
    return x3, y3


def s3(x, y, z, b, modulus):
    xy = mul(x, y, modulus)
    t = xy ^ mul(x, z, modulus) ^ mul(y, z, modulus)
    return mul(t, t, modulus) ^ mul(xy, z, modulus) ^ b


def label(point):
    return None if point is None else point[0]


def root_match(u, v, target, b, modulus):
    if target is None:
        return u == v
    if u is None:
        return v == target
    if v is None:
        return u == target
    return s3(u, v, target, b, modulus) == 0


def ordered(labels):
    return sorted(labels, key=lambda x: -1 if x is None else x)


def run_cell(n, modulus, a, b):
    q = 1 << n
    assert modulus.bit_length() == n + 1
    assert all(mul(x, inverse(x, modulus), modulus) == 1
               for x in range(1, q))
    points = [(x, y) for x in range(q) for y in range(q)
              if mul(y, y, modulus) ^ mul(x, y, modulus) ==
              curve_rhs(x, a, b, modulus)]
    group = [None] + points
    assert all(on_curve(pt, a, b, modulus) for pt in group)
    add_table = {(left, right): add(left, right, a, modulus)
                 for left, right in itertools.product(group, repeat=2)}
    assert all(on_curve(s, a, b, modulus) for s in add_table.values())
    assert all(add_table[None, pt] == pt == add_table[pt, None]
               for pt in group)
    assert all(add_table[pt, (pt[0], pt[1] ^ pt[0])] is None
               for pt in points)
    associativity_checks = 0
    for left, middle, right in itertools.product(group, repeat=3):
        assert add_table[add_table[left, middle], right] == \
               add_table[left, add_table[middle, right]]
        associativity_checks += 1

    xs = sorted({x for x, _ in points})
    lifts = {x: [pt for pt in points if pt[0] == x] for x in xs}
    pair_points = {}
    pair_labels = {}
    pair_mismatches = []
    for x, y in itertools.product(xs, repeat=2):
        sums = {add_table[left, right]
                for left, right in itertools.product(lifts[x], lifts[y])}
        pair_points[x, y] = sums
        direct = {label(pt) for pt in sums}
        polynomial = {z for z in range(q) if s3(x, y, z, b, modulus) == 0}
        if x == y:
            polynomial.add(None)
        if direct != polynomial:
            pair_mismatches.append({"x": x, "y": y,
                                    "direct": ordered(direct),
                                    "polynomial": ordered(polynomial)})
        pair_labels[x, y] = polynomial

    targets = [None] + xs
    checked = true_relations = naive_missed = 0
    tree_mismatches = []
    first_naive_miss = None
    for x1, x2, x3, x4 in itertools.product(xs, repeat=4):
        direct = {label(add_table[left, right])
                  for left, right in itertools.product(
                      pair_points[x1, x2], pair_points[x3, x4])}
        trees = (((x1, x2), (x3, x4)),
                 ((x1, x3), (x2, x4)))
        for target in targets:
            expected = target in direct
            witnesses = []
            affine = []
            for left, right in trees:
                pairs = itertools.product(pair_labels[left], pair_labels[right])
                witnesses.append(any(root_match(u, v, target, b, modulus)
                                     for u, v in pairs))
                pairs = itertools.product(pair_labels[left], pair_labels[right])
                affine.append(any(u is not None and v is not None and
                                  target is not None and
                                  s3(u, v, target, b, modulus) == 0
                                  for u, v in pairs))
            crossed = all(witnesses)
            checked += 1
            true_relations += expected
            if expected and not all(affine):
                naive_missed += 1
                if first_naive_miss is None:
                    first_naive_miss = [x1, x2, x3, x4, target]
            if witnesses != [expected, expected] or crossed != expected:
                tree_mismatches.append([x1, x2, x3, x4, target,
                                        expected, witnesses, crossed])
    return {"n": n, "modulus": modulus, "a": a, "b": b,
            "group_order": len(group), "x_count": len(xs),
            "associativity_checks": associativity_checks,
            "pair_checks": len(pair_points), "checked": checked,
            "true_relations": true_relations,
            "naive_affine_only_missed": naive_missed,
            "first_naive_miss": first_naive_miss,
            "pair_mismatches": pair_mismatches,
            "tree_mismatches": tree_mismatches}


def main():
    out = ROOT / "result.json"
    if out.exists():
        raise FileExistsError(out)
    rows = []
    for n, modulus in FIELDS:
        for a, b in CURVES:
            row = run_cell(n, modulus, a, b)
            rows.append(row)
            print(f"n={n} mod={modulus} a={a} b={b}: "
                  f"{row['checked']} cells, {len(row['pair_mismatches'])} pair "
                  f"and {len(row['tree_mismatches'])} tree mismatches")
    result = {"schema": "binary-crossed-s3-witness/v1", "cells": rows,
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "protocol_sha256": hashlib.sha256((ROOT / "PROTOCOL.md").read_bytes()).hexdigest()}
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    assert all(not row["pair_mismatches"] and not row["tree_mismatches"]
               for row in rows)


if __name__ == "__main__":
    main()
