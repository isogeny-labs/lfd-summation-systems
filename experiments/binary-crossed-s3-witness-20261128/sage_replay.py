#!/usr/bin/env sage -python
"""Independent Sage point-arithmetic replay of the frozen binary S3 cells."""

import hashlib
import itertools
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing


ROOT = Path(__file__).resolve().parent
CELLS = ((3, 0b1011, 0, 1), (3, 0b1011, 1, 1),
         (3, 0b1011, 0, 2), (3, 0b1011, 1, 2),
         (4, 0b10011, 0, 1), (4, 0b10011, 1, 1),
         (4, 0b10011, 0, 2), (4, 0b10011, 1, 2))


def replay_cell(n, modulus, a, b, expected):
    poly = PolynomialRing(GF(2), 'z')
    z = poly.gen()
    field = GF(1 << n, name='u', modulus=sum(
        ((z ** j for j in range(n + 1) if modulus >> j & 1)), poly.zero()))
    u = field.gen()
    values = [sum((u ** j for j in range(n) if i >> j & 1), field(0))
              for i in range(1 << n)]

    def integer(v):
        coeffs = v.polynomial()
        return sum(int(coeffs[j]) << j for j in range(n))

    assert all(integer(values[i]) == i for i in range(1 << n))
    curve = EllipticCurve(field, [field(1), values[a], field(0), field(0), values[b]])
    points = curve.points()
    xs = sorted({integer(pt[0]) for pt in points if not pt.is_zero()})
    lifts = {x: [pt for pt in points if not pt.is_zero() and
                 integer(pt[0]) == x] for x in xs}

    def label(pt):
        return None if pt.is_zero() else integer(pt[0])

    def s3(x, y, t):
        x, y, t = values[x], values[y], values[t]
        return (x*y + x*t + y*t)**2 + x*y*t + values[b]

    def root_match(left, right, target):
        if target is None:
            return left == right
        if left is None:
            return right == target
        if right is None:
            return left == target
        return s3(left, right, target) == 0

    pair_points = {}
    pair_labels = {}
    pair_mismatches = []
    for x, y in itertools.product(xs, repeat=2):
        sums = {left + right for left, right in itertools.product(lifts[x], lifts[y])}
        pair_points[x, y] = sums
        direct = {label(pt) for pt in sums}
        roots = {t for t in range(1 << n) if s3(x, y, t) == 0}
        if x == y:
            roots.add(None)
        if roots != direct:
            pair_mismatches.append([x, y, sorted(direct, key=str),
                                    sorted(roots, key=str)])
        pair_labels[x, y] = roots

    checked = true_relations = naive_missed = 0
    tree_mismatches = []
    targets = [None] + xs
    for x1, x2, x3, x4 in itertools.product(xs, repeat=4):
        direct = {label(left + right) for left, right in itertools.product(
            pair_points[x1, x2], pair_points[x3, x4])}
        trees = (((x1, x2), (x3, x4)), ((x1, x3), (x2, x4)))
        for target in targets:
            answer = target in direct
            matches = []
            affine = []
            for first, second in trees:
                pairs = itertools.product(pair_labels[first], pair_labels[second])
                matches.append(any(root_match(left, right, target)
                                   for left, right in pairs))
                pairs = itertools.product(pair_labels[first], pair_labels[second])
                affine.append(any(left is not None and right is not None and
                                  target is not None and
                                  s3(left, right, target) == 0
                                  for left, right in pairs))
            checked += 1
            true_relations += answer
            naive_missed += answer and not all(affine)
            if matches != [answer, answer] or all(matches) != answer:
                tree_mismatches.append([x1, x2, x3, x4, target,
                                        answer, matches])
    return {'n': n, 'modulus': modulus, 'a': a, 'b': b,
            'group_order': len(points), 'x_count': len(xs),
            'pair_checks': len(pair_points), 'checked': checked,
            'true_relations': true_relations,
            'naive_affine_only_missed': naive_missed,
            'pair_mismatches': pair_mismatches,
            'tree_mismatches': tree_mismatches,
            'matches_original': {
                key: value == expected[key]
                for key, value in [('group_order', len(points)),
                                   ('x_count', len(xs)),
                                   ('pair_checks', len(pair_points)),
                                   ('checked', checked),
                                   ('true_relations', true_relations),
                                   ('naive_affine_only_missed', naive_missed)]}}


def main():
    out = ROOT / 'sage-replay.json'
    if out.exists():
        raise FileExistsError(out)
    original_path = ROOT / 'result.json'
    original = json.loads(original_path.read_text())
    assert len(original['cells']) == len(CELLS)
    rows = []
    for cell, expected in zip(CELLS, original['cells']):
        assert cell == tuple(expected[k] for k in ('n', 'modulus', 'a', 'b'))
        row = replay_cell(*cell, expected)
        rows.append(row)
        print(f"n={cell[0]} a={cell[2]} b={cell[3]}: "
              f"{row['checked']} cells, {len(row['pair_mismatches'])} pair, "
              f"{len(row['tree_mismatches'])} tree mismatches")
    result = {'schema': 'binary-crossed-s3-sage-replay/v1', 'cells': rows,
              'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'original_sha256': hashlib.sha256(original_path.read_bytes()).hexdigest(),
              'plan_sha256': hashlib.sha256((ROOT/'SAGE-REPLAY-PLAN.md').read_bytes()).hexdigest()}
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    assert all(not row['pair_mismatches'] and not row['tree_mismatches'] and
               all(row['matches_original'].values()) for row in rows)


if __name__ == '__main__':
    main()
