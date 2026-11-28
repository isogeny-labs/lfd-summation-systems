#!/usr/bin/env python3
"""Generate frozen, matched O-aware S3 systems and a direct Boolean eliminant."""

import hashlib
import importlib.util
import itertools
import json
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "binary-crossed-s3-witness-20261128" / "run.py"
spec = importlib.util.spec_from_file_location("binary_witness", SOURCE)
ec = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ec)
GRID = ((3, 0b1011, 2), (4, 0b10011, 3), (5, 0b100101, 3))
ARMS = ("tree", "crossed", "direct")


def xor(*polys):
    out = set()
    for poly in polys:
        out.symmetric_difference_update(poly)
    return out


def product(left, right):
    out = set()
    for a in left:
        for b in right:
            term = a | b  # Boolean quotient: v^2 = v.
            if term in out:
                out.remove(term)
            else:
                out.add(term)
    return out


ONE = {0}
ZERO = set()


def vpoly(i):
    return {1 << i}


def const_field(value, n):
    return [set(ONE) if value & (1 << i) else set() for i in range(n)]


def field_xor(*values):
    return [xor(*(v[i] for v in values)) for i in range(len(values[0]))]


def field_mul(left, right, modulus):
    n = len(left)
    cols = [set() for _ in range(2 * n - 1)]
    for i in range(n):
        for j in range(n):
            cols[i + j] = xor(cols[i + j], product(left[i], right[j]))
    for deg in range(2 * n - 2, n - 1, -1):
        top = cols[deg]
        for bit in range(n):
            if modulus & (1 << bit):
                cols[deg - n + bit] = xor(cols[deg - n + bit], top)
    return cols[:n]


def s3_poly(x, y, z, b, modulus):
    xy = field_mul(x, y, modulus)
    xz = field_mul(x, z, modulus)
    yz = field_mul(y, z, modulus)
    t = field_xor(xy, xz, yz)
    return field_xor(field_mul(t, t, modulus),
                     field_mul(xy, z, modulus), const_field(b, len(x)))


def check_s3_poly(n, modulus, b):
    variables = [[vpoly(j * n + i) for i in range(n)] for j in range(3)]
    polys = s3_poly(*variables, b, modulus)
    rng = random.Random(0x533 + n)
    samples = [(0, 0, 0), (1, 1, 1)]
    samples += [tuple(rng.randrange(1 << n) for _ in range(3))
                for _ in range(256)]
    for x, y, z in samples:
        bits = x | (y << n) | (z << (2 * n))
        actual = sum(evaluate(poly, bits) << i for i, poly in enumerate(polys))
        assert actual == ec.s3(x, y, z, b, modulus), (n, x, y, z)


def evaluate(poly, bits):
    return sum(1 for mask in poly if mask & bits == mask) & 1


def anf_of_truth(values, nvars):
    coeff = list(values)
    assert len(coeff) == 1 << nvars
    for bit in range(nvars):
        for mask in range(1 << nvars):
            if mask & (1 << bit):
                coeff[mask] ^= coeff[mask ^ (1 << bit)]
    return {mask for mask, value in enumerate(coeff) if value}


def format_poly(poly, names):
    if not poly:
        return "0"
    terms = []
    for mask in sorted(poly, key=lambda x: (x.bit_count(), x)):
        terms.append("*".join(names[i] for i in range(len(names))
                              if mask & (1 << i)) or "1")
    return "+".join(terms)


def point_data(n, modulus, a, b, k):
    q = 1 << n
    points = [(x, y) for x in range(q) for y in range(q)
              if ec.on_curve((x, y), a, b, modulus)]
    xs = sorted({x for x, _ in points})
    base = [x for x in xs if x < (1 << k)]
    lifts = {x: [p for p in points if p[0] == x] for x in xs}
    pair = {(x, y): {ec.add(p, r, a, modulus)
                     for p, r in itertools.product(lifts[x], lifts[y])}
            for x, y in itertools.product(base, repeat=2)}
    pair_labels = {(x, y): {ec.label(p) for p in pair[x, y]}
                   for x, y in pair}
    targets = []
    for target in xs:
        relation = {quad for quad in itertools.product(base, repeat=4)
                    if any(ec.label(ec.add(u, v, a, modulus)) == target
                           for u in pair[quad[0], quad[1]]
                           for v in pair[quad[2], quad[3]])}
        targets.append((target, relation))
    selected = []
    for desired in (True, False):
        hit = next(((target, relation) for target, relation in targets
                    if bool(relation) == desired), None)
        if hit is not None and hit[0] not in [x[0] for x in selected]:
            selected.append(hit)
    return base, xs, pair_labels, selected


def check_projection(base, target, expected, pair_labels, b, modulus):
    for quad in itertools.product(base, repeat=4):
        pairings = (((0, 1), (2, 3)), ((0, 2), (1, 3)))
        roots = []
        for (i, j), (k, l) in pairings:
            roots.append(any(ec.root_match(u, v, target, b, modulus)
                             for u in pair_labels[quad[i], quad[j]]
                             for v in pair_labels[quad[k], quad[l]]))
        if roots != [quad in expected] * 2:
            raise AssertionError((quad, target, roots, quad in expected))


def build_system(n, k, modulus, b, base, target, expected, arm):
    leaf_count = 4 * k
    names = [f"v{i}" for i in range(leaf_count)]
    leaves = [[vpoly(j * k + i) if i < k else set()
               for i in range(n)] for j in range(4)]
    eqs = []
    allowed = [0 if x in base else 1 for x in range(1 << k)]
    membership = anf_of_truth(allowed, k)
    for j in range(4):
        eqs.append({sum(1 << (j * k + i) for i in range(k)
                         if mask & (1 << i)) for mask in membership})
    if arm == "direct":
        truth = []
        for bits in range(1 << leaf_count):
            quad = tuple((bits >> (j * k)) & ((1 << k) - 1)
                         for j in range(4))
            truth.append(0 if quad in expected else 1)
        eqs.append(anf_of_truth(truth, leaf_count))
    else:
        pairings = (((0, 1), (2, 3)),)
        if arm == "crossed":
            pairings += (((0, 2), (1, 3)),)
        for pairing in pairings:
            nodes = []
            for i, j in pairing:
                start = len(names)
                names += [f"v{start+t}" for t in range(n + 1)]
                x = [vpoly(start + t) for t in range(n)]
                o = vpoly(start + n)
                finite = xor(ONE, o)
                for poly in s3_poly(leaves[i], leaves[j], x, b, modulus):
                    eqs.append(product(finite, poly))
                for bit in range(n):
                    eqs.append(product(o, xor(leaves[i][bit], leaves[j][bit])))
                    eqs.append(product(o, x[bit]))
                nodes.append((x, o))
            (u, ou), (v, ov) = nodes
            fu, fv = xor(ONE, ou), xor(ONE, ov)
            both_finite = product(fu, fv)
            t = const_field(target, n)
            for poly in s3_poly(u, v, t, b, modulus):
                eqs.append(product(both_finite, poly))
            for bit in range(n):
                eqs.append(product(product(ou, fv), xor(v[bit], t[bit])))
                eqs.append(product(product(ov, fu), xor(u[bit], t[bit])))
            eqs.append(product(ou, ov))
    # Boolean field equations must retain the square in the ordinary ring.
    polys = [format_poly(poly, names) for poly in eqs if poly]
    polys += [f"{name}^2+{name}" for name in names]
    return ",".join(names) + "\n2\n" + ",\n".join(polys) + "\n", len(names), len(polys)


def main():
    input_dir = HERE / "inputs"
    input_dir.mkdir(exist_ok=True)
    records = []
    for n, modulus, k in GRID:
        check_s3_poly(n, modulus, 1)
        for a in (0, 1):
            b = 1
            base, xs, pair_labels, selected = point_data(n, modulus, a, b, k)
            for target, expected in selected:
                check_projection(base, target, expected, pair_labels, b, modulus)
                for arm in ARMS:
                    content, nvars, neqs = build_system(
                        n, k, modulus, b, base, target, expected, arm)
                    path = input_dir / f"n{n}_a{a}_b{b}_t{target}_{arm}.ms"
                    path.write_text(content)
                    records.append({"n": n, "modulus": modulus, "k": k,
                                    "a": a, "b": b, "target": target,
                                    "base": base, "all_rational_x": xs,
                                    "relations": len(expected), "arm": arm,
                                    "variables": nvars, "equations": neqs,
                                    "input": str(path.relative_to(HERE)),
                                    "sha256": hashlib.sha256(content.encode()).hexdigest()})
    (HERE / "manifest.json").write_text(json.dumps(records, indent=2) + "\n")
    print(json.dumps({"systems": len(records), "queries": len(records) // 3,
                      "arms": list(ARMS)}, sort_keys=True))


if __name__ == "__main__":
    main()
