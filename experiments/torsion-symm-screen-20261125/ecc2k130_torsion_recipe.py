#!/usr/bin/env python3
"""Exact type-II ONB reciprocal-closure screen for the fixed ECC2K-130 recipe."""
from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from csd.gf2_ghs import GF2n  # noqa: E402

HERE = Path(__file__).resolve().parent
PARAMS = Path("/Volumes/SSD990-2/crypto/docs/ic/params/ecc2k130-fixed.json")
M = 131
N = 2 * M + 1
MASK = (1 << N) - 1


def from_coords(value):
    out = 0
    for i in range(1, M + 1):
        if value >> (i - 1) & 1:
            out |= (1 << i) | (1 << (N - i))
    return out


def normalize(value):
    return value ^ MASK if value & 1 else value


def to_coords(value):
    return (normalize(value) >> 1) & ((1 << M) - 1)


def multiply(a, b):
    out = 0
    while b:
        low = b & -b
        shift = low.bit_length() - 1
        out ^= a if shift == 0 else ((a << shift) | (a >> (N - shift))) & MASK
        b ^= low
    return normalize(out)


ONE = from_coords((1 << M) - 1)


def power(base, exponent):
    value = ONE
    while exponent:
        if exponent & 1:
            value = multiply(value, base)
        base = multiply(base, base)
        exponent >>= 1
    return value


def batch_inverse(values):
    prefixes = []
    product = ONE
    for value in values:
        assert value
        prefixes.append(product)
        product = multiply(product, value)
    inverse = power(product, (1 << M) - 2)
    out = [0] * len(values)
    for i in range(len(values) - 1, -1, -1):
        out[i] = multiply(inverse, prefixes[i])
        inverse = multiply(inverse, values[i])
    assert inverse == ONE
    return out


def polynomial_to_onb(value, images):
    out = 0
    while value:
        low = value & -value
        out ^= images[low.bit_length() - 1]
        value ^= low
    return out


def inverse_linear_map(images):
    pivots = {}
    for i, image in enumerate(images):
        x, coeff = image, 1 << i
        while x:
            j = x.bit_length() - 1
            if j not in pivots:
                pivots[j] = (x, coeff)
                break
            y, d = pivots[j]
            x ^= y
            coeff ^= d
        else:
            raise AssertionError("polynomial-to-ONB map is singular")
    assert len(pivots) == M

    def convert(x):
        coeff = 0
        while x:
            j = x.bit_length() - 1
            y, d = pivots[j]
            x ^= y
            coeff ^= d
        return coeff

    return convert


def main():
    raw = PARAMS.read_bytes()
    params = json.loads(raw)
    assert hashlib.sha256(raw).hexdigest() == "f41972681d260ef4eeb2fea7441f31859f4a121fdd5f2ecc4ff3dee54082aeae"
    curve = params["curve"]
    assert curve["degree"] == M and curve["cofactor"] == 4
    assert params["factor_base"] == {"kind": "hamming-weight", "weight": 2}
    polynomial = int(curve["polynomial"], 16)
    F = GF2n(M, polynomial)
    root = from_coords(int(curve["onb_root"], 16))
    powers = []
    value = ONE
    for i in range(M + 1):
        powers.append(to_coords(value))
        value = multiply(value, root)
    polynomial_value = 0
    for i in range(M + 1):
        if polynomial >> i & 1:
            polynomial_value ^= powers[i]
    assert polynomial_value == 0
    images = powers[:M]
    to_polynomial = inverse_linear_map(images)
    assert all(to_polynomial(image) == 1 << i for i, image in enumerate(images))

    coordinates = [1 << i for i in range(M)]
    coordinates += [(1 << i) | (1 << j) for i in range(M) for j in range(i + 1, M)]
    assert len(coordinates) == M + M * (M - 1) // 2 == 8646
    elements = [from_coords(x) for x in coordinates]
    assert all(multiply(ONE, x) == x for x in elements)
    inverses = batch_inverse(elements)
    counts = {1: Counter(), 2: Counter()}
    first_outside = None
    sample_checks = []
    for i, (coord, value, inverse) in enumerate(zip(coordinates, elements, inverses)):
        assert multiply(value, inverse) == ONE
        reciprocal = to_coords(inverse)
        weight = reciprocal.bit_count()
        counts[coord.bit_count()][weight] += 1
        if weight > 2 and first_outside is None:
            first_outside = {"input_onb": hex(coord), "input_weight": coord.bit_count(),
                             "reciprocal_onb": hex(reciprocal), "reciprocal_weight": weight}
        if i < 16 or i in (M, len(elements) - 1):
            a, b = to_polynomial(coord), to_polynomial(reciprocal)
            assert F.mul(a, b) == 1 and F.inv(a) == b
            sample_checks.append(i)
    assert first_outside is not None
    x = 1
    for order in range(1, 130):
        x = (x * 2) % 131
        assert x != 1
    assert (x * 2) % 131 == 1

    result = {"schema": "torsion-symm-ecc2k130-recipe-v1",
              "protocol": "PROTOCOL_ECC2K130_T2_COMPAT_20261125.md",
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "params_sha256": hashlib.sha256(raw).hexdigest(),
              "field_degree": M, "onb_order": N, "frobenius_order_mod_131": 130,
              "tested_coordinates": len(coordinates),
              "counts_by_input_weight": {str(w): dict(sorted(c.items())) for w, c in counts.items()},
              "reciprocal_weight_at_most_two": sum(c[w] for c in counts.values() for w in (1, 2)),
              "first_outside": first_outside,
              "independent_polynomial_basis_checks": sample_checks}
    (HERE / "ecc2k130-torsion-recipe.json").write_text(json.dumps(result, sort_keys=True, indent=2)+"\n")
    print("recipe candidates", len(coordinates), "reciprocal weight<=2",
          result["reciprocal_weight_at_most_two"], "first outside", first_outside)


if __name__ == "__main__":
    main()
