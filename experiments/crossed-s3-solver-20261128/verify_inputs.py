#!/usr/bin/env python3
"""Check generated msolve equations on all selected direct relations."""

import itertools
import json
from pathlib import Path

import generate

HERE = Path(__file__).resolve().parent


def equations(path):
    lines = path.read_text().splitlines()
    names = lines[0].split(",")
    assert lines[1] == "2"
    expressions = [compile(s.strip().replace("^", "**"), str(path), "eval")
                   for s in "\n".join(lines[2:]).split(",") if s.strip()]
    return names, expressions


def satisfies(names, expressions, bits):
    values = {name: bit for name, bit in zip(names, bits)}
    return all(eval(expr, {"__builtins__": {}}, values) % 2 == 0
               for expr in expressions)


def one_witness(quad, target, pairing, labels, b, modulus):
    (i, j), (k, l) = pairing
    return next(((u, v) for u in labels[quad[i], quad[j]]
                 for v in labels[quad[k], quad[l]]
                 if generate.ec.root_match(u, v, target, b, modulus)), None)


def bits_of(quad, witnesses, k, n):
    bits = [(x >> i) & 1 for x in quad for i in range(k)]
    for pair in witnesses:
        for label in pair:
            bits += [(0 if label is None else label >> i & 1)
                     for i in range(n)]
            bits.append(int(label is None))
    return bits


def main():
    records = json.loads((HERE / "manifest.json").read_text())
    checks = 0
    direct_truth_checks = 0
    for n, modulus, k in generate.GRID:
        for a in (0, 1):
            base, _, labels, selected = generate.point_data(n, modulus, a, 1, k)
            for target, expected in selected:
                rows = {r["arm"]: r for r in records
                        if (r["n"], r["a"], r["target"]) == (n, a, target)}
                compiled = {arm: equations(HERE / row["input"])
                            for arm, row in rows.items()}
                pairings = (((0, 1), (2, 3)), ((0, 2), (1, 3)))
                for quad in expected:
                    witnesses = [one_witness(quad, target, pairing, labels, 1,
                                             modulus) for pairing in pairings]
                    assert all(witness is not None for witness in witnesses)
                    for arm in generate.ARMS:
                        names, expressions = compiled[arm]
                        pairs = witnesses[:1] if arm == "tree" else (
                            witnesses if arm == "crossed" else [])
                        bits = bits_of(quad, pairs, k, n)
                        assert len(bits) == len(names)
                        assert satisfies(names, expressions, bits), (
                            n, a, target, arm, quad, witnesses)
                        checks += 1
                names, expressions = compiled["direct"]
                for quad in itertools.product(range(1 << k), repeat=4):
                    bits = bits_of(quad, [], k, n)
                    assert satisfies(names, expressions, bits) == (quad in expected)
                    direct_truth_checks += 1
    result = {"positive_witness_equation_checks": checks,
              "direct_truth_table_checks": direct_truth_checks,
              "mismatches": 0}
    (HERE / "input-verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
