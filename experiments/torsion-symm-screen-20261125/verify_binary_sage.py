#!/usr/bin/env sage -python
"""Independent dense-space check of the characteristic-two LFD thresholds."""
import json
import hashlib
import sys
from pathlib import Path

from sage.all import GF, VectorSpace, vector, PolynomialRing, Matrix

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from binary_edwards_degree import field, image_generators, product


def contains_one(generators, K, D):
    monoms = [(i, degree-i) for degree in range(D, -1, -1)
              for i in range(degree, -1, -1)]
    index = {m: j for j, m in enumerate(monoms)}
    space = VectorSpace(K, len(monoms))

    def row(poly):
        result = vector(K, len(monoms))
        for monom, coeff in poly.dict().items():
            key = tuple(monom)
            if sum(key) <= D:
                result[index[key]] = coeff
        return result

    basis = space.subspace([row(g) for g in generators if g and g.total_degree() <= D])
    while True:
        new_rows = []
        for v in basis.basis():
            if any(v[index[m]] for m in monoms if sum(m) == D):
                continue
            for axis in (0, 1):
                w = vector(K, len(monoms))
                for m in monoms:
                    if v[index[m]]:
                        shifted = (m[0]+(axis == 0), m[1]+(axis == 1))
                        w[index[shifted]] = v[index[m]]
                new_rows.append(w)
        next_basis = space.subspace(list(basis.basis())+new_rows)
        if next_basis.dimension() == basis.dimension():
            break
        basis = next_basis
    one = vector(K, len(monoms))
    one[index[(0, 0)]] = 1
    return one in basis, basis.dimension()


def main():
    data = json.loads((HERE/"binary-edwards-yield.json").read_text())
    rank_checks = []
    for source in data["cells"]:
        for arm in ("baseline", "colored"):
            rows = source["relation_rows"][arm]
            rank = int(Matrix(GF(source["r"]), rows).rank()) if rows else 0
            assert rank == source["final_rank"][arm]
            rank_checks.append({"n": source["n"], "arm": arm, "sage_rank": rank,
                                "stored_rank": source["final_rank"][arm]})
    c = next(x for x in data["cells"] if x["n"] == 9)
    K, cv = field(9, c["modulus"])
    a = sorted({t for col in c["columns"] if col["color"] == "A" for t in col["t_pair"]})
    b = sorted({t for col in c["columns"] if col["color"] == "B" for t in col["t_pair"]})
    w = sorted(set(a) | set(b))
    target = next(x for x in c["per_target"] if x["baseline_rows"] == 0)
    tr = cv(target["target"][0] ^ target["target"][1])
    cases = []
    for arm, (left, right) in (("baseline", (w, w)), ("colored", (a, b))):
        R = PolynomialRing(K, names=("x", "y"), order="degrevlex")
        x, y = R.gens()
        s = (x*x+x)*(y*y+y)
        direct = [product([cv(v) for v in left], x, K),
                  product([cv(v) for v in right], y, K),
                  x**K.cardinality()-x, y**K.cardinality()-y,
                  (x+y)**2+s*(tr*tr+tr)+tr*tr]
        I, image, _ = image_generators(K, [cv(v) for v in left],
                                       [cv(v) for v in right])
        e, q = I.gens()
        invariant = image+[e*e+q*(tr*tr+tr)+tr*tr]
        for form, generators, D in (("direct", direct, 14 if arm == "baseline" else 8),
                                     ("invariant", invariant, 8 if arm == "baseline" else 5)):
            below, rank_below = contains_one(generators, K, D-1)
            at, rank_at = contains_one(generators, K, D)
            gb = generators[0].parent().ideal(generators).groebner_basis()
            unit = len(gb) == 1 and gb[0] == 1
            item = {"arm": arm, "form": form, "target_scalar": target["scalar"],
                    "degree": D, "below_contains_one": bool(below), "at_contains_one": bool(at),
                    "below_rank": int(rank_below), "at_rank": int(rank_at),
                    "groebner_unit": bool(unit)}
            cases.append(item)
            print(json.dumps(item, sort_keys=True), flush=True)
            assert not below and at and unit
    (HERE/"binary-sage-verification.json").write_text(json.dumps({"schema": "torsion-symm-binary-sage-verification/v1",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "yield_json_sha256": hashlib.sha256((HERE/"binary-edwards-yield.json").read_bytes()).hexdigest(),
        "rank_checks": rank_checks, "cases": cases}, sort_keys=True, indent=2)+"\n")


if __name__ == "__main__":
    main()
