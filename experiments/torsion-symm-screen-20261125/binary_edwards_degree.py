#!/usr/bin/env sage -python
"""Exact scalar-coordinate unit-ideal LFD for matched binary Edwards cells."""
from __future__ import annotations

import hashlib
import argparse
import json
import math
from pathlib import Path
import resource
import sys
import time

from sage.all import GF, PolynomialRing

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
CAP = 12
SECONDS = 60
MAX_COLUMNS = 2000


def field(n, modulus):
    R = PolynomialRing(GF(2), "z")
    z = R.gen()
    poly = sum(((modulus >> i) & 1) * z**i for i in range(n+1))
    K = GF(2**n, "a", modulus=poly)
    a = K.gen()
    powers = [a**i for i in range(n)]
    return K, lambda value: sum((powers[i] for i in range(n) if value >> i & 1), K(0))


def product(values, variable, K):
    out = variable.parent()(1)
    for value in values:
        out *= variable - K(value)
    return out


def polynomial_degree(poly):
    return max((sum(exps) for exps in poly.dict()), default=-1)


def closure_unit_degree(generators, cap=CAP, seconds=SECONDS):
    """Iterated V_D closure; first D with 1 in the span is exact unit LFD."""
    started = time.perf_counter()
    profile = []
    for D in range(cap+1):
        if math.comb(D+2, 2) > MAX_COLUMNS:
            return {"lfd": None, "status": "column-cap", "profile": profile}
        basis = {}
        queue = [{tuple(m): c for m, c in g.dict().items()}
                 for g in generators if polynomial_degree(g) <= D]
        at = 0
        while at < len(queue):
            if time.perf_counter() - started > seconds:
                return {"lfd": None, "status": "time-cap", "profile": profile}
            row = dict(queue[at])
            at += 1
            while row:
                lead = max(row, key=lambda m: (sum(m), m[0]))
                if lead not in basis:
                    if sum(lead) > D:
                        raise AssertionError(f"closure monomial exceeds cap: D={D}, lead={lead}")
                    inv = row[lead]**-1
                    row = {m: c*inv for m, c in row.items()}
                    basis[lead] = row
                    if sum(lead) < D:
                        for index in (0, 1):
                            shifted = {}
                            for m, c in row.items():
                                q = (m[0]+(index == 0), m[1]+(index == 1))
                                if sum(q) <= D:
                                    shifted[q] = c
                            queue.append(shifted)
                    break
                factor = row[lead]
                for m, c in basis[lead].items():
                    new = row.get(m, 0) - factor*c
                    if new:
                        row[m] = new
                    else:
                        row.pop(m, None)
        profile.append({"degree": D, "rank": len(basis), "columns": math.comb(D+2, 2)})
        if len(basis) > math.comb(D+2, 2):
            raise AssertionError(f"closure rank exceeds columns: D={D}, keys={list(basis)}")
        if (0, 0) in basis:
            return {"lfd": D, "status": "exact", "profile": profile,
                    "elapsed_s": time.perf_counter()-started}
    return {"lfd": None, "status": "degree-cap", "profile": profile,
            "elapsed_s": time.perf_counter()-started}


def image_generators(K, values_a, values_b):
    R = PolynomialRing(K, names=("e", "s"), order="degrevlex")
    e, s = R.gens()
    image = sorted({(u+v, (u*u+u)*(v*v+v)) for u in values_a for v in values_b},
                   key=lambda p: (str(p[0]), str(p[1])))
    if not image:
        return R, [R(1)], 0
    J = R.ideal([e-image[0][0], s-image[0][1]])
    for ev, sv in image[1:]:
        J = J.intersection(R.ideal([e-ev, s-sv]))
    return R, list(J.groebner_basis())+[e**K.cardinality()-e, s**K.cardinality()-s], len(image)


def cell(source, cap=CAP):
    n = source["n"]
    K, from_int = field(n, source["modulus"])
    q = 1 << n
    # Verify the Sage field has exactly the repository's polynomial basis.
    from csd.ic_dreg import binary_field
    F = binary_field(n)
    for i in range(q):
        for j in range(q):
            if from_int(F.mul(i, j)) != from_int(i)*from_int(j):
                raise AssertionError("field modulus/representation mismatch")
    ta = sorted({t for col in source["columns"] if col["color"] == "A" for t in col["t_pair"]})
    tb = sorted({t for col in source["columns"] if col["color"] == "B" for t in col["t_pair"]})
    tw = sorted(set(ta) | set(tb))
    bases = {"baseline": (tw, tw), "colored": (ta, tb)}
    targets = [x for x in source["per_target"] if x["baseline_rows"] == 0][:4]
    output = {"n": n, "l": source["l"], "field_order": q,
              "selected_target_scalars": [x["scalar"] for x in targets], "bases": {}, "targets": []}
    for arm, (left, right) in bases.items():
        start = time.perf_counter()
        cpu_start = resource.getrusage(resource.RUSAGE_SELF)
        L, U = [from_int(v) for v in left], [from_int(v) for v in right]
        direct = PolynomialRing(K, names=("x", "y"), order="degrevlex")
        x, y = direct.gens()
        direct_base = [product(L, x, K), product(U, y, K), x**q-x, y**q-y]
        invariant, invariant_base, image_size = image_generators(K, L, U)
        cpu_end = resource.getrusage(resource.RUSAGE_SELF)
        output["bases"][arm] = {"left_t": left, "right_t": right,
                               "image_size": image_size,
                               "direct_generator_degrees": [polynomial_degree(g) for g in direct_base],
                               "invariant_generator_degrees": [polynomial_degree(g) for g in invariant_base],
                               "setup_wall_s": time.perf_counter()-start,
                               "setup_cpu_s": (cpu_end.ru_utime+cpu_end.ru_stime)-
                                              (cpu_start.ru_utime+cpu_start.ru_stime)}
        for target in targets:
            tr = from_int(target["target"][0] ^ target["target"][1])
            a = x*y
            e = x+y
            s = (x*x+x)*(y*y+y)
            direct_f3 = e*e+s*(tr*tr+tr)+tr*tr
            ev, sv = invariant.gens()
            quotient_f3 = ev*ev+sv*(tr*tr+tr)+tr*tr
            for form, gens in (("direct", direct_base+[direct_f3]),
                               ("invariant", invariant_base+[quotient_f3])):
                started = time.perf_counter()
                cpu0 = resource.getrusage(resource.RUSAGE_SELF)
                exact = closure_unit_degree(gens, cap=cap)
                cpu1 = resource.getrusage(resource.RUSAGE_SELF)
                output["targets"].append({"scalar": target["scalar"], "arm": arm,
                                          "form": form,
                                          "generator_degrees": [polynomial_degree(g) for g in gens],
                                          "closure": exact,
                                          "wall_s": time.perf_counter()-started,
                                          "cpu_s": (cpu1.ru_utime+cpu1.ru_stime)-
                                                   (cpu0.ru_utime+cpu0.ru_stime)})
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cap", type=int, default=CAP)
    parser.add_argument("--degrees", type=int, nargs="+", default=[5, 7, 9])
    parser.add_argument("--out", type=Path, default=HERE/"binary-edwards-degree.json")
    args = parser.parse_args()
    yield_data = json.loads((HERE/"binary-edwards-yield.json").read_text())
    source = Path(__file__).read_bytes()
    output = {"schema": "torsion-symm-binary-degree/v1",
              "protocol": "PROTOCOL_BINARY_20261125.md", "source_sha256": hashlib.sha256(source).hexdigest(),
              "yield_source_sha256": yield_data["source_sha256"],
              "degree_cap": args.cap,
              "cells": [cell(x, cap=args.cap) for x in yield_data["cells"] if x["n"] in args.degrees]}
    args.out.write_text(json.dumps(output, sort_keys=True, indent=2)+"\n")
    for item in output["cells"]:
        print(item["n"], [(x["scalar"], x["arm"], x["form"], x["closure"]["lfd"],
                           x["closure"]["status"]) for x in item["targets"]])


if __name__ == "__main__":
    main()
