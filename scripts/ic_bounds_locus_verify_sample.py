#!/usr/bin/env python3
"""Verifier step 1: read the producer's per-target rows (data/ and heldout/), re-derive with
own code the cheap facts for EVERY row (target list completeness, generator, |Z|, internal
consistency of falls/LFD with the stored dimensions), and draw the stratified sample that
ic_bounds_locus_verify.py recomputes.  Seed: random.Random(777 + p)."""
import glob
import json
import os
import random
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ic_bounds_locus_verify import generator, subgroup, is_square, nmono  # own code

ROOT = os.path.dirname(HERE)
EXP = os.path.join(ROOT, "experiments", "ic-bounds-locus")
GENERIC = {6: (13, 1, 605, 111), 7: (15, 1, 829, 115)}


def load(p, k):
    rows = []
    for sub in ("data", "heldout"):
        for fn in sorted(glob.glob(os.path.join(EXP, sub, f"p{p}-k{k}-*.jsonl"))):
            for line in open(fn):
                r = json.loads(line)
                r["_set"] = sub
                rows.append(r)
    return rows


def zcount(p, k, a, b, t, H):
    x1 = H[:, None]
    x2 = H[None, :]
    v = ((x1 - x2) ** 2 * t * t - 2 * ((x1 + x2) * (x1 * x2 + a) + 2 * b) * t + (x1 * x2 - a) ** 2 - 4 * b * (x1 + x2)) % p
    return int(np.sum(v == 0))


def main():
    out = {}
    samples = {}
    for p, k in ((193, 6), (257, 7)):
        rows = load(p, k)
        H = np.array(subgroup(p, k), dtype=np.int64)
        rep = dict(rows=len(rows), by_set={}, gen_mismatch=[], z_mismatch=[], missing_targets=[], extra_targets=[],
                   internal=[], deviating=0)
        curves = sorted({(r["a"], r["b"], r["_set"]) for r in rows})
        for a, b, s in curves:
            want = {t for t in range(p) if is_square(t ** 3 + a * t + b, p)}
            have = {r["t"] for r in rows if (r["a"], r["b"]) == (a, b)}
            if want - have:
                rep["missing_targets"].append([a, b, sorted(want - have)])
            if have - want:
                rep["extra_targets"].append([a, b, sorted(have - want)])
            rep["by_set"].setdefault(s, 0)
            rep["by_set"][s] += len(have)
        dev, gen = [], []
        for r in rows:
            g = generator(p, k, r["a"], r["b"], r["t"])
            if sorted([m, c] for m, c in g.items()) != sorted(r["generator"]):
                rep["gen_mismatch"].append([r["a"], r["b"], r["t"]])
            z = zcount(p, k, r["a"], r["b"], r["t"], H)
            if z != r["Z"]:
                rep["z_mismatch"].append([r["a"], r["b"], r["t"], z, r["Z"]])
            d = {x["D"]: x for x in r["rows"]}
            # internal consistency: falls from stored dims
            falls = [D for D in d if D - 1 in d and d[D]["low"] > d[D - 1]["dim"]]
            if falls != r["falls"]:
                rep["internal"].append(["falls", r["a"], r["b"], r["t"], falls, r["falls"]])
            if r.get("exact_lfd") is not None:
                Dm = max(d)
                if not (d[Dm]["one"] or d[Dm]["codim"] == r["Z"]):
                    rep["internal"].append(["cert", r["a"], r["b"], r["t"]])
                if r["exact_lfd"] != (max(falls) if falls else None):
                    rep["internal"].append(["lfd-vs-fall", r["a"], r["b"], r["t"], r["exact_lfd"], falls])
            sig = (d[3]["dim"], d[3]["low"], d[4]["dim"], d[4]["low"])
            (dev if sig != GENERIC[k] else gen).append((r["a"], r["b"], r["t"], r["_set"], sig))
        rep["deviating"] = len(dev)
        rng = random.Random(777 + p)
        Hs = set(int(h) for h in H)
        # stratified generic sample: per curve, one t in H and one t not in H (when present), + 2 more random
        gs = []
        for a, b, s in curves:
            pool = [x for x in gen if (x[0], x[1]) == (a, b)]
            inH = [x for x in pool if x[2] in Hs]
            outH = [x for x in pool if x[2] not in Hs]
            for P in (inH, outH):
                if P:
                    gs.append(rng.choice(P))
            rest = [x for x in pool if x not in gs]
            gs += rng.sample(rest, min(2, len(rest)))
        samples[f"p{p}-k{k}"] = dict(deviating=dev, generic=gs)
        out[f"p{p}-k{k}"] = rep
    json.dump(out, open(os.path.join(EXP, "verify", "rowchecks.json"), "w"), indent=1)
    json.dump(samples, open(os.path.join(EXP, "verify", "sample.json"), "w"), indent=1)
    for key, rep in out.items():
        print(key, {k: (v if not isinstance(v, list) else len(v)) for k, v in rep.items()})
        print("  sample:", len(samples[key]["deviating"]), "deviating +", len(samples[key]["generic"]), "generic")


if __name__ == "__main__":
    main()
