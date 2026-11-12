#!/usr/bin/env python3
"""Verifier step 5: ordinary-ring cross-check of the chain-quotient closure.
lfd_engine.closure (full polynomial ring in 2k variables, chain equations + reduced generator,
own generator from ic_bounds_locus_verify.generator) versus the verifier's chain closure:
checks dim V_D(ordinary) = C(2k+D, D) - |M_<=D| + dim V_D(chain), and the same for the
low parts (degree <= D-1)."""
import argparse
import json
import os
import sys
import time
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "experiments", "ic-curve-validation", "code"))
import lfd_engine  # noqa: E402
from ic_bounds_locus_verify import generator, chain_closure, nmono  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prime", type=int, default=193)
    ap.add_argument("--k", type=int, default=6)
    ap.add_argument("--targets", required=True)
    ap.add_argument("--degrees", default="3,4")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    p, k = args.prime, args.k
    r = 2 * k
    with open(args.out, "a") as f:
        for tok in args.targets.split(";"):
            a, b, t = map(int, tok.split(","))
            g = generator(p, k, a, b, t)
            inputs, _ = lfd_engine.prime_inputs(sorted(g.items()), k, p)
            for D in map(int, args.degrees.split(",")):
                t0 = time.time()
                o = lfd_engine.closure(inputs, r, D, p)
                c = chain_closure(g, k, D, p)
                chain_part = comb(r + D, D) - nmono(k, D)
                chain_part_low = comb(r + D - 1, D - 1) - nmono(k, D - 1)
                rec = dict(prime=p, k=k, a=a, b=b, t=t, D=D, ordinary_dim=o["dim"], ordinary_low=o["dim_low"],
                           ordinary_one=o["contains_one"], chain_dim=c["dim"], chain_low=c["low"], chain_one=c["one"],
                           dim_ok=o["dim"] == chain_part + c["dim"], low_ok=o["dim_low"] == chain_part_low + c["low"],
                           seconds=round(time.time() - t0, 1))
                f.write(json.dumps(rec) + "\n")
                f.flush()
                print(rec, flush=True)


if __name__ == "__main__":
    main()
