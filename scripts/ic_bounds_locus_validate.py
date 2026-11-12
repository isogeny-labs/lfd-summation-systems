"""Independent re-check of exceptional-locus closures in the ORDINARY ring (bounds item 2).

Uses experiments/ic-curve-validation/code/lfd_engine.py (written from the Lean
definition, imports no repository closure engine): inputs are the 2k chain
equations and the generator in the full polynomial ring F_p[u_(i,j)], all
monomials of degree <= D.  For each requested record the ordinary dimensions
dim V_D and dim(V_D ∩ P_<=D-1) are compared with the chain-quotient values of
scripts/ic_bounds_locus.py through the identity

    ordinary = C(2k+D, D) - |M_<=D| + quotient     (D >= 2).

Usage
  python3 scripts/ic_bounds_locus_validate.py --data FILE.jsonl --targets t1,t2 --degrees 3,4 --out OUT.json
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'experiments/ic-curve-validation/code'))
import lfd_engine as L  # noqa: E402


def squarefree_count(r, D):
    return sum(math.comb(r, d) for d in range(min(r, D) + 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', required=True)
    ap.add_argument('--targets', required=True)
    ap.add_argument('--degrees', default='3,4')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    want = {int(t) for t in args.targets.split(',')}
    degrees = [int(d) for d in args.degrees.split(',')]
    results = []
    for line in open(args.data):
        rec = json.loads(line)
        if rec['t'] not in want:
            continue
        p, k = rec['prime'], rec['k']
        r = 2 * k
        inputs, G = L.prime_inputs(rec['generator'], k, p)
        rows = {row['D']: row for row in rec['rows']}
        for D in degrees:
            t0 = time.time()
            res = L.closure(inputs, r, D, p)
            q = rows[D]
            offset = math.comb(r + D, D) - squarefree_count(r, D)
            offset_low = math.comb(r + D - 1, D - 1) - squarefree_count(r, D - 1)
            item = dict(prime=p, k=k, a=rec['a'], b=rec['b'], t=rec['t'], D=D,
                        ordinary_dim=res['dim'], ordinary_low=res['dim_low'], ordinary_one=res['contains_one'],
                        quotient_dim=q['dim'], quotient_low=q['low'],
                        predicted_ordinary_dim=offset + q['dim'], predicted_ordinary_low=offset_low + q['low'],
                        seconds=round(time.time() - t0, 1))
            item['agree'] = (item['ordinary_dim'] == item['predicted_ordinary_dim']
                             and item['ordinary_low'] == item['predicted_ordinary_low']
                             and item['ordinary_one'] == q['one'])
            results.append(item)
            print(json.dumps(item), flush=True)
    prev = json.loads(Path(args.out).read_text()) if Path(args.out).exists() else []
    Path(args.out).write_text(json.dumps(prev + results, indent=1) + '\n')
    assert all(x['agree'] for x in results)


if __name__ == '__main__':
    main()
