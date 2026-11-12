#!/usr/bin/env python3
"""Compare the producer's exported RREF bases (experiments/ic-bounds-fall/p40961-s7/k*-D4.npz,
k*-D5-low.npz) with bases recomputed by the independent closure of
scripts/ic_bounds_fall_verify.py, as row spaces, and check the recorded file hashes.
Usage: python3 scripts/ic_bounds_fall_verify_npz.py --k 6 7 8 9 10 --k5 6 7
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import ic_bounds_fall_verify as V  # noqa: E402

DIR = ROOT / 'experiments/ic-bounds-fall/p40961-s7'


def load_as_rows(path, G):
    z = np.load(path)
    masks = z['masks'].astype(np.int64)
    k = G.k
    a1, a2 = masks & ((1 << k) - 1), masks >> k
    cols = G.idx[a1 * G.N + a2]
    assert np.all(cols >= 0)
    R = np.zeros((z['basis'].shape[0], G.n))
    R[:, cols] = z['basis'].astype(np.float64)
    return R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--k', type=int, nargs='*', default=[6, 7, 8, 9, 10])
    ap.add_argument('--k5', type=int, nargs='*', default=[])
    ap.add_argument('--out', type=Path, default=ROOT / 'experiments/ic-bounds-fall/verify/npz.json')
    args = ap.parse_args()
    I, g = V.inst_g('p40961')
    p = I['p']
    res = {}
    for k in args.k:
        meta = json.loads((DIR / f'k{k}-D4.json').read_text())
        fsha = hashlib.sha256((DIR / f'k{k}-D4.npz').read_bytes()).hexdigest()
        r = V.closure(g, k, 4, p, cap=4 * k * k + 172 * k + 200, done_cap=600, log=lambda m: None)
        P = load_as_rows(DIR / f'k{k}-D4.npz', r['G'])
        res[f'D4_k{k}'] = dict(file_sha_matches=(fsha == meta['basis_file']['sha256']),
                               compare_producer_vs_own=V.compare(P, r['E'], p))
        print(k, res[f'D4_k{k}'], flush=True)
    for k in args.k5:
        meta = json.loads((DIR / f'k{k}-D5.json').read_text())
        fsha = hashlib.sha256((DIR / f'k{k}-D5-low.npz').read_bytes()).hexdigest()
        r = V.closure(g, k, 5, p, log=lambda m: None)
        own_low = r['E'][r['pw'] <= 4]
        P = load_as_rows(DIR / f'k{k}-D5-low.npz', r['G'])
        res[f'D5low_k{k}'] = dict(file_sha_matches=(fsha == meta['low_basis_file']['sha256']),
                                  compare_producer_vs_own=V.compare(P, own_low, p))
        print(k, res[f'D5low_k{k}'], flush=True)
    args.out.write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
