"""Fit the exceptional curves of the normalized (A, B) plane (bounds item 2).

Reads the plane sweep (scripts/ic_bounds_locus_plane.py, t = 1, every (A, B) in
F_p^2), groups the points by closure signature, and for each non-generic class
finds the polynomials P(A, B) over F_p of least total degree that vanish on all
of its points (null space of the monomial evaluation matrix), reports the
dimension of that space and, when it is one, the factorization of P over F_p
and whether P vanishes ONLY on the class (no generic point on P = 0).

Usage
  python3 scripts/ic_bounds_locus_curve.py --plane 'experiments/ic-bounds-locus/plane/p193-k6-plane-c*.jsonl' --out FILE.json
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from csd import elim  # noqa: E402


def signature(r):
    if r['saturated4'] and not r['one4']:
        return ('V3', r['V3'], 'low3', r['low3'], 'V4', 'saturated (|M_<=4| - |Z|)')
    return ('V3', r['V3'], 'low3', r['low3'], 'V4', r['V4'], 'low4', r['low4'], 'one4', r['one4'])


def monomials(d):
    return [(i, j) for s in range(d + 1) for i in range(s + 1) for j in [s - i]]


def nullspace_mod_p(M, p):
    """Basis of {x : M x = 0} over F_p (columns of M)."""
    rows, cols = M.shape
    A = M.copy() % p
    pivots = []
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if A[i, c]:
                piv = i
                break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        A[r] = A[r] * pow(int(A[r, c]), -1, p) % p
        for i in range(rows):
            if i != r and A[i, c]:
                A[i] = (A[i] - A[i, c] * A[r]) % p
        pivots.append(c)
        r += 1
        if r == rows:
            break
    free = [c for c in range(cols) if c not in pivots]
    basis = []
    for f in free:
        v = np.zeros(cols, dtype=np.int64)
        v[f] = 1
        for i, c in enumerate(pivots):
            v[c] = (-A[i, f]) % p
        basis.append(v)
    return basis


def fit(points, p, max_degree):
    for d in range(1, max_degree + 1):
        mons = monomials(d)
        if len(points) < len(mons) - 1 and d > 1:
            pass
        M = np.array([[pow(A, i, p) * pow(B, j, p) % p for (i, j) in mons] for A, B in points], dtype=np.int64)
        rank = elim.fp_echelon(M, p)[0]
        if rank < len(mons):
            null = nullspace_mod_p(M, p)
            return d, mons, null
    return None, None, None


def poly_str(vec, mons):
    terms = []
    for c, (i, j) in zip(vec, mons):
        if c:
            terms.append(f'{int(c)}*A^{i}*B^{j}')
    return ' + '.join(terms)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--plane', required=True)
    ap.add_argument('--max-degree', type=int, default=24)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    rows = [json.loads(line) for f in sorted(glob.glob(args.plane)) for line in open(f)]
    p = rows[0]['prime']
    k = rows[0]['k']
    seen = {(r['A'], r['B']) for r in rows}
    classes = defaultdict(list)
    for r in rows:
        classes[signature(r)].append((r['A'], r['B']))
    generic = max(classes, key=lambda s: len(classes[s]))
    report = dict(prime=p, k=k, points=len(rows), complete_plane=len(seen) == p * p,
                  generic_signature=list(generic), classes=[])
    for sig, pts in sorted(classes.items(), key=lambda kv: -len(kv[1])):
        item = dict(signature=list(sig), count=len(pts), generic=(sig == generic))
        if sig != generic:
            item['points'] = sorted(pts) if len(pts) <= 400 else None
            d, mons, null = fit(sorted(set(pts)), p, args.max_degree)
            item['least_degree'] = d
            item['nullspace_dimension'] = (len(null) if null is not None else None)
            if null is not None and len(null) == 1:
                vec = null[0]
                item['polynomial'] = poly_str(vec, mons)
                # every plane point on P = 0, by class
                on = Counter()
                for r in rows:
                    val = sum(int(c) * pow(r['A'], i, p) * pow(r['B'], j, p) for c, (i, j) in zip(vec, mons)) % p
                    if val == 0:
                        on[str(signature(r))] += 1
                item['plane_points_on_P_by_class'] = dict(on)
                try:
                    import sympy as sp
                    A, B = sp.symbols('A B')
                    expr = sum(int(c) * A ** i * B ** j for c, (i, j) in zip(vec, mons))
                    item['factorization_mod_p'] = str(sp.factor_list(sp.Poly(expr, A, B, modulus=p)))
                except Exception as error:  # noqa: BLE001
                    item['factorization_mod_p'] = f'unavailable: {error}'
        report['classes'].append(item)
        print(item['count'], sig, item.get('least_degree'), item.get('nullspace_dimension'),
              item.get('factorization_mod_p', '')[:300])
    Path(args.out).write_text(json.dumps(report, indent=1, default=str) + '\n')


if __name__ == '__main__':
    main()
