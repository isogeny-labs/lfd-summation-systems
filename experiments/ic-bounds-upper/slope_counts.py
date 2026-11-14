"""Dimension counts behind the conjectured slope 1/4 (Section 6 of docs/proofs/U1-upper-bounds.md).

(a) Free-extension polynomials.  If V_D(k) is the free extension of a base, its profile satisfies
    P_k(z) = (1+z)^(2k) Q_D(z) mod z^(D+1) with a fixed "virtual depth-0 profile" Q_D.  Q_4 and Q_5
    are read off the certified profiles (experiments/ic-bounds-d5, ic-bounds-fall); the script
    checks that the four certified D = 5 profiles give the same Q_5, tabulates
    f_D(k) = [z^D] (1+z)^(2k) Q_D against |M_<=D(2k)|, and reports the last depth at which
    f_D(k) > |M_<=D| (there the free extension is impossible).

(b) Counting bound of Proposition 5: LFD <= D forces dim I_{D'+1} <= (2k+1) dim I_{D'} for all
    D' >= D.  Exact evaluation ranks on Z-heavy instances (p close to 2^k), and the threshold
    D*(k) = min{D : |M_<=D(2k)| >= 2^k} against floor((k+14)/4).
"""
from __future__ import annotations

import json
import math
import os

import numpy as np

from chain_algebra import Window, Echelon, kernel_basis, solutions

HERE = os.path.dirname(os.path.abspath(__file__))


def M(k, D):
    return sum(math.comb(2 * k, j) for j in range(D + 1))


def Q_of(profile, k):
    D = len(profile) - 1
    inv = [(-1) ** j * math.comb(2 * k + j - 1, j) for j in range(D + 1)]
    return [sum(profile[i] * inv[j - i] for i in range(j + 1)) for j in range(D + 1)]


def f(Q, k, j=None):
    D = len(Q) - 1
    j = D if j is None else j
    return sum(Q[i] * math.comb(2 * k, j - i) for i in range(j + 1))


def part_a():
    prof5 = {10: (0, 0, 2, 138, 3163, 18569), 11: (0, 0, 2, 142, 3441, 25033),
             12: (0, 0, 2, 146, 3727, 32057), 13: (0, 0, 2, 150, 4021, 39657)}
    Q5s = {k: Q_of(pr, k) for k, pr in prof5.items()}
    assert len({tuple(q) for q in Q5s.values()}) == 1
    Q5 = Q5s[10]
    # D = 4: dim 4k^2 + 172k - 571, low_3 = 4k + 87, low_2 = 2 (k = 6..13)
    prof4 = {k: (0, 0, 2, 4 * k + 87, 4 * k * k + 172 * k - 571) for k in range(6, 14)}
    Q4s = {k: Q_of(pr, k) for k, pr in prof4.items()}
    assert len({tuple(q) for q in Q4s.values()}) == 1
    Q4 = Q4s[6]
    Q3 = [0, 0, 1, 1]            # dim V_3 = 2k + 1 = span(g, u g)
    out = dict(Q3=Q3, Q4=Q4, Q5=Q5, table=[])
    for D, Q in ((3, Q3), (4, Q4), (5, Q5)):
        rows = []
        for k in range(2, 17):
            rows.append(dict(k=k, free=f(Q, k), monomials=M(k, D), free_exceeds=f(Q, k) > M(k, D)))
        last = max((r['k'] for r in rows if r['free_exceeds']), default=None)
        out['table'].append(dict(D=D, Q=Q, rows=rows, last_depth_free_impossible=last,
                                 first_depth_free_possible=None if last is None else last + 1,
                                 conjectured_base=4 * D - 10))
    return out


def part_b():
    out = dict(instances=[], thresholds=[])
    for (p, a, b, t, k, dmax) in [(17, 2, 3, 5, 4, 5), (257, 1, 3, 190, 8, 4), (257, 5, 7, 11, 8, 4),
                                  (257, 1, 3, 190, 7, 5), (193, 5, 7, 11, 6, 6), (97, 3, 4, 5, 5, 6)]:
        Z = solutions(p, a, b, t, k)
        rec = dict(p=p, a=a, b=b, t=t, k=k, points=len(Z), degrees=[])
        prev = None
        for D in range(1, dmax + 1):
            w = Window(k, k, D)
            Kb = kernel_basis(w.eval_matrix(Z, p), p)
            dimI = len(Kb)
            row = dict(D=D, monomials=w.n, eval_rank=w.n - dimI, dimI=dimI)
            if prev is not None:
                # generation test I_D = I_{D-1} + sum_v v I_{D-1} and the counting inequality
                pw, pK = prev
                emb = np.zeros((len(pK), w.n), dtype=np.int64)
                emb[:, [w.index[m] for m in pw.mons]] = pK
                E = Echelon(w.n, p)
                E.add(emb)
                for vi in range(len(w.vars)):
                    if E.rank == dimI:
                        break
                    E.add(w.shift(emb, vi))
                row.update(generated_from_below=(E.rank == dimI), span_from_below=E.rank,
                           counting_ok=(dimI <= (2 * k + 1) * len(pK)))
            rec['degrees'].append(row)
            prev = (w, Kb)
        out['instances'].append(rec)
    for k in range(4, 129, 4):
        Dstar = next(D for D in range(1, 2 * k + 1) if M(k, D) >= 2 ** k)
        Dstar2 = next(D for D in range(1, 2 * k + 1) if M(k, D) >= 2 ** (k + 1))
        out['thresholds'].append(dict(k=k, Dstar_Z_2k=Dstar, Dstar_Z_2k1=Dstar2, conjecture=(k + 14) // 4,
                                      selector_bound=k + 2))
    return out


if __name__ == '__main__':
    res = dict(free_extension=part_a(), counting=part_b())
    json.dump(res, open(os.path.join(HERE, 'results-slope-counts.json'), 'w'), indent=1, sort_keys=True)
    for t in res['free_extension']['table']:
        print('D', t['D'], 'Q', t['Q'], 'last depth with f > |M|:', t['last_depth_free_impossible'], 'conjectured base', t['conjectured_base'])
        print('   ', [(r['k'], r['free'], r['monomials']) for r in t['rows'] if 3 <= r['k'] <= 12])
    for inst in res['counting']['instances']:
        print({x: inst[x] for x in ('p', 'k', 'points')})
        for r in inst['degrees']:
            print('   ', r)
    print([(r['k'], r['Dstar_Z_2k'], r['Dstar_Z_2k1'], r['conjecture']) for r in res['counting']['thresholds']])
