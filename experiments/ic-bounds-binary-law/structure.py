#!/usr/bin/env python3
"""Effective size of every two-summand system in table.jsonl.

For each distinct system (k < n; the inverse family is skipped) the Boolean
equations are taken from experiments/ic-curve-sweep/systems.jsonl (k <= 7) or
regenerated with csd.ic_dreg.binary_descent from the stored specification
(k = 8..10; the sha256 system hash is compared with the stored one).  Over GF(2):

  rank      rank of the equation span
  m_eff     rank of the degree-2 parts (number of independent quadratic equations)
  L         dimension of the degree <= 1 part of the span (= rank - m_eff)
  lin_rank  rank of the linear parts of that degree <= 1 subspace
  const     1 lies in the span (unit at degree <= 2 by the inputs alone)
  r_eff     r - lin_rank (variables left after the linear equations)
  bilinear  every degree-2 monomial is a_i b_j (one variable from each summand)

Output: structure.jsonl keyed by system_hash.
"""
import json, os, sys, hashlib

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)


def rank_rows(rows):
    """GF(2) rank of integer bit-vectors; returns (rank, basis dict pivot->row)."""
    basis = {}
    for v in rows:
        while v:
            p = v.bit_length() - 1
            if p in basis:
                v ^= basis[p]
            else:
                basis[p] = v
                break
    return len(basis), basis


def analyse(equations, r):
    k = r // 2
    mons = sorted({m for e in equations for m in e})
    idx = {m: i for i, m in enumerate(mons)}
    # order columns so that degree-2 monomials are the high bits
    order = sorted(mons, key=lambda m: (m.bit_count(), m))
    pos = {m: i for i, m in enumerate(order)}
    vecs = []
    for e in equations:
        v = 0
        for m in e:
            v ^= 1 << pos[m]
        vecs.append(v)
    rank, basis = rank_rows(vecs)
    low_mask = sum(1 << pos[m] for m in order if m.bit_count() <= 1)
    # reduced echelon: rows with pivot in a degree<=1 column span the degree<=1 subspace
    low = [v for p, v in basis.items() if (1 << p) & low_mask]
    L = len(low)
    m_eff = rank - L
    lin = []
    const = False
    for v in low:
        w = 0
        for m in order:
            if (v >> pos[m]) & 1 and m.bit_count() == 1:
                w ^= m
        lin.append(w)
    lin_rank, _ = rank_rows(lin)
    # 1 in span: lowest basis row equal to the constant monomial alone
    if 0 in pos:
        c = 1 << pos[0]
        red = c
        for p in sorted(basis, reverse=True):
            if (red >> p) & 1:
                red ^= basis[p]
        const = red == 0
    A = (1 << k) - 1
    bil = all(m.bit_count() < 2 or (m.bit_count() == 2 and (m & A).bit_count() == 1) for m in mons)
    return dict(rank=rank, m_eff=m_eff, L=L, lin_rank=lin_rank, const=const, r_eff=r - lin_rank,
                bilinear=bil, equations=len(equations))


def main():
    table = [json.loads(l) for l in open(os.path.join(HERE, 'table.jsonl'))]
    want = {}
    for t in table:
        if t['source'] != 'ic-binary-inverse':
            want.setdefault(t['system_hash'], t)
    eqs = {}
    for l in open(os.path.join(ROOT, 'experiments/ic-curve-sweep/systems.jsonl')):
        x = json.loads(l)
        if x['system_hash'] in want:
            eqs[x['system_hash']] = x['equations']
    from csd import ic_dreg as I
    from csd.ic_research_campaign import field
    from csd.ic_arbitrary_bases import canonical_basis
    m2 = {}
    for l in open(os.path.join(ROOT, 'experiments/ic-bounds-binary/results-m2.jsonl')):
        x = json.loads(l)
        m2.setdefault(x['system_hash'], x)
    for h, x in m2.items():
        if h in eqs:
            continue
        F = field(x['n'])
        assert F.modulus == x['modulus']
        E = I.BinaryCurve(F, x['a'], x['b'])
        basis = canonical_basis(x['basis'])
        e = sorted(set(tuple(sorted(row)) for row in I.binary_descent(E, [basis, basis], x['target'][0]) if row))
        system = dict(kind='binary', variables=2 * x['k'], equations=e)
        key = hashlib.sha256(json.dumps(system, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        assert key == h, ('hash mismatch', x['case_id'])
        eqs[h] = [list(r) for r in e]
    out = open(os.path.join(HERE, 'structure.jsonl'), 'w')
    for h, t in want.items():
        s = analyse(eqs[h], 2 * t['k'])
        s['system_hash'] = h
        out.write(json.dumps(s, sort_keys=True) + '\n')
    out.close()
    print('systems', len(want))


if __name__ == '__main__':
    main()
