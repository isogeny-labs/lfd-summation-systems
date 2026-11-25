"""Degree-3 exclusion and degree-4 syzygy witnesses for frozen split cells.

The associated graded chain quotient has u_j^2=0, while auxiliary z_i
remain ordinary variables. FFD here means homogeneous syzygy modulo Koszul
rows in that quotient, matching the convention stated in the paper.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

REPO = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(REPO))
from csd import ic_dreg as I  # noqa: E402
from csd import ic_bounds_split_m as S  # noqa: E402


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def multiply_linear_top(top, mask, p):
    out = {}
    for (term_mask, es), c in top.items():
        if term_mask & mask:
            continue
        key = (term_mask | mask, es)
        out[key] = (out.get(key, 0) + c) % p
    return {k: v for k, v in out.items() if v}


def audit(row):
    p, m, k = row['p'], row['m'], row['k']
    assert p % 2 == 1 and m >= 3 and k >= 2
    E = I.PrimeCurve(p, row['a'], row['b'])
    gens, _ = S.split_generators_m(E, row['target'][0], k, m)
    assert S.system_hash(p, m, k, gens) == row['system_hash']
    degrees = [S.generator_degree(g) for g in gens]
    assert degrees == row['generator_degrees']
    assert degrees == [4] * (m - 2) + [3]
    last = gens[-1]
    top = {key: c for key, c in last.items()
           if key[0].bit_count() + sum(key[1]) == 3}
    assert top  # Degree 3 has just this one generator and a nonzero top.
    assert all(es == (0,) * (m - 3) + (2,) for _, es in top)
    assert all(mask.bit_count() == 1 for mask, _ in top)
    assert len(top) <= 2
    masks = sorted(mask for mask, _ in top)
    if len(masks) == 2:
        a, b = masks
        ca = top[(a, (0,) * (m - 3) + (2,))]
        cb = top[(b, (0,) * (m - 3) + (2,))]
        # (ca*a-cb*b)(ca*a+cb*b)=0 because a^2=b^2=0.
        witness = [(a, ca), (b, -cb % p)]
    else:
        # A sole square-zero chain variable annihilates its own top.
        witness = [(masks[0], 1)]
    image = {}
    for mask, c in witness:
        for key, value in multiply_linear_top(top, mask, p).items():
            image[key] = (image.get(key, 0) + c * value) % p
    assert not any(image.values())
    assert any(c for _, c in witness)
    # No Koszul row can occur in degree 4: the two least input degrees sum to 7.
    assert sorted(degrees)[:2] == [3, 4]
    assert row['status'] == 'exact-unit'
    assert row['exact_lfd'] == row['lfd_lower_bound']
    assert all(rec['verified_closed'] is True for rec in row['profile']
               if not rec['contains_one'])
    return {'p': p, 'm': m, 'k': k, 'seed': row['seed'],
            'system_hash': row['system_hash'], 'ffd_lower_bound': 4,
            'ffd_upper_bound': 4, 'lfd_lower_bound': row['lfd_lower_bound'],
            'lfd_upper_bound': row['exact_lfd'],
            'degree_3_top': [[mask, list(es), c] for (mask, es), c in sorted(top.items())],
            'degree_4_syzygy_coefficient': [[mask, c] for mask, c in witness]}


def main():
    root = REPO / 'experiments/ic-bounds-split-m'
    cells = [(3, 5), (4, 2), (4, 3), (4, 4), (5, 2), (6, 2)]
    files = [root / f'cell-m{m}-k{k}.jsonl' for m, k in cells]
    receipts = []
    table = []
    for (m, k), path in zip(cells, files):
        rows = [json.loads(line) for line in path.read_text().splitlines()]
        assert rows and all((r['m'], r['k']) == (m, k) for r in rows)
        assert len({(r['p'], r['seed']) for r in rows}) == len(rows)
        current = [audit(row) for row in rows]
        receipts += current
        table.append({'m': m, 'k': k, 'systems': len(rows),
                      'fields': dict(sorted(Counter(str(r['p']) for r in rows).items())),
                      'ffd_exact_counts': {'4': len(rows)},
                      'lfd_exact_counts': dict(sorted(Counter(str(r['exact_lfd']) for r in rows).items()))})
    output = {'status': 'passed', 'scope': 'six completed frozen split cells',
              'convention': 'associated graded chain quotient, Koszul trivial module',
              'reason_lower': 'only one degree-3 input, with nonzero top; no degree-3 syzygy',
              'reason_upper': 'explicit nonzero degree-4 annihilator of its two-term top; Koszul degree at least 7',
              'cells': table, 'systems': receipts,
              'input_sha256': {p.name: sha(p) for p in files},
              'source_sha256': {name: sha(REPO/name) for name in
                                ['csd/ic_bounds_split_m.py', 'csd/ic_bounds_m3.py',
                                 'csd/ic_dreg.py']},
              'audit_source_sha256': sha(Path(__file__))}
    dest = HERE / 'audit-split-ffd.json'
    dest.write_text(json.dumps(output, sort_keys=True, indent=2) + '\n')
    print(json.dumps({'status': 'passed', 'cells': len(table), 'systems': len(receipts)}))


if __name__ == '__main__':
    main()
