"""Independent bit-vector replay of the frozen binary homogeneous FFD ranks."""
from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = Path(__file__).resolve().parents[2]
SWEEP = REPO / 'experiments/ic-curve-sweep'


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rank(rows) -> int:
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            old = pivots.get(pivot)
            if old is None:
                pivots[pivot] = row
                break
            row ^= old
    return len(pivots)


def masks(n: int, d: int):
    if d < 0 or d > n:
        return []
    return [sum(1 << j for j in c) for c in combinations(range(n), d)]


def mul(top, monomial):
    return [t | monomial for t in top if not t & monomial]


def audit_row(system, result):
    eqs = sorted({tuple(sorted(set(eq))) for eq in system['equations'] if eq})
    assert all(0 <= m < (1 << system['variables']) for eq in eqs for m in eq)
    if any(eq == (0,) for eq in eqs):
        assert result['status'] == 'constant-input' and result['exact_ffd'] is None
        return
    assert result['status'] == 'exact'
    degrees = [max(m.bit_count() for m in eq) for eq in eqs]
    tops = [tuple(m for m in eq if m.bit_count() == d) for eq, d in zip(eqs, degrees)]
    assert degrees == result['degrees']
    scanned = result['scanned']
    assert [r['degree'] for r in scanned] == list(range(min(degrees), result['exact_ffd'] + 1))
    for record in scanned:
        D = record['degree']
        cols = [(i, m) for i, d in enumerate(degrees) for m in masks(system['variables'], D - d)]
        where = {col: j for j, col in enumerate(cols)}
        outputs = masks(system['variables'], D)
        output_index = {m: j for j, m in enumerate(outputs)}
        images = [sum(1 << output_index[t] for t in mul(tops[i], m)) for i, m in cols]
        image_rank = rank(images)
        trivial = []
        for i, d in enumerate(degrees):
            for mono in masks(system['variables'], D - 2 * d):
                trivial.append(sum(1 << where[(i, t)] for t in mul(tops[i], mono)))
            for j in range(i + 1, len(degrees)):
                for mono in masks(system['variables'], D - d - degrees[j]):
                    left = sum(1 << where[(i, t)] for t in mul(tops[j], mono))
                    right = sum(1 << where[(j, t)] for t in mul(tops[i], mono))
                    trivial.append(left ^ right)
        for vector in trivial:
            image = 0
            while vector:
                bit = vector & -vector
                image ^= images[bit.bit_length() - 1]
                vector ^= bit
            assert image == 0
        trivial_rank = rank(trivial)
        syzygy = len(cols) - image_rank
        assert record == {'degree': D, 'coefficient_columns': len(cols),
                          'image_rank': image_rank, 'syzygy_dimension': syzygy,
                          'trivial_dimension': trivial_rank,
                          'nontrivial_dimension': syzygy - trivial_rank}
        assert (syzygy > trivial_rank) == (D == result['exact_ffd'])


def main():
    receipt_path = HERE / 'binary-ffd.jsonl'
    lines = [json.loads(line) for line in receipt_path.open()]
    header = lines[0]
    assert header['kind'] == 'header' and header['systems'] == 1745
    for name, value in header['input_sha256'].items():
        assert digest(SWEEP / name) == value
    for name, value in header['source_sha256'].items():
        assert digest(REPO / name) == value
    systems = {s['system_hash']: s for s in map(json.loads, (SWEEP / 'systems.jsonl').open())
               if s['kind'] == 'binary'}
    results = {r['system_hash']: r for r in map(json.loads, (SWEEP / 'results.jsonl').open())}
    seen = set()
    for row in lines[1:]:
        key = row['system_hash']
        assert row['kind'] == 'system' and key in systems and key not in seen
        assert row['variables'] == systems[key]['variables']
        audit_row(systems[key], row['result'])
        seen.add(key)
    assert seen == set(systems)
    cells = defaultdict(set)
    for case in map(json.loads, (SWEEP / 'cases.jsonl').open()):
        if case['kind'] in ('binary', 'koblitz'):
            assert case['system_hash'] in systems
            cells[(case['kind'], case['n'], case['m'])].add(case['system_hash'])
    records = {row['system_hash']: row['result'] for row in lines[1:]}
    table = []
    for key, hashes in sorted(cells.items()):
        ffd = Counter(records[h]['exact_ffd'] for h in hashes)
        lfd = Counter(results[h]['exact_lfd'] for h in hashes)
        assert all(results[h]['lfd_lower_bound'] == results[h]['lfd_upper_bound'] ==
                   results[h]['exact_lfd'] for h in hashes)
        table.append({'field_kind': key[0], 'n': key[1], 'm': key[2],
                      'systems': len(hashes),
                      'ffd_exact_counts': {str(k): v for k, v in sorted(ffd.items(), key=lambda x: str(x[0]))},
                      'lfd_exact_counts': {str(k): v for k, v in sorted(lfd.items())}})
    output = {'status': 'passed', 'scope': 'frozen binary and Koblitz m=2 sweep systems',
              'systems': len(seen), 'cells': table,
              'input_sha256': header['input_sha256'],
              'ffd_receipt_sha256': digest(receipt_path),
              'audit_source_sha256': digest(Path(__file__)),
              'rank_method': 'independent integer-bit Gaussian elimination',
              'constant_input_convention': 'FFD not defined; recorded separately'}
    dest = HERE / 'audit-binary-ffd.json'
    dest.write_text(json.dumps(output, sort_keys=True, indent=2) + '\n')
    print(dest, 'systems', len(seen), 'cells', len(table))


if __name__ == '__main__':
    main()
