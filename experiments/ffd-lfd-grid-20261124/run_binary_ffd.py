"""Compute source-bound homogeneous FFD for the frozen binary sweep systems.

This uses the repository's explicitly declared square-zero, Koszul/Frobenius
convention. It keeps every system, including capped and constant-input cases.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
REPO = Path(__file__).resolve().parents[2]
SWEEP = REPO / 'experiments/ic-curve-sweep'
OUTPUT = HERE / 'binary-ffd.jsonl'
sys.path.insert(0, str(REPO))
from csd.ic_binary_ffd import truncated_ffd  # noqa: E402


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    sources = {str(p.relative_to(REPO)): sha(p) for p in (
        REPO / 'csd/ic_binary_ffd.py',
        REPO / 'csd/ic_boolean_closure.py',
        REPO / 'csd/elim.py',
    )}
    inputs = {name: sha(SWEEP / name) for name in
              ('cases.jsonl', 'systems.jsonl', 'results.jsonl')}
    systems = [json.loads(line) for line in (SWEEP / 'systems.jsonl').open()]
    expected = {s['system_hash']: s for s in systems if s['kind'] == 'binary'}
    assert len(expected) == sum(s['kind'] == 'binary' for s in systems)
    seen = set()
    if OUTPUT.exists():
        lines = [json.loads(line) for line in OUTPUT.open()]
        assert lines[0] == {'kind': 'header', 'input_sha256': inputs,
                            'source_sha256': sources, 'systems': len(expected),
                            'convention': 'F2 square-zero Koszul/Frobenius'}
        for record in lines[1:]:
            assert record['kind'] == 'system' and record['system_hash'] in expected
            assert record['system_hash'] not in seen
            seen.add(record['system_hash'])
    else:
        HERE.mkdir(parents=True, exist_ok=True)
        with OUTPUT.open('w') as out:
            out.write(json.dumps({'kind': 'header', 'input_sha256': inputs,
                                  'source_sha256': sources, 'systems': len(expected),
                                  'convention': 'F2 square-zero Koszul/Frobenius'},
                                 sort_keys=True) + '\n')
            out.flush()
    counts = Counter()
    with OUTPUT.open('a') as out:
        for s in systems:
            if s['kind'] != 'binary' or s['system_hash'] in seen:
                continue
            result = truncated_ffd(s['equations'], s['variables'],
                                   max_columns=16384, max_bytes=128_000_000)
            record = {'kind': 'system', 'system_hash': s['system_hash'],
                      'variables': s['variables'], 'result': result}
            out.write(json.dumps(record, sort_keys=True) + '\n')
            out.flush()
            seen.add(s['system_hash'])
            counts[result['status']] += 1
            if len(seen) % 100 == 0:
                print(f'{len(seen)}/{len(expected)} {dict(counts)}', flush=True)
    assert seen == set(expected)
    print(f'complete {len(seen)}; newly computed {dict(counts)}', flush=True)


if __name__ == '__main__':
    main()
