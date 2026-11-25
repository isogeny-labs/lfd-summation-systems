"""Homogeneous FFD receipts for the frozen GF(2^6)..GF(2^11) S3 panel."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE / 'inputs'
REPO = HERE.parents[1]
INPUTS = [ROOT / name for name in (
    'ffd-semaev-n6-sd-input-20261123.json',
    'ffd-semaev-fixed-basis-n7n8-np6-20261123.json',
    'ffd-semaev-fixed-basis-n9n10-np6-20261123.json',
    'ffd-semaev-n11-20261124/fixed-basis-n11-np6.json',
)]
sys.path.insert(0, str(REPO))
from csd.ic_binary_ffd import truncated_ffd  # noqa: E402


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest_json(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def main():
    sources = {name: sha(REPO/name) for name in
               ('csd/ic_binary_ffd.py', 'csd/ic_boolean_closure.py', 'csd/elim.py')}
    files = {str(p.relative_to(ROOT)): sha(p) for p in INPUTS}
    systems = {}
    for path in INPUTS:
        data = json.loads(path.read_text())
        assert data['complete'] is True and data['nprime'] == 6
        for row in data['rows']:
            key = row['n'], row['target_x']
            assert key not in systems
            systems[key] = row
    assert len(systems) == 1977
    outpath = HERE/'semaev-homogeneous-ffd.jsonl'
    assert not outpath.exists()
    with outpath.open('w') as out:
        out.write(json.dumps({'kind': 'header', 'inputs_sha256': files,
                              'sources_sha256': sources, 'cases': len(systems),
                              'convention': 'F2 square-zero Koszul/Frobenius'}, sort_keys=True) + '\n')
        for (n, target), row in sorted(systems.items()):
            result = truncated_ffd(row['equations'], 12,
                                   max_columns=16384, max_bytes=128_000_000)
            out.write(json.dumps({'kind': 'case', 'n': n, 'target_x': target,
                                  'equations_sha256': digest_json(row['equations']),
                                  'result': result}, sort_keys=True) + '\n')
    print(outpath, 'cases', len(systems), 'sha256', sha(outpath))


if __name__ == '__main__':
    main()
