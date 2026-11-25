"""Independent rank replay and LFD binding for the fixed-basis S3 panel."""
from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path

from audit_binary_ffd import audit_row
from run_semaev_ffd import INPUTS, REPO, ROOT, HERE, digest_json, sha

LFD_RECEIPTS = [ROOT/name for name in (
    'ffd-semaev-sd-n6-20261123.jsonl',
    'ffd-semaev-sd-n7n8-20261123.jsonl',
    'ffd-semaev-sd-panel-20261123.jsonl',
    'ffd-semaev-n11-20261124/exact-sd-n11.jsonl',
)]


def main():
    path = HERE/'semaev-homogeneous-ffd.jsonl'
    receipt = [json.loads(line) for line in path.open()]
    header = receipt[0]
    assert header['kind'] == 'header' and header['cases'] == 1977
    assert header['convention'] == 'F2 square-zero Koszul/Frobenius'
    for rel, value in header['inputs_sha256'].items():
        assert sha(ROOT/rel) == value
    for rel, value in header['sources_sha256'].items():
        assert sha(REPO/rel) == value
    systems = {}
    affine = {}
    for input_path in INPUTS:
        for row in json.loads(input_path.read_text())['rows']:
            key = row['n'], row['target_x']
            assert key not in systems
            systems[key] = row
            affine[key] = row.get('static_first_fall_through_3', row.get('ffd_through_4'))
    lfd = {}
    for input_path in LFD_RECEIPTS:
        lines = [json.loads(line) for line in input_path.open()]
        assert lines[0]['kind'] == 'header'
        for row in lines[1:]:
            key = row['n'], row['target_x']
            assert key in systems and key not in lfd and row['kind'] == 'case'
            assert row['equations_sha256'] == digest_json(systems[key]['equations'])
            assert row['sd_lower'] == row['sd_upper'] == row['lfd']
            lfd[key] = row
    assert set(lfd) == set(systems)
    counts = {n: Counter() for n in range(6, 12)}
    seen = set()
    for row in receipt[1:]:
        key = row['n'], row['target_x']
        assert row['kind'] == 'case' and key in systems and key not in seen
        assert row['equations_sha256'] == digest_json(systems[key]['equations'])
        audit_row({'equations': systems[key]['equations'], 'variables': 12}, row['result'])
        assert row['result']['status'] == 'exact'
        ffd = row['result']['exact_ffd']
        counts[row['n']][(ffd, lfd[key]['lfd'])] += 1
        seen.add(key)
    assert seen == set(systems)
    mismatches = sum(ffd != affine[key] for key, ffd in
                     ((key, next(r['result']['exact_ffd'] for r in receipt[1:]
                                 if (r['n'], r['target_x']) == key)) for key in systems))
    output = {
        'status': 'passed', 'cases': len(seen),
        'scope': 'frozen fixed six-coordinate S3 presentation over GF(2^6)..GF(2^11)',
        'ffd_convention': 'homogeneous square-zero Koszul/Frobenius',
        'ffd_lfd_counts_by_n': {str(n): {f'{f},{l}': v for (f, l), v in sorted(c.items())}
                                for n, c in counts.items()},
        'homogeneous_vs_affine_first_fall_mismatches': mismatches,
        'inputs_sha256': header['inputs_sha256'],
        'lfd_receipts_sha256': {str(p.relative_to(ROOT)): sha(p) for p in LFD_RECEIPTS},
        'ffd_receipt_sha256': sha(path),
        'audit_source_sha256': sha(Path(__file__)),
        'rank_method': 'independent integer-bit Gaussian elimination',
    }
    dest = HERE/'audit-semaev-ffd.json'
    dest.write_text(json.dumps(output, sort_keys=True, indent=2) + '\n')
    print(dest, 'cases', len(seen), 'mismatches', mismatches)


if __name__ == '__main__':
    main()
