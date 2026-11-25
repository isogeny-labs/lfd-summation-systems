"""Independent receipt check against frozen points, group law and arm parity."""
from __future__ import annotations

import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from csd.gf2_ghs import GF2n  # noqa: E402
from csd.ic_dreg import BinaryCurve  # noqa: E402

HERE = Path(__file__).resolve().parent


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    manifest_raw = (HERE/'manifest.json').read_bytes()
    result = json.loads((HERE/'results.json').read_text())
    assert result['manifest_sha256'] == sha(manifest_raw)
    assert result['runner_sha256'] == sha((HERE/'run.py').read_bytes())
    manifest = json.loads(manifest_raw)
    assert manifest['parent_sha256'] == sha((HERE.parent/'subspace-remainder-cost-20261125'/'manifest.json').read_bytes())
    assert manifest['protocol_sha256'] == sha((HERE/'PROTOCOL.md').read_bytes())
    assert manifest['generator_sha256'] == sha((HERE/'make_manifest.py').read_bytes())
    for label, filename in [('ic_rootfind_sha256','ic_rootfind.py'),
                            ('rootfind_c_sha256','_rootfind.c'),
                            ('ic_dreg_sha256','ic_dreg.py'),
                            ('gf2_ghs_sha256','gf2_ghs.py')]:
        assert result['environment'][label] == sha((ROOT/'csd'/filename).read_bytes())
    cases = defaultdict(dict)
    counts = defaultdict(int)
    for row in result['receipts']:
        n, k = row['cell']
        c = next(c for c in manifest['cells'] if (c['n'],c['k']) == (n,k))
        target_row = c['targets'][row['index']]
        E = BinaryCurve(GF2n(n,c['modulus']),0,c['b'])
        target = tuple(target_row['target'])
        assert row['kind'] == target_row['kind']
        assert row['status'] == 'complete', (n,k,row['index'],row['arm'],row['status'])
        assert row['charged_wall_seconds'] >= 0 and row['charged_cpu_seconds'] >= 0
        assert row['peak_rss_raw'] > 0
        pairs = []
        for pair in row['pairs']:
            a, b = tuple(pair[0]), tuple(pair[1])
            assert a <= b and a[0] < (1<<k) and b[0] < (1<<k)
            assert E.on_curve(a) and E.on_curve(b) and E.add(a,b) == target
            pairs.append((a,b))
        assert pairs == sorted(set(pairs))
        if 'witness' in target_row:
            assert tuple(sorted(map(tuple,target_row['witness']))) in pairs
        key = (n,k,row['index'])
        assert row['arm'] not in cases[key]
        cases[key][row['arm']] = pairs
        counts[row['arm']] += 1
    assert len(cases) == 15
    for key, arms in cases.items():
        assert set(arms) == {'sat','lrf','image','brute'}, key
        assert all(arms[a] == arms['brute'] for a in arms), key
    assert all(v == 15 for v in counts.values())
    print(json.dumps({'cases':len(cases), 'receipts':len(result['receipts']),
                      'arm_counts':dict(counts), 'parity_mismatches':0},sort_keys=True))


if __name__ == '__main__':
    main()
