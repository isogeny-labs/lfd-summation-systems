"""Freeze the three-summand pilot inputs before any timing."""
from __future__ import annotations

import hashlib
import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from csd.gf2_ghs import GF2n  # noqa: E402
from csd.ic_dreg import BinaryCurve  # noqa: E402

HERE = Path(__file__).resolve().parent
PARENT = HERE / 'manifest.json'
PROTOCOL = HERE / 'M3_PILOT_PROTOCOL.md'
OUT = HERE / 'm3_manifest.json'
CASES = ((23, 11), (29, 13), (29, 14))


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def make():
    parent = json.loads(PARENT.read_text())
    cells = []
    for n, k in CASES:
        item = next(x for x in parent['cells'] if (x['n'], x['k']) == (n, k))
        E = BinaryCurve(GF2n(n, item['modulus']), 0, item['b'])
        fb = []
        for x in range(1 << k):
            p = E.lift_x(x)
            if p is not None:
                fb.append(p)
                if E.neg(p) != p:
                    fb.append(E.neg(p))
        fb.sort()
        assert len(fb) == item['factor_base_points']
        assert digest(json.dumps(fb, separators=(',', ':')).encode()) == item['factor_base_sha256']
        targets = [dict(x) for x in item['targets'][:2]]
        rng = random.Random(20263125 + 100*n + k)
        witness = sorted(rng.choice(fb) for _ in range(3))
        t = E.add(E.add(witness[0], witness[1]), witness[2])
        assert t is not None, 'frozen pilot target must be affine'
        targets.append({'kind': 'planted', 'target': list(t),
                        'witness': [list(p) for p in witness]})
        cells.append({key: item[key] for key in
                      ('n', 'k', 'modulus', 'b', 'order', 'q', 'cofactor', 'G', 'Q',
                       'factor_base_points', 'factor_base_sha256')} | {'targets': targets})
    return {'parent_manifest_sha256': digest(PARENT.read_bytes()),
            'protocol_sha256': digest(PROTOCOL.read_bytes()),
            'source_sha256': digest(Path(__file__).read_bytes()), 'cells': cells}


def main():
    data = make()
    if OUT.exists():
        assert json.loads(OUT.read_text()) == data, 'frozen three-summand manifest changed'
    else:
        OUT.write_text(json.dumps(data, sort_keys=True, indent=2) + '\n')
    for c in data['cells']:
        print(c['n'], c['k'], c['factor_base_points'], len(c['targets']))


if __name__ == '__main__':
    main()
