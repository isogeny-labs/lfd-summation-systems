"""Freeze a bounded high-kernel-dimension panel from the existing campaign."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from csd.gf2_ghs import GF2n  # noqa: E402
from csd.ic_dreg import BinaryCurve  # noqa: E402

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / 'subspace-remainder-cost-20261125' / 'manifest.json'
CELLS = {(23, 11), (29, 13), (29, 14)}
INDICES = (0, 1, 64, 65)


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def main() -> None:
    parent_raw = PARENT.read_bytes()
    parent = json.loads(parent_raw)
    cells = []
    for c in parent['cells']:
        if (c['n'], c['k']) not in CELLS:
            continue
        E = BinaryCurve(GF2n(c['n'], c['modulus']), 0, c['b'])
        fb = sorted((x, y) for x in range(1 << c['k'])
                    for p in [E.lift_x(x)] if p is not None
                    for y in ({p[1], E.neg(p)[1]}))
        assert len(fb) == c['factor_base_points']
        assert digest(json.dumps(fb, separators=(',', ':')).encode()) == c['factor_base_sha256']
        targets = [{'parent_index': i, **c['targets'][i]} for i in INDICES]
        for p in fb:
            doubled = E.add(p, p)
            if doubled is not None and doubled[0] != 0:
                targets.append({'parent_index': None, 'kind': 'diagonal',
                                'target': list(doubled), 'witness': [list(p), list(p)]})
                break
        else:
            raise AssertionError('no nondegenerate diagonal target')
        cells.append({key: c[key] for key in ('n', 'k', 'modulus', 'b',
                      'factor_base_points', 'factor_base_sha256')}
                     | {'targets': targets})
    assert {(c['n'], c['k']) for c in cells} == CELLS
    out = {'parent_sha256': digest(parent_raw),
           'protocol_sha256': digest((HERE/'PROTOCOL.md').read_bytes()),
           'generator_sha256': digest(Path(__file__).read_bytes()),
           'cells': cells}
    raw = json.dumps(out, sort_keys=True, indent=2) + '\n'
    path = HERE/'manifest.json'
    if path.exists():
        assert path.read_text() == raw, 'frozen manifest changed'
    else:
        path.write_text(raw)
    print([(c['n'], c['k'], len(c['targets'])) for c in cells])


if __name__ == '__main__':
    main()
