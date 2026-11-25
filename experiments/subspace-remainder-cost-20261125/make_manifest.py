"""Freeze exact binary curve, factor-base, and target identities before timing."""
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
FIXTURES = ((23, 8, (1 << 23) | (1 << 5) | 1, 894692, 349589, 24),
            (23, 10, (1 << 23) | (1 << 5) | 1, 894692, 349589, 24),
            (23, 11, (1 << 23) | (1 << 5) | 1, 894692, 349589, 24),
            (29, 12, (1 << 29) | (1 << 2) | 1, 326057084, 8388791, 64),
            (29, 13, (1 << 29) | (1 << 2) | 1, 326057084, 8388791, 64),
            (29, 14, (1 << 29) | (1 << 2) | 1, 326057084, 8388791, 64))


def as_json_point(p):
    return None if p is None else [int(p[0]), int(p[1])]


def generator(E, h, q):
    for x in range(1 << E.F.N):
        p = E.lift_x(x)
        if p is None:
            continue
        g = E.mul(h, p)
        if g is not None:
            assert E.mul(q, g) is None
            return g
    raise AssertionError('missing prime-subgroup generator')


def cell(n, k, modulus, b, q, h, cached_g):
    F = GF2n(n, modulus)
    E = BinaryCurve(F, 0, b)
    g = cached_g.get(n)
    if g is None:
        g = generator(E, h, q)
        cached_g[n] = g
    assert E.on_curve(g) and E.mul(q, g) is None
    Q = E.mul(113, g)
    fb = []
    for x in range(1 << k):
        p = E.lift_x(x)
        if p is not None:
            fb.append(p)
            neg = E.neg(p)
            if neg != p:
                fb.append(neg)
    fb.sort()
    assert len(fb) == len(set(fb))
    assert all(E.on_curve(p) and p[0] < (1 << k) for p in fb)
    targets = []
    rng = random.Random(20261125 + 100*n + k)
    for _ in range(64):
        u, v = rng.randrange(q), rng.randrange(q)
        t = E.add(E.mul(u, g), E.mul(v, Q))
        targets.append({'kind': 'random', 'u': u, 'v': v, 'target': as_json_point(t)})
    rng = random.Random(20262125 + 100*n + k)
    for _ in range(16):
        a, bpoint = rng.choice(fb), rng.choice(fb)
        t = E.add(a, bpoint)
        targets.append({'kind': 'planted', 'target': as_json_point(t),
                        'witness': [as_json_point(p) for p in sorted((a, bpoint))]})
    return {'n': n, 'k': k, 'modulus': modulus, 'b': b, 'order': h*q,
            'q': q, 'cofactor': h, 'G': as_json_point(g), 'Q': as_json_point(Q),
            'factor_base_points': len(fb), 'targets': targets}


def main():
    protocol = (HERE/'PROTOCOL.md').read_bytes()
    source = Path(__file__).read_bytes()
    cached_g = {}
    rows = [cell(*args, cached_g) for args in FIXTURES]
    out = {'protocol_sha256': hashlib.sha256(protocol).hexdigest(),
           'source_sha256': hashlib.sha256(source).hexdigest(), 'cells': rows}
    path = HERE/'manifest.json'
    if path.exists():
        assert json.loads(path.read_text()) == out, 'frozen manifest changed'
    else:
        path.write_text(json.dumps(out, sort_keys=True, indent=2)+'\n')
    for row in rows:
        print(row['n'], row['k'], row['factor_base_points'],
              row['G'], len(row['targets']), flush=True)


if __name__ == '__main__':
    main()
