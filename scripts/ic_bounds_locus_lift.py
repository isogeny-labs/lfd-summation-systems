"""Lift the exceptional cubic C1 (V_4 three below generic) from F_p to Q (bounds item 2).

For each prime, the C1 points of the normalized plane (t = 1; V_3 = 2k+1 and
dim V_4 = 602 at k = 6, from scripts/ic_bounds_locus_plane.py output) determine
the unique cubic through them over F_p.  The cubics are normalized, combined by
CRT and rationally reconstructed; the candidate over Q is then reduced modulo
each prime and compared with every fitted cubic, and its zero set is compared
with the whole sampled data (no generic point may lie on it).

Usage: python3 scripts/ic_bounds_locus_lift.py --out FILE.json
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

from ic_bounds_locus_curve import fit  # noqa: E402

FILES = {
    193: ['experiments/ic-bounds-locus/plane/p193-k6-plane-c*.jsonl', 'experiments/ic-bounds-locus/plane/p193-k6-B0.jsonl'],
    257: ['experiments/ic-bounds-locus/plane/p257-k6-lines-A1-20.jsonl', 'experiments/ic-bounds-locus/plane/p257-k6-B0.jsonl'],
    449: ['experiments/ic-bounds-locus/plane/p449-k6-lines-A1-6.jsonl', 'experiments/ic-bounds-locus/plane/p449-k6-B0.jsonl'],
    577: ['experiments/ic-bounds-locus/plane/p577-k6-lines-A1-6.jsonl'],
}
MONS = [(i, j) for s in range(4) for i in range(s + 1) for j in [s - i]]


def rows(p):
    out = []
    for pattern in FILES[p]:
        for f in sorted(glob.glob(str(ROOT / pattern))):
            out += [json.loads(line) for line in open(f)]
    return out


def ratrec(x, m):
    bound = math.isqrt(m // 2)
    r0, r1, s0, s1 = m, x % m, 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    return Fraction(r1, s1)


def crt(residues, moduli):
    x, m = 0, 1
    for r, n in zip(residues, moduli):
        t = ((r - x) * pow(m, -1, n)) % n
        x, m = x + m * t, m * n
    return x % m, m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    cubics = {}
    data = {}
    for p in FILES:
        R = rows(p)
        if not R:
            continue
        data[p] = R
        pts = sorted({(r['A'], r['B']) for r in R if r['V3'] == 13 and r['V4'] == 602})
        if len(pts) < 10:
            continue  # too few points to fix a cubic: this prime is a test prime only
        d, mons, null = fit(pts, p, 3)
        assert d == 3 and len(null) == 1, (p, d, None if null is None else len(null))
        vec = {m: int(c) for m, c in zip(mons, null[0])}
        inv = pow(vec[(3, 0)], -1, p)
        cubics[p] = dict(points=len(pts), coefficients={f'{i},{j}': vec[(i, j)] * inv % p for (i, j) in MONS})
    primes = sorted(cubics)
    lifted = {}
    for key in cubics[primes[0]]['coefficients']:
        x, m = crt([cubics[p]['coefficients'][key] for p in primes], primes)
        lifted[key] = ratrec(x, m)
    report = dict(primes=primes, cubics=cubics, modulus=math.prod(primes),
                  rational_reconstruction={k: (str(v) if v is not None else None) for k, v in lifted.items()})
    if all(v is not None for v in lifted.values()):
        den = math.lcm(*[v.denominator for v in lifted.values()])
        integer = {k: int(v * den) for k, v in lifted.items()}
        g = math.gcd(*integer.values())
        integer = {k: v // g for k, v in integer.items()}
        report['integer_cubic'] = integer
        # leave-one-prime-out consistency: reconstruct from the other primes and compare
        checks = {}
        for p in sorted(data):
            if p not in primes:
                lc = integer['3,0']
                on = {}
                for r in data[p]:
                    v = sum(c * pow(r['A'], int(k.split(',')[0]), p) * pow(r['B'], int(k.split(',')[1]), p)
                            for k, c in integer.items()) % p
                    if v == 0:
                        key = f"V3={r['V3']},V4={r['V4']}"
                        on[key] = on.get(key, 0) + 1
                c1 = sum(1 for r in data[p] if r['V3'] == 13 and r['V4'] == 602)
                checks[p] = dict(test_prime=True, c1_points_sampled=c1, sampled_points_on_cubic_by_class=on)
                continue
            others = [q for q in primes if q != p]
            rec = {}
            for key in integer:
                x, m = crt([cubics[q]['coefficients'][key] for q in others], others)
                rec[key] = ratrec(x, m)
            checks[p] = dict(reconstructed_without_p=all(rec[k] == lifted[k] for k in integer))
            lc = integer['3,0']
            checks[p]['reduces_to_fitted'] = all(
                (integer[k] * pow(lc, -1, p) - cubics[p]['coefficients'][k]) % p == 0 for k in integer)
            # zero set of the integer cubic on all sampled points of this prime, by class
            on = {}
            for r in data[p]:
                v = sum(c * pow(r['A'], int(k.split(',')[0]), p) * pow(r['B'], int(k.split(',')[1]), p)
                        for k, c in integer.items()) % p
                if v == 0:
                    key = f"V3={r['V3']},V4={r['V4']}"
                    on[key] = on.get(key, 0) + 1
            checks[p]['sampled_points_on_cubic_by_class'] = on
        report['checks'] = checks
    Path(args.out).write_text(json.dumps(report, indent=1) + '\n')
    print(json.dumps(report, indent=1)[:4000])


if __name__ == '__main__':
    main()
