"""Item 4 (c): plan the measurement of the depth k_7 at which the two-summand prime LFD reaches 7.

What must be measured.  For a two-summand subgroup instance at depth k with ev(M_<6) onto Z:
  * LFD >= 7 at k  <=>  V_6 is deficient: dim V_6 < |M_<=6| - |Z|  (unit ideal: 1 not in V_6),
    certified only by a verified-closed V_6 (status 'closed-deficient' of csd/ic_bounds_exact10.py);
  * LFD <= 6 at k  <=>  V_6 = I cap M_<=6 (status 'reached'); with V_5 verified deficient at the same
    k the exact LFD is then 6 (strict fall at 6).
So k_7 = least k whose degree-6 closure is deficient; the quantity to measure at each candidate k
is dim V_6 (exact) or, equivalently by the argument in README.md, dim(V_6 cap M_<=5) against
|M_<=5| - |Z|.

This script writes experiments/ic-bounds-exact10/plan.json: per depth the prime (2^k | p-1, the
smallest such below 2^16 when one exists, else 65537), column counts, the memory model of
csd.ic_bounds_exact10.memory_model for the dense split engine at full rank, a time extrapolation
calibrated on the measured mid-scale runs (cols^3 law, single thread, loaded host; an estimate,
not a measurement) and the exact commands.  It runs no closure.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from csd import ic_bounds_exact10 as X
from csd.field import is_prime


def prime_for(k: int) -> int:
    small = [c * (1 << k) + 1 for c in range(1, 64) if c * (1 << k) + 1 < 65536 and is_prime(c * (1 << k) + 1)]
    if small:
        return 12289 if 12289 in small else small[0]
    return next(c * (1 << k) + 1 for c in range(1, 10 ** 6) if is_prime(c * (1 << k) + 1))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--ks', type=int, nargs='+', default=list(range(9, 17)))
    ap.add_argument('--calibration', type=Path, default=ROOT / 'experiments/ic-bounds-exact10/calibration.json',
                    help='measured (cols, seconds) of reached split closures; optional')
    ap.add_argument('--out', type=Path, default=ROOT / 'experiments/ic-bounds-exact10/plan.json')
    args = ap.parse_args()
    alpha = None
    if args.calibration.exists():
        cal = json.loads(args.calibration.read_text())
        pts = [(c['cols'], c['seconds']) for c in cal['runs'] if c['seconds'] > 5]
        if pts:  # least squares of seconds = alpha * cols^3 through the origin
            alpha = sum(s * c ** 3 for c, s in pts) / sum(c ** 6 for c, _ in pts)
    rows = []
    for k in args.ks:
        p = prime_for(k)
        row = dict(k=k, prime=p, word_bytes=2 if p < 65536 else 4)
        for D in (5, 6, 7):
            mem = X.memory_model(k, D, prime=p, combos=2000, tile=1024)
            row[f'D{D}'] = dict(cols=mem['cols'], low_cols=mem['low_cols'], sym_cols=mem['sym_cols'],
                                anti_cols=mem['anti_cols'], parts_gb=round(mem['parts_bytes'] / 1e9, 2),
                                peak_gb_model=round(mem['peak_bytes'] / 1e9, 2),
                                hours_est=(round(alpha * mem['cols'] ** 3 / 3600, 1) if alpha else None))
        tag = f'k{k}-p{p}'
        row['commands'] = {
            'D5_closure_lower_bound': (f'python3 scripts/ic_bounds_exact10_degree6.py --k {k} --degree 5 --prime {p} '
                                       f'--a 1 --b 0 --seed 10 --tag {tag}-unit --ckpt $CKPT/{tag}-D5'),
            'D6_closure': (f'python3 scripts/ic_bounds_exact10_degree6.py --k {k} --degree 6 --prime {p} --a 1 --b 0 '
                           f'--seed 10 --tag {tag}-unit --ckpt $CKPT/{tag}-D6'),
            'D6_proper': (f'python3 scripts/ic_bounds_exact10_degree6.py --k {k} --degree 6 --prime {p} --a 1 --b 0 '
                          f'--seed 11 --want-solutions --tag {tag}-proper --ckpt $CKPT/{tag}-D6-proper'),
        }
        rows.append(row)
        print(json.dumps(dict(k=k, p=p, D6_cols=row['D6']['cols'], D6_peak_gb=row['D6']['peak_gb_model'],
                              D6_hours=row['D6']['hours_est'], D7_cols=row['D7']['cols'],
                              D7_peak_gb=row['D7']['peak_gb_model'])))
    known = dict(
        exact_lfd_by_depth={'k<=5': 4, '6<=k<=9': 5, 'k=10': '>= 6 (three instances, verified-closed V_5)'},
        onset_depths={'LFD 5 (V_4 deficient)': 6, 'LFD 6 (V_5 deficient)': 10},
        dims={'dim V_4 for k>=6 (exact dense, fit)': '4k^2+172k-571', 'dim V_5 at k=10': 18569,
              'dim(V_5 cap M_<=4) at k=10': 3163},
        pkm_magma_dreg={'k=10': 6, 'k=11': 6},
        conjecture=('k_7 lies in 12..17; the onset gaps so far are 6 -> 10; an equal gap gives 14. '
                    'This is a conjecture from two data points, not a bound.'))
    args.out.write_text(json.dumps(dict(known=known, time_alpha_seconds_per_col3=alpha, rows=rows), indent=1) + '\n')


if __name__ == '__main__':
    main()
