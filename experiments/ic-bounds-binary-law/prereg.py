#!/usr/bin/env python3
"""Write the preregistered predictions (predictions.json) BEFORE any closure is run.

For every case of the frozen cells below (case specifications from
scripts/ic_bounds_binary_run.py `specs`, seed 20261110) this script builds the
Boolean equations and the independently enumerated zero set Z (field arithmetic,
csd.ic_bounds_binary.m2_system; no closure, no LFD), computes the structural class
and e = r_eff - m_eff (structure.analyse), and records the LFD predicted by the two
frozen candidate laws of analysis.py:

  A (e only):  affine class LFD = 2 if e <= -1, 3 if e <= 3, 4 if e <= 8, 5 if e >= 9
  A_Z (grid rule, measured |Z|, k >= 8): LFD = 5 if |Z| > 470 else 4
  B (k + e):   affine class LFD = 2 + [k+e >= 3] + [k+e >= 9] + [k+e >= 19]

Only affine-class (x_R != 0, 1 not in the input span) cases enter the A/B scoring.
Const cases (1 in the input span) are predicted LFD <= 2 by the trivial argument
(1 is a combination of degree-2 inputs); order-two cases are recorded without a
scored prediction (exploratory).
"""
import json, os, sys, hashlib, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, 'scripts')); sys.path.insert(0, HERE)
from ic_bounds_binary_run import specs  # noqa: E402
from csd import ic_dreg as I  # noqa: E402
from csd import ic_bounds_binary as Bb  # noqa: E402
from csd.ic_research_campaign import field  # noqa: E402
from structure import analyse  # noqa: E402
from analysis import A_threshold, B_threshold  # noqa: E402

SEED = 20261110
CELLS = [
    # name, n, k, random bases, invariant bases, purpose
    ('P1', 10, 9, 0, 1, 'trace-zero hyperplane: e = 9, k+e = 18 (A: 5, B: 4)'),
    ('P2', 19, 10, 1, 0, 'random base: e = 1, k+e = 11 (A: 3, B: 4)'),
    ('P3', 16, 9, 1, 0, 'random base: e = 2, k+e = 11 (A: 3, B: 4)'),
    ('P4', 19, 8, 1, 0, 'random base: e = -3, k+e = 5 (A: 2, B: 3)'),
    ('P5', 11, 9, 1, 0, 'control, random base: e = 7, k+e = 16 (A: 4, B: 4)'),
]


def main():
    rows = []
    for name, n, k, nr, ni, purpose in CELLS:
        for spec in specs('m2', [n], [k], nr, ni, 4, SEED):
            F = field(n)
            E = I.BinaryCurve(F, spec['a'], spec['b'])
            planted = [tuple(P) for P in spec['planted']] if spec['planted'] else None
            system = Bb.m2_system(E, spec['basis'], spec['base_kind'], tuple(spec['target']),
                                  spec['target_class'], planted)
            st = analyse(system['equations'], 2 * k)
            cls = 'const' if st['const'] else ('order-two' if spec['target'][0] == 0 else 'affine')
            e = st['r_eff'] - st['m_eff']
            Z = len(system['points'])
            row = dict(cell=name, case_id=spec['case_id'], system_hash=system['system_hash'], n=n, k=k,
                       a=spec['a'], b=spec['b'], base_kind=spec['base_kind'], target_class=spec['target_class'],
                       cls=cls, r_eff=st['r_eff'], m_eff=st['m_eff'], e=e, Z=Z)
            if cls == 'affine':
                row.update(pred_A=A_threshold(e), pred_B=B_threshold(k + e),
                           pred_A_Z=(5 if Z > 470 else 4) if k >= 8 else None)
            elif cls == 'const':
                row.update(pred_const='<=2')
            rows.append(row)
            print(json.dumps(row))
    src = {f: hashlib.sha256(open(os.path.join(HERE, f), 'rb').read()).hexdigest()
           for f in ('analysis.py', 'structure.py', 'prereg.py')}
    out = dict(written_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), seed=SEED,
               cells=[dict(name=c[0], n=c[1], k=c[2], random_bases=c[3], invariant_bases=c[4], purpose=c[5])
                      for c in CELLS],
               laws=dict(A='affine: 2 if e<=-1, 3 if e<=3, 4 if e<=8, 5 if e>=9 (e = r_eff - m_eff)',
                         A_Z='affine, k>=8: 5 if |Z|>470 else 4 (measured |Z|)',
                         B='affine: 2 + [k+e>=3] + [k+e>=9] + [k+e>=19]'),
               scoring='per distinct affine system: hit iff certified exact LFD equals the prediction; '
                       'counts of hits/misses per law per cell; cells P1-P4 discriminate A from B, P5 is a control',
               commands=['python3 scripts/ic_bounds_binary_run.py --family m2 --n %d --k %d --random-bases %d '
                         '--invariant-bases %d --seed %d --max-degree 6 --max-columns 14000 --max-bytes 700000000 '
                         '--out experiments/ic-bounds-binary-law/results-prereg.jsonl' % (c[1], c[2], c[3], c[4], SEED)
                         for c in CELLS],
               source_sha256=src, cases=rows)
    path = os.path.join(HERE, 'predictions.json')
    assert not os.path.exists(path), 'predictions.json exists: preregistration is write-once'
    json.dump(out, open(path, 'w'), indent=1)
    print('wrote', path, len(rows), 'cases')


if __name__ == '__main__':
    main()
