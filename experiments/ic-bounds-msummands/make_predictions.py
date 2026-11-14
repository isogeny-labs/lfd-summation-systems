#!/usr/bin/env python3
"""Writes predictions-2026-11-12.json from survey-prime.json and survey-binary.json (generator degrees, |Z|: produced by
`scripts/ic_bounds_msum_run.py survey`, which builds generators and zero sets and runs NO closure).  Write-once."""
import json, os, sys, hashlib
HERE = os.path.dirname(os.path.abspath(__file__))
OUTF = os.path.join(HERE, 'predictions-2026-11-12.json')
if os.path.exists(OUTF):
    sys.exit('predictions file exists: write-once')
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
prime = json.load(open(os.path.join(HERE, 'survey-prime.json')))
binary = json.load(open(os.path.join(HERE, 'survey-binary.json')))
cases = {}
for r in prime:
    mk = r['m'] * r['k']; d = r['deg_g']
    cases[r['name']] = dict(family='prime', m=r['m'], k=r['k'], generator_degree=d, solutions=r['Z'],
        prop1_variables=r['prop1_vars'], kernel_upper_bound=mk + 1,
        pred_a_upper=mk + 1, pred_a_attained=(mk + 1 if d == mk else None), pred_b_d_plus_1=d + 1,
        lfd_lower_bound_if_prop1=(d + 1 if r['prop1_vars'] > 0 else None),
        predicted_lfd=d + 1, predicted_strict_falls=[d + 1])
for r in binary:
    rr = r['r']; d = r['d_in']
    cases[r['name']] = dict(family='binary', m=r['m'], k=r['k'], n=r['n'], variables=rr, input_degree=d, solutions=r['Z'],
        kernel_upper_bound=rr + 1, pred_c_upper=min(rr, d + 2), pred_c_point=min(rr, d + 1),
        point_scored=(d == rr), gpu_cell=(r['name'].startswith('B-m4n9k4')))
doc = dict(
 frozen='2026-11-12', file='predictions-2026-11-12.json', write_once=True,
 frozen_before='every closure for m >= 4 and every m = 3 case not already recorded (the cases below); generators and Z (surveys) were built before, '
               'closures for these cases were not run',
 disclosure_smoke_tests=('Before this freeze two closures were run as smoke tests of the pipeline: prime m = 4, k = 2 and m = 5, k = 2 on '
   'p = 97, instance seed 1 (planted), engine dense, results LFD 9 and 11 with one strict fall at the same degree. Those two instances '
   '(seed 1) are NOT among the cases below (cases use seeds 10..15); the predictions for m = 4, 5 at k = 2 are therefore informed by these '
   'two smoke results and are not blind. Everything else (m = 4 at k = 3, binary cells, GPU cells) is blind.'),
 inputs={'survey-prime.json': sha(os.path.join(HERE, 'survey-prime.json')), 'survey-binary.json': sha(os.path.join(HERE, 'survey-binary.json')),
         'csd/ic_bounds_msum.py': sha(os.path.join(HERE, '../../csd/ic_bounds_msum.py')),
         'scripts/ic_bounds_msum_run.py': sha(os.path.join(HERE, '../../scripts/ic_bounds_msum_run.py'))},
 candidates={
  'a': dict(statement='kernel-checked upper bound LFD <= m k + 1 holds for every prime case, and is attained when the generator has the maximal chain degree d = m k '
            '(as m = 3, k = 2: deg g = 6 = m k, LFD 7 = m k + 1, recorded in experiments/ic-bounds-m3).',
            reason='recorded m = 3 values: 7 at k = 2 (= 3k+1) on 105 of 116 cases (the 11 others have deg g = 5, LFD 6); Proposition 1 (experiments/ic-bounds-m3/README.md) gives LFD >= d + 1 when '
                   'NF(u g) has degree <= d for some chain variable u (automatic when d = m k) and x_i^(2^j) g is not proportional to g; with the kernel bound this forces LFD = m k + 1 when d = m k.',
            per_cell='m=4,k=2: d=8=mk: 9; m=4,k=3: d=12=mk: 13; m=5,k=2: d=10=mk: 11; m=3,k=2: 7. For d < m k (m = 3, k = 3, 4: d = 6, m k = 9, 12) (a) only says <= 10, 13 and does not predict the recorded 7.'),
  'b': dict(statement='LFD = deg(g) + 1 (one strict fall, at d + 1) for the prime direct systems at small depth (as m = 3: d = 6, LFD 7 for k = 2, 3, 4, 5 on 36 + 12 recorded instances). '
            'NOT true for m = 2 (d = 2... recorded LFD 4 at k <= 5, i.e. d + 2); for m = 3 it fails by depth 6 (experiments/ic-bounds-m3-depth: V_7 falls behind the ideal at k = 6, LFD >= 8).',
            reason='recorded m = 3 data', per_cell='LFD = d + 1 = 7 (m=3, all k = 2, 3, 4), 9 (m=4,k=2), 13 (m=4,k=3), 11 (m=5,k=2). On all these cells (a) and (b) agree except m = 3, k = 3, 4: they do NOT '
            'discriminate (d = m k for m >= 4 at k <= 3); the discriminating cell is m = 4, k = 4 (d = 12 < 16): (a) attained would be 17, (b) 13; that cell is in gpu/ and was not run.'),
  'c': dict(statement='binary direct systems: LFD <= min(r, d_in + 2) always (r = m k variables, d_in = maximal Boolean degree of the n equations), and LFD = min(r, d_in + 1) in at least 50 percent of the '
            'cases with d_in = r (scored cases).',
            reason='recorded binary m = 3 direct (experiments/ic-bounds-binary/results-m3-direct.jsonl, 695 rows, n = 5..9, k = 2..5): d_in = 6: LFD in 3..7 (max 7 = d_in + 1; at k = 2, r = d_in = 6, max 6 = r, never r + 1 = 7); '
                   'd_in = 3: max 5 = d_in + 2; at k = 2 with d_in = r = 6: LFD = r in 104 of 145 rows (72 percent), max 6; at k >= 3 with d_in = 6: LFD = d_in + 1 = 7 in 327 of 368 rows (89 percent). m = 2 is outside the rule (d_in = 2, LFD up to 5 at k = 10).',
            per_cell='see cases[...].pred_c_upper / pred_c_point')},
 global_predictions={
  'G1': 'every prime case has an exact certified LFD (no capped case) and it equals d + 1: 36 cases (18 m=3, 6 m=4,k=2, 6 m=4,k=3, 6 m=5,k=2).',
  'G2': 'no prime case exceeds m k + 1; every prime case with prop1_variables > 0 has LFD >= d + 1.',
  'G3': 'every prime case has exactly one strict fall, at D = d + 1 (V_D = 0 below d, V_d = span g, as recorded at m = 3).',
  'G4': 'binary: every case with an exact LFD has LFD <= min(r, d_in + 2); no binary case reaches r + 1.',
  'G5': 'binary: among exact cases with d_in = r, at least 50 percent have LFD = r.',
  'G6': 'm = 3 prime: LFD does not depend on k = 2, 3, 4 (7 everywhere).',
  'G7_gpu_m4k4': 'prime m = 4, k = 4 (65,536 monomials): LFD in [13, 17]; point prediction (b): 13 (the proved lower bound when Proposition 1 applies; prop1 variables: 12 of 16).',
  'G8_gpu_m5k3': 'prime m = 5, k = 3 (32,768 monomials; d = 15 = m k on the two surveyed instances): LFD = 16 forced by Proposition 1 and the kernel bound.',
  'G9_scoring': 'a prediction is scored a hit iff the exact certified value equals predicted_lfd (prime) / satisfies pred_c_upper (binary hard) / equals pred_c_point (binary point); capped or failed cases are scored as misses and listed.'},
 cases=cases)
json.dump(doc, open(OUTF, 'w'), indent=1, sort_keys=True)
print(sha(OUTF), len(cases))
