"""Item 4: mid-scale check of csd/ic_bounds_exact10.py on stored depth-8 (or, with --depth-dir, depth-9) systems.

For each selected system of experiments/ic-chain-depth8/systems.jsonl (index into the file):
  * D = 4: split closure, expected 'closed-deficient' with dim V_4 = 4k^2 + 172k - 571 = 1061
    (experiments/ic-depth-step/growth.json, exact dense values);
  * D = 5: split closure from the generator rows, expected 'reached' with the rank and low
    dimension stored in experiments/ic-chain-depth8/results.jsonl (shift engine certificate);
  * optionally (--save-v5 DIR, first system only) the D = 5 basis in the probe's checkpoint format
    (basis.npy uint16 over SC.monomials(k, 5), state.json), to exercise
    scripts/ic_bounds_exact10_degree6.py at k = 8 (resume path included).
Writes experiments/ic-bounds-exact10/midscale-k8.jsonl.
"""
from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np

from csd import ic_bounds_exact10 as X


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--systems', type=int, nargs='+', default=[0, 14, 6])
    ap.add_argument('--save-v5', type=Path, default=None)
    ap.add_argument('--terms', type=int, default=4)
    ap.add_argument('--seed', type=int, default=20261107)
    ap.add_argument('--depth-dir', type=Path, default=ROOT / 'experiments/ic-chain-depth8')
    ap.add_argument('--out', type=Path, default=ROOT / 'experiments/ic-bounds-exact10/midscale-k8.jsonl')
    args = ap.parse_args()
    d8 = args.depth_dir
    systems = [json.loads(line) for line in open(d8 / 'systems.jsonl')]
    cases = [json.loads(line) for line in open(d8 / 'cases.jsonl')]
    results = {json.loads(line)['system_hash']: json.loads(line) for line in open(d8 / 'results.jsonl')}
    with args.out.open('w') as fo:
        for n, idx in enumerate(args.systems):
            sysd, case = systems[idx], cases[idx]
            assert sysd['system_hash'] == case['system_hash']
            k, p = sysd['k'], sysd['prime']
            g = {int(m): int(c) for m, c in sysd['generators'][0]}
            Z = list(sysd['coordinate_solutions'])
            g0 = sysd['generator_value']
            stored = results[sysd['system_hash']]
            cert = [q for q in stored.get('certificate_checks', []) if q['degree'] == 5][0]
            row = dict(index=idx, system_hash=sysd['system_hash'], k=k, prime=p, a=case['a'], b=case['b'],
                       target=case['target'], base_kind=case['base_kind'], scale=case['base_scale'], points=len(Z),
                       stored=dict(status=stored['status'], exact_lfd=stored['exact_lfd'], rank5=cert['rank'],
                                   target5=cert['target_rank'], low5=cert['low_dimension']))
            t = time.time()
            r4 = X.split_degree(g, k, p, 4, points=Z, generator_value=g0, seed=args.seed + 4, terms=args.terms)
            row['D4'] = {key: r4[key] for key in ('status', 'exact_dimension', 'low_dimension', 'target', 'rounds',
                                                  'verify_passes', 'closed')}
            row['D4']['seconds'] = round(time.time() - t, 1)
            t = time.time()
            r5, sc5 = X.split_degree(g, k, p, 5, points=Z, generator_value=g0, seed=args.seed + 5, terms=args.terms,
                                     keep=True)
            row['D5'] = {key: r5[key] for key in ('status', 'rank', 'exact_dimension', 'low_dimension', 'target',
                                                  'evaluation_onto_lt', 'rounds', 'verify_passes')}
            row['D5']['seconds'] = round(time.time() - t, 1)
            row['D5']['cols'] = X.monomial_count(k, 5)
            # exact LFD: V_4 exact (closed), V_5 reached with ev(M_<5) onto: fall at 5 iff dim V_4 < |M_<5| - |Z|
            fall5 = r4['exact_dimension'] < X.monomial_count(k, 4) - r5['eval_rank_lt']
            row['split_exact_lfd'] = 5 if (r5['status'] == 'reached' and r4['status'] == 'closed-deficient' and fall5) else None
            row['checks'] = dict(D4_dim_fit=r4['exact_dimension'] == 4 * k * k + 172 * k - 571,
                                 D4_closed_deficient=r4['status'] == 'closed-deficient',
                                 D5_rank=r5['exact_dimension'] == cert['rank'],
                                 D5_low=r5['low_dimension'] == cert['low_dimension'],
                                 D5_target=r5['target'] == cert['target_rank'],
                                 lfd=row['split_exact_lfd'] == stored['exact_lfd'])
            row['all_pass'] = all(row['checks'].values())
            row['maxrss_mb'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024
            if args.save_v5 is not None and n == 0:
                args.save_v5.mkdir(parents=True, exist_ok=True)
                rows = np.vstack(list(sc5.basis_full_rows()))
                from csd import elim
                out = []
                rank, _, _ = elim.fp_echelon(rows, p, out)
                assert rank == r5['exact_dimension']
                np.save(args.save_v5 / 'basis.npy', out[0].astype(np.uint16))
                (args.save_v5 / 'state.json').write_text(json.dumps(dict(phase='reached', rank=int(rank))))
                row['saved_v5'] = str(args.save_v5)
            fo.write(json.dumps(row) + '\n')
            fo.flush()
            print(json.dumps(dict(index=idx, points=len(Z), D4=row['D4']['status'], dim4=r4['exact_dimension'],
                                  D5=r5['status'], dim5=r5['exact_dimension'], stored5=cert['rank'],
                                  lfd=row['split_exact_lfd'], pass_=row['all_pass'], s4=row['D4']['seconds'],
                                  s5=row['D5']['seconds'], rss=row['maxrss_mb'])), flush=True)


if __name__ == '__main__':
    main()
