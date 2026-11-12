#!/usr/bin/env python3
"""Verifier driver: stratified subsample of experiments/ic-bounds-m3/cases.jsonl, recomputed with
scripts/ic_bounds_m3_verify_engine.py (independent engine), compared with the producer's rows.jsonl.

  python3 scripts/ic_bounds_m3_verify_run.py select          # -> verify/sample.json (seed 99173)
  python3 scripts/ic_bounds_m3_verify_run.py run --budget-min 20   # resumable -> verify/recheck.jsonl
  python3 scripts/ic_bounds_m3_verify_run.py compare         # -> verify/compare.json
"""
import json, time, random, argparse, sys, hashlib
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ic_bounds_m3_verify_engine as V

ROOT = Path(__file__).resolve().parents[1]
EXP = ROOT / 'experiments' / 'ic-bounds-m3'
VD = EXP / 'verify'
SEED = 99173


def load():
    cases = [json.loads(l) for l in open(EXP / 'cases.jsonl')]
    rows = {json.loads(l)['case_id']: json.loads(l) for l in open(EXP / 'rows.jsonl')}
    return cases, rows


def select():
    cases, rows = load()
    rng = random.Random(SEED)
    jobs = []
    for c in cases:
        k = c['k']; r = rows[c['case_id']]
        special = (r['generator_degree'] != 6 and k >= 2) or r['direct_lfd'] in (0, 3, 6) or \
                  (k == 3 and r['split_free_lfd'] == 5) or (k == 1 and r['split_free_lfd'] == 4)
        encs = []
        # direct: all at k<=3; k=4: all special + stratified 1 per (prime, curve_class)
        if k <= 3: encs.append('direct')
        encs.append('m2')
        if k <= 2: encs.append('split')
        jobs.append(dict(case_id=c['case_id'], encs=encs, special=special))
    # k=3 split: special + one random per (prime, curve_class); k=4 direct: special + one per (prime, class)
    by = {}
    for c in cases:
        by.setdefault((c['k'], c['prime'], c['curve_class']), []).append(c['case_id'])
    jm = {j['case_id']: j for j in jobs}
    for (k, p, cc), ids in sorted(by.items()):
        if k == 3:
            pick = rng.choice(ids); jm[pick]['encs'].append('split')
        if k == 4:
            pick = rng.choice(ids); jm[pick]['encs'].append('direct')
            if rng.random() < 0.3:
                pick2 = rng.choice(ids); jm[pick2]['encs'].append('split')
    for j in jobs:
        c = next(x for x in cases if x['case_id'] == j['case_id'])
        if j['special'] and c['k'] == 3 and 'split' not in j['encs']: j['encs'].append('split')
        if j['special'] and c['k'] == 4 and 'direct' not in j['encs'] and rows[j['case_id']]['generator_degree'] != 6:
            j['encs'].append('direct')
        j['encs'] = sorted(set(j['encs']))
    (VD / 'sample.json').write_text(json.dumps(dict(seed=SEED, jobs=jobs), indent=0))
    from collections import Counter
    cnt = Counter()
    for j in jobs:
        c = next(x for x in cases if x['case_id'] == j['case_id'])
        for e in j['encs']: cnt[(c['k'], e)] += 1
    print(sorted(cnt.items()))


def run(args):
    cases, rows = load()
    cm = {c['case_id']: c for c in cases}
    jobs = json.loads((VD / 'sample.json').read_text())['jobs']
    out = VD / 'recheck.jsonl'
    done = set()
    if out.exists():
        for l in open(out):
            d = json.loads(l); done.add((d['case_id'], d['enc']))
    terms = V.s4_terms()
    t_end = time.time() + 60 * args.budget_min
    todo = [(j['case_id'], e) for j in jobs for e in j['encs'] if (j['case_id'], e) not in done
            and (args.k is None or cm[j['case_id']]['k'] in args.k) and (args.enc is None or e in args.enc)
            and (args.primes is None or cm[j['case_id']]['prime'] in args.primes)]
    # cheap first
    cost = {'m2': 0, 'direct': 1, 'split': 2}
    todo.sort(key=lambda x: (cm[x[0]]['k'], cost[x[1]]))
    with open(out, 'a') as f:
        for cid, e in todo:
            if time.time() > t_end - 30:
                print('budget reached; rerun'); return
            c = cm[cid]
            t0 = time.time()
            try:
                r = V.exact_profile(e, c['prime'], c['a'], c['b'], c['target'][0], c['k'], terms,
                                    deadline=t_end, verify_closed=(c['k'] <= 2 or args.verify_closed))
            except TimeoutError:
                print('timeout on', cid, e); return
            r['case_id'] = cid; r['seconds'] = round(time.time() - t0, 2)
            f.write(json.dumps(r) + '\n'); f.flush()
            print(cid, c['k'], e, r['status'], r.get('lfd'), r['seconds'], flush=True)
    print('done')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['select', 'run'])
    ap.add_argument('--budget-min', type=float, default=20)
    ap.add_argument('--k', type=int, nargs='*', default=None)
    ap.add_argument('--enc', nargs='*', default=None)
    ap.add_argument('--verify-closed', action='store_true')
    ap.add_argument('--primes', type=int, nargs='*', default=None)
    a = ap.parse_args()
    select() if a.cmd == 'select' else run(a)
