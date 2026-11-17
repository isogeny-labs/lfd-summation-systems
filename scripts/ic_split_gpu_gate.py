#!/usr/bin/env python3
"""Equality gate for csd/ic_bounds_split_gpu.py against csd/ic_bounds_split_m.py.

Instances are the certified ones of experiments/ic-bounds-split-m/cell-m*-k*.jsonl (first line of
each cell: prime, curve, target), rebuilt with split_generators_m.  For every (cell, D):

  live      the reference closure (InPlaceClosure, verification pass included) is run here and
            compared field by field with the candidate: dim V_D, low-degree profile, contains 1,
            verified closed;
  recorded  the candidate is compared with the reference's committed record of the same instance
            (used where the reference takes too long for this machine).

Candidate arms: device x split, e.g. cpu:split, cpu:plain, fake-gpu:split, gpu:split.
On cells of at most --space-max columns the spaces are compared too, not only their dimensions:
every candidate basis row, lifted to monomial coordinates, must lie in the reference V_D.

Negative controls (--controls): planted faults (an operator dropped, a wrong gather table, one
input coefficient changed) must make the gate fail on a closure without 1.  The changed coefficient
gives another system with the same dimensions, so only the comparison of spaces sees it; the record
says for each control whether the dimensions alone would have caught it.
Exit status 1 on any mismatch (or on a negative control that passes).
"""
import argparse, json, os, pathlib, sys, time

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from csd import ic_dreg as I                      # noqa: E402
from csd import ic_bounds_split_m as S            # noqa: E402
from csd import ic_bounds_split_gpu as G          # noqa: E402

CELLS = ROOT / 'experiments/ic-bounds-split-m'
FIELDS = ('monomials_le_D', 'dim_V_D', 'low_profile', 'contains_one')


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def instance(m, k, line=0):
    path = CELLS / f'cell-m{m}-k{k}.jsonl'
    if path.exists():
        r = json.loads(path.read_text().splitlines()[line])
        prof = {x['D']: x for x in r['profile']}
    elif (m, k) == (5, 3):          # the capped cell: prime and seed recorded, instance redrawn by the same rule
        sys.path.insert(0, str(ROOT / 'scripts'))
        import run_ic_split_m as RUN
        cap = json.loads((CELLS / 'cell-m5-k3-capped.json').read_text())
        E, inst = RUN.unit_instance(cap['p'], m, k, cap['seed'])
        r = dict(p=cap['p'], a=inst['a'], b=inst['b'], target=[inst['t']])
        prof = {x['D']: dict(x, contains_one=False, low_profile=None) for x in cap['closures_without_one']}
    else:
        raise SystemExit(f'no certified instance for ({m},{k})')
    E = I.PrimeCurve(r['p'], r['a'], r['b'])
    gens, _ = S.split_generators_m(E, r['target'][0], k, m)
    return r['p'], gens, prof


def compare(ref, rec):
    diff = [f for f in FIELDS if ref.get(f) is not None and ref[f] != rec[f]]
    if not rec['contains_one'] and rec['verified_closed'] is not True:
        diff.append('verified_closed')
    if rec['status'] not in ('one', 'closed'):
        diff.append('status:' + rec['status'])
    if rec['guard_hits']:
        diff.append('guard_hits')
    if rec.get('space_equal') is False:
        diff.append('space')
    return diff


def run_candidate(gens, m, k, p, D, arm, ref_gens=None, **kw):
    """Candidate record.  With ref_gens the spaces themselves are compared as well (small cases):
    the reference V_D is rebuilt and every candidate basis row, lifted to monomial coordinates, must
    lie in it; with equal dimensions the two spaces are then equal."""
    device, mode = arm.split(':')
    rec, eng = G.closure(gens, m, k, p, D, device=device, split=(mode == 'split'), **kw)
    if ref_gens is not None and not rec['contains_one']:
        sp = S.SplitSpace(m, k, D, m - 2)
        cl = S.InPlaceClosure(sp, p, 768)
        cl.add_inputs([sp.row(g, p) for g in ref_gens if S.generator_degree(g) <= D])
        cl.saturate()
        inside = all(cl.contains(v) for q, pt in enumerate(eng.parts) for v in eng.co.lift(pt.host_rows(), q, p))
        rec['space_equal'] = bool(inside and not cl.contains_one and cl.rank == rec['dim_V_D'])
    return rec


def perturb(gens, p):
    g0 = dict(gens[0])
    key = sorted(g0)[0]
    g0[key] = (g0[key] + 1) % p or 1
    return [g0] + [dict(g) for g in gens[1:]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--live', default='4,2:3-6;5,2:3-5;4,3:3-5;3,5:3-5;6,2:3-5;4,4:3-5',
                    help='cells and degree ranges compared with a live reference run')
    ap.add_argument('--recorded', default='', help='same syntax, compared with the committed reference record')
    ap.add_argument('--arms', default='cpu:split,cpu:plain,fake-gpu:split,fake-gpu:plain')
    ap.add_argument('--controls', action='store_true', help='also run the negative controls')
    ap.add_argument('--batch', type=int, default=1024)
    ap.add_argument('--space-max', type=int, default=5000, help='compare the spaces themselves up to this many columns')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    arms = a.arms.split(',')
    rows, failures = [], []

    def cells(spec):
        for item in filter(None, spec.split(';')):
            mk, dr = item.split(':')
            m, k = map(int, mk.split(','))
            lo, hi = (map(int, dr.split('-')) if '-' in dr else (int(dr), int(dr)))
            yield m, k, range(lo, hi + 1)

    for mode, spec in (('live', a.live), ('recorded', a.recorded)):
        for m, k, degs in cells(spec):
            p, gens, prof = instance(m, k)
            for D in degs:
                if mode == 'live':
                    t0 = time.time()
                    ref = S.closure(gens, m, k, p, D)
                    ref_s = round(time.time() - t0, 3)
                    if D in prof and any(prof[D].get(f) is not None and prof[D][f] != ref[f] for f in FIELDS):
                        failures.append(f'reference differs from its record m={m} k={k} D={D}')
                else:
                    if D not in prof:
                        raise SystemExit(f'no recorded closure for ({m},{k}) D={D}')
                    ref, ref_s = prof[D], prof[D].get('seconds')
                for arm in arms:
                    rec = run_candidate(gens, m, k, p, D, arm, batch=a.batch,
                                        ref_gens=gens if (mode == 'live' and ref['monomials_le_D'] <= a.space_max) else None)
                    diff = compare(ref, rec)
                    row = dict(mode=mode, m=m, k=k, p=p, D=D, arm=arm, equal=not diff, differences=diff,
                               reference={f: ref.get(f) for f in FIELDS}, reference_seconds=ref_s,
                               reference_verified=ref.get('verified_closed'),
                               candidate={f: rec[f] for f in FIELDS}, verified_closed=rec['verified_closed'],
                               space_equal=rec.get('space_equal'), part_rank=rec['part_rank'], symmetry=rec['symmetry'], seconds=rec['seconds'],
                               peak_device_gib=rec['peak_device_gib'], loadavg=os.getloadavg()[0])
                    rows.append(row)
                    log(f'{mode} m={m} k={k} D={D} {arm}: {"EQUAL" if not diff else "DIFFERENT " + str(diff)} '
                        f'dim={rec["dim_V_D"]} one={rec["contains_one"]} cand={rec["seconds"]}s ref={ref_s}s')
                    if diff:
                        failures.append(f'm={m} k={k} D={D} {arm}: {diff}')
    controls = []
    if a.controls:
        m, k, D = 4, 2, 5
        p, gens, prof = instance(m, k)
        ref = S.closure(gens, m, k, p, D)
        dev = arms[0].split(':')[0]
        tests = [('drop-operator', dict(fault='drop-operator'), gens, f'{dev}:split'),
                 ('wrong-gather', dict(fault='wrong-gather'), gens, f'{dev}:split'),
                 ('wrong-gather-plain', dict(fault='wrong-gather'), gens, f'{dev}:plain'),
                 ('perturbed-input-coefficient', {}, perturb(gens, p), f'{dev}:split'),
                 ('perturbed-input-coefficient-plain', {}, perturb(gens, p), f'{dev}:plain'),
                 ('no-fault (must pass)', {}, gens, f'{dev}:split')]
        for name, kw, gg, arm in tests:
            rec = run_candidate(gg, m, k, p, D, arm, ref_gens=gens, **kw)
            diff = compare(ref, rec)
            dims_only = [d for d in diff if d != 'space']
            want_fail = not name.startswith('no-fault')
            ok = bool(diff) == want_fail
            controls.append(dict(control=name, arm=arm, m=m, k=k, D=D, detected=bool(diff), differences=diff,
                                 detected_by_dimensions_alone=bool(dims_only),
                                 candidate={f: rec[f] for f in FIELDS}, reference={f: ref[f] for f in FIELDS},
                                 symmetry=rec['symmetry'], verified_closed=rec['verified_closed'], as_expected=ok))
            log(f'control {name} {arm}: detected={bool(diff)} {diff} dim={rec["dim_V_D"]} (reference {ref["dim_V_D"]})')
            if not ok:
                failures.append(f'control {name}: detected={bool(diff)}')
    out = dict(status='failed' if failures else 'passed', date=time.strftime('%Y-%m-%d'), arms=arms,
               comparisons=len(rows), equal=sum(r['equal'] for r in rows), failures=failures,
               rows=rows, controls=controls, engine_sha256=G.ENGINE_SHA256)
    if a.out:
        pathlib.Path(a.out).parent.mkdir(parents=True, exist_ok=True)
        pathlib.Path(a.out).write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(dict(status=out['status'], comparisons=len(rows), equal=out['equal'],
                          controls=[(c['control'], c['detected']) for c in controls], failures=failures)))
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
