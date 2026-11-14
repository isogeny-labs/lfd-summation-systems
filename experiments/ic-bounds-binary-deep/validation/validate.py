#!/usr/bin/env python3
"""Validation of csd/ic_bounds_bin_big.py against the existing engine (csd/ic_boolean_closure.py +
csd/ic_local_certificate.py) on systems drawn from the existing certified data.

For every selected system: the instance is rebuilt from its stored specification (basis, curve, target) with the
independent builders of the new driver; the system hash must equal the stored one and |Z| the stored count.  Then
for every degree D up to the certificate degree: dim V_D of the new engine (cpu, and the fake-gpu stand-in with
small tiles and batches) against (i) the stored certified record, (ii) a live run of the existing engine, and the
canonical RREF hash of the new basis against the RREF hash of the existing engine's basis re-reduced
independently.  Also fall profile, LFD, certificate degree, local checks and |Z|.

usage (limits of the workflow: ulimit -v 1500000, one thread):
  python3 validate.py select            write selection.json
  python3 validate.py run [--only i-j]  run the selected systems, append to records.jsonl
  python3 validate.py negative          negative controls  -> negative.json
  python3 validate.py summary           records.jsonl -> VALIDATION.md, summary.json
"""
import json
import os
import random
import sys
import time
import hashlib

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import numpy as np  # noqa: E402
import ic_bounds_bin_big_run as R  # noqa: E402
from csd import ic_bounds_bin_big as EB  # noqa: E402
from csd import ic_boolean_closure as C  # noqa: E402
from csd.ic_local_certificate import local_profile  # noqa: E402

EXP = os.path.join(ROOT, 'experiments')
SEED = 20261112


def load_candidates():
    res = {}
    for l in open(os.path.join(EXP, 'ic-curve-sweep/results.jsonl')):
        x = json.loads(l)
        res[x['system_hash']] = x
    out = []
    for l in open(os.path.join(EXP, 'ic-curve-sweep/cases.jsonl')):
        c = json.loads(l)
        if c['kind'] == 'prime' or c['m'] != 2 or c['system_hash'] not in res:
            continue
        x = res[c['system_hash']]
        if x.get('exact_lfd') is None:
            continue
        lc = x.get('local_checks') or []
        out.append(dict(source='ic-curve-sweep', n=c['n'], k=c['k'], modulus=c['modulus'], a=c['a'], b=c['b'], basis=c['basis'],
                        target=c['target'], target_class=c['target_class'], base_kind=c['base_kind'],
                        curve_class=c['kind'], frobenius_invariant=c['frobenius_invariant'], subfield=c.get('subfield', False),
                        system_hash=c['system_hash'], ref=dict(status=x['status'], lfd=x['exact_lfd'], falls=x['strict_falls'],
                                                              dims=x['quotient_ranks'], Z=x.get('point_count'),
                                                              cert=(lc[-1]['degree'] if lc else x['exact_lfd']))))
    for l in open(os.path.join(EXP, 'ic-bounds-binary/results-m2.jsonl')):
        x = json.loads(l)
        if x.get('exact_lfd') is None:
            continue
        out.append(dict(source='ic-bounds-binary', n=x['n'], k=x['k'], modulus=x['modulus'], a=x['a'], b=x['b'], basis=x['basis'],
                        target=x['target'], target_class=x['target_class'], base_kind=x['base_kind'], curve_class=x['curve_class'],
                        frobenius_invariant=x['frobenius_invariant'], subfield=x.get('subfield', False),
                        system_hash=x['system_hash'], ref=dict(status=x['status'], lfd=x['exact_lfd'], falls=x['strict_falls'],
                                                              dims=x['quotient_ranks'], Z=x['solutions'],
                                                              cert=x.get('certificate_degree', x['exact_lfd']))))
    seen, uniq = set(), []
    for c in out:
        if c['system_hash'] not in seen:
            seen.add(c['system_hash'])
            uniq.append(c)
    return uniq


def spec_of(c, name):
    return dict(name=name, n=c['n'], k=c['k'], base='given', curve='given', a=c['a'], b=c['b'], basis=c['basis'],
                target='given', target_point=c['target'], seed=0)


def classify(inst):
    if inst['structure']['const']:
        return 'const'
    if inst['n_points'] == 0:
        return 'unit'
    return 'order-two' if inst['target'][0] == 0 else 'affine'


def select():
    cands = load_candidates()
    rng = random.Random(SEED)
    rng.shuffle(cands)
    # classify by building instances (cheap for k <= 8); k = 9, 10 only from the stored class labels
    want = []

    def build(c, name):
        inst = R.build_instance(spec_of(c, name))
        assert inst['system_hash'] == c['system_hash'], ('hash', name)
        return inst
    buckets = {}
    for c in cands:
        if c['k'] > 8 or c['k'] < 4:
            continue
        key = (c['k'], c['base_kind'], c['curve_class'])
        buckets.setdefault(key, []).append(c)
    chosen, labels = [], {}
    # (1) k = 4..8: per k, both base kinds, both curve families: scan until each class label is seen
    for k in range(4, 9):
        for bk in ('random-subspace', 'frobenius-invariant'):
            for cc in ('koblitz', 'binary'):
                pool = buckets.get((k, bk, cc), [])
                if not pool:
                    pool = buckets.get((k, bk, 'ordinary'), [])
                got = 0
                for c in pool[:60]:
                    inst = build(c, 'v%d' % len(chosen))
                    lab = classify(inst)
                    cnt = sum(1 for x in chosen if x['k'] == k and x['bk'] == bk and x['label'] == lab and x['cc'] == cc)
                    if cnt == 0 and (lab in ('affine', 'order-two') or k <= 7):
                        chosen.append(dict(c, label=lab, bk=bk, cc=cc))
                        got += 1
                    if got >= (3 if lab == 'affine' else 4) and len({x['label'] for x in chosen if x['k'] == k and x['bk'] == bk and x['cc'] == cc}) >= 2:
                        break
    # (2) classes that need an explicit hunt: const and unit (non-const, |Z| = 0) at k = 5..8
    for lab in ('const', 'unit'):
        for k in (5, 6, 7, 8):
            for c in cands:
                if c['k'] != k:
                    continue
                if lab == 'unit' and c['ref']['status'] != 'exact-unit':
                    continue
                inst = build(c, 'v%d' % len(chosen))
                if classify(inst) == lab and not any(x['label'] == lab and x['k'] == k for x in chosen):
                    chosen.append(dict(c, label=lab, bk=c['base_kind'], cc=c['curve_class']))
                    break
    # (3) k = 9, 10: a few (the existing engine needs 5-100 s each)
    for k, n_take in ((9, 3), (10, 2)):
        pool = [c for c in cands if c['k'] == k and c['source'] == 'ic-bounds-binary']
        for c in pool:
            if n_take == 0:
                break
            if c['ref']['lfd'] <= 4 and c['target_class'] != 'planted-rational-decomposition':
                chosen.append(dict(c, label='n/a', bk=c['base_kind'], cc=c['curve_class']))
                n_take -= 1
    sel = []
    for i, c in enumerate(chosen):
        c = dict(c, name='v%02d' % i)
        sel.append(c)
    json.dump(sel, open(os.path.join(HERE, 'selection.json'), 'w'), indent=1)
    from collections import Counter
    print(len(sel), Counter((x['label'], x['cc']) for x in sel), Counter(x['k'] for x in sel))


# subfield systems: V = F_(2^k) inside F_(2^2k), order-two target (0, 1), b = 1, a in {0, 1}
SUBFIELD = [dict(name='sf%d_%d_a%d' % (2 * k, k, a), n=2 * k, k=k, base='subfield', curve='koblitz%d' % a, target='order-two', seed=0)
            for k in (4, 6, 8) for a in (0, 1)]


def old_hash(eqs, r, D, max_bytes):
    sp = C.boolean_degree_closure(eqs, r, D, max_columns=10 ** 6, max_bytes=max_bytes)
    N = len(sp.masks)
    if sp.contains_one:
        return EB.canonical_hash_unit(r, D, N), N
    EB.use_device('cpu')
    Acc, piv = EB.rref_block(EB.DV.to_dev(sp.basis.copy()), 1 << 14)

    def rows_iter(free):
        bits = EB.unpack(Acc)[:, free] if len(free) else np.zeros((len(piv), 0), dtype=np.uint8)
        for row in bits:
            yield np.packbits(row, bitorder='little').tobytes()
    return EB.canonical_hash(r, D, N, piv, rows_iter), int(len(piv))


def settings(k):
    if k <= 5:
        return dict(batch=3, tile_rows=8, col_block=64, compact_frac=0.05)
    if k <= 8:
        return dict(batch=37, tile_rows=64, col_block=192, compact_frac=0.1, tiles_on='host')
    return dict(batch=200, tile_rows=512, col_block=1024, compact_frac=0.2)


def run_one(inst, ref=None, label=None, source=None, hash_limit_k=8):
    k, r = inst['k'], inst['variables']
    eqs, pts = inst['equations'], inst['points']
    rec = dict(name=inst['name'], source=source, n=inst['n'], k=k, cls=classify(inst), base_kind=inst['spec'].get('base'),
               curve=inst['curve'], frobenius_invariant=inst['frobenius_invariant'], scaled_subfield=inst['scaled_subfield'],
               Z=len(pts), system_hash=inst['system_hash'], e=inst['structure']['e'])
    big = k >= 9
    maxc, maxb = (14000, 700_000_000) if big else (200000, 1_000_000_000)
    t = time.time()
    old = local_profile(eqs, r, pts, max_degree=7, max_columns=maxc, max_bytes=maxb)
    rec['old_seconds'] = round(time.time() - t, 2)
    rec['old'] = dict(status=old['status'], lfd=old['exact_lfd'], falls=old['strict_falls'], dims=old['quotient_ranks'],
                      cert=(old['exact_lfd'] if old['status'] == 'exact-unit' else old.get('local_upper_bound')),
                      local_checks=old['local_checks'])
    if old['exact_lfd'] is None:
        rec['skipped'] = 'existing engine capped: ' + str(old.get('cap_reason'))
        return rec
    Dmax = rec['old']['cert']
    # stored certified record
    if ref:
        rec['stored'] = ref
        rec['stored_equal'] = (ref['lfd'] == old['exact_lfd'] and ref['falls'] == old['strict_falls']
                               and ref['dims'][:len(old['quotient_ranks'])] == old['quotient_ranks']
                               and (ref['Z'] is None or ref['Z'] == len(pts)))
    # old hashes
    old_hashes = []
    if k <= hash_limit_k or (big and Dmax <= 4):
        for D in range(0, Dmax + 1):
            old_hashes.append(old_hash(eqs, r, D, maxb)[0])
    rec['old_hash_checked'] = bool(old_hashes)
    for mode, kw in (('cpu', dict(device='cpu')), ('fake-gpu', dict(device='fake-gpu', **settings(k)))):
        t = time.time()
        new = EB.local_profile_bin(eqs, r, pts, Dmax + 1, kw)
        sec = round(time.time() - t, 2)
        hashes = [x['rref_sha256'] for x in new['per_degree']]
        rec[mode] = dict(status=new['status'], lfd=new['exact_lfd'], cert=new['certificate_degree'], falls=new['strict_falls'],
                         dims=new['quotient_ranks'], local_checks=new['local_checks'], seconds=sec,
                         closed=all(x['closed'] for x in new['per_degree']),
                         repairs=sum(x['repairs'] for x in new['per_degree']),
                         compactions=sum(x['compactions'] for x in new['per_degree']),
                         hashes=[h[:16] for h in hashes], settings=kw)
        rec[mode]['dims_equal_per_D'] = new['quotient_ranks'] == old['quotient_ranks']
        rec[mode]['falls_equal'] = new['strict_falls'] == old['strict_falls']
        rec[mode]['lfd_equal'] = new['exact_lfd'] == old['exact_lfd']
        rec[mode]['cert_equal'] = new['certificate_degree'] == rec['old']['cert']
        rec[mode]['local_checks_equal'] = new['local_checks'] == old['local_checks']
        rec[mode]['hashes_equal_old'] = (hashes[:len(old_hashes)] == old_hashes) if old_hashes else None
        rec[mode]['Z_equal'] = len(pts) == (ref['Z'] if ref and ref['Z'] is not None else len(pts))
        rec[mode]['full_hashes'] = hashes
    rec['cpu_equals_fake_hashes'] = rec['cpu']['full_hashes'] == rec['fake-gpu']['full_hashes']
    for m in ('cpu', 'fake-gpu'):
        rec[m].pop('full_hashes')
    keys = ('dims_equal_per_D', 'falls_equal', 'lfd_equal', 'cert_equal', 'local_checks_equal', 'Z_equal')
    rec['all_ok'] = bool(all(rec[m][kk] for m in ('cpu', 'fake-gpu') for kk in keys)
                         and all(rec[m]['hashes_equal_old'] in (True, None) for m in ('cpu', 'fake-gpu'))
                         and rec['cpu_equals_fake_hashes'] and rec['cpu']['closed'] and rec['fake-gpu']['closed']
                         and rec.get('stored_equal', True))
    return rec


def cmd_run(only=None):
    sel = json.load(open(os.path.join(HERE, 'selection.json')))
    items = [('sel', c) for c in sel] + [('sub', s) for s in SUBFIELD]
    lo, hi = (0, len(items)) if not only else map(int, only.split('-'))
    path = os.path.join(HERE, 'records.jsonl')
    done = set()
    if os.path.exists(path):
        done = {json.loads(l)['name'] for l in open(path) if l.strip()}
    for i, (kind, c) in enumerate(items):
        if not lo <= i < hi:
            continue
        name = c['name']
        if name in done:
            continue
        if kind == 'sel':
            inst = R.build_instance(spec_of(c, name))
            assert inst['system_hash'] == c['system_hash']
            rec = run_one(inst, ref=c['ref'], label=c['label'], source=c['source'])
            rec['label_in_selection'] = c['label']
            rec['stored_curve_class'] = c['curve_class']
            rec['stored_base_kind'] = c['base_kind']
        else:
            inst = R.build_instance(c)
            rec = run_one(inst, source='subfield-construction')
        with open(path, 'a') as f:
            f.write(json.dumps(rec, sort_keys=True) + '\n')
        print(i, name, rec['cls'], 'k=%d' % rec['k'], 'Z=%d' % rec['Z'], 'ok=', rec.get('all_ok'), rec.get('skipped', ''),
              'old %.1fs cpu %.1fs fake %.1fs' % (rec['old_seconds'], rec.get('cpu', {}).get('seconds', 0), rec.get('fake-gpu', {}).get('seconds', 0)), flush=True)


def cmd_negative():
    """(1) withheld shifts of one row: with --no-repair the sweep must fail (closed false); with repair it must
    repair and reach the identical hash.  (2) a basis row deleted from a finished closure must be caught by the
    sweep.  (3) a corrupted checkpoint basis (one bit flipped in a row) must be caught."""
    sel = json.load(open(os.path.join(HERE, 'selection.json')))
    out = []
    pick = [c for c in sel if c['label'] == 'affine' and 5 <= c['k'] <= 7][:3] + [c for c in sel if c['label'] == 'order-two' and c['k'] == 6][:1]
    for c in pick:
        inst = R.build_instance(spec_of(c, c['name']))
        r, eqs, pts = inst['variables'], inst['equations'], inst['points']
        D = c['ref']['cert'] if c['ref']['cert'] else c['ref']['lfd']
        base = EB.Engine(eqs, r, D, device='cpu').run()
        row = dict(name=c['name'], k=c['k'], D=D, clean_hash=base['rref_sha256'], clean_rank=base['rank'], withhold={}, deleted_row={})
        for wh in (0, 1, 2, 5):
            e1 = EB.Engine(eqs, r, D, device='cpu', withhold=wh, repair=False)
            r1 = e1.run()
            e2 = EB.Engine(eqs, r, D, device='fake-gpu', withhold=wh, repair=True, batch=5, tile_rows=8, col_block=64)
            r2 = e2.run()
            row['withhold'][str(wh)] = dict(no_repair_closed=r1['closed'], no_repair_rank=r1['rank'], rank_deficit=base['rank'] - r1['rank'],
                                            repaired_closed=r2['closed'], repairs=r2['repairs'], repaired_hash_equal=r2['rref_sha256'] == base['rref_sha256'])
        # delete a low basis row of a finished closure
        eng = EB.Engine(eqs, r, D, device='cpu')
        eng.run()
        low = np.flatnonzero(eng.row_piv[:eng.nrows] >= eng.ntop)
        victim = int(low[len(low) // 2])
        # remove its contribution: zero the row (it stays a pivot with zero B, i.e. e_piv alone): the space changes
        t = eng.tiles[victim // eng.TR]
        t['a'][victim % eng.TR] ^= t['a'][victim % eng.TR]          # row becomes e_pivot (B = 0): a wrong but valid RREF row
        eng.repair = False
        row['deleted_row'] = dict(victim=victim, sweep_clean_after_corruption=bool(eng.sweep()))
        out.append(row)
        print(json.dumps(row), flush=True)
    json.dump(out, open(os.path.join(HERE, 'negative.json'), 'w'), indent=1)


def cmd_summary():
    recs = [json.loads(l) for l in open(os.path.join(HERE, 'records.jsonl')) if l.strip()]
    from collections import Counter
    lines = ['# Validation of csd/ic_bounds_bin_big.py against the existing engine', '',
             'Generated by validate.py summary from records.jsonl (%d systems).' % len(recs), '',
             '| system | src | n | k | class | base | curve | Z | LFD (old/cpu/fake) | cert D | falls | dims equal (all D) | RREF hash equal | all |',
             '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
    for r_ in recs:
        if r_.get('skipped'):
            lines.append('| %s | %s | %d | %d | %s | | | %d | skipped (%s) | | | | | |' % (r_['name'], r_['source'], r_['n'], r_['k'], r_['cls'], r_['Z'], r_['skipped'][:40]))
            continue
        hc = [r_[m]['hashes_equal_old'] for m in ('cpu', 'fake-gpu')]
        lines.append('| %s | %s | %d | %d | %s | %s | %s | %d | %s/%s/%s | %s | %s | %s/%s | %s | %s |' % (
            r_['name'], (r_['source'] or '')[:12], r_['n'], r_['k'], r_['cls'], r_['base_kind'] if r_['base_kind'] != 'given' else r_.get('stored_base_kind', '')[:6],
            r_['curve']['family'], r_['Z'], r_['old']['lfd'], r_['cpu']['lfd'], r_['fake-gpu']['lfd'], r_['old']['cert'], r_['old']['falls'],
            'yes' if r_['cpu']['dims_equal_per_D'] else 'NO', 'yes' if r_['fake-gpu']['dims_equal_per_D'] else 'NO',
            'yes' if all(h is True for h in hc) else ('n/a' if hc[0] is None else 'NO'), 'yes' if r_['all_ok'] else 'NO'))
    ok = [r_ for r_ in recs if r_.get('all_ok')]
    summ = dict(systems=len(recs), all_ok=len(ok), skipped=sum(1 for r_ in recs if r_.get('skipped')),
                by_class=dict(Counter(r_['cls'] for r_ in recs)), by_k=dict(Counter(r_['k'] for r_ in recs)),
                by_curve=dict(Counter(r_['curve']['family'] for r_ in recs)),
                frobenius_invariant=sum(1 for r_ in recs if r_['frobenius_invariant']),
                subfield=sum(1 for r_ in recs if r_['scaled_subfield']),
                degrees_compared=sum(len(r_['old']['dims']) for r_ in recs if not r_.get('skipped')),
                hash_checked_systems=sum(1 for r_ in recs if r_.get('old_hash_checked')),
                stored_equal=sum(1 for r_ in recs if r_.get('stored_equal')), with_stored=sum(1 for r_ in recs if 'stored' in r_))
    lines += ['', '```', json.dumps(summ, indent=1), '```']
    open(os.path.join(HERE, 'VALIDATION_TABLE.md'), 'w').write('\n'.join(lines) + '\n')
    json.dump(summ, open(os.path.join(HERE, 'summary.json'), 'w'), indent=1)
    print(json.dumps(summ))


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'select':
        select()
    elif cmd == 'run':
        cmd_run(sys.argv[3] if len(sys.argv) > 3 and sys.argv[2] == '--only' else None)
    elif cmd == 'negative':
        cmd_negative()
    elif cmd == 'summary':
        cmd_summary()
