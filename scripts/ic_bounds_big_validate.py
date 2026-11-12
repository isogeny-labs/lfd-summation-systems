"""LOCAL validation of csd/ic_bounds_big.py against the reference engine csd/ic_bounds_d5_verify.py.

The reduced row echelon basis of a subspace is unique for a fixed column order, so the big engine must
reproduce the reference basis hashes exactly; every check below compares rank, profile, part ranks,
contains_one and both basis hashes.  Run with light limits, one thread, scratch under --work:

  ulimit -v 1500000; OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 scripts/ic_bounds_big_validate.py small18 --work DIR

Subcommands: small18, refs (k9, k10 saved reference hashes), p65537, resume (SIGTERM / kill -9), negative
(withheld shifts), stopone, threads, rates.  Each writes experiments/ic-bounds-k14/validation/<name>.json.
"""
import argparse
import shlex
import hashlib
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from csd import ic_bounds_d5_verify as V  # noqa: E402

RUN = ROOT / 'scripts/ic_bounds_big_run.py'
ENG = ROOT / 'csd/ic_bounds_big.py'
OUTDIR = ROOT / 'experiments/ic-bounds-k14/validation'
SMALL = json.load(open(ROOT / 'experiments/ic-bounds-d5/verify/small.json'))
SMALL_INST = {'probe-unit': ('unit-12289', 12289, 1, 0, 10, False), 'proper-13313': ('p13313', 13313, 1, 0, 11, True),
              'random-15361': ('p15361', 15361, 4021, 9917, 12, True)}
SRC = {n: hashlib.sha256(open(p, 'rb').read()).hexdigest() for n, p in
       (('engine', ENG), ('driver', RUN), ('validate', Path(__file__)))}
MODES = {
    'plain': ['--device', 'cpu'],
    'fake-gpu': ['--device', 'fake-gpu'],
    'tiny': ['--device', 'cpu', '--tile-rows', '5', '--batch', '3', '--col-block', '33', '--compact-frac', '0.02'],
    'tiny-fake-gpu': ['--device', 'fake-gpu', '--tile-rows', '5', '--batch', '3', '--col-block', '33',
                      '--compact-frac', '0.02'],
}


def save(name, obj):
    OUTDIR.mkdir(parents=True, exist_ok=True)
    obj = dict(obj, source_sha256=SRC)
    (OUTDIR / f'{name}.json').write_text(json.dumps(obj, indent=1) + '\n')
    print('saved', name, flush=True)


def inst_file(work, tag, k, prime, a, b, seed, want, nonres=False, D=5):
    path = Path(work) / f'{tag}-k{k}.json'
    if not path.exists():
        cmd = [sys.executable, str(RUN), 'instance', '--k', str(k), '--prime', str(prime), '--a', str(a), '--b', str(b),
               '--seed', str(seed), '--D', str(D), '--out', str(path)]
        if want:
            cmd.append('--want')
        if nonres:
            cmd.append('--unit-nonresidue')
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    return str(path)


def load_g(path):
    d = json.load(open(path))
    return d, {tuple(int(x) for x in key.split(',')): v for key, v in d['generator'].items()}


def reference(path, D, gens, split=True):
    d, g = load_g(path)
    C = V.Closure(g, d['k'], d['prime'], D, split=split, gens=gens)
    C.run()
    prof, parts = C.profile()
    return dict(rank=C.rank(), profile=prof, part_profiles=parts, part_ranks=[P.rank for P in C.parts],
                contains_one=C.contains_one(), basis_sha256=[P.sha() for P in C.parts])


def big(path, D, extra, out, env=None, timeout=3000):
    cmd = [sys.executable, str(RUN), 'closure', '--instance', path, '--D', str(D), '--out', str(out)] + extra
    t = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, env=dict(os.environ, **(env or {})), timeout=timeout)
    rec = json.load(open(out)) if os.path.exists(out) else None
    return r.returncode, rec, round(time.time() - t, 1), r


FIELDS = ('rank', 'profile', 'part_profiles', 'part_ranks', 'contains_one', 'basis_sha256')


def same(ref, rec):
    return all(ref[f] == rec[f] for f in FIELDS)


def cmd_small18(a):
    work = Path(a.work)
    work.mkdir(parents=True, exist_ok=True)
    rows = []
    for ref in SMALL['rows']:
        tag, prime, aa, bb, seed, want = SMALL_INST[ref['instance']]
        k, D = ref['k'], ref['D']
        path = inst_file(work, tag, k, prime, aa, bb, seed, want)
        r = reference(path, D, 'g')
        agree_small = (r['rank'] == ref['ref_dense']['dim'] and r['profile'] == ref['ref_dense']['profile']
                       and r['contains_one'] == ref['ref_dense']['contains_one'])
        row = dict(instance=ref['instance'], k=k, D=D, gens='g', ref_rank=r['rank'], ref_profile=r['profile'],
                   ref_hashes=[h[:16] for h in r['basis_sha256']], reference_agrees_with_small_json=agree_small, modes={})
        for mode, flags in MODES.items():
            out = work / 'out.json'
            if out.exists():
                out.unlink()
            rc, rec, sec, _ = big(path, D, flags + ['--gens', 'g'], out)
            ok = rc == 0 and rec is not None and same(r, rec) and rec['closed'] and rec['rref_ok']
            row['modes'][mode] = dict(match=bool(ok), closed=rec and rec['closed'], seconds=sec, rank=rec and rec['rank'])
        row['all_match'] = all(m['match'] for m in row['modes'].values())
        print(json.dumps(dict(instance=ref['instance'], k=k, D=D, all=row['all_match'], small=agree_small)), flush=True)
        rows.append(row)
    save('small18', dict(rows=rows, all_match=all(r['all_match'] for r in rows),
                         all_ref_agree=all(r['reference_agrees_with_small_json'] for r in rows), modes=list(MODES)))


def cmd_refs(a):
    work = Path(a.work)
    rows = []
    for name, tag, k in (('k9-p12289-s10-D5.json', 'unit-12289', 9), ('k10-p12289-s10-D5.json', 'unit-12289', 10)):
        ref = json.load(open(ROOT / 'experiments/ic-bounds-d5/verify' / name))
        path = inst_file(work, tag, k, 12289, 1, 0, 10, False)
        for mode in a.modes.split(','):
            out = work / 'out.json'
            if out.exists():
                out.unlink()
            rc, rec, sec, _ = big(path, 5, MODES[mode] + ['--work-gb', '0.4', '--batch', '256'], out)
            ok = (rc == 0 and rec['rank'] == ref['rank'] and rec['profile'] == ref['profile']
                  and rec['part_profiles'] == ref['part_profiles'] and rec['part_ranks'] == ref['part_ranks']
                  and rec['contains_one'] == ref['contains_one'] and rec['basis_sha256'] == ref['basis_sha256'])
            row = dict(case=name, mode=mode, match=bool(ok), rank=rec['rank'], profile=rec['profile'],
                       closed=rec['closed'], seconds=sec, madds=rec['madds'], madds_dense=rec['madds_dense'],
                       madds_detail=rec['madds_detail'], seconds_detail=rec['seconds_detail'], maxrss_mb=rec['maxrss_mb'],
                       sweep_random=rec['sweep_random'], sweep_full=rec['sweep_full'], ref_seconds=ref['seconds'],
                       ref_madds=ref['madds'], hashes=[h[:16] for h in rec['basis_sha256']])
            print(json.dumps({k: row[k] for k in ('case', 'mode', 'match', 'seconds')}), flush=True)
            rows.append(row)
    save('refs-k9-k10', dict(rows=rows, all_match=all(r['match'] for r in rows)))


def cmd_p65537(a):
    work = Path(a.work)
    rows = []
    for tag, aa, bb, seed, want, nonres in (('unit-65537', 1, 0, 10, False, True), ('proper-65537', 4021, 9917, 12, True, False)):
        for k, D in ((6, 5), (7, 5), (8, 5), (6, 6), (7, 6)):
            path = inst_file(work, tag, k, 65537, aa, bb, seed, want, nonres)
            r = reference(path, D, 'all')
            row = dict(case=f'{tag}-k{k}-D{D}', ref_rank=r['rank'], ref_profile=r['profile'], ref_one=r['contains_one'],
                       ref_hashes=[h[:16] for h in r['basis_sha256']], modes={})
            for mode in ('plain', 'fake-gpu', 'tiny'):
                out = work / 'out.json'
                if out.exists():
                    out.unlink()
                rc, rec, sec, _ = big(path, D, MODES[mode], out)
                row['modes'][mode] = dict(match=bool(rc == 0 and same(r, rec) and rec['closed']), seconds=sec)
                if 'points' in rec:
                    row['modes'][mode].update(points=rec['points'], deficient=rec.get('deficient'),
                                              sound=rec.get('soundness_vanish_on_Z'), eval_rank=rec.get('eval_rank_below_D'))
            row['all_match'] = all(m['match'] for m in row['modes'].values())
            print(json.dumps(dict(case=row['case'], all=row['all_match'], rank=r['rank'], one=r['contains_one'])), flush=True)
            rows.append(row)
    save('p65537', dict(rows=rows, all_match=all(r['all_match'] for r in rows), hash_dtype='uint32 (p >= 65536)'))


def cmd_resume(a):
    work = Path(a.work)
    path = inst_file(work, 'p15361', a.k, 15361, 4021, 9917, 12, True)
    D = a.D
    ref = reference(path, D, 'all')
    common = ['--device', a.device, '--tile-rows', '40', '--batch', '16', '--sweep-chunk', '16', '--log-seconds', '0.5'] + shlex.split(a.extra)
    rows = {}
    # uninterrupted
    out = work / 'u.json'
    out.unlink(missing_ok=True)
    rc, rec, sec, _ = big(path, D, common, out)
    rows['uninterrupted'] = dict(match=same(ref, rec), seconds=sec)
    base_rec = rec
    # SIGTERM then resume
    ck = work / 'ck-term'
    subprocess.run(['rm', '-rf', str(ck)])
    out = work / 'term.json'
    out.unlink(missing_ok=True)
    cmd = [sys.executable, str(RUN), 'closure', '--instance', path, '--D', str(D), '--out', str(out), '--ckpt', str(ck),
           '--ckpt-minutes', '0.05'] + common
    sig_results = []
    for delay in (float(x) for x in a.delays.split(',')):
        p = subprocess.Popen(cmd + (['--resume'] if sig_results else []), stdout=subprocess.PIPE, text=True)
        time.sleep(delay)
        p.send_signal(signal.SIGTERM)
        so, _ = p.communicate(timeout=600)
        sig_results.append(p.returncode)
        if p.returncode == 0:
            break
    p = subprocess.run(cmd + ['--resume'], capture_output=True, text=True)
    rec = json.load(open(out))
    rows['sigterm_resume'] = dict(match=same(ref, rec), same_rounds=rec['rounds'] == base_rec['rounds'],
                                  same_madds=rec['madds'] == base_rec['madds'], exit_codes_of_interrupted_runs=sig_results,
                                  final_exit=p.returncode, resumes=rec['resumes'], closed=rec['closed'], repairs=rec['repairs'],
                                  ckpt_left=sorted(os.listdir(ck)))
    # SIGTERM triggered by the log: once during the generator rows, once in the full sweep
    for stage in ('gens', 'sweep_full'):
        ck = work / f'ck-{stage}'
        subprocess.run(['rm', '-rf', str(ck)])
        out = work / f'{stage}.json'
        out.unlink(missing_ok=True)
        cmd = [sys.executable, str(RUN), 'closure', '--instance', path, '--D', str(D), '--out', str(out), '--ckpt', str(ck)] \
            + [c for c in common if c not in ('--log-seconds', '0.5')] + ['--log-seconds', '0']
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
        sent = False
        for line in p.stdout:
            if not sent and f'"stage": "{stage}"' in line and (stage == 'gens' or '"sweep_pos": 32' in line):
                p.send_signal(signal.SIGTERM)
                sent = True
        p.wait()
        code = p.returncode
        p2 = subprocess.run(cmd + ['--resume'], capture_output=True, text=True)
        rec = json.load(open(out))
        rows[f'sigterm_in_{stage}'] = dict(match=same(ref, rec), same_rounds=rec['rounds'] == base_rec['rounds'], interrupted_exit=code, final_exit=p2.returncode,
                                           resumes=rec['resumes'], closed=rec['closed'], signal_sent=sent)
    # kill -9 during a checkpoint (second checkpoint is slowed down by the test hook)
    ck = work / 'ck-kill'
    subprocess.run(['rm', '-rf', str(ck)])
    out = work / 'kill.json'
    out.unlink(missing_ok=True)
    cmd = [sys.executable, str(RUN), 'closure', '--instance', path, '--D', str(D), '--out', str(out), '--ckpt', str(ck),
           '--ckpt-minutes', '0.03'] + common
    env = dict(os.environ, ICBIG_TEST_CKPT_SLEEP='30', ICBIG_TEST_CKPT_SLEEP_FROM='2')
    p = subprocess.Popen(cmd, env=env, stdout=subprocess.DEVNULL)
    t0 = time.time()
    killed_in_tmp = False
    while time.time() - t0 < 600:
        time.sleep(0.2)
        if os.path.isdir(ck):
            names = os.listdir(ck)
            if any(n.startswith('.tmp-00000002') for n in names) and any(n.startswith('ck-00000001') for n in names):
                if len(os.listdir(os.path.join(ck, [n for n in names if n.startswith('.tmp-00000002')][0]))) > 0:
                    p.send_signal(signal.SIGKILL)
                    killed_in_tmp = True
                    break
        if p.poll() is not None:
            break
    p.wait()
    left = sorted(os.listdir(ck)) if os.path.isdir(ck) else []
    p2 = subprocess.run(cmd + ['--resume'], capture_output=True, text=True)
    rec = json.load(open(out))
    rows['kill9_during_checkpoint'] = dict(match=same(ref, rec), same_rounds=rec['rounds'] == base_rec['rounds'], killed_during_second_checkpoint=killed_in_tmp,
                                           ckpt_dir_after_kill=left, final_exit=p2.returncode, resumes=rec['resumes'],
                                           closed=rec['closed'])
    save(f'resume-k{a.k}-D{D}' + a.tag_suffix, dict(extra=a.extra, k=a.k, D=D, device=a.device, rows=rows, all_match=all(r['match'] and r.get('same_rounds', True) for r in rows.values())))


def cmd_negative(a):
    work = Path(a.work)
    path = inst_file(work, 'p15361', a.k, 15361, 4021, 9917, 12, True)
    ref = reference(path, a.D, a.gens)
    rows = []
    base = ['--gens', a.gens, '--tile-rows', '40', '--batch', '16', '--device', a.device] + shlex.split(a.extra)
    for part, idx, tr in ((0, 0, 6), (0, 1, 6), (0, 2, 6), (0, 3, 6), (0, 5, 6), (0, 8, 6), (0, 13, 6), (0, 21, 6),
                          (0, 40, 6), (1, 0, 6), (1, 2, 6), (1, 7, 6), (1, 30, 6), (0, 0, 0), (1, 0, 0)):
        w = f'{part}:{idx}'
        out = work / 'neg.json'
        out.unlink(missing_ok=True)
        rc, rec, sec, r = big(path, a.D, ['--withhold', w, '--no-repair', '--t-random', str(tr)] + base, out)
        out2 = work / 'neg2.json'
        out2.unlink(missing_ok=True)
        rc2, rec2, sec2, _ = big(path, a.D, ['--withhold', w, '--t-random', str(tr)] + base, out2)
        rows.append(dict(withheld=w, t_random=tr, no_repair_exit=rc, no_repair_closed=rec['closed'],
                         sweep_random=rec['sweep_random'], sweep_full=rec['sweep_full'],
                         hash_wrong_without_repair=rec['basis_sha256'] != ref['basis_sha256'], rank_without_repair=rec['rank'],
                         ref_rank=ref['rank'], repair_exit=rc2, repair_repairs=rec2['repairs'], repair_rows=rec2['repair_rows'],
                         repair_closed=rec2['closed'], repair_hash_match=same(ref, rec2)))
        print(json.dumps({k: rows[-1][k] for k in ('withheld', 't_random', 'no_repair_exit', 'no_repair_closed', 'hash_wrong_without_repair', 'rank_without_repair', 'ref_rank', 'repair_repairs', 'repair_closed', 'repair_hash_match')}), flush=True)
    save(f'negative-withheld-k{a.k}-D{a.D}-gens{a.gens}' + a.tag_suffix, dict(extra=a.extra, k=a.k, D=a.D, gens=a.gens, rows=rows,
         note='withholding a low row is detected whenever the final basis is wrong; rows whose shifts already lie in the span of the others do not change the result and are reported as such',
         detected_all=all(not r['no_repair_closed'] and r['no_repair_exit'] == 3 for r in rows if r['hash_wrong_without_repair']),
         repaired_all=all(r['repair_hash_match'] and r['repair_closed'] for r in rows)))


def cmd_scan(a):
    """withhold each pending low row in turn (part 0 and 1, positions 0..n-1 of the pending list at the moment
    the part first has more than idx pending rows), randomized sweep off, no repair: every withheld row whose
    omission changes the final basis must be caught by the full sweep (exit 3)."""
    work = Path(a.work)
    path = inst_file(work, 'p15361', a.k, 15361, 4021, 9917, 12, True)
    ref = reference(path, a.D, a.gens)
    rows = []
    for part in (0, 1):
        for idx in range(a.n):
            out = work / 'scan.json'
            out.unlink(missing_ok=True)
            rc, rec, sec, _ = big(path, a.D, ['--withhold', f'{part}:{idx}', '--no-repair', '--t-random', '0', '--gens', a.gens,
                                                '--tile-rows', '40', '--batch', '16'], out)
            wrong = rec['basis_sha256'] != ref['basis_sha256']
            rows.append(dict(part=part, idx=idx, exit=rc, hash_wrong=wrong, closed=rec['closed'], rank=rec['rank'],
                             sweep_nonzero=rec['sweep_full']['nonzero_residuals']))
    ess = [r for r in rows if r['hash_wrong']]
    save(f'scan-withheld-k{a.k}-D{a.D}-gens{a.gens}', dict(k=a.k, D=a.D, gens=a.gens, tested=len(rows), essential=len(ess),
         essential_detected=sum(1 for r in ess if r['exit'] == 3 and not r['closed']),
         nonessential_clean=sum(1 for r in rows if not r['hash_wrong'] and r['exit'] == 0 and r['closed']),
         false_alarms=sum(1 for r in rows if not r['hash_wrong'] and r['exit'] != 0),
         undetected_wrong=sum(1 for r in ess if r['exit'] == 0), rows=rows))


def dev_for(work_gb, tile_rows, part_cols, tile_budget_mb):
    """--dev-gb that gives the requested tile budget (fake-gpu: no context slack)."""
    comp = tile_rows * part_cols * 4 / 2 ** 20
    return round((work_gb * 1024 + comp + tile_budget_mb) / 1024, 6)


def cmd_spill(a):
    """closures with a device budget small enough to force spilling (fake-gpu enforces the budget): hashes must
    equal the reference / saved reference records, and the spill counters must be nonzero."""
    work = Path(a.work)
    rows = []
    saved = {k: json.load(open(ROOT / 'experiments/ic-bounds-d5/verify' / f'k{k}-p12289-s10-D5.json')) for k in (9, 10)}
    cases = [('unit-12289', 12289, 1, 0, 10, False, False, 9, 5, 128, 6331, 55, 0.2),
             ('unit-12289', 12289, 1, 0, 10, False, False, 10, 5, 128, 10878, 150, 0.25),
             ('p15361', 15361, 4021, 9917, 12, True, False, 7, 6, 64, 3300, 13, 0.2),
             ('proper-65537', 65537, 4021, 9917, 12, True, False, 8, 5, 64, 3450, 14, 0.2),
             ('unit-65537', 65537, 1, 0, 10, False, True, 8, 5, 64, 3450, 14, 0.2)]
    for tag, prime, aa, bb, seed, want, nonres, k, D, tr, pc, tb, wg in cases:
        path = inst_file(work, tag, k, prime, aa, bb, seed, want, nonres)
        ref = None
        if tag == 'unit-12289' and k in saved:
            sv = saved[k]
            ref = dict(rank=sv['rank'], profile=sv['profile'], part_profiles=sv['part_profiles'], part_ranks=sv['part_ranks'],
                       contains_one=sv['contains_one'], basis_sha256=sv['basis_sha256'])
        else:
            ref = reference(path, D, 'all')
        dg = dev_for(wg, tr, pc, tb)
        out = work / 'spill.json'
        out.unlink(missing_ok=True)
        rc, rec, sec, r = big(path, D, ['--device', 'fake-gpu', '--tile-rows', str(tr), '--work-gb', str(wg), '--dev-gb', str(dg),
                                          '--batch', '256'], out)
        sp = rec['spill'] if rec else {}
        ok = bool(rc == 0 and same(ref, rec) and rec['closed'] and sp.get('spills', 0) > 0 and sp.get('fetches', 0) > 0)
        row = dict(case=f'{tag}-k{k}-D{D}', dev_gb=dg, work_gb=wg, tile_rows=tr, tile_budget_mb=tb, match=ok, rank=rec['rank'], closed=rec['closed'],
                   spill=sp, seconds=sec, dev_peak_live_mb=rec['dev_peak_live_mb'], dev_limit_mb=(sp.get('dev_limit_bytes') or 0) >> 20,
                   contains_one=rec['contains_one'], stderr_tail=r.stderr[-300:] if rc else '')
        print(json.dumps({k_: row[k_] for k_ in ('case', 'match', 'rank', 'seconds', 'dev_peak_live_mb', 'dev_limit_mb')}), flush=True)
        rows.append(row)
        save('spill', dict(rows=rows, all_match=all(r_['match'] for r_ in rows)))


def cmd_nosplit(a):
    work = Path(a.work)
    rows = []
    for tag, prime, aa, bb, seed, want in (('p15361', 15361, 4021, 9917, 12, True), ('unit-12289', 12289, 1, 0, 10, False)):
        for k, D in ((6, 3), (6, 4), (6, 5), (7, 4)):
            path = inst_file(work, tag, k, prime, aa, bb, seed, want)
            r = reference(path, D, 'all', split=False)
            row = dict(case=f'{tag}-k{k}-D{D}', ref_rank=r['rank'], modes={})
            for mode in ('plain', 'tiny-fake-gpu'):
                out = work / 'ns.json'
                out.unlink(missing_ok=True)
                rc, rec, sec, _ = big(path, D, MODES[mode] + ['--nosplit'], out)
                row['modes'][mode] = dict(match=bool(rc == 0 and same(r, rec) and rec['closed']), seconds=sec)
            rows.append(row)
            print(json.dumps(row), flush=True)
    save('nosplit', dict(rows=rows, all_match=all(m['match'] for r in rows for m in r['modes'].values())))


def cmd_rates(a):
    """closure + sweep at one thread on the unit p = 12289 instance for (D, k) pairs, with and without the
    reference engine's number (its madds counts n x rank x free: the dense method)."""
    work = Path(a.work)
    rows = []
    for dk in a.cases.split(','):
        D, k = (int(x) for x in dk.split(':'))
        path = inst_file(work, 'unit-12289', k, 12289, 1, 0, 10, False)
        out = work / f'rate-{D}-{k}.json'
        out.unlink(missing_ok=True)
        rc, rec, sec, _ = big(path, D, ['--work-gb', '0.4', '--batch', '256', '--log-seconds', '1000'], out)
        refname = ROOT / f'experiments/ic-bounds-d5/verify/k{k}-p12289-s10-D{D}.json'
        ref = json.load(open(refname)) if refname.exists() else None
        row = dict(D=D, k=k, columns=rec['columns'], rank=rec['rank'], contains_one=rec['contains_one'], seconds=rec['seconds'],
                   madds_closure_phase=rec['madds_closure_phase'], madds_total=rec['madds'], madds_dense=rec['madds_dense'],
                   madds_detail=rec['madds_detail'], seconds_detail=rec['seconds_detail'], rounds=rec['rounds'],
                   maxrss_mb=rec['maxrss_mb'], closed=rec['closed'], peak_entries=rec['peak_entries'],
                   ref_madds_dense=ref and ref['madds'], ref_seconds=ref and ref['seconds'],
                   ref_hash_match=ref and ref['basis_sha256'] == rec['basis_sha256'])
        d = rec['madds_detail']
        sd = rec['seconds_detail']
        row['rate_reduce_madds_per_s'] = d['reduce'] / sd['reduce'] if sd['reduce'] else None
        row['rate_update_madds_per_s'] = d['update'] / sd['update'] if sd['update'] else None
        row['rate_rref_madds_per_s'] = d['rref'] / sd['rref'] if sd['rref'] else None
        row['sparse_over_dense_reduce'] = d['reduce'] / d['reduce_dense'] if d['reduce_dense'] else None
        row['sparse_over_dense_total'] = row['madds_closure_phase'] / (d['reduce_dense'] + d['update'] + d['rref']) if d['reduce_dense'] else None
        print(json.dumps(row), flush=True)
        rows.append(row)
        save(f'rates-{a.tag}', dict(rows=rows, note='one thread (OPENBLAS_NUM_THREADS=1), engine defaults except --work-gb 0.4 --batch 256'))


def cmd_stopone(a):
    work = Path(a.work)
    rows = []
    for k, D in ((6, 5), (7, 5), (6, 6)):
        path = inst_file(work, 'unit-12289', k, 12289, 1, 0, 10, False)
        ref = reference(path, D, 'all')
        out = work / 's1.json'
        out.unlink(missing_ok=True)
        rc, rec, sec, _ = big(path, D, ['--stop-on-one', '--device', a.device], out)
        out2 = work / 's2.json'
        out2.unlink(missing_ok=True)
        rc2, rec2, sec2, _ = big(path, D, ['--device', a.device], out2)
        rows.append(dict(k=k, D=D, ref_contains_one=ref['contains_one'], exit_stop=rc, status=rec['status'],
                         contains_one=rec['contains_one'], closure_complete=rec['closure_complete'], closed=rec['closed'],
                         rank_at_stop=rec['rank'], full_rank=ref['rank'], full_run_contains_one=rec2['contains_one'],
                         full_run_matches=same(ref, rec2)))
    save('stop-on-one', dict(rows=rows, passed=all(r['contains_one'] and not r['closure_complete'] and r['exit_stop'] == 2
                                                    and r['full_run_matches'] for r in rows)))


def cmd_threads(a):
    work = Path(a.work)
    path = inst_file(work, 'p15361', 7, 15361, 4021, 9917, 12, True)
    ref = reference(path, 5, 'all')
    rows = []
    for th in (1, 2):
        out = work / 'th.json'
        out.unlink(missing_ok=True)
        env = {k: v for k, v in os.environ.items() if k not in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS')}
        cmd = [sys.executable, str(RUN), 'closure', '--instance', path, '--D', '5', '--out', str(out), '--threads', str(th)]
        subprocess.run(cmd, env=env, capture_output=True)
        rec = json.load(open(out))
        rows.append(dict(threads=th, blas_threads_env=rec['blas_threads_env'], match=same(ref, rec), seconds=rec['seconds']))
    save('threads', dict(rows=rows, all_match=all(r['match'] for r in rows), default_threads_note='len(os.sched_getaffinity(0)) = %d here' % len(os.sched_getaffinity(0))))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('cmd')
    ap.add_argument('--work', required=True)
    ap.add_argument('--k', type=int, default=7)
    ap.add_argument('--D', type=int, default=5)
    ap.add_argument('--device', default='cpu')
    ap.add_argument('--delays', default='2,2')
    ap.add_argument('--modes', default='plain')
    ap.add_argument('--gens', default='all')
    ap.add_argument('--cases', default='5:6')
    ap.add_argument('--n', type=int, default=40)
    ap.add_argument('--extra', default='')
    ap.add_argument('--tag', default='x')
    ap.add_argument('--tag-suffix', default='')
    a = ap.parse_args()
    globals()['cmd_' + a.cmd.replace('-', '_')](a)


if __name__ == '__main__':
    main()

