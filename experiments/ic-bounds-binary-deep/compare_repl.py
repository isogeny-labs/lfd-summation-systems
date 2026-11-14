#!/usr/bin/env python3
"""Compare a replication of the binary deep runs with the first run, record by record (REPLICATION.md).

For every system of the replication and every degree D of its ladder it compares, against the first run's
fetched/<name>.json.D<D>.json: columns, rank, codim, profile, contains_one, closed, unit, system_hash and the
canonical basis hash rref_sha256; and on the final record: exact_lfd, certificate_degree, strict_falls,
quotient_ranks, local_checks, n_points, points_sha256, status.  It also reports the settings that differ
(device, batch, tile_rows, col_block) so that an equal hash is known to come from a different execution.
usage: compare_repl.py --run NAME=DIR [--run NAME=DIR ...] [--tar FILE ...] [--first DIR] [--out FILE]
  --tar FILE   unpack a repl-<name>.tgz of the pod wrapper into fetched-repl/ first (members are checked)
"""
import argparse
import glob
import json
import os
import tarfile

HERE = os.path.dirname(os.path.abspath(__file__))
DEG_KEYS = ['columns', 'rank', 'codim', 'profile', 'contains_one', 'closed', 'unit', 'system_hash', 'rref_sha256']
FIN_KEYS = ['exact_lfd', 'certificate_degree', 'strict_falls', 'quotient_ranks', 'local_checks', 'n_points',
            'points_sha256', 'status', 'system_hash']
SET_KEYS = ['device', 'batch', 'tile_rows', 'col_block', 'tiles']


def load(path):
    with open(path) as f:
        return json.load(f)


def unpack(tar_path, dest):
    with tarfile.open(tar_path, 'r:gz') as tar:
        for m in tar.getmembers():
            name = m.name
            assert not name.startswith('/') and '..' not in name.split('/'), name
            assert m.isfile() or m.isdir(), name
        tar.extractall(dest, filter='data')


def compare(name, rdir, first):
    systems, problems = [], []
    finals = sorted(p for p in glob.glob(os.path.join(rdir, '*.json'))
                    if '.json.D' not in p and not p.endswith('.status.json')
                    and os.path.basename(p) not in ('chain-status.json',) and not os.path.basename(p).startswith('selftest'))
    if not finals:
        # an empty or wrong directory must not read as a pass
        problems.append('no system records in %s (point --run at the out/ directory of the replication)' % rdir)
    for path in finals:
        sysname = os.path.basename(path)[:-5]
        fpath = os.path.join(first, sysname + '.json')
        if not os.path.exists(fpath):
            problems.append('%s: no record in the first run' % sysname)
            continue
        a, b = load(fpath), load(path)
        row = dict(system=sysname, final_equal=True, degrees=0, degrees_equal=0, hashes_equal=0, proper_hashes_equal=0,
                   proper_degrees=0, differing=[], settings=[])
        for key in FIN_KEYS:
            if a.get(key) != b.get(key):
                row['final_equal'] = False
                row['differing'].append('final.' + key)
        cert = b.get('certificate_degree')
        if cert is None:
            row['differing'].append('no certificate degree in the replication')
            cert = -1
        for d in range(0, cert + 1):
            pa, pb = fpath + '.D%d.json' % d, path + '.D%d.json' % d
            if not (os.path.exists(pa) and os.path.exists(pb)):
                row['differing'].append('D%d: record missing (%s)' % (d, 'first' if not os.path.exists(pa) else 'replication'))
                row['degrees'] += 1
                continue
            da, db = load(pa), load(pb)
            row['degrees'] += 1
            bad = [k for k in DEG_KEYS if da.get(k) != db.get(k)]
            if not bad:
                row['degrees_equal'] += 1
            same = da.get('rref_sha256') == db.get('rref_sha256') and da.get('rref_sha256') is not None
            row['hashes_equal'] += int(same)
            # a basis of the whole space (or of zero) has a hash that says nothing: count the proper subspaces apart
            proper = isinstance(db.get('rank'), int) and isinstance(db.get('columns'), int) and 0 < db['rank'] < db['columns']
            row['proper_degrees'] += int(proper)
            row['proper_hashes_equal'] += int(proper and same)
            row['differing'] += ['D%d.%s' % (d, k) for k in bad]
            if d == cert:
                row['settings'] = [dict(key=k, first=da.get(k), replication=db.get(k)) for k in SET_KEYS if da.get(k) != db.get(k)]
                row['rref_sha256_at_certificate'] = db.get('rref_sha256')
                row['columns_at_certificate'] = db.get('columns')
                row['seconds'] = dict(first=a.get('seconds'), replication=b.get('seconds'))
        row['equal'] = row['final_equal'] and row['degrees'] == row['degrees_equal'] and not row['differing']
        if not row['equal']:
            problems.append('%s: %s' % (sysname, ', '.join(row['differing'])))
        systems.append(row)
    return dict(run=name, directory=os.path.relpath(rdir, HERE), systems=len(systems),
                equal=sum(r['equal'] for r in systems), degrees=sum(r['degrees'] for r in systems),
                degrees_equal=sum(r['degrees_equal'] for r in systems), hashes_equal=sum(r['hashes_equal'] for r in systems),
                proper_degrees=sum(r['proper_degrees'] for r in systems),
                proper_hashes_equal=sum(r['proper_hashes_equal'] for r in systems),
                with_different_settings=sum(1 for r in systems if r['settings']), problems=problems, rows=systems)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', action='append', default=[], help='NAME=DIR of a replication out directory')
    ap.add_argument('--tar', action='append', default=[])
    ap.add_argument('--first', default=os.path.join(HERE, 'fetched'))
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    for t in a.tar:
        unpack(t, os.path.join(HERE, 'fetched-repl'))
    results = []
    for spec in a.run:
        name, rdir = spec.split('=', 1)
        results.append(compare(name, rdir, a.first))
    for r in results:
        print('%s: %d systems, %d equal; %d degree records, %d equal, %d basis hashes equal (%d of %d on proper subspaces); %d systems run with different settings; %d problems'
              % (r['run'], r['systems'], r['equal'], r['degrees'], r['degrees_equal'], r['hashes_equal'], r['proper_hashes_equal'],
                 r['proper_degrees'], r['with_different_settings'], len(r['problems'])))
        for p in r['problems']:
            print('  PROBLEM', p)
    if a.out:
        json.dump(dict(first=os.path.relpath(a.first, HERE), runs=results), open(a.out, 'w'), indent=1)
    if not results or any(r['problems'] for r in results):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
