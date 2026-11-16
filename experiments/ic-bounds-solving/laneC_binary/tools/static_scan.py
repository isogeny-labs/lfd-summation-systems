#!/usr/bin/env python3
"""Static Boolean Macaulay scan (lane C, solving-degree study 20261115).

For a binary m=2 system (n Boolean equations of degree 2 in r=2k variables) build the bit-packed
Macaulay matrix Mac_D (every equation times every squarefree monomial of degree <= D-2, field
equations built in: csd.ic_dreg.boolean_macaulay of the repository), echelonise it over GF(2)
IN PLACE (tools/gf2rank.c, OpenMP), and record per degree D:

  rank, deficiency = cols - |Z| - rank   (dim of (I cap M_<=D) / Mac_D; 0 iff row space = I cap M_<=D),
  ev   = rank of the evaluation of M_<D on Z (needs ev = |Z|),
  one  = 1 in the row space (unit ideals), lin = affine forms in the row space.
d_full = least D with deficiency 0 and ev = |Z|  (proper ideal)   /  least D with 1 in row space (unit),
exactly the definition of scripts/ic_bounds_binary_related.py.  Z is recomputed here by brute force
over all 2^r assignments (independent of the repository's point enumeration) and compared with the stored |Z|.

Sources of instances:
  results-m2.jsonl   (k = 8..10, rebuilt with scripts/ic_bounds_binary_validate.rebuild, system hash checked)
  ic-curve-sweep systems.jsonl (k = 2..7, equations stored, found by system hash from ic-bounds-binary-law/table.jsonl)
"""
import argparse, ctypes, json, os, sys, time, resource, hashlib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = '/home/user/ml-cryptanalysis'
sys.path.insert(0, REPO); sys.path.insert(0, os.path.join(REPO, 'scripts'))
from csd import ic_dreg as I  # noqa: E402

SO = os.environ.get('GF2RANK_SO', os.path.join(HERE, 'gf2rank.so'))


def lib():
    if not os.path.exists(SO):
        os.system(f'gcc -O3 -march=native -fopenmp -shared -fPIC -o {SO} {os.path.join(HERE, "gf2rank.c")}')
    L = ctypes.CDLL(SO)
    L.gf2_rank_inplace.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_long, ctypes.c_long, ctypes.c_void_p, ctypes.c_int]
    L.gf2_rank_inplace.restype = ctypes.c_long
    return L


def echelon(M, cols, threads):
    piv = np.zeros(min(M.shape[0], cols), dtype=np.int64)
    r = lib().gf2_rank_inplace(M.ctypes.data, M.shape[0], M.shape[1], cols, piv.ctypes.data, threads)
    return r, piv[:r]


def brute_points(eqs, r):
    x = np.arange(1 << r, dtype=np.uint32)
    ok = np.ones(1 << r, dtype=bool)
    for e in eqs:
        v = np.zeros(1 << r, dtype=np.uint8)
        for m in e:
            m = np.uint32(m)
            v ^= ((x & m) == m).astype(np.uint8)
        ok &= (v == 0)
    return [int(p) for p in x[ok]]


def eval_rank(points, r, D, threads):
    """rank of evaluation of squarefree monomials of degree < D on the points."""
    masks = I._masks_upto(r, D - 1)
    pts = np.asarray(points, dtype=np.uint64)
    E = ((pts[:, None] & masks[None, :]) == masks[None, :])
    P = np.packbits(E, axis=1, bitorder='little')
    pad = (-P.shape[1]) % 8
    if pad:
        P = np.pad(P, ((0, 0), (0, pad)))
    W = np.ascontiguousarray(P).view(np.uint64)
    rk, _ = echelon(W, len(masks), threads)
    return rk


def scan(eqs, r, points, max_cols, max_gb, max_degree, threads, time_cap):
    Z = len(points)
    d_in = max(m.bit_count() for e in eqs for m in e)
    out = dict(d_full=None, d_full_lb=None, stop_reason=None, per_degree=[], d_lin=None)
    t_start = time.time()
    for D in range(d_in, max_degree + 1):
        cols = I.n_cols(r, D)
        nrows = sum(I.n_cols(r, D - max(m.bit_count() for m in e)) for e in eqs)
        gb = nrows * ((cols + 63) // 64) * 8 / 1e9
        if cols > max_cols or gb > max_gb:
            out['d_full_lb'] = D; out['stop_reason'] = f'cap at D={D}: cols={cols} rows={nrows} matrix={gb:.2f}GB'
            return out
        if time.time() - t_start > time_cap:
            out['d_full_lb'] = D; out['stop_reason'] = f'time cap before D={D}'
            return out
        t0 = time.time()
        M, cols = I.boolean_macaulay(eqs, r, D)
        rank, piv = echelon(M, cols, threads)
        del M
        piv_set = set(int(c) for c in piv)
        lin = sum(1 for c in piv_set if c >= cols - (r + 1))
        one = (cols - 1) in piv_set
        rec = dict(D=D, rows=nrows, columns=cols, rank=rank, lin=lin, one=one,
                   deficiency=(cols - Z - rank) if Z else None, seconds=round(time.time() - t0, 1))
        if Z == 0:
            full = one
        else:
            full = False
            if rank == cols - Z:
                ev = eval_rank(points, r, D, threads)
                rec['ev'] = ev
                full = (ev == Z)
        out['per_degree'].append(rec)
        if full:
            out['d_full'] = D
            return out
    out['d_full_lb'] = max_degree + 1; out['stop_reason'] = f'max_degree {max_degree}'
    return out


def load_instances(ks):
    T = {}
    for l in open(os.path.join(REPO, 'experiments/ic-bounds-binary-law/table.jsonl')):
        t = json.loads(l)
        if t['source'] == 'ic-binary-inverse' or t['k'] not in ks:
            continue
        T.setdefault(t['system_hash'], t)
    sweep = {}
    need = {h for h, t in T.items() if t['source'] == 'ic-curve-sweep'}
    for l in open(os.path.join(REPO, 'experiments/ic-curve-sweep/systems.jsonl')):
        x = json.loads(l)
        if x['system_hash'] in need and x['kind'] == 'binary':
            sweep[x['system_hash']] = x['equations']
    res = {}
    for l in open(os.path.join(REPO, 'experiments/ic-bounds-binary/results-m2.jsonl')):
        x = json.loads(l)
        res.setdefault(x['system_hash'], x)
    inst = []
    for h, t in T.items():
        if t['source'] == 'ic-curve-sweep':
            inst.append(dict(system_hash=h, k=t['k'], n=t['n'], lfd=t['lfd'], cert=t['cert_degree'], Z=t['Z'], eqs=sweep[h], src='sweep'))
        else:
            inst.append(dict(system_hash=h, k=t['k'], n=t['n'], lfd=t['lfd'], cert=t['cert_degree'], Z=t['Z'], row=res[h], src='bounds'))
    return inst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--k', type=int, nargs='+', required=True)
    ap.add_argument('--per-group', type=int, default=2, help='systems per (k, n, LFD, unit?) group, by hash order')
    ap.add_argument('--max-cols', type=int, default=70000)
    ap.add_argument('--max-gb', type=float, default=2.8)
    ap.add_argument('--max-degree', type=int, default=12)
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--time-cap', type=float, default=3600, help='seconds of scanning per system before the next degree is refused')
    ap.add_argument('--hashes', nargs='*', help='only these system hashes')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    done = set()
    if os.path.exists(a.out):
        done = {json.loads(l)['system_hash'] for l in open(a.out) if l.strip()}
    inst = load_instances(set(a.k))
    groups = {}
    for x in sorted(inst, key=lambda x: x['system_hash']):
        groups.setdefault((x['k'], x['n'], x['lfd'], x['Z'] == 0), []).append(x)
    sel = [x for g in groups.values() for x in g[:a.per_group]]
    if a.hashes:
        sel = [x for x in inst if x['system_hash'] in set(a.hashes)]
    sel.sort(key=lambda x: (x['k'], x['n'], x['lfd'], x['system_hash']))
    sel = [x for x in sel if x['system_hash'] not in done]
    print(len(sel), 'systems', flush=True)
    for i, x in enumerate(sel):
        if x['src'] == 'bounds':
            from ic_bounds_binary_validate import rebuild
            sysm = rebuild(x['row']); assert sysm['system_hash'] == x['system_hash']
            eqs, r = sysm['equations'], sysm['variables']
        else:
            eqs, r = x['eqs'], 2 * x['k']
            canon = dict(kind='binary', variables=r, equations=eqs)
        points = brute_points(eqs, r)
        assert len(points) == x['Z'], (x['system_hash'], len(points), x['Z'])
        t0 = time.time()
        res = scan(eqs, r, points, a.max_cols, a.max_gb, a.max_degree, a.threads, a.time_cap)
        rec = dict(system_hash=x['system_hash'], k=x['k'], n=x['n'], r=r, Z=len(points), lfd=x['lfd'], cert_degree=x['cert'],
                   src=x['src'], seconds=round(time.time() - t0, 1), peak_rss_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024, **res)
        if rec['d_full'] is not None and rec['cert_degree'] is not None:
            assert rec['cert_degree'] <= rec['d_full'], rec
        open(a.out, 'a').write(json.dumps(rec, sort_keys=True) + '\n')
        print(f"[{i+1}/{len(sel)}] k={x['k']} n={x['n']} lfd={x['lfd']} D*={x['cert']} Z={len(points)} d_full={rec['d_full']} lb={rec['d_full_lb']} {rec['seconds']}s", flush=True)


if __name__ == '__main__':
    main()
