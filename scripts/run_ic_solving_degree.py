"""Solving degree of the paper's exact-LFD systems (experiments/ic-solving-degree/PROTOCOL.md).

sd_sigma = max(LFD, max.GB.deg_sigma) (Caminata--Gorla Thm 3.1).  Each system is tagged unit / count / cube
(docs/proofs/N10, Cor. 4) and, when it is proper with a solution set of at most --max-points points, its exact
max.GB.deg is computed for DRL and grlex by Buchberger--Moeller on the squarefree staircase
(csd/ic_solving_degree.py).
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from csd import ic_solving_degree as SD  # noqa: E402

EXP = ROOT / 'experiments'


def jl(path):
    return [json.loads(line) for line in open(path) if line.strip()]


def boolean_zero_set(equations, r):
    """Zero set of Boolean polynomials (lists of squarefree monomial masks) on {0,1}^r by a subset (zeta) transform."""
    assert r <= 24
    T = np.zeros(1 << r, dtype=np.uint64)
    words = []
    for j, eq in enumerate(equations):
        w, bit = divmod(j, 64)
        while len(words) <= w:
            words.append(np.zeros(1 << r, dtype=np.uint64))
        for mask in eq:
            words[w][int(mask)] ^= np.uint64(1) << np.uint64(bit)
    alive = np.ones(1 << r, dtype=bool)
    for W in words:
        for i in range(r):
            V = W.reshape(-1, 2, 1 << i)
            V[:, 1, :] ^= V[:, 0, :]
        alive &= (W == 0)
    del T
    return [int(z) for z in np.flatnonzero(alive)]


# ---------------------------------------------------------------------------- adapters
# each yields dict(set, id, family, m, k, r, p, lfd, npoints, points=PointsObject|None, extra)

def sweep_like(set_name, directory):
    systems = jl(EXP / directory / 'systems.jsonl')
    results = {r['system_hash']: r for r in jl(EXP / directory / 'results.jsonl')}
    for s in systems:
        res = results[s['system_hash']]
        lfd = res.get('exact_lfd')
        if lfd is None:
            continue
        sol = s['coordinate_solutions']
        r = s['variables']
        if s['kind'] == 'prime':
            k = s['k']
            exps = [((z & ((1 << k) - 1)), z >> k) for z in sol]
            pts = SD.PrimeChainPoints(exps, 2, k, s['generator_value'], s['prime']) if sol else None
            yield dict(set=set_name, id=s['system_hash'][:16], family='prime', m=2, k=k, r=r, p=s['prime'],
                       lfd=lfd, npoints=len(sol), points=pts, system=s)
        else:
            pts = SD.BooleanPoints(sol, r) if sol else None
            yield dict(set=set_name, id=s['system_hash'][:16], family='binary', m=2, k=None, r=r, p=2,
                       lfd=lfd, npoints=len(sol), points=pts, system=s)


def msummands():
    for f in sorted(glob.glob(str(EXP / 'ic-bounds-msummands/records/*.json'))):
        rec = json.load(open(f))
        if rec.get('exact_lfd') in (None, 'None'):
            continue
        m, k = int(rec['m']), int(rec['k'])
        n = rec.get('solutions')
        n = int(n) if n not in (None, 'None') else None
        yield dict(set='S3', id=rec['case'], family=rec['family'], m=m, k=k, r=m * k, p=int(rec['p']) if 'p' in rec else 2,
                   lfd=int(rec['exact_lfd']), npoints=n, points=None, record=f)
    for f in sorted(glob.glob(str(EXP / 'ic-bounds-msummands/pod-records/*/out/*.json'))):
        if '.D' in Path(f).name or 'status' in Path(f).name:
            continue
        rec = json.load(open(f))
        st = str(rec.get('status', ''))
        if not st.startswith('exact') or rec.get('exact_lfd') in (None, 'None'):
            continue
        import re
        mm = re.search(r'm(\d+)(?:n\d+)?k(\d+)', Path(f).name)
        m, k = int(mm.group(1)), int(mm.group(2))
        n = rec.get('solutions', rec.get('n_points'))
        n = int(n) if n not in (None, 'None') else None
        fam = rec.get('family', 'binary' if Path(f).name.startswith('G-') else 'prime')
        yield dict(set='S3pod', id=rec.get('case', Path(f).stem), family=fam, m=m, k=k, r=m * k,
                   p=int(rec['p']) if 'p' in rec else 2, lfd=int(rec['exact_lfd']), npoints=n, points=None, record=f)


def msum():
    for rec in jl(EXP / 'ic-bounds-msum/results.jsonl'):
        if rec.get('lfd_exact') is None:
            continue
        m, k = rec['m'], rec['k']
        yield dict(set='S4', id='m%dk%dp%da%db%dt%d' % (m, k, rec['p'], rec['a'], rec['b'], rec['t']), family='prime', m=m,
                   k=k, r=m * k, p=rec['p'], lfd=rec['lfd_exact'], npoints=rec['zeros'], points=None)


def binary_results():
    """S5: systems rebuilt from the stored parameters (the records keep only the equation count); the rebuilt
    system hash must equal the record's, and Z is the builder's field-level enumeration."""
    sys.path.insert(0, str(EXP / 'ic-curve-validation' / 'code'))
    from csd import ic_dreg as I
    from csd import ic_bounds_binary as B
    from csd.ic_research_campaign import field
    for name in ['results-m2.jsonl', 'results-m3-direct.jsonl', 'results-m3-split.jsonl']:
        for rec in jl(EXP / 'ic-bounds-binary' / name):
            if rec.get('exact_lfd') is None:
                continue
            r = rec['variables']
            m = 2 if 'm2' in name else 3
            row = dict(set='S5', id='%s:%s' % (name[8:-6], rec['case_id']), family='binary', m=m, k=rec.get('k'), r=r,
                       p=2, lfd=rec['exact_lfd'], npoints=rec.get('solutions'), points=None)
            if rec.get('solutions') and rec['solutions'] <= MAX_POINTS[0] and not TAGS_ONLY[0]:
                F = field(rec['n'])
                assert F.modulus == rec['modulus']
                E = I.BinaryCurve(F, rec['a'], rec['b'])
                target = tuple(rec['target'])
                planted = [tuple(P) for P in rec['planted']] if rec.get('planted') else None
                if rec['family'] == 'm2':
                    system = B.m2_system(E, rec['basis'], rec['base_kind'], target, rec['target_class'], planted)
                elif rec['family'] == 'm3-direct':
                    system = B.m3_direct_system(E, rec['basis'], target)
                else:
                    system = B.m3_split_system(E, rec['basis'], target)
                assert system['system_hash'] == rec['system_hash'], rec['case_id']
                assert len(system['points']) == rec['solutions']
                row['equations'] = system['equations']
                row['field_points'] = system['points']
            yield row


def binary_deep():
    outcomes = json.load(open(EXP / 'ic-bounds-binary-deep/outcomes-2026-11-12.json'))
    lfd = {row['name']: row for row in outcomes['rows']}
    for f in sorted(glob.glob(str(EXP / 'ic-bounds-binary-deep/instances/*.json'))):
        inst = json.load(open(f))
        row = lfd.get(inst['name'])
        if row is None or row.get('lfd') is None or not str(row.get('status', '')).startswith('exact'):
            continue
        r = inst['variables']
        pts = SD.BooleanPoints(inst['points'], r) if inst['points'] else None
        yield dict(set='S6', id=inst['name'], family='binary', m=2, k=inst.get('k'), r=r, p=2, lfd=row['lfd'],
                   npoints=len(inst['points']), points=pts)


# exact LFD of the depth 10-16 instances, from experiments/ic-bounds-k14/RESULTS.md (None: only a bound there)
DEEP_LFD = {'proper-p114689-s22-k13': 6, 'proper-p114689-s22-k14': 7, 'proper-p65537-s12-k13': 6,
            'unit-p65537-a1-s10-k10': 6, 'unit-p65537-a1-s10-k11': 6, 'unit-p65537-a1-s10-k12': 6,
            'unit-p65537-a1-s10-k13': 6, 'unit-p65537-a1-s10-k14': 7, 'unit-p163841-a3-s21-k13': 6,
            'unit-p163841-a3-s21-k14': 7}


def deep_prime():
    from csd.ic_curve_sweep import subgroup_generator_value
    for f in sorted(glob.glob(str(EXP / 'ic-bounds-k14/instances/*.json'))):
        inst = json.load(open(f))
        name = Path(f).stem
        Z = inst['Z'] if isinstance(inst.get('Z'), list) else []
        p, k = inst['prime'], inst['k']
        pts = None
        if Z and len(Z) <= MAX_POINTS[0] and not TAGS_ONLY[0]:
            pts = SD.PrimeChainPoints([tuple(z) for z in Z], 2, k, subgroup_generator_value(p, k), p)
        yield dict(set='S7', id=name, family='prime', m=2, k=k, r=2 * k, p=p, lfd=DEEP_LFD.get(name),
                   npoints=inst['points'] if isinstance(inst.get('points'), int) else len(Z), points=pts)


MAX_POINTS = [3000]
TAGS_ONLY = [False]

SOURCES = dict(S1=lambda: sweep_like('S1', 'ic-curve-sweep'),
               S2=lambda: (x for d in ('ic-chain-depth8', 'ic-chain-depth9') for x in sweep_like('S2', d)),
               S3=msummands, S4=msum, S5=binary_results, S6=binary_deep, S7=deep_prime)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--sets', default='S1,S2,S3,S4,S5,S6,S7')
    ap.add_argument('--max-points', type=int, default=3000)
    ap.add_argument('--tags-only', action='store_true')
    ap.add_argument('--out', default=str(EXP / 'ic-solving-degree' / 'rows.jsonl'))
    a = ap.parse_args()
    MAX_POINTS[0], TAGS_ONLY[0] = a.max_points, a.tags_only
    rows = []
    t0 = time.time()
    for set_name in a.sets.split(','):
        for inst in SOURCES[set_name]():
            lfd, n, r = inst['lfd'], inst['npoints'], inst['r']
            tag, sd = SD.resolve(lfd, n, r) if lfd is not None else (None, None)
            row = dict(set=inst['set'], id=inst['id'], family=inst['family'], m=inst['m'], k=inst['k'], r=r, p=inst['p'],
                       lfd=lfd, points=n, tag=tag, count_bound=(max(2, (n + 1).bit_length() - 1) if n else None))
            pts = inst.get('points')
            if inst.get('field_points') is not None:
                zs = boolean_zero_set(inst['equations'], r)
                row['zeta_equals_field_enumeration'] = zs == sorted(inst['field_points'])
                pts = SD.BooleanPoints(inst['field_points'], r)
            elif pts is None and inst.get('equations') is not None and n and n <= a.max_points and not a.tags_only:
                zs = boolean_zero_set(inst['equations'], r)
                assert len(zs) == n, (inst['id'], len(zs), n)
                if inst.get('points_sha256'):
                    row['points_sha256_ok'] = hashlib.sha256(json.dumps(zs).encode()).hexdigest() == inst['points_sha256']
                pts = SD.BooleanPoints(zs, r)
            if pts is not None and n and n <= a.max_points and not a.tags_only:
                for order in ('drl', 'grlex'):
                    st = SD.staircase(pts, n, r, inst['p'], order)
                    row['gb_' + order] = st['max_gb_degree']
                    row['hf_' + order] = st['standard_by_degree']
                    row['gens_' + order] = st['generators_by_degree']
                assert row['hf_drl'] == row['hf_grlex'], 'affine Hilbert function must not depend on the order'
            if lfd is not None:
                gb = row.get('gb_drl')
                if tag is not None:
                    row['sd_drl'] = lfd
                    if gb is not None and gb > lfd:
                        row['tag_violation'] = True       # would refute N10 Cor. 4; reported, never dropped
                elif gb is not None:
                    row['sd_drl'] = max(lfd, gb)
                    row['sd_grlex'] = max(lfd, row['gb_grlex'])
                    row['tag'] = 'computed'
            rows.append(row)
        print(json.dumps(dict(set=set_name, rows=sum(1 for x in rows if x['set'].startswith(set_name)),
                              seconds=round(time.time() - t0, 1))), flush=True)
    with open(a.out, 'w') as f:
        for row in rows:
            f.write(json.dumps(row, sort_keys=True) + '\n')
    print('wrote', a.out, len(rows))


if __name__ == '__main__':
    main()
