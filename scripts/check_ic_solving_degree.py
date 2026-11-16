"""Independent checks C1-C4 of experiments/ic-solving-degree/PROTOCOL.md.

  export   write the Sage input (C1/C2 sample) to experiments/ic-solving-degree/sage-in.json
  compare  read the Sage output; C1 (GB vs Buchberger-Moeller), C2 (sd from an independent dense
           mutant-iterated Macaulay closure, CG Def. 1.1), C3 (brute-force Z), C4 (|Z| <= 2^{k+1});
           writes experiments/ic-solving-degree/checks.json
"""
from __future__ import annotations

import argparse
import itertools
import json
import random
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from csd import ic_solving_degree as SD  # noqa: E402
from csd.elim import fp_echelon  # noqa: E402

EXP = ROOT / 'experiments'
OUT = EXP / 'ic-solving-degree'
SEED = 20261114


def jl(path):
    return [json.loads(line) for line in open(path) if line.strip()]


def load_sweep_like(directory):
    systems = jl(EXP / directory / 'systems.jsonl')
    results = {r['system_hash']: r for r in jl(EXP / directory / 'results.jsonl')}
    return [(s, results[s['system_hash']]) for s in systems]


def mask_exp(mask, r):
    return [(mask >> v) & 1 for v in range(r)]


def input_polys(s):
    """Generators + chain equations (prime) or Boolean equations (binary) as [[exps, coef], ...]."""
    r = s['variables']
    polys = []
    if s['kind'] == 'prime':
        p, k = s['prime'], s['k']
        for g in s['generators']:
            polys.append([[mask_exp(m, r), c % p] for m, c in g])
        for i in range(r // k):
            for j in range(k):
                v = i * k + j
                sq = [0] * r
                sq[v] = 2
                if j < k - 1:
                    polys.append([[sq, 1], [mask_exp(1 << (v + 1), r), p - 1]])
                else:
                    polys.append([[sq, 1], [[0] * r, p - 1]])
    else:
        for eq in s['equations']:
            polys.append([[mask_exp(m, r), 1] for m in eq])
        for v in range(r):
            sq = [0] * r
            sq[v] = 2
            polys.append([[sq, 1], [mask_exp(1 << v, r), 1]])
    return polys


def sample():
    rng = random.Random(SEED)
    s1 = load_sweep_like('ic-curve-sweep')
    proper = [x for x in s1 if x[0]['coordinate_solutions']]
    c1 = rng.sample(proper, 150)
    small = [x for x in s1 if x[0]['variables'] <= 8]
    c2 = rng.sample([x for x in small if x[0]['kind'] == 'prime'], 20) + \
        rng.sample([x for x in small if x[0]['kind'] != 'prime'], 20)
    s2 = load_sweep_like('ic-chain-depth8') + load_sweep_like('ic-chain-depth9')
    return c1, c2, s2


def extras():
    """Systems outside the seeded sample that C1 must also cover: every system with max.GB.deg > LFD
    (rows.jsonl) and two deep instances.  Returns [(id, p, r, polys, points_object, npoints)]."""
    rows = jl(OUT / 'rows.jsonl')
    over = {r['id'] for r in rows if r.get('gb_drl') is not None and r['lfd'] is not None and r['gb_drl'] > r['lfd']}
    out = []
    for s, res in load_sweep_like('ic-curve-sweep'):
        if s['system_hash'][:16] in over:
            out.append((s['system_hash'], s['prime'] if s['kind'] == 'prime' else 2, s['variables'], input_polys(s),
                        points_of(s), len(s['coordinate_solutions'])))
    sys.path.insert(0, str(ROOT / 'scripts'))
    import run_ic_solving_degree as RUN
    for row in RUN.binary_results():
        if row['id'] in over:
            fake = dict(kind='binary', variables=row['r'], equations=row['equations'])
            out.append((row['id'], 2, row['r'], input_polys(fake), SD.BooleanPoints(row['field_points'], row['r']),
                        len(row['field_points'])))
    from csd.ic_curve_sweep import subgroup_generator_value
    for name in ['proper-p114689-s22-k13', 'proper-p65537-s12-k13']:
        inst = json.load(open(EXP / 'ic-bounds-k14/instances' / (name + '.json')))
        p, k = inst['prime'], inst['k']
        gen = []
        for key, c in inst['generator'].items():
            a1, a2 = map(int, key.split(','))
            gen.append([a1 | (a2 << k), c])
        fake = dict(kind='prime', variables=2 * k, prime=p, k=k, generators=[gen])
        out.append((name, p, 2 * k, input_polys(fake),
                    SD.PrimeChainPoints([tuple(z) for z in inst['Z']], 2, k, subgroup_generator_value(p, k), p),
                    len(inst['Z'])))
    return out


def export(_a):
    c1, c2, s2 = sample()
    seen, items = set(), []
    for eid, p, r, polys, _, _ in extras():
        seen.add(eid)
        items.append(dict(id=eid, p=p, r=r, polys=polys))
    for s, res in c1 + c2 + s2:
        if s['system_hash'] in seen:
            continue
        seen.add(s['system_hash'])
        items.append(dict(id=s['system_hash'], p=s['prime'] if s['kind'] == 'prime' else 2, r=s['variables'],
                          polys=input_polys(s)))
    (OUT / 'sage-in.json').write_text(json.dumps(items))
    print('exported', len(items))


def points_of(s):
    sol = s['coordinate_solutions']
    r = s['variables']
    if s['kind'] == 'prime':
        k = s['k']
        return SD.PrimeChainPoints([(z & ((1 << k) - 1), z >> k) for z in sol], 2, k, s['generator_value'], s['prime'])
    return SD.BooleanPoints(sol, r)


# ------------------------------------------------------------------ C2: dense closure, CG Definition 1.1
def monomials(r, d):
    # columns in decreasing degree, so that the echelon rows with pivot in degree < d span the
    # elements of the row space of degree < d (the rows the mutant step must multiply)
    out = []
    for deg in range(d, -1, -1):
        for c in itertools.combinations_with_replacement(range(r), deg):
            e = [0] * r
            for v in c:
                e[v] += 1
            out.append(tuple(e))
    return out


class Closure:
    """V_{F,d}: row space of the degree-d Macaulay matrix after the mutant iteration (CG Section 1)."""

    def __init__(self, polys, r, p, d):
        self.r, self.p, self.d = r, p, d
        self.mons = monomials(r, d)
        self.col = {m: i for i, m in enumerate(self.mons)}
        rows = []
        for f in polys:
            fd = max(sum(e) for e, _ in f)
            for mu in self.mons:
                if sum(mu) + fd <= d:
                    rows.append(self.vec([(tuple(a + b for a, b in zip(e, mu)), c) for e, c in f]))
        self.basis = self.echelon(np.array(rows, dtype=np.int64)) if rows else np.zeros((0, len(self.mons)), np.int64)
        while True:
            new = []
            for row in self.basis:
                nz = np.flatnonzero(row)
                deg = max(sum(self.mons[i]) for i in nz)
                if deg >= d:
                    continue
                for v in range(r):
                    shifted = np.zeros(len(self.mons), dtype=np.int64)
                    for i in nz:
                        e = list(self.mons[i])
                        e[v] += 1
                        shifted[self.col[tuple(e)]] = row[i]
                    new.append(shifted)
            if not new:
                break
            before = len(self.basis)
            self.basis = self.echelon(np.vstack([self.basis, np.array(new)]))
            if len(self.basis) == before:
                break

    def vec(self, terms):
        v = np.zeros(len(self.mons), dtype=np.int64)
        for e, c in terms:
            v[self.col[tuple(e)]] = (v[self.col[tuple(e)]] + c) % self.p
        return v

    def echelon(self, M):
        out = []
        rank, _, _ = fp_echelon(M % self.p, self.p, out)
        return out[0][:rank] if rank else np.zeros((0, M.shape[1]), np.int64)

    def contains(self, terms):
        if max(sum(e) for e, _ in terms) > self.d:
            return False
        v = self.vec(terms)
        rank, _, _ = fp_echelon(np.vstack([self.basis, v[None, :]]) % self.p, self.p)
        return rank == len(self.basis)


def solving_degree_direct(polys, basis, r, p, dmax=9):
    for d in range(0, dmax + 1):
        C = Closure(polys, r, p, d)
        if all(C.contains([(tuple(e), c) for e, c in g]) for g in basis):
            return d
    return None


# ------------------------------------------------------------------ compare
def compare(a):
    c1, c2, s2 = sample()
    sage = {}
    if (OUT / 'sage-out.json').exists():
        sage.update({x['id']: x for x in json.load(open(OUT / 'sage-out.json'))})
    for f in [OUT / 'sage-out.jsonl', *sorted(OUT.glob('sage-shard*-out.jsonl')), *sorted(OUT.glob('sage-bin*-out.jsonl'))]:
        if not f.exists():
            continue
        for x in jl(f):
            if not x.get('timeout'):
                sage[x['id']] = x
    report = dict(C1=[], C2=[], C3=[], C4={}, not_run=[])
    # C1
    for s, res in c1 + s2:
        g = sage.get(s['system_hash'])
        if g is None:
            report['not_run'].append(dict(check='C1', id=s['system_hash'][:16], r=s['variables']))
            continue
        r = s['variables']
        p = s['prime'] if s['kind'] == 'prime' else 2
        n = len(s['coordinate_solutions'])
        st = SD.staircase(points_of(s), n, r, p, 'drl')
        lead = sorted([mask_exp(m, r) for m in st['squarefree_generators']] +
                      [[2 if w == v else 0 for w in range(r)] for v in
                       (x.bit_length() - 1 for x in st['standard_variables'])])
        report['C1'].append(dict(id=s['system_hash'][:16], points=n, dimension=g['dimension'],
                                 gb_bm=st['max_gb_degree'], gb_sage=g['max_degree'],
                                 leading_equal=lead == g['leading'],
                                 ok=(g['dimension'] == n and st['max_gb_degree'] == g['max_degree'] and lead == g['leading'])))
    # C1 extras: every system with max.GB.deg > LFD, and two deep instances
    for eid, p, r, polys, pts, n in extras():
        g = sage.get(eid)
        if g is None:
            report['not_run'].append(dict(check='C1', id=str(eid)[:24], r=r))
            continue
        st = SD.staircase(pts, n, r, p, 'drl')
        lead = sorted([mask_exp(m, r) for m in st['squarefree_generators']] +
                      [[2 if w == v else 0 for w in range(r)] for v in
                       (x.bit_length() - 1 for x in st['standard_variables'])])
        report['C1'].append(dict(id=str(eid)[:24], points=n, dimension=g['dimension'], gb_bm=st['max_gb_degree'],
                                 gb_sage=g['max_degree'], leading_equal=lead == g['leading'], extra=True,
                                 ok=(g['dimension'] == n and st['max_gb_degree'] == g['max_degree'] and lead == g['leading'])))
    # C2
    for s, res in c2:
        g = sage.get(s['system_hash'])
        if g is None or 'basis' not in g:
            report['not_run'].append(dict(check='C2', id=s['system_hash'][:16], r=s['variables']))
            continue
        r = s['variables']
        p = s['prime'] if s['kind'] == 'prime' else 2
        sd = solving_degree_direct(input_polys(s), g['basis'], r, p)
        predicted = max(res['exact_lfd'], g['max_degree'])
        report['C2'].append(dict(id=s['system_hash'][:16], kind=s['kind'], r=r, lfd=res['exact_lfd'],
                                 gb=g['max_degree'], sd_direct=sd, predicted=predicted, ok=sd == predicted))
    # C3: brute-force Z from the generator on the chain grid (prime S1)
    rng = random.Random(SEED + 3)
    s1 = load_sweep_like('ic-curve-sweep')
    for s, res in rng.sample([x for x in s1 if x[0]['kind'] == 'prime'], 100):
        k, p, g0 = s['k'], s['prime'], s['generator_value']
        N = 1 << k
        zs = []
        for E1 in range(N):
            for E2 in range(N):
                val = 0
                for mask, c in s['generators'][0]:
                    a1, a2 = mask & (N - 1), mask >> k
                    val += c * pow(g0, (a1 * E1 + a2 * E2) % N, p)
                if val % p == 0:
                    zs.append(E1 | (E2 << k))
        report['C3'].append(dict(id=s['system_hash'][:16], ok=sorted(zs) == sorted(s['coordinate_solutions'])))
    # C4: |Z| <= 2^{k+1} on two-summand systems (rows.jsonl)
    rows = jl(OUT / 'rows.jsonl')
    viol, checked = [], 0
    for r_ in rows:
        if r_['m'] == 2 and r_['k'] and r_['points'] is not None and r_['set'] in ('S1', 'S2', 'S5', 'S6', 'S7') \
                and (r_['family'] == 'prime' or r_['r'] == 2 * r_['k']):
            checked += 1
            if r_['points'] > 2 ** (r_['k'] + 1):
                viol.append(r_['id'])
    report['C4'] = dict(checked=checked, violations=viol, ok=not viol)
    summary = {c: (dict(n=len(v), ok=sum(x['ok'] for x in v)) if isinstance(v, list) else v)
               for c, v in report.items() if c != 'not_run'}
    summary['not_run'] = report['not_run']
    report['summary'] = summary
    (OUT / 'checks.json').write_text(json.dumps(report, indent=1, sort_keys=True))
    print(json.dumps(summary, indent=1))


# ------------------------------------------------------------------ C1b: PolyBoRi on every computed binary system
def polybori_export(_a):
    """polybori-all-in.json: every binary system of S1 and S5 with a computed max.GB.deg (graded-lex)."""
    sys.path.insert(0, str(ROOT / 'scripts'))
    import run_ic_solving_degree as RUN
    rows = {r['id']: r for r in jl(OUT / 'rows.jsonl')}
    out = []
    for row in RUN.binary_results():
        if rows.get(row['id'], {}).get('gb_grlex') is not None:
            out.append(dict(id=row['id'], r=row['r'], equations=row['equations'], points=row['field_points']))
    for s, res in load_sweep_like('ic-curve-sweep'):
        if s['kind'] != 'prime' and rows.get(s['system_hash'][:16], {}).get('gb_grlex') is not None:
            out.append(dict(id=s['system_hash'], r=s['variables'], equations=s['equations'],
                            points=s['coordinate_solutions']))
    (OUT / 'polybori-all-in.json').write_text(json.dumps(out))
    print('exported', len(out))


def polybori_compare(_a):
    """Compare polybori-all-out.jsonl (sage polybori_gb.sage) with the grlex staircase and rows.jsonl.
    The Boolean ring leaves the field equations implicit: in K[u] the reduced basis also contains
    u_v^2 + u_v for every standard variable u_v, which has degree 2."""
    inp = {s['id']: s for s in json.load(open(OUT / 'polybori-all-in.json'))}
    rows = {r['id']: r for r in jl(OUT / 'rows.jsonl')}
    n = ok = 0
    bad = []
    for o in jl(OUT / 'polybori-all-out.jsonl'):
        s = inp[o['id']]
        n += 1
        st = SD.staircase(SD.BooleanPoints(s['points'], s['r']), len(s['points']), s['r'], 2, order='grlex')
        rid = o['id'] if o['id'].startswith('m') else o['id'][:16]
        pb = max(o['max_degree'], 2 if st['standard_variables'] else 0)
        lead = sorted(st['squarefree_generators']) == o['leading_masks']
        good = st['max_gb_degree'] == pb == rows[rid]['gb_grlex'] and lead
        ok += good
        if not good:
            bad.append(dict(id=o['id'][:24], bm=st['max_gb_degree'], polybori=o['max_degree'],
                            polybori_with_squares=pb, recorded=rows[rid]['gb_grlex'], leading_equal=lead))
    rep = dict(check='C1b: PolyBoRi reduced deglex basis (plus u_v^2+u_v for standard variables) vs the '
                     'Buchberger-Moeller grlex staircase (max degree, leading monomials) and the recorded gb_grlex; '
                     'every computed binary system of S1 and S5', n=n, ok=ok, mismatches=bad)
    (OUT / 'checks-polybori-all.json').write_text(json.dumps(rep, indent=1))
    print(n, ok, bad[:5])


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('cmd', choices=['export', 'compare', 'polybori-export', 'polybori-compare'])
    a = ap.parse_args()
    {'export': export, 'compare': compare, 'polybori-export': polybori_export,
     'polybori-compare': polybori_compare}[a.cmd](a)
