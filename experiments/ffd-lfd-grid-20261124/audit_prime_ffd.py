"""Bind each frozen prime sweep system to the kernel-checked FFD theorem."""
from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = Path(__file__).resolve().parents[2]
SWEEP = REPO / 'experiments/ic-curve-sweep'


def sha(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generator_from_case(c):
    p, k, a, b, t = c['prime'], c['k'], c['a'], c['b'], c['target'][0]
    scale = c['base_scale']
    N = 1 << k
    coeff = {
        (2, 0): t*t, (0, 2): t*t, (1, 1): -2*t*t-2*a,
        (2, 1): -2*t, (1, 2): -2*t,
        (1, 0): -2*a*t-4*b, (0, 1): -2*a*t-4*b,
        (0, 0): a*a-4*b*t, (2, 2): 1,
    }
    output = {}
    for (e1, e2), value in coeff.items():
        mask = (e1 % N) | ((e2 % N) << k)
        output[mask] = (output.get(mask, 0) + value * pow(scale, e1+e2, p)) % p
    return {mask: value for mask, value in output.items() if value}


def main():
    paths = ('cases.jsonl', 'systems.jsonl', 'results.jsonl')
    input_hashes = {name: sha(SWEEP / name) for name in paths}
    theorem_sources = {name: sha(REPO / name) for name in (
        'formal/ICFallDegrees/PrimeFieldConclusions.lean',
        'formal/ICFallDegrees/PrimeM2General.lean',
        'formal/ICFallDegrees/GeneralQuadraticClassification.lean')}
    systems = {s['system_hash']: s for s in map(json.loads, (SWEEP/'systems.jsonl').open())
               if s['kind'] == 'prime'}
    results = {r['system_hash']: r for r in map(json.loads, (SWEEP/'results.jsonl').open())}
    cases = defaultdict(list)
    for c in map(json.loads, (SWEEP/'cases.jsonl').open()):
        if c['kind'] == 'prime':
            assert c['system_hash'] in systems
            cases[c['system_hash']].append(c)
    assert set(cases) == set(systems) and len(systems) == 1260
    cells = defaultdict(set)
    ffd = {}
    for key, system in systems.items():
        assert system['m'] == 2 and system['k'] >= 2 and system['prime'] >= 5
        assert system['variables'] == 2 * system['k']
        saved = dict((int(mask), int(value)) for mask, value in system['generators'][0])
        assert len(system['generators']) == 1
        for case in cases[key]:
            assert case['m'] == 2 and case['prime'] == system['prime']
            assert case['k'] == system['k']
            assert generator_from_case(case) == saved
            p, a, t, scale = case['prime'], case['a'], case['target'][0], case['base_scale']
            assert scale % p != 0
            # S3(scale*u,scale*v,t;a,b) = scale^4 S3(u,v,t/scale;a/scale^2,b/scale^3).
            # The normalized input differs only by this nonzero input scalar.
            degree = 3 if (3*t*t+a) % p == 0 else 4
            assert ffd.get(key, degree) == degree
            ffd[key] = degree
            cells[(p, system['m'])].add(key)
        r = results[key]
        assert r['lfd_lower_bound'] == r['lfd_upper_bound'] == r['exact_lfd']
        # An independent finite consistency check against the closure profile.
        assert r['strict_falls'] and r['strict_falls'][0] == ffd[key]
    table = []
    for (p, m), keys in sorted(cells.items()):
        fcount = Counter(ffd[k] for k in keys)
        lcount = Counter(results[k]['exact_lfd'] for k in keys)
        table.append({'field_kind': 'prime', 'p': p, 'n': 1, 'm': m,
                      'systems': len(keys),
                      'ffd_exact_counts': {str(k): v for k, v in sorted(fcount.items())},
                      'lfd_exact_counts': {str(k): v for k, v in sorted(lcount.items())}})
    out = {'status': 'passed', 'scope': 'frozen prime m=2 sweep systems',
           'systems': len(systems), 'cells': table,
           'input_sha256': input_hashes, 'theorem_source_sha256': theorem_sources,
           'theorem': 'ICFallDegrees.PrimeFieldConclusions.ordinary_ffd',
           'theorem_condition': 'FFD=3 iff 3*t^2+a=0 mod p; otherwise FFD=4',
           'coset_transport': 'nonzero equal scale on both free abscissas; divide S3 by scale^4',
           'audit_source_sha256': sha(Path(__file__))}
    dest = HERE/'audit-prime-ffd.json'
    dest.write_text(json.dumps(out, sort_keys=True, indent=2) + '\n')
    print(dest, 'systems', len(systems), 'cells', len(table))


if __name__ == '__main__':
    main()
