#!/usr/bin/env python3
"""Compare the certified split-presentation values of experiments/ic-bounds-m3
with the stage B readings of the independent Rust instrument (crypto repository),
copied to experiments/crossrepo-pkm-stage-b/source/.

It checks the copied files against their recorded digests, extracts both
tables, and requires agreement on every depth at which both exist. It runs
no closure. Output: experiments/crossrepo-pkm-stage-b/comparison.json.
"""
import collections
import hashlib
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
HERE = ROOT / 'experiments/crossrepo-pkm-stage-b'
M3 = ROOT / 'experiments/ic-bounds-m3'


def certified_split():
    """Depth -> histogram of the exact LFD of the free split encoding, unit ideals only."""
    cases = {}
    for line in (M3 / 'cases.jsonl').read_text().splitlines():
        case = json.loads(line)
        cases[case['case_id']] = case
    seen, table, proper = set(), collections.defaultdict(collections.Counter), collections.Counter()
    for path in sorted(M3.glob('results*.jsonl')):
        for line in path.read_text().splitlines():
            row = json.loads(line)
            if row['encoding'] != 'split' or row['case_id'] in seen:
                continue
            seen.add(row['case_id'])
            free = row['results']['free']
            if free['exact_lfd'] is None:
                continue
            depth = cases[row['case_id']]['k']
            if free['certificate']['kind'] == 'unit':
                table[depth][free['exact_lfd']] += 1
            else:
                proper[depth] += 1
    return table, proper


def stage_b():
    """(m, depth) -> set of per-draw readings of the Kummer tower arm, both primes."""
    score = json.loads((HERE / 'source/score.json').read_text())
    table = collections.defaultdict(list)
    for cell in score['cells']:
        if cell['kind'] == 'kummer' and cell['control'] == 'tower':
            table[(cell['m'], cell['t'])] += [(cell['p'], cell['d_F'], cell['ffd'], cell['D'])]
    return table


def main():
    failures = []
    source = json.loads((HERE / 'SOURCE.json').read_text())
    for name, digest in source['sha256'].items():
        if hashlib.sha256((HERE / 'source' / name).read_bytes()).hexdigest() != digest:
            failures.append(f'digest mismatch: {name}')
    certified, proper = certified_split()
    readings = stage_b()
    overlap = []
    for depth in sorted(certified):
        cells = readings.get((3, depth))
        if not cells:
            continue
        values = set(certified[depth])
        draws = [int(v) for _, d_f, _, _ in cells for v in d_f.split() if v.isdigit()]
        censored = [v for _, d_f, _, _ in cells for v in d_f.split() if not v.isdigit()]
        agree = len(values) == 1 and set(draws) == values
        overlap.append(dict(m=3, depth=depth, certified_unit_histogram=dict(certified[depth]),
                            certified_proper_cases=proper[depth], stage_b_draws=draws,
                            stage_b_censored=censored, agree=agree))
        if not agree:
            failures.append(f'm=3 depth {depth}: certified {dict(certified[depth])}, stage B {draws}')
    if [row['depth'] for row in overlap] != [2, 3, 4]:
        failures.append('expected overlap at depths 2, 3 and 4')
    first_fall = sorted({v for cells in readings.values() for _, _, ffd, _ in cells for v in ffd.split()})
    if first_fall != ['4']:
        failures.append(f'first fall degree not constant 4: {first_fall}')
    beyond = {f'm={m},depth={t}': [dict(p=p, d_F=d, D=big) for p, d, _, big in cells]
              for (m, t), cells in sorted(readings.items()) if (m, t) not in {(3, 2), (3, 3), (3, 4)}}
    amendment = []
    for name in ('M4-t4-k0-c7.jsonl', 'M4-t4-k1-c7.jsonl'):
        row = json.loads((HERE / 'source' / name).read_text().splitlines()[0])
        ones = [c['c'] for c in row['closures'] if c['one']]
        amendment.append(dict(file=name, m=4, depth=4, first_degree_with_one=min(ones) if ones else None,
                              closures_without_one=[c['c'] for c in row['closures'] if not c['one']]))
        if not ones or min(ones) != 7:
            failures.append(f'{name}: expected 1 first at degree 7')
    out = dict(status='failed' if failures else 'passed',
               scope='comparison of recorded values; no closure was run; the stage B readings carry no certificate',
               overlap=overlap, stage_b_first_fall_degrees=first_fall,
               stage_b_only=beyond, four_summands_depth_four=amendment, failures=failures)
    (HERE / 'comparison.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({k: out[k] for k in ('status', 'overlap', 'four_summands_depth_four', 'failures')}, indent=2))
    return bool(failures)


if __name__ == '__main__':
    sys.exit(main())
