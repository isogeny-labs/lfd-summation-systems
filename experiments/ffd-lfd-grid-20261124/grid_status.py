"""Reconcile the frozen FFD/LFD manifests with certificate receipts."""
from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bounds(counts, *, nullable=False):
    values = [int(k) for k in counts if k != 'None']
    assert values and sum(counts.values()) > 0
    if not nullable:
        assert 'None' not in counts
    return {'lower': min(values), 'upper': max(values),
            'undefined': counts.get('None', 0), 'histogram': counts}


def verified_f4_upper(p, d8):
    names = [f'f4-p{p}-{kind}.json' for kind in ('result', 'trace', 'replay')]
    present = [(HERE / name).exists() for name in names]
    if not any(present):
        return None
    assert all(present), f'incomplete F4 receipt set for {p}'
    producer, trace, replay = (read(name) for name in names)
    expected_sources = {
        'f4_min/src/main.rs': '06a47329aa402a02e70589c1c30a5465223a6f17b2a6d2dbb5df5813ec1476b7',
        'f4_min/src/f4_fp_tower.rs': 'a4bbf68f23912cbda4e955f83debdd5af2641be05a0aba0f23ce4dbb3bb6e6bb',
        'f4-linux-x86_64': '60e1fff70e9c67b1dd0ec778d12e5e8566a29e9d30e1d45bd86ea4771f623a11',
    }
    assert all(sha(HERE / name) == digest for name, digest in expected_sources.items())
    input_digest = sha(d8)
    for receipt in (producer, replay):
        assert (receipt['p'], receipt['m'], receipt['k'], receipt['degree_cap']) == (p, 5, 3, 8)
        assert receipt['input_sha256'] == input_digest
        assert receipt['inconsistent'] is True
        assert receipt['solving_degree_max'] == 8
        assert receipt['timed_out'] is False and receipt['oversize'] is False
        assert receipt['replay_failure'] is None
    assert producer['replay_mode'] is False and replay['replay_mode'] is True
    assert trace['input_sha256'] == input_digest
    assert (trace['p'], trace['m'], trace['k'], trace['inconsistent'],
            trace['solving_degree_max']) == (p, 5, 3, True, 8)
    assert producer['learned_steps'] == len(trace['steps']) == replay['steps']
    trace_digest = sha(HERE / names[1])
    assert producer['trace_sha256'] == replay['replayed_trace_sha256'] == trace_digest
    return {'upper': 8, 'status': 'certified', 'method': 'exact F4 trace replay',
            'producer_sha256': sha(HERE / names[0]),
            'trace_sha256': trace_digest,
            'replay_sha256': sha(HERE / names[2]),
            'source_sha256': expected_sources}


def main():
    binary = read('audit-binary-ffd.json')
    prime = read('audit-prime-ffd.json')
    semaev = read('audit-semaev-ffd.json')
    split = read('audit-split-ffd.json')
    split8 = read('audit-split-m5k3-ffd.json')
    assert all(r['status'] == 'passed' for r in (binary, prime, semaev, split, split8))
    log = HERE / 'lfd-sweep-replay-20261124.log'
    replay = json.loads(log.read_text().splitlines()[-1])
    assert replay == {'status': 'passed: cases regenerated and every exact result replayed',
                      'cases': 3163, 'systems': 3005}

    curve_cells = []
    for audit in (binary, prime):
        for c in audit['cells']:
            assert sum(c['lfd_exact_counts'].values()) == c['systems']
            assert sum(c['ffd_exact_counts'].values()) == c['systems']
            curve_cells.append({'field_kind': c['field_kind'],
                                'field': f"GF({c['p']})" if 'p' in c else f"GF(2^{c['n']})",
                                'n': c['n'], 'm': c['m'], 'systems': c['systems'],
                                'ffd': bounds(c['ffd_exact_counts'], nullable=True),
                                'lfd': bounds(c['lfd_exact_counts'])})
    assert len(curve_cells) == 21
    assert binary['systems'] + prime['systems'] == 3005
    # Ten identical binary systems are assigned to both the binary and
    # Koblitz field cells; one of them is a constant-input system.
    assert sum(c['systems'] for c in curve_cells) == 3015
    assert sum(c['ffd']['undefined'] for c in curve_cells) == 65
    binary_rows = [json.loads(s) for s in (HERE/'binary-ffd.jsonl').read_text().splitlines()[1:]]
    assert len(binary_rows) == binary['systems']
    assert sum(r['result']['status'] == 'constant-input' for r in binary_rows) == 64

    semaev_cells = []
    for n, joint in sorted(semaev['ffd_lfd_counts_by_n'].items(), key=lambda x: int(x[0])):
        fc, lc = Counter(), Counter()
        for pair, count in joint.items():
            f, l = map(int, pair.split(','))
            fc[str(f)] += count
            lc[str(l)] += count
        semaev_cells.append({'field_kind': 'binary-fixed-six-coordinate-S3',
                             'field': f'GF(2^{n})', 'n': int(n), 'm': 2,
                             'systems': sum(joint.values()),
                             'ffd': bounds(dict(fc)), 'lfd': bounds(dict(lc)),
                             'ffd_lfd_joint': joint})
    assert len(semaev_cells) == 6
    assert sum(c['systems'] for c in semaev_cells) == semaev['cases'] == 1977

    split_cells = []
    for c in split['cells']:
        assert sum(c['fields'].values()) == c['systems']
        assert sum(c['lfd_exact_counts'].values()) == c['systems']
        for field, count in sorted(c['fields'].items()):
            # The per-system exact LFD and FFD receipts are retained in split['systems'].
            systems = [s for s in split['systems'] if
                       (s['m'],s['k'],s['p']) == (c['m'],c['k'],int(field))]
            assert len(systems) == count
            assert all(s['ffd_lower_bound'] == s['ffd_upper_bound'] == 4 for s in systems)
            assert all(s['lfd_lower_bound'] == s['lfd_upper_bound'] for s in systems)
            lc = Counter(str(s['lfd_upper_bound']) for s in systems)
            split_cells.append({'field_kind': 'prime-split-chain', 'field': f'GF({field})',
                                'n': 1, 'm': c['m'], 'k': c['k'], 'systems': count,
                                'ffd': bounds({'4': count}), 'lfd': bounds(dict(lc))})
    assert len(split_cells) == 15 and len(split['systems']) == 24

    assert len(split8['cells']) == 2
    for c in split8['cells']:
        p = c['p']
        assert (c['m'],c['k'],c['ffd_lower_bound'],c['ffd_upper_bound']) == (5,3,4,4)
        d7 = HERE / f'split-m5-k3-p{p}-D7.txt'
        d8 = HERE / f'split-m5-k3-p{p}-D8.txt'
        producer = read(f'split-m5-k3-D7-p{p}-producer.json')
        verifier = read(f'split-m5-k3-D7-p{p}-replay.json')
        assert d7.read_text().splitlines()[1:] == d8.read_text().splitlines()[1:]
        assert producer['input_sha256'] == verifier['input_sha256']
        assert producer['original_input_sha256'] == sha(d7)
        assert (producer['status'],producer['p'],producer['m'],producer['k'],producer['D'],
                producer['contains_one']) == ('closed',p,5,3,7,False)
        assert (verifier['status'],verifier['closed_check'],verifier['contains_one']) == (
            'passed',True,False)
        upper = verified_f4_upper(p, d8)
        lfd = {'lower': 8, 'upper': None, 'status': 'open',
               'D7_input_sha256': sha(d7), 'D7_replay_sha256':
               sha(HERE / f'split-m5-k3-D7-p{p}-replay.json')}
        if upper is not None:
            lfd.update(upper)
        split_cells.append({'field_kind': 'prime-split-chain', 'field': f'GF({p})',
                            'n': 1, 'm': 5, 'k': 3, 'systems': 1,
                            'ffd': bounds({'4': 1}),
                            'lfd': lfd,
                            'D8_input_sha256': sha(d8)})
    assert len(split_cells) == 17
    missing = [{'field': c['field'], 'n': c['n'], 'm': c['m'], 'k': c['k'],
                'needed': 'verified LFD upper bound'}
               for c in split_cells if c['lfd']['upper'] is None]
    assert 0 <= len(missing) <= 2
    out = {'status': 'complete' if not missing else 'incomplete',
           'scope': 'frozen curve sweep, fixed-basis Semaev panel, and frozen split protocol; publication depth14/depth18 continuations tracked separately',
           'frozen_manifest_cells_fully_certified_where_FFD_defined': not missing,
           'tracks': {'frozen_curve_sweep': curve_cells,
                      'fixed_basis_semaev': semaev_cells,
                      'frozen_split_protocol': split_cells},
           'missing': missing,
           'constant_input_FFD_undefined_distinct': 64,
           'constant_input_cell_assignments': 65,
           'curve_sweep_distinct_systems': 3005,
           'curve_sweep_cell_assignments': 3015,
           'input_receipt_sha256': {name: sha(HERE/name) for name in (
               'audit-binary-ffd.json','audit-prime-ffd.json','audit-semaev-ffd.json',
               'audit-split-ffd.json','audit-split-m5k3-ffd.json',
               'lfd-sweep-replay-20261124.log')},
           'audit_source_sha256': sha(Path(__file__))}
    dest = HERE / 'grid-status.json'
    dest.write_text(json.dumps(out, sort_keys=True, indent=2) + '\n')
    print(json.dumps({'status': out['status'], 'cells': {k:len(v) for k,v in out['tracks'].items()},
                      'missing': missing, 'constant_input_FFD_undefined_distinct': 64}))


if __name__ == '__main__':
    main()
