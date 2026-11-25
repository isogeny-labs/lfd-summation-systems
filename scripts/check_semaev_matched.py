#!/usr/bin/env python3
"""Replay exact controls and verify the saved matched-benchmark record structure.

Does not reproduce timing or RSS. Source hashes bind the measured implementation.
"""
from collections import Counter
import contextlib
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_semaev_matched as RUN


def main():
    out = RUN.OUT
    env = json.loads((out/'environment.json').read_text())
    for path, digest in env['source_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest, path
    for path, key in [('PROTOCOL.md','protocol_sha256'), ('instances.json','instances_sha256')]:
        assert hashlib.sha256((out/path).read_bytes()).hexdigest() == env[key], path
    with tempfile.TemporaryDirectory() as tmp:
        target = Path(tmp)
        with contextlib.redirect_stdout(io.StringIO()):
            RUN.audit(target)
            regenerated = RUN.instances(target)
        assert json.loads(json.dumps(regenerated)) == json.loads((out/'instances.json').read_text())
        assert json.loads((target/'correctness.json').read_text()) == json.loads((out/'correctness.json').read_text())
        runs = [json.loads(s) for s in (out/'runs.jsonl').read_text().splitlines()]
        expected = {(c['id'], f, rep, False) for c in regenerated['cases']
                    for f in ('original','split') for rep in range(3)}
        expected |= {(c['id'], f, None, True) for c in regenerated['cases']
                     for f in ('original','split')}
        keys = [(r['case_id'],r['formulation'],r['repetition'],r['diagnostic']) for r in runs]
        assert len(keys) == len(set(keys)) and set(keys) == expected
        oracle_counts = {c['id']:len(c['oracle']) for c in regenerated['cases']}
        for r in runs:
            if r['status']=='ok' and not r['diagnostic']:
                assert r['projected_solutions'] == oracle_counts[r['case_id']]
                for key in ('setup_seconds','gb_seconds','recovery_seconds','peak_rss_kib'):
                    assert r[key] >= 0
            if r['diagnostic'] and r['status']=='ok':
                assert all(r['diagnostics'][k] is None for k in
                           ('matrix_rows','matrix_columns','matrix_rank','matrix_sparsity'))
        (target/'instances.json').write_bytes((out/'instances.json').read_bytes())
        (target/'runs.jsonl').write_bytes((out/'runs.jsonl').read_bytes())
        with contextlib.redirect_stdout(io.StringIO()):
            RUN.summarize(target)
        assert json.loads((target/'summary.json').read_text()) == json.loads((out/'summary.json').read_text())
    print(json.dumps(dict(status='passed', records=len(runs),
                          statuses=dict(Counter(r['status'] for r in runs)),
                          scope='exact audit, case selection, provenance, structure and summary; not timing replay')))


if __name__ == '__main__':
    main()
