#!/usr/bin/env python3
"""Summarize the verification-pass closure re-runs in experiments/ic-bounds-binary/verify/.

  python3 scripts/ic_bounds_binary_verify_summary.py > experiments/ic-bounds-binary/verify/summary.json
"""
import collections
import hashlib
import json
import os

V = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments', 'ic-bounds-binary', 'verify')
FLAGS = ['hash_match', 'points_match', 'lfd_match', 'dstar_match', 'dims_match']


def main():
    out = {}
    for fam in ('m2', 'm3-direct', 'm3-split'):
        p = os.path.join(V, f'{fam}.jsonl')
        if not os.path.exists(p):
            continue
        rows = [json.loads(l) for l in open(p)]
        stored = [json.loads(l) for l in open(os.path.join(V, '..', f'results-{fam}.jsonl'))]
        d = dict(file_sha256=hashlib.sha256(open(p, 'rb').read()).hexdigest(), cases_checked=len(rows),
                 cases_stored=len(stored), closures_run=sum(r.get('closure_from') is None for r in rows),
                 distinct_hashes=len({r['my_system_hash'] for r in rows}))
        for f in FLAGS:
            d[f'{f}_false'] = [r['case_id'] for r in rows if not r.get(f)]
        d['solutions_mismatch'] = [r['case_id'] for r in rows if r['solutions'] != r['stored_solutions']]
        d['unit_mismatch'] = [r['case_id'] for r in rows if r['unit'] != (r['stored_status'] == 'exact-unit')]
        d['falls_mismatch'] = [r['case_id'] for r in rows if r['falls'] != [x for x in (r['stored_falls'] or [])
                                                                             if x <= (r['dstar'] or 99)]]
        nv = [r for r in rows if 'naive_match' in r and r.get('closure_from') is None]
        d['naive_checked'] = len(nv)
        d['naive_false'] = [r['case_id'] for r in nv if not r['naive_match']]
        ex = [r for r in rows if r.get('closure_from') is None and len(r['dims']) - 1 > (r['dstar'] or 0)]
        d['extra_degree_scanned'] = len(ex)
        d['fall_after_certificate'] = [r['case_id'] for r in rows if r.get('fall_after_certificate')]
        per = {}
        for k in sorted({r['k'] for r in rows}):
            rr = [r for r in rows if r['k'] == k]
            per[k] = dict(cases=len(rr), max_lfd=max(r['lfd'] for r in rr),
                          hist=dict(sorted(collections.Counter(r['lfd'] for r in rr).items())),
                          dstar_minus_lfd=dict(collections.Counter(r['dstar'] - r['lfd'] for r in rr)),
                          max_dstar=max(r['dstar'] for r in rr),
                          seconds=round(sum(r['total_seconds'] for r in rr if r.get('closure_from') is None), 1),
                          peak_rss_mb=max(r['peak_rss_mb'] for r in rr))
        d['per_k'] = per
        d['missing_case_ids'] = len({r['case_id'] for r in stored} - {r['case_id'] for r in rows})
        out[fam] = d
    print(json.dumps(out, indent=1, sort_keys=True))


if __name__ == '__main__':
    main()
