"""Summarize experiments/ic-bounds-freeext/records/*.json into table.md and summary.json.

Per record (one instance, one degree D): for every computed base depth k0, whether
W_{k0+1}(V_D(k0)) is closed (deterministic sweep), contains the generator rows, equals
V_D(k0+1); the least valid base depth (closed AND equal AND 1 notin V_D(k0) AND g in V_D(k0));
the first depth where V_D(k) != I_k cap M_<=D (only depths with 2^k | p-1); and whether the
profile recursion p_j(k+1) = p_j(k) + 2 p_{j-1}(k) + p_{j-2}(k) holds from that depth on."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'experiments/ic-bounds-freeext'


def summarize(r):
    ks = sorted(int(k) for k in r['closures'])
    cl = r['closures']
    ext = r.get('extensions', {})
    rows = []
    for k in ks:
        c = cl[str(k)]
        e = ext.get(str(k))
        rows.append(dict(k=k, rank=c['rank'], profile=c['profile'], one=c['contains_one'], ideal_dim=c['ideal_dim'],
                         deficient=c['deficient'], points=c['points'],
                         W_closed=None if e is None else e['W_closed'],
                         W_gen_outside=None if e is None else e['W']['generator_rows_outside'],
                         W_sweep_bad=None if e is None else e['W']['sweep_full']['nonzero_residuals'],
                         W_equals_next=None if e is None else e.get('equals_V_next', e['W'].get('equals_reference_V_next')),
                         W_equals_by=None if e is None else (e.get('equals_V_next_by') or
                                                             ('recorded V_D(k+1) SHA-256' if 'equals_reference_V_next' in e['W'] else None)),
                         recursion_to_next=None if e is None else e.get('recursion_holds'),
                         hyp=None if e is None else e['kernel_hypotheses_hold']))
    valid = [x['k'] for x in rows if x['W_closed'] and x['W_equals_next'] and not x['one']]
    closed_hyp = [x['k'] for x in rows if x['hyp']]
    fall = next((x['k'] for x in rows if x['deficient']), None)
    rec_from = None
    if fall is not None:
        rr = [x['recursion_to_next'] for x in rows if x['k'] >= fall and x['recursion_to_next'] is not None]
        rec_from = all(rr) if rr else None
    wit = {k: (v['found'][0]['e'] if v['found'] else None) for k, v in r.get('witness', {}).items()}
    return dict(name=None, p=r['p'], a=r['a'], b=r['b'], seed=r['seed'], D=r['D'], split=r['split'],
                X=r['target_info'].get('target', [r['target_info'].get('X')])[0], exceptional=r['exceptional'],
                depths=ks, least_valid_base=min(valid) if valid else None,
                least_base_kernel_hypotheses=min(closed_hyp) if closed_hyp else None,
                first_deficient_depth=fall, recursion_from_fall=rec_from, witnesses=wit, rows=rows)


def main():
    recs = sorted((OUT / 'records').glob('*.json'))
    out = []
    lines = ['| record | D | X | exceptional | least base k0 (closed & equal) | least k0 with kernel hypotheses | first deficient depth | recursion from there | per-k0: closed / equal |',
             '|---|---|---|---|---|---|---|---|---|']
    for f in recs:
        r = json.loads(f.read_text())
        s = summarize(r)
        s['name'] = f.stem
        out.append(s)
        exc = ','.join(k for k, v in s['exceptional'].items() if v) or '-'
        per = '; '.join(f"{x['k']}: {'C' if x['W_closed'] else 'nc'}/{'=' if x['W_equals_next'] else ('?' if x['W_equals_next'] is None else '≠')}"
                        for x in s['rows'] if x['W_closed'] is not None)
        lines.append(f"| {f.stem} | {s['D']} | {s['X']} | {exc} | {s['least_valid_base']} | {s['least_base_kernel_hypotheses']} | "
                     f"{s['first_deficient_depth']} | {s['recursion_from_fall']} | {per} |")
    (OUT / 'summary.json').write_text(json.dumps(out, indent=1) + '\n')
    (OUT / 'table.md').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    sys.exit(main())
