"""Item 4: summarize experiments/ic-bounds-exact10 (validation, profiles, mid-scale, runner tests).

Writes summary.json (agreement counts, per-row table), calibration.json (cols, seconds of reached
split closures, used by scripts/ic_bounds_exact10_plan.py), hashes.json (sha256 of the sources and
data files) and cost.svg / cost.png (memory model of the dense split engine against depth, with the
mid-scale measurements).
"""
from __future__ import annotations

import glob
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
OUT = ROOT / 'experiments/ic-bounds-exact10'


def load(pattern):
    rows = []
    for path in sorted(glob.glob(str(OUT / pattern))):
        rows += [json.loads(line) for line in open(path) if line.strip()]
    return rows


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# block-engine file each validation batch imported (csd/ic_chain_block_closure.py was edited by the
# main session at 22:36 and again at 23:31:40 on 2026-11-07; the second edit changed only the default
# float64 tile of _modmat, an exact computation either way)
BLOCK_VERSIONS = {
    'k56': '865032b874b2e88e692d913c747f8add27e5ca5d36593b9325d4231c3c748db9',
    'k7a': '865032b874b2e88e692d913c747f8add27e5ca5d36593b9325d4231c3c748db9',
    'k7b': '865032b874b2e88e692d913c747f8add27e5ca5d36593b9325d4231c3c748db9',
    'k7c': '865032b874b2e88e692d913c747f8add27e5ca5d36593b9325d4231c3c748db9',
    'k7d': '865032b874b2e88e692d913c747f8add27e5ca5d36593b9325d4231c3c748db9',
    'k7e': 'de5e39705f88b6d975a05db056208b6481b220f00d5bfcb395efb7cb3a6329ca',
    'k56-v2': 'de5e39705f88b6d975a05db056208b6481b220f00d5bfcb395efb7cb3a6329ca',
}


def by_batch(prefix):
    out = {}
    for path in sorted(glob.glob(str(OUT / f'{prefix}-*.jsonl'))):
        tag = Path(path).stem[len(prefix) + 1:]
        rows = [json.loads(line) for line in open(path) if line.strip()]
        ok = sum(1 for r in rows if (r['all_pass'] if 'all_pass' in r else (r['lfd_agree'] and r['falls_agree'])))
        out[tag] = dict(rows=len(rows), passed=ok, block_engine_sha256=BLOCK_VERSIONS.get(tag))
    return out


def main():
    # headline counts use one run per (instance, k, D): the k56-v2 repeat is reported per batch only
    val = [r for path in sorted(glob.glob(str(OUT / 'validation-*.jsonl'))) if not path.endswith('-v2.jsonl')
           for r in (json.loads(line) for line in open(path) if line.strip())]
    prof = [r for path in sorted(glob.glob(str(OUT / 'profiles-*.jsonl'))) if not path.endswith('-v2.jsonl')
            for r in (json.loads(line) for line in open(path) if line.strip())]
    mid = load('midscale-k*.jsonl')
    checks = {}
    table = []
    for r in val:
        for key, ok in r['checks'].items():
            c = checks.setdefault(key, [0, 0])
            c[0] += bool(ok)
            c[1] += 1
        table.append(dict(instance=r['instance'], k=r['k'], D=r['degree'], points=r['points'], cols=r['cols'],
                          ref_source=r['ref_source'], ref_dim=r['ref_dim'], ref_one=r['ref_contains_one'],
                          block_one=r.get('block_engine', {}).get('contains_one'),
                          block_rank=r.get('block_engine', {}).get('rank'),
                          split_status=r['split_scratch']['status'], split_dim=r['split_scratch']['exact_dimension'],
                          split_low=r['split_scratch']['low_dimension'],
                          prev_dim=r.get('split_from_prev_sparse', {}).get('exact_dimension'),
                          neg_missing=r.get('negative_control', {}).get('missing'), all_pass=r['all_pass']))
    one_false = [t for t in table if t['ref_one'] is False]
    one_false_unit = [t for t in one_false if t['points'] == 0]
    summary = dict(
        validation_rows=len(val), validation_all_pass=sum(1 for r in val if r['all_pass']),
        validation_rows_without_reference=sum(1 for r in val if r['ref_dim'] is None),
        checks={k: dict(passed=v[0], total=v[1]) for k, v in sorted(checks.items())},
        rows_with_1_not_in_VD=len(one_false), unit_rows_with_1_not_in_VD=len(one_false_unit),
        proper_rows=sum(1 for t in table if t['points']), unit_rows=sum(1 for t in table if not t['points']),
        block_engine_rows=sum(1 for t in table if t['block_one'] is not None),
        block_engine_undecided_rows=sum(1 for t in table if t['block_one'] is False),
        block_engine_stall_equals_dense=sum(1 for t in table if t['block_one'] is False and t['block_rank'] == t['ref_dim']),
        profiles=len(prof), profiles_lfd_agree=sum(1 for r in prof if r['lfd_agree']),
        profiles_falls_agree=sum(1 for r in prof if r['falls_agree']),
        profile_lfd_by_k={str(k): sorted({r['split']['exact_lfd'] for r in prof if r['k'] == k}) for k in
                          sorted({r['k'] for r in prof})},
        midscale_rows=len(mid), midscale_all_pass=sum(1 for r in mid if r['all_pass']),
        validation_batches=by_batch('validation'), profile_batches=by_batch('profiles'),
        table=table)
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    runs = []
    for r in val:
        s = r['split_scratch']  # noqa: shadows nothing used below
        if s['status'] == 'reached' and s['seconds'] > 0:
            runs.append(dict(source='validation-scratch', k=r['k'], D=r['degree'], cols=r['cols'], seconds=s['seconds']))
    for r in mid:
        if r['D5']['status'] == 'reached':
            runs.append(dict(source='midscale-D5', k=r['k'], D=5, cols=r['D5']['cols'], seconds=r['D5']['seconds']))
    for path in sorted(glob.glob(str(OUT / 'degree*-k*.json'))):
        d = json.loads(open(path).read())
        if d['result']['status'] == 'reached':
            total = d.get('seconds_total')
            log = OUT / 'logs' / f"runner-k{d['instance']['k']}-D{d['degree']}.log"
            if total is None and log.exists():  # runs recorded before seconds_total: sum the segments' 'done' lines
                total = sum(json.loads(line)['elapsed'] for line in open(log)
                            if line.startswith('{') and '"stage": "done"' in line)
            runs.append(dict(source=Path(path).name, k=d['instance']['k'], D=d['degree'], cols=d['result']['cols'],
                             seconds=total or d['seconds']))
    (OUT / 'calibration.json').write_text(json.dumps(dict(
        note='wall seconds of reached split closures on a shared 4-core host (load about 4), single BLAS thread; '
             'used only for an order-of-magnitude extrapolation', runs=runs), indent=1) + '\n')
    files = sorted(glob.glob(str(ROOT / 'csd/ic_bounds_exact10*.py')) + glob.glob(str(ROOT / 'scripts/ic_bounds_exact10*.py'))
                   + [str(ROOT / 'csd/ic_chain_block_closure.py'), str(ROOT / 'csd/ic_chain_shift_closure.py'),
                      str(ROOT / 'csd/ic_fall_certificates.py'), str(ROOT / 'csd/ic_curve_sweep.py'),
                      str(ROOT / 'csd/ic_dreg.py'), str(ROOT / 'csd/elim.py')])
    data = sorted(p for p in glob.glob(str(OUT / '*.json*')) if not p.endswith('hashes.json'))
    (OUT / 'hashes.json').write_text(json.dumps(dict(
        source_sha256={str(Path(f).relative_to(ROOT)): sha(f) for f in files},
        data_sha256={str(Path(f).relative_to(ROOT)): sha(f) for f in data}), indent=1) + '\n')
    plan = OUT / 'plan.json'
    if plan.exists():
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        pl = json.loads(plan.read_text())
        ks = [r['k'] for r in pl['rows']]
        fig, ax = plt.subplots(figsize=(7, 4.2))
        for D, color in ((5, '#4C78A8'), (6, '#F58518'), (7, '#54A24B')):
            ax.plot(ks, [r[f'D{D}']['peak_gb_model'] for r in pl['rows']], 'o-', color=color,
                    label=f'degree {D}: dense split engine, model peak')
        ax.axhline(15, color='#888', ls='--', lw=1)
        ax.text(ks[0], 17, '15 GB host', color='#555', fontsize=8)
        ax.axhline(512, color='#888', ls=':', lw=1)
        ax.text(ks[0], 580, '512 GB host', color='#555', fontsize=8)
        ax.set_yscale('log')
        ax.set_xlabel('chain depth k')
        ax.set_ylabel('peak memory, GB (model)')
        ax.set_title('Exact degree-D closure cost (two summands); LFD 7 needs a deficient V_6')
        ax.legend(fontsize=8, loc='upper left')
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(OUT / 'cost.svg')
        fig.savefig(OUT / 'cost.png', dpi=130)
    print(json.dumps({k: v for k, v in summary.items() if k != 'table'}, indent=1))


if __name__ == '__main__':
    main()
