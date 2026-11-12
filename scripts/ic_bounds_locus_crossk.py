"""Compare the exceptional classes at depths k = 6 and k = 7 on the same normalized
points (A, B) over p = 257 (bounds item 2).  For every k = 7 deviating target of the
main sweep and every non-generic point of the p = 257, k = 6 line sweep, compute
the dense V_3, V_4 (modulo chains, no zero set) at both depths for (A, B, t = 1).

Usage: python3 scripts/ic_bounds_locus_crossk.py --out FILE.json
"""
import argparse, glob, json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'scripts'))
from ic_bounds_locus_plane import s3_generator, closure_row  # noqa: E402
from ic_bounds_locus_classify import load, summarize, deviations  # noqa: E402

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out', required=True); args = ap.parse_args()
    p = 257
    pts = {}
    for r in load(p, 7, ROOT / 'experiments/ic-bounds-locus/data'):
        s = summarize(r)
        if deviations(s, 7) and r['t']:
            it = pow(r['t'], -1, p)
            pts[(r['a'] * it * it % p, r['b'] * it ** 3 % p)] = dict(source='k7-main', a=r['a'], b=r['b'], t=r['t'])
    for f in glob.glob(str(ROOT / 'experiments/ic-bounds-locus/plane/p257-k6-lines-*.jsonl')):
        for line in open(f):
            r = json.loads(line)
            if (r['V3'], r['V4'], r['low4']) != (13, 605, 111):
                pts.setdefault((r['A'], r['B']), dict(source='k6-lines'))
    out = []
    for (A, B), src in sorted(pts.items()):
        rec = dict(A=A, B=B, **src)
        for k in (6, 7):
            row = closure_row(p, k, s3_generator(p, k, A, B, 1), certify=False, with_zeros=True)
            rec[f'k{k}'] = row
        out.append(rec)
        print(A, B, src['source'], [(rec[f'k{k}']['V3'], rec[f'k{k}']['V4'], rec[f'k{k}']['low4'], rec[f'k{k}']['Z'], rec[f'k{k}']['saturated4']) for k in (6, 7)], flush=True)
    Path(args.out).write_text(json.dumps(out, indent=1) + '\n')

if __name__ == '__main__':
    main()
