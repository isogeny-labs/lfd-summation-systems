#!/usr/bin/env python3
"""data/sd_vs_k_m3.svg : sd (via delta+maxwt g^-1; brute force at k=2,3 marked), d_lf, 3k+1 (K1), delta+3k."""
import json, os, glob
D = os.path.join(os.path.dirname(__file__), '..', 'data')
S = json.load(open(D + '/summary.json'))
rows = []
for f in glob.glob(D + '/scan_k*.jsonl'):
    for l in open(f):
        r = json.loads(l)
        if r['unit'] and not r['locus5']: rows.append(r)
ks = sorted({r['k'] for r in rows})
sd_by_k = {k: (min(r['predicted_sd'] for r in rows if r['k'] == k), max(r['predicted_sd'] for r in rows if r['k'] == k)) for k in ks}
bf = {}
for b in S['bruteforce']:
    if not b['locus5']: bf.setdefault(b['k'], []).append(b['sd_bruteforce'])
for f in glob.glob(D + '/bruteforce_k3_*.json'):
    b = json.load(open(f))
    if b.get('sd_bruteforce') and not b['locus5']: bf.setdefault(3, []).append(b['sd_bruteforce'])
W, H, L, R, T, B = 760, 480, 70, 30, 40, 60
xs = lambda k: L + (k - 1.5) / (8.5 - 1.5) * (W - L - R)
ys = lambda v: H - B - v / 32 * (H - T - B)
o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" font-family="sans-serif" font-size="12">',
     f'<rect width="{W}" height="{H}" fill="#ffffff"/>',
     f'<text x="{W/2}" y="22" text-anchor="middle" font-size="14">m = 3 prime chains: static solving degree sd vs depth k  [pcc:unverified]</text>']
for v in range(0, 33, 4):
    o.append(f'<line x1="{L}" x2="{W-R}" y1="{ys(v)}" y2="{ys(v)}" stroke="#ddd"/><text x="{L-8}" y="{ys(v)+4}" text-anchor="end">{v}</text>')
for k in range(2, 9):
    o.append(f'<line x1="{xs(k)}" x2="{xs(k)}" y1="{T}" y2="{H-B}" stroke="#eee"/><text x="{xs(k)}" y="{H-B+18}" text-anchor="middle">{k}</text>')
o.append(f'<text x="{W/2}" y="{H-12}" text-anchor="middle">depth k</text><text x="18" y="{H/2}" transform="rotate(-90 18 {H/2})" text-anchor="middle">degree</text>')
def line(pts, col, dash=''):
    d = ' '.join(f'{xs(k):.1f},{ys(v):.1f}' for k, v in pts)
    return f'<polyline points="{d}" fill="none" stroke="{col}" stroke-width="2" {dash}/>'
K = list(range(2, 9))
o.append(line([(k, 3 * k + 6) for k in K], '#1b6ca8', 'stroke-dasharray="6 4"'))
o.append(line([(k, 3 * k + 1) for k in K], '#c0392b', 'stroke-dasharray="2 4"'))
o.append(line([(k, 7) for k in range(2, 6)], '#2e8b57'))
o.append(line([(k, sd_by_k[k][1]) for k in ks], '#222'))
for k in ks:
    o.append(f'<circle cx="{xs(k)}" cy="{ys(sd_by_k[k][1])}" r="4" fill="#222"/>')
    if sd_by_k[k][0] != sd_by_k[k][1]:
        o.append(f'<path d="M{xs(k)},{ys(sd_by_k[k][0])-5} l5,5 l-5,5 l-5,-5 z" fill="#fff" stroke="#222" stroke-width="1.5"/>')
    for v in set(bf.get(k, [])):
        o.append(f'<rect x="{xs(k)-7}" y="{ys(v)-7}" width="14" height="14" fill="none" stroke="#e08a00" stroke-width="2"/>')
for k in range(2, 6):
    o.append(f'<circle cx="{xs(k)}" cy="{ys(7)}" r="3" fill="#2e8b57"/>')
o.append(f'<path d="M{xs(6)},{ys(8)-6} l6,10 l-12,0 z" fill="#2e8b57"/><text x="{xs(6)+9}" y="{ys(8)+3}" fill="#2e8b57">&gt;= 8 (one non-unit instance)</text>')
o.append(f'<polygon points="{xs(2)},{ys(7)-9} {xs(2)-7},{ys(7)+5} {xs(2)+7},{ys(7)+5}" fill="none" stroke="#2e8b57" stroke-width="2"/>')
lg = [('#222', 'sd = delta + maxwt(g^-1) (max over scan = 3k+6; open diamond = exceptions, 3k+5)'), ('#e08a00', 'sd by ordinary-ring brute force (k=2, k=3 if run)'),
      ('#1b6ca8', 'delta + 3k  (delta = 6)'), ('#c0392b', '3k + 1  (K1 kernel bound, mk+1)'), ('#2e8b57', 'd_lf = 7 (recorded k=2..5; measured by me at k=2,3; k=6 only >= 8)')]
for i, (c, t) in enumerate(lg):
    y = T + 18 + 17 * i
    o.append(f'<rect x="{L+12}" y="{y-9}" width="12" height="12" fill="{c}"/><text x="{L+30}" y="{y+1}">{t}</text>')
o.append('</svg>')
open(D + '/sd_vs_k_m3.svg', 'w').write('\n'.join(o))
print('sd per k', sd_by_k, 'bruteforce', bf)
