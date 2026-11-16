"""Driver for tasks 1-4.  For (k, p): draw random (a,b,t), get sd via g^{-1} (sd_inverse; None = g has a zero
on mu_N x mu_N, i.e. non-unit), optionally brute force sd (k <= 6), optionally d_lf/d_ff.
Usage: python3 run_scan.py k p n_unit max_tries brute(0/1) dlf(none|ord|quot) n_dlf seed
Writes data/scan_k{k}_p{p}_s{seed}.json.  [pcc:unverified]"""
import sys, os, json, time, random
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import sd_inverse as SI, static_bruteforce as SB, lf_quotient as Q

k, p, n_unit, max_tries, brute, dlf, n_dlf, seed = sys.argv[1:9]
k, p, n_unit, max_tries, brute, n_dlf, seed = int(k), int(p), int(n_unit), int(max_tries), int(brute), int(n_dlf), int(seed)
rng = random.Random(seed * 7919 + p * 31 + k)
rows, tries, t_start = [], 0, time.time()
t_inv = t_bf = t_dlf = 0.0
while len(rows) < n_unit and tries < max_tries:
    tries += 1
    a, b = random_curve(p, rng); t = rng.randrange(p)
    t0 = time.time(); r = SI.sd_inverse(a, b, t, p, k); t_inv += time.time() - t0
    if r is None: continue
    row = dict(p=p, a=a, b=b, t=t, k=k, sd_inverse=r['sd_inverse'], maxwt=r['maxwt'], zero_top_coeff=r['zero_top_coeff'],
               zeros_per_weight=r['zeros_per_weight'], nzeros=r['nzeros'], inverse_verified=r['inverse_verified'],
               sd_bruteforce=None, d_lf=None, d_ff=None, d_lf_source=None)
    if brute:
        t0 = time.time(); sd, dims = SB.sd_nf(a, b, t, p, k); t_bf += time.time() - t0
        row['sd_bruteforce'] = sd
        row['P1_match'] = (sd == r['sd_inverse'])
    if dlf != 'none' and len(rows) < n_dlf:
        t0 = time.time()
        if dlf == 'ord':
            o = SB.lfd_per_degree(a, b, t, p, k, Dmax=k + 3)
            row.update(d_lf=o['d_lf'], d_ff=o['d_ff'], d_lf_source='lfd_engine (ordinary ring, exact)')
        else:
            q = Q.profile(a, b, t, p, k, Dmax=k + 3)
            row.update(d_lf=q['d_lf'], d_ff=q['d_ff'], d_lf_source='quotient closure (lf_quotient.py, exact)')
        t_dlf += time.time() - t0
    rows.append(row)
out = dict(k=k, p=p, seed=seed, tries=tries, n_unit=len(rows), unit_rate=len(rows) / max(1, tries),
           wall_sec=round(time.time() - t_start, 1), sec_inverse=round(t_inv, 1), sec_bruteforce=round(t_bf, 1), sec_dlf=round(t_dlf, 1),
           rows=rows)
json.dump(out, open(os.path.join(os.path.dirname(__file__), '..', 'data', 'scan_k%d_p%d_s%d.json' % (k, p, seed)), 'w'), indent=1)
sds = [r['sd_inverse'] for r in rows]
print('k=%d p=%d unit %d/%d  sd=2k+2: %d  other: %s  wall %.1fs  (inv %.1f bf %.1f dlf %.1f)' % (
    k, p, len(rows), tries, sum(s == 2 * k + 2 for s in sds), sorted(set(s for s in sds if s != 2 * k + 2)), time.time() - t_start, t_inv, t_bf, t_dlf))
bad = [r for r in rows if r.get('P1_match') is False]
print('P1 mismatches:', len(bad), 'of', sum(1 for r in rows if r['sd_bruteforce'] is not None))
