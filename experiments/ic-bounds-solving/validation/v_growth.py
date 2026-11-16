#!/usr/bin/env python3
"""Validator check for the growth statements beyond the brute-force range, through my own g^{-1} (Theorem-backed):
  * NTT self-check against the O(N^3) matmul DFT at small N;
  * m = 3: lane B scan rows at k = 3..8 (recompute delta, maxwt, c_top flag, sd = delta + maxwt) and fresh random unit
    instances of my own at p = 65537, k = 3..6;
  * m = 2: lane A rows at k = 11, 12, 13 recomputed (sd = 2 + maxwt, c_top flag, nzeros).
"""
import glob, json, os, random, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vcommon as V

LB = '/home/user/ml-cryptanalysis/experiments/ic-bounds-solving/laneB_prime_m3/data'
LA = '/home/user/ml-cryptanalysis/experiments/ic-bounds-solving/laneA_prime_m2/data/all_instances.json'


def selfcheck():
    rng = np.random.default_rng(1)
    out = []
    for (p, N) in [(97, 32), (257, 64), (12289, 128), (65537, 256)]:
        w = V.root_of_unity(p, N)
        X = rng.integers(0, p, size=(N, N), dtype=np.int64)
        A = V.dft_axes(X, w, p, method='matmul')
        B = V.dft_axes(X, w, p, method='ntt')
        back = V.idft_axes(B, w, p)
        out.append(dict(p=p, N=N, ntt_eq_matmul=bool((A == B).all()), roundtrip=bool((back == X).all())))
    return out


def main():
    rep = dict(selfcheck=selfcheck())
    print('selfcheck', rep['selfcheck'], flush=True)
    rng = random.Random(20261127)
    # ---- m = 3, lane B rows
    files = {3: 'scan_k3_p65537.jsonl', 4: 'scan_k4_p65537.jsonl', 5: 'scan_k5_p65537.jsonl', 6: 'scan_k6_p65537.jsonl',
             7: 'scan_k7_p8390273.jsonl', 8: 'scan_k8_p67109633.jsonl'}
    m3 = []
    for k, f in files.items():
        rows = [json.loads(l) for l in open(os.path.join(LB, f))]
        rows = [r for r in rows if r.get('unit')]
        pick = rng.sample(rows, min(3 if k <= 7 else 2, len(rows)))
        for r in pick:
            t0 = time.time()
            info, g, ginv = V.generator_and_inverse(3, k, r['p'], r['a'], r['b'], r['t'])
            rec = dict(k=k, p=r['p'], a=r['a'], b=r['b'], t=r['t'],
                       laneB=dict(delta=r.get('delta'), maxwt=r.get('maxwt_ginv'), top_zero=r.get('top_coeff_zero'), sd=r.get('predicted_sd')),
                       mine=dict(unit=info['unit'], delta=info.get('delta'), maxwt=info.get('maxwt'), top_zero=info.get('c_top_zero'),
                                 sd=info.get('sd_predicted'), conv=info.get('convolution_ok'), charsum=info.get('c_top_charsum_agrees')),
                       sec=round(time.time() - t0, 1))
            rec['agree'] = (info['unit'] and info['delta'] == r.get('delta') and info['maxwt'] == r.get('maxwt_ginv')
                            and info['c_top_zero'] == r.get('top_coeff_zero') and info['sd_predicted'] == r.get('predicted_sd'))
            rec['sd_eq_3k_plus_6'] = (info.get('sd_predicted') == 3 * k + 6)
            m3.append(rec)
            print('m3 laneB', rec, flush=True)
    # lane B's k = 3 brute-force instance
    bf = json.load(open(os.path.join(LB, 'bruteforce_k3_0.json')))
    info, g, ginv = V.generator_and_inverse(3, 3, bf['p'], bf['a'], bf['b'], bf['t'])
    rep['laneB_bruteforce_k3_instance'] = dict(laneB=dict(delta=bf['delta'], maxwt=bf['maxwt_ginv'], sd_bruteforce=bf['sd_bruteforce']),
                                              mine=dict(delta=info['delta'], maxwt=info['maxwt'], sd=info['sd_predicted']),
                                              agree=(info['delta'] == bf['delta'] and info['maxwt'] == bf['maxwt_ginv'] and info['sd_predicted'] == bf['sd_bruteforce']))
    print('laneB k3 bf instance', rep['laneB_bruteforce_k3_instance'], flush=True)
    # fresh m = 3 unit instances at p = 65537
    fresh = []
    for k in [3, 4, 5, 6]:
        p = 65537
        got = 0
        tries = 0
        while got < 3 and tries < 400:
            tries += 1
            a = rng.randrange(1, p)
            b = rng.randrange(1, p)
            if (4 * a ** 3 + 27 * b * b) % p == 0:
                continue
            t = rng.randrange(p)
            info, g, ginv = V.generator_and_inverse(3, k, p, a, b, t)
            if not info['unit']:
                continue
            got += 1
            fresh.append(dict(k=k, p=p, a=a, b=b, t=t, delta=info['delta'], maxwt=info['maxwt'], sd=info['sd_predicted'],
                              top_zero=info['c_top_zero'], conv=info['convolution_ok'], charsum=info['c_top_charsum_agrees'],
                              sd_eq_3k_plus_6=(info['sd_predicted'] == 3 * k + 6), on_locus=((5 * t ** 3 + 3 * a * t + 6 * b) % p == 0)))
            print('m3 fresh', fresh[-1], flush=True)
    rep['m3_laneB_rows'] = m3
    rep['m3_fresh'] = fresh
    # ---- m = 2, lane A rows at k = 11, 12, 13
    rows = json.load(open(LA))
    m2 = []
    for k in [11, 12, 13]:
        cand = [r for r in rows if r['k'] == k]
        for r in rng.sample(cand, min(2, len(cand))):
            t0 = time.time()
            info, g, ginv = V.generator_and_inverse(2, k, r['p'], r['a'], r['b'], r['t'])
            rec = dict(k=k, p=r['p'], a=r['a'], b=r['b'], t=r['t'],
                       laneA=dict(sd=r['sd_inverse'], maxwt=r['maxwt'], top_zero=r['zero_top_coeff'], nzeros=r['nzeros']),
                       mine=dict(unit=info['unit'], delta=info.get('delta'), maxwt=info.get('maxwt'), sd=info.get('sd_predicted'),
                                 top_zero=info.get('c_top_zero'), nzeros=info.get('nzeros_ginv'), conv=info.get('convolution_ok')),
                       sec=round(time.time() - t0, 1))
            rec['agree'] = (info['unit'] and info['maxwt'] == r['maxwt'] and info['sd_predicted'] == r['sd_inverse']
                            and info['c_top_zero'] == r['zero_top_coeff'] and info['nzeros_ginv'] == r['nzeros'] and info['convolution_ok'])
            m2.append(rec)
            print('m2 laneA', rec, flush=True)
    rep['m2_laneA_rows_k11_13'] = m2
    rep['summary'] = dict(m3_laneB_agree=sum(r['agree'] for r in m3), m3_laneB_n=len(m3),
                          m3_fresh_n=len(fresh), m3_fresh_sd_3k6=sum(r['sd_eq_3k_plus_6'] for r in fresh),
                          m2_k11_13_agree=sum(r['agree'] for r in m2), m2_k11_13_n=len(m2))
    json.dump(rep, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out', 'growth.json'), 'w'), indent=1)
    print(json.dumps(rep['summary']))


if __name__ == '__main__':
    main()
