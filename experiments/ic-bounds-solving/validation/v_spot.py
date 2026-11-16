#!/usr/bin/env python3
"""Validator check 4: spot-check lane A's data/all_instances.json.

(a) 30 stratified random rows (k = 3..10) recomputed with my own g^{-1} (unit flag, maxwt, sd = 2 + maxwt,
    top coefficient zero flag, number of zero coefficients, convolution check); rows with k <= 4 also get my
    own ordinary-ring static brute force (v_static.static_sd) compared with the row's sd_bruteforce / sd_inverse.
(b) aggregate recounts from the rows: total rows, exceptions (sd != 2k+2), zero-top-coefficient rows, their
    coincidence, the k = 3, p = 73 exception fraction, duplicates, brute-force rows and mismatches, and the
    per-(k, p) P2 test (generic fraction >= 1 - 3/p, falsifier < 1 - 10/p).
"""
import argparse, collections, json, os, random, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vcommon import generator_and_inverse
from v_static import static_sd

STRATA = {3: 6, 4: 6, 5: 4, 6: 3, 7: 3, 8: 3, 9: 2, 10: 3}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--laneA', default='/home/user/ml-cryptanalysis/experiments/ic-bounds-solving/laneA_prime_m2/data/all_instances.json')
    ap.add_argument('--seed', type=int, default=20261126)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    rows = json.load(open(a.laneA))
    rng = random.Random(a.seed)
    report = {}
    # ---- (b) aggregates
    agg = {}
    agg['rows'] = len(rows)
    exc = [r for r in rows if r['sd_inverse'] != 2 * r['k'] + 2]
    agg['exceptions'] = len(exc)
    agg['exceptions_by_k_p'] = {f"{k},{p}": c for (k, p), c in sorted(collections.Counter((r['k'], r['p']) for r in exc).items())}
    agg['exceptions_all_sd_2k_plus_1'] = all(r['sd_inverse'] == 2 * r['k'] + 1 for r in exc)
    agg['zero_top_coeff_rows'] = sum(1 for r in rows if r['zero_top_coeff'])
    agg['zero_top_coeff_iff_exception'] = all((r['zero_top_coeff']) == (r['sd_inverse'] != 2 * r['k'] + 2) for r in rows)
    agg['sd_inverse_eq_2_plus_maxwt_all'] = all(r['sd_inverse'] == 2 + r['maxwt'] for r in rows)
    k3p73 = [r for r in rows if r['k'] == 3 and r['p'] == 73]
    e73 = [r for r in k3p73 if r['sd_inverse'] != 8]
    agg['k3_p73'] = dict(rows=len(k3p73), exceptions=len(e73), fraction=len(e73) / len(k3p73), one_over_p=1 / 73)
    bf = [r for r in rows if r.get('sd_bruteforce') is not None]
    agg['bruteforce_rows'] = len(bf)
    agg['bruteforce_by_k'] = dict(sorted(collections.Counter(r['k'] for r in bf).items()))
    agg['bruteforce_mismatches'] = sum(1 for r in bf if r['sd_bruteforce'] != r['sd_inverse'])
    agg['bruteforce_exceptions_checked'] = sum(1 for r in bf if r['sd_inverse'] != 2 * r['k'] + 2)
    agg['inverse_verified_false'] = sum(1 for r in rows if not r.get('inverse_verified'))
    agg['duplicate_instance_keys'] = len(rows) - len({(r['p'], r['a'], r['b'], r['t'], r['k']) for r in rows})
    agg['rows_by_k'] = dict(sorted(collections.Counter(r['k'] for r in rows).items()))
    agg['nzeros_eq_zero_implies_generic'] = all((r['nzeros'] > 0) or (not r['zero_top_coeff']) for r in rows)
    cells = collections.defaultdict(lambda: [0, 0])
    for r in rows:
        c = cells[(r['k'], r['p'])]
        c[0] += 1
        c[1] += (r['sd_inverse'] == 2 * r['k'] + 2)
    p2 = []
    for (k, p), (n, gen) in sorted(cells.items()):
        frac = gen / n
        p2.append(dict(k=k, p=p, n=n, generic=gen, fraction=round(frac, 5), bound_1_minus_3_over_p=round(1 - 3 / p, 5),
                       falsifier_1_minus_10_over_p=round(1 - 10 / p, 5), passes_1_minus_3_over_p=frac >= 1 - 3 / p,
                       below_falsifier=frac < 1 - 10 / p))
    agg['p2_cells'] = p2
    agg['p2_cells_failing_1_minus_3_over_p'] = [c for c in p2 if not c['passes_1_minus_3_over_p']]
    agg['p2_cells_below_falsifier'] = [c for c in p2 if c['below_falsifier']]
    report['aggregates'] = agg
    # ---- (a) 30 rows
    picks = []
    for k, cnt in STRATA.items():
        cand = [r for r in rows if r['k'] == k]
        # ensure at least one exception is included where exceptions exist at that k (k = 3, 4)
        ex = [r for r in cand if r['sd_inverse'] != 2 * k + 2]
        chosen = []
        if ex:
            chosen.append(rng.choice(ex))
        rest = [r for r in cand if r not in chosen]
        chosen += rng.sample(rest, cnt - len(chosen))
        picks += chosen
    checks = []
    for r in picks:
        t0 = time.time()
        info, g, ginv = generator_and_inverse(2, r['k'], r['p'], r['a'], r['b'], r['t'])
        c = dict(k=r['k'], p=r['p'], a=r['a'], b=r['b'], t=r['t'], laneA=dict(sd_inverse=r['sd_inverse'], maxwt=r['maxwt'],
                 zero_top_coeff=r['zero_top_coeff'], nzeros=r['nzeros'], sd_bruteforce=r.get('sd_bruteforce')),
                 mine=dict(unit=info['unit'], delta=info.get('delta'), maxwt=info.get('maxwt'), sd_predicted=info.get('sd_predicted'),
                           c_top_zero=info.get('c_top_zero'), nzeros=info.get('nzeros_ginv'), convolution_ok=info.get('convolution_ok'),
                           c_top_charsum_agrees=info.get('c_top_charsum_agrees')))
        ok = (info['unit'] and info['maxwt'] == r['maxwt'] and info['sd_predicted'] == r['sd_inverse']
              and info['c_top_zero'] == r['zero_top_coeff'] and info['nzeros_ginv'] == r['nzeros'] and info['convolution_ok'])
        c['agree_inverse'] = bool(ok)
        if r['k'] <= 4:
            sd, prof = static_sd(2, r['k'], r['p'], g, info['delta'])
            c['mine']['sd_static'] = sd
            c['agree_static_vs_laneA_bruteforce'] = (r.get('sd_bruteforce') is None) or (sd == r['sd_bruteforce'])
            c['agree_static_vs_inverse'] = (sd == r['sd_inverse'])
        c['seconds'] = round(time.time() - t0, 2)
        checks.append(c)
        print(f"k={r['k']} p={r['p']} a={r['a']} b={r['b']} t={r['t']}: laneA sd={r['sd_inverse']} maxwt={r['maxwt']} top0={r['zero_top_coeff']} nz={r['nzeros']} bf={r.get('sd_bruteforce')} | "
              f"mine maxwt={info.get('maxwt')} sd={info.get('sd_predicted')} top0={info.get('c_top_zero')} nz={info.get('nzeros_ginv')} "
              f"static={c['mine'].get('sd_static')} agree={c['agree_inverse']} {c['seconds']}s", flush=True)
    report['spot_rows'] = checks
    report['spot_summary'] = dict(n=len(checks), agree_inverse=sum(c['agree_inverse'] for c in checks),
                                  static_checked=sum(1 for c in checks if 'sd_static' in c['mine']),
                                  static_agree_inverse=sum(1 for c in checks if c.get('agree_static_vs_inverse')),
                                  static_agree_laneA_bf=sum(1 for c in checks if c.get('agree_static_vs_laneA_bruteforce')))
    json.dump(report, open(a.out, 'w'), indent=1)
    print(json.dumps(report['spot_summary']))
    print('aggregates:', json.dumps({k: v for k, v in agg.items() if k not in ('p2_cells',)}, indent=None)[:3000])


if __name__ == '__main__':
    main()
