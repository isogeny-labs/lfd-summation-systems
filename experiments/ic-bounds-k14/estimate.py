"""ESTIMATE (not a measurement) of the multiply-adds the big engine performs, extrapolated from the local
runs in validation/rates-*.json.

Model: M(D, c) = M(D, c_last) * (c / c_last)^beta with c the number of columns; beta is the local log-log
slope between consecutive measured depths k.  The band uses the smallest and the largest of the last three
local slopes.  A backtest predicts the measured D = 5 points k = 10, 11 from k <= 9 only.
Assumptions that matter: same instance family (unit curve, p = 12289, 1 in V_D for k <= 9 so the
rank is full there, which is NOT the regime of the target run), one thread, closure phase only (the
certificate sweep costs about as much as the closure, see rates), engine defaults.  Rates for the GPU time
are assumptions, not measurements.
"""
import json
import math
import sys
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
rows = []
for f in sorted((HERE / 'validation').glob('rates-*.json')):
    rows += json.load(open(f))['rows']
data = {}
for r in rows:
    data.setdefault(r['D'], {})[r['k']] = (r['columns'], r['madds_closure_phase'], r['contains_one'])


def cols(k, D):
    return sum(comb(2 * k, w) for w in range(D + 1))


def slopes(D, upto=None):
    ks = sorted(k for k in data[D] if upto is None or k <= upto)
    out = []
    for a, b in zip(ks, ks[1:]):
        (c0, m0, _), (c1, m1, _) = data[D][a], data[D][b]
        out.append((b, math.log(m1 / m0) / math.log(c1 / c0)))
    return out


def extrapolate(D, k, upto=None):
    ks = sorted(kk for kk in data[D] if upto is None or kk <= upto)
    kl = ks[-1]
    c_last, m_last, _ = data[D][kl]
    sl = [s for _, s in slopes(D, upto)][-3:]
    c = cols(k, D)
    return [m_last * (c / c_last) ** b for b in (min(sl), sum(sl) / len(sl), max(sl))], sl


out = dict(measured={D: {k: dict(columns=v[0], madds_closure=v[1], contains_one=v[2]) for k, v in sorted(d.items())}
                     for D, d in data.items()}, slopes={D: slopes(D) for D in data})
back = {}
for k in sorted(data[5]):
    if k >= 10:
        (lo, mid, hi), sl = extrapolate(5, k, upto=9)
        back[k] = dict(predicted_low_mid_high=[lo, mid, hi], measured=data[5][k][1], ratio_measured_over_mid=data[5][k][1] / mid)
out['backtest_D5_from_k<=9'] = back
rate_lo, rate_hi = 6e12, 12e12      # ASSUMPTION: sustained exact-FP64 multiply-adds/s on one B200 (about 30-60% of the 20 T peak)
# certificate sweep: measured total / closure-phase madds on the non-full-rank runs (D = 5, k = 10, 11)
sweep_factor = [r['madds_total'] / r['madds_closure_phase'] for r in rows if r['rank'] < r['columns']]
out['sweep_factor_measured'] = sweep_factor
sf = sum(sweep_factor) / len(sweep_factor) if sweep_factor else 1.0
nonfull = [r for r in rows if r['rank'] < r['columns']]
blended = sum(r['madds_total'] for r in nonfull) / sum(r['seconds'] for r in nonfull)   # madds/s, one thread, all kernels, closure + sweep
out['blended_single_thread_madds_per_s'] = blended
tab = []
for D in (5, 6):
    for k in range(max(data[D]) + 1, 15):
        (lo, mid, hi), sl = extrapolate(D, k)
        tab.append(dict(D=D, k=k, columns=cols(k, D), madds_low=lo, madds_mid=mid, madds_high=hi, slopes_used=sl,
                        madds_with_sweep_mid=mid * sf,
                        gemm_bound_hours_low=lo * sf / rate_hi / 3600, gemm_bound_hours_high=hi * sf / rate_lo / 3600,
                        single_thread_days_at_blended_rate=mid * sf / blended / 86400))
out['estimate'] = tab
out['gpu_rate_assumption_madds_per_s'] = [rate_lo, rate_hi]
(HERE / 'validation' / 'estimate.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps(out['backtest_D5_from_k<=9'], indent=1))
for t in tab:
    print(f"D={t['D']} k={t['k']} c={t['columns']:>7} madds {t['madds_low']:.2e} .. {t['madds_mid']:.2e} .. {t['madds_high']:.2e}  "
          f"single-thread days {t['single_thread_days_at_blended_rate']:.2f} with sweep (mid) {t['madds_with_sweep_mid']:.2e}  gemm-bound h {t['gemm_bound_hours_low']:.3f} .. {t['gemm_bound_hours_high']:.3f}")
