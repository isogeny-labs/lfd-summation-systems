#!/usr/bin/env python3
"""g^{-1} in F_p[(Z/2^k)^3]: evaluate on H^3 (NTT), invert pointwise, inverse NTT; support statistics."""
import sys, os, time, json, random
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from build_m3 import *

_WT = {}
def weight_array(k):
    if k not in _WT:
        N = 1 << k
        pc = popcount_table(N).astype(np.int8)
        _WT[k] = pc[:, None, None] + pc[None, :, None] + pc[None, None, :]
    return _WT[k]

def analyse(p, a, b, t, k, check_conv=False):
    t0 = time.time()
    N = 1 << k
    g, delta = build_g(p, a, b, t, k)
    w = root_N(p, N)
    v = dft3(g, w, p)
    nzero = int((v == 0).sum())
    rec = dict(p=p, k=k, a=a, b=b, t=t, delta=delta, n_zeros_on_H3=nzero, unit=(nzero == 0),
               locus5=((5 * t ** 3 + 3 * a * t + 6 * b) % p == 0))
    if nzero:
        rec['seconds'] = round(time.time() - t0, 3)
        return rec
    gi = idft3(modinv_arr(v, p), w, p)
    if check_conv:
        prod = np.zeros_like(g)
        for idx in np.argwhere(g):
            e = tuple(int(x) for x in idx)
            prod = (prod + g[e] * np.roll(gi, e, axis=(0, 1, 2))) % p
        one = np.zeros_like(g); one[0, 0, 0] = 1
        rec['inverse_verified_by_convolution'] = bool((prod == one).all())
    W = weight_array(k)
    nzmask = gi != 0
    maxwt = int(W[nzmask].max())
    zeros_per_w = {}
    for wt in range(3 * k + 1):
        sel = W == wt
        zeros_per_w[wt] = [int((~nzmask[sel]).sum()), int(sel.sum())]
    rec.update(maxwt_ginv=maxwt, top_coeff_zero=bool(gi[N - 1, N - 1, N - 1] == 0),
               top_coeff=int(gi[N - 1, N - 1, N - 1]), zeros_per_weight=zeros_per_w,
               predicted_sd=delta + maxwt, generic_sd=delta + 3 * k, seconds=round(time.time() - t0, 3))
    return rec

if __name__ == '__main__':
    p, k, a, b, t = map(int, sys.argv[1:6])
    print(json.dumps(analyse(p, a, b, t, k, check_conv=(k <= 3)), indent=1))
