"""Statistics for benchmark records: medians and IQR, bootstrap intervals, paired ratios, log2-linear fits, noise floor.

Timings have heavy right tails, so everything is a median or a geometric mean; intervals come from a seeded
percentile bootstrap over instances (the unit of replication), never over repetitions of one instance.
"""
from __future__ import annotations

import math

import numpy as np


def summ(xs):
    a = np.asarray(sorted(xs), dtype=float)
    if a.size == 0:
        return dict(n=0)
    return dict(n=int(a.size), median=float(np.median(a)), q1=float(np.percentile(a, 25)),
                q3=float(np.percentile(a, 75)), min=float(a[0]), max=float(a[-1]))


def boot_ci(xs, fn=np.median, B=2000, alpha=0.05, seed=0):
    a = np.asarray(xs, dtype=float)
    if a.size == 0:
        return (math.nan, math.nan)
    if a.size == 1:
        return (float(a[0]), float(a[0]))
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, a.size, size=(B, a.size))
    vals = np.array([fn(a[i]) for i in idx])
    return (float(np.percentile(vals, 100 * alpha / 2)), float(np.percentile(vals, 100 * (1 - alpha / 2))))


def geomean(xs):
    a = np.asarray(xs, dtype=float)
    return float(np.exp(np.mean(np.log(a)))) if a.size else math.nan


def paired_ratio(num: dict, den: dict, B=2000, seed=0):
    """Geometric mean of num[key] / den[key] over the common keys, with a bootstrap interval over keys."""
    keys = sorted(set(num) & set(den))
    r = np.array([num[k] / den[k] for k in keys], dtype=float)
    if r.size == 0:
        return dict(n=0)
    lo, hi = boot_ci(np.log(r), np.mean, B=B, seed=seed)
    return dict(n=int(r.size), gmean=float(np.exp(np.log(r).mean())), lo=float(np.exp(lo)), hi=float(np.exp(hi)),
                min=float(r.min()), max=float(r.max()))


def fit_log2(xs, ys, B=2000, seed=0):
    """y ~ c * 2^(s x): OLS of log2 y on x.  Returns slope s (doublings per unit of x) with a bootstrap interval,
    intercept, R^2 and the number of points.  A slope of 1 means y doubles per step of x."""
    x = np.asarray(xs, dtype=float)
    y = np.log2(np.asarray(ys, dtype=float))
    if x.size < 3 or np.ptp(x) == 0:
        return dict(n=int(x.size))

    def ols(ix):
        xx, yy = x[ix], y[ix]
        if np.ptp(xx) == 0:
            return math.nan, math.nan
        s, c = np.polyfit(xx, yy, 1)
        return s, c

    s, c = ols(np.arange(x.size))
    pred = s * x + c
    ss_res, ss_tot = float(((y - pred) ** 2).sum()), float(((y - y.mean()) ** 2).sum())
    rng = np.random.default_rng(seed)
    slopes = [ols(rng.integers(0, x.size, x.size))[0] for _ in range(B)]
    slopes = [v for v in slopes if not math.isnan(v)]
    return dict(n=int(x.size), slope=float(s), lo=float(np.percentile(slopes, 2.5)), hi=float(np.percentile(slopes, 97.5)),
                intercept=float(c), r2=(1 - ss_res / ss_tot) if ss_tot else math.nan)


def noise_floor(pairs):
    """pairs: list of (t_A, t_B) of the same solver on the same instance run as two arms.  Returns the spread of
    B/A: its median, 5th and 95th percentile.  Differences between solvers smaller than this are not resolved."""
    r = np.array([b / a for a, b in pairs if a > 0 and b > 0], dtype=float)
    if r.size == 0:
        return dict(n=0)
    return dict(n=int(r.size), median=float(np.median(r)), p05=float(np.percentile(r, 5)), p95=float(np.percentile(r, 95)))
