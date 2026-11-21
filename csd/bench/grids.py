"""Named parameter grids.  Each returns a list of Spec; seeds are 0 .. seeds-1."""
from __future__ import annotations

from .spec import Spec


def _seeds(seeds, f):
    return [f(s) for s in range(seeds)]


def smoke(seeds=2):
    """Seconds: every solver that applies, tiny sizes, planted targets (non-empty cross-check)."""
    out = []
    for k in (4, 5, 6):
        out += _seeds(seeds, lambda s: Spec(n=3 * k - 2, k=k, m=3, seed=s, plant=True))
    out += _seeds(seeds, lambda s: Spec(n=13, k=4, m=4, seed=s, plant=True))
    return out


def m3(seeds=3, kmin=5, kmax=15, plant=False):
    """m = 3, V = <1, w, ..>, n = 3k - 2: the scaling table of experiments/ic-rootfind/."""
    return [Spec(n=3 * k - 2, k=k, m=3, seed=s, plant=plant) for k in range(kmin, kmax + 1) for s in range(seeds)]


def m4(seeds=3, kmin=4, kmax=9, plant=False):
    return [Spec(n=4 * k - 3, k=k, m=4, seed=s, plant=plant) for k in range(kmin, kmax + 1) for s in range(seeds)]


def basis(seeds=3, n=31, k=10):
    """The factor base decides: <1..w^(j-1)> plus k - j random vectors."""
    out = []
    for j in (k, k - 1, k - 2, 0):
        out += [Spec(n=n, k=k, m=3, seed=s, basis="standard" if j == k else "mixed", j=None if j == k else j)
                for s in range(seeds)]
    return out


def regime(seeds=3, k=12, nmin=27, nmax=45):
    """Fixed k, n from below 3k (free unknowns) to above it (no solvable systems)."""
    return [Spec(n=n, k=k, m=3, seed=s) for n in range(nmin, nmax + 1, 2) for s in range(seeds)]


def calib(seeds=3):
    """Instances of experiments/ic-rootfind/scaling-m3.jsonl that enumeration can still check (k <= 10)."""
    return m3(seeds=seeds, kmin=7, kmax=10)


GRIDS = dict(smoke=smoke, m3=m3, m4=m4, basis=basis, regime=regime, calib=calib)
