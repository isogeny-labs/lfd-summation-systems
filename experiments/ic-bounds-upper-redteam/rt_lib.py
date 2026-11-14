"""Red-team builder for two-summand prime chain systems in the ORDINARY polynomial ring.
Own code; only the closure / exact-LFD engine is imported from the independent validator
(experiments/ic-curve-validation/code/lfd_engine.py).  No import of experiments/ic-bounds-upper."""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'ic-curve-validation', 'code'))
import lfd_engine as E


def s3_coeffs(p, a, b, t, s1=1, s2=1):
    """Coefficients of S_3(s1*x1, s2*x2, t) as {(e1, e2): c}, written out from
    S_3 = (x1-x2)^2 t^2 - 2((x1+x2)(x1 x2 + a) + 2b) t + (x1 x2 - a)^2 - 4b(x1+x2)."""
    base = {
        (2, 2): 1,
        (2, 1): -2 * t, (1, 2): -2 * t,
        (2, 0): t * t, (0, 2): t * t,
        (1, 1): -2 * t * t - 2 * a,
        (1, 0): -2 * a * t - 4 * b, (0, 1): -2 * a * t - 4 * b,
        (0, 0): a * a - 4 * b * t,
    }
    return {(e1, e2): c * pow(s1, e1, p) * pow(s2, e2, p) % p for (e1, e2), c in base.items()}


def s3_eval(p, a, b, t, x1, x2):
    return ((x1 - x2) ** 2 * t * t - 2 * ((x1 + x2) * (x1 * x2 + a) + 2 * b) * t
            + (x1 * x2 - a) ** 2 - 4 * b * (x1 + x2)) % p


def mu(p, N):
    return [x for x in range(1, p) if pow(x, N, p) == 1]


def build(p, a, b, t, k1, k2, s1=1, s2=1, raw=False):
    """Returns (inputs, points, r).  Variables: u_{1,j} = index j, u_{2,j} = index k1 + j."""
    r = k1 + k2
    ks, offs = (k1, k2), (0, k1)
    inputs = []
    for i in range(2):
        for j in range(ks[i]):
            sq = [0] * r
            sq[offs[i] + j] = 2
            f = {tuple(sq): 1}
            if j + 1 < ks[i]:
                nx = [0] * r
                nx[offs[i] + j + 1] = 1
                f[tuple(nx)] = p - 1
            else:
                f[tuple([0] * r)] = p - 1
            inputs.append(f)
    G = {}
    for (e1, e2), c in s3_coeffs(p, a, b, t, s1, s2).items():
        m = [0] * r
        for i, e in ((0, e1), (1, e2)):
            if e == 1:
                m[offs[i]] += 1
            elif e == 2:
                if raw:
                    m[offs[i]] += 2
                elif ks[i] >= 2:
                    m[offs[i] + 1] += 1
                # ks[i] == 1: x^2 = 1, nothing to add
        m = tuple(m)
        G[m] = (G.get(m, 0) + c) % p
    G = {m: c for m, c in G.items() if c}
    inputs.append(G)
    H1, H2 = mu(p, 1 << k1), mu(p, 1 << k2)
    assert len(H1) == 1 << k1 and len(H2) == 1 << k2
    pts = []
    for x1 in H1:
        for x2 in H2:
            if s3_eval(p, a, b, t, s1 * x1 % p, s2 * x2 % p) == 0:
                pts.append(tuple([pow(x1, 1 << j, p) for j in range(k1)] + [pow(x2, 1 << j, p) for j in range(k2)]))
    return inputs, pts, r


def exact(p, a, b, t, k1, k2, s1=1, s2=1, raw=False, Dmax=None):
    inputs, pts, r = build(p, a, b, t, k1, k2, s1, s2, raw)
    bound = min(k1, k2) + (4 if raw else 2)
    if Dmax is None:
        Dmax = bound + 1
    res = E.exact_lfd(inputs, r, p, pts, Dmax)
    return res, bound, len(pts)


def nonsingular(p, a, b):
    return (4 * a ** 3 + 27 * b * b) % p != 0


def coset_reps(p, k):
    """One representative per coset of mu_{2^k} in F_p^*."""
    H = set(mu(p, 1 << k))
    seen, reps = set(), []
    for c in range(1, p):
        if c not in seen:
            reps.append(c)
            seen |= {c * h % p for h in H}
    return reps
