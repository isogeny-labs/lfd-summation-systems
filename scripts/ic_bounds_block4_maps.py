"""Tests A, C, N of docs/proofs/N9-block-four.md: candidate depth-raising maps between closures.

  doubling  OUT D k r_max          Test A.  Is g^(r) = g(x1^(2^r), x2^(2^r)) in V_D(k)?  And the dimension of
                                   L_r = {Q in span(1, u_{i,j}) : Q g^(r) in V_D(k)} (weight-<=1 multipliers),
                                   plus whether x1 - beta is in L_r for a beta of order not dividing 2^k.
  curve     OUT r_max npts         Test A'.  Values of S_3(x1^(2^r), x2^(2^r), t) at points (x1, x2) of the
                                   curve S_3(x1, x2, t) = 0 over F_p (count of nonzero values).
  coset     OUT D k1 k2 r          Test C.  Coset instances g_s(y, x2) = g(s y, x2), s over representatives of
                                   mu_{2^(k1+r)} / mu_{2^k1}: is 1 in V_D[g_s](k1, k2)?
  norm      OUT D k1 k2 r          Test N.  h_r(X, x2) = prod_{zeta in mu_{2^r}} g(zeta x1, x2) (a polynomial in
                                   X = x1^(2^r)): weight of h_r, is 1 in V_D[h_r](k1, k2)?
  normentry OUT Dmax k1 k2 r       Test N.  least E <= Dmax with h_r(x1^(2^r), x2) in V_E[g](k1, k2).

Unit instance of experiments/ic-bounds-k14 (p = 65537); every closure is followed by a deterministic full
sweep (field sweep_failures, 0 = closed).
"""
import json
import sys
import time

import numpy as np

import ic_bounds_block4_lib as B4

INST = B4.UNIT
P = INST['p']
G = B4.generator(INST)


def emit(out, rec):
    rec = dict(rec, inst=INST)
    print(json.dumps(rec), flush=True)
    with open(out, 'a') as fh:
        fh.write(json.dumps(rec) + '\n')


def frob(g, r):
    return {(e1 << r, e2 << r): c for (e1, e2), c in g.items()}


def pmul(f, h, p):
    out = {}
    for (a1, a2), c in f.items():
        for (b1, b2), d in h.items():
            key = (a1 + b1, a2 + b2)
            out[key] = (out.get(key, 0) + c * d) % p
    return {e: v for e, v in out.items() if v}


def norm_x1(g, r, p):
    """prod over zeta in mu_{2^r} of g(zeta x1, x2), as a dict in (X = x1^(2^r), x2)."""
    z = B4.root_of_unity(p, r) if r > 0 else 1
    h = {(0, 0): 1}
    for e in range(1 << r):
        zz = pow(z, e, p)
        h = pmul(h, {(e1, e2): c * pow(zz, e1, p) % p for (e1, e2), c in g.items()}, p)
    assert all(e1 % (1 << r) == 0 for e1, _ in h), 'norm not a polynomial in x1^(2^r)'
    return {(e1 >> r, e2): c for (e1, e2), c in h.items()}


def wt_of(d, k1, k2):
    return max(B4.popcount(e1 % (1 << k1)) + B4.popcount(e2 % (1 << k2)) for e1, e2 in d)


def reduce_dict(d, k1, k2, p):
    out = {}
    for (e1, e2), c in d.items():
        key = (e1 % (1 << k1), e2 % (1 << k2))
        out[key] = (out.get(key, 0) + c) % p
    return {e: v for e, v in out.items() if v}


def closure(k1, k2, D, gen):
    C = B4.FastClosure(k1, k2, D, P, [gen])
    C.run()
    return C, C.sweep_full()


def cmd_doubling(out, D, k, rmax):
    t0 = time.time()
    C, bad = closure(k, k, D, G)
    one = C.contains_one()
    beta = 3  # 3 is a primitive root mod 65537: order 2^16, not dividing 2^k for k < 16
    for r in range(1, rmax + 1):
        gr = reduce_dict(frob(G, r), k, k, P)
        inn = C.contains(gr)
        # multipliers Q in span(1, u_{i,j}); Q g^(r) has weight <= wt(g^(r)) + 1
        mults = [{(0, 0): 1}] + [{(1 << j, 0): 1} for j in range(k)] + [{(0, 1 << j): 1} for j in range(k)]
        rows, ok = [], True
        for m in mults:
            try:
                rows.append(C.L.vec(reduce_dict(pmul(m, gr, P), k, k, P), P))
            except KeyError:
                ok = False
                break
        dimL = lin = None
        if ok:
            R = C.P.residual(np.array(rows))          # residuals of Q_m g^(r) modulo V_D(k)
            # L_r = left kernel of R: dim = #mults - rank(R)
            _, piv = B4.rref(R.astype(np.int64).T, P) if R.any() else (None, [])
            rank = len(piv)
            dimL = len(mults) - rank
            lin = C.contains(reduce_dict(pmul({(1, 0): 1, (0, 0): (-beta) % P}, gr, P), k, k, P))
        emit(out, dict(test='doubling', D=D, k=k, r=r, V_dim=C.dim(), V_one=one, sweep_failures=bad,
                       g_r_weight=wt_of(gr, k, k), g_r_in_V=inn, dim_L_weight1_multipliers=dimL,
                       multipliers=len(mults), x1_minus_beta_times_g_r_in_V=lin, beta=beta,
                       seconds=round(time.time() - t0, 2)))


def cmd_curve(out, rmax, npts):
    pts = []
    rng = np.random.default_rng(7)
    co = G
    while len(pts) < npts:
        x1 = int(rng.integers(1, P))
        c = [0, 0, 0]
        for (e1, e2), v in co.items():
            c[e2] = (c[e2] + v * pow(x1, e1, P)) % P
        if c[2] == 0:
            continue
        disc = (c[1] * c[1] - 4 * c[2] * c[0]) % P
        if pow(disc, (P - 1) // 2, P) != 1:
            continue
        s = next(y for y in range(P) if y * y % P == disc) if disc else 0
        x2 = (-c[1] + s) * pow(2 * c[2], P - 2, P) % P
        if x2:
            pts.append((x1, x2))
    def ev(d, x1, x2):
        return sum(v * pow(x1, e1, P) * pow(x2, e2, P) for (e1, e2), v in d.items()) % P
    assert all(ev(G, a, b) == 0 for a, b in pts)
    for r in range(1, rmax + 1):
        gr = frob(G, r)
        nz = sum(1 for a, b in pts if ev(gr, a, b) != 0)
        emit(out, dict(test='curve', r=r, points=len(pts), nonzero_values_of_g_r_on_curve=nz))
    # weight-<=1 multipliers Q = c + sum c_{ij} x_i^(2^j), j < J: does a nonzero Q vanish at all the points?
    # J = 16: all weight-1 functions on F_p^* (x^(2^16) = 1 there); J = 15: exponents < 2^15 only.
    for J in (16, 15):
        funcs = [lambda a, b: 1] + [lambda a, b, j=j: pow(a, 1 << j, P) for j in range(J)] + \
                [lambda a, b, j=j: pow(b, 1 << j, P) for j in range(J)]
        M = np.array([[f(a, b) for f in funcs] for a, b in pts], dtype=np.int64)
        _, piv = B4.rref(M, P)
        emit(out, dict(test='curve-weight1', J=J, functions=len(funcs), points=len(pts), evaluation_rank=len(piv),
                       note='J = 16: the one relation is x1^(2^15) + x2^(2^15) = 0 on C(F_p) (Legendre symbols; '
                            'the unit property of this instance), a relation on F_p-points only'))

def cmd_coset(out, D, k1, k2, r):
    t0 = time.time()
    z = B4.root_of_unity(P, k1 + r)        # generator of mu_{2^(k1+r)}; s = z^c, c < 2^r, represent the cosets
    res = []
    for c in range(1 << r):
        s = pow(z, c, P)
        gs = {(e1, e2): v * pow(s, e1, P) % P for (e1, e2), v in G.items()}
        C, bad = closure(k1, k2, D, gs)
        res.append(dict(c=c, one=C.contains_one(), dim=C.dim(), profile=C.profile(), sweep_failures=bad,
                        solutions=coset_solutions(gs, k1, k2)))
    emit(out, dict(test='coset', D=D, k1=k1, k2=k2, r=r, cosets=res, all_one=all(x['one'] for x in res),
                   seconds=round(time.time() - t0, 2)))


def coset_solutions(gs, k1, k2):
    h1 = B4.root_of_unity(P, k1); h2 = B4.root_of_unity(P, k2)
    H2 = np.array([pow(h2, e, P) for e in range(1 << k2)], dtype=np.int64)
    pw2 = [np.ones_like(H2), H2, H2 * H2 % P]
    cnt = 0
    for e in range(1 << k1):
        x1 = pow(h1, e, P)
        v = np.zeros_like(H2)
        for (e1, e2), c in gs.items():
            v = (v + (c * pow(x1, e1, P) % P) * pw2[e2]) % P
        cnt += int(np.count_nonzero(v == 0))
    return cnt


def cmd_norm(out, D, k1, k2, r):
    t0 = time.time()
    h = norm_x1(G, r, P)
    hr = reduce_dict(h, k1, k2, P)
    w = wt_of(hr, k1, k2)
    if w > D:
        emit(out, dict(test='norm', D=D, k1=k1, k2=k2, r=r, h_weight=w, h_terms=len(hr), one=False,
                       note='weight of h_r exceeds D: h_r not an input at this budget'))
        return
    C, bad = closure(k1, k2, D, hr)
    emit(out, dict(test='norm', D=D, k1=k1, k2=k2, r=r, h_weight=w, h_terms=len(hr), one=C.contains_one(),
                   dim=C.dim(), profile=C.profile(), sweep_failures=bad, seconds=round(time.time() - t0, 2)))


def cmd_normentry(out, Dmax, k1, k2, r):
    t0 = time.time()
    h = norm_x1(G, r, P)
    lifted = reduce_dict({(e1 << r, e2): c for (e1, e2), c in h.items()}, k1, k2, P)
    w = wt_of(lifted, k1, k2)
    entry = None
    for E in range(max(2, w), Dmax + 1):
        C, bad = closure(k1, k2, E, G)
        if C.contains(lifted):
            entry = dict(E=E, V_one=C.contains_one(), sweep_failures=bad)
            break
    emit(out, dict(test='normentry', Dmax=Dmax, k1=k1, k2=k2, r=r, lifted_weight=w, entry=entry,
                   seconds=round(time.time() - t0, 2)))


if __name__ == '__main__':
    cmd = sys.argv[1]
    a = [int(x) for x in sys.argv[3:]]
    {'doubling': cmd_doubling, 'curve': cmd_curve, 'coset': cmd_coset, 'norm': cmd_norm,
     'normentry': cmd_normentry}[cmd](sys.argv[2], *a)
