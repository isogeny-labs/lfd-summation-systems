#!/usr/bin/env python3
"""Which permutations of summands and auxiliaries are symmetries of the split presentation?

Candidates: tau = (pi, rho), pi a permutation of the m summands (blocks of the mask), rho a
permutation of the m - 2 auxiliaries.  These are exactly the degree-preserving monomial
permutations that permute the shift variables block-wise.  For each candidate, on a certified
instance:
  inputs   span(tau g_1, .., tau g_{m-1}) = span(g_1, .., g_{m-1})   (rank test, exact mod p)
  closure  tau V_D = V_D for the closed V_D without 1 (every permuted basis row reduces to zero)
A candidate that stabilises the span of the inputs stabilises every V_D; the second test asks
whether V_D has symmetries the inputs do not show.
"""
import itertools, json, pathlib, sys
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'scripts'))
from csd import ic_bounds_split_m as S, ic_bounds_split_gpu as G   # noqa: E402
from csd import ic_bounds_big as B                                  # noqa: E402
import ic_split_gpu_gate as GATE                                    # noqa: E402


def perm_index(sp, pi, rho):
    k, N = sp.k, (1 << sp.k) - 1
    out = np.empty(len(sp), dtype=np.int64)
    for i, (mask, es) in enumerate(sp.mons):
        new = 0
        for blk in range(sp.m):
            new |= ((mask >> (blk * k)) & N) << (pi[blk] * k)
        nes = [0] * len(es)
        for j, e in enumerate(es):
            nes[rho[j]] = e
        out[i] = sp.index[(new, tuple(nes))]
    return out


def rank(M, p):
    return G.host_rref_transform(M, p)[1]


def main():
    out = []
    for (m, k, D) in [(4, 2, 5), (5, 2, 5), (4, 3, 5)]:
        p, gens, _ = GATE.instance(m, k)
        rec, eng = G.closure(gens, m, k, p, D, device='cpu', split=False)
        assert rec['status'] == 'closed' and rec['verified_closed']
        sp = eng.sp
        rows = np.asarray([sp.row(g, p) for g in gens], dtype=np.int64)
        basis = eng.parts[0].host_rows()
        r0 = rank(rows, p)
        sym_inputs, sym_closure, n = [], [], 0
        for pi in itertools.permutations(range(m)):
            for rho in itertools.permutations(range(m - 2)):
                n += 1
                idx = perm_index(sp, pi, rho)
                moved = np.zeros_like(rows); moved[:, idx] = rows
                a = rank(np.concatenate([rows, moved]), p) == r0
                mb = np.zeros_like(basis); mb[:, idx] = basis
                R = eng.parts[0].reduce(B.DV.to_dev(mb.astype(np.float64)))
                c = not bool(np.any(B.DV.to_host(R)))
                if a: sym_inputs.append((pi, rho))
                if c: sym_closure.append((pi, rho))
        line = dict(m=m, k=k, D=D, p=p, candidates=n, dim_V_D=rec['dim_V_D'],
                    stabilise_inputs=sym_inputs, stabilise_V_D=sym_closure)
        out.append(line)
        print(json.dumps(line), flush=True)
    (ROOT / 'experiments/gpu-split-m/raw/symmetry_group.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
