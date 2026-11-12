#!/usr/bin/env python3
"""Item 3: relate the ordinary-ring dimension of V_D to the repository's 'modulo the chains' dimension.

Relation tested:   dim V_D (ordinary)  =  dim( ker NF on P_{<=D} )  +  dim NF(V_D)
  * NF  = reduction by the chain rewriting u_{i,j}^2 -> u_{i,j+1}, u_{i,k-1}^2 -> 1 (own code, monomial -> squarefree
          monomial, degree non-increasing, linear).  Chain equations are a Groebner basis (pure powers in distinct
          variables) so  ker NF = chain ideal C.
  * ker NF cap P_{<=D} has dimension |M_{<=D}| - |squarefree monomials of degree <= D|  ("nonstandard count").
  * C cap P_{<=D} is contained in V_D because the chain equations are inputs and V^chain_D (closure of the chain
    equations alone) already fills it: checked as dim V^chain_D == |M_{<=D}| - rank(NF|P_{<=D}).
  * 'modulo the chains' dimension  q_D := dim NF(V_D) = dim V_D - dim C_{<=D}.
Reproduces: q_3 = 2k+1 (growth.json dim_V3), q_4 (growth.json dim_V4) where affordable.
"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, '/home/user/ml-cryptanalysis')
import numpy as np
import lfd_engine as L
from csd import ic_dreg as I   # allowed helpers only: PrimeCurve, pkm_generator

ROOT = '/home/user/ml-cryptanalysis'


def nf_matrix(sp, k):
    """T[c, c'] = 1 iff NF(monomial c) = monomial c' (columns of the same space, squarefree targets)."""
    r = 2 * k
    T = np.zeros((sp.n, sp.n), dtype=np.float64)
    for c, m in enumerate(sp.cols):
        e = list(m)
        for i in range(2):
            for j in range(k):
                v = i * k + j
                while e[v] >= 2:           # u^2 -> next (or 1 at the end of the chain)
                    e[v] -= 2
                    if j + 1 < k:
                        e[v + 1] += 1
        T[c, sp.idx[tuple(e)]] = 1
    return T


def chain_only_inputs(k, p):
    ins, _ = L.prime_inputs([], k, p)
    return ins[:-1]  # drop the (empty) generator


def profile(gen_items, k, p, Dlist, label):
    r = 2 * k
    inputs, G = L.prime_inputs(gen_items, k, p)
    pts = L.prime_points(G, k, p)
    chain = chain_only_inputs(k, p)
    rows = []
    for D in Dlist:
        t0 = time.time()
        sp = L.Space(r, D)
        T = nf_matrix(sp, k)
        std = int(np.sum(np.all(np.array(sp.cols) <= 1, axis=1)))
        nonstd = sp.n - std
        rank_nf = L.rank_mod(T.T, p)               # image of NF on P_{<=D} = span of squarefree monomials
        kerdim = sp.n - rank_nf
        Vc = L.closure(chain, r, D, p, space=sp)    # chain equations alone
        V = L.closure(inputs, r, D, p, space=sp)
        # NF(V_D): rank of  E @ T  (mod p)
        img = np.mod(V['V'].E @ T, p) if V['V'].rk else np.zeros((0, sp.n))
        q = L.rank_mod(img, p) if img.shape[0] else 0
        # low part V_D cap P_{<=D-1}
        lowmask = V['V'].piv >= sp.low_start
        imglow = img[lowmask]
        qlow = L.rank_mod(imglow, p) if imglow.shape[0] else 0
        # |Z| and I_D
        Mfull = L.eval_matrix(pts, sp, p) if pts else np.zeros((0, sp.n), dtype=np.int64)
        rk_eval = L.rank_mod(Mfull, p) if pts else 0
        dimI = sp.n - rk_eval
        # NF(I_D): dimension of the ideal mod chains = std monomials minus evaluation rank on std monomials
        stdcols = np.flatnonzero(np.all(np.array(sp.cols) <= 1, axis=1))
        rk_eval_std = L.rank_mod(Mfull[:, stdcols], p) if pts else 0
        row = dict(label=label, k=k, p=p, D=D, monomials=sp.n, squarefree=std, nonstandard=nonstd,
                   dim_ker_NF=kerdim, dim_Vchain=Vc['dim'], dim_V_ordinary=V['dim'],
                   dim_NF_V=q, dim_NF_V_low=qlow, dim_I_ordinary=dimI, points=len(pts),
                   relation_holds=(V['dim'] == kerdim + q), chain_fills_ker=(Vc['dim'] == kerdim),
                   one_in_V=V['contains_one'], seconds=round(time.time() - t0, 1))
        rows.append(row)
        print(json.dumps(row), flush=True)
    return rows


if __name__ == '__main__':
    mode = sys.argv[1]
    out = []
    if mode == 'sweep':
        # two sweep systems with k = 3 (one proper, one unit), hashes given as prefixes
        S = {}
        for l in open(ROOT + '/experiments/ic-curve-sweep/systems.jsonl'):
            s = json.loads(l)
            S[s['system_hash']] = s
        R = {}
        for l in open(ROOT + '/experiments/ic-curve-sweep/results.jsonl'):
            r = json.loads(l)
            R[r['system_hash']] = r
        for pref in sys.argv[2:]:
            h = [x for x in S if x.startswith(pref)][0]
            s = S[h]
            gen = [(int(m), int(c)) for m, c in s['generators'][0]]
            rows = profile(gen, s['k'], s['prime'], [1, 2, 3, 4], h[:12])
            print('repo quotient_ranks for', h[:12], R[h]['quotient_ranks'], 'status', R[h]['status'])
            out.append(dict(hash=h, rows=rows, repo_quotient_ranks=R[h]['quotient_ranks']))
    else:
        g = json.load(open(ROOT + '/experiments/ic-depth-step/growth.json'))
        kmax3 = int(sys.argv[2])
        kmax4 = int(sys.argv[3])
        for rr in g['rows']:
            k = rr['k']
            if k > max(kmax3, kmax4):
                continue
            p, a, b, x = rr['prime'], rr['a'], rr['b'], rr['target_x']
            E = I.PrimeCurve(p, a, b)
            G = I.pkm_generator(E, x, k)
            Dl = [3] + ([4] if k <= kmax4 else [])
            rows = profile(sorted(G.items()), k, p, Dl, 'growth p=%d a=%d b=%d x=%d' % (p, a, b, x))
            q3 = [r_ for r_ in rows if r_['D'] == 3][0]
            res = dict(growth=rr, mine_dim_NF_V3=q3['dim_NF_V'], match_dimV3=q3['dim_NF_V'] == rr['dim_V3'],
                       mine_points=q3['points'], match_points=q3['points'] == rr['solutions'])
            if k <= kmax4:
                q4 = [r_ for r_ in rows if r_['D'] == 4][0]
                res.update(mine_dim_NF_V4=q4['dim_NF_V'], match_dimV4=q4['dim_NF_V'] == rr['dim_V4'])
            print('GROWTH', json.dumps(res), flush=True)
            out.append(res)
    json.dump(out, open(os.path.join(HERE, 'item3_%s.json' % mode), 'w'), indent=1)
