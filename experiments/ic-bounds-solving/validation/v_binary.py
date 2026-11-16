#!/usr/bin/env python3
"""Validator check 5: own Boolean static Macaulay d_full for proper binary m = 2 systems (lane C definition):

  B = F_2[x_1..x_r]/(x_i^2 + x_i); equations f_1..f_n (lists of monomial masks, XOR-sum).
  Mac_D (Boolean form) = span{ NF(mu f_l) : mu squarefree, deg mu <= D - deg f_l }, NF(mu f_l) = XOR over monomials
  m of f_l of (mu | m) (Boolean product = union of supports, repeated monomials cancel mod 2).
  Columns: squarefree monomials of degree <= D.  Z: brute force over all 2^r assignments.
  d_full = least D with  rank(Mac_D) = cols_D - |Z|  and  rank(eval of monomials of degree < D on Z) = |Z|.
  (Then Mac_D = I cap M_<=D because Mac_D is inside I cap M_<=D, whose dimension is cols_D - |Z| once the
  evaluation map on Z is onto.)  We continue two more degrees to observe whether the condition stays true.
GF(2) linear algebra: own Python-int bitset elimination (pivot on the highest set bit).
"""
import argparse, json, os, sys, time
from itertools import combinations

SYSTEMS = '/home/user/ml-cryptanalysis/experiments/ic-curve-sweep/systems.jsonl'
TABLE = '/home/user/ml-cryptanalysis/experiments/ic-bounds-solving/laneC_binary/data/binary_table.jsonl'
SCAN = '/home/user/ml-cryptanalysis/experiments/ic-bounds-solving/laneC_binary/data/static_scan_k2_7.jsonl'


def gf2_rank(rows):
    piv = {}
    for row in rows:
        while row:
            h = row.bit_length() - 1
            pr = piv.get(h)
            if pr is None:
                piv[h] = row
                break
            row ^= pr
    return len(piv), piv


def squarefree_monomials(r, dmax):
    out = []
    for d in range(dmax + 1):
        for comb in combinations(range(r), d):
            m = 0
            for i in comb:
                m |= 1 << i
            out.append(m)
    return out


def brute_Z(eqs, r):
    pts = []
    for x in range(1 << r):
        ok = True
        for e in eqs:
            v = 0
            for m in e:
                if (x & m) == m:
                    v ^= 1
            if v:
                ok = False
                break
        if ok:
            pts.append(x)
    return pts


def d_full_scan(eqs, r, Z, Dmax):
    degs = [max(bin(m).count('1') for m in e) for e in eqs]
    prof = []
    d_full = None
    for D in range(max(degs), Dmax + 1):
        t0 = time.time()
        cols = squarefree_monomials(r, D)
        cidx = {m: i for i, m in enumerate(cols)}
        rows = set()
        for e, de in zip(eqs, degs):
            for mu in squarefree_monomials(r, D - de):
                row = 0
                for m in e:
                    row ^= 1 << cidx[mu | m]
                if row:
                    rows.add(row)
        rank, piv = gf2_rank(list(rows))
        one_in = False
        # 1 in row space: reduce e_1 (column of the constant monomial, index 0) by the pivots
        v = 1 << cidx[0]
        while v:
            h = v.bit_length() - 1
            if h in piv:
                v ^= piv[h]
            else:
                break
        one_in = (v == 0)
        low = [m for m in cols if bin(m).count('1') < D]
        ev_rows = []
        for z in Z:
            rr = 0
            for j, m in enumerate(low):
                if (z & m) == m:
                    rr |= 1 << j
            ev_rows.append(rr)
        ev_rank, _ = gf2_rank(ev_rows)
        full = (rank == len(cols) - len(Z)) and (ev_rank == len(Z))
        rec = dict(D=D, columns=len(cols), rows_distinct=len(rows), rank=rank, deficiency=len(cols) - len(Z) - rank,
                   ev_rank_below_D=ev_rank, Z=len(Z), full=full, one_in_rowspace=one_in, sec=round(time.time() - t0, 2))
        prof.append(rec)
        if full and d_full is None:
            d_full = D
        if d_full is not None and D >= d_full + 2:
            break
    return d_full, prof


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hashes', nargs='+', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    want = set(a.hashes)
    sysd = {}
    for l in open(SYSTEMS):
        x = json.loads(l)
        if x.get('kind') == 'binary' and x['system_hash'] in want:
            sysd[x['system_hash']] = x
    table = {}
    for l in open(TABLE):
        x = json.loads(l)
        if x['system_hash'] in want:
            table.setdefault(x['system_hash'], x)
    scan = {}
    for l in open(SCAN):
        x = json.loads(l)
        if x['system_hash'] in want:
            scan[x['system_hash']] = x
    out = []
    for h in a.hashes:
        s = sysd[h]
        eqs = [[int(m) for m in e] for e in s['equations']]
        r = s['variables']
        t0 = time.time()
        Z = brute_Z(eqs, r)
        tb = table.get(h, {})
        sc = scan.get(h, {})
        d_full, prof = d_full_scan(eqs, r, Z, Dmax=min(r, (sc.get('d_full') or 6) + 3))
        rec = dict(system_hash=h, r=r, k=r // 2, n=len(eqs), Z=len(Z), laneC=dict(Z=tb.get('Z'), d_full=tb.get('d_full'), lfd=tb.get('lfd'),
                   dstar=tb.get('dstar'), per_degree=[(d['D'], d['columns'], d['rank'], d.get('deficiency'), d.get('ev')) for d in sc.get('per_degree', [])]),
                   mine=dict(d_full=d_full, per_degree=[(d['D'], d['columns'], d['rank'], d['deficiency'], d['ev_rank_below_D'], d['full'], d['one_in_rowspace']) for d in prof]),
                   agree_Z=(len(Z) == tb.get('Z')), agree_d_full=(d_full == tb.get('d_full')),
                   ranks_agree=all(any(md[0] == ld[0] and md[1] == ld[1] and md[2] == ld[2] for ld in [(d['D'], d['columns'], d['rank']) for d in sc.get('per_degree', [])])
                                   for md in [(d['D'], d['columns'], d['rank']) for d in prof if d['D'] <= (sc.get('d_full') or 0)]),
                   stays_full_after=all(d['full'] for d in prof if d_full is not None and d['D'] >= d_full),
                   ge_lfd=(d_full is not None and tb.get('lfd') is not None and d_full >= tb['lfd']),
                   seconds=round(time.time() - t0, 1))
        out.append(rec)
        print(f"{h[:12]} k={rec['k']} n={rec['n']} Z={len(Z)} (laneC Z={tb.get('Z')}) my d_full={d_full} laneC d_full={tb.get('d_full')} lfd={tb.get('lfd')} "
              f"D*={tb.get('dstar')} agree={rec['agree_d_full']} ranks_agree={rec['ranks_agree']} stays_full={rec['stays_full_after']} {rec['seconds']}s", flush=True)
        print('   mine   :', rec['mine']['per_degree'])
        print('   laneC  :', rec['laneC']['per_degree'])
    json.dump(out, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
