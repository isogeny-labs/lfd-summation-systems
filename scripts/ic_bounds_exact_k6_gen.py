#!/usr/bin/env python3
"""Kernel certificates for an exact last fall degree at depth 6 (two-summand prime chain).

For one instance (p, a, b, t) with 64 | p - 1 and a unit ideal, this script builds the two
certificates checked by the Lean modules `ICFallDegrees.ExactK6Cert` / `ExactK6LowerSound` /
`ExactK6UpperSound`:

  lower : a reduced echelon basis of the degree-4 quotient closure V_4(6) (the group algebra
          F_p[(Z/64)^2] modulo the chains), in the `N5Cert.Data` format with complement codes
          m = 4095 - (e0 + 64 e1); its rows of pivot weight <= 3 as sparse lists.  The checker
          tests that all 12 translates of every such row lie in the row space, that the
          generator does and that the constant is a free column:  1 is not in V_4.
  upper : a derivation of 1 in the degree-5 quotient closure.  Row 0 is the generator; row n is
          a combination  sum c * tau_s(row i)  with i < n and wt(s) + deg(row i) <= 5; the last
          row is 1.  Rows are packed in one natural with 128 x 128 slots (slot e0 + 128 e1), so
          tau_s is one left shift and the wrap is folded by the final scan.

Everything the Lean checker computes is recomputed here entry by entry (`ref_check_lower`,
`ref_check_upper`: the Python reference checker), and the Lean data, kernel and instance modules
are emitted.

    python3 scripts/ic_bounds_exact_k6_gen.py --params 1153,139,1094,185 --suffix P1153 \
        --out-dir experiments/ic-bounds-exact-k6/p1153-a139-b1094-t185

The linear algebra is self-contained (numpy, int64, exact mod p); no repository engine is
imported.  Resident memory stays below 1 GB.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

N = 64
NN = N * N
WT = np.array([bin(c % N).count('1') + bin(c // N).count('1') for c in range(NN)])
MUS = [[mu for mu in range(NN) if WT[mu] == w] for w in range(13)]
ORDER = sorted(range(NN), key=lambda c: (-WT[c], c))   # heaviest first


def popcount(n: int) -> int:
    return bin(n).count('1')


# ----------------------------------------------------------------------------- group algebra

def generator(p: int, a: int, b: int, t: int) -> np.ndarray:
    """Normal form of S_3(x1, x2, t) modulo the chains: exponent pairs (e0, e1) in {0,1,2}^2,
    index e0 + 64 e1.  Same polynomial as `PrimeM2General.normalized`."""
    g = np.zeros(NN, dtype=np.int64)
    co = {(2, 2): 1, (2, 1): -2 * t, (1, 2): -2 * t, (2, 0): t * t, (0, 2): t * t,
          (1, 1): -2 * (t * t + a), (1, 0): -2 * a * t - 4 * b, (0, 1): -2 * a * t - 4 * b,
          (0, 0): a * a - 4 * b * t}
    for (e0, e1), v in co.items():
        g[e0 + N * e1] = v % p
    return g


def generator_entries(p: int, a: int, b: int, t: int) -> list[tuple[int, int, int]]:
    g = generator(p, a, b, t)
    out = [(e0, e1, int(g[e0 + N * e1])) for e1 in range(3) for e0 in range(3)]
    assert all(v != 0 for _, _, v in out), 'a generator coefficient vanishes: instance not supported'
    return out


def translate(v: np.ndarray, mu: int) -> np.ndarray:
    m0, m1 = mu % N, mu // N
    return np.roll(np.roll(v.reshape(N, N), m1, axis=0), m0, axis=1).reshape(-1)


def count_zeros(p: int, a: int, b: int, t: int) -> int:
    H = sorted({pow(x, (p - 1) // N, p) for x in range(1, p)})
    assert len(H) == N
    return sum(1 for u in H for w in H
               if ((u - w) ** 2 * t * t - 2 * ((u + w) * (u * w + a) + 2 * b) * t + (u * w - a) ** 2
                   - 4 * b * (u + w)) % p == 0)


class Echelon:
    """Reduced row echelon basis over F_p on a fixed column list, built by insertion; optionally
    every basis row is tracked as a combination of the inserted independent vectors."""

    def __init__(self, p: int, cols: np.ndarray, track: bool):
        self.p, self.cols, self.m, self.track = p, cols, len(cols), track
        self.R = np.zeros((0, self.m), dtype=np.int64)
        self.piv: list[int] = []
        self.T = np.zeros((0, 0), dtype=np.int64)
        self.lab: list = []

    def reduce(self, v):
        p = self.p
        if not self.piv:
            return v % p, None
        co = v[self.piv] % p
        nz = np.nonzero(co)[0]
        if len(nz) == 0:
            return v % p, co
        return (v - co[nz] @ self.R[nz]) % p, co

    def contains(self, v) -> bool:
        return not self.reduce(v)[0].any()

    def insert(self, v, label=None) -> bool:
        p = self.p
        r, co = self.reduce(v)
        nzr = np.nonzero(r)[0]
        if len(nzr) == 0:
            return False
        c = int(nzr[0])
        iv = pow(int(r[c]), p - 2, p)
        n = len(self.piv)
        if self.track:
            tv = np.zeros(n + 1, dtype=np.int64)
            tv[n] = 1
            Tn = np.zeros((n + 1, n + 1), dtype=np.int64)
            Tn[:n, :n] = self.T
            if n:
                nz = np.nonzero(co)[0]
                tv[:n] = (-(co[nz] @ self.T[nz])) % p
            tv = tv * iv % p
            self.lab.append(label)
        r = r * iv % p
        if n:
            col = self.R[:, c].copy()
            nz = np.nonzero(col)[0]
            if len(nz):
                self.R[nz] = (self.R[nz] - np.outer(col[nz], r)) % p
                if self.track:
                    Tn[nz] = (Tn[nz] - np.outer(col[nz], tv)) % p
        self.R = np.vstack([self.R, r])
        self.piv.append(c)
        if self.track:
            Tn[n] = tv
            self.T = Tn
        return True


# ----------------------------------------------------------------------------- lower: V_4(6)

def closure4(p: int, g: np.ndarray, log=print):
    """Reduced echelon basis of the degree-4 quotient closure (columns of weight <= 4, heaviest
    first).  Closedness is not trusted from here: the reference checker re-tests it."""
    cols = np.array([c for c in ORDER if WT[c] <= 4])
    E = Echelon(p, cols, track=False)
    for w in range(3):
        for mu in MUS[w]:
            E.insert(translate(g, mu)[cols])
    gens = [1 << j for j in range(6)] + [N << j for j in range(6)]
    done: set[bytes] = set()
    while True:
        grew = False
        wcol = WT[cols]
        for i in [i for i, c in enumerate(E.piv) if wcol[c] <= 3]:
            row = E.R[i].copy()
            key = row.tobytes()
            if key in done:
                continue
            done.add(key)
            full = np.zeros(NN, dtype=np.int64)
            full[cols] = row
            for mu in gens:
                grew |= E.insert(translate(full, mu)[cols])
        if not grew:
            break
    prof = [int((WT[cols][E.piv] <= e).sum()) for e in range(5)]
    log(f'[lower] dim V_4 = {len(E.piv)}, profile dim(V_4 ∩ M_<=j) = {prof}')
    return E.R, E.piv, cols, prof


def slot_width(bound: int) -> int:
    """Smallest multiple of 4 with bound < 2^w."""
    w = bound.bit_length()
    return (w + 3) // 4 * 4


def build_lower(p: int, a: int, b: int, t: int, log=print) -> dict:
    g = generator(p, a, b, t)
    R, piv, cols, prof = closure4(p, g, log)
    nrows, ncols = R.shape
    code = lambda c: 4095 - int(c)
    w = slot_width(512 * p * p)
    piv_set = set(piv)
    free = [j for j in range(ncols) if j not in piv_set]
    piv_codes = [code(cols[j]) for j in piv]
    # free column 0 is a dummy with code 0 (the monomial of weight 12, never a column): N5Cert's
    # table cannot address free column 0 in an entry, and `globalOK` wants code 0 free.
    free_codes = [0] + [code(cols[j]) for j in free]
    assert nrows < 1024 and 0 not in piv_codes and len(set(free_codes)) == len(free_codes)
    tbl = 0
    for i, m in enumerate(piv_codes):
        tbl |= (i + 1) << (16 * m)
    for f, m in enumerate(free_codes):
        tbl |= (1024 + f) << (16 * m)
    piv_enc = sum(m << (16 * i) for i, m in enumerate(piv_codes))
    free_enc = sum(m << (16 * f) for f, m in enumerate(free_codes))
    packed, sparse = [], []
    for i in range(nrows):
        Rw, entries = 0, []
        for f, j in enumerate(free, start=1):
            v = int(R[i, j])
            if v:
                Rw |= v << (w * f)
                entries.append((code(cols[j]), v))
        packed.append(Rw)
        sparse.append(entries)
    low_idx = [i for i, m in enumerate(piv_codes) if popcount(m) >= 9]
    low = [(i, sparse[i]) for i in low_idx]
    gen = [(code(e0 + N * e1), v) for e0, e1, v in generator_entries(p, a, b, t)]
    return dict(p=p, w=w, nrows=nrows, nfree=len(free_codes), piv_codes=piv_codes,
                free_codes=free_codes, tbl=tbl, piv_enc=piv_enc, free_enc=free_enc, packed=packed,
                low=low, gen=gen, profile=prof,
                free_entries=sum(len(s) for s in sparse))


def tree_path(i: int):
    """Path to index i under N5Cert.Tree.get: 0 = root; odd -> left with i//2; even -> right with
    i//2 - 1."""
    path = []
    while i != 0:
        if i % 2 == 1:
            path.append('L')
            i = i // 2
        else:
            path.append('R')
            i = i // 2 - 1
    return path


def build_tree(values: list[int]):
    root: dict = {}
    for i, v in enumerate(values):
        node = root
        for step in tree_path(i):
            node = node.setdefault(step, {})
        node['v'] = v
        node['i'] = i
    return root


def tree_get(node, i: int) -> int:
    for step in tree_path(i):
        node = node.get(step)
        if node is None:
            return 0
    return node.get('v', 0)


def shift_code(i: int, j: int, m: int) -> int:
    """ExactK6Cert.shiftCode."""
    c = max(4095 - m, 0)
    if i == 0:
        return max(4095 - ((c % 64 + 2 ** j) % 64 + 64 * (c // 64)), 0)
    return max(4095 - (c % 64 + 64 * ((c // 64 + 2 ** j) % 64)), 0)


def ref_check_lower(cert: dict, log=print) -> dict:
    """Python mirror of the Lean checks of the lower half (N5Cert.globalOK/rowOK/freeOK/lowIdxOK/
    member and ExactK6Cert.boxOK/constOK/lowIdxOK/lowRowOK)."""
    p, w, nrows, nfree = cert['p'], cert['w'], cert['nrows'], cert['nfree']
    B = 1 << w
    tree = build_tree(cert['packed'])
    tbl_get = lambda m: (cert['tbl'] >> (16 * m)) & 65535
    piv = lambda i: (cert['piv_enc'] >> (16 * i)) & 65535
    free_code = lambda f: (cert['free_enc'] >> (16 * f)) & 65535
    row = lambda i: tree_get(tree, i)
    stats = dict(member_tests=0, ops=0)

    def entry_ok(m, x, free_only=False):
        if not (0 < x < p and m < 16384):
            return False
        V = tbl_get(m)
        if V == 0:
            return False
        v = V - 1
        if not free_only and v < 1024:
            return v < nrows and piv(v) == m
        return 1024 <= v and v + 1 - 1024 < nfree and free_code(v + 1 - 1024) == m

    def accum(u):
        acc = 0
        for m, x in u:
            V = tbl_get(m)
            if V == 0:
                continue
            v = V - 1
            acc += x * row(v) if v < 1024 else (p - x) << (w * (v + 1 - 1024))
        return acc

    def slots_zero(Rw):
        for _ in range(nfree):
            if Rw % B % p != 0:
                return False
            Rw //= B
        return Rw == 0

    def member(u):
        stats['member_tests'] += 1
        stats['ops'] += 50 * len(u) + 4 * nfree
        return len(u) <= 512 and all(entry_ok(m, x) for m, x in u) and slots_zero(accum(u))

    def pack_free(L):
        return sum(x << (w * (tbl_get(m) - 1024)) for m, x in L)

    assert 512 * p * p < B and tbl_get(0) >= 1024 and p >= 3, 'globalOK'
    for i in range(nrows):
        assert tbl_get(piv(i)) == i + 1 and piv(i) < 16384, f'rowOK {i}'
        Rw = row(i)
        for _ in range(nfree):
            assert Rw % B < p, f'rowOK slots {i}'
            Rw //= B
        assert Rw == 0
    for f in range(nfree):
        assert tbl_get(free_code(f)) >= 1024 and free_code(f) < 16384, f'freeOK {f}'
    assert not [i for i in range(nrows) if popcount(piv(i)) <= 3], 'N5 lowIdxOK with low = []'
    assert member(cert['gen']), 'generator not in the row space'
    assert all(piv(i) < 4096 for i in range(nrows)) and all(free_code(f) < 4096 for f in range(nfree)), 'boxOK'
    assert tbl_get(4095) >= 1024, 'constOK'
    assert [i for i in range(nrows) if popcount(piv(i)) >= 9] == [i for i, _ in cert['low']], 'lowIdxOK'
    max_len = 0
    for idx, entries in cert['low']:
        assert idx < nrows and len(entries) <= 256
        assert all(entry_ok(m, x, free_only=True) for m, x in entries)
        assert pack_free(entries) == row(idx)
        full = [(piv(idx), 1)] + entries
        for j in range(6):
            for i in (0, 1):
                u = [(shift_code(i, j, m), x) for m, x in full]
                max_len = max(max_len, len(u))
                assert member(u), f'low row {idx}: translate ({i},{j}) not in the row space'
    log(f'[lower] reference check passed: {nrows} rows, {nfree} free columns (1 dummy), '
        f'{len(cert["low"])} low rows, {stats["member_tests"]} membership tests, '
        f'longest test vector {max_len}, estimated kernel operations {stats["ops"]:.3g}')
    return stats


# ----------------------------------------------------------------------------- upper: derivation of 1

def derive_one(p: int, g: np.ndarray, log=print):
    """Derivation of 1 in the degree-5 quotient closure.

    Generators h_0 = g, h_1, ...: each new h_n is an element of lower degree than expected in the
    span of the admissible translates { tau_mu h_i : wt(mu) <= 5 - deg h_i } of the earlier ones.
    Per round only the new elements of the least degree are kept (their translates reach
    furthest).  Returns (vector, degree, terms) per generator, terms = (coefficient, i, mu)."""
    D = 5
    cols = np.array([c for c in ORDER if WT[c] <= D])
    wcol = WT[cols]
    G = [(g, 2, None)]
    Wsp = Echelon(p, cols, track=True)
    Kn = {e: Echelon(p, cols, track=False) for e in (2, 3, 4)}

    def add_translates(i):
        h, d, _ = G[i]
        for wt in range(0, D - d + 1):
            for mu in MUS[wt]:
                v = translate(h, mu)[cols]
                Wsp.insert(v, (i, mu))
                for e in (2, 3, 4):
                    if d + wt <= e:
                        Kn[e].insert(v)

    def terms_of(i):
        return [(int(Wsp.T[i, j]),) + Wsp.lab[j] for j in np.nonzero(Wsp.T[i])[0]]

    def full(v):
        out = np.zeros(NN, dtype=np.int64)
        out[cols] = v
        return out

    add_translates(0)
    rnd = 0
    while True:
        rnd += 1
        pw = wcol[Wsp.piv]
        prof = [int((pw <= e).sum()) for e in range(6)]
        log(f'[upper] round {rnd}: {len(G)} generators, rank {len(Wsp.piv)}, profile {prof}')
        if prof[0] > 0:
            i = int(np.nonzero(pw == 0)[0][0])
            G.append((full(Wsp.R[i]), 0, terms_of(i)))
            break
        new = []
        for e in (1, 2, 3, 4):
            for i in np.nonzero(pw == e)[0]:
                v = Wsp.R[i]
                if e >= 2 and Kn[e].contains(v):
                    continue
                new.append((full(v), e, terms_of(i)))
                if e >= 2:
                    Kn[e].insert(v)
            if new:
                break
        if not new:
            raise SystemExit('derivation stuck: 1 is not in the degree-5 closure of this instance')
        for x in new:
            G.append(x)
            add_translates(len(G) - 1)
    # keep only the generators the final row depends on
    need = {len(G) - 1}
    stack = [len(G) - 1]
    while stack:
        n = stack.pop()
        for (_, i, _) in (G[n][2] or []):
            if i not in need:
                need.add(i)
                stack.append(i)
    order = sorted(need)
    assert order[0] == 0
    ren = {n: k for k, n in enumerate(order)}
    out = []
    for n in order:
        vec, d, terms = G[n]
        out.append((vec, d, None if terms is None else sorted((c, ren[i], mu) for c, i, mu in terms)))
    log(f'[upper] derivation: {len(out)} rows (of {len(G)} generators found), '
        f'{sum(len(tm or []) for _, _, tm in out)} terms, '
        f'rows by degree {[sum(1 for _, d, _ in out if d == e) for e in range(5)]}')
    return out


def build_upper(p: int, a: int, b: int, t: int, log=print) -> dict:
    g = generator(p, a, b, t)
    rows = derive_one(p, g, log)
    worst = max(sum(c for c, _, _ in tm) for _, _, tm in rows[1:])
    w = slot_width(worst * p)
    packed, degs, terms = [], [], []
    for vec, d, tm in rows:
        H = 0
        for c in np.nonzero(vec)[0]:
            c = int(c)
            H |= int(vec[c]) << (w * (c % N + 128 * (c // N)))
        packed.append(H)
        degs.append(d)
        terms.append(None if tm is None else [(c, mu % N + 128 * (mu // N), i) for c, i, mu in tm])
    gen = [(e0 + 128 * e1, v) for e0, e1, v in generator_entries(p, a, b, t)]
    assert len(degs) < 4096
    return dict(p=p, w=w, packed=packed, degs=degs, terms=terms, gen=gen,
                deg_enc=sum(d << (4 * n) for n, d in enumerate(degs)))


def ref_check_upper(cert: dict, log=print) -> dict:
    """Python mirror of ExactK6Cert.uGlobalOK / genOK / rowCheck and of the final test."""
    p, w = cert['p'], cert['w']
    B = 1 << w
    tree = build_tree(cert['packed'])
    row = lambda n: tree_get(tree, n)
    deg = lambda n: (cert['deg_enc'] >> (4 * n)) & 15
    stats = dict(terms=0, big_ops=0, small_ops=0, max_terms=0)

    def row_scan(n, acc):
        """scanBlocks/scanSlots, slot by slot."""
        H, d = row(n), deg(n)
        lo, hi = acc % B ** 8192, acc // B ** 8192
        M128, M64 = B ** 128, B ** 64
        for e1 in range(64):
            xs = [lo % M128 % M64, lo % M128 // M64, hi % M128 % M64, hi % M128 // M64]
            h = H % M128
            for e0 in range(64):
                hs = h % B
                if (sum(x % B for x in xs) + (p - 1) * hs) % p != 0:
                    return False
                if not hs < p:
                    return False
                if not (hs == 0 or popcount(e0) + popcount(e1) <= d):
                    return False
                xs = [x // B for x in xs]
                h //= B
            if h != 0:
                return False
            lo, hi, H = lo // M128, hi // M128, H // M128
        stats['small_ops'] += 64 * 64 * 25
        return hi == 0 and H == 0

    assert p >= 2 and w >= 1, 'uGlobalOK'
    gen = cert['gen']
    assert row(0) == sum(x << (w * c) for c, x in gen) and sum(x for _, x in gen) < B, 'genOK pack'
    assert all(c < 16384 for c, _ in gen) and row_scan(0, row(0)), 'genOK scan'
    nrows = len(cert['packed'])
    for n in range(1, nrows):
        terms = cert['terms'][n]
        for c, s, i in terms:
            assert i < n and s % 128 < 64 and s // 128 < 64, f'termOK row {n}'
            assert popcount(s % 128) + popcount(s // 128) + deg(i) <= 5, f'degree budget row {n}'
        assert sum(c for c, _, _ in terms) * p < B, f'slot bound row {n}'
        acc = sum(c * (row(i) << (w * s)) for c, s, i in terms)
        assert row_scan(n, acc), f'row {n}: the combination does not fold to the row'
        stats['terms'] += len(terms)
        stats['big_ops'] += 3 * len(terms)
        stats['max_terms'] = max(stats['max_terms'], len(terms))
    assert row(nrows - 1) == 1, 'the last row is not 1'
    log(f'[upper] reference check passed: {nrows} rows, {stats["terms"]} terms '
        f'(largest row {stats["max_terms"]}), slot width {w}, '
        f'{stats["big_ops"]} big-number operations, about {stats["small_ops"]:.3g} scan operations')
    return stats


# ----------------------------------------------------------------------------- Lean emission

def hx(n: int) -> str:
    return hex(n)


def emit_tree(L: list[str], values: list[int], prefix: str, split_depth: int, named: bool = False) -> str:
    """Emit a Tree holding values[i] at Tree.get index i; subtrees below `split_depth` become
    their own definitions (keeps every definition small).  With `named`, every value is its own
    definition `<prefix>v<i>`.  Returns the root expression."""
    root = build_tree(values)
    if named:
        for i, v in enumerate(values):
            L.append(f'def {prefix}v{i} : Nat := {hx(v)}')
        L.append('')

    def val(node) -> str:
        if named and 'i' in node:
            return f"{prefix}v{node['i']}"
        return hx(node.get('v', 0))

    def whole(node) -> str:
        if node is None:
            return '.leaf'
        return f"(.node {whole(node.get('L'))} {val(node)} {whole(node.get('R'))})"

    def emit(node, name, depth):
        if node is None:
            return '.leaf'
        if depth == 0:
            L.append(f'def {name} : Tree := {whole(node)}')
            L.append('')
            return name
        l = emit(node.get('L'), name + 'l', depth - 1)
        r = emit(node.get('R'), name + 'r', depth - 1)
        return f"(.node {l} {val(node)} {r})"

    return emit(root, prefix, split_depth)


def lean_data(lower: dict, upper: dict, inst: dict, suffix: str, low_chunk: int, term_chunk: int) -> str:
    L = []
    L.append('import ICFallDegrees.ExactK6Cert')
    L.append('')
    L.append('/-! Generated by `scripts/ic_bounds_exact_k6_gen.py`; do not edit by hand.')
    L.append(f'Instance `{inst["tag"]}`: p = {inst["p"]}, a = {inst["a"]}, b = {inst["b"]}, t = {inst["t"]}, depth 6.')
    L.append(f'Lower half: reduced echelon basis of V_4(6), {lower["nrows"]} rows, {lower["nfree"]} free columns '
             f'(one dummy), {len(lower["low"])} rows of pivot weight ≤ 3, slot width {lower["w"]}.')
    L.append(f'Upper half: derivation of 1 in the degree-5 quotient closure, {len(upper["packed"])} rows, '
             f'{sum(len(t or []) for t in upper["terms"])} terms, slot width {upper["w"]}.')
    L.append('Encodings: `ICFallDegrees.ExactK6Cert`. -/')
    L.append('')
    L.append(f'namespace ICFallDegrees.ExactK6Data{suffix}')
    L.append('open ICFallDegrees.N5Cert ICFallDegrees.ExactK6Cert')
    L.append('')
    L.append(f'def p : Nat := {inst["p"]}')
    L.append(f'def a : Nat := {inst["a"]}')
    L.append(f'def b : Nat := {inst["b"]}')
    L.append(f'def t : Nat := {inst["t"]}')
    L.append('')
    L.append('/-! ## Lower half -/')
    L.append('')
    L.append(f'def pivEnc : Nat := {hx(lower["piv_enc"])}')
    L.append(f'def freeEnc : Nat := {hx(lower["free_enc"])}')
    L.append(f'def tbl : Nat := {hx(lower["tbl"])}')
    L.append('def lgen : List Entry := [' + ', '.join(f'({m}, {x})' for m, x in lower['gen']) + ']')
    L.append('')
    top = emit_tree(L, lower['packed'], 'ltree', 3)
    L.append(f'def lrows : Tree := {top}')
    L.append('')
    L.append('/-- The echelon data; `low = []` because no complement code of a column has popcount ≤ 3. -/')
    L.append('def ldata : Data :=')
    L.append(f'  {{ p := p, w := {lower["w"]}, nRows := {lower["nrows"]}, nFree := {lower["nfree"]}, pivEnc := pivEnc, '
             'freeEnc := freeEnc, tbl := tbl, rows := lrows, low := [], gen := lgen }')
    L.append('')
    low = lower['low']
    names = []
    for k in range(0, len(low), low_chunk):
        name = f'low{k // low_chunk}'
        names.append(name)
        L.append(f'def {name} : List (Nat × List Entry) := [')
        body = []
        for idx, entries in low[k:k + low_chunk]:
            body.append(f'  ({idx}, [' + ', '.join(f'({m}, {x})' for m, x in entries) + '])')
        L.append(',\n'.join(body))
        L.append(']')
        L.append('')
    L.append('def lowChunks : List (List (Nat × List Entry)) := [' + ', '.join(names) + ']')
    L.append('')
    L.append('/-- The rows of pivot weight ≤ 3 as sparse lists. -/')
    L.append('def low : List (Nat × List Entry) := lowChunks.flatten')
    L.append('')
    L.append('/-! ## Upper half -/')
    L.append('')
    top = emit_tree(L, upper['packed'], 'utree', 0, named=True)
    L.append(f'def urows : Tree := {top}')
    L.append('')
    L.append('def ugen : List Entry := [' + ', '.join(f'({c}, {x})' for c, x in upper['gen']) + ']')
    L.append('')
    L.append('def udata : UData :=')
    L.append(f'  {{ p := p, w := {upper["w"]}, rows := urows, degEnc := {hx(upper["deg_enc"])}, gen := ugen }}')
    L.append('')
    for n, terms in enumerate(upper['terms']):
        if terms is None:
            continue
        parts = []
        for k in range(0, len(terms), term_chunk):
            name = f'terms{n}_{k // term_chunk}'
            parts.append(name)
            L.append(f'def {name} : List Term := [' +
                     ', '.join(f'({c}, {s}, {i})' for c, s, i in terms[k:k + term_chunk]) + ']')
        L.append(f'def terms{n} : List Term := ' + ' ++ '.join(parts))
        L.append('')
    L.append('/-- The terms of row `n`. -/')
    L.append('def termsOf : Nat → List Term')
    for n, terms in enumerate(upper['terms']):
        if terms is not None:
            L.append(f'  | {n} => terms{n}')
    L.append('  | _ => []')
    L.append('')
    L.append(f'end ICFallDegrees.ExactK6Data{suffix}')
    return '\n'.join(L) + '\n'


def lean_kernel_lower(lower: dict, inst: dict, suffix: str, row_chunk: int, nlow_chunks: int) -> str:
    nrows = lower['nrows']
    L = []
    L.append(f'import ICFallDegrees.ExactK6Data{suffix}')
    L.append('import ICFallDegrees.ExactK6LowerSound')
    L.append('')
    L.append('/-! Generated by `scripts/ic_bounds_exact_k6_gen.py`; do not edit by hand.')
    L.append(f'Kernel checks of the lower certificate (`1 ∉ V_4`) of instance `{inst["tag"]}`.  Every theorem')
    L.append('closed by `decide +kernel` is the Lean kernel evaluating the checker on the data; the checks are')
    L.append('chunked so that one kernel evaluation stays small. -/')
    L.append('')
    L.append(f'namespace ICFallDegrees.ExactK6KernelLower{suffix}')
    L.append(f'open ICFallDegrees.N5Cert ICFallDegrees.N5CertSound ICFallDegrees.ExactK6Data{suffix}')
    L.append('')
    L.append('-- One kernel evaluation at a time.')
    L.append('set_option Elab.async false')
    L.append('')
    L.append('theorem global_ok : globalOK ldata = true := by decide +kernel')
    L.append('theorem lowIdx_nil : N5Cert.lowIdxOK ldata = true := by decide +kernel')
    L.append('theorem gen_ok : member ldata ldata.gen = true := by decide +kernel')
    L.append('theorem free_all : allBelow (freeOK ldata) ldata.nFree = true := by decide +kernel')
    L.append('theorem box_ok : ExactK6Cert.boxOK ldata = true := by decide +kernel')
    L.append('theorem const_ok : ExactK6Cert.constOK ldata = true := by decide +kernel')
    L.append('theorem lowIdx_ok : ExactK6Cert.lowIdxOK ldata low = true := by decide +kernel')
    L.append('')
    bounds = list(range(0, nrows, row_chunk)) + [nrows]
    for i in range(len(bounds) - 1):
        lo, hi = bounds[i], bounds[i + 1]
        L.append(f'theorem rows_{i} : rangeAll (rowOK ldata) {lo} {hi - lo} = true := by decide +kernel')
    L.append('')
    L.append('theorem rows_all : ∀ i < ldata.nRows, rowOK ldata i = true := by')
    L.append('  intro i hi')
    L.append(f'  have hi\' : i < {nrows} := hi')
    for i in range(len(bounds) - 1):
        hi = bounds[i + 1]
        if i < len(bounds) - 2:
            L.append(f'  rcases Nat.lt_or_ge i {hi} with h | h')
            L.append(f'  · exact rangeAll_sound _ _ _ rows_{i} i (by omega) (by omega)')
        else:
            L.append(f'  exact rangeAll_sound _ _ _ rows_{i} i (by omega) (by omega)')
    L.append('')
    for i in range(nlow_chunks):
        L.append(f'theorem lowc_{i} : low{i}.all (ExactK6Cert.lowRowOK ldata) = true := by decide +kernel')
    L.append('')
    L.append('theorem low_all : ∀ r ∈ low, ExactK6Cert.lowRowOK ldata r = true := by')
    L.append('  have H : lowChunks.all (fun l => l.all (ExactK6Cert.lowRowOK ldata)) = true := by')
    names = ', '.join(f'lowc_{i}' for i in range(nlow_chunks))
    L.append(f'    simp only [lowChunks, List.all_cons, List.all_nil, {names}, Bool.true_and]')
    L.append('  intro r hr')
    L.append('  obtain ⟨l, hl, hrl⟩ := List.mem_flatten.mp hr')
    L.append('  exact List.all_eq_true.mp (List.all_eq_true.mp H l hl) r hrl')
    L.append('')
    L.append('/-- The echelon data is valid (`low = []`: the depth-7 specific checks are vacuous). -/')
    L.append('theorem valid : Valid ldata where')
    L.append('  global := global_ok')
    L.append('  rows := rows_all')
    L.append('  free := allBelow_sound _ _ free_all')
    L.append('  lowIdx := lowIdx_nil')
    L.append('  low := fun r hr => absurd hr List.not_mem_nil')
    L.append('  gen := gen_ok')
    L.append('  nRows_lt := by decide')
    L.append('')
    L.append(f'end ICFallDegrees.ExactK6KernelLower{suffix}')
    return '\n'.join(L) + '\n'


def lean_kernel_upper(upper: dict, inst: dict, suffix: str) -> str:
    nrows = len(upper['packed'])
    L = []
    L.append(f'import ICFallDegrees.ExactK6Data{suffix}')
    L.append('import ICFallDegrees.ExactK6UpperSound')
    L.append('')
    L.append('/-! Generated by `scripts/ic_bounds_exact_k6_gen.py`; do not edit by hand.')
    L.append(f'Kernel checks of the upper certificate (`1 ∈ V_5`) of instance `{inst["tag"]}`: one')
    L.append('`decide +kernel` theorem per row of the derivation. -/')
    L.append('')
    L.append(f'namespace ICFallDegrees.ExactK6KernelUpper{suffix}')
    L.append(f'open ICFallDegrees.N5Cert ICFallDegrees.ExactK6Cert ICFallDegrees.ExactK6Data{suffix}')
    L.append('')
    L.append('-- One kernel evaluation at a time.')
    L.append('set_option Elab.async false')
    L.append('')
    L.append('theorem uglobal_ok : uGlobalOK udata = true := by decide +kernel')
    L.append('theorem ugen_ok : genOK udata = true := by decide +kernel')
    L.append(f'theorem ufinal : udata.row {nrows - 1} = 1 := by decide +kernel')
    L.append('')
    for n in range(1, nrows):
        L.append(f'theorem urow_{n} : rowCheck udata {n} terms{n} = true := by decide +kernel')
    L.append('')
    L.append(f'theorem urows_all : ∀ n, 1 ≤ n → n ≤ {nrows - 1} → rowCheck udata n (termsOf n) = true := by')
    L.append('  intro n h1 h2')
    L.append('  interval_cases n')
    for n in range(1, nrows):
        L.append(f'  · exact urow_{n}')
    L.append('')
    L.append(f'end ICFallDegrees.ExactK6KernelUpper{suffix}')
    return '\n'.join(L) + '\n'


def lean_instance(lower: dict, upper: dict, inst: dict, suffix: str) -> str:
    p, a, b, t = inst['p'], inst['a'], inst['b'], inst['t']
    nrows = len(upper['packed'])
    c1, c2, c3 = (-2 * (t * t + a)) % p, (-2 * t) % p, (t * t) % p
    c4, c5 = (-2 * a * t - 4 * b) % p, (a * a - 4 * b * t) % p
    ucodes = [c for c, _ in upper['gen']]
    lcodes = [m for m, _ in lower['gen']]
    L = []
    L.append(f'import ICFallDegrees.ExactK6KernelLower{suffix}')
    L.append(f'import ICFallDegrees.ExactK6KernelUpper{suffix}')
    L.append('import Mathlib.Tactic.NormNum.Prime')
    L.append('')
    L.append('/-! Generated by `scripts/ic_bounds_exact_k6_gen.py`; do not edit by hand.')
    L.append('')
    L.append(f'# Exact last fall degree 5 at depth 6: p = {p}, a = {a}, b = {b}, t = {t}')
    L.append('')
    L.append('For the normalized two-summand chain presentation `PrimeM2General.inputs 6 _ a b t` over')
    L.append(f'`ZMod {p}`, the kernel-checked certificates `ExactK6KernelLower{suffix}` (a containing space of the')
    L.append(f'degree-4 closure that excludes `1`) and `ExactK6KernelUpper{suffix}` (a derivation of `1` in the')
    L.append('degree-5 quotient closure) give `¬ InClosure F 4 1` and `InClosure F 5 1`, hence')
    L.append('`IsLastFallDegree F 5` by `exactLFD_of_first_one`. -/')
    L.append('')
    L.append(f'namespace ICFallDegrees.ExactK6Instance{suffix}')
    L.append('')
    L.append(f'open ICFallDegrees ICFallDegrees.N5Cert ICFallDegrees.N5CertSound ICFallDegrees.ExactK6Data{suffix}')
    L.append('open ICFallDegrees.ChainTranslation ICFallDegrees.ChainQuotient ICFallDegrees.PrimeM2General')
    L.append('open ICFallDegrees.GeneralQuadraticTop ICFallDegrees.ChainClosureBridge MvPolynomial')
    L.append('')
    L.append('noncomputable section')
    L.append('')
    L.append('/-- The prime field of the instance. -/')
    L.append(f'abbrev K := ZMod {p}')
    L.append('')
    L.append(f'instance : Fact (Nat.Prime {p}) := ⟨by norm_num⟩')
    L.append('')
    L.append(f'theorem p_cast : ((p : ℕ) : K) = 0 := ZMod.natCast_self {p}')
    L.append('')
    L.append('theorem two_ne : (2 : K) ≠ 0 := by decide')
    L.append('')
    L.append('theorem hk : 2 ≤ 6 := by norm_num')
    L.append('')
    L.append('theorem var00 : var 2 6 0 ⟨0, by omega⟩ = active 6 hk 0 := Fin.ext (by simp [var_val, active])')
    L.append('theorem var01 : var 2 6 0 ⟨1, by omega⟩ = active 6 hk 1 := Fin.ext (by simp [var_val, active])')
    L.append('theorem var10 : var 2 6 1 ⟨0, by omega⟩ = active 6 hk 2 := Fin.ext rfl')
    L.append('theorem var11 : var 2 6 1 ⟨1, by omega⟩ = active 6 hk 3 := Fin.ext rfl')
    L.append('')
    L.append('/-- The generator list of the upper certificate is the normalized generator. -/')
    L.append('theorem ugen_eq : (ugen.map fun e => (e.2 : K) • ExactK6Upper.M (K := K) e.1).sum =')
    L.append(f'    normalized 6 hk ({a} : K) ({b} : K) ({t} : K) := by')
    L.append('  simp only [ugen, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, ExactK6Upper.M]')
    rws = ', '.join(f'digits_two_bits hk (col {c}) (by decide) (by decide)' for c in ucodes)
    L.append(f'  rw [{rws}]')
    L.append('  simp only [col_zero, col_one, var00, var01, var10, var11]')
    L.append('  norm_num [bit]')
    L.append('  simp only [normalized, top, tail, smul_eq_C_mul, monomial_single_add,')
    L.append('    ← X_pow_eq_monomial, pow_one, mul_one]')
    L.append(f'  have h1 : (-2 * ({t} ^ 2 + {a}) : K) = {c1} := by decide')
    L.append(f'  have h2 : (-2 * {t} : K) = {c2} := by decide')
    L.append(f'  have h3 : ({t} ^ 2 : K) = {c3} := by decide')
    L.append(f'  have h4 : (-2 * {a} * {t} - 4 * {b} : K) = {c4} := by decide')
    L.append(f'  have h5 : ({a} ^ 2 - 4 * {b} * {t} : K) = {c5} := by decide')
    L.append('  rw [h1, h2, h3, h4, h5]')
    L.append('  ring')
    L.append('')
    L.append('/-- The generator list of the lower certificate is the same polynomial. -/')
    L.append('theorem lgen_eq : ExactK6Lower.Ψ K (toLo K ldata.gen) =')
    L.append(f'    normalized 6 hk ({a} : K) ({b} : K) ({t} : K) := by')
    L.append('  rw [← ugen_eq, ExactK6Lower.Ψ_toLo]')
    L.append('  have hdec : ∀ m c : ℕ, (4095 - m) % 64 = c % 128 → (4095 - m) / 64 = c / 128 →')
    L.append('      ExactK6Lower.dec (col m) = col c := by')
    L.append('    intro m c h0 h1')
    L.append('    rw [ExactK6Lower.dec_col]')
    L.append('    funext i')
    L.append('    fin_cases i')
    L.append('    · exact h0')
    L.append('    · exact h1')
    L.append('  simp only [ldata, lgen, ugen, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil,')
    L.append('    ExactK6Upper.M]')
    rws = ', '.join(f'hdec {m} {c} (by decide) (by decide)' for m, c in zip(lcodes, ucodes))
    L.append(f'  rw [{rws}]')
    L.append('')
    L.append(f'/-- **Lower half.** `1` is not in the degree-4 closure. -/')
    L.append(f'theorem not_one_in_closure_four :')
    L.append(f'    ¬ InClosure (Set.range (inputs 6 hk ({a} : K) ({b} : K) ({t} : K))) 4 1 :=')
    L.append(f'  ExactK6Lower.not_one ldata ExactK6KernelLower{suffix}.valid p_cast two_ne low')
    L.append(f'    ExactK6KernelLower{suffix}.box_ok ExactK6KernelLower{suffix}.const_ok')
    L.append(f'    ExactK6KernelLower{suffix}.lowIdx_ok ExactK6KernelLower{suffix}.low_all _ _ _ lgen_eq')
    L.append('')
    L.append('/-- **Upper half.** `1` is in the degree-5 closure. -/')
    L.append('theorem one_in_closure_five :')
    L.append(f'    InClosure (Set.range (inputs 6 hk ({a} : K) ({b} : K) ({t} : K))) 5 1 := by')
    L.append('  have zero := reduced_chain_ideal_trivial (K := K) 2 6 two_ne')
    L.append('  have included := chainInputs_two_subset 6 hk')
    L.append(f'    ({a} : K) ({b} : K) ({t} : K)')
    L.append('  have gdeg : ∀ g ∈ chainInputs K 2 6, g.totalDegree ≤ 5 := by')
    L.append('    rintro g ⟨v, rfl⟩')
    L.append('    exact (chainPolynomial_totalDegree_le 2 6 v).trans (by norm_num)')
    L.append('  rw [one_mem_iff (nonzero := chainInputs_nonzero 2 6) zero _ 5 included gdeg')
    L.append('    (chain_one_reduced (K := K) 2 6)]')
    L.append('  have hin : InQuotClosure (chainInputs K 2 6) (chainInputs_nonzero 2 6)')
    L.append(f'      (Set.range (inputs 6 hk ({a} : K) ({b} : K) ({t} : K))) 5')
    L.append('      (udata.gen.map fun e => (e.2 : K) • ExactK6Upper.M (K := K) e.1).sum := by')
    L.append('    have e : udata.gen = ugen := rfl')
    L.append('    rw [e, ugen_eq, ← ChainModelMap.normalForm_normalized two_ne 6 hk]')
    L.append('    refine InQuotClosure.input ⟨Fin.last _, ?_⟩ ?_')
    L.append('    · simp only [inputs, Fin.lastCases_last]')
    L.append(f'    · have h := (degree_and_top 6 hk ({a} : K) ({b} : K) ({t} : K) (Fin.last _)).1')
    L.append('      simp only [inputs, Fin.lastCases_last] at h')
    L.append('      rw [h]')
    L.append('      norm_num')
    L.append('  have key := ExactK6Upper.derivation_sound two_ne')
    L.append(f'    (Set.range (inputs 6 hk ({a} : K) ({b} : K) ({t} : K))) udata p_cast')
    L.append(f'    ExactK6KernelUpper{suffix}.uglobal_ok ExactK6KernelUpper{suffix}.ugen_ok hin {nrows - 1} termsOf')
    L.append(f'  exact key ExactK6KernelUpper{suffix}.urows_all ExactK6KernelUpper{suffix}.ufinal')
    L.append('')
    L.append(f'/-- **Exact last fall degree 5 at depth 6** for the two-summand prime chain presentation over')
    L.append(f'`𝔽_{p}` with `a = {a}`, `b = {b}`, `t = {t}`. -/')
    L.append('theorem exact_lfd_five :')
    L.append(f'    IsLastFallDegree (Set.range (inputs 6 hk ({a} : K) ({b} : K) ({t} : K))) 5 :=')
    L.append('  exactLFD_of_first_one (by norm_num) one_in_closure_five not_one_in_closure_four')
    L.append('')
    L.append('end')
    L.append('')
    L.append(f'end ICFallDegrees.ExactK6Instance{suffix}')
    return '\n'.join(L) + '\n'


# ----------------------------------------------------------------------------- main

def sha256_text(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--params', required=True, help='p,a,b,t')
    ap.add_argument('--suffix', required=True, help='suffix of the Lean module names, e.g. P1153')
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--lean-dir', default='formal/ICFallDegrees')
    ap.add_argument('--low-chunk', type=int, default=4, help='low rows per kernel theorem')
    ap.add_argument('--row-chunk', type=int, default=150, help='row bound checks per kernel theorem')
    ap.add_argument('--term-chunk', type=int, default=200, help='terms per list literal')
    ap.add_argument('--no-lean', action='store_true', help='build and reference-check only')
    args = ap.parse_args()
    p, a, b, t = (int(x) for x in args.params.split(','))
    assert (p - 1) % 64 == 0 and all(p % q for q in range(2, int(p ** 0.5) + 1)), 'need a prime with 64 | p - 1'
    tag = f'p{p}-a{a}-b{b}-t{t}'
    inst = dict(p=p, a=a, b=b, t=t, tag=tag)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    t0 = time.time()
    zeros = count_zeros(p, a, b, t)
    disc = (4 * a ** 3 + 27 * b * b) % p
    rhs = (t ** 3 + a * t + b) % p
    log(f'instance {tag}: scalar chain solutions |Z| = {zeros}, discriminant {"non" if disc else ""}zero, '
        f't is {"" if pow(rhs, (p - 1) // 2, p) == 1 else "not "}the abscissa of a rational point, '
        f'3t^2 + a = {(3 * t * t + a) % p}')
    assert zeros == 0, 'the ideal is not the unit ideal'
    lower = build_lower(p, a, b, t, log)
    lstats = ref_check_lower(lower, log)
    upper = build_upper(p, a, b, t, log)
    ustats = ref_check_upper(upper, log)
    log(f'built and reference-checked in {time.time() - t0:.0f} s')
    summary = dict(instance=inst, scalar_solutions=zeros,
                   lower=dict(rows=lower['nrows'], free_columns=lower['nfree'], slot_width=lower['w'],
                              profile=lower['profile'], low_rows=len(lower['low']),
                              free_entries=lower['free_entries'], **lstats),
                   upper=dict(rows=len(upper['packed']), slot_width=upper['w'], degrees=upper['degs'],
                              terms_per_row=[len(tm or []) for tm in upper['terms']], **ustats))
    cert = dict(instance=inst,
                lower=dict(p=p, w=lower['w'], nrows=lower['nrows'], nfree=lower['nfree'],
                           piv_codes=lower['piv_codes'], free_codes=lower['free_codes'],
                           gen=lower['gen'], low=lower['low'],
                           packed_sha256=sha256_text(','.join(hx(x) for x in lower['packed']))),
                upper=dict(p=p, w=upper['w'], degs=upper['degs'], gen=upper['gen'], terms=upper['terms'],
                           packed_sha256=sha256_text(','.join(hx(x) for x in upper['packed']))))
    if not args.no_lean:
        lean_dir = Path(args.lean_dir)
        nlow_chunks = (len(lower['low']) + args.low_chunk - 1) // args.low_chunk
        files = {
            f'ExactK6Data{args.suffix}.lean': lean_data(lower, upper, inst, args.suffix, args.low_chunk, args.term_chunk),
            f'ExactK6KernelLower{args.suffix}.lean': lean_kernel_lower(lower, inst, args.suffix, args.row_chunk, nlow_chunks),
            f'ExactK6KernelUpper{args.suffix}.lean': lean_kernel_upper(upper, inst, args.suffix),
            f'ExactK6Instance{args.suffix}.lean': lean_instance(lower, upper, inst, args.suffix),
        }
        summary['lean_sha256'] = {}
        for name, text in files.items():
            (lean_dir / name).write_text(text)
            summary['lean_sha256'][name] = sha256_text(text)
            log(f'wrote {lean_dir / name} ({len(text) / 1e6:.2f} MB)')
    (out / 'certificate.json').write_text(json.dumps(cert, separators=(',', ':')))
    (out / 'summary.json').write_text(json.dumps(summary, indent=1))
    (out / 'generator.log').write_text('\n'.join(lines) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
