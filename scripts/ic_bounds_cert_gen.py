#!/usr/bin/env python3
"""Generate the kernel-checkable N5 base certificate (B1, B2, G) for one instance.

Reads the exported RREF basis of V_4(7) (experiments/ic-bounds-fall/<tag>/k7-D4.npz) or
recomputes it with csd.ic_bounds_fall, builds the certificate data in the encoding the Lean
checker `ICFallDegrees.N5Cert` decodes, re-runs the whole check in Python (a reference
implementation of the Lean checker, entry by entry), and writes

  experiments/ic-bounds-cert/<tag>/certificate.json     the data (hex-encoded big naturals)
  formal/ICFallDegrees/N5CertData.lean                   the same data as Lean literals

Encoding.  A monomial of the depth-7 two-summand algebra is the code m = l0 + 128*l1
(l0, l1 < 128).  Column slots: the RREF has nRows pivot columns and nFree free columns.
  pivEnc  : 16-bit slots, slot idx = pivot code of row idx
  freeEnc : 16-bit slots, slot f = code of free column f
  tbl     : 16-bit slots indexed by code m < 16384: 0 (not a column), idx+1 (pivot of row idx),
            or 1024 + f (free column f)
  rows    : a binary tree (N5Cert.Tree, heap indexing of Tree.get) whose entry idx is row idx
            packed: free column f carries the coefficient in bits [w f, w(f+1)), 512 p^2 < 2^w
  low     : (idx, [(m, c), ...]) for every row whose pivot has weight <= 3: its free-column
            entries as a sparse list, c in [1, p)
  gen     : the generator's entries (m, c) with c in [1, p)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

K0 = 7
D = 4
N0 = 1 << K0
def slot_width(p: int) -> int:
    """Smallest multiple of 4 with 512 p^2 < 2^w (the accumulator bound of N5Cert.member)."""
    w = 4
    while 512 * p * p >= 1 << w:
        w += 4
    return w


def popcount(n: int) -> int:
    return bin(n).count('1')


def code_weight(m: int) -> int:
    return popcount(m % N0) + popcount(m // N0)


def load_basis(npz: Path):
    d = np.load(npz)
    B = d['basis'].astype(np.int64)
    pivots = [int(c) for c in d['pivots']]
    masks = [int(m) for m in d['masks']]
    p = int(d['prime'][0])
    assert int(d['k'][0]) == K0 and int(d['D'][0]) == D
    return B, pivots, masks, p


def generator_entries(inst: dict, p: int) -> list[tuple[int, int]]:
    out = []
    for e1, e2, c in inst['generator_pairs']:
        c %= p
        if c:
            out.append((e1 + N0 * e2, c))
    return sorted(out)


def s3_coefficients(p: int, a: int, b: int, t: int) -> dict[tuple[int, int], int]:
    """Normalized generator g = NF(S_3(x1, x2, t)) as {(e1, e2): coeff}, from the formula in
    docs/proofs/N5-reduction.md Section 1 (PrimeM2General.normalized)."""
    c: dict[tuple[int, int], int] = {}

    def add(e, v):
        c[e] = (c.get(e, 0) + v) % p

    add((2, 2), 1)
    add((1, 1), -2 * (t * t + a))
    add((1, 2), -2 * t)
    add((2, 1), -2 * t)
    add((2, 0), t * t)
    add((0, 2), t * t)
    add((1, 0), -(2 * a * t + 4 * b))
    add((0, 1), -(2 * a * t + 4 * b))
    add((0, 0), a * a - 4 * b * t)
    return {e: v for e, v in c.items() if v}


def build(B: np.ndarray, pivots: list[int], masks: list[int], p: int, gen: list[tuple[int, int]]):
    nrows, ncols = B.shape
    w = slot_width(p)
    piv_set = set(pivots)
    free = [c for c in range(ncols) if c not in piv_set]
    piv_codes = [masks[c] for c in pivots]
    free_codes = [masks[c] for c in free]
    assert len(set(piv_codes)) == nrows and not (set(piv_codes) & set(free_codes))
    assert all(0 <= m < 16384 for m in masks)
    assert nrows < 1024 and len(free_codes) < 16384 - 1024
    tbl = 0
    for i, m in enumerate(piv_codes):
        tbl |= (i + 1) << (16 * m)
    for f, m in enumerate(free_codes):
        tbl |= (1024 + f) << (16 * m)
    piv_enc = 0
    for i, m in enumerate(piv_codes):
        piv_enc |= m << (16 * i)
    free_enc = 0
    for f, m in enumerate(free_codes):
        free_enc |= m << (16 * f)
    packed = []
    sparse = []
    for i in range(nrows):
        row = B[i] % p
        assert row[pivots[i]] == 1
        for j, c in enumerate(pivots):
            assert row[c] == (1 if j == i else 0)
        R = 0
        entries = []
        for f, c in enumerate(free):
            v = int(row[c])
            if v:
                R |= v << (w * f)
                entries.append((masks[c], v))
        packed.append(R)
        sparse.append(entries)
    low_idx = [i for i, m in enumerate(piv_codes) if code_weight(m) <= 3]
    low = [(i, sparse[i]) for i in low_idx]
    return dict(p=p, w=w, nrows=nrows, nfree=len(free_codes), piv_codes=piv_codes, free_codes=free_codes,
                tbl=tbl, piv_enc=piv_enc, free_enc=free_enc, packed=packed, sparse=sparse,
                low_idx=low_idx, low=low, gen=gen)


# ----------------------------------------------------------------------------- tree layout

def tree_index_of_path(i: int):
    """Path (list of 'L'/'R') to index i under N5Cert.Tree.get: 0 = root; i odd -> left with i//2;
    i even > 0 -> right with i//2 - 1."""
    path = []
    while i:
        if i % 2 == 1:
            path.append('L'); i //= 2
        else:
            path.append('R'); i = i // 2 - 1
    return path


def build_tree(values: list[int]):
    """Nested dict tree holding values[i] at the Tree.get position of i."""
    root = {}
    for i, v in enumerate(values):
        node = root
        for step in tree_index_of_path(i):
            node = node.setdefault(step, {})
        node['v'] = v
    return root


def tree_to_lean(node, hx) -> str:
    if node is None:
        return '.leaf'
    v = node.get('v', 0)
    return f"(.node {tree_to_lean(node.get('L'), hx)} {hx(v)} {tree_to_lean(node.get('R'), hx)})"


def tree_get(node, i: int) -> int:
    for step in tree_index_of_path(i):
        node = node.get(step)
        if node is None:
            return 0
    return node.get('v', 0)


# ----------------------------------------------------------------------------- reference checker

def ref_check(cert: dict, log=print) -> dict:
    """Mirror of the Lean checker (N5Cert), on the certificate data; raises on any failure."""
    p, w, nrows, nfree = cert['p'], cert['w'], cert['nrows'], cert['nfree']
    tbl, piv_enc, free_enc = cert['tbl'], cert['piv_enc'], cert['free_enc']
    Bm = 1 << w
    tree = build_tree(cert['packed'])

    def tbl_get(m):
        return (tbl >> (16 * m)) & 0xFFFF

    def piv_get(i):
        return (piv_enc >> (16 * i)) & 0xFFFF

    def free_get(f):
        return (free_enc >> (16 * f)) & 0xFFFF

    def row(i):
        return tree_get(tree, i)

    assert 512 * p * p < Bm and tbl_get(0) >= 1024 and p >= 3, 'globalOK'
    for i in range(nrows):
        assert tbl_get(piv_get(i)) == i + 1 and piv_get(i) < 16384, ('rowOK', i)
        R = row(i)
        for f in range(nfree):
            assert R % Bm < p, ('rowOK slot', i, f)
            R //= Bm
        assert R == 0
    for f in range(nfree):
        assert tbl_get(free_get(f)) >= 1024 and free_get(f) < 16384, ('freeOK', f)
    low_idx = [i for i in range(nrows) if code_weight(piv_get(i)) <= 3]
    assert low_idx == [i for i, _ in cert['low']], 'lowIdxOK'

    def entry_ok(m, x, free_only=False):
        if not (0 < x < p and m < 16384):
            return False
        v = tbl_get(m)
        if v == 0:
            return False
        v -= 1
        if v < 1024:
            return (not free_only) and v < nrows and piv_get(v) == m
        f = v + 1 - 1024
        return f < nfree and free_get(f) == m

    def accum(u):
        R = 0
        for m, x in u:
            v = tbl_get(m)
            if v == 0:
                continue
            v -= 1
            if v < 1024:
                R += x * row(v)
            else:
                R += (p - x) << (w * (v + 1 - 1024))
        return R

    def member(u):
        if len(u) > 512 or not all(entry_ok(m, x) for m, x in u):
            return False
        R = accum(u)
        for f in range(nfree):
            if (R % Bm) % p:
                return False
            R //= Bm
        return R == 0

    def pack_free(L):
        R = 0
        for m, c in L:
            v = tbl_get(m)
            if v:
                R += c << (w * (v - 1024))
        return R

    def shift(entries, i, j, carry):
        out = []
        for m, c in entries:
            l0, l1 = m % 128, m // 128
            if i == 0:
                s = l0 + (1 << j)
                if s < 128:
                    if not carry:
                        out.append((s + 128 * l1, c))
                elif carry:
                    out.append((s - 128 + 128 * l1, c))
            else:
                s = l1 + (1 << j)
                if s < 128:
                    if not carry:
                        out.append((l0 + 128 * s, c))
                elif carry:
                    out.append((l0 + 128 * (s - 128), c))
        return out

    assert member(cert['gen']), 'generator not in S'
    t0 = time.time()
    n_checks = 0
    max_len = 0
    for idx, entries in cert['low']:
        assert idx < nrows and len(entries) <= 256 and all(entry_ok(m, c, True) for m, c in entries)
        assert pack_free(entries) == row(idx), ('low row packing', idx)
        full = [(piv_get(idx), 1)] + entries
        for j in range(K0):
            for blk in range(2):
                for carry in (False, True):
                    u = shift(full, blk, j, carry)
                    max_len = max(max_len, len(u))
                    assert member(u), ('shift', idx, blk, j, carry)
                    n_checks += 1
    log(f'reference checker: {n_checks} shift-membership checks passed in {time.time() - t0:.1f}s; '
        f'longest shifted vector {max_len}')
    return dict(low_rows=len(cert['low']), shift_checks=n_checks, max_shift_len=max_len,
                low_entries=sum(len(e) for _, e in cert['low']))


# ----------------------------------------------------------------------------- Lean emission

def lean_data(cert: dict, tag: str, inst: dict, low_chunk: int) -> str:
    def hx(n: int) -> str:
        return '0x%x' % n if n else '0'
    L = []
    L.append('import ICFallDegrees.N5Cert')
    L.append('')
    L.append('/-! Generated by `scripts/ic_bounds_cert_gen.py`; do not edit by hand.')
    L.append(f'Instance `{tag}`: p = {cert["p"]}, a = {inst["a"]}, b = {inst["b"]}, t = {inst["t"]}.')
    L.append(f'Reduced echelon basis of V_4(7): {cert["nrows"]} rows, {cert["nfree"]} free columns, '
             f'{sum(len(r) for r in cert["sparse"])} free-column entries; {len(cert["low"])} rows with pivot '
             f'weight ≤ 3; slot width {cert["w"]}.')
    L.append('Encoding: `ICFallDegrees.N5Cert`. -/')
    L.append('')
    L.append('namespace ICFallDegrees.N5CertData')
    L.append('open ICFallDegrees.N5Cert')
    L.append('')
    L.append(f'def p : Nat := {cert["p"]}')
    L.append(f'def a : Nat := {inst["a"]}')
    L.append(f'def b : Nat := {inst["b"]}')
    L.append(f'def t : Nat := {inst["t"]}')
    L.append(f'def w : Nat := {cert["w"]}')
    L.append(f'def nRows : Nat := {cert["nrows"]}')
    L.append(f'def nFree : Nat := {cert["nfree"]}')
    L.append(f'def pivEnc : Nat := {hx(cert["piv_enc"])}')
    L.append(f'def freeEnc : Nat := {hx(cert["free_enc"])}')
    L.append(f'def tbl : Nat := {hx(cert["tbl"])}')
    gen = ', '.join(f'({m}, {c})' for m, c in cert['gen'])
    L.append(f'def gen : List Entry := [{gen}]')
    L.append('')
    # packed rows as a tree literal, built from subtrees to keep each term small
    tree = build_tree(cert['packed'])
    L.append('/-- Packed rows (free-column slots of `w` bits), as a `Tree` indexed by `Tree.get`. -/')
    # emit the two depth-1 subtrees ... recursively name subtrees at depth 3
    names = []
    def emit(node, name, depth):
        if node is None:
            return '.leaf'
        if depth == 0:
            L.append(f'def {name} : Tree := {tree_to_lean(node, hx)}')
            L.append('')
            names.append(name)
            return name
        l = emit(node.get('L'), name + 'l', depth - 1)
        r = emit(node.get('R'), name + 'r', depth - 1)
        return f"(.node {l} {hx(node.get('v', 0))} {r})"
    top = emit(tree, 'tree', 3)
    L.append(f'def rows : Tree := {top}')
    L.append('')
    # low rows in chunks
    chunk_names = []
    low = cert['low']
    for s in range(0, len(low), low_chunk):
        name = f'low{s // low_chunk}'
        chunk_names.append(name)
        L.append(f'def {name} : List (Nat × List Entry) := [')
        parts = []
        for idx, entries in low[s:s + low_chunk]:
            es = ', '.join(f'({m}, {c})' for m, c in entries)
            parts.append(f'  ({idx}, [{es}])')
        L.append(',\n'.join(parts))
        L.append(']')
        L.append('')
    L.append('def lowChunks : List (List (Nat × List Entry)) := [' + ', '.join(chunk_names) + ']')
    L.append('')
    L.append('def low : List (Nat × List Entry) := lowChunks.flatten')
    L.append('')
    L.append('def data : Data :=')
    L.append('  { p := p, w := w, nRows := nRows, nFree := nFree, pivEnc := pivEnc, freeEnc := freeEnc,')
    L.append('    tbl := tbl, rows := rows, low := low, gen := gen }')
    L.append('')
    L.append('end ICFallDegrees.N5CertData')
    return '\n'.join(L) + '\n'


def lean_kernel(cert: dict, tag: str, suffix: str, row_chunk: int) -> str:
    """The per-instance kernel-check module: every check as a `decide +kernel` theorem (split into
    chunks so that each kernel evaluation stays within memory), and `valid : Valid data`."""
    nrows = cert['nrows']
    nlow_chunks = (len(cert['low']) + cert['low_chunk'] - 1) // cert['low_chunk']
    L = []
    L.append(f'import ICFallDegrees.N5CertData{suffix}')
    L.append('import ICFallDegrees.N5CertSound')
    L.append('')
    L.append('/-! Generated by `scripts/ic_bounds_cert_gen.py`; do not edit by hand.')
    L.append(f'Kernel checks of the N5 base certificate of instance `{tag}` (data in `N5CertData{suffix}`).')
    L.append('Every theorem below is closed by `decide +kernel`: the Lean kernel evaluates the checker of')
    L.append('`ICFallDegrees.N5Cert` on the data.  The checks are split into chunks so that each kernel')
    L.append('evaluation stays within a few GB of memory.  `valid` assembles them. -/')
    L.append('')
    L.append(f'namespace ICFallDegrees.N5CertKernel{suffix}')
    L.append(f'open ICFallDegrees.N5Cert ICFallDegrees.N5CertSound ICFallDegrees.N5CertData{suffix}')
    L.append('')
    L.append('-- One kernel evaluation at a time: each holds a few GB of kernel cache.')
    L.append('set_option Elab.async false')
    L.append('')
    L.append('theorem global_ok : globalOK data = true := by decide +kernel')
    L.append('theorem lowIdx_ok : lowIdxOK data = true := by decide +kernel')
    L.append('theorem gen_ok : member data data.gen = true := by decide +kernel')
    L.append('theorem free_all : allBelow (freeOK data) data.nFree = true := by decide +kernel')
    L.append('')
    bounds = list(range(0, nrows, row_chunk)) + [nrows]
    for i in range(len(bounds) - 1):
        lo, hi = bounds[i], bounds[i + 1]
        L.append(f'theorem rows_{i} : rangeAll (rowOK data) {lo} {hi - lo} = true := by decide +kernel')
    L.append('')
    L.append('theorem rows_all : ∀ i < data.nRows, rowOK data i = true := by')
    L.append('  intro i hi')
    L.append(f'  have hi\' : i < {nrows} := hi')
    for i in range(len(bounds) - 1):
        lo, hi = bounds[i], bounds[i + 1]
        if i < len(bounds) - 2:
            L.append(f'  rcases Nat.lt_or_ge i {hi} with h | h')
            L.append(f'  · exact rangeAll_sound _ _ _ rows_{i} i (by omega) (by omega)')
        else:
            L.append(f'  exact rangeAll_sound _ _ _ rows_{i} i (by omega) (by omega)')
    L.append('')
    for i in range(nlow_chunks):
        L.append(f'theorem lowc_{i} : low{i}.all (lowRowOK data) = true := by decide +kernel')
    L.append('')
    L.append('theorem low_all : ∀ r ∈ data.low, lowRowOK data r = true := by')
    L.append('  have H : lowChunks.all (fun l => l.all (lowRowOK data)) = true := by')
    names = ', '.join(f'lowc_{i}' for i in range(nlow_chunks))
    L.append(f'    simp only [lowChunks, List.all_cons, List.all_nil, {names}, Bool.true_and]')
    L.append('  intro r hr')
    L.append('  obtain ⟨l, hl, hrl⟩ := List.mem_flatten.mp hr')
    L.append('  exact List.all_eq_true.mp (List.all_eq_true.mp H l hl) r hrl')
    L.append('')
    L.append('/-- All kernel checks pass: the certificate data is valid. -/')
    L.append('theorem valid : Valid data where')
    L.append('  global := global_ok')
    L.append('  rows := rows_all')
    L.append('  free := allBelow_sound _ _ free_all')
    L.append('  lowIdx := lowIdx_ok')
    L.append('  low := low_all')
    L.append('  gen := gen_ok')
    L.append('  nRows_lt := by decide')
    L.append('')
    L.append(f'end ICFallDegrees.N5CertKernel{suffix}')
    return '\n'.join(L) + '\n'


def lean_instance(cert: dict, tag: str, suffix: str, inst: dict) -> str:
    """The per-instance theorem module: (B1), (B2), (G) from `valid`, and the unconditional
    instance theorems `not_one_in_closure` / `lfd_ge_five` via `ChainModelMap`."""
    p, a, b, t = cert['p'], inst['a'], inst['b'], inst['t']
    g = dict(cert['gen'])
    # coefficient identities of the normalized generator, in the order of its codes
    c00, c10, c20, c01, c11, c21, c02, c12 = (g.get(m, 0) for m in (0, 1, 2, 128, 129, 130, 256, 257))
    assert g.get(258) == 1 and c10 == c01 and c20 == c02 and c21 == c12
    L = []
    L.append(f'import ICFallDegrees.N5CertKernel{suffix}')
    L.append('import Mathlib.Tactic.NormNum.Prime')
    L.append('')
    L.append('/-! Generated by `scripts/ic_bounds_cert_gen.py`; do not edit by hand.')
    L.append('')
    L.append(f'# N5 made unconditional for the instance `{tag}`: p = {p}, a = {a}, b = {b}, t = {t}')
    L.append('')
    L.append('The three per-instance hypotheses of `ChainModelMap.not_one_of_base` / `n5_lfd_ge_five` —')
    L.append('(B1) `1 ∉ S`, (B2) `W_8(S)` closed, (G) the generator lies in `S` — are discharged from the')
    L.append(f'kernel-checked certificate `N5CertKernel{suffix}.valid`, with `S` the row space of the')
    L.append('certified reduced echelon basis of `V_4(7)`.  Conclusion: for every depth `k ≥ 8`, `1` is not in')
    L.append('the degree-4 closure of the two-summand presentation (`not_one_in_closure`), and, if the')
    L.append('ideal is the unit ideal, no `e < 5` is an LFD bound (`lfd_ge_five`): `LFD ≥ 5` at every `k ≥ 8`. -/')
    L.append('')
    L.append(f'namespace ICFallDegrees.N5CertInstance{suffix}')
    L.append('')
    L.append(f'open ICFallDegrees ICFallDegrees.N5Cert ICFallDegrees.N5CertSound ICFallDegrees.N5CertData{suffix}')
    L.append(f'open ICFallDegrees.N5CertKernel{suffix} ICFallDegrees.DepthReduction ICFallDegrees.ChainModelMap')
    L.append('open ICFallDegrees.ChainTranslation ICFallDegrees.ChainQuotient ICFallDegrees.PrimeM2General')
    L.append('open ICFallDegrees.GeneralQuadraticTop MvPolynomial')
    L.append('')
    L.append('noncomputable section')
    L.append('')
    L.append(f'/-- The prime field of the instance. -/')
    L.append(f'abbrev K := ZMod {p}')
    L.append('')
    L.append(f'instance : Fact (Nat.Prime {p}) := ⟨by norm_num⟩')
    L.append('')
    L.append(f'theorem p_cast : ((data.p : ℕ) : K) = 0 := ZMod.natCast_self {p}')
    L.append('')
    L.append('theorem two_ne : (2 : K) ≠ 0 := by decide')
    L.append('')
    L.append('/-- (B2): the depth-8 free extension of `S` is closed. -/')
    L.append('theorem base : Closed 7 8 4 (W 7 4 (S K data) 8) := closed valid p_cast')
    L.append('')
    L.append('/-- (B1): the constant is not in `S`. -/')
    L.append('theorem one : Finsupp.single (0 : Fin 2 → ℕ) (1 : K) ∉ S K data := one_not_mem valid')
    L.append('')
    L.append('theorem var00 {k : ℕ} (hk : 2 ≤ k) : var 2 k 0 ⟨0, by omega⟩ = active k hk 0 :=')
    L.append('  Fin.ext (by simp [var_val, active])')
    L.append('theorem var01 {k : ℕ} (hk : 2 ≤ k) : var 2 k 0 ⟨1, by omega⟩ = active k hk 1 :=')
    L.append('  Fin.ext (by simp [var_val, active])')
    L.append('theorem var10 {k : ℕ} (hk : 2 ≤ k) : var 2 k 1 ⟨0, by omega⟩ = active k hk 2 :=')
    L.append('  Fin.ext (by simp [var_val, active])')
    L.append('theorem var11 {k : ℕ} (hk : 2 ≤ k) : var 2 k 1 ⟨1, by omega⟩ = active k hk 3 :=')
    L.append('  Fin.ext (by simp [var_val, active]; omega)')
    L.append('')
    L.append('/-- The generator list of the certificate is the normalized generator, at every depth. -/')
    L.append('theorem Φ_gen (k : ℕ) (hk2 : 2 ≤ k) :')
    L.append('    Φ (K := K) 7 k (embed 0 (toLo K gen)) = normalized k hk2 (a : K) (b : K) (t : K) := by')
    L.append('  simp only [gen, toLo_cons, toLo_nil, map_add, map_smul, map_zero, Φ_embed_eCode]')
    L.append('  rw [digits_two_bits hk2 (col 0) (by decide) (by decide), digits_two_bits hk2 (col 1) (by decide) (by decide),')
    L.append('    digits_two_bits hk2 (col 2) (by decide) (by decide), digits_two_bits hk2 (col 128) (by decide) (by decide),')
    L.append('    digits_two_bits hk2 (col 129) (by decide) (by decide), digits_two_bits hk2 (col 130) (by decide) (by decide),')
    L.append('    digits_two_bits hk2 (col 256) (by decide) (by decide), digits_two_bits hk2 (col 257) (by decide) (by decide),')
    L.append('    digits_two_bits hk2 (col 258) (by decide) (by decide)]')
    L.append('  simp only [col_zero, col_one, var00 hk2, var01 hk2, var10 hk2, var11 hk2]')
    L.append('  norm_num [bit]')
    L.append('  simp only [normalized, top, tail, a, b, t, Nat.cast_ofNat, smul_eq_C_mul, monomial_single_add,')
    L.append('    ← X_pow_eq_monomial, pow_one, mul_one]')
    L.append(f'  have h1 : (-2 * ({t} ^ 2 + {a}) : K) = {c11} := by decide')
    L.append(f'  have h2 : (-2 * {t} : K) = {c21} := by decide')
    L.append(f'  have h3 : ({t} ^ 2 : K) = {c20} := by decide')
    L.append(f'  have h4 : (-2 * {a} * {t} - 4 * {b} : K) = {c10} := by decide')
    L.append(f'  have h5 : ({a} ^ 2 - 4 * {b} * {t} : K) = {c00} := by decide')
    L.append('  rw [h1, h2, h3, h4, h5]')
    L.append('  ring')
    L.append('')
    L.append('theorem gen_support (l : Fin 2 → ℕ) (hl : l ∈ (toLo K gen).support) : LowBox 7 l ∧ lowWt l ≤ 4 := by')
    L.append('  have := toLo_support_subset _ hl')
    L.append('  simp only [gen, List.map_cons, List.map_nil, List.toFinset_cons, List.toFinset_nil,')
    L.append('    Finset.mem_insert, Finset.mem_singleton, Finset.not_mem_empty, or_false] at this')
    L.append('  rcases this with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl <;>')
    L.append('    exact ⟨lowBox_col (by norm_num), by rw [lowWt_col (by norm_num)]; decide⟩')
    L.append('')
    L.append('/-- (G): the generator lies in `S`. -/')
    L.append('theorem generatorInBase (k : ℕ) (hk2 : 2 ≤ k) :')
    L.append('    GeneratorInBase\' 7 4 (S K data) k hk2 (a : K) (b : K) (t : K) := fun _ =>')
    L.append('  ⟨toLo K gen, gen_mem valid p_cast, fun l hl => (gen_support l hl).1, fun l hl => (gen_support l hl).2,')
    L.append('    Φ_gen k hk2⟩')
    L.append('')
    L.append(f'/-- **N5 for the instance, unconditional.** For every depth `k ≥ 8`, `1` is not in the degree-4')
    L.append(f'closure of the two-summand chain presentation over `𝔽_{p}` with `a = {a}`, `b = {b}`, `t = {t}`. -/')
    L.append('theorem not_one_in_closure (k : ℕ) (hk : 8 ≤ k) (hk2 : 2 ≤ k) :')
    L.append(f'    ¬ InClosure (Set.range (inputs k hk2 ({a} : K) ({b} : K) ({t} : K))) 4 1 :=')
    L.append('  ChainModelMap.not_one_of_base 7 4 (S K data) base one k hk hk2 two_ne _ _ _ (generatorInBase k hk2)')
    L.append('')
    L.append('/-- If moreover the ideal is the unit ideal at depth `k`, then no `e < 5` is an LFD bound:')
    L.append('`LFD ≥ 5` at every depth `k ≥ 8`. -/')
    L.append('theorem lfd_ge_five (k : ℕ) (hk : 8 ≤ k) (hk2 : 2 ≤ k)')
    L.append(f'    (unit : (1 : FiniteQuotient.P (2 * k) K) ∈ Ideal.span (Set.range (inputs k hk2 ({a} : K) ({b} : K) ({t} : K)))) :')
    L.append(f'    ∀ e < 5, ¬ HasLFDBound (Set.range (inputs k hk2 ({a} : K) ({b} : K) ({t} : K))) e :=')
    L.append('  ChainModelMap.n5_lfd_ge_five (S K data) base one k hk hk2 two_ne _ _ _ (generatorInBase k hk2) unit')
    L.append('')
    L.append('end')
    L.append('')
    L.append(f'end ICFallDegrees.N5CertInstance{suffix}')
    return '\n'.join(L) + '\n'


def compute_basis(p: int, a: int, b: int, t: int):
    """Recompute the reduced echelon basis of V_4(7) with the repository engine."""
    from csd import ic_bounds_fall as F
    pairs = sorted(((e1, e2, c) for (e1, e2), c in s3_coefficients(p, a, b, t).items()))
    gm = {e1 | (e2 << K0): c for e1, e2, c in pairs}
    grid = F.make_grid(K0, D)
    r = F.closure(gm, grid, p, verify=True, reduced=True)
    assert r['closed'] and not r['contains_one']
    return r['basis'].astype(np.int64), r['pivots'], grid.masks, pairs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tag', default='p40961-s7')
    ap.add_argument('--npz', default=None, help='exported basis (default: experiments/ic-bounds-fall/<tag>/k7-D4.npz)')
    ap.add_argument('--instance', default=None)
    ap.add_argument('--params', default=None, help='p,a,b,t: recompute the basis for this instance instead of loading it')
    ap.add_argument('--lean-out', default=None)
    ap.add_argument('--out-dir', default=None)
    ap.add_argument('--low-chunk', type=int, default=4, help='low rows per data chunk (one kernel theorem each)')
    ap.add_argument('--row-chunk', type=int, default=100, help='rows per rowOK kernel theorem')
    ap.add_argument('--suffix', default='', help='module suffix: N5CertData<suffix>, N5CertKernel<suffix>')
    args = ap.parse_args()
    if args.params:
        p, a, b, t = (int(x) for x in args.params.split(','))
        B, pivots, masks, pairs = compute_basis(p, a, b, t)
        inst = dict(prime=p, a=a, b=b, t=t, generator_pairs=pairs, source='recomputed by csd.ic_bounds_fall.closure')
        tag = args.tag if args.tag != 'p40961-s7' else f'p{p}-a{a}-b{b}-t{t}'
        npz = None
    else:
        tag = args.tag
        fall = ROOT / 'experiments/ic-bounds-fall' / tag
        npz = Path(args.npz) if args.npz else fall / 'k7-D4.npz'
        inst_path = Path(args.instance) if args.instance else fall / 'instance.json'
        inst = json.loads(inst_path.read_text())
        inst['t'] = inst['target'][0]
        B, pivots, masks, p = load_basis(npz)
        assert p == inst['prime']
    gen = generator_entries(inst, p)
    formula = s3_coefficients(p, inst['a'], inst['b'], inst['t'])
    assert sorted((e1 + N0 * e2, c) for (e1, e2), c in formula.items()) == gen, 'generator differs from the normalized S_3 formula'
    cert = build(B, pivots, masks, p, gen)
    stats = ref_check(cert)
    out_dir = Path(args.out_dir) if args.out_dir else ROOT / 'experiments/ic-bounds-cert' / tag
    out_dir.mkdir(parents=True, exist_ok=True)
    cert['low_chunk'] = args.low_chunk
    lean = lean_data(cert, tag, inst, args.low_chunk)
    lean = lean.replace('namespace ICFallDegrees.N5CertData', 'namespace ICFallDegrees.N5CertData' + args.suffix)
    lean = lean.replace('end ICFallDegrees.N5CertData', 'end ICFallDegrees.N5CertData' + args.suffix)
    lean_out = Path(args.lean_out) if args.lean_out else ROOT / f'formal/ICFallDegrees/N5CertData{args.suffix}.lean'
    lean_out.write_text(lean)
    kernel = lean_kernel(cert, tag, args.suffix, args.row_chunk)
    kernel_out = lean_out.with_name(f'N5CertKernel{args.suffix}.lean')
    kernel_out.write_text(kernel)
    instance_lean = lean_instance(cert, tag, args.suffix, inst)
    instance_out = lean_out.with_name(f'N5CertInstance{args.suffix}.lean')
    instance_out.write_text(instance_lean)
    data = dict(tag=tag, prime=p, a=inst['a'], b=inst['b'], t=inst['t'], k0=K0, D=D,
                basis_source=(str(npz.relative_to(ROOT)) if npz else inst['source']),
                npz_sha256=(hashlib.sha256(npz.read_bytes()).hexdigest() if npz else None),
                instance_hash=inst.get('instance_hash'), nrows=cert['nrows'], nfree=cert['nfree'],
                entries=sum(len(r) for r in cert['sparse']), low_rows=len(cert['low']), checks=stats,
                slot_width=cert['w'], low_chunk=args.low_chunk, generator=cert['gen'],
                piv_codes=cert['piv_codes'], free_codes=cert['free_codes'],
                low=[(i, e) for i, e in cert['low']],
                packed_rows_hex=['%x' % n for n in cert['packed']],
                lean_data_file=(str(lean_out.relative_to(ROOT)) if lean_out.is_relative_to(ROOT) else str(lean_out)),
                lean_data_sha256=hashlib.sha256(lean.encode()).hexdigest(),
                lean_kernel_file=(str(kernel_out.relative_to(ROOT)) if kernel_out.is_relative_to(ROOT) else str(kernel_out)),
                lean_kernel_sha256=hashlib.sha256(kernel.encode()).hexdigest(),
                lean_instance_file=(str(instance_out.relative_to(ROOT)) if instance_out.is_relative_to(ROOT) else str(instance_out)),
                lean_instance_sha256=hashlib.sha256(instance_lean.encode()).hexdigest())
    (out_dir / 'certificate.json').write_text(json.dumps(data, indent=1))
    print(json.dumps({k: v for k, v in data.items() if k not in ('packed_rows_hex', 'low', 'piv_codes', 'free_codes')}, indent=1))
    print('lean data:', lean_out, len(lean), 'bytes; kernel file:', kernel_out)


if __name__ == '__main__':
    main()
