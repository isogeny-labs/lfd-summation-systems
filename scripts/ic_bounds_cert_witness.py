#!/usr/bin/env python3
"""N7 witness for the proper-ideal lower bound: a low ideal element of weight <= 4 outside V_4(7).

Input: the N5 base certificate experiments/ic-bounds-cert/<tag>/certificate.json (reduced echelon
basis of S = V_4(7), generator g; written by scripts/ic_bounds_cert_gen.py).

The script looks for an element  u = h * g  of the ideal, where the cofactor
h = sum_e c_e x0^e0 x1^e1 has exponents so small that no monomial of any shifted generator
leaves the low box [0,128)^2 (so the same product is an ideal element at every depth k >= 8,
with no carry out of the low window), such that

  (i)  every monomial of u has binary weight wt(l0) + wt(l1) <= 4, and
  (ii) u is NOT in S (its echelon residual is nonzero in some free column).

Search, in this order, all exact linear algebra over F_p:
  1. single shifts x0^e0 x1^e1 * g with e0, e1 <= 125;
  2. the spaces U_E = span{ x0^e0 x1^e1 g : e0, e1 < E } for E = 6, 14, 30, 62 and then the
     least E (by bisection; U_E is increasing in E) with  U_E cap M_{<=4}  not inside S.
The witness is taken from the least such E, made sparse by a second elimination on the
low-weight columns.

Output:
  experiments/ic-bounds-cert/<tag>/witness.json            search record, cofactor, product
  formal/ICFallDegrees/N5CertWitnessData<P>.lean          the data as Lean literals
  formal/ICFallDegrees/N5CertWitnessInstance<P>.lean      kernel checks and the instance theorem

The Lean checker is ICFallDegrees.N5CertWitness (prodOK, outsideOK); `ref_check` below mirrors
it step by step (same merge, same fuel, same accumulator) and must pass before anything is
written.

  python3 scripts/ic_bounds_cert_witness.py --tag p40961-s7 --suffix P40961
  python3 scripts/ic_bounds_cert_witness.py --tag p11-a2-b5-t7 --suffix P11
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FUEL = 4096


def pc(n: int) -> int:
    return bin(n).count('1')


def code_weight(m: int) -> int:
    return pc(m % 128) + pc(m // 128)


class Base:
    """The certified echelon basis, read from certificate.json."""

    def __init__(self, cert: dict):
        self.p = cert['prime']
        self.w = cert['slot_width']
        self.B = 1 << self.w
        self.piv = cert['piv_codes']
        self.free = cert['free_codes']
        self.rows = [int(h, 16) for h in cert['packed_rows_hex']]
        self.gen = [tuple(x) for x in cert['generator']]
        self.pidx = {m: i for i, m in enumerate(self.piv)}
        self.fidx = {m: f for f, m in enumerate(self.free)}
        assert len(set(self.free)) == len(self.free) and not set(self.free) & set(self.piv)

    def accum(self, u):
        """N5Cert.accum."""
        R = 0
        for m, x in u:
            if m in self.pidx:
                R += x * self.rows[self.pidx[m]]
            else:
                R += (self.p - x) << (self.w * self.fidx[m])
        return R

    def residual(self, u):
        """Free columns f where the echelon residual of u is nonzero: [(f, slot mod p)]."""
        R = self.accum(u)
        out = []
        for f in range(len(self.free)):
            s = (R % self.B) % self.p
            if s:
                out.append((f, s))
            R //= self.B
        assert R == 0
        return out


def rref_mod(M: np.ndarray, p: int, cols) -> tuple[np.ndarray, list[int]]:
    """Row reduce M modulo p with pivots in the given columns (in that order)."""
    M %= p
    r = 0
    piv = []
    for c in cols:
        if r == M.shape[0]:
            break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0:
            continue
        i = r + int(nz[0])
        if i != r:
            M[[r, i]] = M[[i, r]]
        M[r] = (M[r] * pow(int(M[r, c]), p - 2, p)) % p
        col = M[:, c].copy()
        col[r] = 0
        idx = np.nonzero(col)[0]
        if len(idx):
            M[idx] = (M[idx] - np.outer(col[idx], M[r])) % p
        piv.append(c)
        r += 1
    return M, piv


def single_shifts(base: Base) -> dict:
    adm = out = 0
    for e0 in range(126):
        for e1 in range(126):
            u = [(m + e0 + 128 * e1, c) for m, c in base.gen]
            if max(code_weight(m) for m, _ in u) > 4:
                continue
            adm += 1
            out += bool(base.residual(u))
    return dict(admissible=adm, outside_S=out)


def box_kernel(base: Base, E: int, sparse: bool = False):
    """U_E cap M_{<=4}: its dimension and the basis vectors that are outside S.

    Returns (dim, [(terms, u, residual)]) with terms = [((e0, e1), c)] the cofactor and
    u = [(code, coeff)] the product."""
    p = base.p
    shifts = [(e0, e1) for e1 in range(E) for e0 in range(E)]
    mons = sorted({(m % 128 + e0, m // 128 + e1) for e0, e1 in shifts for m, _ in base.gen})
    assert all(x < 128 and y < 128 for x, y in mons)
    high = [mm for mm in mons if pc(mm[0]) + pc(mm[1]) > 4]
    low = [mm for mm in mons if pc(mm[0]) + pc(mm[1]) <= 4]
    cidx = {mm: i for i, mm in enumerate(high + low)}
    nh, nl, ns = len(high), len(low), len(shifts)
    M = np.zeros((ns, nh + nl + ns), dtype=np.int64)
    for r, (e0, e1) in enumerate(shifts):
        for m, c in base.gen:
            M[r, cidx[(m % 128 + e0, m // 128 + e1)]] = c
        M[r, nh + nl + r] = 1
    M, pv = rref_mod(M, p, range(nh))
    Kr = M[len(pv):].copy()
    del M
    assert not Kr[:, :nh].any()
    if sparse:
        # second elimination on the low-weight columns, heaviest exponents last, so that the
        # rows are as short as the kernel allows
        Kr, _ = rref_mod(Kr, p, range(nh, nh + nl))
    out = []
    for r in range(Kr.shape[0]):
        u = [(low[j][0] + 128 * low[j][1], int(Kr[r, nh + j])) for j in np.nonzero(Kr[r, nh:nh + nl])[0]]
        rs = base.residual(u)
        if rs:
            terms = [(shifts[j], int(Kr[r, nh + nl + j])) for j in np.nonzero(Kr[r, nh + nl:])[0]]
            out.append((terms, u, rs))
    return Kr.shape[0], out


# ----------------------------------------------------------------------------- reference checker

def shift_scale(t, L):
    e0, e1, c = t
    return [(m + e0 + 128 * e1, c * x) for m, x in L]


def add_p(n, A, B):
    """N5CertWitness.addP (iterative transcription of the same case analysis)."""
    out = []
    i = j = 0
    while True:
        if n == 0:
            return out + A[i:] + B[j:]
        if i == len(A):
            return out + B[j:]
        if j == len(B):
            return out + A[i:]
        n -= 1
        a, b = A[i], B[j]
        if a[0] < b[0]:
            out.append(a); i += 1
        elif b[0] < a[0]:
            out.append(b); j += 1
        else:
            out.append((a[0], a[1] + b[1])); i += 1; j += 1


def mul_acc(F, g, terms):
    acc = []
    for t in reversed(terms):
        acc = add_p(F, shift_scale(t, g), acc)
    return acc


def normalize(p, L):
    return [(m, x % p) for m, x in L if x % p]


def chunk_terms(terms, size):
    return [terms[i:i + size] for i in range(0, len(terms), size)]


def ref_check(base: Base, chunks, parts, u, slot) -> dict:
    """Mirror of prodOK (one per chunk), sumOK and outsideOK."""
    p = base.p
    for terms, part in zip(chunks, parts):
        for e0, e1, _ in terms:
            assert all(m % 128 + e0 < 128 and m // 128 + e1 < 128 for m, _ in base.gen), 'termOK'
        assert normalize(p, mul_acc(FUEL, base.gen, terms)) == part, 'prodOK'
    acc = []
    for q in reversed(parts):
        acc = add_p(FUEL, q, acc)
    assert normalize(p, acc) == u, 'sumOK'
    for m, x in u:
        assert 0 < x < p and m < 16384 and (m in base.pidx or m in base.fidx), 'entryOK'
        assert pc(m) <= 4, 'weight'
    cs = sum(x for _, x in u)
    assert cs * p < base.B, 'slot bound'
    assert slot < len(base.free)
    assert all(f == slot or base.free[f] != base.free[slot] for f in range(len(base.free)))
    R = base.accum(u)
    val = (R >> (base.w * slot)) % base.B
    assert val % p != 0, 'residual slot'
    return dict(coef_sum=cs, slot_bound=cs * p, slot_value_mod_p=val % p)


def independent_check(base: Base, terms, u):
    """The product recomputed with dictionaries (no merge, no packing)."""
    p = base.p
    prod: dict[tuple[int, int], int] = {}
    for e0, e1, c in terms:
        for m, x in base.gen:
            key = (m % 128 + e0, m // 128 + e1)
            assert key[0] < 128 and key[1] < 128
            prod[key] = (prod.get(key, 0) + c * x) % p
    prod = {k: v for k, v in prod.items() if v}
    assert prod == {(m % 128, m // 128): x for m, x in u}
    assert all(pc(a) + pc(b) <= 4 for a, b in prod)


# ----------------------------------------------------------------------------- Lean emission

def chunked(items, per_line):
    lines = []
    for i in range(0, len(items), per_line):
        lines.append('  ' + ', '.join(items[i:i + per_line]))
    return ',\n'.join(lines)


def lean_data(suffix: str, tag: str, base: Base, chunks, parts, u, slot, E) -> str:
    nterms = sum(len(c) for c in chunks)
    L = [
        'import ICFallDegrees.N5CertWitness',
        '',
        '/-! Generated by `scripts/ic_bounds_cert_witness.py`; do not edit by hand.',
        f'Instance `{tag}` (p = {base.p}): the N7 witness.  `terms_i` are the chunks of a cofactor `h`',
        f'with exponents `< {E}` ({nterms} terms `(e₀, e₁, c)` in {len(chunks)} chunks), `part_i` the product of',
        f'chunk `i` with the generator, `u = h · g` the whole product ({len(u)} monomials, all of weight',
        '≤ 4), and `slot` a free column of the echelon basis of `V_4(7)` where the residual of `u` is',
        'nonzero.  Encoding: `ICFallDegrees.N5CertWitness`. -/',
        '',
        f'namespace ICFallDegrees.N5CertWitnessData{suffix}',
        'open ICFallDegrees.N5Cert ICFallDegrees.N5CertWitness',
        '',
        f'def fuel : Nat := {FUEL}',
        f'def slot : Nat := {slot}',
        '',
    ]
    for i, (terms, part) in enumerate(zip(chunks, parts)):
        L += [f'def terms_{i} : List Term := [', chunked([f'({a}, {b}, {c})' for a, b, c in terms], 8), ']', '',
              f'def part_{i} : List Entry := [', chunked([f'({m}, {x})' for m, x in part], 10), ']', '']
    L += ['def parts : List (List Entry) := [' + ', '.join(f'part_{i}' for i in range(len(parts))) + ']', '',
          'def u : List Entry := [', chunked([f'({m}, {x})' for m, x in u], 10), ']', '',
          f'end ICFallDegrees.N5CertWitnessData{suffix}', '']
    return '\n'.join(L)


def lean_instance(suffix: str, tag: str, p: int, a: int, b: int, t: int, nchunks: int) -> str:
    prods = '\n'.join(f'theorem prod_{i} : prodOK data fuel terms_{i} part_{i} = true := by decide +kernel'
                      for i in range(nchunks))
    allp = ''.join(f'all_cons (N5CertWitness.witness_lowIdeal (d := data) p_cast prod_{i})\n    (' for i in range(nchunks))
    allp += 'all_nil' + ')' * nchunks
    return f'''import ICFallDegrees.N5CertInstance{suffix}
import ICFallDegrees.N5CertWitnessData{suffix}

/-! Generated by `scripts/ic_bounds_cert_witness.py`; do not edit by hand.

# N7 for the instance `{tag}`: p = {p}, a = {a}, b = {b}, t = {t}

`LFD ≥ 5` at every depth `k ≥ 8`, **with no unit-ideal hypothesis**: the kernel checks below
(a chunked product check and a non-membership check)
certify a low ideal element `u = h · g` of weight `≤ 4` that is not in the certified base
`S = V_4(7)`, and `ShiftWitness.lfd_gt_of_lowIdeal` turns it, together with the N5 base
certificate `N5CertKernel{suffix}.valid` (imported, not re-checked), into the lower bound. -/

namespace ICFallDegrees.N5CertWitnessInstance{suffix}

open ICFallDegrees ICFallDegrees.N5Cert ICFallDegrees.N5CertSound ICFallDegrees.N5CertData{suffix}
open ICFallDegrees.N5CertKernel{suffix} ICFallDegrees.N5CertInstance{suffix}
open ICFallDegrees.N5CertWitness ICFallDegrees.N5CertWitnessData{suffix}
open ICFallDegrees.DepthReduction ICFallDegrees.ChainModelMap ICFallDegrees.ShiftWitness
open ICFallDegrees.PrimeM2General

set_option Elab.async false

/-! The product checks, one per chunk of the cofactor: `part_i` is `(∑ terms_i) · gen` modulo
`p`, computed inside the low box. -/
{prods}

/-- The sum check: `u` is the sum of the partial products modulo `p`. -/
theorem sum_ok : sumOK data.p fuel parts u = true := by decide +kernel

/-- The non-membership check: weight `≤ 4` and a nonzero echelon residual in column `slot`. -/
theorem outside_ok : outsideOK data u slot = true := by decide +kernel

noncomputable section

theorem char_p (n : ℕ) (h : (n : K) = 0) : data.p ∣ n :=
  (ZMod.natCast_zmod_eq_zero_iff_dvd n {p}).mp h

theorem parts_lowIdeal : ∀ q ∈ parts, LowIdeal 7 (toLo K gen) (toLo K q) :=
  {allp}

/-- The witness is a low ideal element of the generator. -/
theorem witness_lowIdeal : LowIdeal 7 (toLo K gen) (toLo K u) :=
  sum_lowIdeal p_cast parts_lowIdeal sum_ok

/-- The witness has weight `≤ 4`. -/
theorem witness_weight (l : Fin 2 → ℕ) (hl : l ∈ (toLo K u).support) : lowWt l ≤ 4 :=
  (support_of_outsideOK (d := data) outside_ok l hl).2

/-- The witness is not in the base `S = V_4(7)`. -/
theorem witness_not_mem : toLo K u ∉ S K data :=
  not_mem_of_outsideOK valid p_cast char_p outside_ok

/-- **N7 for the instance, unconditional.**  For every depth `k ≥ 8`, no `e < 5` is an LFD bound
of the two-summand chain presentation over `𝔽_{p}` with `a = {a}`, `b = {b}`, `t = {t}`:
`LFD ≥ 5`.  No hypothesis on the ideal. -/
theorem lfd_ge_five (k : ℕ) (hk : 8 ≤ k) (hk2 : 2 ≤ k) :
    ∀ e < 5, ¬ HasLFDBound (Set.range (inputs k hk2 ({a} : K) ({b} : K) ({t} : K))) e := by
  intro e he
  exact lfd_gt_of_lowIdeal 7 4 (S K data) base k hk hk2 two_ne _ _ _ (toLo K gen)
    (gen_mem valid p_cast) (fun l hl => (gen_support l hl).1) (fun l hl => (gen_support l hl).2)
    (Φ_gen k hk2) (toLo K u) witness_lowIdeal witness_weight witness_not_mem e (by omega)

/-- The same statement with the depth bound only (`2 ≤ k` follows from `8 ≤ k`). -/
theorem lfd_ge_five' :
    ∀ k (hk : 8 ≤ k), ∀ e < 5,
      ¬ HasLFDBound (Set.range (inputs k (by omega) ({a} : K) ({b} : K) ({t} : K))) e :=
  fun k hk => lfd_ge_five k hk (by omega)

end

end ICFallDegrees.N5CertWitnessInstance{suffix}
'''


# ----------------------------------------------------------------------------- main

def search(base: Base, record: dict, t0: float):
    record['single_shifts'] = single_shifts(base)
    print('single shifts x0^e0 x1^e1 g, e <= 125, weight <= 4:', record['single_shifts'], flush=True)

    boxes = {}

    def probe(E):
        if E not in boxes:
            dim, out = box_kernel(base, E)
            boxes[E] = dict(E=E, shifts=E * E, dim_low_weight=dim, basis_vectors_outside_S=len(out))
            print(f'  box E={E}: dim(U_E cap M_<=4) = {dim}, basis vectors outside S: {len(out)} '
                  f'[{time.time() - t0:.0f}s]', flush=True)
        return boxes[E]['basis_vectors_outside_S'] > 0

    for E in (6, 14, 30):
        probe(E)
    assert probe(62), 'no witness with cofactor exponents < 62'
    lo, hi = 30, 62          # invariant: no witness at lo, a witness at hi
    assert not probe(lo)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if probe(mid):
            hi = mid
        else:
            lo = mid
    E = hi
    record['boxes'] = [boxes[e] for e in sorted(boxes)]
    record['least_E'] = E
    print(f'least E with a witness: {E}', flush=True)

    dim, out = box_kernel(base, E, sparse=True)
    assert out
    terms_raw, u, rs = min(out, key=lambda x: (len(x[0]), len(x[1])))
    # terms sorted by increasing shift code: mulAcc then merges every new term at the front
    terms = sorted(((e0, e1, c) for (e0, e1), c in terms_raw), key=lambda x: x[0] + 128 * x[1])
    u = sorted(u)
    record['search_seconds'] = round(time.time() - t0, 1)
    return E, terms, u, rs, rs[0][0]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--tag', default='p40961-s7')
    ap.add_argument('--suffix', default='P40961')
    ap.add_argument('--no-lean', action='store_true')
    ap.add_argument('--chunk', type=int, default=96, help='cofactor terms per kernel check')
    ap.add_argument('--from-record', action='store_true',
                    help='reuse the cofactor of an existing witness.json (skip the search)')
    args = ap.parse_args()

    cdir = ROOT / 'experiments' / 'ic-bounds-cert' / args.tag
    cert_path = cdir / 'certificate.json'
    cert = json.loads(cert_path.read_text())
    base = Base(cert)
    p = base.p
    t0 = time.time()
    record = dict(tag=args.tag, prime=p, a=cert['a'], b=cert['b'], t=cert['t'], k0=7, D=4,
                  certificate_sha256=hashlib.sha256(cert_path.read_bytes()).hexdigest())

    assert not base.residual(base.gen), 'generator not in S'
    old = json.loads((cdir / 'witness.json').read_text()) if args.from_record else None
    if old is not None:
        assert old['certificate_sha256'] == record['certificate_sha256']
        for key in ('single_shifts', 'boxes', 'least_E'):
            record[key] = old[key]
        record['search_seconds'] = old.get('search_seconds', old.get('seconds'))
        E = old['least_E']
        terms = [tuple(x) for x in old['witness']['terms']]
        u = [tuple(x) for x in old['witness']['u']]
        rs = base.residual(u)
        slot = rs[0][0]
    else:
        E, terms, u, rs, slot = search(base, record, t0)
    independent_check(base, terms, u)
    chunks = chunk_terms(terms, args.chunk)
    parts = [normalize(p, mul_acc(FUEL, base.gen, c)) for c in chunks]
    chk = ref_check(base, chunks, parts, u, slot)
    record['witness'] = dict(
        cofactor_terms=len(terms), cofactor_max_exponent=[max(t[0] for t in terms), max(t[1] for t in terms)],
        product_monomials=len(u), product_max_exponent=[max(m % 128 for m, _ in u), max(m // 128 for m, _ in u)],
        residual_free_columns=len(rs), slot=slot, slot_free_code=base.free[slot], fuel=FUEL, checks=chk,
        chunk_size=args.chunk, chunks=len(chunks), part_sizes=[len(q) for q in parts],
        terms=[list(t) for t in terms], u=[list(e) for e in u])
    print(f'witness: {len(terms)} cofactor terms, {len(u)} monomials, residual in {len(rs)} free columns, '
          f'slot {slot} (code {base.free[slot]}); checks {chk}', flush=True)

    if not args.no_lean:
        fdir = ROOT / 'formal' / 'ICFallDegrees'
        files = {
            f'N5CertWitnessData{args.suffix}.lean': lean_data(args.suffix, args.tag, base, chunks, parts, u, slot, E),
            f'N5CertWitnessInstance{args.suffix}.lean': lean_instance(args.suffix, args.tag, p, cert['a'],
                                                                      cert['b'], cert['t'], len(chunks)),
        }
        record['lean'] = {}
        for name, text in files.items():
            (fdir / name).write_text(text)
            record['lean'][f'formal/ICFallDegrees/{name}'] = hashlib.sha256(text.encode()).hexdigest()
            print('wrote', fdir / name, len(text), 'bytes')
    (cdir / 'witness.json').write_text(json.dumps(record, indent=1) + '\n')
    print('wrote', cdir / 'witness.json')


if __name__ == '__main__':
    main()
