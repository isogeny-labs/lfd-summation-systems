"""Structure of the fall spaces of the two-summand prime chain closure (bounds program, item 5).

Reads (or recomputes) the degree-4 closures V_4(k) of one instance for k = 6..kmax and tests the
k-uniform description

    W_k := Free_k^{(k0)}(V_4(k0)) = (+)_{h} tau_h ( V_4(k0) ∩ M_{<= 4 - wt h} ),

the sum over the exponent pairs h supported on levels >= k0 with wt h <= 2 (translations that touch
no level below k0 and hence never carry into the low window).  For k0 = 7:

  * W_k contains every admissible generator row tau_mu g (mu = mu_low + mu_high),
  * 1 ∈ W_k  iff  1 ∈ V_4(7)  (the h = 0 component of 1 is 1),
  * W_k is closed under every admissible variable product for all k >= 8 as soon as W_8 is
    (closure is checked monomial-locally; the only k-dependence is a carry into level k0 or a
    wrap at the top level, both of which stay inside the high patterns),
  * dim W_k = dim V_4(7) + 2(k-7) dim(V_4(7) ∩ M_<=3) + C(2k-14, 2) dim(V_4(7) ∩ M_<=2)
            = 829 + 230 (k-7) + 2 C(2k-14, 2) = 4k^2 + 172k - 571  on generic instances.

The script verifies, for every k in range: V_4(k) == W_k as row spaces, the exclusion certificate
on W_k (contains generator rows, closed, 1 not in W_k) computed independently of V_4(k), and the
stability V_4(k) ∩ M_{levels <= 6} == V_4(7).  It also records why k0 = 6 fails (wrap-around at
level 5: 42 extra dimensions in V_4(6)), the level-set Moebius decomposition of the fall space
V_4 ∩ M_<=3 and of V_4, the birth-round histograms, the explicit weight-2 fall element w (symbolic
formula over Q(a, b, t) checked on the instance), the multiplier weight filtration for a unit-ideal
instance, and the pivot-set stability of V_4(7) across instances.

Everything is exact for the instances and depths computed; nothing here is a bound for other k
or other instances unless stated as a written proof (the reduction lemma above).

Example:
  ulimit -v 3000000; python3 scripts/ic_bounds_fall_analyze.py --prime 40961 --seed 7 --kmax 10
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from csd import elim  # noqa: E402
from csd import ic_bounds_fall as F  # noqa: E402


def load_or_compute(out: Path, pairs, p: int, k: int, product_cap: int, log) -> tuple[dict, F.Grid]:
    grid = F.make_grid(k, 4)
    path = out / f'k{k}-D4.npz'
    if path.exists():
        z = np.load(path)
        assert int(z['prime'][0]) == p and int(z['k'][0]) == k and int(z['D'][0]) == 4
        res = dict(basis=z['basis'].astype(np.int64), pivots=[int(c) for c in z['pivots']],
                   pivot_weights=[int(w) for w in z['pivot_weights']], birth=[int(b) for b in z['birth']],
                   rank=int(z['basis'].shape[0]), history=None, rounds=int(max(z['birth'])), closed=None,
                   contains_one=(len(grid.masks) - 1) in set(int(c) for c in z['pivots']))
        log(f'k={k}: loaded {path.name} ({res["rank"]} rows)')
        return res, grid
    t0 = time.time()
    res = F.closure(F.encode_generator(pairs, k), grid, p, product_cap=product_cap)
    log(f'k={k}: computed V_4, rank {res["rank"]}, closed {res["closed"]}, {time.time() - t0:.1f}s')
    return res, grid


def restrict_to_levels(result: dict, grid: F.Grid, grid_small: F.Grid, p: int) -> np.ndarray:
    """Basis of V ∩ M_{levels < k_small}, written in the coordinates of grid_small."""
    LOW = (1 << grid_small.k) - 1
    keep = F.level_mask_columns(grid, LOW)
    dim, basis = F.intersect_columns(result['basis'], keep, p, want_basis=True)
    emb = np.zeros((dim, len(grid_small.masks)), dtype=np.int64)
    for c in np.flatnonzero(keep):
        a1, a2 = grid.pair(grid.masks[c])
        emb[:, grid_small.column(a1, a2)] = basis[:, c]
    return emb


def moebius_table(basis: np.ndarray, grid: F.Grid, p: int) -> dict:
    dims = {L: F.intersect_columns(basis, F.level_mask_columns(grid, L), p) for L in range(1 << grid.k)}
    e = F.exact_level_counts(dims)
    pieces = {','.join(str(j) for j in range(grid.k) if L >> j & 1): v for L, v in sorted(e.items()) if v}
    # sum of the level-restricted subspaces: is the space the (direct) sum of its pieces?
    maximal = [L for L in e if e[L] and not any(e[M] and M != L and (M & L) == L for M in e)]
    stacked = np.vstack([F.intersect_columns(basis, F.level_mask_columns(grid, L), p, want_basis=True)[1] for L in maximal]) \
        if maximal else np.zeros((0, basis.shape[1]), dtype=np.int64)
    span = elim.fp_echelon(stacked, p)[0] if stacked.shape[0] else 0
    return dict(exact_level_counts=pieces, negative_counts=int(sum(1 for v in e.values() if v < 0)),
                sum_of_counts=int(sum(e.values())), dim=int(basis.shape[0]),
                span_of_level_pieces=int(span), decomposes_into_level_pieces=bool(span == basis.shape[0]))


def symbolic_w():
    """The unique (mod g) combination w = Q g, Q supported on {0,1,2}^2, of weight <= 2, over Q(a,b,t)."""
    import sympy as sp
    a, b, t, x1, x2 = sp.symbols('a b t x1 x2')
    S3 = (x1 - x2) ** 2 * t ** 2 - 2 * ((x1 + x2) * (x1 * x2 + a) + 2 * b) * t + (x1 * x2 - a) ** 2 - 4 * b * (x1 + x2)
    qs = {(i, j): sp.Symbol(f'q{i}{j}') for i in range(3) for j in range(3)}
    Q = sum(qs[i, j] * x1 ** i * x2 ** j for (i, j) in qs)
    P = sp.Poly(sp.expand(Q * S3), x1, x2)
    bad = [(3, 1), (3, 2), (3, 4), (1, 3), (2, 3), (4, 3), (3, 3)]  # exponent pairs in {0..4}^2 of weight >= 3
    eqs = [P.coeff_monomial(x1 ** e1 * x2 ** e2) for e1, e2 in bad]
    sol = sp.linsolve(eqs, list(qs.values()))
    (vec,) = sol
    free = sorted(set().union(*[v.free_symbols for v in vec]) & set(qs.values()), key=str)
    return dict(S3=str(S3), killed_monomials=bad, solution=[str(sp.factor(v)) for v in vec],
                unknown_order=[f'q{i}{j}' for (i, j) in qs], free_parameters=[str(f) for f in free],
                note='two free parameters: q00 (the g component) and one more; the second basis vector is w')


def numeric_w_check(pairs, inst, p: int, grid: F.Grid, res: dict) -> dict:
    """Solve the same 7 linear conditions numerically mod p and compare with V_4 ∩ M_<=2."""
    N = grid.N
    g = F.encode_generator(pairs, grid.k)
    hs = [(i, j) for i in range(3) for j in range(3)]
    bad = [(3, 1), (3, 2), (3, 4), (1, 3), (2, 3), (4, 3), (3, 3)]
    M = np.zeros((len(bad), len(hs)), dtype=np.int64)
    for col, (h1, h2) in enumerate(hs):
        for mask, c in g.items():
            s1, s2 = mask & (N - 1), mask >> grid.k
            if (h1 + s1, h2 + s2) in bad:
                M[bad.index((h1 + s1, h2 + s2)), col] = (M[bad.index((h1 + s1, h2 + s2)), col] + c) % p
    out = []
    rank, piv, _ = elim.fp_echelon(M, p, out)
    R = F.rref(out[0], list(piv), p)
    free = [j for j in range(len(hs)) if j not in piv]
    kernel = []
    for f in free:
        v = np.zeros(len(hs), dtype=np.int64)
        v[f] = 1
        for i, c in enumerate(piv):
            v[c] = (-R[i, f]) % p
        kernel.append(v)
    rows = []
    for v in kernel:
        row = np.zeros(len(grid.masks), dtype=np.int64)
        for col, (h1, h2) in enumerate(hs):
            if v[col]:
                for mask, c in g.items():
                    s1, s2 = mask & (N - 1), mask >> grid.k
                    cc = grid.column(h1 + s1, h2 + s2)
                    row[cc] = (row[cc] + v[col] * c) % p
        rows.append(row)
    rows = np.asarray(rows)
    low2 = F.low_part(res, grid, 2)
    return dict(kernel_dim=len(kernel), max_weight=int(max(grid.weights[np.flatnonzero(r)].max() for r in rows)),
                equals_V4_cap_M2=F.same_row_space(rows, low2, p), free_unknowns=[f'q{hs[f][0]}{hs[f][1]}' for f in free])


def w_formula_check(pairs, inst, p: int, grid: F.Grid, res: dict) -> dict:
    """Evaluate the symbolic multiplier Q_w (q22 = 1, q00 = 0) at (a, b, t) mod p, form w = Q_w g and check
    weight <= 2 and w ∈ V_4.  Denominator a + 3 t^2 must be nonzero mod p."""
    a, b, t = inst['a'] % p, inst['b'] % p, inst['target'][0] % p
    den = (a + 3 * t * t) % p
    if den == 0:
        return dict(skipped='a + 3t^2 = 0 mod p')
    inv = pow(den, -1, p)
    q = {(2, 2): 1, (1, 2): 2 * t % p, (2, 1): 2 * t % p, (1, 1): 2 * (a + 5 * t * t) % p,
         (2, 0): (-t * (11 * a * t + 8 * b + 17 * t ** 3)) * inv % p,
         (0, 2): (-t * (11 * a * t + 8 * b + 17 * t ** 3)) * inv % p,
         (1, 0): 2 * (5 * a * a * t + 2 * a * b + 15 * a * t ** 3 - 2 * b * t * t + 16 * t ** 5) * inv % p,
         (0, 1): 2 * (5 * a * a * t + 2 * a * b + 15 * a * t ** 3 - 2 * b * t * t + 16 * t ** 5) * inv % p}
    g = F.encode_generator(pairs, grid.k)
    N = grid.N
    row = np.zeros(len(grid.masks), dtype=np.int64)
    for (h1, h2), qc in q.items():
        for mask, c in g.items():
            s1, s2 = mask & (N - 1), mask >> grid.k
            col = grid.column(h1 + s1, h2 + s2)
            row[col] = (row[col] + qc * c) % p
    nz = np.flatnonzero(row)
    weight = int(grid.weights[nz].max())
    out = []
    rank, piv, _ = elim.fp_echelon(res['basis'], p, out)
    in_V4 = elim.fp_echelon(np.vstack([res['basis'], row[None, :]]), p)[0] == rank
    low2 = F.low_part(res, grid, 2)
    gen_row = F.generator_rows(g, grid, p, max_shift_weight=0)[0]
    spans = F.same_row_space(np.vstack([row[None, :], gen_row]), low2, p)
    return dict(weight=weight, in_V4=bool(in_V4), support=int(nz.size), g_and_w_span_V4_cap_M2=spans['equal'],
                denominator_a_plus_3t2=den)


def multiplier_filtration(pairs, p: int, k: int, res: dict, grid: F.Grid, max_rows: int = 7000) -> dict:
    """dim(Q ∩ B_w) for Q = V_4 g^{-1} (unit ideal) and B_w the Dirac masses of weight <= w;
    the same for the fall part; and kappa(w) = dim{q in B_w : wt(q g) <= 4}."""
    g = F.encode_generator(pairs, k)
    N = 1 << k
    Q = F.multipliers(res['basis'], g, grid, p).reshape(res['rank'], -1)
    wt = F.weight_table(k).reshape(-1)
    low = [i for i, w in enumerate(res['pivot_weights']) if w <= 3]
    f, f3, kappa4, kappa3 = [], [], [], []
    for w in range(2 * k + 1):
        other = np.flatnonzero(wt > w)
        f.append(int(res['rank'] - (elim.fp_echelon(Q[:, other], p)[0] if other.size else 0)))
        f3.append(int(len(low) - (elim.fp_echelon(Q[low][:, other], p)[0] if other.size else 0)))
    for bound, store in ((4, kappa4), (3, kappa3)):
        for w in range(2 * k + 1):
            hs = [h for h in range(N * N) if wt[h] <= w]
            if len(hs) > max_rows:
                break  # kappa is recorded only while the Dirac matrix fits the memory budget
            cols: dict[int, int] = {}
            entries = []
            for h in hs:
                h1, h2 = divmod(h, N)
                row = {}
                for mask, c in g.items():
                    s1, s2 = mask & (N - 1), mask >> k
                    a1, a2 = (h1 + s1) % N, (h2 + s2) % N
                    if int(a1).bit_count() + int(a2).bit_count() > bound:
                        col = cols.setdefault(a1 * N + a2, len(cols))
                        row[col] = (row.get(col, 0) + c) % p
                entries.append(row)
            M = np.zeros((len(hs), max(1, len(cols))), dtype=np.int64)
            for i, row in enumerate(entries):
                for j, c in row.items():
                    M[i, j] = c
            store.append(int(len(hs) - elim.fp_echelon(M, p)[0]))
    return dict(dim_Q_cap_Bw=f, dim_Qfall_cap_Bw=f3, kappa4=kappa4, kappa3=kappa3,
                kappa_note=f'kappa recorded for w = 0..{len(kappa4) - 1} (Dirac matrix capped at {max_rows} rows)',
                max_multiplier_weight=int(next(w for w in range(2 * k + 1) if f[w] == res['rank'])),
                note='Q ∩ B_w = kappa4(w) for small w means: every weight-<=4 combination of translates of g by '
                     'weight <= w lies in V_4; the first w where they differ is where V_4 stops being "everything reachable"')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--prime', type=int, default=40961)
    ap.add_argument('--seed', type=int, default=7)
    ap.add_argument('--kmax', type=int, default=10)
    ap.add_argument('--k0', type=int, default=7)
    ap.add_argument('--product-cap', type=int, default=20000)
    ap.add_argument('--moebius-kmax', type=int, default=8)
    ap.add_argument('--unit-prime', type=int, default=1073153, help='prime for the unit-ideal multiplier analysis (0 to skip)')
    ap.add_argument('--unit-seed', type=int, default=1)
    ap.add_argument('--unit-kmax', type=int, default=7)
    ap.add_argument('--extra-instances', type=str, default='12289:3,12289:j1728',
                    help='comma list prime:seed or prime:j1728 (y^2=x^3+x with random target) for cross-instance checks')
    ap.add_argument('--out', type=Path, default=ROOT / 'experiments/ic-bounds-fall')
    ap.add_argument('--analysis-name', default='analysis.json')
    ap.add_argument('--skip-moebius', action='store_true')
    ap.add_argument('--skip-w', action='store_true')
    args = ap.parse_args()
    p, k0 = args.prime, args.k0
    out = args.out / f'p{p}-s{args.seed}'
    out.mkdir(parents=True, exist_ok=True)
    report: dict = dict(prime=p, seed=args.seed, k0=k0, kmax=args.kmax)
    logs = []

    def log(msg):
        print(msg, flush=True)
        logs.append(msg)

    inst_path = out / 'instance.json'
    if inst_path.exists():
        inst = json.loads(inst_path.read_text())
    else:
        inst = F.random_instance(p, 6, args.seed)
        inst['instance_hash'] = F.instance_hash(inst)
        inst_path.write_text(json.dumps(inst, indent=1))
    pairs = inst['generator_pairs']
    report['instance'] = {a: inst[a] for a in ('prime', 'a', 'b', 'target', 'target_order', 'curve_order', 'instance_hash')}

    results, grids = {}, {}
    for k in range(6, args.kmax + 1):
        results[k], grids[k] = load_or_compute(out, pairs, p, k, args.product_cap if k <= 10 else 6000, log)
    report['dimensions'] = {k: dict(dim_V4=results[k]['rank'], low_profile=F.low_profile(results[k], grids[k]),
                                    formula=4 * k * k + 172 * k - 571, formula_low=4 * k + 87,
                                    birth_histogram=dict(Counter(results[k]['birth'])),
                                    low_birth_histogram=dict(Counter(b for b, w in zip(results[k]['birth'], results[k]['pivot_weights']) if w <= 3)))
                            for k in results}

    # --- free extension from k0 ---------------------------------------------------------
    base, gbase = results[k0], grids[k0]
    base_profile = F.low_profile(base, gbase)
    report['base'] = dict(k0=k0, dim=base['rank'], low_profile=base_profile, contains_one=base['contains_one'],
                          basis_sha256=F.row_digest(base['basis']))
    free_checks = {}
    for k in range(k0, args.kmax + 1):
        g = F.encode_generator(pairs, k)
        rows, info = F.free_extension(base, gbase, grids[k], p)
        cmp = F.same_row_space(rows, results[k]['basis'], p)
        t0 = time.time()
        cert = F.check_containing_space(rows, g, grids[k], p, product_cap=6000)
        predicted = base_profile[4] + 2 * (k - k0) * base_profile[3] + math.comb(2 * (k - k0), 2) * base_profile[2] \
            + math.comb(2 * (k - k0), 3) * base_profile[1] + math.comb(2 * (k - k0), 4) * base_profile[0]
        stable = None
        if k > k0:
            emb = restrict_to_levels(results[k], grids[k], gbase, p)
            stable = F.same_row_space(emb, base['basis'], p)
        free_checks[k] = dict(rows=int(rows.shape[0]), predicted_dim=predicted, rows_by_high_weight=info['rows_by_high_weight'],
                              equals_V4=cmp, certificate=cert, stable_low_part_equals_base=stable,
                              certificate_seconds=round(time.time() - t0, 1))
        log(f'k={k}: Free_k^({k0})(V_4({k0})) dim {cmp["rank_a"]} (predicted {predicted}), V_4 dim {cmp["rank_b"]}, '
            f'equal {cmp["equal"]}, certificate {cert["certifies_one_not_in_VD"]}, stable {None if stable is None else stable["equal"]}')
    report['free_extension_from_k0'] = free_checks

    # --- why k0 = 6 fails ----------------------------------------------------------------
    if 6 in results and k0 > 6:
        g6, grid6 = results[6], grids[6]
        fails = {}
        for k in range(7, min(args.kmax, 9) + 1):
            rows, info = F.free_extension(g6, grid6, grids[k], p)
            cmp = F.same_row_space(rows, results[k]['basis'], p)
            emb = restrict_to_levels(results[k], grids[k], grid6, p)
            st = F.same_row_space(emb, g6['basis'], p)
            cert = F.check_containing_space(rows, F.encode_generator(pairs, k), grids[k], p, product_cap=6000)
            fails[k] = dict(free_dim=cmp['rank_a'], V4_dim=cmp['rank_b'], sum_dim=cmp['rank_sum'],
                            dims_only_in_free=cmp['rank_sum'] - cmp['rank_b'], dims_only_in_V4=cmp['rank_sum'] - cmp['rank_a'],
                            stable_low_part_dim=st['rank_a'], stable_low_part_in_V4_6=st['a_in_b'],
                            V4_6_minus_stable=st['rank_b'] - st['rank_a'], certificate=cert)
            log(f'k0=6 -> k={k}: free dim {cmp["rank_a"]}, V_4 dim {cmp["rank_b"]}, sum {cmp["rank_sum"]}; '
                f'V_4(k) ∩ M_lev<=5 has dim {st["rank_a"]} ⊂ V_4(6) ({st["a_in_b"]}), closed {cert["closed"]}')
        report['free_extension_from_k0_6_fails'] = fails
        # which weight-<=3 elements of V_4(6) touch level 5 (wrap-prone) versus V_4(7)
        for k, name in ((6, 'V4_6'), (7, 'V4_7')):
            if k in results:
                gr = grids[k]
                top = np.asarray([bool((gr.pair(m)[0] >> (k - 1)) & 1 or (gr.pair(m)[1] >> (k - 1)) & 1) for m in gr.masks])
                low = F.low_part(results[k], gr, 3)
                report[f'{name}_fall_space_rank_on_top_level_columns'] = int(elim.fp_echelon(low[:, top], p)[0])

    # --- level-set decomposition ----------------------------------------------------------
    moeb = {}
    for k in ([] if args.skip_moebius else range(6, min(args.kmax, args.moebius_kmax) + 1)):
        low = F.low_part(results[k], grids[k], 3)
        moeb[k] = dict(fall_space=moebius_table(low, grids[k], p), V4=moebius_table(results[k]['basis'], grids[k], p))
        log(f'k={k}: fall-space level pieces {moeb[k]["fall_space"]["exact_level_counts"]} '
            f'(decomposes: {moeb[k]["fall_space"]["decomposes_into_level_pieces"]}); V_4 negative counts {moeb[k]["V4"]["negative_counts"]}')
    report['level_decomposition'] = moeb
    # the fall space as base pieces plus free translates of g and w
    fall_structure = {}
    for k in range(k0 + 1, args.kmax + 1):
        gk = grids[k]
        g = F.encode_generator(pairs, k)
        low_base = F.low_part(base, gbase, 3)
        low2 = F.low_part(base, gbase, 2)
        parts = [F.embed_rows(low_base, gbase, gk)]
        for h in F.high_patterns(k, k0, 1)[1:]:
            parts.append(F.embed_rows(low2, gbase, gk, h))
        rows = np.vstack(parts)
        cmp = F.same_row_space(rows, F.low_part(results[k], gk, 3), p)
        fall_structure[k] = dict(predicted=int(low_base.shape[0] + 2 * (k - k0) * low2.shape[0]), compare=cmp)
    report['fall_space_structure'] = dict(description='F_4(k) = F_4(k0) (+) sum_{j>=k0, i} tau_{2^j e_i} span(g, w)',
                                          checks=fall_structure)

    # --- explicit w ----------------------------------------------------------------------
    if not args.skip_w:
        try:
            report['w_symbolic'] = symbolic_w()
        except Exception as exc:  # pragma: no cover
            report['w_symbolic'] = dict(error=repr(exc))
        report['w_numeric'] = {k: numeric_w_check(pairs, inst, p, grids[k], results[k]) for k in (6, 7) if k in results}
        report['w_formula_on_instance'] = {k: w_formula_check(pairs, inst, p, grids[k], results[k]) for k in (6, 7) if k in results}
        log(f'w check: {report["w_numeric"]}; formula check {report["w_formula_on_instance"]}')

    # --- multiplier filtration on a unit-ideal instance ---------------------------------
    if args.unit_prime:
        up = args.unit_prime
        uinst = F.random_instance(up, 6, args.unit_seed)
        mf = dict(instance={a: uinst[a] for a in ('prime', 'a', 'b', 'target', 'target_order')}, solutions_k6=len(uinst['points_at_k']))
        for k in range(6, args.unit_kmax + 1):
            gr = F.make_grid(k, 4)
            r = F.closure(F.encode_generator(uinst['generator_pairs'], k), gr, up)
            try:
                mf[k] = multiplier_filtration(uinst['generator_pairs'], up, k, r, gr)
                mf[k]['dim_V4'] = r['rank']
                log(f'unit instance k={k}: dim {r["rank"]}, multiplier filtration {mf[k]["dim_Q_cap_Bw"]}, max weight {mf[k]["max_multiplier_weight"]}')
            except ValueError as exc:
                mf[k] = dict(skipped=str(exc))
        report['multiplier_filtration'] = mf

    # --- cross-instance checks ------------------------------------------------------------
    cross = {}
    base_pivots = set(base['pivots'])
    for spec in [s for s in args.extra_instances.split(',') if s]:
        pr, sd = spec.split(':')
        pr = int(pr)
        if sd == 'j1728':
            import random as _r
            E = F.I.PrimeCurve(pr, 1, 0)
            rng = _r.Random(10)
            R = F.S.random_rational_target(E, rng, pr)
            pairs2 = sorted(((e1, e2, c) for (e1, e2), c in F.I.s3_prime_coeffs(E, R[0]).items()))
            desc = dict(prime=pr, a=1, b=0, target=list(R), note='j=1728 curve y^2=x^3+x, random target (depth-10 probe family)')
        else:
            inst2 = F.random_instance(pr, 6, int(sd))
            pairs2 = inst2['generator_pairs']
            desc = {a: inst2[a] for a in ('prime', 'a', 'b', 'target', 'target_order')}
        res2 = {k: F.closure(F.encode_generator(pairs2, k), F.make_grid(k, 4), pr) for k in range(6, 10)}
        b2, gb2 = res2[k0], F.make_grid(k0, 4)
        checks = {}
        for k in range(k0 + 1, 10):
            gk = F.make_grid(k, 4)
            rows, info = F.free_extension(b2, gb2, gk, pr)
            cmp = F.same_row_space(rows, res2[k]['basis'], pr)
            cert = F.check_containing_space(rows, F.encode_generator(pairs2, k), gk, pr, product_cap=6000)
            checks[k] = dict(dim_V4=res2[k]['rank'], formula=4 * k * k + 172 * k - 571, equals_free=cmp['equal'],
                             certificate=cert['certifies_one_not_in_VD'])
        cross[spec] = dict(instance=desc, dims={k: res2[k]['rank'] for k in res2},
                           low_profiles={k: F.low_profile(res2[k], F.make_grid(k, 4)) for k in res2},
                           base_pivot_set_equal=(set(b2['pivots']) == base_pivots),
                           base_pivot_symmetric_difference=len(set(b2['pivots']) ^ base_pivots), free_checks=checks)
        log(f'instance {spec}: dims {cross[spec]["dims"]}, pivot set equal {cross[spec]["base_pivot_set_equal"]}, free checks {[(k, c["equals_free"], c["certificate"]) for k, c in checks.items()]}')
    report['cross_instances'] = cross
    report['log'] = logs
    (out / args.analysis_name).write_text(json.dumps(report, indent=1, default=str))
    log(f'wrote {out / args.analysis_name}')


if __name__ == '__main__':
    main()
