# Item 1: exact LFD of three-summand (m = 3) prime subgroup-chain systems

> **Note (2026-11-12).** Some runs listed below as not run were run later and are reported in `paper/curve-lfd-bounds.tex`, including the `k = 6` degree-7 closure (`experiments/ic-bounds-m3-depth/`) and the three-to-eight-summand runs (`experiments/ic-bounds-msum/`). The lists below record the state when this file was written.

Status: complete for k = 1, 2, 3, 4. 461 seeded cases (5 primes x 6 curves x up to 4 targets x
4 depths); every case has an exact certified LFD in all three encodings (direct m = 3, split m = 3
with a free auxiliary, m = 2 reference): 1383 exact values, no capped case. Seed 20261107.

## Systems

Curve `E: y^2 = x^3 + a x + b` over `F_p`, `p` prime with `2^4 | p - 1`; `H` the subgroup of order
`2^k` of `F_p^*`, `H = {g0^e}` with `g0 = csd.ic_curve_sweep.subgroup_generator_value(p, k)`.
Each summand `x_i` (`i = 1, 2, 3`) is written through the chain `u_{i,0} = x_i`,
`u_{i,j}^2 = u_{i,j+1}`, `u_{i,k-1}^2 = 1`. Modulo the chains a monomial is a squarefree mask (block
`i` holds the bits of the exponent of `x_i` mod `2^k`); multiplication by `u_{i,j}` adds `2^j` to
that exponent mod `2^k`. The target `R = (t, y_R)` is fixed.

* **direct**: inputs = the `3k` chain equations and `g = NF(S_4(x1, x2, x3, t))`, where
  `S_4(x1,x2,x3,x4) = Res_X(S_3(x1,x2,X), S_3(x3,x4,X))` and
  `S_3(u,v,w) = (u-v)^2 w^2 - 2((u+v)(uv+a)+2b) w + (uv-a)^2 - 4b(u+v)`. `S_4` was derived here
  (resultant of two quadratics, `csd/ic_bounds_m3.py:s4_polynomial`) and checked: it vanishes at the
  x-coordinates of 347 random rational decompositions `P1 + P2 + P3 = R` over six primes (17 to
  1000003), is symmetric in its four arguments, has degree 4 in each variable and total degree 12,
  125 terms generically (`selfcheck.json`, `symbolic.json`). Every rational decomposition with
  `x_i in H` found in the sweep lies in the scalar zero set `Z` (asserted per case).
  Input degree = degree of the normal form: `deg NF(S_4) <= 6` (exponents 0..4 reduce to masks of
  weight <= 2 per block; weight 2 only for exponent 3), `<= 3` at `k = 1`.
* **split** (left chain): inputs = the `3k` chains and
  `g1 = NF(S_3(x1, x2, z))`, `g2 = NF(S_3(z, x3, t))`, with ONE auxiliary coordinate `z`
  (the x-coordinate of `P1 + P2`). `z` is an ordinary polynomial variable: no chain, no subgroup
  restriction, never reduced (this is NOT `csd.ic_research_campaign.prime_split`, which writes the
  auxiliary as a chain variable in `H`). Ring `F_p[u_{1..3,*}, z]`, `3k + 1` variables;
  `deg g1 = 4`, `deg g2 = 3` for `k >= 2`. Two ideals are certified from the same closures:
  - `I_free = (chains, g1, g2)` - `z` ranges over the algebraic closure (the encoding used for the
    headline numbers);
  - `I_field = (chains, g1, g2, z^p - z)` - `z` ranges over the full field `F_p`. The field
    equation has degree `p`, so `V_D(F_field) = V_D(F_free)` for every `D < p`; all closures here
    have `D <= 12 < 17 <= p`.
  Why this encoding: it is the standard decomposition `S_4 = Res_z(S_3(x1,x2,z), S_3(z,x3,t))`
  with the intermediate sum kept as a free coordinate over the full field, the way relation
  collection with `m >= 3` splits the summation polynomial; restricting `z` to a chain would
  change the system (relations whose partial sum is outside `H` would be lost).
* **m2** (reference): `{chains, NF(S_3(x1, x2, t))}` at the same `(p, curve, t, k)`, computed with
  the existing `csd.ic_curve_sweep.prime_local_profile`.

## Certificates (what makes a value exact)

`V_D` is computed exactly modulo the chains (equal to the ordinary closure for `D >= 2`; the
chains are a Groebner basis with pairwise coprime leading terms `u_{i,j}^2`):

* `direct`, `k <= 3`: `csd.ic_fall_certificates.chain_degree_closure(m=3)`, every closure not
  containing 1 re-checked with `check_closed_exclusion_certificate` (containing, multiplication
  stable, excludes 1).
* `split`, `k <= 3`: `csd.ic_bounds_m3.mixed_closure` (dense; the closure is verified closed: all
  variable shifts of all rows of degree `< D` re-reduce to zero).
* `k = 4` (both encodings) and any dense budget overflow: `csd.ic_bounds_m3_batched.batched_closure`
  (work queue + reduced echelon basis, exact float64 products with `n (p-1)^2 < 2^53`; final
  verification pass: every variable shift of every basis row of degree `< D` reduces to zero).
  When the exact `dim(I cap M_<=D)` is known (evaluation rank of `M_<=D` on the quotient algebra)
  the saturation stops as soon as the span reaches it: span `subset V_D subset I cap M_<=D` with
  equal dimensions gives `V_D` exactly (`closure_check = rank-equals-dim(I cap M_<=D)` in the
  results). It reproduces `chain_degree_closure` / `mixed_closure` rank and degree profile on 180
  (system, D) pairs and 20 full exact profiles, k = 1..3, both encodings (`selfcheck.json`).

Exact value:

* **unit**: `exact_lfd` = least `D` with `1 in V_D` (all smaller closures exclude 1).
* **proper (local-to-global)**: at the first `D` where the monomials of degree `< D` span `R/I`
  and `dim V_D/chains = |M_<=D| - dim R/I`, the exact LFD is the last strict fall
  (`dim(V_D cap M_<=D-1) > dim V_{D-1}`) at or below `D`. `dim R/I`: direct `|Z|` with `Z`
  enumerated over `H^3`; split `sum_{x in H^3} deg gcd(g1(x,z), g2(x,z))` (and with `z^p - z` for
  `I_field`), spanning checked by the rank of the images in `(+)_x F_p[z]/(h_x)`.
* **capped**: only the last strict fall of exact closures is reported as a lower bound; a stall is
  never used as a bound.

**Proof of the local-to-global criterion (general form; no radicality or rationality needed).**
Let `A = R/I` be finite dimensional, let the images of `M_<D` span `A`, and let
`dim V_D = |M_<=D| - dim A`. Since `V_D subset I cap M_<=D` and `M_<=D -> A` is onto (it contains
`M_<D`), `dim(I cap M_<=D) = |M_<=D| - dim A`, so `V_D = I cap M_<=D`. Induction on `D' >= D`
with hypothesis `V_D' = I cap M_<=D'`: every monomial `mu` of degree `<= D'` is congruent mod
`V_D'` to some `rho(mu)` in `span M_<D` (because `M_<D` spans `A` and the difference lies in
`I cap M_<=D'`). A monomial of degree `D'+1` is `x mu` with `deg mu = D'`;
`x mu - x rho(mu) = x (mu - rho(mu))` is an admissible product in `V_{D'+1}`, and `x rho(mu)` has
degree `<= D <= D'`, so it is congruent mod `V_D' subset V_{D'+1}` to an element of `span M_<D`.
Hence every `f in I cap M_<=D'+1` is congruent mod `V_{D'+1}` to some `rho in span M_<D`, which
lies in `I cap M_<D subset V_D subset V_{D'+1}`; so `f in V_{D'+1}`. Therefore no strict fall
occurs above `D`, every ideal element `h` lies in `V_{max(D, deg h)}`, and the least such bound is
the last strict fall (for `e < D`: `I cap M_<=e = V_D cap M_<=e`, which equals `V_e` exactly when
there is no fall in `(e, D]`).

## Proven statements (written proofs; not Lean-checked)

**Proposition 1 (carry lower bound; proof).** Let `F = chains + {g}` with `g` a chain normal form of
degree `d >= 2` used as an input of degree `d`, over any number `m` of blocks. Suppose a chain
variable `u = u_{i,j}` satisfies (a) `NF(u g)` has degree `<= d`, and (b) `g`, as a function on
`H^m`, is not supported inside one level set `{x_i^(2^j) = c}`. Then `LFD(F) >= d + 1`.

*Proof.* Let `W = (C cap M_<=d) + F_p g`, `C` the chain ideal. The chains form a Groebner basis for
any degree order (leading terms `u_{i,j}^2`, pairwise coprime), so every `c in C` has a standard
representation and its top-degree form lies in the monomial ideal `(u_{i,j}^2)`: it contains no
squarefree monomial. The top form of `g` consists of squarefree monomials. Hence `c + lambda g`
with `lambda != 0` has degree exactly `d` when it lies in `M_<=d`, so `W cap M_<d = C cap M_<d`,
whose variable multiples stay in `C cap M_<=d subset W`. `W` contains every input of degree `<= d`
(chains, `g`), so `V_d subset W` by minimality. Now `h = NF(u g)` lies in `I` and has degree
`<= d` by (a). If `h in W`, then `h - lambda g in C` for some `lambda`, i.e. (both being normal
forms) `h = lambda g`, i.e. `x_i^(2^j) g = lambda g` on `H^m`, contradicting (b). So `h` is not in
`V_d`, hence not in any `V_e` with `e <= d` (`V_e subset V_{e+1}`), and the least `D` with
`h in V_{max(D, deg h)}` is at least `d + 1`.

**Corollary 2 (direct m = 3; proof plus symbolic identity).** For every depth `k >= 3`, every prime
`p > 5` with `2^k | p - 1`, every curve and every target with `5 t^3 + 3 a t + 6 b != 0 (mod p)`, the
direct system `{chains, NF(S_4(x1,x2,x3,t))}` has `LFD >= 7`. For `k = 2` the same holds whenever
`x1 g` is not proportional to `g` on `H^3` (checked per instance).

*Proof.* Exponents of `S_4` are `0..4` in each variable; for `k >= 2` the only exponent of chain
weight 2 is `3`, so the unique possible degree-6 mask is `(3,3,3)` and its coefficient is the
coefficient of `x1^3 x2^3 x3^3` in `S_4`, which is `-8 (5 t^3 + 3 a t + 6 b)` (sympy over
`Z[a,b,t]`, `symbolic.json`; re-checked per case in `rows.jsonl`, field `c333`). So `d = 6`.
(a) with `u = u_{1,0}`: the mask `(3,3,3)` goes to `(4,3,3)` (weight 5; `(0,3,3)`, weight 4, at
`k = 2`), every other mask has weight `<= 5` and gains at most 1. (b) for `k >= 3`: the coefficient of
`x3^4` in `S_4(x1,x2,x3,t)` is `S_3(x1,x2,t)^2` (`symbolic.json`). Fix `x1 in H`, `x1 != t`.
`S_3(x1, x2, t)` has degree 2 in `x2` with leading coefficient `(x1 - t)^2 != 0`, so it vanishes for at
most 2 values `x2`; for the other `x2`, `S_4(x1,x2,.,t)` has degree 4 and at most 4 roots. So the
slice `{x1} x H^2` holds at most `2*2^k + 4*2^k = 6*2^k < 4^k` zeros of `g` when `k >= 3`: `g` is
nonzero somewhere on every slice `x1 != t`, i.e. on at least `2^k - 1 >= 7` distinct slices, so it
is not supported in one level set `{x1 = c}`. Proposition 1 gives `LFD >= 7`.

Together with the kernel-checked `LFD <= mk + 1 = 3k + 1`, at `k = 2` this pins `LFD = 7` for every
instance satisfying the hypotheses (the exact runs confirm it); for `k >= 3` it is a lower bound
valid at every depth.

**Remark 3 (field-equation split; proof).** If `dim R/I_free > dim R/I_field` (some auxiliary root `z`
outside `F_p`, or a multiple root), then `I_free` is strictly inside `I_field`, so `z^p - z` is not in `I_free`.
Where the monomials of degree `< D0` span `R/I_free` (every certified case), `z^p - z = rho + f` with
`f in I_free` and `deg rho < D0 < p`; `rho` lies in `I_field \ I_free`, and it is missing from every
`V_e(F_field)`, `e < p`, because those coincide with the closures of `F_free`, which lie in
`I_free`. Hence `LFD(F_field) >= p`: with the auxiliary written over the full field through its
field equation the LFD is governed by `p`, not by `k`. When the two dimensions agree the ideals
coincide and so do the closures and LFDs.

## Results (exact, certified)

Grid (`cases.jsonl`, `skips.json`): primes 17, 97, 257, 641, 7681; per prime the curves
`(0,1)` and `(0,b)` (j = 0), `(1,0)` and `(a,0)` (j = 1728), two random `(a,b)`; targets: the
order-two point (first root of `x^3+ax+b`, 27 of 30 curves have one), two random rational points
(fixed across k), and a planted `P1 + P2 + P3` with `x_i in H` (one per curve and k; 10 skipped at
k = 1 where no point has `x in H`). Every planted target has a rational decomposition inside `Z`
(113 of 113). Note `H = F_17^*` at p = 17, k = 4.

Exact LFD histograms (key: LFD value, value: number of cases):

| k | direct m=3 | direct, `deg NF(S_4) = 6` | split m=3, `I_free` | split, `I_field` | m=2 reference | direct - m2 (paired) | split - m2 (paired) |
|---|---|---|---|---|---|---|---|
| 1 | 4: 100, 3: 5, 0: 6 | (deg <= 3 at k = 1) | 5: 83, 4: 28 | exact 4/5: 100; `>= p`: 11 | 3: 111 | +1: 100, 0: 5, -3: 6 | +2: 83, +1: 28 |
| 2 | 7: 105, 6: 11 | 7: 105 | 5: 116 | exact 5: 87; `>= p`: 29 | 4: 116 | +3: 105, +2: 11 | +1: 116 |
| 3 | 7: 117 | 7: 107 | 6: 111, 5: 6 | exact 6: 80; `>= p`: 37 | 4: 117 | +3: 117 | +2: 111, +1: 6 |
| 4 | 7: 117 | 7: 106 | 6: 117 | exact 6: 69; `>= p`: 48 | 4: 117 | +3: 117 | +2: 117 |

Kernel-checked upper bounds for comparison: `LFD <= 3k+1` (direct; 4, 7, 10, 13) and `2k+1` (m = 2).
Unit ideals among the direct systems: 74, 63, 51, 27 at k = 1..4 (the rest are proper with
local-to-global certificates); m = 2: 103, 94, 83, 66. `|Z|` (direct) median / max: 0/6, 0/28,
6/122, 18/488 at k = 1..4.

Degenerate k = 1 rows: six direct LFD-0 cases are `(a,0)` curves with the target `(0,0)`: for
`a = 1` `NF(S_4)` is the nonzero constant 64 (a constant input, LFD 0); for `p = 641, a = 446` the
degree-2 generator never produces a strict fall (6 solutions). Five k = 1 rows have
`deg NF = 2` and LFD 3.

Dimensions of `V_D` (quotient by the chains; the ordinary dimension is
`C(r+D, D) - |M_<=D| + q`, stored per row as `dim_V_D`). Min / median / max over all cases:

| encoding | k | D | `|M_<=D|` | `dim V_D / chains` min / median / max |
|---|---|---|---|---|
| direct | 2 | 6 / 7 | 64 / 64 | 1 / 1 / 64 ; 38 / 64 / 64 |
| direct | 3 | 6 / 7 | 466 / 502 | 1 / 1 / 10 ; 380 / 496 / 502 |
| direct | 4 | 6 / 7 | 2510 / 3302 | 1 / 1 / 13 ; 2814 / 3284 / 3302 |
| split | 2 | 4 / 5 | 129 / 192 | 25 / 25 / 47 ; 161 / 192 / 192 |
| split | 3 | 5 / 6 | 825 / 1291 | 246 / 248 / 797 ; 1162 / 1288 / 1291 |
| split | 4 | 5 / 6 | 2772 / 5282 | 346 / 348 / 394 ; 4750 / 5262 / 5282 |

Full table in `table.md`. Representative profile (p = 97, curve (0,1), planted target, k = 4,
target (92,78)): direct, quotient ranks `q_D` = 0,0,0,0,0,0,1,3278 for D = 0..7 (ordinary
`dim V_D` = 0,0,12,156,1026,4602,16055,50364), single strict fall at 7, `|Z| = 24`, certified at
D = 7 -> LFD 7; split `I_free`: q = 0,0,0,1,27,348,5246 (ordinary 0,0,12,169,1221,6144,27096),
falls at 4,5,6, `dim A = 36`, certified at 6 -> LFD 6 (`I_field` has `dim A = 20 < 36`, so
`LFD(F_field) >= p`); m = 2 at the same target: LFD 4.

### Reading (labels: exact = certified computation; empirical = pattern; proof = written above)

* exact: every one of the 461 direct m = 3 systems has an exact LFD; for k = 2, 3, 4 the value is 7
  on every case except 11 k = 2 cases with `deg NF(S_4) = 5` (LFD 6). At k = 2 the value 7 equals the
  kernel upper bound `3k+1 = 7`: the bound is attained.
* proof (Corollary 2): `LFD >= 7` for the direct system at every `k >= 3` whenever
  `5t^3 + 3at + 6b != 0`; the hypotheses hold in 107 + 106 cases at k = 3, 4 and the instance check
  in 105 cases at k = 2. The exceptional locus `5t^3 + 3at + 6b = 0` (it contains `(0,0)` on j = 1728)
  gives `deg NF = 5`; there the exact value is 6 at k = 2 and 7 at k = 3, 4 (10 and 11 cases).
* empirical: direct LFD = `deg NF(S_4) + 1` for every case with `deg NF = 6` (k = 2, 3, 4) and every
  `deg NF = 3` case at k = 1; the upper side (`LFD <= 7`) is not proved for general k. The
  degree-7 certificate used here cannot work at every depth: it needs `|M_<7| >= |Z|`, where
  `|M_<7|` grows polynomially (about `(3k)^6/720`) while `|Z|` is about `4 * 8^k / p` (observed medians
  above), so for primes close to `2^k` the certificate, or the value 7, must change at some depth.
  This is not a lower bound beyond 7; no claim is made about where it changes.
* exact: split encoding with the auxiliary free: LFD 5 at k = 1 (mostly), 5 at k = 2, 6 at k = 3
  and 4 (k = 3 exceptions with LFD 5: six targets with `y_R = 0`). Compared with the direct
  encoding at the same k it is 1 above at k = 1, 2 below at k = 2 and 1 below at k = 3, 4, at the
  price of one more variable (`3k+1`) and a larger certifying space (`|M_<=6| = 5282` mixed
  monomials vs `|M_<=7| = 3302` for the direct system at k = 4).
* proof (Remark 3) + exact: with the field equation `z^p - z` (auxiliary restricted to `F_p`), the
  LFD equals the free one when no auxiliary root leaves `F_p` (100, 87, 80, 69 cases at k = 1..4) and
  is `>= p` otherwise (11, 29, 37, 48 cases): the restriction to the full field moves the LFD to the
  scale of `p`.
* comparison with m = 2 at the same (p, curve, target, k): m = 2 is exactly 3 at k = 1 and 4 at
  k = 2, 3, 4 on all cases (consistent with the existing sweep `experiments/ic-curve-sweep`, where every
  subgroup-base system has LFD 4 at k = 2..5). Paired gap direct - m2 = +3 at k = 2 (105 cases), 3, 4 (all cases); split - m2 = +1 at
  k = 2 and +2 at k = 3, 4. All three encodings are constant in k over k = 2..4 (direct 7, m2 4)
  or step once (split 5 -> 6 at k = 3).

### Independent validation

`scripts/ic_bounds_m3_validate.py` recomputes the LFD in the ORDINARY ring (chains as explicit
inputs, full monomial space, no chain-quotient shortcut) with the independent validator engine
`experiments/ic-curve-validation/code/lfd_engine.py` on every k = 1 case and every 4th k = 2 case:
direct 140 / 140 agree; split 122 / 122 agree (split checked where `I_free = I_field`, so the zero
set is rational and the ideal radical, which that engine's certificate assumes). `validation.json`.

## Not done / limits

* Depth k >= 5 is not run (see heavy runs below); nothing here bounds the m = 3 LFD for k >= 5
  except the proven `LFD >= 7` (Corollary 2) and the kernel `LFD <= 3k+1`.
* The cases are a seeded structured grid, not a population sample; no claim about relation yield,
  relation rank or ECDLP cost is made.
* Propositions 1-3 are written proofs, not Lean-checked.

## Heavy runs not executed

* k = 5 direct (r = 15 chain variables, `|M_<=7| = 16384`): batched engine, one instance about
  2-4 h and 5-6 GB (`B` is 16384 x 16384 float64 = 2.1 GB plus products), e.g.
  `python3 -c "from csd import ic_dreg as I, ic_bounds_m3 as M; print(M.exact_profile('direct', I.PrimeCurve(97,0,1), 92, 5, engine='batched', max_columns=17000)['results'])"`
  (target x 92 on (0,1) mod 97; uses p = 97, 2^5 | 96).
* k = 5 split (`|M_<=6| = 17548` mixed monomials): similar, 6 GB, 3-5 h per instance.

## Reproduce

```
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 scripts/ic_bounds_m3_selfcheck.py        # S_4 checks, engine cross-check -> selfcheck.json
python3 scripts/ic_bounds_m3_symbolic.py         # symbolic identities -> symbolic.json
python3 scripts/ic_bounds_m3_run.py generate     # seeded grid (seed 20261107) -> cases.jsonl, skips.json
python3 scripts/ic_bounds_m3_run.py run --k 1 2 --budget-min 20
python3 scripts/ic_bounds_m3_run.py run --k 3 --budget-min 22            # rerun until "done"
python3 scripts/ic_bounds_m3_run.py run --k 4 --encodings direct m2 --results results-k4.jsonl --budget-min 22
python3 scripts/ic_bounds_m3_run.py select-k4-split      # k4-split-selection.json (35 cases)
python3 scripts/ic_bounds_m3_run.py run --k 4 --encodings split --select k4-split-selection.json --results results-k4-split.jsonl --budget-min 21
python3 scripts/ic_bounds_m3_run.py select-k4-split2     # k4-split-selection2.json (15 cases)
python3 scripts/ic_bounds_m3_run.py run --k 4 --encodings split --select k4-split-selection2.json --results results-k4-split.jsonl --budget-min 21
python3 scripts/ic_bounds_m3_run.py run --k 4 --encodings split --primes 17 97 257 --results results-k4-split-a.jsonl --budget-min 22
python3 scripts/ic_bounds_m3_run.py run --k 4 --encodings split --primes 641 7681 --results results-k4-split-b.jsonl --budget-min 22
# (every run command is resumable: rerun it until it prints "done")
python3 scripts/ic_bounds_m3_summarize.py        # summary.json, rows.jsonl, table.md, lfd_vs_k.png/svg
python3 scripts/ic_bounds_m3_validate.py         # independent ordinary-ring re-check -> validation.json
```

Run conditions: 4-core host shared with other jobs, `ulimit -v 3000000`, single-threaded BLAS;
medians 0.6 s (direct) / 7 s (split) per case at k = 3, 36 s / 69 s at k = 4; peak RSS about
0.9 GB. Every run command was split into chunks of at most 22 minutes.

## Files (sha256, first 16 hex digits)

| file | sha256 |
|---|---|
| cases.jsonl | 342cc28785eb2772 |
| results.jsonl (k = 1..3, all encodings) | 14c5eb015f2758fb |
| results-k4.jsonl (k = 4 direct, m2) | ffa0e6e27f1c3b7d |
| results-k4-split.jsonl / -a / -b (k = 4 split) | 36d839b4a1d21e93 / 40aaeb638383fba1 / 976f02bd65d6e4a4 |
| rows.jsonl (merged per case) | 2fed8fa19ed47b5d |
| selfcheck.json / symbolic.json / validation.json | 2a9b78abc01eda2b / e74507177ed480c2 / bcfd3c79fd0cf9bb |
| csd/ic_bounds_m3.py / csd/ic_bounds_m3_batched.py | see `sha256sum` (code may be re-hashed at review) |

Each result line records `generator_hash` (sha256 of the normalized generator(s), p, k), the
profile per degree (`monomials_le_D`, `quotient_rank`, `dim_V_D` ordinary, evaluation ranks,
`closure_check`, `engine`) and the certificate; `cases.jsonl` records p, a, b, curve class, target,
target class, planted points, k and the subgroup generator `g0`.
