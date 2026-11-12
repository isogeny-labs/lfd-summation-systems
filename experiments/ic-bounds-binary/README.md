# Item 6: binary ordinary and Koblitz systems at larger k, three summands, and what a cost statement needs

Status (2026-11-08): 1871 seeded cases, 1726 distinct systems, **every case has an exact certified
LFD** (1355 proper ideals certified by the local-to-global rank criterion, 516 unit ideals by
`1 in V_D`; no capped case). Seed `20261108`. Total closure time 5.3 CPU-hours, peak 1.28 GB per
process. No ECDLP running-time or speed claim is made anywhere in this directory.

| family | what | k | cases | exact | LFD histogram | max LFD | written `min(k,l)+2` | kernel `k+l+1` / `r+1` |
|---|---|---|---|---|---|---|---|---|
| m2 | two summands, `S_3`, `r = 2k` | 8 | 192 | 192 | 2:3, 3:70, 4:119 | **4** | 10 | 17 |
| m2 | | 9 | 192 | 192 | 2:16, 3:30, 4:146 | **4** | 11 | 19 |
| m2 | | 10 | 184 | 184 | 2:17, 3:4, 4:133, 5:30 | **5** | 12 | 21 |
| m3-direct | three summands, `S_4` descended, `r = 3k` | 2 | 201 | 201 | 3:43, 4:12, 5:42, 6:104 | **6** | (m=2 only) | 7 |
| m3-direct | | 3 | 231 | 231 | 4:54, 5:6, 6:30, 7:141 | **7** | | 10 |
| m3-direct | | 4 | 255 | 255 | 4:23, 5:41, 6:11, 7:180 | **7** | | 13 |
| m3-direct | | 5 | 8 | 8 | 5:2, 7:6 | **7** | | 16 |
| m3-split | three summands, two `S_3`, auxiliary `z` in `F_(2^n)`, `r = 3k+n` | 2 | 201 | 201 | 3:80, 4:121 | **4** | | 3k+n+1 = 16 (max) |
| m3-split | | 3 | 231 | 231 | 3:14, 4:213, 5:4 | **5** | | 19 (max) |
| m3-split | | 4 | 176 | 176 | 4:166, 5:10 | **5** | | 21 (max) |

Full tables (by `n`, curve class, base kind, target class, with `|Z|`, Macaulay sizes, closure
seconds and memory): [`tables.md`](tables.md); machine-readable: [`summary.json`](summary.json);
figures: [`lfd_vs_k.png`](lfd_vs_k.png) (`.svg`, `.pdf`), [`cost_proxy.png`](cost_proxy.png).

## Systems

* Field `F_(2^n)`, modulus `csd.ic_dreg.binary_modulus(n)` (stored per row as `modulus`).
  m2: `n = 11, 12, 13, 14`, `k = 8, 9, 10` (all `k < n`). m3: `n = 5, ..., 9`.
* Curves `y^2 + xy = x^3 + a x^2 + b`: Koblitz `(0,1)`, `(1,1)` and two seeded random ordinary
  curves per `n` (`b` not in `{0,1}`; `csd.ic_bounds_binary.curves_for`, values in every row and
  listed in `tables.md`), e.g. `n = 11`: `(304,1583)`, `(1509,529)`.
* Coordinate spaces `V` of dimension `k` (canonical reduced basis stored per row):
  `frobenius-invariant` (`csd.ic_curve_sweep.invariant_candidates`, i.e. images `g(sigma)` of
  `frobenius_invariant_subspaces`, at most two per `(n,k)`) and `random-subspace` (seeded, verified
  not Frobenius-invariant). Invariant subspaces of dimension 8, 9, 10 exist only where `x^n - 1`
  has a divisor of degree `n - k`: `n = 11` only `k = 10` (the trace-zero hyperplane
  `{z^2 + z}`, see P4), `n = 12` and `n = 14` all three, `n = 13` none.
* Targets: order two `(0, sqrt b)`; two random rational points; one planted target `P + Q` with
  `x(P), x(Q) in V`. **Caveat for m = 3:** the planted target is a planted *two*-term sum, so it does
  not guarantee a three-term decomposition (97 of 149 m3-direct and 83 of 127 m3-split planted
  systems have `|Z| > 0`).
* Families (`csd/ic_bounds_binary.py`):
  - `m2`: Weil descent of `S_3(x1, x2, x_R)`, `x1, x2 in V` (`csd.ic_curve_sweep.prepare_binary`),
    `n` Boolean equations of degree 2 in `2k` variables;
  - `m3-direct`: Weil descent of `S_4(x1, x2, x3, x_R)` with three copies of `V`
    (`csd.ic_dreg.binary_descent`), `n` equations of Boolean degree 6 (`S_4` has degree 4 in each
    variable and `x^3` has Boolean degree 2) in `3k` variables;
  - `m3-split`: `S_3(x1, x2, z) = 0 = S_3(z, x3, x_R)` with an auxiliary coordinate `z` ranging over
    all of `F_(2^n)` (`csd.ic_research_campaign.binary_split`), `2n` equations of degree 3 in
    `3k + n` variables.
* Seeds: every random choice is drawn from `random.Random(sha256("20261108:..."))` keyed by
  `(n)`, `(n,k)` or `(family,n,k,a,b,basis)`; each row stores the full specification
  (`modulus, a, b, basis, target, target_class, planted`), its `case_id`, the `system_hash`
  (sha256 of the canonical equation list) and `points_sha256`.

## Certificates (what makes each value exact)

`csd.ic_local_certificate.local_profile`, closure modulo the field equations `x_i^2 = x_i`
(`csd.ic_boolean_closure`), every closure re-checked by `check_containing_space` (contains the
inputs, stable under every variable shift):

* unit ideal: `exact_lfd` = least `D` with `1 in V_D`;
* proper ideal: at the first `D*` with evaluation of `M_<D*` onto `Z` and
  `codim V_D* = |Z|` (both GF(2) ranks), the exact LFD is the last strict fall
  (`formal/ICFallDegrees/LocalReductionBound`; the inequality used is recorded per row in
  `local_checks`, and `runner` asserts `codim = |Z|`).

`Z` is enumerated independently of the Boolean equations: m2 by `csd.ic_dreg.binary_solutions`
(field arithmetic, solving the F_2-affine equation in `x2` for each `x1`) and asserted equal to the
descent's `boolean_points`; m3-direct by brute force over `V^3` with `s4_binary_value`; m3-split by
the scalar `S_3` chain over all `z`. `local_profile` re-checks that every point vanishes on every
equation.

### Independent checks

| file | what | result |
|---|---|---|
| `validation-hashes.json` | every row rebuilt from its stored specification (fresh `Z` cross-check included), `system_hash` compared | 1871 / 1871 match |
| `validation-m3.json` | engine A: `experiments/ic-curve-validation/code/lfd_engine.py` (ordinary ring `GF(2)[z]`, field equations as degree-2 inputs, no shared closure code), brute-force `Z` over all `2^r` assignments; 10 m3-direct (`r = 6, 9`) + 10 m3-split (`r = 11, 12`) rows | 20 / 20 equal LFD |
| `validation-m3.json`, `validation-m3direct-mxl.json` | engine B: `csd.ic_dreg._mxl_closure` (MutantXL closure of the full Boolean Macaulay matrix to a fixed point), rank of `V_D` compared at every `D <= D*` | 11 / 11 rows agree at every checked degree (m3-direct `k = 2..5` incl. a `k = 5`, LFD 7 row through `D = 7`) |
| `validation-m2.json` | engine B on two rows per `k` for m2, including two `k = 10`, LFD-5 rows (`n = 11, |Z| = 499`; `n = 12, |Z| = 552`) through `D = 5` | 6 / 6 agree at every degree |

Engine A is dense and cannot reach the m2 sizes (`C(r+D, D)` columns); engine B shares the GF(2)
echelon C kernel (`csd/_elim.c`) but not the closure logic.

## Results

### (a) Two summands, `k = 8, 9, 10`

* **exact.** Max LFD is 4 at `k = 8, 9` and 5 at `k = 10`. Every value is far below the written
  bound `min(k,l)+2 = k+2` (10, 11, 12) and the kernel bound `k+l+1 = 2k+1` (17, 19, 21). The
  distance remaining to the written bound is 6, 7, 7.
* **exact.** These are the first certified binary two-summand systems with LFD 5 in the
  repository (the k <= 7 sweep `experiments/ic-curve-sweep` has max 4). 30 systems, all at `k = 10`:
  `n = 11` (all 20 proper `k = 10` systems) and `n = 12` (10 systems).
* **empirical.** At `k = 10` the value is separated by the number of solutions with no overlap in
  this grid: LFD 5 for `|Z|` in `[472, 1056]` (30 systems), LFD <= 4 for `|Z|` in `[46, 352]`
  (137 proper systems). At `k = 8, 9` there is no such separation (LFD 3 and 4 overlap in `|Z|`).
* **empirical.** `|Z|` tracks `2^(2k-n)`: mean `|Z|` over random bases and non-order-two targets is
  31.8, 124.2, 508.2 for `n = 11` (`2^(2k-n)` = 32, 128, 512), 16.3, 68.8, 248.2 for `n = 12`,
  7.1, 33.2, 124.6 for `n = 13`, 4.8, 16.3, 64.8 for `n = 14`. The trace-zero hyperplane at
  `n = 11` doubles it (`|Z|` 991 to 1056) or empties it (order-two target: unit ideal, LFD 2).
  Rational unordered decompositions are about 0.45 of `|Z|` at every `k`.
* **exact.** Order-two targets: LFD histogram 3:46, 4:2 (`k = 8`); 3:30, 4:18 (`k = 9`);
  3:4, 4:36, 5:6 (`k = 10`). Unit ideals: LFD 2 at `k = 9, 10`, but 2, 3 or 4 at `k = 8`.
* **exact.** No curve-class or base-kind effect on the maximum: Koblitz and random ordinary curves,
  random and Frobenius-invariant bases all reach 4 at `k = 8, 9` and 5 at `k = 10`; the
  Frobenius-invariant bases enter the LFD-5 regime through larger `|Z|`, not otherwise
  (descriptive only; the grid is not a controlled comparison).

### (b) Three summands

* **exact.** Direct `S_4` systems: max LFD 6 at `k = 2` (`r = 6`, kernel `r+1 = 7`), 7 at
  `k = 3, 4, 5` (`r = 9, 12, 15`). Each exceeds the two-summand formula `k+2` (which is proved only
  for `m = 2`; no contradiction). LFD 7 occurs for unit and proper ideals alike.
* **empirical.** The direct maximum stays at 7 from `k = 3` to `k = 5`; the `k = 5` evidence is
  8 systems on one curve (`(1,1)`) and one random basis, so it is a pattern, not a plateau claim.
* **exact.** Split systems (auxiliary full-field coordinate) have max LFD 4 at `k = 2` and 5 at
  `k = 3, 4` while having `n` more variables. On the 161 shared order-two targets (same `n, k`,
  curve, basis) LFD(direct) > LFD(split) in 54 pairs, equal in 103, smaller in 4.
* Coverage gap: m3-split `k = 4` has `n = 5, 6, 7` complete but only 1 case at `n = 8` and none at
  `n = 9` (199 s and 1.28 GB per `n = 8` case); m3-direct `k = 5` only `n = 7, 9`, curve `(1,1)`, one
  random basis (8 cases, 50 to 350 s and up to 1.25 GB each). See "Runs not done".

### (c) Cost-related quantities per system

Each row of `results-*.jsonl` stores `solutions = |Z|`, `input_degree`, `exact_lfd`,
`certificate_degree D*`, `macaulay_at_certificate {degree, rows, columns}`,
`cost_degree = max(LFD, input degree)` with `macaulay_at_cost_degree`,
`dim_closure_at_certificate` and `codim_closure_at_certificate (= |Z|)`. Columns are
squarefree monomials `|M_<=D| = sum_{i<=D} C(r,i)`; rows are the standard Boolean Macaulay
rows `sum_e |M_<=D-deg e|`. `related-*.jsonl` adds the static (non-iterated) Macaulay degrees
`d_lin` (all affine-linear forms vanishing on `Z` in the row space) and `d_full` (the row space
alone meets the certificate), from `scripts/ic_bounds_binary_related.py` (m3-direct every row,
column cap 32000, `k = 5` cap 17000; m2 and m3-split the first row per
(family, n, k, curve, base, target) group, cap 15000).

Per-`(n,k)` table in `tables.md` ("Cost-related quantities"). Headline counts:

* `D* - LFD` is 0 in 1817 of 1871 rows and 1 in the other 54; `D* <= 5` for every m2 and m3-split row,
  `D* <= 7` for every m3-direct row. Largest certificate matrices: m2 `k = 10`,
  `D* = 5`: 21,700 columns, up to 16,212 standard rows (`2^14.4` columns against `2^20` for the full
  Boolean space); m3-direct `k = 5`, `D* = 7`: 16,384 columns (`2^14` against `2^15`).
* HKYY degree `max(LFD, |Z|)` (see below): it is `>= r` in 1041 of the 1355 proper systems (m2 403
  of 525, every proper m2 system at `n = 11`; m3-direct 406 of 449; m3-split 232 of 381), i.e. the HKYY route is charged the whole space `M_<=r`
  (`2^r` columns, the size of exhaustive coordinate enumeration) in most of the grid.
* Static Macaulay: m3-direct `d_full` is 5 to 9 at `k = 2`, 7 to 10 at `k = 3`, 9 to 12 (`= r`) at
  `k = 4` against `D* <= 7`; for m2 the static scan never met the certificate below the cap
  (`d_full > 6` at `k = 8`, `> 5` at `k = 9`, `> 4` at `k = 10`) in 111 of the 117 sampled rows; the
  other 6 are unit ideals (five with `1` already in the degree-2 Macaulay row space, one `k = 8`
  row with `d_full = 5` against `D* = 4`). `d_lin` is not a solving
  degree here: with `|Z| >= r + 1` points the target linear space is often `0` and `d_lin = 2`.

## Proofs used above (written in full; none is Lean-checked)

**P1 (counting lower bound for the certificate degree).** If evaluation `span M_<D -> F_2^Z` is
onto then `|M_<D| >= |Z|`. Hence `D* >= min{D : sum_{i<D} C(r,i) >= |Z|}`. (Recorded per row as
`counting_certificate_lower`; it is never larger than `D*`, asserted in the summary; it equals
`D*` in 144 proper rows: 89 m2, 43 m3-direct `k = 2`, 12 m3-split.)

**P2 (Z is read off `V_D*` by linear algebra of size `|Z|`).** Let `R = F_2[x]/(x_i^2 + x_i)` and
`I` the ideal of the inputs. `R` is the ring of all functions `F_2^r -> F_2`, so `I` is the set of
functions vanishing on `Z` and evaluation gives `R/I = F_2^Z`. Assume at `D` that the image of `V_D`
in `R` equals `I cap M_<=D` and that `M_<D` evaluates onto `F_2^Z` (the certificate). Choose
`B subset M_<D` whose evaluations form a basis of `F_2^Z`. For `b in B`, `x_i b` (normal form) has
degree `<= D`, and there are unique `c_{b'}` with `x_i b - sum c_{b'} b'` vanishing on `Z`, i.e.
lying in `I cap M_<=D = V_D`; since `span B cap V_D = 0` and `dim M_<=D / V_D = |Z| = |B|`, the
`c_{b'}` are obtained by reducing `x_i b` and `B` modulo an echelon basis of `V_D`. This gives the
matrix `T_i` of multiplication by `x_i` on `R/I` in the basis `B`. Under `R/I = F_2^Z`, `T_i` is
diagonal with entries `z_i` (`z in Z`), so the joint eigenvalue vectors of the commuting
`T_1, ..., T_r` are exactly the points of `Z` (e.g. split by `ker T_1`, `ker (T_1 + 1)`
recursively, `O(r |Z|^3)` field operations). So once `V_D*` is known, `Z` costs an echelon form
of `|M_<=D*|` columns plus `r` matrices of size `|Z|`.

**P3 (summation polynomial depends on `b` only; order-two target).** For
`y^2 + xy = x^3 + a x^2 + b`, `S_3(x1,x2,x3) = (x1x2 + x1x3 + x2x3)^2 + x1x2x3 + b`
(`csd.ic_research_campaign.s3_boolean`, `binary_descent` for m = 2), and `binary_descent` builds
`S_4` from the same `A, B, C` blocks with constants `b` only. So the Koblitz curves `(0,1)` and
`(1,1)` give identical equations for the same `x_R` (this is why 1871 cases are 1726 systems).
At `x_R = 0`: `S_3(x1,x2,0) = (x1x2)^2 + b = (x1x2 + sqrt b)^2`, so
`Z = {(x1,x2) in V^2 : x1 x2 = sqrt b}`.

**P4 (the only 10-dimensional Frobenius-invariant subspace of `F_(2^11)` is `{z^2 + z}`).**
`F_(2^n)` is a cyclic `F_2[sigma]`-module `= F_2[x]/(x^n - 1)`; for `n = 11`,
`x^11 - 1 = (x + 1) g(x)` with `g` irreducible of degree 10 (2 has order 10 mod 11). A
10-dimensional invariant subspace is the image of a degree-1 divisor, i.e. of `sigma + 1`, which is
`{z^2 + z}` (kernel `F_2`, dimension 10). There is none of dimension 8 or 9; at `n = 13`
(2 has order 12 mod 13) there is none of dimension 8, 9, 10.

**P5 (`D* <= d_full`).** Every row `mu f` of the Boolean Macaulay matrix at `D` is obtained from
the input `f` by `deg mu` single-variable multiplications, each of an element of degree `< D`
(normal form degree grows by at most one), so the row space lies in `V_D subset I cap M_<=D`. If
the row space already has codimension `|Z|` and `M_<D` is onto `F_2^Z`, the same holds for `V_D`,
so the closure certificate holds at `D`. (Asserted for every row in `related-*.jsonl`.)

## What a cost statement needs (no ECDLP claim is made)

Solving cost depends on the LFD together with the number of solutions (Huang-Kosters-Yang-Yeo,
Prop. 2.8, as summarized in `docs/reviews-curve-lfd-novelty-20261105.md`: solving needs
`V_{max(d_F, e)}` with `e` the number of solutions). For these Boolean systems:

1. **The HKYY degree is useless in most of the grid.** `max(LFD, |Z|) >= r` in 1041 of 1355
   proper systems, and `M_<=r` is the whole `2^r`-dimensional Boolean space, the size of
   enumerating the coordinates. Since `|Z|` is near `2^(2k-n)` (empirical), any `n <= 2k - log2(2k)`
   puts the two-summand systems in this regime.
2. **The measured certificate degree is small and bounded only by measurement.** `D*` is LFD or
   LFD + 1 and at most 5 (m2, m3-split) or 7 (m3-direct) here; P2 shows `Z` is then a linear
   algebra problem of size `|M_<=D*|` and `|Z|`. A cost statement needs a proved bound on `D*`
   (or on LFD plus a degree at which `M_<D` spans `R/I`) as a function of `k, n`; the available
   proved bounds are the kernel `2k + 1` (`r + 1`) and the written `k + 2` for m = 2, against
   measured 4 to 5. For m = 3 there is no written bound and the measured value (7 at `r = 9..15`)
   is already above the m = 2 formula.
3. **The closure, not one Macaulay matrix, is what reaches `D*`.** The static matrix needs
   `d_full` up to `r` (m3-direct `k = 4`) or more than the cap (m2), so a cost statement must
   charge the iterated closure (rounds, rows added), which this item records only as wall seconds
   (`closure_seconds`) and memory (`peak_rss_mb`), not as an operation count.
4. **A stopping rule that does not know `|Z|`.** The certificate tests `codim V_D = |Z|` with `Z`
   enumerated independently. A solver must detect `D*` from `V_D` alone (for instance a
   border-basis commutation test on `T_i`); no such test is established or run here.
5. **Relation yield and count of systems.** About 0.45 `|Z|` rational unordered decompositions per
   system (empirical), the factor base size `2^k`, and the number of targets to process are outside
   this item; together with 2 to 4 they are what an index-calculus cost would need.

## Runs not done (resource limits)

* m2 `n = 11, k = 10`: done for the trace-zero basis and the first random basis on curves
  `(0,1)`, `(1,1)`, `(304,1583)` (24 cases); 24 cases left (second random basis, curve
  `(1509,529)`). Each LFD-5 system took 157 to 702 s and 1.2 GB.
* m3-split `k = 4`, `n = 8, 9`: 79 cases left (`n = 9`, `D = 5`: 27,896 columns, may exceed the
  0.9 GB byte cap).
* m3-direct `k = 5` full grid (`n = 6..9`, four curves, invariant and random bases): 152 cases left.
* m2 `k >= 11`: `D = 5` already has 35,443 columns at `k = 11`; the `n <= 2k - 9` (LFD-5) regime
  would need about 3 GB and 30 min per system and `D = 6` would need tens of GB.

Exact commands, memory and time estimates are in the structured result of this item.

## Reproduce

All commands with `ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`.
Runs are resumable (rows already in `--out` are skipped) and were executed in chunks of at most
24 minutes (`timeout 1440`, `--submit-seconds 900..1100` so started cases finish).

```sh
# (a) m = 2
python3 scripts/ic_bounds_binary_run.py --family m2 --n 14 13 12 --k 8 9 10 --workers 2 \
    --submit-seconds 960 --out experiments/ic-bounds-binary/results-m2.jsonl
python3 scripts/ic_bounds_binary_run.py --family m2 --n 11 --k 8 9 --workers 2 \
    --out experiments/ic-bounds-binary/results-m2.jsonl
for c in 0 1 2; do python3 scripts/ic_bounds_binary_run.py --family m2 --n 11 --k 10 --random-bases 1 \
    --curve-indices $c --workers 2 --submit-seconds 960 --out experiments/ic-bounds-binary/results-m2.jsonl; done
# (b) m = 3
python3 scripts/ic_bounds_binary_run.py --family m3-direct --n 5 6 7 8 9 --k 2 3 4 --max-degree 14 \
    --out experiments/ic-bounds-binary/results-m3-direct.jsonl
python3 scripts/ic_bounds_binary_run.py --family m3-direct --n 7 9 --k 5 --random-bases 1 --invariant-bases 0 \
    --curve-indices 1 --max-degree 14 --submit-seconds 960 --out experiments/ic-bounds-binary/results-m3-direct.jsonl
python3 scripts/ic_bounds_binary_run.py --family m3-split --n 5 6 7 --k 2 3 4 --max-degree 14 \
    --out experiments/ic-bounds-binary/results-m3-split.jsonl
python3 scripts/ic_bounds_binary_run.py --family m3-split --n 8 9 --k 2 3 --max-degree 14 \
    --out experiments/ic-bounds-binary/results-m3-split.jsonl   # plus the single n=8, k=4 case
# (c) static Macaulay degrees
python3 scripts/ic_bounds_binary_related.py --results experiments/ic-bounds-binary/results-m3-direct.jsonl \
    --out experiments/ic-bounds-binary/related-m3-direct.jsonl            # k=5 rows: --max-columns 17000
python3 scripts/ic_bounds_binary_related.py --results experiments/ic-bounds-binary/results-m2.jsonl \
    --per-group 1 --max-columns 15000 --out experiments/ic-bounds-binary/related-m2.jsonl
python3 scripts/ic_bounds_binary_related.py --results experiments/ic-bounds-binary/results-m3-split.jsonl \
    --per-group 1 --max-columns 15000 --out experiments/ic-bounds-binary/related-m3-split.jsonl
# validation
python3 scripts/ic_bounds_binary_validate.py --results experiments/ic-bounds-binary/results-m3-direct.jsonl \
    experiments/ic-bounds-binary/results-m3-split.jsonl --per-family 10 --max-vars 12 --mxl-per-k 1 \
    --mxl-max-columns 8000 --out experiments/ic-bounds-binary/validation-m3.json
python3 scripts/ic_bounds_binary_validate.py --results experiments/ic-bounds-binary/results-m3-direct.jsonl \
    --per-family 0 --mxl-per-k 1 --mxl-max-columns 17000 --out experiments/ic-bounds-binary/validation-m3direct-mxl.json
python3 scripts/ic_bounds_binary_validate.py --results experiments/ic-bounds-binary/results-m2.jsonl \
    --per-family 0 --mxl-per-k 2 --mxl-max-columns 22000 --out experiments/ic-bounds-binary/validation-m2.json
python3 scripts/ic_bounds_binary_validate.py --per-family 0 --hash-all --out experiments/ic-bounds-binary/validation-hashes.json
# tables and figures
python3 scripts/ic_bounds_binary_summarize.py
```

The m3-split `n = 8, 9`, `k = 4` command of the original grid was stopped after one case; the
m2 `n = 12` run was resumed once after a chunk timeout. Resuming never re-runs a stored `case_id`.

Input-file hashes (sha256):

```
3902fa9f4e2aae74493eb5cf44835329e28f262ba0231dd29d1675e40fb5ca7d  results-m2.jsonl
c28a7c4ebc0e6d275616fc8f7172ad29bc82cfe4f877b2f4ac5a7c85bfcb1e87  results-m3-direct.jsonl
8d76db26375729a00c7c58865a0771f48e5dc5a74b721e3c0a0385ed98586a2d  results-m3-split.jsonl
7c05334bf8f3cb615b4d506000f18aa4877930a94733db151be1a91c11cfa93b  related-m2.jsonl
5e801ef902a7793d0eabf0da1011699bf011a1d56efc2b050ca5e7b71957b048  related-m3-direct.jsonl
a19ff23cf2f9cbd40b58b284bc58834091d2a272937ee6bdac7557a026ac5190  related-m3-split.jsonl
```

ecdb: not queried. The curves are seeded small-field curves generated here, and no curve-level
fact (group order, isogeny class) is used; the computation is of the descended equations only.
