# Three-summand direct chain system: last fall degree against depth

> **Note (2026-11-12).** Some runs listed below as not run were run later and are reported in `paper/curve-lfd-bounds.tex`, including the `k = 6` degree-7 closure reported in the paper (Section 7, "Three summands"). The lists below record the state when this file was written.

Status (2026-11-10): exact LFD at depth k = 5 on 12 instances (8 primes with 32 | p - 1), closure
profiles for D = 5..8 at k = 2..5 on the same instances, the free-extension recursion test, and
the two-summand reference on the same instances. The depth-6, degree-7 closure has run on one instance (see "Depth 6"): V_7 falls behind
the ideal at k = 6, so LFD >= 8 there. Labels: **exact** (certified
computation, command and certificate named), **proof** (written out here or in
`experiments/ic-bounds-m3/README.md`), **empirical**, **conjecture**.

## System

Curve `E: y^2 = x^3 + a x + b` over `F_p`, `H` the subgroup of order `2^k` of `F_p^*`. Inputs:
the `3k` chain equations `u_{i,j}^2 = u_{i,j+1}`, `u_{i,k-1}^2 = 1` (`u_{i,0} = x_i`) and
`g = NF(S_4(x1, x2, x3, t))`, with `S_4` and the normal form exactly as in
`experiments/ic-bounds-m3` (`csd.ic_bounds_m3.direct_generator`). Every target here satisfies
`5t^3 + 3at + 6b != 0`, so `deg g = 6` (checked per row: `generator_degree = 6`).
`V_D` is the iterated degree-D closure, represented modulo the chains (equal to the ordinary
closure for `D >= 2`), `M_<=j` the squarefree chain masks of degree `<= j`,
`p_j = dim(V_D cap M_<=j)`, `Z` the zero set of `g` on `H^3` (`R/I = F_p^Z`; the chain ideal is
radical), `dim(I cap M_<=j) = |M_<=j| - rank(ev: M_<=j -> F_p^Z)`.

## Engine and certificates

`scripts/ic_bounds_m3_depth.py` (`CompactClosure`). Same object and same column order as
`csd.ic_bounds_m3_batched.batched_closure`, but the basis is kept in reduced row echelon form
with only its non-pivot columns stored (uint16 residues, row blocks), so a reduction costs
`rank x #non-pivot columns`; float64 products, every accumulated sum below `n (p-1)^2 < 2^53`.

* Every row added is the input or a variable shift of a span element of degree `< D`, so the span
  is inside `V_D`.
* A run closes either by a **full closedness sweep** (every variable shift of every final basis row
  of degree `< D` reduces to zero, the input is in the span: then the span equals `V_D`;
  `closure_check = closed-sweep`), by `1` becoming a pivot (`unit`), or by the **dimension
  sandwich** (rank reaches the exact `dim(I cap M_<=D)`, and span `subset V_D subset I cap M_<=D`;
  `closure_check = rank=dim(I cap M_<=D)`), as in `experiments/ic-bounds-m3`.
* Exact LFD (the rule of `experiments/ic-bounds-m3/README.md`, criterion proved there): `V_D = 0`
  for `D <= 5` (no input of degree `<= D`); unit: least `D` with `1 in V_D`; proper: at the first
  `D` with `rank ev(M_<D) = |Z|` and `dim V_D = |M_<=D| - |Z|`, LFD = last strict fall
  (`dim(V_D cap M_<=D-1) > dim V_{D-1}`) at or below `D`. Once this holds at `D`, it holds at every
  `D' >= D` (`V_D' = I cap M_<=D'`), so no degree above `D` falls behind.

Engine checks (exact):
* `selfcheck-k23.json`, `selfcheck-k4.json`: rank and `dim(V_D cap M_<=D-1)` agree with
  `csd.ic_bounds_m3.exact_profile` (`chain_degree_closure` at k <= 3, `batched_closure` at k = 4)
  for every D of the reference profile, with and without the sandwich: 112 (instance, k, D, mode)
  comparisons on 3 instances at k = 2, 3 and 2 instances at k = 4, all equal.
* Block storage re-check (chunk 100, batch 200, several blocks and compactions) reproduces 8 k = 3, 4
  profiles; a run without the sandwich at k = 4 closes by the full sweep with the same profile.
* `nosandwich-k5.jsonl`: k = 5, D = 7 without the sandwich (full closedness sweep) on two instances;
  compared in `table.md` (last section).

## Instances

Seeded (`make_instance`, seed string `m3depth-<p>-<seed>-<planted>`): a random nonsingular curve and
either a random rational target (`-r`) or a planted target `P1 + P2 + P3` with `x(P_i)` in the
order-4 subgroup (`-pl`), which decomposes over the factor base at every depth `k >= 2`. The same
`(p, a, b, t)` is used at every depth.

| instance | p | a | b | target (t, y) | 5t^3+3at+6b | \|Z\| at k = 2, 3, 4, 5 |
|---|---|---|---|---|---|---|
| p97-r | 97 | 44 | 26 | (51, 73) | 67 | 0, 0, 45, 469 |
| p193-r | 193 | 73 | 50 | (153, 115) | 26 | 6, 12, 69, 544 |
| p257-r | 257 | 217 | 48 | (86, 50) | 153 | 0, 6, 24, 195 |
| p7681-r | 7681 | 1876 | 2678 | (6223, 7322) | 3423 | 0, 0, 0, 12 |
| p97-pl | 97 | 45 | 4 | (94, 91) | 66 | 1, 11, 83, 383 |
| p193-pl | 193 | 123 | 170 | (21, 137) | 69 | 3, 15, 52, 295 |
| p257-pl | 257 | 8 | 50 | (16, 26) | 90 | 10, 22, 58, 268 |
| p7681-pl | 7681 | 4554 | 578 | (4298, 3246) | 3052 | 10, 22, 46, 103 |
| p353-r | 353 | 219 | 74 | (123, 73) | 26 | 0, 3, 15, 141 |
| p449-pl | 449 | 331 | 31 | (382, 304) | 442 | 10, 22, 79, 229 |
| p577-r | 577 | 21 | 468 | (444, 360) | 313 | 0, 0, 6, 39 |
| p641-pl | 641 | 408 | 325 | (29, 3) | 423 | 6, 6, 36, 129 |

## Results

### Task 1: exact LFD at k = 5 (exact)

| instance | \|Z\| | profile of V_7, p_0..p_7 | dim(I cap M_<=7) | certificate | LFD |
|---|---|---|---|---|---|
| p97-r | 469 | 0 0 0 125 1475 4475 9480 15915 | 15915 | local-to-global at 7 (rank ev(M_<7) = 469 = \|Z\|) | **7** |
| p193-r | 544 | 0 0 0 48 1397 4400 9405 15840 | 15840 | local-to-global at 7 | **7** |
| p257-r | 195 | 0 0 0 381 1746 4749 9754 16189 | 16189 | local-to-global at 7 | **7** |
| p7681-r | 12 | 0 6 109 564 1929 4932 9937 16372 | 16372 | local-to-global at 7 | **7** |
| p97-pl | 383 | 0 0 0 193 1558 4561 9566 16001 | 16001 | local-to-global at 7 | **7** |
| p193-pl | 295 | 0 0 0 281 1646 4649 9654 16089 | 16089 | local-to-global at 7 | **7** |
| p257-pl | 268 | 0 0 0 326 1676 4676 9681 16116 | 16116 | local-to-global at 7 | **7** |
| p7681-pl | 103 | 0 0 66 491 1841 4841 9846 16281 | 16281 | local-to-global at 7 | **7** |
| p353-r | 141 | 0 0 5 435 1800 4803 9808 16243 | 16243 | local-to-global at 7 | **7** |
| p449-pl | 229 | 0 0 0 365 1715 4715 9720 16155 | 16155 | local-to-global at 7 | **7** |
| p577-r | 39 | 0 0 82 537 1902 4905 9910 16345 | 16345 | local-to-global at 7 | **7** |
| p641-pl | 129 | 0 0 9 447 1812 4815 9820 16255 | 16255 | local-to-global at 7 | **7** |

`|M_<=7| = 16384` (15 chain variables). `V_6 = span(g)` (rank 1) on every instance, `V_7 cap M_<=6`
has dimension `>= 9405`, so the only strict fall is at 7; `V_7 = I cap M_<=7` and `M_<7` spans
`R/I`. Answer: **still 7 at k = 5, on all 12 instances (8 primes with 32 | p - 1, generic targets)**.
On p7681-r and p97-r the k = 5, D = 7 closure was also run without the sandwich and closed by the full
sweep with the identical profile (`nosandwich-k5.jsonl`).
With the proof `LFD >= 7` (Corollary 2 of `experiments/ic-bounds-m3`, k >= 3, generic targets) the
lower side is a proof at every depth; the value at k = 5 is exact only for these 12 systems.

All k = 2..4 values on the same instances are 7 as well (36 systems, `table.md`; unit ideal where
`|Z| = 0`, 9 of them). k = 5 runs: 3.5 to 7.5 min each, single thread, peak RSS 0.8 GB (1.3 GB for the first
run with batch 1024).

### Task 2: closure profiles by degree and depth (exact)

Full table (every instance, k = 2..5; D = 6, 7 everywhere, D = 8 at k <= 4 on the first 8 instances): `table.md`,
machine-readable `closures.jsonl`, `summary.json`.

* **D = 5**: `V_5 = 0` modulo the chains for every generic target at every depth (the only input of
  degree `<= 5` would be `g`, of degree 6), so `p_j = 0`, while `dim(I cap M_<=5)` is 53..63,
  360..382, 1503..1586, 4400..4932 at k = 2..5 (column `ideal_dim_le_j[5]` of the D = 7 rows).
  This is immediate from the definition (no input of degree `<= 5`, and modulo the chains the
  closure of the chain equations alone is 0); no D = 5 runs were needed.
* **D = 6**: `V_6 = span(g)`, profile `(0,0,0,0,0,0,1)` at every k = 2..5 on every instance (exact),
  while `dim(I cap M_<=6) = |M_<=6| - rank ev(M_<=6)` is 54..64, 444..466, 2427..2510,
  9405..9937 at k = 2..5.
* **D = 7**: `V_7 = I cap M_<=7` at every k = 2..5 on every instance (exact), so
  `p_j(k) = |M_<=j| - rank ev(M_<=j -> F_p^Z)`: the profile is fixed by the zero set, not by a
  free extension of the previous depth.
* **D = 8**: computed at k = 2..4 on the first 8 instances (= `I cap M_<=8`); at every other
  computed (instance, k), including k = 5, `V_8 = I cap M_<=8` follows from the certificate at 7.

**Three-block recursion `p_j(k+1) = sum_w C(3,w) p_{j-w}(k)`** (`table.md`, section
"Recursion"): 
* D = 6: holds at every step k = 2 -> 3 -> 4 -> 5, but vacuously (both sides are `span(g)`,
  `p = (0,...,0,1)`).
* D = 7, 8: holds only when both depths are unit ideals (then `V_D = M_<=D` and `|M_<=j|` obeys the
  recursion by Vandermonde): p97-r 2 -> 3, p7681-r 2 -> 3 -> 4, p577-r 2 -> 3. It fails on every step
  where either depth has `Z` nonempty (45 of 52 D = 7, 8 steps; the other 7 are the unit-ideal
  steps), already at `j <= 2`. Example p193-r, D = 7,
  k = 4 -> 5: predicted `0 0 14 272 1457 4396 9397 15832`, actual `0 0 0 48 1397 4400 9405 15840`.

So at k <= 5 there is **no depth from which the free-extension recursion holds non-trivially, for
any D**: the two-summand mechanism (closure = free extension of the one below, above a base depth)
needs the closure to have fallen behind the ideal, and for three summands no degree `>= 7` has fallen
behind at k <= 5.

### Task 3: where each degree falls behind

| D | first depth with V_D != I cap M_<=D (all 12 instances) | label |
|---|---|---|
| 5 | k = 2 (V_5 = 0, I cap M_<=5 != 0) | exact |
| 6 | k = 2 (V_6 = span(g)) | exact |
| 7 | not at k <= 5; k = 6: see "Depth 6" | exact for k <= 5 |
| 8 | not at k <= 5 (follows from D = 7) | exact for k <= 5 |

### Depth 6 (k = 6, D = 7): run on one instance, V_7 falls behind

Outcome (2026-11-10, **exact**, one instance, one engine): on p7681-r (p = 7681, a = 1876, b = 2678,
t = 6223, 48 solutions) the degree-7 closure at k = 6 has dimension 19,958 with profile
(0, 0, 0, 0, 0, 3, 2375, 19958), `closure_check = closed-sweep`, against 62,956 for the ideal part
(63,004 columns, evaluation rank 48), and 1 is not in it. This is the second reading below:
V_7 falls behind at k = 6 and LFD(6) >= 8 on this instance. The degree-6 closure at k = 6 is span(g)
(rank 1). The run took 1381 s and 5.7 GB on a rented machine's CPU (the local container kept
restarting); rows are in `closures.jsonl` (instance p7681-r, k = 6, D = 6 and 7), and `table.md` /
`summary.json` were regenerated. The text below is the plan as written before the run.


`|M_<=7| = 63004` at 18 chain variables; the compact engine needs about 2 GB of stored residues at
its peak plus about 1.5 GB of work arrays, so it must run under the heavy lock. It was queued at
03:31 on instance p7681-r (`|Z|` small, so the sandwich target is cheap to compute) and was still
waiting for the lock when this report was written (the lock holder runs a five-job queue). The queued
command (detached, log in the scratchpad `m3depth/k6.log`) appends rows to `closures-k6.jsonl`;
then `python3 scripts/ic_bounds_m3_depth_analyze.py` reads them (`table.md`, `summary.json`).
Reading of the outcome, decided in advance:
* `closure_check = rank=dim(I cap M_<=D)` and `eval_rank_le_j[6] = |Z|`: LFD(6) = 7 on p7681-r
  (exact, local-to-global), and D = 7 has not fallen behind at k = 6, unlike two summands where the
  LFD-degree closure falls behind exactly at k = 6.
* `closure_check = closed-sweep` with `rank < dim(I cap M_<=7)`: V_7 falls behind at k = 6 and
  LFD(6) >= 8 on p7681-r (exact lower bound; the analyzer reports `lower-bound-deficient-closure`).
  The exact value would need D = 8 at k = 6 (`|M_<=8| = 106762`), out of reach here. The recursion
  could then be tested only from k = 6 to 7, which is also out of reach (`|M_<=7| = 198440`).
Expected cost: 2 to 4 hours with 2 threads (k = 5 took 4 to 8 minutes single threaded; the work grows
roughly like `n^3`).

### Predictions for k = 6, 7, 8

No extrapolation is supported by these data. The recursion does not hold non-trivially at any
depth, so it cannot carry the profile forward; the only proven statements are `LFD >= 7` for every
k >= 3 (proof, generic targets) and `LFD <= 3k + 1` = 19, 22, 25 (kernel). Labelled guesses:
* **conjecture (weak)**: LFD(6) = 7 for generic targets at primes with 64 | p - 1. Support: V_7 =
  I cap M_<=7 with room to spare at k = 2..5 on all 12 instances (the evaluation of `M_<=5` already
  has full rank `|Z|` at k = 5 on all 12 instances), and `|M_<7| = 31180` exceeds the expected
  `|Z|` (about `4 * 8^6 / p`, at most about 5400 for p >= 193). Against: for two summands the
  LFD-degree closure falls behind at k = 6 even on unit-ideal instances (`|Z| = 0`), so the
  `|Z|` count is not what decides it. The queued run decides this for one instance.
* k = 7, 8: **no prediction**. The degree-7 certificate needs `|M_<7| >= |Z|`; with `|M_<7|` =
  82160 and 190051 at 21 and 24 variables and `|Z|` about `4 * 8^k / p`, it fails for the smallest
  admissible primes at k = 8 (p = 257: `|Z|` about 2.6e5). That is a limit of the certificate, not a
  lower bound on the LFD.

## Comparison with two summands (same instances, same factor base)

`m2.jsonl` (`scripts/ic_bounds_m3_depth_m2.py`, engine `csd.ic_curve_sweep.prime_local_profile`,
system `{2k chains, NF(S_3(x1, x2, t))}`, same `(p, a, b, t)`): exact LFD 4 at k = 2..5 on all 8
instances, and at k = 6 (9 instances with 64 | p - 1) 5 on 8 of them and 4 on p449-pl (k = 7 is capped
by the column budget of that engine). Except p449-pl at k = 6 this matches the two-summand
conjecture `LFD = floor((k+14)/4)` (stated for generic systems; exact on the ic-bounds-fall data,
k = 2..13); p449-pl is a counterexample to that formula as a statement about every instance.

| k (same H, same factor base) | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|
| m = 2, exact | 4 (12/12) | 4 (12/12) | 4 (12/12) | 4 (12/12) | 5 (8/9), 4 (1/9) |
| m = 3 direct, exact | 7 (12/12) | 7 (12/12) | 7 (12/12) | 7 (12/12) | >= 8 on one instance (p7681-r) |
| m = 3, variables 3k / m = 2 variables 2k | 6/4 | 9/6 | 12/8 | 15/10 | 18/12 |

* At equal factor base (equal k) the three-summand system sits 3 degrees above the two-summand one
  at k = 2..5 on every paired instance (exact).
* At equal number of chain variables (3k = 2k'): m = 3 at k = 2, 4 (6, 12 variables) has LFD 7,
  m = 2 at k' = 3, 6 has LFD 4, 5 (exact on these instances).
* Mechanism difference (exact on these data): for m = 2 the input has degree 2 and the fall is at
  `deg + 2 = 4`; for m = 3 the input has degree 6 and the fall is at `deg + 1 = 7`. For m = 2 the
  closure at the LFD falls behind the ideal from k = 6 (degree 4, giving LFD 5); for m = 3 the
  degree-7 closure equals the ideal through k = 5.

What the data do **not** establish: any value at k >= 6 beyond what the depth-6 run below shows;
any statement for other curves, primes or targets than the 12 instances (k = 5), the 36 systems here at k = 2..4 and the
234 systems of `experiments/ic-bounds-m3` at k = 3, 4; any relation between LFD and the cost of relation collection or of an ECDLP
attack (no cost claim is made; m = 3 also uses more variables and a factor base with three times
as many summands per relation).

## Not done / limits

* D = 8 at k = 6 (`|M_<=8| = 106762` at 18 variables) and any k = 7 closure at D = 7
  (`|M_<=7| = 198440`) do not fit this host (15 GB RAM).
* Only the dimension profile is tested against the recursion, not the pivot-level structure of a
  free extension.
* No Lean statement is produced here.

## Reproduce

```
cd /home/user/ml-cryptanalysis
ulimit -v 1500000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 scripts/ic_bounds_m3_depth.py selfcheck --instances p97-r p257-pl p7681-r --k 2 3
python3 scripts/ic_bounds_m3_depth.py selfcheck --instances p97-r p257-pl --k 4
python3 scripts/ic_bounds_m3_depth.py run --instances p97-r p193-r p257-r p7681-r p97-pl p193-pl p257-pl p7681-pl --k 2 3 4 --D 6 7 8
python3 scripts/ic_bounds_m3_depth.py run --instances p7681-r p257-r p97-r p193-r p97-pl p193-pl p257-pl p7681-pl --k 5 --D 6 7 --batch 384 --chunk 2048
python3 scripts/ic_bounds_m3_depth.py run --instances p7681-r p97-r --k 5 --D 7 --batch 384 --chunk 2048 --no-target --out experiments/ic-bounds-m3-depth/nosandwich-k5.jsonl
python3 scripts/ic_bounds_m3_depth.py run --instances p353-r p449-pl p577-r p641-pl --k 2 3 4 5 --D 6 7 --batch 384 --chunk 2048
python3 scripts/ic_bounds_m3_depth_m2.py --k 2 3 4 5 6 7
# depth 6 (heavy lock, about 2.5 GB resident):
flock /tmp/lfd-locks/heavy.lock sh -c 'ulimit -v 7000000; OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 python3 scripts/ic_bounds_m3_depth.py run --instances p7681-r --k 6 --D 6 7 --batch 512 --chunk 1024 --out experiments/ic-bounds-m3-depth/closures-k6.jsonl'
python3 scripts/ic_bounds_m3_depth_analyze.py    # summary.json, table.md
```

## Files

| file | role | sha256 (16) |
|---|---|---|
| `scripts/ic_bounds_m3_depth.py` | engine, instances, run, selfcheck | 971d41877d0f2863 |
| `scripts/ic_bounds_m3_depth_analyze.py` | exact LFD, recursion, falling-behind tables | 34c949a15d59150e |
| `scripts/ic_bounds_m3_depth_m2.py` | two-summand reference on the same instances | ce82f7fdcc524c67 |
| `experiments/ic-bounds-m3-depth/closures.jsonl` | one row per (instance, k, D): profile, ideal dims, eval ranks, |Z|, closure check | d5d66179ac855c8c |
| `experiments/ic-bounds-m3-depth/nosandwich-k5.jsonl` | k = 5, D = 7 closed by the full sweep (no sandwich) | 90467d215ac95cea |
| `experiments/ic-bounds-m3-depth/m2.jsonl` | m = 2 exact LFD rows | 2e7f357a456d734c |
| `experiments/ic-bounds-m3-depth/selfcheck-k23.json` | engine vs exact_profile, k = 2, 3 | 0173cabc596638bc |
| `experiments/ic-bounds-m3-depth/selfcheck-k4.json` | engine vs exact_profile, k = 4 | eb477a7b6aefbb70 |
| `experiments/ic-bounds-m3-depth/summary.json` | analyzer output | e3c6cf44cdf15af8 |
| `experiments/ic-bounds-m3-depth/table.md` | analyzer output (all tables) | 34922055078a4428 |
| `closures.jsonl` (k = 6 rows) | the depth-6 run on p7681-r, D = 6 and 7 | exact |
