# Last fall degree (LFD) of the direct m-summand systems, m = 3, 4, 5

Status 2026-11-12. Labels: **exact** (certified computation, command and certificate named), **proof** (written out in
`experiments/ic-bounds-m3/README.md`, not Lean-checked), **empirical**, **unverified**. Raw facts; interpretation only where
a prediction is scored. Code: `csd/ic_bounds_msum.py` (generators by evaluation, zero sets, ladders),
`csd/ic_bounds_msum_checks.py` (validation), `scripts/ic_bounds_msum_run.py` (driver: `survey`, `instance`, `validate`, `lfd`,
`lfd-file`, `plan`). Frozen predictions: `predictions-2026-11-12.json`, sha256
`deb4b304a065be94de56c051bec585f15a2ae7ece34c077266a6b390228eaad7`.

## Systems

**Prime chains.** `E: y^2 = x^3 + a x + b / F_p`, `2^k | p-1`, `H = <g0>` of order `2^k` (`g0 = subgroup_generator_value(p, k)`),
variables `u_{i,j}` (i < m, j < k) with the chain equations, `x_i = u_{i,0}`; input `g = NF(S_{m+1}(x_1..x_m, t))` modulo the
chains (exponent vector `a_i` mod `2^k`, chain mask = binary digits of `a_i` in block i, degree = number of one bits), the same
object as `csd.ic_bounds_m3.direct_generator` (m = 3) and `csd.ic_curve_sweep.prime_coset_generator` (m = 2). Closure `V_D`,
LFD rule, certificates: exactly those of `experiments/ic-bounds-m3/README.md` (unit ideal: least D with 1 in V_D; proper ideal:
first D with rank ev(M_<D on Z) = |Z| and dim V_D = |M_<=D| - |Z|, LFD = last strict fall at or below D). Z = zero set of the
evaluated `S_{m+1}` on the grid `H^m` (enumeration, independent of the closure). Engine: `CompactClosure` of
`scripts/ic_bounds_m3_depth.py` (generic through `MixedSpace(m, k, D)`), and for at most 1,500 columns also the dense
`csd.ic_fall_certificates.chain_degree_closure` (records say which).

**Binary / Koblitz.** `E: y^2 + xy = x^3 + a x^2 + b / F_(2^n)`, `V` an F_2-subspace of dimension k (basis canonical reduced),
`x_i = sum_j c_{i,j} v_j`, Boolean variable index `i*k + j`; the n coordinate functions of `S_{m+1}(x_1..x_m, x_R)` as Boolean
polynomials (field equations built in), i.e. `binary_descent` for m = 2, 3. Certificate and LFD rule:
`csd.ic_local_certificate.local_profile` (python, r <= 12) and `csd.ic_bounds_bin_big.local_profile_bin` (bit-packed, cpu, r >= 13):
same definition, cross-checked below.

**Generators by evaluation (new).** `S_n(y_1..y_n) = Res_X(S_{j+1}(y_1..y_j, X), S_{n-j+1}(y_{j+1}..y_n, X))`, numbers
substituted for every variable except X. Each factor is a polynomial in X of FIXED formal degree `2^(j-1)` (padded coefficient list;
`S_{j+1}(ys, X)` for j >= 3 by interpolation at `2^(j-1)+1` points) and the resultant is the determinant of the Sylvester matrix of
that formal size, so the value is a polynomial in the remaining variables. Prime: values on the `2^(mk)` grid points, then the
inverse DFT over F_p (size `2^k` per axis) gives the coefficients on exponent vectors (the chain quotient is the algebra of functions on
the grid when `2^k | p-1`). Binary: truth table on the `2^(mk)` assignments, Moebius transform (all n coordinates at once on uint64 words).
Defaults: split j = 2 for S_4, S_5; j = n/2 for S_6.

## Validation (exact statements; files in `validation/`)

| id | what was compared (exact statement) | n | result | file |
|---|---|---|---|---|
| V1 | recursion value vs `csd.ic_bounds_m3.s4_value` (m = 3) at 1200 random points over p = 17, 97, 1000003 (759 degenerate: equal coordinates, zeros, all four equal), closed-form fast path vs generic Sylvester determinant: equal at every point; ratio to `s4_value` is the single constant 1 for every prime. S_5, S_6, S_7: invariant under permutation of the arguments and under the choice of split j (exact equality, p = 97, 360 points); degree in each variable exactly 2^(n-2) (2, 4, 8, 16 for n = 3..6, finite differences, top coefficient may vanish at one random point so 4 trials each); vanishing at the x-coordinates of 959 random rational decompositions P_1 + ... + P_m = R, m = 2..5, p = 97, 1009, 7681 (0 nonzero) | 1200 + | ok | `summation-prime.json` |
| V1b | binary: recursion = `s4_binary_value` at 900 points (n = 5, 6, 8; a, b in Koblitz and random), fast = Sylvester; permutation/split invariance of S_5, S_6; vanishing at random decompositions m = 2..5 | 9 curves | ok | `summation-binary.json` |
| V2 | prime generators by grid + inverse DFT vs existing symbolic builders: m = 2 `prime_coset_generator` (20 instances, k = 2, 3, 4, 5), m = 3 `direct_generator` (18 instances, k = 2, 3, 4): same support, ratio {1} (every coefficient equal); sha256 of the m = 3 generator equal to the recorded `generator_hash` of `experiments/ic-bounds-m3-depth/closures.jsonl` for 48 (instance, k) pairs (12 instances, k = 2..5) | 86 | ok | `prime-generators.json` |
| V3 | binary: Moebius-transform equations vs `binary_descent` for m = 2 (90) and m = 3 (90): the canonical equation lists (sets of sorted monomial-mask tuples = the unique Boolean normal form) are EQUAL, not only their spans (n = 5..9, k = 2..4, Koblitz and random curves, random bases, two targets); Z equal to `binary_solutions`; `system_hash` and |Z| equal to the recorded rows of `experiments/ic-bounds-binary/results-m3-direct.jsonl` on 24 rows (k = 2..5) | 204 | ok | `binary-generators.json` |
| V3b | binary m = 4, 5 (no existing builder): truth tables for different splits and generic Sylvester vs fast path are identical arrays; planted decomposition in Z; rational decompositions (group law) subset of Z | 9 | ok | `binary-m45-consistency.json` |
| V3c | prime m = 3, 4, 5, all 30 case instances: rational decompositions found by the group law are a subset of Z (Z is larger: y outside F_p) | 30 | ok | `prime-Z-vs-group-law.json` |
| V4 | engines: dense `chain_degree_closure` vs `CompactClosure` without and with the dimension sandwich: per-D quotient rank and low profile identical (m = 2..5, 14 systems); binary `local_profile` (python) vs `local_profile_bin` (bit-packed, cpu): LFD, falls, quotient ranks identical (10 systems, m = 3, 4, 5) | 24 | ok | `engines-prime.json`, `engines-binary.json` |
| C1 | positive control m = 2 prime, k = 2..5: generator-list hash (= recorded `system_hash`), Z, exact LFD, status, strict falls equal to `experiments/ic-curve-sweep` (8 per k, unit and proper mixed: 3 unit + 5 proper at k = 2, 3, 2 unit + 6 proper at k = 4, 5; LFD 4 everywhere) | 32 | ok | `control-m2-prime.json` |
| C2 | positive control m = 2 binary, k = 4..7: system_hash, |Z|, LFD, status, falls equal to the recorded rows | 32 | ok | `control-m2-binary.json` |
| C3 | positive control m = 3 prime, the 12 instances of `experiments/ic-bounds-m3-depth` at k = 2, 3, 4 (new evaluation-based generator): LFD 7, |Z|, `V_D` profile and rank at every recorded D >= 6 equal to `closures.jsonl` | 36 | ok | `control-m3-prime.json` |
| C4 | positive control m = 3 binary, 28 recorded rows k = 2..5 (unit and proper, d_in 3 and 6, LFD 3..7): system_hash, |Z|, LFD, status, falls, certificate degree equal | 28 | ok | `control-m3-binary.json` |
| N1 | negative control prime (6 systems, 17 trials): one grid value altered. A zero of Z turned nonzero (5 trials, g' used with the TRUE Z): detected in 5 of 5 (assertion `unit ideal with nonempty zero set`, or V_D not inside the true ideal at D = LFD, or status capped); a nonzero turned zero (6): status capped in 6 of 6 (no certificate against the true Z); nonzero shifted (6): the zero set is unchanged so the ideal is the same: detected in 2 of 6 (LFD 10 instead of 7 at m = 3, k = 3 and a capped status elsewhere). The CONSISTENT altered system (g', Z') is a different system: its LFD equals the original in 14 of 17 trials (m = 4, 5 at k = 2 are forced to 9, 11 by d = m k, so a wrong generator does not change the value there) and differs in 3 (10 instead of 7 at m = 3, k = 3) | 17 | detected where it must | `negative-prime.json` |
| N2 | negative control binary (4 systems, 8 trials): one truth-table bit flipped at a zero of Z: certificate refuses (point does not vanish) 4 of 4; flipped at a nonzero point: result differs in 1 of 4 (the others are a different system with the same LFD) | 8 | ok | `negative-binary.json` |

What is NOT covered: k >= 5 prime against the recorded k = 5 closures (generator hash only, V2), binary m = 4, 5 against any pre-existing record (none exists), the GPU engine (not run).
Re-derivations of N1 counts are in the JSON (`detected`, `consistent_altered`).

## Predictions and what happened

Predictions file `predictions-2026-11-12.json` (sha256 `deb4b304a065be94de56c051bec585f15a2ae7ece34c077266a6b390228eaad7`) was written after `survey` (generators, degrees, |Z|; no closure) and before any closure of the
cases below. **Disclosure:** before the freeze two smoke-test closures ran (prime m = 4, k = 2 and m = 5, k = 2, p = 97, instance seed 1, results LFD 9 and 11); those instances are not among the cases (seeds 10..15), but the
k = 2 predictions for m = 4, 5 are informed by them. After the freeze only additions were made to `csd/ic_bounds_msum.py` (`binary_degree_and_zeros`) and the run script (`instance`, `lfd-file`, `plan`, `build_binary_light`); current sha256 of
`csd/ic_bounds_msum.py` a1bc04e7..., `csd/ic_bounds_msum_checks.py` 70fdd2ab..., `scripts/ic_bounds_msum_run.py` 9f8aa0a8... (the frozen file records the pre-addition hashes of the first and last).

Scoring (`score.py` -> `scores.json`, per-case rows there):

| prediction | scope | hit | miss | note |
|---|---|---|---|---|
| G1/(b) prime LFD = d + 1 (equal to m k + 1 here; a coincides) | 36 prime cases (18 m = 3, 6 each m = 4 k = 2, 3 and m = 5 k = 2) | 36 | 0 | values 7 (m = 3, k = 2, 3, 4), 9, 13, 11 |
| G2 LFD <= m k + 1; LFD >= d + 1 where Proposition 1 applies | 36 | 36 | 0 | Prop. 1 variables > 0 in all 36 |
| G3 exactly one strict fall, at d + 1 | 36 | 36 | 0 | |
| G6 m = 3 LFD independent of k = 2, 3, 4 | 18 | 18 | 0 | |
| G4/(c) binary LFD <= min(r, d_in + 2); never r + 1 | 48 exact binary cases | 48 | 0 | 4 cases not run (r = 16: `B-m4n9k4-*`, column budget), not scored |
| G5/(c) point: LFD = min(r, d_in + 1), scored where d_in = r | 32 cases | 29 | 3 | misses: `B-m4n6k2-koblitz0-ran-ord` (7 vs 8), `B-m4n6k2-koblitz1-inv-ran` (7 vs 8), `B-m4n7k2-ordA-ran-pla` (7 vs 8); 29 / 32 = 91 percent >= the 50 percent threshold |
| G7 m = 4, k = 4 in [13, 17], point 13; G8 m = 5, k = 3 = 16 | GPU cells | not run | | `gpu/` |

(a) and (b) agree on every prime cell that was run, so the cells do not discriminate them (d = m k for m >= 4, k <= 3; for m = 3 d = 6 < m k at k = 3, 4 and (a) is only an upper bound there). The prime values at d = m k
(9, 13, 11, and 7 at m = 3, k = 2) are also forced by Proposition 1 (proof, hypotheses checked per instance) together with the kernel-checked upper bound m k + 1 (taken from `experiments/ic-bounds-m3/README.md`; not re-checked here): the runs confirm
the implementation and the hypotheses, they are not an independent measurement of an unknown number. A binary cell with an exact value not predicted by the point rule: m = 3 order-two targets (d_in = 3, LFD 4 or 5; not scored).

## Results (every case that was run; `results.md` has the same tables, `records/*.json` one file per case)

### Prime chains (all cases)

| m | k | case | p | target | deg g | terms | \|Z\| | class | LFD | cert. degree | strict falls | kernel bound | status | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3 | 2 | P-m3k2-p193-r | 193 | random-rational | 6 | 61 | 0 | unit | 7 | 7 | [7] | 7 | exact-unit | 0.01 |
| 3 | 2 | P-m3k2-p257-pl | 257 | planted-H4 | 6 | 64 | 1 | proper | 7 | 7 | [7] | 7 | exact-proper-local | 0.02 |
| 3 | 2 | P-m3k2-p449-r | 449 | random-rational | 6 | 64 | 0 | unit | 7 | 7 | [7] | 7 | exact-unit | 0.01 |
| 3 | 2 | P-m3k2-p641-pl | 641 | planted-H4 | 6 | 64 | 3 | proper | 7 | 7 | [7] | 7 | exact-proper-local | 0.02 |
| 3 | 2 | P-m3k2-p7681-r | 7681 | random-rational | 6 | 64 | 0 | unit | 7 | 7 | [7] | 7 | exact-unit | 0.01 |
| 3 | 2 | P-m3k2-p97-pl | 97 | planted-H4 | 6 | 64 | 10 | proper | 7 | 7 | [7] | 7 | exact-proper-local | 0.02 |
| 3 | 3 | P-m3k3-p193-r | 193 | random-rational | 6 | 125 | 0 | unit | 7 | 7 | [7] | 10 | exact-unit | 0.32 |
| 3 | 3 | P-m3k3-p257-pl | 257 | planted-H4 | 6 | 125 | 22 | proper | 7 | 7 | [7] | 10 | exact-proper-local | 0.68 |
| 3 | 3 | P-m3k3-p449-r | 449 | random-rational | 6 | 125 | 0 | unit | 7 | 7 | [7] | 10 | exact-unit | 0.31 |
| 3 | 3 | P-m3k3-p641-pl | 641 | planted-H4 | 6 | 125 | 6 | proper | 7 | 7 | [7] | 10 | exact-proper-local | 0.78 |
| 3 | 3 | P-m3k3-p7681-r | 7681 | random-rational | 6 | 125 | 0 | unit | 7 | 7 | [7] | 10 | exact-unit | 0.34 |
| 3 | 3 | P-m3k3-p97-pl | 97 | planted-H4 | 6 | 122 | 3 | proper | 7 | 7 | [7] | 10 | exact-proper-local | 0.81 |
| 3 | 4 | P-m3k4-p193-r | 193 | random-rational | 6 | 125 | 12 | proper | 7 | 7 | [7] | 13 | exact-proper-local | 4.61 |
| 3 | 4 | P-m3k4-p257-pl | 257 | planted-H4 | 6 | 125 | 67 | proper | 7 | 7 | [7] | 13 | exact-proper-local | 5.29 |
| 3 | 4 | P-m3k4-p449-r | 449 | random-rational | 6 | 125 | 12 | proper | 7 | 7 | [7] | 13 | exact-proper-local | 4.63 |
| 3 | 4 | P-m3k4-p641-pl | 641 | planted-H4 | 6 | 125 | 7 | proper | 7 | 7 | [7] | 13 | exact-proper-local | 4.25 |
| 3 | 4 | P-m3k4-p7681-r | 7681 | random-rational | 6 | 125 | 0 | unit | 7 | 7 | [7] | 13 | exact-unit | 4.75 |
| 3 | 4 | P-m3k4-p97-pl | 97 | planted-H4 | 6 | 125 | 56 | proper | 7 | 7 | [7] | 13 | exact-proper-local | 4.62 |
| 4 | 2 | P-m4k2-p193-r | 193 | random-rational | 8 | 256 | 0 | unit | 9 | 9 | [9] | 9 | exact-unit | 0.07 |
| 4 | 2 | P-m4k2-p257-pl | 257 | planted-H4 | 8 | 256 | 19 | proper | 9 | 9 | [9] | 9 | exact-proper-local | 0.14 |
| 4 | 2 | P-m4k2-p449-r | 449 | random-rational | 8 | 255 | 0 | unit | 9 | 9 | [9] | 9 | exact-unit | 0.11 |
| 4 | 2 | P-m4k2-p641-pl | 641 | planted-H4 | 8 | 256 | 32 | proper | 9 | 9 | [9] | 9 | exact-proper-local | 0.15 |
| 4 | 2 | P-m4k2-p7681-r | 7681 | random-rational | 8 | 256 | 0 | unit | 9 | 9 | [9] | 9 | exact-unit | 0.07 |
| 4 | 2 | P-m4k2-p97-pl | 97 | planted-H4 | 8 | 256 | 4 | proper | 9 | 9 | [9] | 9 | exact-proper-local | 0.16 |
| 4 | 3 | P-m4k3-p193-r | 193 | random-rational | 12 | 4084 | 0 | unit | 13 | 13 | [13] | 13 | exact-unit | 13.32 |
| 4 | 3 | P-m4k3-p257-pl | 257 | planted-H4 | 12 | 4078 | 92 | proper | 13 | 13 | [13] | 13 | exact-proper-local | 14.45 |
| 4 | 3 | P-m4k3-p449-r | 449 | random-rational | 12 | 4096 | 6 | proper | 13 | 13 | [13] | 13 | exact-proper-local | 13.55 |
| 4 | 3 | P-m4k3-p641-pl | 641 | planted-H4 | 12 | 4096 | 43 | proper | 13 | 13 | [13] | 13 | exact-proper-local | 14.59 |
| 4 | 3 | P-m4k3-p7681-r | 7681 | random-rational | 12 | 4096 | 0 | unit | 13 | 13 | [13] | 13 | exact-unit | 14.97 |
| 4 | 3 | P-m4k3-p97-pl | 97 | planted-H4 | 12 | 4036 | 152 | proper | 13 | 13 | [13] | 13 | exact-proper-local | 17.83 |
| 5 | 2 | P-m5k2-p193-r | 193 | random-rational | 10 | 1024 | 0 | unit | 11 | 11 | [11] | 11 | exact-unit | 3.79 |
| 5 | 2 | P-m5k2-p257-pl | 257 | planted-H4 | 10 | 1024 | 136 | proper | 11 | 11 | [11] | 11 | exact-proper-local | 6.63 |
| 5 | 2 | P-m5k2-p449-r | 449 | random-rational | 10 | 1024 | 0 | unit | 11 | 11 | [11] | 11 | exact-unit | 3.87 |
| 5 | 2 | P-m5k2-p641-pl | 641 | planted-H4 | 10 | 1024 | 75 | proper | 11 | 11 | [11] | 11 | exact-proper-local | 7.2 |
| 5 | 2 | P-m5k2-p7681-r | 7681 | random-rational | 10 | 1024 | 0 | unit | 11 | 11 | [11] | 11 | exact-unit | 3.28 |
| 5 | 2 | P-m5k2-p97-pl | 97 | planted-H4 | 10 | 1024 | 136 | proper | 11 | 11 | [11] | 11 | exact-proper-local | 6.64 |

### Binary (all cases)

| m | n | k | r | case | curve | base | target | d_in | eqs | \|Z\| | class | LFD | cert. degree | strict falls | r+1 | status | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3 | 6 | 2 | 6 | B-m3n6k2-koblitz0-ran-ord | koblitz0 | random | order-two | 3 | 6 | 10 | proper | 4 | 4 | [4] | 7 | exact-proper-local | 0.01 |
| 3 | 6 | 2 | 6 | B-m3n6k2-koblitz1-inv-ran | koblitz1 | invariant | random | 6 | 6 | 0 | unit | 6 | 6 | [4, 5, 6] | 7 | exact-unit | 0.03 |
| 3 | 6 | 2 | 6 | B-m3n6k2-ordA-ran-pla | ordA | random | planted | 6 | 6 | 10 | proper | 6 | 6 | [5, 6] | 7 | exact-proper-local | 0.02 |
| 3 | 6 | 2 | 6 | B-m3n6k2-ordB-inv-ran | ordB | invariant | random | 6 | 6 | 3 | proper | 6 | 6 | [5, 6] | 7 | exact-proper-local | 0.02 |
| 3 | 8 | 3 | 9 | B-m3n8k3-koblitz0-ran-ord | koblitz0 | random | order-two | 3 | 8 | 22 | proper | 4 | 4 | [3, 4] | 10 | exact-proper-local | 0.05 |
| 3 | 8 | 3 | 9 | B-m3n8k3-koblitz1-inv-ran | koblitz1 | invariant | random | 6 | 8 | 0 | unit | 6 | 6 | [6] | 10 | exact-unit | 0.1 |
| 3 | 8 | 3 | 9 | B-m3n8k3-ordA-ran-pla | ordA | random | planted | 6 | 8 | 22 | proper | 7 | 7 | [7] | 10 | exact-proper-local | 0.23 |
| 3 | 8 | 3 | 9 | B-m3n8k3-ordB-inv-ran | ordB | invariant | random | 6 | 8 | 0 | unit | 6 | 6 | [6] | 10 | exact-unit | 0.06 |
| 3 | 9 | 4 | 12 | B-m3n9k4-koblitz0-ran-ord | koblitz0 | random | order-two | 3 | 8 | 52 | proper | 4 | 5 | [4] | 13 | exact-proper-local | 0.85 |
| 3 | 9 | 4 | 12 | B-m3n9k4-koblitz1-inv-ran | koblitz1 | random(no in | random | 6 | 9 | 12 | proper | 7 | 7 | [7] | 13 | exact-proper-local | 7.11 |
| 3 | 9 | 4 | 12 | B-m3n9k4-ordA-ran-pla | ordA | random | planted | 6 | 9 | 12 | proper | 7 | 7 | [7] | 13 | exact-proper-local | 6.82 |
| 3 | 9 | 4 | 12 | B-m3n9k4-ordB-inv-ran | ordB | random(no in | random | 6 | 9 | 9 | proper | 7 | 7 | [7] | 13 | exact-proper-local | 6.39 |
| 3 | 8 | 5 | 15 | B-m3n8k5-koblitz0-ran-ord | koblitz0 | random | order-two | 3 | 8 | 220 | proper | 5 | 6 | [3, 4, 5] | 16 | exact-proper-local | 22.38 |
| 3 | 8 | 5 | 15 | B-m3n8k5-koblitz1-inv-ran | koblitz1 | invariant | random | 6 | 8 | 0 | unit | 7 | 7 | [7] | 16 | exact-unit | 57.89 |
| 3 | 8 | 5 | 15 | B-m3n8k5-ordA-ran-pla | ordA | random | planted | 6 | 8 | 193 | proper | 7 | 7 | [7] | 16 | exact-proper-local | 61.68 |
| 3 | 8 | 5 | 15 | B-m3n8k5-ordB-inv-ran | ordB | invariant | random | 6 | 8 | 482 | proper | 7 | 7 | [7] | 16 | exact-proper-local | 66.38 |
| 4 | 6 | 2 | 8 | B-m4n6k2-koblitz0-ran-ord | koblitz0 | random | order-two | 8 | 6 | 19 | proper | 7 | 7 | [7] | 9 | exact-proper-local | 0.06 |
| 4 | 6 | 2 | 8 | B-m4n6k2-koblitz1-inv-ran | koblitz1 | invariant | random | 8 | 4 | 206 | proper | 7 | 7 | [7] | 9 | exact-proper-local | 0.04 |
| 4 | 6 | 2 | 8 | B-m4n6k2-ordA-ran-pla | ordA | random | planted | 8 | 6 | 19 | proper | 8 | 8 | [8] | 9 | exact-proper-local | 0.06 |
| 4 | 6 | 2 | 8 | B-m4n6k2-ordB-inv-ran | ordB | invariant | random | 8 | 6 | 0 | unit | 8 | 8 | [7, 8] | 9 | exact-unit | 0.08 |
| 4 | 7 | 2 | 8 | B-m4n7k2-koblitz0-ran-ord | koblitz0 | random | order-two | 8 | 7 | 0 | unit | 8 | 8 | [7, 8] | 9 | exact-unit | 0.07 |
| 4 | 7 | 2 | 8 | B-m4n7k2-koblitz1-inv-ran | koblitz1 | random(no in | random | 8 | 7 | 0 | unit | 8 | 8 | [8] | 9 | exact-unit | 0.03 |
| 4 | 7 | 2 | 8 | B-m4n7k2-ordA-ran-pla | ordA | random | planted | 8 | 7 | 19 | proper | 7 | 7 | [7] | 9 | exact-proper-local | 0.06 |
| 4 | 7 | 2 | 8 | B-m4n7k2-ordB-inv-ran | ordB | random(no in | random | 8 | 7 | 0 | unit | 8 | 8 | [7, 8] | 9 | exact-unit | 0.08 |
| 4 | 8 | 3 | 12 | B-m4n8k3-koblitz0-ran-ord | koblitz0 | random | order-two | 8 | 8 | 12 | proper | 9 | 9 | [9] | 13 | exact-proper-local | 13.7 |
| 4 | 8 | 3 | 12 | B-m4n8k3-koblitz1-inv-ran | koblitz1 | invariant | random | 12 | 8 | 0 | unit | 12 | 12 | [11, 12] | 13 | exact-unit | 24.99 |
| 4 | 8 | 3 | 12 | B-m4n8k3-ordA-ran-pla | ordA | random | planted | 12 | 8 | 140 | proper | 12 | 12 | [11, 12] | 13 | exact-proper-local | 32.72 |
| 4 | 8 | 3 | 12 | B-m4n8k3-ordB-inv-ran | ordB | invariant | random | 12 | 8 | 0 | unit | 12 | 12 | [11, 12] | 13 | exact-unit | 25.46 |
| 4 | 9 | 3 | 12 | B-m4n9k3-koblitz0-ran-ord | koblitz0 | random | order-two | 8 | 9 | 12 | proper | 9 | 9 | [9] | 13 | exact-proper-local | 14.54 |
| 4 | 9 | 3 | 12 | B-m4n9k3-koblitz1-inv-ran | koblitz1 | invariant | random | 12 | 9 | 0 | unit | 12 | 12 | [11, 12] | 13 | exact-unit | 23.75 |
| 4 | 9 | 3 | 12 | B-m4n9k3-ordA-ran-pla | ordA | random | planted | 12 | 9 | 12 | proper | 12 | 12 | [12] | 13 | exact-proper-local | 17.68 |
| 4 | 9 | 3 | 12 | B-m4n9k3-ordB-inv-ran | ordB | invariant | random | 12 | 9 | 80 | proper | 12 | 12 | [11, 12] | 13 | exact-proper-local | 36.06 |
| 4 | 9 | 4 | 16 | B-m4n9k4-koblitz0-ran-ord | | | | 8 | | 66 | | | | | | not-run: estimated columns 58651 at D=10 exceed --max-columns 40000 | |
| 4 | 9 | 4 | 16 | B-m4n9k4-koblitz1-inv-ran | | | | 12 | | 115 | | | | | | not-run: estimated columns 65519 at D=14 exceed --max-columns 40000 | |
| 4 | 9 | 4 | 16 | B-m4n9k4-ordA-ran-pla | | | | 12 | | 312 | | | | | | not-run: estimated columns 65519 at D=14 exceed --max-columns 40000 | |
| 4 | 9 | 4 | 16 | B-m4n9k4-ordB-inv-ran | | | | 12 | | 12 | | | | | | not-run: estimated columns 65519 at D=14 exceed --max-columns 40000 | |
| 5 | 6 | 2 | 10 | B-m5n6k2-koblitz0-ran-ord | koblitz0 | random | order-two | 10 | 6 | 211 | proper | 10 | 10 | [9, 10] | 11 | exact-proper-local | 1.34 |
| 5 | 6 | 2 | 10 | B-m5n6k2-koblitz1-inv-ran | koblitz1 | invariant | random | 10 | 6 | 0 | unit | 10 | 10 | [9, 10] | 11 | exact-unit | 1.16 |
| 5 | 6 | 2 | 10 | B-m5n6k2-ordA-ran-pla | ordA | random | planted | 10 | 6 | 136 | proper | 10 | 10 | [10] | 11 | exact-proper-local | 0.89 |
| 5 | 6 | 2 | 10 | B-m5n6k2-ordB-inv-ran | ordB | invariant | random | 10 | 6 | 136 | proper | 10 | 10 | [9, 10] | 11 | exact-proper-local | 1.7 |
| 5 | 7 | 2 | 10 | B-m5n7k2-koblitz0-ran-ord | koblitz0 | random | order-two | 10 | 7 | 136 | proper | 10 | 10 | [9, 10] | 11 | exact-proper-local | 1.58 |
| 5 | 7 | 2 | 10 | B-m5n7k2-koblitz1-inv-ran | koblitz1 | random(no in | random | 10 | 7 | 0 | unit | 10 | 10 | [9, 10] | 11 | exact-unit | 1.13 |
| 5 | 7 | 2 | 10 | B-m5n7k2-ordA-ran-pla | ordA | random | planted | 9 | 7 | 136 | proper | 9 | 9 | [9] | 11 | exact-proper-local | 0.83 |
| 5 | 7 | 2 | 10 | B-m5n7k2-ordB-inv-ran | ordB | random(no in | random | 10 | 7 | 0 | unit | 10 | 10 | [9, 10] | 11 | exact-unit | 1.44 |
| 5 | 8 | 3 | 15 | B-m5n8k3-koblitz0-ran-ord | koblitz0 | random | order-two | 15 | 8 | 1126 | proper | 15 | 15 | [14, 15] | 16 | exact-proper-local | 958.17 |
| 5 | 8 | 3 | 15 | B-m5n8k3-koblitz1-inv-ran | koblitz1 | invariant | random | 15 | 8 | 0 | unit | 15 | 15 | [14, 15] | 16 | exact-unit | 660.71 |
| 5 | 8 | 3 | 15 | B-m5n8k3-ordA-ran-pla | ordA | random | planted | 15 | 8 | 215 | proper | 15 | 15 | [14, 15] | 16 | exact-proper-local | 951.52 |
| 5 | 8 | 3 | 15 | B-m5n8k3-ordB-inv-ran | ordB | invariant | random | 15 | 8 | 0 | unit | 15 | 15 | [14, 15] | 16 | exact-unit | 737.03 |
| 5 | 9 | 3 | 15 | B-m5n9k3-koblitz0-ran-ord | koblitz0 | random | order-two | 15 | 9 | 801 | proper | 15 | 15 | [14, 15] | 16 | exact-proper-local | 898.96 |
| 5 | 9 | 3 | 15 | B-m5n9k3-koblitz1-inv-ran | koblitz1 | invariant | random | 15 | 9 | 0 | unit | 15 | 15 | [15] | 16 | exact-unit | 294.6 |
| 5 | 9 | 3 | 15 | B-m5n9k3-ordA-ran-pla | ordA | random | planted | 15 | 9 | 71 | proper | 15 | 15 | [14, 15] | 16 | exact-proper-local | 884.9 |
| 5 | 9 | 3 | 15 | B-m5n9k3-ordB-inv-ran | ordB | invariant | random | 15 | 9 | 0 | unit | 15 | 15 | [14, 15] | 16 | exact-unit | 689.63 |

### Histogram of exact LFD per (family, m, k)

| family | m | k | LFD: count |
|---|---|---|---|
| binary | 3 | 2 | 4: 1, 6: 3 |
| binary | 3 | 3 | 4: 1, 6: 2, 7: 1 |
| binary | 3 | 4 | 4: 1, 7: 3 |
| binary | 3 | 5 | 5: 1, 7: 3 |
| binary | 4 | 2 | 7: 3, 8: 5 |
| binary | 4 | 3 | 12: 6, 9: 2 |
| binary | 4 | 4 | not-run: estimated c: 4 |
| binary | 5 | 2 | 10: 7, 9: 1 |
| binary | 5 | 3 | 15: 8 |
| prime | 3 | 2 | 7: 6 |
| prime | 3 | 3 | 7: 6 |
| prime | 3 | 4 | 7: 6 |
| prime | 4 | 2 | 9: 6 |
| prime | 4 | 3 | 13: 6 |
| prime | 5 | 2 | 11: 6 |


## Limits

* Prime cells: k <= 4 (m = 3), k <= 3 (m = 4), k = 2 (m = 5); 6 instances per cell (3 planted-H4 targets, 3 random rational), p in {97, 193, 257, 449, 641, 7681}; instance seeds 10..15. Planted targets at k = 2 lie in the subgroup generated by the factor-base points,
  so |Z| is large (up to 136 of 1024 at m = 5). Not run: prime m = 4, k = 4 (65,536 monomials), m = 5, k = 3 (32,768), m = 3, k >= 5 beyond the recorded closures; binary r >= 16 (`B-m4n9k4-*`: 4 cases) and r = 18..24: see `gpu/PLAN.md`
  (the prime m-block GPU engine does not exist: `csd/ic_bounds_big.py` is two-block; what is missing is listed there; every GPU time is an extrapolation, unverified).
* The grid is a seeded structured sample, not a population sample. Binary instances: one curve family member per class (Koblitz (0,1), (1,1), two random ordinary curves), random or Frobenius-invariant bases (where an invariant space of that dimension exists; otherwise a random base, recorded in `base_kind`),
  targets order-two / random / planted (planted m-sum). No claim about relation yield or ECDLP cost.
* Z for m >= 4 is the zero set of the evaluated summation polynomial on the grid (independent of the closure, validated by symmetry, split independence, degrees, vanishing at actual decompositions, and group-law subset checks), not an independent algebraic enumeration of all (including non-rational) decompositions.
* The kernel upper bound m k + 1 and Proposition 1 are used as stated in the repository (written proofs / kernel-checked statement); nothing here re-checks them in Lean.
* The binary r = 15 cells (m = 5, k = 3; m = 3, k = 5) ran on the bit-packed engine on the cpu device (RSS 1.9 GB, under the heavy lock, up to 16 min per cell); cpu throughput was about 10 GFLOP/s (3 threads), below the earlier 26 GFLOP/s calibration.
* Pre-freeze smoke tests: see above. `B-m4n9k4-*` and `GP-*` are not results.
* This is not a claim about index-calculus speed: LFD is one proxy for the degree reached; the success criterion in AGENTS.md (iterative improvement of the baseline) is not addressed by these measurements.

## Commands to rerun (from the repository root; memory limits as in AGENT_RULES)

```
ulimit -v 7000000; export OMP_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3     # heavy lock for anything above 1.2 GB resident
python3 scripts/ic_bounds_msum_run.py survey --out experiments/ic-bounds-msummands/survey-prime.json --family prime     # generators, degrees, Z (no closure)
python3 scripts/ic_bounds_msum_run.py survey --family binary --out experiments/ic-bounds-msummands/survey-binary.json
python3 experiments/ic-bounds-msummands/make_predictions.py                    # write-once
python3 scripts/ic_bounds_msum_run.py validate                                 # all checks -> validation/ (about 11 min; control_m3_binary 7 min)
python3 scripts/ic_bounds_msum_run.py lfd --family prime --engine dense        # 36 prime cases (4 min)
python3 scripts/ic_bounds_msum_run.py lfd --family binary --max-columns 40000 --budget-min 115   # 52 binary cases, resumable (about 2 h: the 8 + 8 r = 15 m = 5 cells take 11-16 min each)
python3 experiments/ic-bounds-msummands/score.py; python3 experiments/ic-bounds-msummands/make_tables.py
python3 scripts/ic_bounds_msum_run.py plan                                     # gpu/ (about 10 min: builds / evaluates the r = 16..20 instances)
```
Cell names, seeds and primes are in `prime_cells()`, `binary_cells()`, `GPU_BIN`, `GPU_PRIME` of the run script. Logs of this session were in the scratchpad (`msum/*.log`), not kept.

## Files

`predictions-2026-11-12.json` (+ `make_predictions.py`, `survey-prime.json`, `survey-binary.json`), `records/` (one JSON per case: instance parameters, `instance_hash`, `generator_hash` (prime) or `system_hash` (binary), `points_sha256`, profile per D, certificate),
`validation/`, `scores.json`, `results.md`, `gpu/` (`PLAN.md`, `cost.json`, `cells.json`, `chain_msum.py`, `specs/`).


## Pod run (2026-11-12; plan `gpu/RUN.md`, written before the run)

One H200 pod (`icb-msum`, `mv5hxv2uqelul3`), created 06:15:51Z, terminated by the copy loop at 09:08:59Z.
Records: `fetched/` and `fetched-p2/` (archives as copied), unpacked into `pod-records/<chain>/`. Scoring:
`python3 score_pod.py` → `pod-scores.json`, `pod-table.md`.

Operational record:
- First job (`multi_msum.py`): the binary chain ran as planned; the four prime chains failed at once with
  `FileNotFoundError: /workspace/code/csd/_elim.c` (the code bundle held only `csd/*.py`). Their failed records are in
  `fetched/msum-pm*.tgz`.
- Second job (`gpu/rerun_prime.py`, bundle with `csd/*.c`): the same prime systems. Stopped by a TERM signal at
  about 09:05Z after the last m = 4, k = 4 instance finished and had been copied off; at that moment only the second
  m = 3, k = 6 instance was running (started about 08:25Z), and the first had ended at its 2-hour limit at D = 8.
  The second m = 3, k = 6 instance has no record. Not run at all (plan): binary r = 20 except one cell, r = 24, prime
  m = 3, k = 7.
- Binary degree caps were raised by one after the freeze (RUN.md); G-m3n10k6-koblitz0-ran-ord is capped even at the
  raised cap 6.
- Identity: prime records carry the generator hash and equal the frozen one in `gpu/cost.json`. Binary instances were
  rebuilt on the pod and their system hash compared with the frozen one before each run where a frozen hash exists
  (r = 16, 18); for the r = 20, 21 cells no instance was built before the run, so there is no frozen hash to compare.

| case | family | m | k | deg g or d_in | r | \|Z\| | status | LFD | strict falls | s | hash equal | predictions |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| G-m4n9k4-koblitz0-ran-ord | binary | 4 | 4 | 8 | 16 | 66 | exact-proper-local | 9 | [9] | 17.47 | True | upper True, point 9 True |
| G-m4n9k4-koblitz1-inv-ran | binary | 4 | 4 | 12 | 16 | 115 | exact-proper-local | 13 | [13] | 15.88 | True | upper True, point 13 True |
| G-m4n9k4-ordA-ran-pla | binary | 4 | 4 | 12 | 16 | 312 | exact-proper-local | 13 | [13] | 15.96 | True | upper True, point 13 True |
| G-m3n10k6-koblitz0-ran-ord | binary | 3 | 6 | 3 | 18 | 412 | capped | - | [4, 5] | 11.98 | True |  |
| G-m3n10k6-koblitz1-inv-ran | binary | 3 | 6 | 6 | 18 | 1227 | exact-proper-local | 7 | [7] | 17.01 | True | upper True, point 7 True |
| G-m3n10k6-ordA-ran-pla | binary | 3 | 6 | 6 | 18 | 375 | exact-proper-local | 7 | [7] | 14.28 | True | upper True, point 7 True |
| G-m4n10k5-koblitz0-ran-ord | binary | 4 | 5 | 8 | 20 | 2086 | exact-proper-local | 9 | [9] | 1244.46 | no frozen hash | upper True, point 9 True |
| G-m3n11k7-koblitz0-ran-ord | binary | 3 | 7 | 6 | 21 | None | exact-proper-local | 6 | [3, 4, 5, 6] | 584.9 | no frozen hash | upper True, point 7 False |
| G-m3n11k7-koblitz1-inv-ran | binary | 3 | 7 | 6 | 21 | None | exact-proper-local | 7 | [7] | 154.17 | no frozen hash | upper True, point 7 True |
| G-m3n11k7-ordA-ran-pla | binary | 3 | 7 | 6 | 21 | None | exact-proper-local | 7 | [7] | 154.83 | no frozen hash | upper True, point 7 True |
| GP-m3k6-p7681-r | prime | 3 | 6 | None |  | None | timeout | - | None | None | None |  |
| GP-m4k4-p97-pl | prime | 4 | 4 | 12 |  | 416 | exact-proper-local | 13 | [13] | 3097.9 | True | G7 range True, (b)=13 True |
| GP-m4k4-p193-r | prime | 4 | 4 | 12 |  | 316 | exact-proper-local | 13 | [13] | 3010.0 | True | G7 range True, (b)=13 True |
| GP-m4k4-p257-pl | prime | 4 | 4 | 12 |  | 553 | exact-proper-local | 13 | [13] | 3187.1 | True | G7 range True, (b)=13 True |
| GP-m4k4-p449-r | prime | 4 | 4 | 12 |  | 60 | exact-proper-local | 13 | [13] | 2919.8 | True | G7 range True, (b)=13 True |
| GP-m4k4-p641-pl | prime | 4 | 4 | 12 |  | 163 | exact-proper-local | 13 | [13] | 2905.1 | True | G7 range True, (b)=13 True |
| GP-m4k4-p7681-r | prime | 4 | 4 | 12 |  | 0 | exact-unit | 13 | [13] | 2894.0 | True | G7 range True, (b)=13 True |
| GP-m5k3-p97-pl | prime | 5 | 3 | 15 |  | 96 | exact-proper-local | 16 | [16] | 707.4 | True | G8 (16) True |
| GP-m5k3-p193-r | prime | 5 | 3 | 15 |  | 305 | exact-proper-local | 16 | [16] | 790.3 | True | G8 (16) True |
| GP-m5k3-p257-pl | prime | 5 | 3 | 15 |  | 1 | exact-proper-local | 16 | [16] | 700.8 | True | G8 (16) True |
| GP-m5k3-p449-r | prime | 5 | 3 | 15 |  | 0 | exact-unit | 16 | [16] | 705.6 | True | G8 (16) True |
| GP-m5k3-p641-pl | prime | 5 | 3 | 15 |  | 736 | exact-proper-local | 16 | [16] | 1222.9 | True | G8 (16) True |
| GP-m5k3-p7681-r | prime | 5 | 3 | 15 |  | 0 | exact-unit | 16 | [16] | 712.3 | True | G8 (16) True |


### Scoring of the frozen predictions on the pod cells

| prediction | cases | hits |
|---|---|---|
| G7: prime m = 4, k = 4, LFD in [13, 17] | 6 | 6 |
| G7 point prediction (b): LFD = d + 1 = 13 | 6 | 6 (candidate (a), "mk + 1 = 17 attained", fails on all 6) |
| G8: prime m = 5, k = 3, LFD = 16 (forced) | 6 | 6 |
| G4: binary LFD <= min(r, d_in + 2) | 9 exact | 9 |
| binary point rule LFD = min(r, d_in + 1) | 9 exact | 8 (miss: G-m3n11k7-koblitz0-ran-ord, LFD 6 against 7) |
| G9: capped or failed cases count as misses | G-m3n10k6-koblitz0-ran-ord capped (LFD >= 5); GP-m3k6-p7681-r timeout (LFD >= 8, already recorded) | |

Reading (empirical): on every prime cell computed, including the one where deg g < mk (m = 4, k = 4), the last fall
degree is deg g + 1 with a single strict fall, on proper and unit ideals alike. The binary point rule d_in + 1 holds on
8 of 9 new exact cases; the miss has an order-two target, where the degree falls stepwise (3, 4, 5, 6). For the three r = 21 cells
d_in is the theoretical bound 6 (the equations were not built before the run), so their point-rule scores use that
bound; for an order-two target the true d_in is smaller (3 at r = 18), so the scoring of that miss is itself uncertain.
