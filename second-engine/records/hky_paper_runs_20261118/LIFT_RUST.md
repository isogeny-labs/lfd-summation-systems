# Degree lift V_D -> V_{D+1} n M_<=D in Rust (`tools/hky_lift`)

Date: 2026-11-18. Worktree `wt/paper-lift` (sparse checkout of isogeny-labs/crypto), local commits only.
Labels: **proof** (argument written here or in LIFT.md), **exact** (computed, record named), **estimate**.
Every number is `[pcc:unverified]`: no pcc claim was filed; nothing here is a leaderboard or benchmark row.

## 0. Summary

- A standalone Rust crate `tools/hky_lift/` (own `Cargo.toml`/`Cargo.lock`, not in the root workspace,
  deps rayon, sha2, serde_json) implements the theorem of
  `lfd-summation-systems/experiments/ic-bounds-k14/LIFT.md` from its statement and proof, as a **third**
  implementation: in ordinary monomial coordinates with **no sigma split** (the theorem with the trivial
  symmetry group) and with V_D taken in its reduced echelon form in the input file's own column order
  (so the independent implementation's decomposition V_D = V_P (+) V_N has V_N = 0 here). It reads the
  engine's `hky.basis/v1` export directly.
- **Exact agreement everywhere it was compared** (section 4): 41 lift runs against direct closures computed
  by this crate (dim V_{D+1}, profile, 1-membership, canonical hash of V_{D+1} n M_<=D, sigma-part
  profiles), all equal to the recorded oracle values of the paper repository as well; 11 lifts on the
  **engine's own exports** (`pkm_tower_pilot --export-basis`, column order 1), where in addition the
  engine's V_D equals this crate's direct closure V_D (same canonical hash); the paper's unit instance
  `unit-p65537-a1-s10-k10`, D = 5 -> 6 (1 in V_6, dim 60,460, the full space) and the proper k = 10
  instance (60,446, 1 not in V_6) from the engine's export, field by field against the recorded records.
- **Negative controls** (section 5): every corrupted basis (raw byte flips, a dropped row, a perturbed
  coefficient; 29 files) is refused (trailer SHA-256, structural checks, random closedness check or the
  generator check); with the refusal overridden, a perturbed proper-instance basis makes the lift report
  1 in V_5 where the true answer is 1 not in V_5 (so the input checks are load-bearing); the test hook
  "no fall shifts" changes the answer (3473 -> 2844, 1 lost; the same 2844 the Python lift's hook gave).
- **No disagreement** between the theorem, the two Python implementations and this one was found.
- k = 14 has **not** been run (it needs the engine's closed V_6(14) export, ~8 h of closure); the command,
  the expected values and a field-by-field verifier are in section 6; projected cost in section 7.

## 1. Build, binary, commands

```bash
export CARGO_TARGET_DIR=/home/user/crypto-wt/paper-lift/target
cd tools/hky_lift && cargo build --release && cargo test --release
# frozen binary (no target-cpu flag; AVX2 kernel chosen at run time; glibc >= 2.34):
#   research/hky_paper_runs_20261118/lift_results/bin/hky_lift   (git-ignored)   sha256 in section 8
hky_lift lift --basis V_D.hkybasis [--expect recorded.json] [--out out.json]   # the lift
hky_lift verify --got out.json --expected recorded.json                        # field-by-field
hky_lift closure --p P --a A --b B --X X --k K --D D [--export F --order 0|1 --scramble N]
hky_lift validate --p .. --k .. --D .. [--vd-basis engine.basis] [--parts]     # lift vs direct closure
hky_lift corrupt --in F --out G --mode drop-row:I|perturb:I|flip-byte:OFF      # negative controls
```
Exit codes: 0 ok; 2 input refused; 3 a field differs from the expectation.
`lift` options: `--include-low` (shift every row, makes the result independent of the closedness of V_D,
LIFT.md section 10), `--batch-tops`, `--work-gb`, `--cache-gb` (follower cache), `--fall-every`,
`--no-verify-falls`, `--no-hash`, `--parts` (sigma-part profiles; small cases only), `--trials`,
`--no-refuse` and the test hooks `--test-no-fall-shifts`, `--test-drop-shift S`.

## 2. Input format adopted

Spec: `/home/user/crypto-wt/paper-engine/research/hky_paper_runs_20261118/BASIS_FORMAT.md`, sha256
`05c5d4e591f06c419829083cd93a3fbc51e6e2602c6802ba33a705a2904819e6`, copied unchanged to
`research/hky_paper_runs_20261118/BASIS_FORMAT.adopted.md`. The reader (`basis.rs`, `read_hky_basis`) makes
every check the spec lists (magic, version, m = 2, mask_encoding, n_cols = binomial sum, the column list:
distinct masks, weight <= D, degree non-increasing, last column = 1, and exactly the named order 0 or 1;
pivots strictly increasing; entries strictly ascending and after the pivot; coefficients in [1, p); the
profile and the contains_one flag recomputed from the pivots; the trailing SHA-256), in one streaming pass,
keeping only the byte offset of each row (rows are re-read with `pread` when needed: the k = 14 file is
estimated at ~90 GB). The engine's sidecar `instance` block (a, b, target[0] = X) gives the generator; the
lift then checks that g and random combinations of the m g (weight m <= D - 2) lie in the stored space and
that random combinations of the low rows stay in it under all 2k shifts, and **refuses** otherwise.

The lift works in the file's own column order (both orders are degree non-increasing, which is all the
theorem needs). The echelon is **not** assumed reduced: `Rref::from_echelon` fully reduces it by blocked
back-substitution (rows from the last pivot up; a block of 2048 rows in parallel against the finished rows
below it, then sequentially inside the block). Memory: rank x #free u32 (section 7).

## 3. Theorem step -> code

All in `tools/hky_lift/src/lift.rs` unless named.

| LIFT.md | here |
|---|---|
| L = span(M_<=D), column order | `layout.rs` `Layout` (file's column list); shift u_{i,j}: `Mon::shift` (add 2^j to block i mod 2^k) |
| V_D closed, reduced echelon basis, rows e_x + B[x] | `basis.rs` `Rref::from_echelon` |
| "regular" weight-D monomial | a pivot of V_D of weight D (unsplit: no orbit condition) |
| special / regular weight-(D+1) monomial | `Lift::new`: Sp = { u n : n free of weight D, weight(u n) = D + 1 } = tops with a free divisor; regular tops = the rest |
| T = F^Sp (+) F^Fr, red | `Tgt` = Reg / T unit / Piv (-B[c] on the free block); `Lift::tgt` |
| Lemma 2 (only rows of pivot weight D) | `run`: `regx`; `--include-low` adds the other rows |
| Lemma 4, followers R_u, "no regular coordinate" | `follower_inner` (asserts no `Tgt::Reg` in a follower row); `rho_neg_n` gives red(-B[c]) mod K |
| Lemma 4/5 generator count: each regular top has exactly D + 1 generators, each divisor a pivot | `run`: combinatorial assertion over all regular tops and the count identity sum_x #{u : u x regular top} = (D+1) |Reg| (both asserted before any product) |
| Lemma 5 (differences gamma_i - gamma_1; direct generators) | `run`: grouped units per batch of regular tops, differences; direct units (u x special or carried to weight <= D) with eterm |
| K in reduced echelon form, Sp first | `KSpace` (rows full width on T, Sp columns first); every stored vector kept reduced modulo K on the live columns, synced when K grows (`Fol::sync`) |
| falls = rows of K with pivot in the Fr block; fixed point (Lemma 1, closure invariant) | `fall_round` / `fall_fixed_point`: each fall row shifted once at its state (all 2k shifts), until none pending |
| explicit check Phi R_u c K (independent impl.'s step (4)) | `verify_falls`: all final fall rows shifted again, residuals must be 0 (`verify_falls.nonzero` in the record) |
| statement 2: dim V_{D+1} = dim V_D + |Reg| + dim K | `finish`: `dim_VD1` |
| statement 3: 1 in V_{D+1} | `finish`: 1 is a pivot of V_D, or its Fr column is a pivot of K (last column, so its K row is the unit vector) |
| statement 1/4: V_{D+1} n L = V_D (+) F, pivots(V_D) u pivots(F) | `finish`: `profile` from the union of the pivots; canonical hash `low_basis_sha256` of the reduced echelon basis of V_{D+1} n L (V_D rows reduced by the fall rows, streamed in chunks) |

**Proof notes for the unsplit form (proof).** Lemma 3 is not needed. Lemma 4: for a pivot x of weight D,
l_x = e_x + B[x] with B[x] on free columns only (RREF), so u l_x has top coordinates u x (coefficient 1) and
u f for free f of weight D; u f of weight D + 1 means f is u f minus one bit, i.e. a free divisor, so u f is
special. The generators of a regular top t are u l_{t / u} for its D + 1 bits; a weight-(D+1) monomial is
reached from weight D only without carry, i.e. from a divisor. The fall invariant: a K row with pivot q in
the Fr block changes only by multiples of new rows with pivot after q, which are again in the Fr block, hence
new (pending) falls.

## 4. Validation (expected vs got)

All runs below with the frozen binary (section 8) on the 4-vCPU development host, shared with other jobs
(1-minute load 13 to 26 throughout: the seconds are **not** benchmarks). Records:
`research/hky_paper_runs_20261118/lift_results/` (`run_all.sh` regenerates every file; `tables.sh` this table).

**(a) Lift vs direct closure, both computed here** (`validate_small.jsonl`, 41 runs). V_D is computed by
this crate's direct closure (`closure.rs`, written from the definition of V_n with "gens all", certified
by an exhaustive sweep: every shift of every final low row reduces to zero), written as an `hky.basis/v1`
file **scrambled** to a non-reduced echelon (each row gets random multiples of 3 rows with later pivots),
read back by the lift; the oracle is the direct closure V_{D+1} in M_<=D+1. Compared: dim V_{D+1}, profile
(dim V_{D+1} n M_<=j, j <= D + 1), 1 in V_{D+1}, the canonical hash of the reduced echelon basis of
V_{D+1} n M_<=D (the direct closure's rows of pivot weight <= D: LIFT.md section 5 "the hash from the direct
closure"), and with `parts` the sigma-part profiles (rank of (1 +- sigma) applied to V_{D+1} n M_<=j).
"recorded oracle" compares dim, 1 and profile with `lfd-summation-systems/experiments/ic-bounds-k14-independent/validation.json`
(the paper's direct closures), and the part dimensions with its `dim_low_parts`.
Variants: `order 1` (the engine's sigma-grouped column order), `include-low`, `fall-every 0` (falls only
at the end), tiny batches with a follower cache of 100 kB (evictions and rebuilds), unscrambled input.

**(b) Lift on the engine's own exports** (`validate_engine.jsonl`, 11 runs; `pkm_tower_pilot
--instance unit-p786433-a1-s10-k{k}.json --D D --export-basis`, column order 1, engine binary sha256
`d31a012b…`). Same comparisons, plus: the engine's V_D equals this crate's direct closure V_D (same
canonical hash in column order 1) in every case.

| case | flags | dim V_D | direct V_{D+1}: dim, 1 in, profile | lift: dim, 1 in | profile, hash of V_{D+1} n M_<=D, sigma parts | recorded oracle (paper repo) | dim T / dim K / dim F / rounds | lift s / direct s | all equal |
|---|---|---|---|---|---|---|---|---|---|
| unit-65537-k3-D3to4 | parts  | 7 | 57, true, 1 7 22 42 57 | 57, true | true, true, true 23+19 | true (recorded parts 23+19) | 49 / 49 / 35 / 5 | 0.01 / 0 | **true** |
| unit-65537-k4-D3to4 | parts  | 9 | 163, true, 1 9 37 93 163 | 163, true | true, true, true 49+44 | true (recorded parts 49+44) | 153 / 153 / 84 / 7 | 0.01 / 0 | **true** |
| unit-65537-k5-D3to4 | parts  | 11 | 386, true, 1 11 56 176 386 | 386, true | true, true, true 91+85 | true (recorded parts 91+85) | 374 / 374 / 165 / 9 | 0.03 / 0.01 | **true** |
| unit-65537-k6-D4to5 | parts  | 605 | 1586, true, 1 13 79 299 794 1586 | 1586, true | true, true, true 408+386 | true (recorded parts 408+386) | 197 / 197 / 189 / 2 | 0.11 / 0.63 | **true** |
| unit-65537-k7-D4to5 | parts  | 829 | 3473, true, 1 15 106 470 1471 3473 | 3473, true | true, true, true 750+721 | true (recorded parts 750+721) | 1790 / 1790 / 642 / 4 | 1.64 / 2.56 | **true** |
| unit-65537-k8-D4to5 | parts  | 1061 | 6885, true, 1 17 137 697 2517 6885 | 6885, true | true, true, true 1277+1240 | true (recorded parts 1277+1240) | 4900 / 4900 / 1456 / 6 | 12.31 / 11.19 | **true** |
| unit-12289-k4-D3to4 |   | 9 | 163, true, 1 9 37 93 163 | 163, true | true, true | true (recorded parts 49+44) | 153 / 153 / 84 / 7 | 0.01 / 0 | **true** |
| unit-12289-k6-D4to5 |   | 605 | 1586, true, 1 13 79 299 794 1586 | 1586, true | true, true | true (recorded parts 408+386) | 197 / 197 / 189 / 2 | 0.06 / 0.32 | **true** |
| unit-12289-k7-D4to5 |   | 829 | 3473, true, 1 15 106 470 1471 3473 | 3473, true | true, true | true (recorded parts 750+721) | 1790 / 1790 / 642 / 4 | 1.47 / 2.17 | **true** |
| unit-163841-k5-D3to4 |   | 11 | 386, true, 1 11 56 176 386 | 386, true | true, true | true (recorded parts 91+85) | 374 / 374 / 165 / 9 | 0.02 / 0.01 | **true** |
| unit-163841-k6-D4to5 |   | 605 | 1586, true, 1 13 79 299 794 1586 | 1586, true | true, true | true (recorded parts 408+386) | 197 / 197 / 189 / 2 | 0.06 / 0.36 | **true** |
| unit-163841-k7-D4to5 |   | 829 | 3473, true, 1 15 106 470 1471 3473 | 3473, true | true, true | true (recorded parts 750+721) | 1790 / 1790 / 642 / 4 | 1.64 / 2.11 | **true** |
| unit-65537-k2-D2to3 | parts  | 1 | 5, false, 0 0 1 5 | 5, false | true, true, true 1+0 | true (recorded parts 1+0) | 14 / 4 / 0 / 0 | 0 / 0 | **true** |
| unit-65537-k3-D2to3 |   | 1 | 7, false, 0 0 1 7 | 7, false | true, true | true (recorded parts 1+0) | 41 / 6 / 0 / 0 | 0 / 0 | **true** |
| unit-65537-k4-D2to3 |   | 1 | 9, false, 0 0 1 9 | 9, false | true, true | true (recorded parts 1+0) | 92 / 8 / 0 / 0 | 0 / 0 | **true** |
| unit-65537-k6-D3to4 | parts  | 13 | 605, false, 0 0 2 111 605 | 605, false | true, true, true 58+53 | true (recorded parts 58+53) | 780 / 591 / 98 / 9 | 0.09 / 0.04 | **true** |
| unit-65537-k7-D3to4 |   | 15 | 829, false, 0 0 2 115 829 | 829, false | true, true | true (recorded parts 60+55) | 1455 / 813 / 100 / 9 | 0.48 / 0.14 | **true** |
| unit-65537-k8-D3to4 |   | 17 | 1061, false, 0 0 2 119 1061 | 1061, false | true, true | true (recorded parts 62+57) | 2499 / 1043 / 102 / 9 | 1.29 / 0.29 | **true** |
| unit-12289-k6-D3to4 |   | 13 | 605, false, 0 0 2 111 605 | 605, false | true, true | true (recorded parts 58+53) | 780 / 591 / 98 / 9 | 0.11 / 0.04 | **true** |
| unit-163841-k7-D3to4 |   | 15 | 829, false, 0 0 2 115 829 | 829, false | true, true | true (recorded parts 60+55) | 1455 / 813 / 100 / 9 | 0.34 / 0.09 | **true** |
| proper-65537-X2-k6-D3to4 | parts  | 13 | 605, false, 0 0 2 111 605 | 605, false | true, true, true 58+53 | true (recorded parts 58+53) | 780 / 591 / 98 / 9 | 0.09 / 0.04 | **true** |
| proper-65537-X2-k6-D4to5 | parts  | 605 | 1584, false, 0 11 77 297 792 1584 | 1584, false | true, true, true 407+385 | true (recorded parts 407+385) | 197 / 195 / 187 / 2 | 0.12 / 0.7 | **true** |
| proper-65537-X2-k7-D3to4 |   | 15 | 829, false, 0 0 2 115 829 | 829, false | true, true | true (recorded parts 60+55) | 1455 / 813 / 100 / 9 | 0.41 / 0.11 | **true** |
| proper-65537-X2-k7-D4to5 | parts  | 829 | 3471, false, 0 13 104 468 1469 3471 | 3471, false | true, true, true 749+720 | true (recorded parts 749+720) | 1790 / 1788 / 640 / 4 | 1.11 / 2.08 | **true** |
| proper-65537-X2-k8-D3to4 |   | 17 | 1061, false, 0 0 2 119 1061 | 1061, false | true, true | true (recorded parts 62+57) | 2499 / 1043 / 102 / 9 | 1.29 / 0.3 | **true** |
| proper-65537-X2-k8-D4to5 | parts  | 1061 | 6881, false, 0 13 133 693 2513 6881 | 6881, false | true, true, true 1275+1238 | true (recorded parts 1275+1238) | 4900 / 4896 / 1452 / 6 | 11.83 / 12.78 | **true** |
| proper-15361-k6-D4to5 |   | 605 | 1584, false, 0 11 77 297 792 1584 | 1584, false | true, true | true (recorded parts 407+385) | 197 / 195 / 187 / 2 | 0.06 / 0.79 | **true** |
| proper-15361-k7-D4to5 |   | 829 | 3469, false, 0 11 102 466 1467 3469 | 3469, false | true, true | true (recorded parts 748+719) | 1790 / 1786 / 638 / 4 | 0.95 / 1.76 | **true** |
| proper-15361-k7-D3to4 |   | 15 | 829, false, 0 0 2 115 829 | 829, false | true, true | true (recorded parts 60+55) | 1455 / 813 / 100 / 9 | 0.36 / 0.1 | **true** |
| order1-unit-65537-k7-D4to5 | order 1 parts  | 829 | 3473, true, 1 15 106 470 1471 3473 | 3473, true | true, true, true 750+721 | true | 1790 / 1790 / 642 / 4 | 1.58 / 2.08 | **true** |
| order1-proper-65537-X2-k7-D4to5 | order 1  | 829 | 3471, false, 0 13 104 468 1469 3471 | 3471, false | true, true | true | 1790 / 1788 / 640 / 4 | 1.7 / 2.29 | **true** |
| order1-unit-65537-k7-D3to4 | order 1  | 15 | 829, false, 0 0 2 115 829 | 829, false | true, true | true | 1455 / 813 / 100 / 9 | 0.48 / 0.16 | **true** |
| low-unit-65537-k7-D4to5 | include-low  | 829 | 3473, true, 1 15 106 470 1471 3473 | 3473, true | true, true | true | 1790 / 1790 / 642 / 4 | 1.53 / 1.97 | **true** |
| low-proper-65537-X2-k8-D4to5 | include-low  | 1061 | 6881, false, 0 13 133 693 2513 6881 | 6881, false | true, true | true | 4900 / 4896 / 1452 / 6 | 10.33 / 13.84 | **true** |
| endfalls-unit-65537-k8-D3to4 | fall-every 0  | 17 | 1061, false, 0 0 2 119 1061 | 1061, false | true, true | true | 2499 / 1043 / 102 / 9 | 1.11 / 0.39 | **true** |
| endfalls-unit-65537-k7-D4to5 | fall-every 0  | 829 | 3473, true, 1 15 106 470 1471 3473 | 3473, true | true, true | true | 1790 / 1790 / 642 / 4 | 1.56 / 2.13 | **true** |
| tiny-unit-65537-k7-D4to5 | batch-tops 7 cache-gb 0.0001  | 829 | 3473, true, 1 15 106 470 1471 3473 | 3473, true | true, true | true | 1790 / 1790 / 642 / 5 | 1.05 / 2.25 | **true** |
| tiny-proper-15361-k7-D4to5 | batch-tops 5 cache-gb 0.0001  | 829 | 3469, false, 0 11 102 466 1467 3469 | 3469, false | true, true | true | 1790 / 1786 / 638 / 5 | 3.02 / 2.29 | **true** |
| noscramble-unit-65537-k6-D4to5 | scramble 0  | 605 | 1586, true, 1 13 79 299 794 1586 | 1586, true | true, true | true | 197 / 197 / 189 / 2 | 0.09 / 0.39 | **true** |
| unit-65537-k9-D4to5 |   | 1301 | 12616, true, 1 19 172 988 4048 12616 | 12616, true | true, true | true (recorded parts 2047+2001) | 10321 / 10321 / 2747 / 8 | 60.09 / 65.56 | **true** |
| proper-65537-X1424-k9-D4to5 |   | 1301 | 12614, false, 0 17 170 986 4046 12614 | 12614, false | true, true | true (recorded parts 2046+2000) | 10321 / 10319 / 2745 / 8 | 60.84 / 55.48 | **true** |
| engine-unit-p786433-k4-D3to4 |  engine export (order 1); engine V_D = our V_D: true | 9 | 163, true, 1 9 37 93 163 | 163, true | true, true, true 49+44 | - | 153 / 153 / 84 / 7 | 0.01 / 0.01 | **true** |
| engine-unit-p786433-k5-D3to4 |  engine export (order 1); engine V_D = our V_D: true | 11 | 386, true, 1 11 56 176 386 | 386, true | true, true, true 91+85 | - | 374 / 374 / 165 / 9 | 0.03 / 0.03 | **true** |
| engine-unit-p786433-k6-D3to4 |  engine export (order 1); engine V_D = our V_D: true | 13 | 605, false, 0 0 2 111 605 | 605, false | true, true, true 58+53 | - | 780 / 591 / 98 / 9 | 0.1 / 0.04 | **true** |
| engine-unit-p786433-k6-D4to5 |  engine export (order 1); engine V_D = our V_D: true | 605 | 1586, true, 1 13 79 299 794 1586 | 1586, true | true, true, true 408+386 | - | 197 / 197 / 189 / 2 | 0.18 / 0.46 | **true** |
| engine-unit-p786433-k6-D5to6 |  engine export (order 1); engine V_D = our V_D: true | 1586 | 2510, true, 1 13 79 299 794 1586 2510 | 2510, true | true, true, true 804+782 | - | 0 / 0 / 0 / 0 | 0.1 / 1.22 | **true** |
| engine-unit-p786433-k7-D3to4 |  engine export (order 1); engine V_D = our V_D: true | 15 | 829, false, 0 0 2 115 829 | 829, false | true, true, true 60+55 | - | 1455 / 813 / 100 / 9 | 0.38 / 0.11 | **true** |
| engine-unit-p786433-k7-D4to5 |  engine export (order 1); engine V_D = our V_D: true | 829 | 3473, true, 1 15 106 470 1471 3473 | 3473, true | true, true, true 750+721 | - | 1790 / 1790 / 642 / 4 | 1.54 / 1.97 | **true** |
| engine-unit-p786433-k7-D5to6 |  engine export (order 1); engine V_D = our V_D: true | 3473 | 6476, true, 1 15 106 470 1471 3473 6476 | 6476, true | true, true, true 1751+1722 | - | 0 / 0 / 0 / 0 | 0.22 / 7.4 | **true** |
| engine-unit-p786433-k8-D3to4 |  engine export (order 1); engine V_D = our V_D: true | 17 | 1061, false, 0 0 2 119 1061 | 1061, false | true, true, true 62+57 | - | 2499 / 1043 / 102 / 9 | 0.92 / 0.25 | **true** |
| engine-unit-p786433-k8-D4to5 |  engine export (order 1); engine V_D = our V_D: true | 1061 | 6885, true, 1 17 137 697 2517 6885 | 6885, true | true, true, true 1277+1240 | - | 4900 / 4900 / 1456 / 6 | 7.31 / 9.08 | **true** |
| engine-unit-p786433-k9-D4to5 |  engine export (order 1); engine V_D = our V_D: true | 1301 | 12616, true, 1 19 172 988 4048 12616 | 12616, true | true, true, true 2047+2001 | - | 10321 / 10321 / 2747 / 8 | 49.47 / 57.61 | **true** |


In the "recorded oracle" column, true = dim, 1 and profile equal to validation.json; the recorded part
dimensions are printed for reference and agree with the sigma-part sums of the previous column wherever
those were computed (`parts`); 52 of 52 rows `all equal`. Engine rows have no recorded oracle (p = 786433
instances are new) and are compared with the direct closure only.

**(c) The paper's k = 10 cells, D = 5 -> 6, from the engine's export, field by field** (`lift-*-D5to6.json`;
`hky_lift lift --basis … --expect …`). Expected values: unit = the recorded on-card lift
`experiments/ic-bounds-k14/runs/lift-unit-p65537-k10-D5to6-gpu.json`; proper = LIFT.md section 5
(2026-11-11 addition) / validation.json oracle (direct closure), in `lift_results/expected-proper-p65537-s12-k10-D5to6.json`.

| field | unit-p65537-a1-s10-k10: expected | got | proper-p65537-s12-k10 (X = 1424): expected | got |
|---|---|---|---|---|
| dim V_5 | 18,569 | 18,569 | 18,569 | 18,569 |
| 1 in V_6 | true | **true** | false | **false** |
| dim V_6 | 60,460 | 60,460 | 60,446 | 60,446 |
| profile | 1 21 211 1351 6196 21700 60460 | equal | 0 7 197 1337 6182 21686 60446 | equal |
| dim F = dim(V_6 n M_<=5) - dim V_5 | 3,131 (sum of the recorded per-part falls) | 3,131 | 3,117 | 3,117 |
| fixed point | true | true | true | true |
| verify_falls (all final falls re-shifted, nonzero residuals) | - | K = T (skipped) | - | 62,340 rows, 0 nonzero |
| input checks (g, random m g; 60 shifted random low combinations) | - | pass | - | pass |
| this crate's T: dim T / dim K / special tops / regular tops | (split: [1965,1954] / [1965,1954]) | 4,199 / 4,199 / 1,068 / 37,692 | | 4,199 / 4,185 / 1,068 / 37,692 |
| madds (products + syncs + K) / follower syncs | 6.50e11 (GPU record, split) | 4.46e11 / 2.71e11 | | 4.54e11 / 2.72e11 |
| wall (shared 4 vCPU), peak RSS | 8.4 s on the H200 | 109 s, 1.8 GB | | 109 s, 1.8 GB |

`dim T`, `dim K`, the special/regular counts are **not** comparable with the split implementations (they
count orbit columns per part; here monomials, and "special" is relative to the monomial echelon form).
`falls` is comparable as dim F.


## 5. Negative controls

`lift_results/negative_controls.jsonl` (`negative_controls.sh`): V_4 of k = 7 (unit: 1 enters at D = 5;
proper X = 2: 1 never enters), exported scrambled, then corrupted with `hky_lift corrupt`. `drop-row:I` and
`perturb:I` rewrite header, profile and trailer, so the file is structurally valid and only the space
changed; `flip-byte` flips one raw bit and leaves the trailer. Every corrupted file is **refused** (22 of 22,
plus 7 further raw flips in `negative_sha_flips.jsonl`: 5 caught by the trailing SHA-256 or truncation,
2 by structural checks first). With the refusal overridden (`--no-refuse`, for the record only):

* a dropped row (input a proper subspace of V_D; refused because it is not closed, or because g or some
  m g is no longer in it) gives the same V_{D+1} (dimension and 1) in all 10 runs: the iteration recovered
  the missing row here, which the theorem does not promise (its hypothesis is a closed V_D), so the refusal
  stays;
* a perturbed coefficient (input NOT inside V_D) makes the **proper** instance report 1 in V_5, dim 3,473,
  where the true answer is 1 not in V_5, dim 3,471: the lift is only as sound as its input, and the input
  checks are what catch this. (The checks are random, 3 trials x 2k shifts plus g and 3 random
  combinations of the m g; they caught all 20 changed spaces here.)

Test hooks (input clean): `--test-no-fall-shifts` changes the answer for both instances (dim 2,844 instead
of 3,473 / 3,471, and 1 lost for the unit instance) — the same 2,844 the Python lift's hook reported for
unit k = 7 (LIFT.md section 5); `--test-drop-shift 0` and `13` leave it unchanged (the dropped generators
are redundant here, as in 14 of the Python lift's 16 runs).

| control | exit | outcome | dim V_D read | dim V_{D+1} | 1 in V_{D+1} | random closedness: shifted combos leaving | generator check: not in |
|---|---|---|---|---|---|---|---|
| unit clean | 0 | done | 829 | 3473 | true | 0 | 0 |
| unit clean, include-low | 0 | done | 829 | 3473 | true | 0 | 0 |
| unit flip-byte:2000 (raw) | 2 | refused: mask 0x20a9 has weight 5 > D | - | - | - | - | - |
| unit drop-row:0 | 2 | refused: the stored space is not closed: shifted random low elements leave it | - | - | - | 3 | 0 |
| unit drop-row:0 --no-refuse | 0 | done | 828 | 3473 | true | 3 | 0 |
| unit perturb:0 | 2 | refused: the stored space is not closed: shifted random low elements leave it | - | - | - | 3 | 0 |
| unit perturb:0 --no-refuse | 0 | done | 829 | 3473 | true | 3 | 0 |
| unit drop-row:100 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 12 | 3 |
| unit drop-row:100 --no-refuse | 0 | done | 828 | 3473 | true | 12 | 3 |
| unit perturb:100 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 12 | 3 |
| unit perturb:100 --no-refuse | 0 | done | 829 | 3473 | true | 12 | 3 |
| unit drop-row:400 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 36 | 3 |
| unit drop-row:400 --no-refuse | 0 | done | 828 | 3473 | true | 36 | 3 |
| unit perturb:400 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 36 | 3 |
| unit perturb:400 --no-refuse | 0 | done | 829 | 3473 | true | 36 | 3 |
| unit drop-row:700 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 42 | 3 |
| unit drop-row:700 --no-refuse | 0 | done | 828 | 3473 | true | 42 | 3 |
| unit perturb:700 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 42 | 3 |
| unit perturb:700 --no-refuse | 0 | done | 829 | 3473 | true | 42 | 3 |
| unit drop-row:828 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 42 | 4 |
| unit drop-row:828 --no-refuse | 0 | done | 828 | 3473 | true | 42 | 4 |
| unit perturb:828 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 42 | 4 |
| unit perturb:828 --no-refuse | 0 | done | 829 | 3473 | true | 42 | 4 |
| unit hook: no fall shifts | 0 | done | 829 | 2844 | - | 0 | 0 |
| unit hook: drop shift 0 | 0 | done | 829 | 3473 | true | 0 | 0 |
| unit hook: drop shift 13 | 0 | done | 829 | 3473 | true | 0 | 0 |
| proper clean | 0 | done | 829 | 3471 | - | 0 | 0 |
| proper clean, include-low | 0 | done | 829 | 3471 | - | 0 | 0 |
| proper flip-byte:2000 (raw) | 2 | refused: mask 0x20a9 has weight 5 > D | - | - | - | - | - |
| proper drop-row:0 | 2 | refused: the stored space is not closed: shifted random low elements leave it | - | - | - | 3 | 0 |
| proper drop-row:0 --no-refuse | 0 | done | 828 | 3471 | - | 3 | 0 |
| proper perturb:0 | 2 | refused: the stored space is not closed: shifted random low elements leave it | - | - | - | 3 | 0 |
| proper perturb:0 --no-refuse | 0 | done | 829 | 3473 | true | 3 | 0 |
| proper drop-row:100 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 12 | 3 |
| proper drop-row:100 --no-refuse | 0 | done | 828 | 3471 | - | 12 | 3 |
| proper perturb:100 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 12 | 3 |
| proper perturb:100 --no-refuse | 0 | done | 829 | 3473 | true | 12 | 3 |
| proper drop-row:400 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 36 | 3 |
| proper drop-row:400 --no-refuse | 0 | done | 828 | 3471 | - | 36 | 3 |
| proper perturb:400 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 36 | 3 |
| proper perturb:400 --no-refuse | 0 | done | 829 | 3473 | true | 36 | 3 |
| proper drop-row:700 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 42 | 3 |
| proper drop-row:700 --no-refuse | 0 | done | 828 | 3471 | - | 42 | 3 |
| proper perturb:700 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 42 | 3 |
| proper perturb:700 --no-refuse | 0 | done | 829 | 3473 | true | 42 | 3 |
| proper drop-row:828 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 42 | 4 |
| proper drop-row:828 --no-refuse | 0 | done | 828 | 3471 | - | 42 | 4 |
| proper perturb:828 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 42 | 4 |
| proper perturb:828 --no-refuse | 0 | done | 829 | 3473 | true | 42 | 4 |
| proper hook: no fall shifts | 0 | done | 829 | 2844 | - | 0 | 0 |
| proper hook: drop shift 0 | 0 | done | 829 | 3471 | - | 0 | 0 |
| proper hook: drop shift 13 | 0 | done | 829 | 3471 | - | 0 | 0 |



## 6. k = 14 reproduction (pod; not run here)

Input: the engine's `hky.basis/v1` export of the **closed** V_6(14) with its sidecar (the sidecar's
`instance` block names a, b, X). The engine run is the other worker's (`pkm_tower_pilot --instance
<lfd-summation-systems>/experiments/ic-bounds-k14/instances/<inst>.json --D 6 --export-basis
/workspace/vd/<inst>-D6.basis --out …`, about 8 h at k = 14; the file is estimated at ~90 GB: the k = 10
export has 1,486 stored entries per row on average, 220.9 MB for 18,569 rows).

```bash
# on the pod, with the frozen binary (section 8) and the recorded lift JSON of the paper repository
tools/hky_lift/scripts/pod_k14.sh lift_results/bin/hky_lift /workspace/vd/unit-p65537-a1-s10-k14-D6.basis \
    experiments/ic-bounds-k14/runs/lift-unit-p65537-k14-D6to7.json /workspace/out
#   CACHE_GB (follower cache, default 120) and BATCH_TOPS (default 8192) are environment knobs;
#   exit 0 = every recorded field equal, 3 = a field differs, 2 = the input was refused.
# the comparison alone, on any host:
hky_lift verify --got /workspace/out/lift-unit-p65537-k14-D6to7.rust.json \
    --expected experiments/ic-bounds-k14/runs/lift-unit-p65537-k14-D6to7.json
```

Expected values (from the recorded lift JSON; `verify` compares them field by field; `falls` and
dim(V_7 n M_<=6) are the sums of the recorded per-part values):

| recorded lift | p | instance | dim V_6 | 1 in V_7 | dim V_7 | profile of V_7 | dim F | dim(V_7 n M_<=6) |
|---|---|---|---|---|---|---|---|---|
| `lift-unit-p65537-k14-D6to7.json` | 65537 | unit-p65537-a1-s10-k14 | 448,922 | true | 1,683,218 (full space) | 1 29 407 3683 24158 122438 499178 1683218 | 50,256 | 499,178 |
| `lift-unit-p163841-k14-D6to7.json` | 163841 | unit-p163841-a3-s21-k14 | 448,922 | true | 1,683,218 (full space) | same | 50,256 | 499,178 |
| `lift-proper-p114689-k14-D6to7.json` | 114689 | proper-p114689-s22-k14 (2,370 solutions) | 448,922 | false | 1,680,848 (= the ideal bound) | 0 0 2 1313 21788 120068 496808 1680848 | 47,886 | 496,808 |

Also `fixed_point: true`. Not comparable (split-only quantities): `dim_T`, `dim_K`, per-part hashes,
`n_single`; the part profiles are comparable only with `--parts`, which is for small cases.
The V_6 input itself should match the recorded closure: rank 448,922, profile
0 0 2 154 4451 73587 448922 (`runs/unit-p65537-k14-D6.json`), which the reader re-derives from the pivots.


## 7. Memory and time

**Memory (estimate; unsplit, u32 storage).** n = #columns of M_<=6, r = rank V_6, nf = n - r free columns,
nT = |Sp| + nf (|Sp| = special tops; at k = 14 the split records give 5,519 special orbits per part, about
11,000 monomials, and the bound |Sp| <= 22 x #free weight-6 monomials = 30,910).

| | k = 14 | k = 15 | k = 16 |
|---|---|---|---|
| n (V_6 columns) | 499,178 | 768,212 | 1,149,017 |
| r (k = 14 recorded; k = 15, 16 at ~90 %) | 448,922 | ~691,000 | ~1,034,000 |
| nf | 50,256 | ~76,800 | ~115,000 |
| B = reduced echelon on the free columns, r x nf x 4 B | **90 GB** | ~212 GB | ~476 GB |
| K (full width, <= nT^2 x 4 B; nT ~ 1.22 nf) | <= 15 GB | <= 35 GB | <= 79 GB |
| one follower R_u (nf x live T columns x 4 B, shrinks as K grows) | <= 12.3 GB | <= 29 GB | <= 64 GB |
| all 2k followers (cache them all: no rebuilds) | 345 GB | 870 GB | 2 TB |
| batch work (`--work-gb`, default 2) | 2 GB | 2 GB | 2 GB |
| input file on disk (pread; only row offsets in RAM) | ~90 GB | ~210 GB | ~470 GB |
| **host RAM: B + K + cache of ~10 followers + work** | **~230 GB** | ~540 GB | ~1.2 TB |

The minimum is B + K + one follower + work (k = 14: ~120 GB; followers are then rebuilt on every use).
The sigma split of the Python implementations would halve B and quarter K and the followers (orbit
coordinates, each part of half the size); it is **not** implemented here (it was left out on purpose, to
keep this a separate implementation of the theorem), so this port needs about twice the Python lift's host
memory (44 GB RSS recorded at k = 14 with the basis on the host).

**Time (estimate).** Measured here: k = 10, D = 5 -> 6: 4.46e11 (unit) and 4.54e11 (proper)
multiply-adds, of which 2.7e11 are follower syncs; the Python split lift recorded 6.50e11 for the unit cell.
At k = 14 the Python lift recorded 7.15e14 (unit p = 65537), 7.30e14 (unit p = 163841) and 1.33e15
(proper p = 114689) multiply-adds. Two scalings bracket this port: the k = 10 ratio (0.69) and the
unsplit-vs-split size factor (each product has 2 x the rows, columns and inner dimension of a part, two
parts: 4 x). Hence **unit: 5e14 to 3e15, proper: 9e14 to 5e15 multiply-adds**, plus the
back-substitution of the input echelon (rows x stored entries on pivot columns x nf; 2-3 s at k = 10;
<= ~5e14 at k = 14 if the stored entries per row scale with nf).
Kernel rate: 2.2 to 2.9 Gmadd/s per thread for the AVX2 u64 tile kernel (`hky_lift bench`, measured on
this host under a load of ~20 on 4 vCPU, so a lower bound). At ~2.5 Gmadd/s per core and 60 % parallel
efficiency on 64 cores (~100 Gmadd/s): **k = 14 unit ~1.5 to 8 h, proper ~2.5 to 14 h** on CPU (the Python
lift took 751 s and 4,679 s on an H200). k = 15 and 16: the products scale like rows x nf x nT, about
(1.54)^3 = 3.6 x and (2.3)^3 = 12 x the k = 14 cost if the regime (few special tops, K filling early)
persists — **estimate only**; nobody has measured V_6 at k = 15, 16.


## 8. Provenance

| item | value |
|---|---|
| source commit of the frozen binary | `ab6ba13436` (branch `wt/paper-lift`; `tools/hky_lift` unchanged since) |
| frozen binary | `research/hky_paper_runs_20261118/lift_results/bin/hky_lift` (git-ignored), sha256 `6824a8bdd8dd326b18c3e52fa1431c80dcf23582f9c52b3d6f375de51b3b1eaa`; rustc 1.97.0, `cargo build --release`, opt-level 3, no target-cpu flag (AVX2 tile kernel selected at run time), dynamically linked against glibc (needs GLIBC_2.34) |
| input spec | `BASIS_FORMAT.adopted.md` = paper-engine `research/hky_paper_runs_20261118/BASIS_FORMAT.md`, sha256 `05c5d4e5…19e6` |
| engine binary used for the exports | paper-engine `target/release/examples/pkm_tower_pilot` (copied to scratch), sha256 `d31a012b5e5d938b1b7ccb80b227ea6cf4322375b6a1170f1a587482204137bc` (= the `pilot_gate` sha recorded in the engine's gate log) |
| records | `lift_results/`: `validate_small.jsonl`, `validate_engine.jsonl`, `negative_controls.jsonl`, `negative_sha_flips.jsonl`, `lift-*-D5to6.json`, the engine sidecars of the k = 10 exports, `binaries.sha256`, `run_all.sh`, `tables.sh` |
| tests | `cargo test --release` (modular gemm against i128, both kernels via `HKY_LIFT_MICRO_KERNEL=1`; shift and column orders) |
| not done | k = 14, 15, 16 runs (need the engine's V_6 exports and a large-memory host); the sigma split; a GPU path; a resume/checkpoint layer for the lift (a k = 14 lift is restarted from the input file) |

**Limits.** The result is exact only if the input is the closed V_D: the lift checks g, random m g and
random shifted low combinations (refusing on any failure) but cannot prove closedness; `--include-low`
removes the dependence on closedness (LIFT.md section 10). No independent certificate of V_{D+1} is produced;
correctness rests on the theorem (whose unsplit form is argued in section 3), the in-run assertions
(no regular coordinate in a follower row, every divisor of a regular top a pivot, D + 1 generators each,
`verify_falls`), and the comparisons above. The lift is slower than the direct closure where almost every
top is special (tiny V_D; e.g. k = 9, D = 4 -> 5), as for the Python lift.

