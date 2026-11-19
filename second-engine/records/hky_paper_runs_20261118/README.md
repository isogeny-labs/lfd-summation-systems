# Paper runs on the HKY engine: exact-instance loader, basis export, checkpoint/resume, depth-16 instance

Date 2026-11-18. Worktree `wt/paper-engine` (sparse checkout of isogeny-labs/crypto), local commits only.
Purpose: a second, independent engine for the degree-6 closures of
`lfd-summation-systems/paper/curve-lfd-bounds.tex` (Limitations: "no second engine has computed the
degree-6 closure at k = 14"; "no exact value is known beyond 14"), run on pods with an exported,
re-checkable basis.

Every number here is `[pcc:unverified]` (no pcc claim filed); labels: **measured** (this machine,
this directory's scripts), **given** (task card or earlier records), **estimate**.
Not available here: Conductor, research-kb, ecdb, AWS/isolab, PMU. No pod was started.

## 1. The system correspondence (derived, then checked coefficient by coefficient)

**Paper** (`curve-lfd-bounds.tex`, normalized two-summand chain presentation; `csd/ic_dreg.py`
`pkm_generator`, `s3_prime_coeffs`; `csd/ic_bounds_d5.py`). Ring
R = F_p[X_1, X_2]/(X_1^{2^k} − 1, X_2^{2^k} − 1), presented by chain variables u_{i,j} = X_i^{2^j}
(j = 0..k−1) with u_{i,j}² = u_{i,j+1}, u_{i,k−1}² = 1; a square-free monomial is a mask with bit
(i−1)·k + j for u_{i,j}, i.e. X_1^{e_1}X_2^{e_2} ↦ e_1 | e_2 << k (multiplication by u_{i,j} adds
2^j to e_i mod 2^k). The factor base is the subgroup itself (`subgroup presentation, c1 = c2 = 1`):
x_i = X_i. The single generator is g = NF(S_3(x_1, x_2, X_R)), X_R the target abscissa, chain degree
2. V_D is the least subspace of M_≤D containing NF(μ g) (wt μ + 2 ≤ D) and closed under h ↦ NF(u·h)
for wt h < D; profile[w] = dim(V_D ∩ M_≤w), w = 0..D.

**Engine** (`hky_fp_tower.rs`, `V_{F,c}`; pilot Kummer tower). Ring
F_p[y]/(y_j² − y_{j+1}, y_{t−1}² − c) per block, x = scale · y_0, generator NF(S_3(x_1, x_2, x_R))
by the tower's square rules. `V_{F,c}` = least subspace of R_≤c containing F ∩ R_≤c and closed
under x·h for deg h < c (module doc). `dims[d]` = dim(V ∩ R_≤d) from pivot degrees (columns sorted by
degree, highest first). When 1 becomes a pivot the engine stops and reports V = M_≤c, which is exact:
1 has degree 0 < c, and multiplying by distinct chain variables one at a time (no carry) reaches every
square-free monomial of degree ≤ c.

**Correspondence.** `--instance` builds the Kummer tower with **scale = 1** and **top constant
c = 1**, t = k, m = 2. Then y_{i,j} ≡ u_{i,j} with engine variable index i·k + j = the paper's mask
bit, x_i = y_{i,0} = X_i, and the engine's generator is exactly g: **no rescaling at all.** (The
"diagonal rescaling" of `hky_engine_speed_20261117/README.md` concerns the pilot's random cosets
x = g·y, scale ≠ 1; for the paper's subgroup instances scale is 1.) The target enters only as the
constant third argument of S_3. Definitions of V_D, profile and contains_one coincide (the span of
{NF(μ g)} is the span of iterated variable products of g).

**Checks that run on every `--instance` invocation** (`examples/pkm_tower_paper/paper.rs`):
(i) the pilot's NF(S_3) equals a direct transcription of `s3_prime_coeffs` + `pkm_generator`
(mask ↦ coefficient, all 9 terms); (ii) it equals the instance file's `generator` field;
(iii) the tower is checked on all 2^k subgroup points; (iv) the target lies on the curve, p is prime.
The run aborts if any fails (except (ii) under an explicit `--override-*` control).

## 2. Build (harness) and binaries

The root crate does not compile; `harness/` compiles the needed modules by path (approach of
`research/rand_vs_struct_fb_20261118/harness/`), plus `sha2`.

```bash
export CARGO_TARGET_DIR=/home/user/crypto-wt/paper-engine/target
cd research/hky_paper_runs_20261118/harness
HKY_ENGINE_COMMIT=$(git rev-parse HEAD) cargo build --release --offline \
    --example pkm_tower_pilot --example hky_basis_check --example hky_paper_instance
```

No `target-cpu`/feature flags (no `.cargo/config`, empty `RUSTFLAGS`); SIMD is chosen at run time.
`ldd`: libgcc_s, libm, libc only; highest GLIBC symbol version 2.34 (≤ 2.35). Release target dir
≈ 200 MB. Frozen binaries are in `bin/` (git-ignored); hashes in §7.

Files:
- `src/cryptanalysis/hky_closure_io.rs` — checkpoint/resume (`ClosureState`), stop on SIGTERM,
  `hky.basis/v1` writer. Hooked into `closure_fast2` at three delimited points (`// ── io hook n of 3`).
- `examples/pkm_tower_paper/paper.rs` — the `--instance` mode of `pkm_tower_pilot`.
- `tools/hky_basis_check.rs` — independent reader/checker of `hky.basis/v1`.
- `tools/hky_paper_instance.rs` — the paper's unit-instance recipe in Rust (§6).
- `gate/` — gate, negative controls, checkpoint test, unit-ness test (scripts and their logs/records).
- `pod/` — job scripts (§8), `PREDICTIONS.md` — predicted wall/memory/disk, written before any pod run.
- `BASIS_FORMAT.md` — the export format.

## 3. Gate: recorded closures reproduced exactly

<!-- GATE:begin -->
Expected = `lfd-summation-systems/experiments/ic-bounds-k14/runs/<cell>.json` (`columns`, `profile` = dim(V_D ∩ M_≤w) for
w = 0..D, `contains_one`); got = this engine's `--instance` record; `match` requires all three equal and the closure
closed. Binaries (§8): `d31a012b` = `bin/pilot_gate` (pre-merge engine); `52b47153` = first provisional build (same
sources, commit-stamped); `0dfbffad` = the merged engine (wt/engine-opt c51685b6a9 merged at 9fb60e9c13); **`39acdeac` =
the final pod binary** (9fb60e9c13 + non-fatal checkpoint/export failures, stamped 450d5b295b). Muladds are equal across all
binaries on every shared cell (the optimizer changes no counter). Scripts: `gate/run_gate.sh`, `gate/final_gate.sh`,
`gate/table.sh`; records under `gate/raw*/`. Walls: 4 shared vCPU; the pre-merge rows ran at load 15–27, the merged/final
rows at load ≈ 4.

| cell (instance : D) | columns | expected profile | got profile | 1 ∈ V_D exp / got | match | muladds | wall (s, shared 4 vCPU) | binary |
|:--|--:|:--|:--|:-:|:-:|--:|--:|:--|
| unit-p65537-a1-s10-k10 : 5 | 21700 | 0,0,2,138,3163,18569 | 0,0,2,138,3163,18569 | false / false | **true** | 16024605579 | 43 | d31a012b |
| unit-p65537-a1-s10-k11 : 5 | 35443 | 0,0,2,142,3441,25033 | 0,0,2,142,3441,25033 | false / false | **true** | 16340811395 | 40 | d31a012b |
| unit-p65537-a1-s10-k12 : 5 | 55455 | 0,0,2,146,3727,32057 | 0,0,2,146,3727,32057 | false / false | **true** | 17441038634 | 45 | d31a012b |
| unit-p65537-a1-s10-k10 : 6 | 60460 | 1,21,211,1351,6196,21700,60460 | 1,21,211,1351,6196,21700,60460 | true / true | **true** | 43706615635 | 106 | d31a012b |
| unit-p65537-a1-s10-k13 : 5 | 83682 | 0,0,2,150,4021,39657 | 0,0,2,150,4021,39657 | false / false | **true** | 17993100441 | 55 | d31a012b |
| unit-p163841-a3-s21-k13 : 5 | 83682 | 0,0,2,150,4021,39657 | 0,0,2,150,4021,39657 | false / false | **true** | 17993415519 | 56 | d31a012b |
| unit-p65537-a1-s10-k14 : 5 | 122438 | 0,0,2,154,4323,47849 | 0,0,2,154,4323,47849 | false / false | **true** | 18299157464 | 61 | d31a012b |
| unit-p163841-a3-s21-k14 : 5 | 122438 | 0,0,2,154,4323,47849 | 0,0,2,154,4323,47849 | false / false | **true** | 18299505466 | 53 | d31a012b |
| unit-p65537-a1-s10-k11 : 6 | 110056 | 1,23,254,1794,9109,35443,110056 | 1,23,254,1794,9109,35443,110056 | true / true | **true** | 321522714727 | 1188 | d31a012b |
| unit-p65537-a1-s10-k15 : 5 | 174437 | 0,0,2,158,4633,56649 | 0,0,2,158,4633,56649 | false / false | **true** | 19233087929 | 72 | d31a012b |
| unit-p65537-a1-s10-k10 : 5 | 21700 | 0,0,2,138,3163,18569 | 0,0,2,138,3163,18569 | false / false | **true** | 16024605579 | 27 | 52b47153 |
| unit-p65537-a1-s10-k10 : 6 | 60460 | 1,21,211,1351,6196,21700,60460 | 1,21,211,1351,6196,21700,60460 | true / true | **true** | 43706615635 | 93 | 52b47153 |
| unit-p65537-a1-s10-k14 : 5 | 122438 | 0,0,2,154,4323,47849 | 0,0,2,154,4323,47849 | false / false | **true** | 18299157464 | 48 | 52b47153 |
| unit-p163841-a3-s21-k14 : 5 | 122438 | 0,0,2,154,4323,47849 | 0,0,2,154,4323,47849 | false / false | **true** | 18299505466 | 49 | 52b47153 |
| unit-p163841-a3-s21-k13 : 5 | 83682 | 0,0,2,150,4021,39657 | 0,0,2,150,4021,39657 | false / false | **true** | 17993415519 | 41 | 52b47153 |
| unit-p65537-a1-s10-k10 : 5 | 21700 | 0,0,2,138,3163,18569 | 0,0,2,138,3163,18569 | false / false | **true** | 16024605579 | 5 | 0dfbffad |
| unit-p65537-a1-s10-k11 : 5 | 35443 | 0,0,2,142,3441,25033 | 0,0,2,142,3441,25033 | false / false | **true** | 16340811395 | 6 | 0dfbffad |
| unit-p65537-a1-s10-k12 : 5 | 55455 | 0,0,2,146,3727,32057 | 0,0,2,146,3727,32057 | false / false | **true** | 17441038634 | 8 | 0dfbffad |
| unit-p65537-a1-s10-k10 : 6 | 60460 | 1,21,211,1351,6196,21700,60460 | 1,21,211,1351,6196,21700,60460 | true / true | **true** | 43706615635 | 16 | 0dfbffad |
| unit-p65537-a1-s10-k13 : 5 | 83682 | 0,0,2,150,4021,39657 | 0,0,2,150,4021,39657 | false / false | **true** | 17993100441 | 9 | 0dfbffad |
| unit-p163841-a3-s21-k13 : 5 | 83682 | 0,0,2,150,4021,39657 | 0,0,2,150,4021,39657 | false / false | **true** | 17993415519 | 10 | 0dfbffad |
| unit-p65537-a1-s10-k14 : 5 | 122438 | 0,0,2,154,4323,47849 | 0,0,2,154,4323,47849 | false / false | **true** | 18299157464 | 11 | 0dfbffad |
| unit-p163841-a3-s21-k14 : 5 | 122438 | 0,0,2,154,4323,47849 | 0,0,2,154,4323,47849 | false / false | **true** | 18299505466 | 9 | 0dfbffad |
| unit-p65537-a1-s10-k15 : 5 | 174437 | 0,0,2,158,4633,56649 | 0,0,2,158,4633,56649 | false / false | **true** | 19233087929 | 9 | 0dfbffad |
| unit-p65537-a1-s10-k11 : 6 | 110056 | 1,23,254,1794,9109,35443,110056 | 1,23,254,1794,9109,35443,110056 | true / true | **true** | 321522714727 | 110 | 0dfbffad |
| unit-p65537-a1-s10-k12 : 6 | 190051 | 1,25,301,2325,12951,55455,190051 | 1,25,301,2325,12951,55455,190051 | true / true | **true** | 2104737035939 | 586 | 0dfbffad |
| unit-p65537-a1-s10-k10 : 5 | 21700 | 0,0,2,138,3163,18569 | 0,0,2,138,3163,18569 | false / false | **true** | 16024605579 | 6 | 39acdeac |
| unit-p65537-a1-s10-k11 : 5 | 35443 | 0,0,2,142,3441,25033 | 0,0,2,142,3441,25033 | false / false | **true** | 16340811395 | 5 | 39acdeac |
| unit-p65537-a1-s10-k12 : 5 | 55455 | 0,0,2,146,3727,32057 | 0,0,2,146,3727,32057 | false / false | **true** | 17441038634 | 7 | 39acdeac |
| unit-p65537-a1-s10-k10 : 6 | 60460 | 1,21,211,1351,6196,21700,60460 | 1,21,211,1351,6196,21700,60460 | true / true | **true** | 43706615635 | 15 | 39acdeac |
| unit-p65537-a1-s10-k13 : 5 | 83682 | 0,0,2,150,4021,39657 | 0,0,2,150,4021,39657 | false / false | **true** | 17993100441 | 9 | 39acdeac |
| unit-p163841-a3-s21-k13 : 5 | 83682 | 0,0,2,150,4021,39657 | 0,0,2,150,4021,39657 | false / false | **true** | 17993415519 | 8 | 39acdeac |
| unit-p65537-a1-s10-k14 : 5 | 122438 | 0,0,2,154,4323,47849 | 0,0,2,154,4323,47849 | false / false | **true** | 18299157464 | 9 | 39acdeac |
| unit-p163841-a3-s21-k14 : 5 | 122438 | 0,0,2,154,4323,47849 | 0,0,2,154,4323,47849 | false / false | **true** | 18299505466 | 10 | 39acdeac |
| unit-p65537-a1-s10-k15 : 5 | 174437 | 0,0,2,158,4633,56649 | 0,0,2,158,4633,56649 | false / false | **true** | 19233087929 | 10 | 39acdeac |
| unit-p65537-a1-s10-k11 : 6 | 110056 | 1,23,254,1794,9109,35443,110056 | 1,23,254,1794,9109,35443,110056 | true / true | **true** | 321522714727 | 103 | 39acdeac |
| unit-p65537-a1-s10-k12 : 6 | 190051 | 1,25,301,2325,12951,55455,190051 | 1,25,301,2325,12951,55455,190051 | true / true | **true** | 2104737035939 | 550 | 39acdeac |

The profile definition is the same on both sides: the records' `profile` ends at `rank` and is cumulative by degree, as
the engine's `dims` (pivot degrees, cumulative). **Every row matches exactly; no discrepancy anywhere.**

Not reproduced locally: unit-p65537 and unit-p163841 k = 13, D = 6 (1 enters at full rank 313,912). One attempt with the merged
engine through `pod/run_job.sh` and a 14.5 GB memory limit (this box has 16 GB) reached 12.4 GB RSS after 37 min, the memory
guard sent SIGTERM, and the clean checkpoint then filled the disk (3.4 GB free) and the engine panicked; `job.json` recorded
`censored_memory_guard_soft` (`gate/raw-merged/k13/`). That exposed the disk-full panic, fixed in 450d5b295b (a failed
checkpoint or export is now reported and the closure continues). k = 13, D = 6 needs > 13 GB (the hky2 pod run peaked at
10.4 GB with sparser rows), so it belongs on a pod, as do the k = 14 D = 6 cells (J1/J2). An earlier pre-merge k = 12, D = 6
attempt died the same way (a 3.4 GB checkpoint at 10 min filled the disk); the merged and final binaries completed it without
checkpoints.
<!-- GATE:end -->

## 4. Negative controls

At D = 5 the profile is generic, so a control there does not discriminate (measured, `gate/negative_controls.log`: X = 4272, X = 4270 and a = 2 all give (0,0,2,138,3163,18569) at k = 10 — as the paper itself reports for random curves at k = 14). Discriminating controls (measured):

| control (k = 10, instance unit-p65537) | D | profile | 1 ∈ V_D | vs recorded (1,21,211,1351,6196,21700,60460), 1 ∈ V_6 |
|:--|:-:|:--|:-:|:--|
| target X = 11932 = x(P_1 + P_2), x(P_i) ∈ H (36 scalar solutions) | 6 | (0,0,175,1315,6160,21664,60424) | no | **differs** (codim 36 = number of solutions: V_6 = I ∩ M_≤6) |
| target X = 4272 (perturbed) | 6 | (1,21,…,60460) | yes | same: a generic X is a unit instance too (expected solutions ≈ 2^20/p ≈ 16 is not a guarantee either way) |
| a = 2 (wrong curve) | 6 | (1,21,…,60460) | yes | same, as above |

Basis-level controls (`gate/negative_controls_basis.log`): the exported basis of the true k = 10, D = 5 instance is rejected by the independent checker when the generator is built with X = 4272 or with a = 2 (`generator in span: false`); the basis exported for X = 4272 has the same profile but a different SHA-256 (086ff0a5… vs 5fa47095…) and is rejected with the true generator. The generator cross-check against the instance file fails on every override (`generator_matches_file: false`). So a wrong a or target cannot pass: either the profile/contains_one changes, or the generator check and the basis check fail.

## 5. Basis export and round trip

`--export-basis PATH` writes `hky.basis/v1` (BASIS_FORMAT.md) at the end of a closed closure and a
sidecar `PATH.json`. The writer lifts the σ-split parts to original monomial coordinates (orbit-grouped
column order 1; `+`/`−` rows leading at the same orbit are replaced by their half-sum and
half-difference, so pivots stay distinct and monic), asserts the profile recomputed from its pivots
equals the closure's, writes atomically (tmp + fsync + rename) and hashes the file read back.

`hky_basis_check --basis PATH [--closedness [--exact]]` re-reads the file with its own parser and
checks: trailing SHA-256 and the sidecar's hashes; the column list (Σ C(2k, d) columns, distinct,
degree ≤ D, the declared order rule); rows (pivots strictly increasing, entries strictly after the
pivot, coefficients in [1, p)); profile and contains_one from pivots. `--closedness` checks that the
generator NF(S_3) (built by the checker itself from p, a, b, X) lies in the span (exact reduction) and
that every product u·h, h a basis row of degree < D and u any of the 2k chain variables (the checker's
own product: exponent addition mod 2^k, not the engine's square rules), lies in the span: `--exact`
reduces every product; the default tests each product against 4 random annihilators of the span
(a product outside the span passes with probability ≤ p^−4). Together: W ⊇ V_D, and with the
engine's W ⊆ V_D, W = V_D.

Measured (`gate/unitness_786433.log`, `gate/negative_controls_basis.log`): every closure without 1 at k = 4..10 (p = 786433) and k = 10, D = 5 (p = 65537) exported, re-read and checked: all `ok`, closedness exact (k ≤ 9: every product reduced) or by annihilators (k = 10: 63,260 products, D = 6 with 1: 434,000). Power of the detector (`--drop-row`, k = 10, D = 5): removing one row (index 18568 / 18500 / 17000 / 5000) makes 12,244 / 8,676 / 5,055 / 152 products fail and the check reject. Sizes: 14.3 bytes of file per final stored engine word; k = 10, D = 5: 220.9 MB; k = 8..11, D = 4: 0.8–1.4 MB (`gate/basis_size_scaling.log`).

## 6. The depth-16 unit instance

**Prime.** The paper's unit-ness argument (Lemma `lem:unit`) needs a primitive 2^{k+1}-th root of unity,
i.e. 2^{k+1} | p − 1, so that the order-2^k subgroup H consists of squares; then a non-residue target
abscissa X with Y ≠ 0 makes S_3(x_1, x_2, X) ≠ 0 on H × H and the ideal is the unit ideal.
For k = 16 this needs 2^17 | p − 1. p = 65537 = 2^16 + 1 fails (2^17 ∤ 2^16). 196609 = 3·2^16 + 1 fails
twice: it is not prime (196609 = 7 · 28087, `factor`), and even so 196608 = 3·2^16 is not divisible by
2^17. Scanning c·2^17 + 1 for c = 1, 2, …: 131073 = 3·43691, 262145 = 5·13·37·109, 393217 = 11·35747,
524289 = 3·174763, 655361 = 7·251·373, **786433 = 3·2^18 + 1 is prime** — the smallest prime with
2^17 | p − 1 (it even has 2^18 | p − 1, so it hosts unit instances up to k = 17). It is below 2^31
(engine limit), of the same size as the paper's primes, and is the pilot's own default prime.
(2013265921 = 15·2^27 + 1 would also do; the smaller prime keeps the instance close to the paper's.)

**Instance.** `tools/hky_paper_instance.rs` reimplements `scripts/ic_bounds_big_run.py cmd_instance
--unit-nonresidue` (want = false): CPython's MT19937 seeded by `random.Random(seed)`, `randrange(1, p)`
by rejection on `getrandbits(bit_length(p − 1))`, `PrimeCurve.lift_x` with the paper's Tonelli–Shanks
`sqrt_mod`, first draw whose abscissa is a non-residue. Self-test (measured): it reproduces both paper
targets exactly — p = 65537, a = 1, seed 10 → (4271, 11720), 1 attempt; p = 163841, a = 3, seed 21 →
(73730, 28291), 1 attempt — and the regenerated instance files for unit-p65537 k = 10, 14 are identical
to the paper's in every field (k = 13/15 and p163841 differ only in the optional eval-rank fields the
paper wrote for some files). The scalar-solution count is enumerated (quadratic in x_2 for each
x_1 ∈ H); positive control: X = x(P_1 + P_2) with x(P_i) ∈ H gives 4 solutions at k = 6.

`instances/unit-p786433-a1-s10-k16.json`: p = 786433, y² = x³ + x (a = 1, b = 0, j = 1728), seed 10,
target **R = (485051, 497175)** (3 attempts), points 0. Its `unitness` block, every item checked
explicitly by the tool: 2^17 | p − 1 ✓, the order-2^16 subgroup generator 346395 is a square ✓
(H ⊆ squares), 485051 is a non-residue ✓, 497175² = 485051³ + 485051 ✓, Y ≠ 0 ✓, b = 0 ✓, scalar
solutions enumerated over H × H: 0 ✓ ⇒ unit ideal (Lemma `lem:unit`). The same (a, b, R) at every k
(the draw does not depend on k): files for k = 4..16.

**Unit-ness by running** (small depth, same prime and target; expected last fall degree
⌊(k+14)/4⌋ = 4 (k ≤ 5), 5 (k = 6..9), 6 (k = 10..13)):

| k | D | expected | profile | 1 ∈ V_D | basis check |
|--:|--:|:--|:--|:-:|:--|
| 4 | 3 / 4 | LFD 4 | (0,0,1,9) / (1,9,37,93,163) | no / yes | ok exact / ok exact |
| 5 | 3 / 4 | LFD 4 | (0,0,1,11) / (1,11,56,176,386) | no / yes | ok / ok |
| 6 | 4 / 5 | LFD 5 | (0,0,2,111,605) / full 1586 | no / yes | ok / ok |
| 7 | 4 / 5 | LFD 5 | (0,0,2,115,829) / full 3473 | no / yes | ok / ok |
| 8 | 4 / 5 | LFD 5 | (0,0,2,119,1061) / full 6885 | no / yes | ok / ok |
| 9 | 4 / 5 | LFD 5 | (0,0,2,123,1301) / full 12616 | no / yes | ok / ok |
| 10 | 5 / 6 | LFD 6 | (0,0,2,138,3163,18569) / full 60460 | no / yes | ok annihilator / ok |
| 11 | 5 / 6 | LFD 6 | (0,0,2,142,3441,25033) / full 110056 | no / yes | ok annihilator / ok (merged engine, `gate/unitness_786433_merged.log`) |

1 enters exactly at ⌊(k+14)/4⌋ at every k = 4..11 on the new prime (k = 10 also re-run on the merged engine: same basis SHA-256 ad097a63…) (measured, `gate/unitness_786433.log`), and the closures without 1 have the same profiles as the paper's p = 65537 instance (k = 10, D = 5: (0,0,2,138,3163,18569)), as expected for this generic invariant. This is the evidence that the recipe gives unit ideals on p = 786433 before it is trusted at k = 16 (where the lemma is the proof).

## 7. Checkpoint/resume

`--checkpoint-dir DIR --checkpoint-every MINUTES` (and `--checkpoint-every-batches N` for tests) write
`DIR/state.ckpt` at a batch boundary: the whole `ClosureState` (both parts' echelon rows, the queue
buckets in order, rows_reduced, contains_one, the yield/interreduction schedule counters, every stats
counter, elapsed time), a fingerprint (SHA-256 of ring, generator, D, column split and every schedule
option) and a trailing SHA-256; written as `.tmp`, fsynced, renamed, directory fsynced, plus
`DIR/state.json` (batch, rank, queue, sha). `--resume DIR` loads it (refusing a wrong fingerprint or
hash) and continues; the engine has no RNG and its schedule depends only on this state, so the
continuation performs the same batches. SIGTERM/SIGINT/SIGUSR1 make the engine write a checkpoint at
the next batch boundary and exit with a closure censored `"signal"` (exit code 3).

Measured (`gate/ckpt-k10D5*.log`, unit-p65537 k = 10, D = 5, 522 batches, 16 interreductions): an uninterrupted run against runs killed with **SIGKILL** and resumed, at three points — batch 150, batch 300 (`ckpt-k10D5-mid`), and after 22 s of wall at batch 475 = 30 % of the muladds, past several interreductions (`ckpt-k10D5-t22`) — and runs stopped with **SIGTERM** (clean checkpoint, exit 3) at batches 145, 308, 488 and resumed. Every resumed run: same columns, profile, contains_one, **muladds (16,024,605,579), rows_reduced (31,631)**, and the exported basis is **byte-identical** (content SHA-256 5fa47095…, file SHA-256 16108d90…). Checkpoint size at batch 488: 83.7 MB, written in 0.4 s.

After the engine-opt merge (rows stored in epoch coordinates; the checkpoint now also saves each part's `np_cols` and resume rebuilds with `Echelon::with_epoch`; the export reads rows through `Echelon::row_full`; fingerprint v2 also covers `np_epochs` and the sweep options) the same test on the merged binary (`gate/ckpt-merged-k10D5-t5.log`: SIGKILL at batch 475 = 30 % of the muladds, SIGTERM at batch 510 = 66 %) and on the final binary (`gate/ckpt-final-k10D5-t3.log`) gives the same muladds and a basis byte-identical to the pre-merge engine's (content SHA-256 5fa47095…): the optimizer changed neither the closure nor its exported basis. Checkpoint or export write failures (disk full) no longer stop the engine: tested with a checkpoint path on `/dev/full` (one failure recorded in `io.checkpoint_failures`, the next checkpoint written, the result and basis unchanged). Size: a checkpoint is up to ≈ 1 × RSS (k = 12, D = 6: > 3.4 GB at 10 min), and the previous one stays on disk while the next is written.

## 8. Frozen binary and pod jobs

The owner's optimizer pass (`wt/engine-opt`) will replace the engine; **the final pod binary is built
after that merge and must pass this gate again** — **done**: the merge is 9fb60e9c13 and its binary passed the full gate (`gate/run_gate.sh <new binary> <dir>` gave
`match: true` on every row; `pod/BINARY.sha256` now names it). The binaries below are the gate
binaries of this report.

| binary | SHA-256 | built from | use |
|:--|:--|:--|:--|
| `bin/pilot_gate` | d31a012b5e5d938b1b7ccb80b227ea6cf4322375b6a1170f1a587482204137bc | source of e4d5a6ec08 (stamp says 50262dc257: built before that commit, same pilot sources) | the gate rows of §3 |
| `bin/pkm_tower_pilot` first provisional build (pre-merge) | 52b47153eb732c85c17a693c8f85bd481e66b5fce90b752ab7d51217c7048d04 | e295580642 (stamped) | re-gated on a subset (§3); superseded |
| `bin/pkm_tower_pilot_merged` | 0dfbffadf7bf4c18762cc53b2aac28f51e3e6de7bfb183bbd34511d011903c08 | 9fb60e9c13 = merge of wt/engine-opt c51685b6a9 (stamped) | full gate (§3), checkpoint test (§7); superseded by the final binary |
| `bin/pkm_tower_pilot_final` (**pod binary**, `pod/BINARY.sha256`; copy it to `/workspace/bin/pkm_tower_pilot`) | 39acdeac81689a49fa98a9b71e9924bdc29418d80247bece72e372aa41cac74f | 450d5b295b (merge + non-fatal checkpoint/export failures; stamped) | full gate (§3), checkpoint test (§7), pods |

All: `cargo build --release`, no target-cpu flags, runtime SIMD (avx512/avx2/portable), dynamically linked against libgcc_s/libm/libc only, max GLIBC_2.34.

Pod scripts (`pod/`, bash + coreutils only): `run_job.sh <job> <instance> <D>` verifies the binary's
SHA-256 against `pod/BINARY.sha256`, logs host facts (CPU model, nproc, MemTotal, cgroup limit, disk,
kernel, SIMD flags, glibc, `ldd`), runs with `--checkpoint-every 30` into `/workspace/out/<job>/ckpt`,
resumes automatically when `ckpt/state.ckpt` exists, skips a job whose `job.json` already says ok,
samples VmRSS/VmHWM every 60 s into `rss.log`, and guards memory: RSS above 85 % of the effective limit
(min of cgroup limit and MemTotal) → SIGTERM (checkpoint + clean censored exit), above 97 % → SIGKILL;
either way `job.json` records `status: censored_memory_guard_*`, and an exit by SIGKILL from outside is
recorded as `censored_sigkill_likely_oom` — never a silent failure. Outputs under
`/workspace/out/<job>/`: `record.json` (engine record), `job.json`, `basis.bin` + `basis.bin.json`,
`rss.log`, `engine.err`, `job.log`. Tested locally: normal completion, the soft→hard guard with an
artificial 120 MB limit (censored record, then a plain re-run completed with the reference basis
SHA-256), and skip-if-done (`gate/raw/podscript-guard-test/`).

| job | command (pod) |
|:--|:--|
| J1 | `pod/J1.sh` = `run_job.sh J1-unit-p65537-a1-s10-k14-D6 /workspace/instances/unit-p65537-a1-s10-k14.json 6` |
| J2 | `pod/J2.sh` = `run_job.sh J2-unit-p163841-a3-s21-k14-D6 /workspace/instances/unit-p163841-a3-s21-k14.json 6` |
| J3 | `pod/J3.sh` = `run_job.sh J3-unit-p65537-a1-s10-k15-D6 /workspace/instances/unit-p65537-a1-s10-k15.json 6` |
| J4 | `pod/J4.sh` = `run_job.sh J4-unit-p786433-a1-s10-k16-D6 /workspace/instances/unit-p786433-a1-s10-k16.json 6` |

Pod layout: `/workspace/bin/pkm_tower_pilot` (the frozen binary), `/workspace/instances/*.json` (the
paper's files for J1–J3, `instances/unit-p786433-a1-s10-k16.json` for J4), this `pod/` directory.
Predicted wall/memory/disk with bands: `PREDICTIONS.md` (merged engine: J1/J2 ≈ 3.3 h, 31 GB; J3 ≈ 19 h, 75 GB;
J4 ≈ 4.3 d, 170 GB, on 32 vCPU; see the revision section).

## 9. Notes and limits

- **Basis size.** An echelon basis of V_6(14) in original coordinates is inherently large: the
  engine's rows are interreduced, so they are dense on the ~50k non-pivot columns; measured 14.3 bytes
  of `hky.basis/v1` per final stored engine word (k = 10, D = 5: 220.9 MB for rank 18,569), predicting
  ≈ 50 GB at k = 14, D = 6 (the GPU engine's free-column RREF is of the same order:
  r × (c − r) ≈ 449k × 50k entries). A compact v2 (free-column dense rows, u16/u32 coefficients) would
  cut this 2–4×; not done here because v1 is the frozen interface.
- The D = 5 profile alone is not a discriminating observable (it is the same for perturbed targets,
  wrong a, and — per the paper — random curves); the gate therefore also compares contains_one and the
  D = 6 profiles, and the negative controls use the D = 6 cell and the basis checker.
- Local walls are on 4 vCPU shared with four other workers (load 15–27); use muladds for comparisons.
