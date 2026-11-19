# Skipping products the closure provably does not need: the lead-cover criterion for hky3 (J4)

Date 2026-11-18. Worktree `wt/paper-engine` (sparse checkout), local commits only. Engine
`src/cryptanalysis/hky_fp_tower_fast2.rs` (`--engine hky3`), driver `examples/pkm_tower_pilot.rs` (main and
`--instance` modes), build `research/hky_paper_runs_20261118/harness/`. **Every number in this file is
`[pcc:unverified]`** (no pcc claim filed; no leaderboard or `benchmarks/README.md` row is claimed). Labels:
**measured** (this machine, a script of this directory), **given** (earlier READMEs/records), **estimate**.
Not reachable from this sandbox: research-kb, the paper collection, ecdb, AWS/isolab, a PMU; no pod was used.
The machine: 4 vCPU (Xeon 2.1 GHz), 15 GB RAM, ~3 GB free disk, shared (load 2–8 during the runs below), so
walls are noisy; multiply-adds (`muladds`) and product counts are deterministic and are the primary evidence.

## 0. Summary

- **A criterion was proved and implemented** (`PROOF.md`): fix an order κ on the variables; a queued product
  x_b·r is not formed when r's leading column is the leading column of an element of
  Σ_{κ(a)<κ(b)} π T_a(V ∩ M_{≤d−1}) (d = deg r), witnessed by a product the engine has formed or computed. The
  proof shows the final span is still exactly V_D (truncation at D, chain normal form, σ-split, degree falls,
  versions/interreduction all covered); it uses only commutation of the variable operators, deg(x·h) ≤ deg h + 1,
  and the echelon's leads. It is **not** a probabilistic test and needs no closedness check to be sound.
  Flag `--hky3-skip on` (default **off**), witness mode `--hky3-skip-witness off|leads|span` (default `leads`).
- **Gate PASS with the skip on** (§4): all **110 reference cells EQUAL**; the closedness check (computed
  without the skip) on **314 closures without 1: 708,997 reductions, 0 non-zero**; faults 1–5 still detected;
  negative controls of the skip: **fault 6** (equal order) and **fault 8** (ignore the order) detected by both
  equality (22 differences each) and closedness (232 / 24,560 non-zero); **fault 7** (covers from leads of any
  degree) was *not* detected on this set (reported, §4.3). **Paper gate: 11/11 recorded closures match**
  (columns, profile, contains_one), closedness 1,637,168 reductions, 0 non-zero. Independent `hky_basis_check`
  on skip-exported bases: k = 8, 9 exact closedness, k = 10 annihilators: all `ok`. 20 unit tests pass,
  including the reference comparisons with the skip in their option grid and `variable_products_commute`.
- **Measured effect** (§3): it depends on the closure's shape.
  - **J4-shaped proxy** (k = 14, D = 6, random curve, no 1, budget 900 s each, same binary, traces):
    multiply-adds to reach rank 325,000: off 1.55e12, **leads 1.02e12 (1.52×)**, span 0.87e12 (1.79×); engine
    seconds 772 / 467 / 352 (**1.65×** / 2.19×). Off reached rank 330,317, leads 343,135, span 346,680 in 900 s.
    The ratio grows with the rank (1.2× at 300k, 1.5× at 325k).
  - Planted t = 10, c = 6 (full no-1 closure, 60,460 columns): muladds 3.64e11 → **1.37e11 (2.66×, leads)**,
    0.80e11 + 0.30e11 witness (span); wall 131–135 s → **55–57 s (2.4×)**, span 66 s (interleaved, 2 reps).
  - Closures where 1 enters at the last pivot and the paper's D = 5 cells: ≈ 1× (t = 11 c = 6: wall −6 %,
    muladds +9 %; paper k = 10…12, D = 6: −9 % to +9 % muladds; D = 5 cells +0.3 % to +5 %): what is skipped
    there is cheap.
  - Products skipped: 46 % of product items (t = 10 planted), 31 % (k = 14 leads, 900 s), 57 % (span), 16–21 %
    (paper D = 5 cells), 23 % (t = 11 c = 6).
- **Other levers for J4** (§5): checkpoints are now written **packed** (each row's columns and coefficients in
  the fewest bytes): measured −49 % (p = 65537, k = 10) and −36 % (p = 786433, k = 11) for the same state;
  resume exact (SIGKILL and SIGTERM, same muladds, byte-identical basis). Expected −25…−30 % at k = 16 (3-byte
  columns and coefficients): two checkpoints ≈ 2 × 0.75 × RSS fit the 320 GB disk with `EXPORT=0`. The skip
  does not change the final space, so the basis and RSS at the end are not reduced; RSS at matched time was
  6.2 GB (leads) vs 8.7 GB (off) on the k = 14 proxy. `VERIFY=1` still doubles memory: not for J4.
- **Revised J4 prediction** (32 vCPU, written before any pod data with this binary, §6): **wall 230,000 s
  (2.7 d), band 100,000 – 950,000 s**; **peak RSS 160 GB [100 – 260]**; checkpoint ≈ 120 GB [80 – 200], two on
  disk while writing ≈ 240 GB; basis export off. Binary `bin/pilot_skip` sha256
  **36a047ff81bc0ea9c9ce5470902acbf73e90f19d4a43cc7711e6b0eaa10700e7** (commit 43d5c40d1f, stamped), pod
  script `research/hky_paper_runs_20261118/pod/J4.sh` (flags `--hky3-skip on --hky3-skip-witness leads`).

## 1. Why skipping, and why earlier rounds did not

The opt round measured that reductions to zero are 51–88 % of the multiply-adds of a no-1 closure and that
their number, products − rank, can only fall if products are not formed; it did not try, because a heuristic
skip would have needed the closedness check (≈ 2.5× the closure) to certify each run. The criterion here is
proved, so a run with the skip is certified by the proof (and the implementation by the gate), not by a check.

## 2. The criterion and the implementation

Statement and proof: `PROOF.md`. Implementation (`hky_fp_tower_fast2.rs`):

- `kappa`: the order on the formed variables (`--hky3-skip-order 1`, default: decreasing index, so that the
  variable served first within a row's products, the last one pushed, has κ = 0).
- `cover[part][column]`: the least κ known to cover the column. Sources: every formed product of a row of
  degree d − 1 ≤ D − 2 whose lead has degree d (`ItemOut::covers`); with `--hky3-skip-witness leads`
  (default) also, when such a row is inserted, the leads of the top-degree parts of its products by every
  variable (no reduction); with `span` the labelled span of those top parts (`Witness`, `insert_batch`; more
  covers, much more work, more memory).
- The skip: when a batch is drawn, an item (r, b) with `cover[part][lead r] < κ(b)` is dropped and counted
  (`skipped_items`); the batch is filled with the next items. Nothing else changes.
- State: `cover` is part of `ClosureState` and of the checkpoint (only when the skip is on; the fingerprint
  gains `skip/v1`, the order and the witness options, so a skip run never resumes a no-skip checkpoint and
  vice versa; with the skip off the fingerprint and schedule are unchanged: same muladds as the gated pod
  binary at the same batch, measured).
- New counters: `queued_items`, `skipped_items`, `covered_cols`, `cover_hist`, `rank_by_deg`,
  `witness_muladds`, `witness_s`, `coverable_items`/`coverable_ops` (with the skip off: what it would skip),
  `ops_by_src`/`prods_by_src`. Faults 6, 7, 8: negative controls (PROOF.md §2.2).

## 3. A/B — measured

### 3.1 J4-shaped: k = 14, D = 6, random curve (no 1), 900 s budget each, traces (`raw/traj_t14c6/`)

Main mode, `cell1.sh` cell t = 14 c = 6 (the hky2 pod's k = 14 D = 6 shape: 1 ∉ V_6, ~90 % rank at the end),
`bin/pilot_skip`, arms run one after another on the otherwise idle machine (load 3.4–4.0), `traj.sh 12500`:

| rank reached | off: muladds / engine s | leads | span |
|--:|:--|:--|:--|
| 250,000 | 8.40e10 / 58 | 8.72e10 / 58 | 6.51e10 / 39 |
| 275,000 | 1.89e11 / 127 | 1.64e11 / 115 | 9.36e10 / 54 |
| 300,000 | 5.09e11 / 338 | 4.24e11 / 227 | 4.05e11 / 194 |
| 312,500 | 9.12e11 / 521 | 7.72e11 / 360 | 7.91e11 / 321 |
| 325,000 | **1.55e12 / 772** | **1.02e12 / 467** | 0.87e12 / 352 |
| 337,500 | (not reached: 330,317 at 900 s) | 1.53e12 / 707 | 1.17e12 / 477 |

At the stop: products 1,375,453 / 951,765 / 512,575; zero 936,996 / 492,290 / 142,052; items skipped 0 /
250,054 of 810,852 / 474,301 of 838,166; interreduction muladds 3.6e11 / 0.73e11 / 2.3e11 (the interreduction
schedule moves with the batches; t = 11 below shows it can also go the other way); peak RSS 8.70 / 6.23 /
8.62 GB; span witness 1.0e10 muladds, 9.2 s.

A second probe on the paper's unit instance itself (`raw/k14_unit_probe/`, `--instance unit-p65537 k14 D6`,
900 s, a development build with the same criterion, span witness; the two arms ran beside other jobs, so only
the counts are comparable): skip reached rank 344,420 with 57 % of items skipped vs 322,907 without; at the
stop 86 % of the degree-5 pivots were covered (`cover_hist`).

### 3.2 Planted t = 10, c = 6 (full no-1 closure), interleaved, 2 reps (`raw/ab_t10c6_planted/`)

| arm | wall s (2 reps) | CPU s | peak RSS MB | muladds | products | skipped / queued items |
|:--|:--|--:|--:|--:|--:|:--|
| off | 130.9, 135.1 | 232 | 604–612 | 3.641e11 | 433,961 | 0 / 216,980 |
| leads | **55.2, 57.5** | 97 | 548–572 | **1.369e11** | 236,369 | 98,796 / 216,980 |
| span | 66.7, 65.7 | 116 | 569–582 | 0.803e11 (+0.297e11 witness, 33 s) | 121,191 | 156,385 / 216,980 |

dims identical in every run ([0,19,209,1349,6194,21698,60458]); load 3.8–7.5 (a k = 14 probe ran beside).

### 3.3 Closures where 1 enters at the end, and the paper's cells

- t = 11 c = 6 (1 at the last pivot): §3.4. Products 515,081 → 331,373 (leads) → 164,825 (span, `raw/dev/`),
  muladds 3.22e11 → 2.95e11 … 3.51e11 depending on where the interreductions fall. The span witness costs
  1.4e11–2.7e11 muladds here (development build, `raw/dev/r11_*`: 163–257 s vs 118 s) and is slower overall.
- Paper gate cells (`raw/paper-gate-skip/rows.jsonl`, leads) against the gated pod binary's muladds: k = 10 D 6
  4.37e10 → 3.99e10; k = 11 D 6 3.215e11 → 3.513e11; k = 12 D 6 2.105e12 → 2.057e12 (560 s vs 550 s); D = 5,
  k = 10…15: +0.3 % … +5 %.

### 3.4 t = 11 c = 6 interleaved (`raw/ab_t11c6/`)

| arm | wall s (2 reps) | CPU s | peak RSS MB | muladds | interreduce muladds | products | skipped / queued |
|:--|:--|--:|--:|--:|--:|--:|:--|
| off | 105.4, 106.9 | 355 | 1,334–1,415 | 3.215e11 | 1.50e10 | 515,081 | 0 / 389,873 |
| leads | **100.7, 98.8** | 334–344 | 1,333–1,401 | 3.513e11 | 1.13e11 | 331,373 | 92,524 / 389,873 |
| leads, `--hky3-interreduce 0.5` (1 run, `raw/ir_t11/`) | 214.9 (load ≈ 4, beside a gate) | 405 | — | 2.945e11 | 0.63e10 | 337,351 | 98,504 |

dims identical ([1,23,254,1794,9109,35443,110056], 1 at the last pivot). Wall −6 %, muladds +9 % (the
interreductions), −8 % with a later interreduction threshold: on this shape the skip is about neutral.

### 3.5 Where the skipped work was

With the skip off the engine counts what the criterion would skip (`coverable_ops`; covers from formed
products only): k = 14 unit probe at 900 s: 213,491 items, 2.27e11 of 1.51e12 muladds (15 %). The rest of the
saving in §3.1 is schedule: fewer products per pivot means fewer reductions late in the closure, where each
costs 10⁶–10⁷ muladds. The skipped products are mostly those of degree-(D − 1) rows whose leads are reached
by products of degree-(D − 2) rows; rows that arose from degree falls whose lead no lower product reaches
(the "mutants") are never skipped — on the paper's D = 5 cells these are most of the degree-4 rows
(k = 10: 1,682 of 3,025 degree-4 leads are in the span witness at all), which is why D = LFD − 1 at small k
gains nothing.

## 4. Gate and soundness — measured

### 4.1 110 reference cells (`run_gate.sh raw/gate-skip --hky3-skip on --hky3-verify-closed on`, `raw/gate-skip.log`)

Reference = the definitional `--engine hky` of `bin/pilot_base` (sha256 a8ff8af5…, records of the opt round
reused), candidate = `bin/pilot_skip` (36a047ff…), compared by `bin/hky_engine_equality` (7a695e94…).

| set | cells | result | closedness (closures without 1 / reductions / non-zero) | items skipped |
|:--|--:|:--|:--|:--|
| m2-tower-t1-8 | 16 | EQUAL | 36 / 19,700 / 0 | 512 of 4,916 |
| m2-tower-planted | 12 | EQUAL | 48 / 165,408 / 0 | 14,146 of 41,340 |
| m2-null-t2-7 | 18 | EQUAL | 58 / 123,646 / 0 | 77,186 of 123,588 |
| m2-naive-t2-4 | 6 | EQUAL | 17 / 109 / 0 | 10 of 82 |
| m3-tower-t2-3 | 12 | EQUAL | 34 / 10,739 / 0 | 956 of 4,617 |
| m4-tower-t2-3 | 2 | EQUAL | 7 / 12,725 / 0 | 2,742 of 7,014 |
| m5 / m6 / m7 (t = 2) | 2+2+2 | EQUAL | 23 / 189,879 / 0 | 65,922 of 133,165 |
| p0-kummer-isogeny (p = 786433) | 36 | EQUAL | 86 / 184,078 / 0 | 17,394 of 45,998 |
| censor-max-cols | 2 | EQUAL | 5 / 2,713 / 0 | 60 of 677 |
| **total** | **110** | **EQUAL** | **314 / 708,997 / 0** | |
| neg-fault1–5 (as before) | | detected: 60 / 40 / 40 / 40 / 36 differences | non-zero 2,104 / 44 / 96 / 220 / 168 | |
| **neg-fault6** (skip at equal order) | | **detected: 22 differences** | 48 / 165,208 / **232** | 18,764 |
| **neg-fault7** (covers from leads of any degree) | | **not detected** (12 cells equal) | 48 / 165,408 / 0 | 15,882 |
| **neg-fault8** (skip whenever covered) | | **detected: 22 differences** | 48 / 165,208 / **24,560** | 25,054 |

**GATE PASS.** (The m5–m7 rows are summed here; `raw/gate-skip.log` has each.)

### 4.2 Paper gate (`paper_gate.sh bin/pilot_skip raw/paper-gate-skip`, `raw/paper-gate-skip.log`)

The 11 recorded closures of `research/hky_paper_runs_20261118/gate/final_gate.sh` (k = 10…15 D = 5, both
primes at k = 13, 14; k = 10, 11, 12 D = 6) with `--hky3-skip on --hky3-verify-closed on`: **11/11 match**
(columns, profile, contains_one, closed); closedness on the 8 closures without 1: 1,637,168 reductions,
0 non-zero; 361,502 of 1,681,623 items skipped. k = 13 D = 6 does not fit this machine (as before).

### 4.3 Negative controls: what they show

Faults 6 and 8 skip products the proof does not cover and are caught by both the equality gate and the
closedness check, so the checks can see an unsound skip of this kind. Fault 7 (taking covers from products whose
lead fell to a lower degree) was not caught on the gate set: either such covers are harmless on these systems
or the gate is too small to show it; the proof does not cover it and it is not used. Unit test
`the_skip_keeps_the_space` also requires fault 8 to be detected.

### 4.4 Other checks

- Checkpoint/resume with the skip on (`raw/ckpt-skip-k10D5*.log`; `gate/ckpt_test.sh` with `EXTRA_FLAGS`):
  SIGKILL and SIGTERM at batches 150–175 and at 2–3 s (batches 408–425): resumed runs have the same columns,
  profile, muladds (16,075,815,269), rows_reduced, and a byte-identical exported basis.
- Independent basis checker on bases exported with the skip (`raw/basis_check_skip.log`): k = 8, 9 (p = 786433,
  D = 4) exact closedness, 1,904 / 2,214 products, 0 failing; k = 10 D = 5 (p = 65537) annihilator mode,
  63,260 products, 0 failing; generator in span in all three.
- Unit tests (`unit/`, `cargo test --release --offline`): 20 pass, among them `variable_products_commute`
  (x_a x_b = x_b x_a and the degree bound on every monomial, t = 3, Kummer and general rules),
  `the_skip_keeps_the_space` (planted/random, t = 3, 4, split/unsplit, c = 3…6, both orders, all witness modes:
  same dims and contains_one as without the skip, closedness 0 non-zero, something skipped, fault 8 seen),
  and the v2 reference comparisons, whose 256-variant option grid now has the skip on in half of the
  variants (all three witness modes, both orders, witness caps 0–2).

## 5. Memory, checkpoints, disk for J4

- **Packed checkpoints** (all runs, not only with the skip; old checkpoints are still read): each stored row
  is written with its columns and coefficients in the fewest bytes the row needs (tags 2 and 3). Same state,
  old vs new binary: k = 10 D 5 p = 65537 batch 500: 83,173,924 → 42,016,046 B (−49 %); k = 11 D 5
  p = 786433 batch 500: 37,057,947 → 23,561,838 B (−36 %). At k = 16 columns need 3 bytes, so expect
  −25…−30 % (estimate): ≈ 0.7–0.75 × RSS per checkpoint instead of ≈ 1×.
- **The skip and memory**: same final space, so the final rows (most of the RSS) are not reduced; fewer
  products in flight and a different interreduction schedule move the peak: 6.23 vs 8.70 GB at the 900 s stop
  on k = 14 (leads vs off, the leads arm at a higher rank), 548–572 vs 604–612 MB on planted t = 10. The span
  witness keeps its own rows (+40–60 % RSS at t = 11: 1.92–2.17 vs 1.35 GB): not recommended for J4.
- **Disk at J4** (320 GB): basis export off (`EXPORT=0` in `J4.sh`: ≈ 270 GB alone); two checkpoints while
  one is written ≈ 2 × 120 GB = 240 GB (estimate) — fits with ≈ 80 GB headroom, where the earlier binary
  (≈ 2 × 170 GB) did not.
- **VERIFY at J4**: the closedness check copies the echelon (≈ 2× RSS): does not fit 256 GB. Not needed for
  soundness of the skip (proof + gate); `hky_basis_check` on a smaller export remains the independent check.

## 6. Revised J4 prediction (written before any pod data with this binary)

Base: the merged engine's prediction (`research/hky_paper_runs_20261118/PREDICTIONS.md`, given): 370,000 s
[170,000 – 950,000], 170 GB [110 – 260], 32 vCPU. Factor from §3.1: 1.5–1.65× at the highest matched rank of the
k = 14 proxy, growing with the rank; the J4 closure spends most of its work in its last 20 % of rank, beyond
what the proxy reached; the paper's own small cells gain ≈ 1×. Central factor 1.6, band 1.0–2.2.

| J4 (unit-p786433 k = 16, D = 6), 32 vCPU | central | band |
|:--|:--|:--|
| wall | **230,000 s (2.7 d)** | 100,000 – 950,000 s |
| peak RSS | **160 GB** | 100 – 260 GB |
| checkpoint (each; two on disk while writing) | **120 GB** (≈ 240 GB peak disk) | 80 – 200 GB |
| basis export | off (`EXPORT=0`) | — |
| expected result | closed, 1 ∉ V_6 (paper: 1 ∉ V_6(k) for k ≥ 14); dims unknown | — |

A run outside the band is a finding to explain (schedule, interreduction timing, host), not a failure.
The skip's own counters (`skipped_items`, `queued_items`, `covered_cols`, `cover_hist`) are in `record.json`.

## 7. Binary, build, files

Frozen binary (git-ignored): `research/hky_engine_skip_20261118/bin/pilot_skip`, sha256
`36a047ff81bc0ea9c9ce5470902acbf73e90f19d4a43cc7711e6b0eaa10700e7`, built from commit `43d5c40d1f` (stamped
`HKY_ENGINE_COMMIT`), `cargo build --release --offline`, no target-cpu flags, empty `RUSTFLAGS`, no
`.cargo/config`; runtime SIMD; links libgcc_s, libm, libc only; highest GLIBC symbol 2.34.

```bash
export CARGO_TARGET_DIR=/home/user/crypto-wt/paper-engine/target
cd research/hky_paper_runs_20261118/harness
HKY_ENGINE_COMMIT=$(git rev-parse HEAD) cargo build --release --offline --example pkm_tower_pilot
cp $CARGO_TARGET_DIR/release/examples/pkm_tower_pilot ../../hky_engine_skip_20261118/bin/pilot_skip
cd ../../hky_engine_skip_20261118/unit && CARGO_TARGET_DIR=$CARGO_TARGET_DIR/unit cargo test --release --offline
cd .. && ./run_gate.sh raw/gate-X --hky3-skip on --hky3-verify-closed on    # needs bin/pilot_base, bin/hky_engine_equality
./paper_gate.sh bin/pilot_skip raw/paper-gate-X
./traj.sh 12500 raw/traj_t14c6/*.trace.jsonl
```

Pod (J4): copy `bin/pilot_skip` to `/workspace/bin/pkm_tower_pilot_skip`, run
`research/hky_paper_runs_20261118/pod/J4.sh` (it checks `pod/BINARY_J4.sha256`, passes
`ENGINE_FLAGS="--hky3-skip on --hky3-skip-witness leads"` through `run_job.sh`, `EXPORT=0`, job name
`J4s-…` so the earlier no-skip J4 checkpoint is never resumed by it). J1–J3 keep `pod/BINARY.sha256`.

Files: `PROOF.md`, `README.md`, `run_gate.sh`, `paper_gate.sh`, `traj.sh`, `unit/` (test crate), `raw/`
(gate logs and gzipped records, A/B runs, traces, probes; `raw/dev/` holds development-build runs quoted
in §3.3).
