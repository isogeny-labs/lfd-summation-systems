# Provenance of `second-engine/`

Everything in this directory was written in the repository `isogeny-labs/crypto`
(branch `claude/lfd-ffd-paper-eprint-status-qq8fit`, PR #4735) and moved here on
19 November 2026 by owner decision ("paper material only"). The source of every
file is commit **`7014d5dd2a008f07132eaa1e916dce7f5f73f05c`** of that repository.

## The runs and the binaries that made them

| Runs | Engine source (crypto commit) | Binary SHA-256 |
|---|---|---|
| J1, J2, J3 (`records/hky_paper_runs_20261118/results/`), the paper gate before them | `450d5b295b` | `39acdeac81689a49fa98a9b71e9924bdc29418d80247bece72e372aa41cac74f` (`pod/BINARY.sha256`) |
| J4s (skip criterion on) | `43d5c40d1f` | `36a047ff81bc0ea9c9ce5470902acbf73e90f19d4a43cc7711e6b0eaa10700e7` (`pod/BINARY_J4.sha256`) |
| J1 lift, J2 lift | `tools/hky_lift` (commits `7ca693ea52`, `ab6ba13436`) | `6824a8bdd8dd326b18c3e52fa1431c80dcf23582f9c52b3d6f375de51b3b1eaa` |
| `hky_pod_runs_20261117` (p = 2013265921, older engine `hky_fp_tower_fast.rs`) | `89a2a70` | `47d6c281…` (`records/hky_engine_speed_20261117/README.md`) |

The pod runs used the full crypto crate (built through the stub harness
`research/hky_paper_runs_20261118/harness/`, not copied here). The binaries were
not kept; the hashes above are the recorded ones.

## What was copied unchanged

Byte-identical to commit `7014d5dd2a` (SHA-256 compared after the copy):

- `engine/src/cryptanalysis/hky_fp_tower.rs`, `hky_fp_tower_fast.rs`,
  `hky_fp_tower_fast2.rs`, `hky_closure_io.rs` (crypto `src/cryptanalysis/`);
- `engine/examples/pkm_tower_paper/paper.rs` (crypto `examples/pkm_tower_paper/`);
- `records/hky_paper_runs_20261118/tools/hky_basis_check.rs`, `hky_paper_instance.rs`;
- `lift/` (crypto `tools/hky_lift/`: sources, `Cargo.toml`, `Cargo.lock`, `scripts/`).
  Built here with rustc 1.97.0 the release binary has SHA-256 `6824a8bd…`, the
  hash of the binary that ran the J1 and J2 lifts.
- `records/hky_paper_runs_20261118/`, `records/hky_engine_opt_20261118/`,
  `records/hky_engine_skip_20261118/`: the crypto `research/` directories of the same
  names, every file except the build harnesses (`harness/`, `unit/`), which are stub
  crates of the crypto crate and do not apply here;
- `records/hky_pod_runs_20261117/`: only the files the paper cites (`RESULTS.md`,
  `LAUNCH.md`, and twelve files of `results/`);
- `records/hky_engine_speed_20261117/README.md`: the description of the engine
  version of the 2026-11-17 runs and its frozen binary.

Paths inside the copied records still read `research/…`, `tools/hky_lift`,
`src/cryptanalysis/…` and `/home/user/lfd-summation-systems/…`: they are the paths of
the crypto repository at the time of the runs. Here `research/<name>` is
`second-engine/records/<name>`, `tools/hky_lift` is `second-engine/lift`, and the
engine files are under `second-engine/engine/`.

## What changed, and why

Only the crate around the unchanged engine files is new:

1. **`engine/src/cryptanalysis/f4_fp_tower.rs`** (new, 26 lines). The engine uses one
   type from the rest of the crypto crate, `f4_fp_tower::SquareRule`; it is copied
   verbatim, with the field list of `f4_fp_tower::TowerRing` (read by the loader).
   The F4 engine of that file is not extracted.
2. **`engine/src/lib.rs`, `engine/src/cryptanalysis/mod.rs`** (new). The library keeps
   the name `crypto_lib` and the module path `crypto_lib::cryptanalysis::*`, so the
   engine files and `paper.rs` compile without a change of path.
3. **`engine/examples/pkm_tower_pilot.rs`** (new, 567 lines, from 4,236). The
   crypto pilot is the generator of the systems; `paper.rs` (its `--instance` mode)
   uses its items through `use super::*`. This file holds those items, copied
   verbatim from the pilot at `7014d5dd2a`, except at the places marked `EXTRACT:`:
   - `PrimeCurve` (from `ic_boundary`) is a struct of the three fields `p, a, b` the
     loader sets and `summation` reads;
   - `CoverMap` (`--cover-map`) is an empty enum: the loader sets `cover: None`, and
     the cover branch of `summation` is cut;
   - `Pol::degree` and `Pol::to_poly` (F4 input) are cut;
   - `system` keeps the `tower` control only (asserted): the pilot's code for it is
     unchanged; the early returns of `naive` and `randbase`, the `null` draws, the
     `multiples` equations and the `ladder` quadrics are cut (the paper instances
     never reach them; the RNG passed in is not drawn from on the tower path, in the
     pilot either);
   - `main` dispatches `--instance` to `paper::run` as the pilot does and refuses
     every other mode;
   - `#![allow(dead_code)]`: the copied enums keep variants (`Kind`, `Step`,
     `Control`, `Engine`) that the instance path does not construct.
   Kept verbatim: `mulm, addm, subm, powm, invm`, `Pol` (the rest), `s3`, `Kind`,
   `Step`, `Tower` with `equations`, `square_rules`, `check`, `nonsingular`,
   `Target`, `Control`, `summation` (the rest), `tower_ring`, `raw_terms`, `Engine`,
   `is_prime_u64`, and `#![recursion_limit = "256"]`.
4. **`engine/Cargo.toml`** (new). Dependencies: `rayon`, `sha2`, `serde_json`, `libc`
   (the stop handlers of `hky_closure_io`), `rand` (the pilot's `StdRng` argument of
   `system`). The crypto crate's other dependencies are not needed.
5. **`engine/gate/extract_gate.sh`** (new): the gate of the extract; it calls the
   record's own `gate/run_gate.sh` and `gate/ckpt_test.sh` unchanged.

Not extracted: the pilot's other modes (F4, `f4_fp_tower`, `sig_fp_tower`,
`f4_tower_generic`, the field backends, the Dickson and isogeny towers, the
controls and the `--wide` and `--aux-search` paths). The 2026-11-17 runs on
p = 2013265921 went through the pilot's general mode (Kummer tower with coset
scale g), so their generator is not in the extract; their engine
(`hky_fp_tower_fast.rs`, `--engine hky2` in `paper.rs`) is.

## Engine source between the runs and `7014d5dd2a`

- From `43d5c40d1f` (J4s) to `7014d5dd2a`: `#[cfg_attr(…, allow(dead_code))]` on two
  SIMD variants (a merge of crypto `main`, `73566fe220`), `rustfmt` (`eda588ab34`) and
  clippy fixes (`4fdce86001`: an unused `powm` and the never-called sequential
  `Witness::insert` with its scratch field removed; `x % p == 0` written
  `x.is_multiple_of(p)`; `div_ceil`, `unwrap_or`, `is_none_or`). No change of logic.
- From `450d5b295b` (J1–J3) to `43d5c40d1f`: the lead-cover skip criterion
  (`--hky3-skip`, default off), packed checkpoints, and their records and proof
  (`records/hky_engine_skip_20261118/`).

So the extract is not the source of either run binary; the local gate
(`gate-20261119/`, summarised in `README.md`) ties it to them: on the eleven
recorded closures it reproduces, its multiply-add and reduced-row counts equal
those of the recorded gates of both run binaries, cell by cell.
