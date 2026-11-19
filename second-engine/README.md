# The second closure engine: records, lift, standalone extract

The paper (`paper/curve-lfd-bounds.tex`, Section "Depths eleven to sixteen" and
Appendix "Data and artifacts") reports closures of a second engine and a third
degree lift, both written in the repository `isogeny-labs/crypto`. This directory
holds what the paper needs from that work, moved here on 19 November 2026 (owner
decision: paper material only). The runs themselves were made with binaries built
from crypto commits; [`PROVENANCE.md`](PROVENANCE.md) lists the commits, the binary
hashes, what was copied unchanged and every change made in the extract.

Every result here is a corroboration, not a certificate in the sense of the paper,
and its records label it `[pcc:unverified]`.

## Layout

| Path | Content |
|---|---|
| `records/hky_paper_runs_20261118/` | the paper runs J1–J4s (`results/`), their instances (`instances/`, the depth-4..16 unit instances over p = 786433), the pod launch log (`pod/LAUNCH.md`), the predictions written before the runs (`PREDICTIONS.md`), the basis format (`BASIS_FORMAT.md`), the lift design (`LIFT_RUST.md`), the gates of the run binaries (`gate/`), the lift validation (`lift_results/`), and two Rust tools (`tools/hky_basis_check.rs`, `tools/hky_paper_instance.rs`) |
| `records/hky_engine_opt_20261118/` | the optimisation round of the engine and its closedness check (`--hky3-verify-closed`), gate and A/B records |
| `records/hky_engine_skip_20261118/` | the lead-cover skip criterion (`--hky3-skip`, run J4s): `PROOF.md`, gate and A/B records |
| `records/hky_pod_runs_20261117/` | the files the paper cites from the 2026-11-17 runs (p = 2013265921, older engine `hky_fp_tower_fast.rs`) |
| `records/hky_engine_speed_20261117/README.md` | the engine version of those runs and its frozen binary |
| `lift/` | the Rust degree lift `hky_lift` (crypto `tools/hky_lift`, unchanged) |
| `engine/` | standalone Rust crate: the engine files (unchanged), the paper-instance loader, the basis checker and the instance generator; `gate/extract_gate.sh` |
| `gate-20261119/` | records of the local gate of the extract (below) |

Paths inside the records are those of the crypto repository: `research/<name>` is
`records/<name>` here, and `tools/hky_lift` is `lift/`.

## Build

Rust 1.89 or later (checked with 1.97.0); one target directory serves both crates:

```sh
export CARGO_TARGET_DIR=$PWD/target          # from second-engine/
(cd engine && cargo build --release --examples)
(cd lift && cargo build --release)
```

This gives `target/release/examples/pkm_tower_pilot` (the engine, `--instance` mode
only), `hky_basis_check`, `hky_paper_instance`, and `target/release/hky_lift`.
`target/` is not committed. The crates are not built in CI (CI runs the Python
self-tests and the Lean build); the gate below is run by hand.

One closure of the paper, for instance degree 6 at depth 12 on the unit instance
over p = 65537 (about 5 GB of memory and 10 minutes on 4 cores):

```sh
target/release/examples/pkm_tower_pilot \
  --instance ../experiments/ic-bounds-k14/instances/unit-p65537-a1-s10-k12.json --D 6 \
  --hky3-verify-closed on --out k12-D6.json
```

Flags (from `engine/examples/pkm_tower_paper/paper.rs`): `--export-basis PATH`
(an `hky.basis/v1` file, `records/hky_paper_runs_20261118/BASIS_FORMAT.md`),
`--checkpoint-dir DIR`, `--checkpoint-every MINUTES`, `--resume DIR`,
`--hky3-verify-closed on`, `--hky3-skip on`. The record's `profile`, `columns`,
`contains_one` and `closed` fields are the ones compared with
`experiments/ic-bounds-k14/runs/*.json`.

## Gate of the extract (19 November 2026)

`engine/gate/extract_gate.sh gate-20261119 <scratch>` on a 4-core, 16 GB host
(`gate-20261119/host.txt`, binary hashes in `gate-20261119/binaries.sha256`):

1. **Closures.** The eleven recorded closures of the paper that fit in 16 GB, from
   `experiments/ic-bounds-k14/runs/`, recomputed with `--hky3-verify-closed on`, and
   again with `--hky3-skip on --hky3-verify-closed on` (the configuration of J4s):

   | cell (instance, D) | recorded profile, 1 ∈ V_D | match | muladds = run-binary gate | wall s / peak MB | skip on: match, skipped/queued, muladds = J4s-binary gate |
   |---|---|:-:|:-:|---|---|
   | `unit-p65537-a1-s10-k10`, 5 | (0,0,2,138,3163,18569), false | true | yes | 22 / 229 | true, 5093/31630, yes |
   | `unit-p65537-a1-s10-k11`, 5 | (0,0,2,142,3441,25033), false | true | yes | 25 / 380 | true, 6553/37851, yes |
   | `unit-p65537-a1-s10-k12`, 5 | (0,0,2,146,3727,32057), false | true | yes | 29 / 470 | true, 8200/44724, yes |
   | `unit-p65537-a1-s10-k13`, 5 | (0,0,2,150,4021,39657), false | true | yes | 34 / 592 | true, 10052/52273, yes |
   | `unit-p163841-a3-s21-k13`, 5 | (0,0,2,150,4021,39657), false | true | yes | 41 / 592 | true, 10052/52273, yes |
   | `unit-p65537-a1-s10-k14`, 5 | (0,0,2,154,4323,47849), false | true | yes | 50 / 660 | true, 12114/60522, yes |
   | `unit-p163841-a3-s21-k14`, 5 | (0,0,2,154,4323,47849), false | true | yes | 49 / 660 | true, 12114/60522, yes |
   | `unit-p65537-a1-s10-k15`, 5 | (0,0,2,158,4633,56649), false | true | yes | 55 / 749 | true, 14374/69495, yes |
   | `unit-p65537-a1-s10-k10`, 6 | (1,21,211,1351,6196,21700,60460), true | true | yes | 27 / 562 | true, 49993/217000, yes |
   | `unit-p65537-a1-s10-k11`, 6 | (1,23,254,1794,9109,35443,110056), true | true | yes | 127 / 1323 | true, 92524/389873, yes |
   | `unit-p65537-a1-s10-k12`, 6 | (1,25,301,2325,12951,55455,190051), true | true | yes | 662 / 4947 | true, 140433/665460, yes |

   "match" is the record's own test (`records/hky_paper_runs_20261118/gate/run_gate.sh`):
   equal columns, profile and `contains_one`, and the closure closed. The multiply-add
   and reduced-row counts equal, cell by cell, those of the recorded gates of the run
   binaries (`records/hky_paper_runs_20261118/gate/raw-final/`, binary `39acdeac…`;
   `records/hky_engine_skip_20261118/raw/paper-gate-skip/`, binary `36a047ff…`). The
   engine's closedness re-check reduced 1,637,168 products, all to zero (in both passes;
   at D = 6 the space is all of M_≤6 and there is nothing to check).
2. **Independent checker.** `hky_basis_check --closedness` on the exported bases of
   k = 10, D = 5 and D = 6: `ok`, profile and rank as recorded, generator in the span,
   closed under the 63,260 and 434,000 chain-variable products (`basis-check/`).
3. **Checkpoint and resume** (`ckpt-k10D5.log`, the record's `ckpt_test.sh`): a run
   killed with SIGKILL after 7 checkpoints and one stopped with SIGTERM after 9, each
   resumed, give the reference closure (dims, profile, muladds, reduced rows, content
   hash) and a byte-identical basis file. The exported basis at k = 10, D = 5 (file
   SHA-256 `16108d90…`, content `5fa47095…`) is the file the run binary exported in its
   own checkpoint test (`records/hky_paper_runs_20261118/gate/raw-final/ckpt-k10D5-t3/ref.json`).
4. **Negative controls** at k = 10, D = 6 (`negctl/`): the planted target X = 11932
   (36 scalar solutions) gives profile (0,0,175,1315,6160,21664,60424) without 1, and a
   perturbed X = 4272 and a wrong a = 2 give unit instances again; all three equal the
   record (`records/hky_paper_runs_20261118/gate/negative_controls_d6.log`), muladds included.
5. **Unit tests** of the engine files: 20 of 20 pass (`unit-tests.log`,
   `cargo test --release --lib`).

The lift crate's validation suite (`lift/scripts/validate_small.sh`, lift against
direct closure and the recorded values of `experiments/ic-bounds-k14-independent`)
gave 41 of 41 cases equal (`gate-20261119/lift/`), as recorded for the pod binary in
`records/hky_paper_runs_20261118/lift_results/validate_small.log`. The lift built
here has SHA-256 `6824a8bd…`, the hash of the binary of the J1 and J2 lifts.

Not reproduced here: the depth-13 degree-6 closures and everything at depth 14 and
above (70 to 235 GB of memory on the pods), and the 2026-11-17 runs, whose
generator is the crypto pilot's general mode (not extracted; see `PROVENANCE.md`).
