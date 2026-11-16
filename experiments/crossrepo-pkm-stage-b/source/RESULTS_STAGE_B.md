# RESULTS: PKM prime fall-degree study, stage B (m = 3–6 ladders) on a reserved RunPod pod (`pfd-b-rp`)

Every number is `[pcc:unverified]`. No pcc claim exists for any of them.

Evidence classes:
- **measured**: read from the rows under `results/RPB/6m9ewd49m2ivbn/`;
- **derived**: computed from those rows by the frozen scorer `prime_fall_degree_tools score` (`TABLE.md`, `score.json`), or
  as stated in the row;
- **cited**: a repository document says so.

## Verdict (the §6 rules as the frozen scorer applies them)

- **The first fall degree is flat.** The median `ffd` is 4 on every Kummer cell (both primes, `m = 3–6`, every `t`) and
  5 on every 2-isogeny cell (measured).
- **The last fall degree grows.** At `P1`, the Kummer `m = 3` ladder rises from 5 to 7 over `t = 2–5`. H-LFD reads
  **"grows"** for the Kummer `m = 3` ladder and for the `m = 3` null. The `m = 4` ladders read "mixed", because their last
  point is censored. The `m = 5` ladders read "insufficient".
- **H-FFD reads "FFD assumption fails"** on the `m = 3` and `m = 5` Kummer ladders: median `D − ffd` rises by at least 2.
  It reads "consistent" on the `m = 4` ladders and on every 2-isogeny ladder.
- **H-gap.** On every cell with both values, `d_F = D`: 41 draws, 0 with `d_F < D`. The F4 solving degree is intrinsic
  here, as at `m = 2` (stage A: 0 of 22).
- **H-null at `m = 3`.**
  - On every completed `N` the null's median equals the Kummer median (5, 6, 6, 7), which reads "shape-driven".
  - At the largest `N` both arms are censored (`≥ 7`). The rule text does not say how to treat censored medians, so the
    reading at the two largest common `N` is **unresolved**, as in stage A.
- **H-m** at fixed `t = 2, 3` reads "mixed" on all six ladders, because the largest `m` is censored. Stage C's reading
  ("grows in m") is not changed.

## Gates

| gate | result |
|:--|:--|
| G1 | pass on the pod (`gates/g1.log`). |
| G2 (`d_F ≤ D`) | 41 draws, 0 violations (derived). |
| G3 | `--jobs 2` run: **not completed (tool race)**. It panicked after 122 comparisons, all of which agreed (`gates/g3.jsonl`). Re-run with one worker per Amendment 21 on the session host: **123 capped runs, 0 disagreements, 0 timeouts** (`gates_rerun_a21/`). |
| G4 | not part of this job. |

## Exact values, derived from the gates

The table prints a censored closure as `≥ c` or `> c` even when `D` closes it. By G2, `d_F ≤ D`, and the closure
certifies `1 ∉ V_{c}`. So:
- Kummer `m = 5`, `t = 3`, both primes: closure complete to `c = 7` without `1`, and `D = 8`, so **`d_F = 8`**.
- Kummer `m = 6`, `t = 2`, both primes: closure `≥ 7` and `D = 7`, so **`d_F = 7`**.
- Kummer `m = 4`, `t = 4`, `P1`: `≥ 7` and `D = 7`, so **`d_F = 7`**. The Amendment 19 closure `M4-t4-k{0,1}-c7` confirms
  this exactly, on both draws (`d_F = 7`).

## The growth table: Kummer, `P1`, lower median of the two draws

Stage A gives `m = 2`; this job and Amendment 19 give `m ≥ 3`.

| `m` \ `t` | 2 | 3 | 4 | 5 | 6 | `ffd` |
|--:|--:|--:|--:|--:|--:|--:|
| 2 | 4 | 4 | 4 | 4 | 5 | 4 (theorem D3) |
| 3 | 5 | 6 | 6 | 7 | ≥ 7 | 4 |
| 4 | 6 | 7 | 7 | ≥ 7 (A19, `c = 6` complete) | | 4 |
| 5 | 6 | 8 | ≥ 8 | ≥ 7 (column limit) | | 4 |
| 6 | 7 | | | | | 4 |

The `m = 2` ladder continues: `d_F = 6` at `t = 10–13` (stage A, `RESULTS_T12_BRACKET.md`, and the Amendment 18 cell
below). Its `t = 14` cell is pending.

Censoring at the ladder ends:
- the registered `max_entries 3·10⁸`: `m = 3`, `t = 6`; `m = 4`, `t = 4`; `m = 6`, `t = 2`;
- `max_cols 2·10⁶`: `m = 5`, `t = 5`;
- the 1800 s wall budget, one draw only: `P0` Kummer `m = 3`, `t = 5`, draw 0. Its pair reads 7.

## Execution record (measured)

- **Pod.** `6m9ewd49m2ivbn`, RunPod CPU, 16 vCPU / 128 GB, US$1.04/h. Created 2026-11-14 22:24Z; job 22:43:18Z to
  08:24:43Z (9.7 h); deleted 2026-11-15 ~11:00Z.
  - `host.txt`: AMD EPYC 9575F, `nproc 16`, `rayon 12`. `MemTotal` shows the host's 527 GB, not the pod's allocation.
- **Binaries.** `binary.sha256` lists the Amendment 17 rebuild: `pkm_tower_pilot 8d49270b…`,
  `prime_fall_degree_tools 17e0445a…`, `hky_fp_tower_check b2cc55fd…`, `pfd_grid_zeros 96283361…`.
  - `host.txt`'s embedded `build-info.json` is the frozen image's descriptive file and lists the image's hashes. It does
    not describe these binaries.
- **Stages.** Every `hky-B` and `d-B` stage ended `exit 0` (`progress.txt`); `g3` ended `exit 101` (Amendment 21).
- **Results tarball.** SHA-256 in `results/RPB/results.tar.gz.sha256`.
- **Route.** Amendments 16 and 17. There is no S3 resume checkpoint, and the raw files are in Git, not S3 (no
  credentials in the container).

## Other cells read in the same period (Amendments 18 and 19)

- **`m = 2`, `t = 13`, draw 0, cap 6.**
  - Measured: `inconsistent true`, not oversize, not timed out, 12,372 s; `grid_zeros 0`.
  - Derived: cap 5 did not refute (`RESULTS_T13_BRACKET.md`), so **`d_F = 6`**. Both `t = 13` draws now read 6, so
    Amendment 7 D.5 reads **"holds at t = 13"**.
- **`m = 2`, `t = 14`, cap 6, both draws.**
  - Measured: `oversize true` at `--max-dense 12800000000` after about 5,100 s, `degree_reached 6`, not refuting.
  - Amendment 14's `-m112` re-run (`--max-dense 27900000000`) is running on 256 GB pods.
- **`m = 4`, Amendment 8 / 13 cells.**
  - Kummer `t = 4`, `c = 7`: `d_F = 7` on both draws, with 1.33·10⁹ entries at `c = 7`. This is the prediction
    `d_F(16) = 7` of Amendment 8.
  - Null `t = 4`, `c = 7`: `d_F = 7` on both draws (1.68·10⁹ entries).
  - Kummer and null `t = 5`, `c = 6`: closure complete to `c = 6` without `1` on both draws, so `d_F ≥ 7`. That is 129,040
    columns and 182,700 entries.
  - Amendment 8's reading needs `t = 5`, `c = 7`, which is not yet run.

## What this does not show

- No exponent, security claim, `S` or leaderboard row (§8).
- A flat `ffd` with a growing `d_F` is evidence on the first-fall-degree assumption for these presentations. It does not
  bound any solver's cost.
- The `m ≥ 4` ladders end censored after 2–3 points, so no slope in `t` is measured for `m ≥ 4`.
