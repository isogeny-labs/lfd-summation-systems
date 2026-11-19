# Predicted wall time, memory and disk for the paper pod jobs

Written 2026-11-18, **before any pod run** of these jobs, from (a) the 2026-11-17 pod data
(`research/hky_pod_runs_20261117/RESULTS.md`, `results/t1[34]_hky2.jsonl`: `hky2`, 32 shared vCPU,
random curves, p = 2013265921) and (b) this directory's local gate (`hky3`, 4 shared vCPU at load
15–27, so local walls are not used as absolute rates; muladds and stored words are deterministic).
Every number is an estimate `[pcc:unverified]`. The binary these predictions are for is
**provisional** (see README §7): the owner's optimizer pass will replace it, which changes the
constant (the prediction should then be rescaled by its measured muladd/wall ratio), not the
closures.

## Inputs

| quantity | value | source |
|:--|--:|:--|
| hky2, k = 14, D = 6, full closure (1 ∉ V_6): wall | 28,747 s | pod, given |
| hky2, k = 14, D = 6: muladds / stored words / peak RSS | 1.637e14 / 2.657e9 / 23.1 GB | pod record |
| hky2, k = 13, D = 6 (1 enters at full rank): wall / muladds / words / RSS | 2,446 s / 2.394e13 / 1.072e9 / 10.4 GB | pod |
| hky3 / hky2 wall at t = 12 | ≈ 0.55 | speed2 README, given |
| hky3 / hky2 muladds at k = 14, D = 5 | 1.830e10 / 2.576e10 = 0.71 | local gate / pod |
| hky3 / hky2 stored words at k = 14, D = 5 | 35.3e6 / 26.2e6 = 1.35 | local gate / pod |
| basis bytes per final stored word (lift ×2, 8 B/entry) | 14.3 | local, k = 10, D = 5 (220.9 MB / 15.45e6 words) |
| checkpoint bytes per stored word | ≈ 5.4 | local, k = 10, D = 5 |
| per-depth growth of D = 6 muladds (13 → 14, hky2) | ×6.8 on ×1.59 columns: exponent ≈ 4.1 | pod |
| per-depth growth of stored words (13 → 14) | ×2.5 on ×1.59 columns: exponent ≈ 2.0 | pod |

Columns of M_≤6: k = 14: 499,178; k = 15: 768,212 (×1.54); k = 16: 1,149,017 (×1.50).

## Predictions (32-vCPU AMD EPYC, hky3 provisional binary)

| job | instance | wall: central [band] | peak RSS: central [band] | basis file | checkpoint (×2 while writing) | expected result |
|:--|:--|:--|:--|:--|:--|:--|
| J1 | unit-p65537 k = 14, D = 6 | **18,000 s (5.0 h)** [11,000 – 32,000 s] | **31 GB** [23 – 45 GB] | ≈ 50 GB [35 – 80] | ≈ 20 GB | closed, 1 ∉ V_6, profile (0,0,2,154,4451,73587,448922) |
| J2 | unit-p163841 k = 14, D = 6 | **18,000 s (5.0 h)** [11,000 – 32,000 s] | **31 GB** [23 – 45 GB] | ≈ 50 GB [35 – 80] | ≈ 20 GB | as J1 (recorded: identical profile) |
| J3 | unit-p65537 k = 15, D = 6 | **106,000 s (29 h)** [50,000 – 250,000 s] | **75 GB** [50 – 110 GB] | ≈ 120 GB [80 – 200] | ≈ 45 GB | closed, 1 ∉ V_6 (paper: Thm restrict); dim unknown |
| J4 | unit-p786433 k = 16, D = 6 | **560,000 s (6.5 d)** [250,000 – 1,400,000 s] | **170 GB** [110 – 260 GB] | ≈ 270 GB [170 – 450] | ≈ 100 GB | closed, 1 ∉ V_6 expected (paper: "1 ∉ V_6(k) for k ≥ 14", for its coefficients) |

How the centrals are formed:

- **J1/J2 wall**: hky2 pod wall × hky3/hky2 ratio; the two ratios available (0.55 wall at t = 12,
  0.71 muladds at k = 14 D = 5) bracket 15,800 – 20,400 s; central 18,000 s. Band: ×0.6 – ×1.8 of
  central (thread-count dependence of the schedule, shared-vCPU neighbours, the D = 6 rank-at-exit
  difference between the random curve and the unit instance, export and checkpoint I/O ≈ 2 %).
- **J1/J2 memory**: hky2 peak 23.1 GB × 1.35 (hky3 stores denser, interreduced rows) = 31 GB.
- **J3/J4 wall**: J1 × (column ratio)^4.1 per depth step: ×5.9 (14 → 15), ×5.3 (15 → 16). Band
  ×0.5 – ×2.4 (one measured step only; the exponent between steps is not known to better than ±1).
- **J3/J4 memory**: ×(column ratio)^2.0 per step: ×2.4, ×2.25.
- **Disk**: basis ≈ 1.6 × peak RSS (14.3 B per stored word, ≈ 8.7 B RSS per word); checkpoint
  ≈ 0.6 × peak RSS, and the old one stays until the new one is renamed over it.

## What this means for the pods

- J1, J2: a 32-vCPU host with **≥ 64 GB RAM** and **≥ 150 GB** free on `/workspace` (basis 50 GB,
  checkpoints 2 × 20 GB, margin). About 5 h each.
- J3: **≥ 128 GB RAM**, **≥ 300 GB** disk; about 1–3 days, so checkpoints matter (30-minute cadence:
  ≈ 45 GB per write, ≈ 1.5 min at 0.5 GB/s, ≈ 5 % overhead — consider `CKPT_MIN=60`).
- J4: **≥ 256 GB RAM** is the central need and the band reaches past it; disk **≥ 600 GB** with the
  export, or run with `EXPORT=0` if only dims/profile/contains_one are needed. On a 32-vCPU host
  the central is 6.5 days. The memory guard turns an overrun into a clean censored record with a
  checkpoint, never a silent OOM kill.

These are the numbers to score the runs against: a run outside its band is a finding to explain
(schedule, host, or a wrong extrapolation exponent), not a failure of the instrument.

## Revision for the merged engine (written 2026-11-18, still before any pod run)

The pod binary is now the merged engine (`wt/engine-opt` c51685b6a9 merged at 9fb60e9c13, sha256
`0dfbffad…`). Its counters are identical to the gate binary's on every gate cell (same muladds), so
the closures and the scaling exponents above are unchanged; only the constant moves. The optimizer
round measured 1.35–1.7× faster no-1 closures (given). Local evidence (measured, 4 vCPU, now lightly
loaded, load ≈ 3.5): k = 11, D = 6 111 s and 1.41 GB, k = 12, D = 6 587 s and 5.04 GB (2.1e12
muladds, ≈ 3.6e9 muladds/s on 4 vCPU).

| job | wall: central [band] | peak RSS | basis | checkpoint |
|:--|:--|:--|:--|:--|
| J1, J2 (k = 14) | **12,000 s (3.3 h)** [7,000 – 24,000 s] | 31 GB [23 – 45] | ≈ 50 GB | ≈ 30 GB |
| J3 (k = 15) | **70,000 s (19 h)** [33,000 – 170,000 s] | 75 GB [50 – 110] | ≈ 120 GB | ≈ 75 GB |
| J4 (k = 16) | **370,000 s (4.3 d)** [170,000 – 950,000 s] | 170 GB [110 – 260] | ≈ 270 GB | ≈ 170 GB |

Checkpoint size is revised up: at k = 12, D = 6 one checkpoint passed 3.4 GB (it filled this
machine's disk, `README.md` §3) while the closure's peak RSS was 5.0 GB, so a checkpoint is up to
≈ 1 × RSS, and while one is written the previous one is still on disk (2 × RSS). Disk needs: J1/J2
≥ 150 GB, J3 ≥ 350 GB, J4 ≥ 700 GB with the export (or `EXPORT=0`). `--hky3-verify-closed on`
(optional, re-checks closedness at the end) doubles peak memory: fits J1–J3 on a 256 GB host, not J4.

## Revision for the skip binary (J4 only; written 2026-11-18, before any pod data with it)

J4 is relaunched with `bin/pilot_skip` (sha256 36a047ff…, commit 43d5c40d1f; `pod/J4.sh`, flags
`--hky3-skip on --hky3-skip-witness leads`, `EXPORT=0`): the lead-cover criterion of
`research/hky_engine_skip_20261118/` (proof, gate PASS) skips products proved to lie in V_6, and checkpoints
are packed. Measured on the k = 14, D = 6 no-1 proxy: 1.5–1.65× fewer multiply-adds / engine seconds to the
same rank (growing with the rank); ≈ 1× on closures where 1 enters at the end. Every number `[pcc:unverified]`.

| job | wall: central [band] | peak RSS | checkpoint (two on disk while writing) | basis |
|:--|:--|:--|:--|:--|
| J4 (k = 16), skip binary | **230,000 s (2.7 d)** [100,000 – 950,000 s] | 160 GB [100 – 260] | ≈ 120 GB each, ≈ 240 GB peak [80 – 200 each] | off |

J1–J3 are unchanged (they keep `pod/BINARY.sha256`).
