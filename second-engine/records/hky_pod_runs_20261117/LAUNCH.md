# Launch record, 2026-11-17: hky2 at depth 13 and 14, optlevel on rented hosts

Owner-authorized ("Do 1-2", this session), environment RunPod (isolab-eks, S3, Conductor and
research-kb are not reachable from the session container, so no resume checkpoint in S3, no
bench.record and no pcc claim exist for these runs; every number they produce is `[pcc:unverified]`).

| item | CPU pod | GPU pod |
|---|---|---|
| name, id | `hky2-cpu-t14`, `gz5n46wfmmg86r` | `optlevel-gpu`, `qkscc4xn0nnj2o` |
| hardware | cpu5m, 32 vCPU, 256 GB, EUR-IS-1; host is a 128-thread AMD EPYC (AVX-512), shared: load average 24 to 65 at start | NVIDIA RTX 2000 Ada 16 GB, 6 vCPU, 31 GB, EU-RO-1; host AMD EPYC 7443 |
| image | `runpod/base:1.0.2-ubuntu2404` | `runpod/pytorch:1.0.2-cu1281-torch280-ubuntu2404` (nvcc 12.8, ncu present, nsys absent) |
| price, cap | US$2.08/h, 86,400 s, worst case US$49.92 | US$0.24/h, 7,200 s, worst case US$0.48 |
| created | 15:5x UTC, REST v2, `pod.py` (the lfd-summation-systems launcher with a CPU branch added) | same |
| pushed (SHA-256 verified by the pod) | `pilot_candidate` 47d6c281… (the frozen binary of `research/hky_engine_speed_20261117`, 7,666,096 bytes), `optlevel.tar.gz` 84d2bb71…, `job.py` 82899a74…, `cpu_job.sh` 16b28418… | `optlevel.tar.gz` 84d2bb71…, `job.py` 82899a74…, `gpu_job.sh` 30fc9cdc… |
| job | `cpu`: optlevel probe, controls, hky2 t = 10 under `measure`; then t = 13 draws 0 and 1 (c <= 6); then t = 14 draw 0 (c <= 6, budget 79,200 s) | `gpu`: optlevel probe, `gpu probe`, `gpu calibrate`, CPU controls |

## Predictions written before any result (from `research/hky_engine_speed_20261117/README.md`)

- t = 13, draw 0: the same dims as on the 4-vCPU container (V_5: 39,657 of 83,682, no 1; 1 enters at c = 6),
  and the same `muladds` (2.4217e13) if the schedule is thread-count independent; this is not known.
  Draw 1: 1 enters at c = 6 (expected last fall degree 6 at t = 13).
- t = 14, c = 6: closure without 1, dimension 448,922 of 499,178 (value of the other repository's
  engine). A closure that contains 1 at c = 6, or another dimension, contradicts that engine.
  Cost estimate 1.2e14 to 4e14 muladds; may be off by 3x or more.
- optlevel: `perf` is not expected to work in a container (GPU host reports `perf_event_paranoid` = 4);
  if so the T-pmu path stays unvalidated and that is the finding.

## Known confounders

- The CPU host is shared and was loaded at start, so the "quiet host" condition for the optlevel
  full-parallel points is not met; the records carry the load average.
- Wall times on this host are not comparable with the container's or the earlier pods' beyond a factor.

## m-sweep (owner request 16:0x UTC: "run m=4-9 … immediately to give more data")

| item | value |
|---|---|
| pods | `msweep-a` `o3cc7qp0n24apx` (EU-RO-1) and `msweep-b` `edjws4kb7e711m` (EUR-IS-1): cpu5m, 32 vCPU, 256 GB, US$2.08/h each, cap 50,400 s, worst case US$29.12 each; both hosts shared and loaded at start (load average 36 to 73) |
| instrument | the same frozen `pilot_candidate` (47d6c281…), `--engine hky2`, Kummer tower (split presentation), P1 = 2013265921, `--random 2`, `--hky-cmax 9`, `--budget 2400` s per closure, `--hky-max-cols 3000000` |
| cells | `msweep-a`: m = 4, 5, 6; `msweep-b`: m = 7, 8, 9; each at `--t-min 2 --max-t 2`, which emits depth 2 and depth 3 (N = m x depth) |
| started | 16:11:40 UTC, `msweep_job.sh` 7d740003… |

Checks and predictions written before the sweep's new cells:

- Calibration against certified values (lfd-summation-systems `experiments/ic-bounds-split-m`): (m, k) = (4,2) 6, (5,2) 6, (4,3) 7, (6,2) 7; (5,3) lower bound 8. Local run of this binary before launch: (4,2) 6, (4,3) 7, (5,2) 6, (5,3) no 1 through degree 7 (censored at 8 by a 200 s budget), (7,2) no 1 through degree 7.
- New cells: (7,2), (8,2), (9,2) and depth 3 for m >= 5. No prediction is stated for their values; the only reading rule is that "1 enters at c" is an exact value for this instrument and a censored closure is a lower bound of the censored degree.
- This instrument has no closure verification pass of the kind the lfd-summation-systems engine has; values new to this sweep are single-instrument and `[pcc:unverified]`.

## Changes after launch (16:15 to 16:30 UTC)

- **More pods for the m-sweep** (owner asked whether a GPU could speed it up; the sparse engine has no GPU path, so
  the sweep was spread over more CPU pods instead): `msweep-m6` `gzyyvzfhwssxi2`, `msweep-m8` `t4rz3s0otppli2`,
  `msweep-m9` `kb81a8op5pd0ah`: cpu3m (cpu5m at 32 vCPU was sold out), 32 vCPU, 256 GB, US$1.76/h each, cap 28,800 s,
  worst case US$14.08 each; started 16:25 UTC with `msweep_job.sh "<m>" "2"`. `msweep-a` and `msweep-b` continue with
  their lists, so m = 6, 8, 9 may be computed twice on different hosts.
- **Error, depth 13/14:** stages B and C of `cpu_job.sh` exited at once with 127: the image has no `/usr/bin/time`.
  No closure was run. Restarted 16:26 UTC as job `bc` (`cpu_job_bc.sh`, no `time`; peak RSS sampled from
  `/proc/<pid>/status` each minute; t = 14 budget 72,000 s to stay under the pod cap).
- **Error, GPU job 2:** `tools/profile.sh ncu` exited 2, "needs gawk"; the triad could not be scored because the GPU
  calibration has no class `triad`. Job `gpu3` (`gpu_job3.sh`) installs gawk and scores the triad as `fma_f64` at
  24 bytes per operation. Its last line archives the wrong directory; the files are fetched one by one.

## 16:34 UTC

- GPU pod `qkscc4xn0nnj2o` terminated after job `gpu3` (GET returns 404). `tools/profile.sh ncu` exits 3 on it:
  `ERR_NVGPUCTRPERM`, the driver withholds GPU performance counters from container users. Records of the three
  GPU kernels are T-gpu-timing (`results/gpu3_*`).
- m = 9 and m = 8 at depth 2/3 were censored by `--hky-max-cols 3000000` (3,981,645 columns at (9,2) D = 8;
  3,947,416 at (8,3) D = 7). Reruns `sweep2` on `msweep-m9` (16:28) and `msweep-m8` (16:34) with a 12,000,000 cap,
  output tag `-cols12M`. `msweep_job.sh` now takes the cap and a tag as arguments.

## 17:02 UTC

- `msweep-b` finished (m = 7, 8, 9 at the 3,000,000 cap) and was terminated (`edjws4kb7e711m`, GET 404); its records are
  `results/b_*.jsonl`. `msweep-m6` (`gzyyvzfhwssxi2`) was terminated as a duplicate of the m = 6 run on `msweep-a`.
- `msweep-m8` `sweep2` finished: (8,3) closed without 1 at D = 7 (3,947,416 columns, 184 s and 181 s); D = 8 has
  13,962,982 columns, over the 12M cap. `sweep3` started 17:02 with `"8" "3" 20000000` (depth 3 and 4, cap 20M).

## 17:31 UTC

- `msweep-m8` (`t4rz3s0otppli2`) terminated (GET 404): (8,3) D = 8 at 13.96M columns was not going to finish in its budget.
- `msweep-a`: job `sweep` stopped with TERM during the second target of (6,3) D = 8 (the first hit the 2,400 s
  budget); `sweep6` runs (6,3) on one target with an 18,000 s budget and a 1,200,000-column cap. The cap also
  limits the depth-4 cell the driver emits afterwards.

## 17:56 UTC

- PR #4568 merged at 17:48 with the records up to 17:35; later records go in a follow-up PR.
- `msweep-m9` (`kb81a8op5pd0ah`) terminated (GET 404) after target 1 closed D = 8 without 1; its D = 9 closure was not waited for.
- Depth 13 finished (both targets, exit 0); depth 14 started 17:46 UTC on the same pod.
- Still running: `hky2-cpu-t14` (depth 14) and `msweep-a` (`sweep6`, (6,3) D = 8, 18,000 s budget).

## 18:42 UTC

- PR #4585 merged at 17:57. Both remaining runs still going: depth 14 c = 6 (56 min in, peak RSS 12.4 GB);
  (6,3) D = 8 (74 min in, RSS 73.7 GB, host load average 78 to 218).

## 19:29 UTC

- Depth 14 c = 6: 1 h 43 min in, peak RSS 14.8 GB. (6,3) D = 8: 2 h 01 min in, peak RSS 90.4 GB. No new closure line.

## 20:21 UTC

- Depth 14 c = 6: 2 h 34 min in, peak RSS 19.8 GB. (6,3) D = 8: 2 h 52 min in, peak RSS 109.4 GB, host load average 126 to 195. No new closure line.
## Note, 2026-11-17 (later the same day): the Python glue was removed

Nothing above is rewritten; the hashes in the table are those of the files as
pushed, which are this directory at commit `833a35d21e`.

- `tools/optlevel` is now a Rust crate. The job scripts call a prebuilt binary
  pushed as `/workspace/in/optlevel` (`cargo build --release --manifest-path
  tools/optlevel/Cargo.toml`; one x86-64 binary runs on any pod, because its
  calibration kernels select their instruction set at run time) instead of
  `python3 tools/optlevel/optlevel.py`. `optlevel.tar.gz` is still needed by the
  GPU jobs for `tools/profile.sh` and `kernels/gpu_calib.cu`.
- The inline `python3` JSON summaries in `msweep_job.sh`, `cpu_job.sh`,
  `gpu_job2.sh` and `gpu_job3.sh` are now `jq` filters or `optlevel summary`.
  The pods need `jq` (`cpu_job.sh` installs it with its other packages).
- `job.py` (SHA-256 `82899a74…`, five lines: run one shell script and pass its
  exit code through) is deleted from the tree. It existed only because the
  control server inside `pod.py` accepts jobs whose first word is `python3`.
  Until `pod.py` is replaced, the same launch is
  `pod.py run NAME JOB -- python3 -c 'import os,sys; os.execv("/bin/bash", ["bash", "/workspace/code/" + sys.argv[1], *sys.argv[2:]])' cpu_job.sh`:
  that one `exec` is the only Python left in a job, and it is `pod.py`'s
  requirement, not the job's.
- `pod.py` is unchanged apart from a note on its second line: it is a copy of
  the lfd-summation-systems launcher and is to be replaced, not extended.
- No pod was started, changed or contacted for this edit. The edited scripts
  were checked with `bash -n` and their `jq` filters against the committed
  result files (`results/b_m7-T2.jsonl` reproduces the summary lines in
  `results/b_msweep.log`); they have not run on a pod.

## 21:18 UTC

- Depth 14 c = 6: 3 h 31 min in, peak RSS 21.2 GB. (6,3) D = 8: 3 h 50 min of its 5 h budget. No new closure line.

## 22:35 UTC

- (6,3) D = 8 on `msweep-a`: censored by `max_entries` (1.4e10) after 14,465 s; the closure did not decide. (6,3)
  stays "at least 8". Pod `o3cc7qp0n24apx` terminated (GET 404).
- Depth 14 c = 6 still running (4 h 49 min, peak RSS 22.0 GB). Only pod `gz5n46wfmmg86r` left from this batch.

## 23:36 UTC

- Depth 14 c = 6 still running (5 h 50 min, peak RSS 22.4 GB). Only pod left. Within the pre-run estimate
  (4 to 13.5 h on 16 cores; this host is 32 shared vCPU under heavy neighbour load).

## 00:48 UTC (2026-11-18)

- Depth 14 c = 6 still running (7 h 02 min, peak RSS 22.8 GB, barely growing — deep in the dense elimination).
  Within the 4-13.5 h estimate. Only pod left.

## 02:11 UTC (2026-11-18) — depth 14 DONE, cleanup

- t = 14, m = 2, c = 6: finished 01:45:30 UTC, 28,747 s, peak RSS 23.1 GB. 1 does NOT enter; dim 448,922 of
  499,178. Matches the CuPy GPU engine exactly (independent cross-engine confirmation of the depth-14 value).
- Last pod `gz5n46wfmmg86r` terminated (GET 404). All batch pods down. Pod proxy hosts removed from the
  environment allowlist.
