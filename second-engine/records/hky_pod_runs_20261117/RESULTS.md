# Results so far (2026-11-17, 17:56 UTC) — all `[pcc:unverified]`, single instrument unless stated

Instrument: frozen `pilot_candidate` (sha256 47d6c281…), `--engine hky2`, prime 2013265921, Kummer tower
(split presentation). "value" = least degree c with 1 in V_c on every target run; a closure that ends without 1
gives "at least c + 1". Records: `results/`. Hosts were shared and loaded; times are indicative only.

## Two summands, depth 13 (pod `hky2-cpu-t14`, 32 vCPU)

| target | c = 5 | c = 6 | wall of c = 6 | muladds of c = 6 | peak RSS |
|---|---|---|--:|--:|--:|
| draw 0 | dim 39,657 of 83,682, no 1 | 1 enters (313,912 columns) | 2,446 s | 2.394e13 | 10.4 GB (whole run so far) |
| draw 1 | dim 39,657 of 83,682, no 1 | 1 enters (313,912 columns) | 2,367 s | 2.394e13 | 11.2 GB (whole run) |

Draw 0 repeats the container run (4 vCPU: 11,775 s, 2.422e13 muladds): same dims, same answer. The muladd counts
differ by 1.1 %, so the schedule depends on the thread count. Both targets give last fall degree 6 at depth 13.

Depth 14 (target 0, c <= 6) DONE: c = 5 no 1 (dim 47,849 of 122,438); c = 6 **no 1, dim 448,922 of 499,178**, 28,747 s (~8 h), peak RSS 23.1 GB. This **independently reproduces the lfd-summation-systems CuPy GPU engine** (expected: no 1, dim 448,922) on a completely different engine (Rust sparse hky2). So d_F(t=14) > 6; with the independent degree-lift (1 in V_7) the last fall degree at depth 14 is 7 — consistent with the conjecture LFD(k)=floor((k+14)/4).

## m-sweep (two targets per cell unless stated)

| m | depth 2 | depth 3 |
|--:|:--|:--|
| 4 | 6 (certified 6) | 7 (certified 7) |
| 5 | 6 (certified 6) | 8, at 226,325 columns, 257 s and 259 s (certified: lower bound 8) |
| 6 | 7 (certified 7), 61,183 columns, 10 s | at least 8: D = 7 closed without 1 (440,772 columns, 2 s); D = 8 (1,136,799 columns) censored by the 1.4e10-entry cap after 14,465 s on one target, so still undecided |
| 7 | 8, at 471,041 columns, 935 s and 934 s | at least 8: D = 8 (4,384,467 columns) not run (column cap) |
| 8 | 8, at 1,466,345 columns, about 200 s; same from two hosts | at least 8: D = 7 closed without 1 (3,947,416 columns, 180 s); D = 8 (13,962,982 columns) abandoned |
| 9 | at least 9: D = 8 closed without 1 on both targets (3,981,645 columns, dimension 1,047,766, 639 s each); D = 9 (11,075,610 columns) hit the 2,400 s budget on target 0 and was not waited for on target 1 | not emitted by the driver |

Calibration: the five cells with certified values (lfd-summation-systems `experiments/ic-bounds-split-m`) are
reproduced. New cells have no closure verification pass and no second instrument.

## optlevel

- CPU (both pods): `perf_event_paranoid` = 4, no hardware counters: T-timing only. Controls on the 32-vCPU pod
  (`results/optlevel_cpu.tar.gz`) and on the GPU pod's host (`results/optlevel_gpu.tar.gz`).
- GPU (RTX 2000 Ada): `gpu_calib.cu` builds and runs; three single-class kernels score 0.994, 0.994 and 1.000 of
  the calibrated ceiling (T-gpu-timing, `results/gpu3_*`). This is a pipeline check: the kernels are the calibration
  kernels themselves. `ncu` is refused (`ERR_NVGPUCTRPERM`), so T-gpu-ncu is not validated.
- Tool finding: on the GPU pod the tool reports "all 48 cores" although the container is allotted 6 vCPUs.
