# Predictions for the production closures V_6 at depth k = 14, 15, 16 (m = 2, no 1), 32 vCPU

Written 2026-11-18, **before any pod data of this round** (no pod was used; `pod/scale_bench.sh` is prepared, not
run). Everything here is an **estimate**, `[pcc:unverified]`. "old" = hky3 at the base commit `c2afb6b334`
(`bin/pilot_base`), "new" = hky3 of this round (`bin/pilot_opt2`, epoch coordinates + sweep, defaults).

## Inputs

| input | value | kind |
|:--|:--|:--|
| columns of V_6 | k = 14: 499,178; k = 15: 768,212; k = 16: 1,149,017 | given |
| hky2 (v1) at k = 14, pod 32 vCPU | 1.637e14 muladds, 28,747 s (5.69e9 muladd/s), peak RSS 23.1 GB, rank 448,922 (90 %), entries 2.50e9 | given (`hky_pod_runs_20261117/results/t14_progress.jsonl`) |
| hky2 at k = 13 c = 6 (1 enters at the last pivot) | 2.394e13 muladds, 2,446 s (9.79e9 muladd/s), 11.2 GB | given |
| hky3/hky2 muladds on full-rank m = 2 closures | 0.47 (t = 10 c = 6), 0.60 (t = 11 c = 6: 3.215e11 / 5.39e11) | measured (t = 11 here) / given |
| muladds exponent in the columns | 3.5 .. 4.2 (v1/v2 fits 3.2..3.7 on n₊; k = 13 → 14 gives 4.15) | estimate |
| new/old wall per muladd (same muladds) on production-shaped cells | 1.70× wall, 2.17× CPU (planted t = 10 c = 6, full no-1 closure), 1.35× (t = 11 c = 6, 110,056 columns, 99 % rank), 1.47× (t = 10 c = 5, 85 % rank); 1.16..1.35× on low-rank c = 5 proxies | measured, 4 shared vCPU |
| stored bytes new/old | ≈ 1.0 (t = 14 c = 5: 283 vs ≈ 282 MB; t = 11 c = 6: 1.19 vs ≈ 1.2 GB), RSS 1.06× at t = 11 c = 6 (sweep accumulators 81 MB) | measured |

The muladd count is **identical** old vs new (deterministic, unit test and every A/B run), so only the rate differs.

## Model

- old muladds M₁₄ = 0.70 × 1.637e14 = 1.15e14 (band 0.60..0.85 → 0.98e14..1.39e14); M_k = M₁₄ (n_k/n₁₄)^e,
  e = 3.85 (band 3.5..4.2): ×5.25 (4.5..6.1) at k = 15, ×24.8 (18.5..33) at k = 16.
- old rate at k = 14 on 32 vCPU: 7e9 muladd/s (band 5.7e9 = hky2's measured rate .. 9e9; hky3 fills more lanes than
  hky2: batch up to 256, parallel second pass). It fell 42 % from k = 13 to k = 14 for hky2 (memory and the
  2 MB-per-thread accumulator leaving cache); assume ×0.8 per depth step (band 0.6..1.0).
- new rate = old rate × g, g = 1.5 at k = 14 (band 1.2..2.5: lower end the local gain, upper end if the pod is
  bandwidth-bound and the sweep's one-read-per-batch of each pivot row and the L2-sized accumulator pay as argued in
  README §3), degrading ×0.9 per depth step (band 0.75..1.0) instead of ×0.8: the reduction no longer streams each
  reducer row per product nor scatters into an accumulator of n columns.
- memory old: 23 GB at k = 14 (band 18..28), × (n_k/n₁₄)^(1.9..2.0) (hky2's entries grew as n^1.90 from k = 13 to 14):
  k = 15 ×2.27..2.37, k = 16 ×4.9..5.3. New = old × 0.95..1.10, plus the sweep's accumulators (capped at 2 GiB;
  estimate 0.3..1 GB).

## Predicted cost of the V_6 closure on one 32-vCPU host

| k | columns | muladds (both) | old wall | new wall | old peak RSS | new peak RSS |
|--:|--:|--:|--:|--:|--:|--:|
| 14 | 499,178 | 1.15e14 (0.98e14 .. 1.39e14) | 16,400 s = 4.6 h (3.0 .. 6.8 h) | 10,900 s = 3.0 h (1.2 .. 5.6 h) | 23 GB (18 .. 28) | 24 GB (18 .. 32) |
| 15 | 768,212 | 6.0e14 (4.4e14 .. 8.5e14) | 107,000 s = 30 h (14 .. 69 h) | 63,500 s = 18 h (5.4 .. 46 h) | 53 GB (41 .. 66) | 55 GB (40 .. 74) |
| 16 | 1,149,017 | 2.85e15 (1.8e15 .. 4.6e15) | 636,000 s = 177 h (56 .. 622 h) | 335,000 s = 93 h (22 .. 333 h) | 117 GB (88 .. 148) | 121 GB (85 .. 165) |

The bands are wide because three factors are uncertain at once (hky3's muladd ratio at k = 14, the exponent, the
pod rate); the **ratio** new/old (0.6 central, 0.4 .. 0.85) is the narrower prediction. The k = 16 old central
(177 h) sits inside the task's 100–280 h estimate.

`--hky3-verify-closed on` (the closedness check) costs a second echelon copy in the part's columns (≈ ×2 memory: does
**not** fit at k = 16 on 256 GB, fits at k = 14 and 15) and measured 2.5× the closure's time at t = 10 c = 5; it is
for the gate cells, not for the production runs.

## Refutation criteria (to be scored from the pod)

- `scale_bench.sh` on the 32-vCPU pod: if new/old wall at 32 threads on the main cell is outside 0.4 .. 0.9, the rate
  model g is refuted; if the counters differ between arms or thread counts, the engine is refuted (stop).
- k = 14 production: refuted if new wall is outside 1.2 .. 5.6 h, or peak RSS outside 18 .. 32 GB, or if 1 enters
  (the expected answer is no 1, dim 448,922, as the hky2 pod run and the GPU engine found).
