# VERIFY: item 4, "exact value at depth 10" (independent check)

Verifier: an independent session, 2026-11-08. I wrote my own closure for this check and did not use the producer's
closure code (`csd/ic_bounds_exact10.py` and `scripts/ic_bounds_exact10_*.py`). The producer's
scripts were run only to test that the README commands work.

The producer claims no exact LFD at depth 10, and this check does not establish one either. What
was checked is the validation grid, the mid-scale reproductions, the k = 8, D = 6 runner result,
the two proofs, the instance construction, and the cost and memory arithmetic.

## Independent engine

`scripts/ic_bounds_exact10_verify.py` (sha256 `a6f29e0c…2dd9`) is written from the definition. Its
only repository import is `csd/ic_dreg.py`, for `PrimeCurve` and the `S_3` coefficients.

- **Columns.** The squarefree chain masks `e1 | e2<<k` of weight `<= D`, sorted by weight
  (descending), then by mask.
- **Rows.** Kept in full RREF with leftmost pivots, stored as uint16 and multiplied in float64 tiles
  (exact, since `r p^2 < 2^53`). With that ordering, the rows whose pivot has weight `< D` span
  `V cap M_<D`.
- **Shifts.** Multiplying by `u_{i,j}` adds `2^j` to block `i` mod `2^k`.
- **`closed`.** A deterministic worklist: every inserted low row is shifted by all `2k` variables,
  and every shift is reduced and inserted, until the list is empty. The result is exactly V_D.
- **`reached`.** rank = `|M_<=D| - rank ev(M_<=D)`. Every inserted row is a genuine element of V_D,
  so the rank cannot exceed the target. Random sparse mixing (6 terms) is used only to reach the
  target faster. A stall falls back to the deterministic worklist, so no stall is ever reported as
  a bound.
- **Instance data.** The generator is `NF(S_3(c t1, c t2, x))`. The subgroup generator is the first
  `h^((p-1)/2^k)` of exact order `2^k`. The points are found by vectorised evaluation.

The generator and the point set match the stored depth-8 and depth-9 systems exactly (both
`generator_matches_stored` and `points_match_stored` are true).

**Self-test (exact).** On 12 systems at k = 3, 4 (p in {257, 7681, 12289}), my LFD and falls match
the independent ordinary-ring validator `experiments/ic-curve-validation/code/lfd_engine.py`
(`verify/selftest.json`). One of the 12 is a proper ideal; all have LFD 4.

## Verdicts

| # | claim (short) | label | verdict | evidence from this check |
|---|---|---|---|---|
| 1 | block engine agrees with dense on "1 in V_D", 72/72 | exact | **confirmed** | My closure gives dim V_D, dim(V_D cap M_<D) and "1 in V_D" on all 72 grid rows. They equal the producer's `ref_*` values and the block engine's recorded `contains_one` in 72/72 (and 48/48 in the v2 repeat): `verify/compare.json`. I re-executed the current block engine (sha `de5e3970…`, as in `hashes.json`) on 5 undecided rows, fed with my own V_{D-1}. It gave `contains_one=False` and sym+anti = my dim V_D on 5/5 (`verify/block-rerun.json`). Counts confirmed: 33 unit rows, 39 proper rows, 45 rows with 1 not in V_D, 6 of them unit (k = 6, 7, D = 4). |
| 2 | in the 45 undecided rows, the block-engine stall rank equals dense dim V_D (observed only) | empirical | **confirmed** | 45/45 against my exact dims. The label is right: a block-engine stall is not a certificate. |
| 3 | split engine matches dense on 72/72 from scratch and 72/72 from V_{D-1} | exact | **confirmed** | The recorded `exact_dimension`, `low_dimension`, `contains_one`, `target` and `status` equal my values on 72/72 rows, and `split_from_prev_sparse` equals them on 72/72. At k = 7, D = 6, where the producer's reference was local-to-global and not a direct closure, my direct closure reaches the same dimension on all 8 instances. |
| 4 | negative control rejected in 56/56 applicable rows | exact | **confirmed** (weak control) | 56/56 applicable rows. In every one, `expected_missing` agrees with my exact dims, and `detected == expected_missing`. Rerunning the README `validate` command (k = 5, two instances, into the scratchpad) re-executed the control: 6/6 detected. Caveat: the control space is far smaller than V_D (for example rank 470 against 636), so it does not test sensitivity to a deficit of a few dimensions. |
| 5 | split_profile reproduces exact LFD and falls, 24/24 (LFD 4 at k = 5; LFD 5 at k = 6, 7; falls [4] / [4,5]) | exact | **confirmed** | My profiles over D = 2..6 give the same LFD and falls on all 24 (instance, k) pairs (`verify/grid-k56.jsonl`, `verify/grid-k7.jsonl`). The certificate is `reached` with ev(M_<D) onto Z, plus exact closed dims below. |
| 6 | mid-scale: k = 8 systems 0, 14, 6 and k = 9 system 0 | exact | **confirmed** | k = 8: dim V_4 = 1061 (closed, deficient: targets 2265, 2517, 2435). V_5 reaches 6633, 6885 (1 in V_5, unit) and 6803, with low dims 2265, 2517, 2435. k = 9 system 0: dim V_4 = 1301 (closed, target 4024). V_5 reaches 12592 with low dim 4024. All four have exact LFD 5. 1061 and 1301 equal 4k^2+172k-571, which is a fit for k = 6..10 and not a law. |
| 7 | runner at k = 8, D = 6 (p = 257, x = 190, 252 points) reaches dim V_6 = 14641 | exact | **confirmed** (value) | My closure, started from my own V_5, reaches 14641 = 14893 - 252, with ev(M_<6) onto Z (rank 252), in 1098 s. The runner's interrupt-and-resume behaviour was not re-executed (only the log was read). The V_5 sha `636d184a…` was not compared, because my basis has a different representation. |
| 8 | soundness of the split certificates | proof | **confirmed** | I checked the argument: sigma is an automorphism fixing g (c1 = c2, S_3 symmetric); `(1+sigma)W (+) (1-sigma)W = W + sigma W` for p odd; the target follows from the split semisimple quotient; closed plus below-target gives LFD >= D+1; reached plus onto gives LFD <= D. `SplitClosure.verify_pass` does reduce all 2k shifts of every low row of both parts, and closure is declared only after a pass returns 0. The numbers agree with an engine that does not split. |
| 9 | (c): LFD >= 7 iff f_6(k) < delta_5(k); delta_5(10) = 3131 | proof | **confirmed** | The reduction (V_6 = I cap M_<=6 iff V_6 cap M_<=5 = I cap M_<=5, given ev(M_<5) onto) is correct as written: write m = x_v m', pick q in M_<5, and x_v(m'-q) is in V_6. 21700 - 18569 = 3131; dim V_5 = 18569 comes from the stored, verified-closed depth-10 V_5 (`experiments/ic-chain-depth10/verification.json`). I also checked ev(M_<5) onto Z for the two proper depth-10 instances (98 and 60 points, rank 98 and 60; targets 21602 and 21640 as stored: `verify/k10-onto.json`). f_5(10) = 1614 and delta_4(10) = 4647 check. |
| 10 | original depth-10 runner set-up > 15 GB; four OOM kills near 22:34 | exact | **confirmed** (attribution unverifiable) | Shape arithmetic: 18569 x 60460 x 8 B = 8.98 GB; V_5 int64 = 3.22 GB; 30318 sym columns give a 4.50 GB projection. dmesg shows four cgroup OOM kills at 11096-11145 s uptime, with total-vm 16.6 GB and anon-rss 12.65-13.35 GB. dmesg does not name the script. **Stale statement:** the README says the patched run (pid 11870) "was still going at the time of writing", but dmesg shows pid 11870 OOM-killed at 14457 s uptime (about 23:31), with anon-rss 8.85 GB. So the 22:36 patch did not by itself keep the run alive. pid 16280 is the current run. |
| 11 | dense split engine cost (6.6/17/45/113/520 GB; time brackets) | empirical | **confirmed**, with a caveat on the cols^2 column | `scripts/ic_bounds_exact10_plan.py --out <scratch>` reproduces `plan.json` exactly. Column counts check: 60460, 110056, 190051, 313912, 499178. The uint16 two-part storage is 3.66, 12.1, 36.1, 98.5 GB, and uint32 at k = 14 is 498 GB. Prime facts (exact): 40961 is the only prime = 1 mod 2^13 below 2^16, and there is none = 1 mod 2^14 (the next is 65537). The cols^3 hours reproduce by least squares on `calibration.json` (9.8, 59, 306, 1377, 5537 h). The cols^2 hours (1.8, 6, 18, 49, 123) do **not** reproduce from the file: least squares gives 2.2, 7.5, 22, 61, 153, and the fit method is not documented. All of these extrapolate 4-34x beyond the calibrated column range (at most 14893). The label "empirical, order of magnitude" is appropriate. |
| 12 | k=10 instances from dry run (x = 9362, 0 points; x = 7412, 98 points; x = 7776, 60 points) | exact | **confirmed** | My own target draw (`random.Random(seed)`, `randrange(1, p)` until x is a curve abscissa, redrawn when solutions are wanted) and point count give the same x and point counts, each on attempt 1 (`verify/dryrun-k10.json`). The README dry-run commands also run and print the logged lines. The generator sha256 values were not compared, because the hash format is the producer's. |
| 13 | k_7 in 12..17 (14 under equal gaps) | conjecture | **confirmed** as labelled | It is labelled as a conjecture and is not used as a bound. The PKM Magma report (d_reg 6 at k = 10, 11) was not checked here. |

## README commands run

- `scripts/ic_bounds_exact10_validate.py --ks 5 --instances rand-p257-x58 proper-p12289-seed43 --out <scratch> --tag=-verify`
  ran, and all 6 rows and 2 profiles pass. It was redirected to the scratchpad so that the
  producer's files stay unchanged.
- The three `scripts/ic_bounds_exact10_degree6.py --dry-run --k 10 …` commands ran.
- `scripts/ic_bounds_exact10_plan.py` ran (`--out` to the scratchpad) and its output is identical
  to `plan.json`.
- `ic_bounds_exact10_summarize.py` was **not** run, because it rewrites existing files
  (`summary.json`, `cost.svg`, ...).
- The heavy runs were not run.

## Not checked

- Any k = 10, D = 6 closure (heavy runs; producer's estimate 7 GB and about 10 h).
- The runner's checkpoint and resume.
- The k = 11..13 plans beyond arithmetic.
- The main session's runner (pid 16280).

## Reproduce

All runs used `ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`.
The commands and output hashes are in `verify/README.md`.
