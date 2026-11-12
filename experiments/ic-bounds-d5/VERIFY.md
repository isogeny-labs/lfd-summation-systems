# Independent verification of item 3 (degree-5 closure growth)

Verifier run of 2026-11-08 on the 4-core / 15 GB host, alongside the long degree-6 job (PID 16280, not
touched). Every Python run used `ulimit -v 3000000` and single-thread BLAS. The longest single run took
19.8 min, peak RSS was 1.77 GB, and the new files total under 0.2 MB. No existing file was edited.

Evidence labels follow the task: **exact** (certified computation, certificate named), **proof**
(argument written out here), **empirical**, **conjecture**. Nothing here is kernel-checked.

## Method and independence

- **My own closure engine.** `csd/ic_bounds_d5_verify.py` (sha256 `c3d123ed…`) is the engine, and
  `scripts/ic_bounds_d5_verify.py` (sha256 `1c94b5bd…`) is the driver.
  - The engine imports nothing from `csd/ic_bounds_d5.py`, from the producer's scripts, or from
    `ic_fall_certificates`, `ic_chain_shift_closure` or `ic_chain_block_closure`.
  - It writes S_3 itself and enumerates Z itself, vectorised over `x2`.
  - It reuses `ic_dreg.PrimeCurve` and `ic_curve_sweep.random_rational_target`, only to redraw the
    recorded targets from `random.Random(seed)`. It uses `ic_dreg.s3_prime_coeffs` only as an assertion
    against my own expansion of S_3.
- **What the engine computes.** It works in the group algebra `F_p[(Z/2^k)^2]`, with weight = popcount
  and `u_{i,j}` = adding `2^j` to `a_i` mod `2^k`.
  - Columns are ordered by weight, descending.
  - Each part is kept in RREF as (pivot columns, coefficients on the free columns). Storage is float32
    and products are exact float64, with every partial sum below 2^53.
  - It has an optional summand-swap split with block-1 shifts and both projections. A `--nosplit` mode
    uses the full coordinates and all 2k shifts.
  - Generator rows are either `NF(m g)` for `wt m <= D - 2`, or `g` alone (`--gens g`).
- **Certificates.** Each closure carries three:
  - **(a) Processing invariant.** Every low pivot is processed exactly once. A low row is changed only by
    adding multiples of rows inserted later, and those rows are also low. So the final low rows lie in
    the span of the states that were shifted.
  - **(b) Randomized sweep.** Take 6 random combinations of the low rows of each part, apply every shift
    and projection, and check that they reduce to 0. A space that is not closed passes with probability
    at most `4k p^-6`.
  - **(c) Full deterministic sweep.** Every shift of every final low row must reduce to 0. This is the
    certificate I rely on, and it completed and passed in every run below.
- **Soundness check on proper ideals.** Random elements of the final space vanish on Z. Separately, the
  monomials of weight < D evaluate onto Z (`rank = |Z|`).
- **Premise check.** The task's premise is that working modulo the chains equals the ordinary closure
  (D >= 2). I checked it directly: the image under NF of the ordinary-ring closure computed by
  `experiments/ic-curve-validation/code/lfd_engine.py` has the same dimension as my chain closure, and
  the same verdict on whether 1 is in V_D. That holds in 11 cases: k = 4 (unit and p = 13313, D = 3, 4, 5)
  and k = 5 (unit, D = 3, 4, 5; p = 15361, D = 4, 5). See `verify/ordinary.json`.

## Results

| check | result | file |
|---|---|---|
| 18 closures at k = 6, 7 (three families, D = 3, 4, 5) | My split engine and my no-split engine, both started from `g` alone and both full-swept, give exactly the validation.json values for every one: dimension, profile, whether 1 is in V_D, and \|Z\|. | `verify/small.json` |
| depth-8 file, 22 systems (20 at k = 8, 2 at k = 7) | All 22 reproduce: dim V_4 (1061, or 829 at k = 7), the D = 5 rank equals the saved target, and the saved low dimension matches. All sweeps pass, soundness holds and evaluation is onto. At D = 5 these are all "reached", i.e. V_5 = I cap M_<=5. | `verify/depth8.json` |
| k = 9, p = 12289, seed 10 | profile [1, 19, 172, 988, 4048, 12616], 1 in V_5, equal to the producer's k9 record | `verify/k9-p12289-s10-D5.json` |
| **k = 10, three instances** | p = 12289 seed 10 (\|Z\| = 0), p = 13313 seed 11 (\|Z\| = 98, target (7412, 6781)) and p = 15361, a = 4021, b = 9917, seed 12 (\|Z\| = 60, target (7776, 4768)). Each gives dim V_5 = 18,569 with profile [0, 0, 2, 138, 3163, 18569], part ranks 9311 + 9258, 1 not in V_5 and the full sweep passing. All three are deficient (targets 21700 / 21602 / 21640), so **LFD >= 6 (exact)**. The per-weight part pivot counts [0,0,2,70,1536,7703] / [0,0,0,66,1489,7703] equal the README's. | `verify/k10-*.json` |
| **k = 11, unit**: p = 12289, a = 1, b = 0, seed 10, target (9362, 1126), \|Z\| = 0 | dim V_5 = **25,033**, profile [0, 0, 2, 142, 3441, 25033], part ranks 12546 + 12487, 1 not in V_5. The full sweep checked 3441 low rows and passed, so **LFD >= 6 (exact)**. Closure 932 s, total 1151 s, peak RSS 1.77 GB. | `verify/k11-p12289-s10-D5.json`, `verify/logs/` |
| **k = 11, proper**: p = 12289, a = 4021, b = 9917, seed 12, target (4408, 4441), \|Z\| = 348 | Identical dimension, profile and part ranks. Target 35,095, deficiency 10,062; the space vanishes on Z, evaluation is onto, and the full sweep passed, so **LFD >= 6 (exact)**. Total 1190 s. | `verify/k11-p12289-s12-D5.json` |
| dim V_4 at k = 6..13, both instances per depth (16 closures) | Every one equals 4k^2 + 172k - 571 (605, 829, 1061, 1301, 1549, 1805, 2069, 2341), and dim(V_4 cap M_<=3) = 4k + 87 throughout. All full-swept; soundness holds on the 15 instances with Z nonempty. | `verify/d4.json` |
| prepared heavy-run instances (k = 12, 13) | Targets and \|Z\| = 1352, 1399, 1484, 1582 reproduce, evaluation is onto, and 40961 is prime with 2^13 dividing 40960. | `verify/instance-*.json` |
| Producer's README commands | The k = 6, D = 4 probe runs (shift and split, output sent to the scratchpad) reproduce the committed records' `basis_sha256` byte for byte. `ic_bounds_d5_validate.py --depths 6 --depth8-limit 1` runs and agrees. `ic_bounds_d5_analyze.py`, run in a shadow root on copies of runs/ and d4-check.json, regenerates summary.json and table.md byte for byte. | (scratch) |
| Pause and resume | `ic_bounds_d5_probe.py --k 10 --max-minutes 1` exits 3 after round 4 (rank 8616). Rerunning the same command gives `basis_sha256 3a377d8f81d607cc…`, identical to the committed k = 10 record. The producer's own resumed-run record (`runs/k10-…-resume.json`, named in `logs/k10-resume.log`) is not in `runs/`; this rerun replaces it as evidence. | (scratch; checkpoint deleted) |
| `verification.json` hashes | 43 files checked. All match, except `csd/ic_chain_block_closure.py`, which is the change the producer disclosed. | — |

## Proof checks

- **Split certificate (claim 5).** The argument is correct.
  - `V_D` is sigma-invariant, because g = NF(S_3) is symmetric and `sigma(u_{1,j} h) = u_{2,j} sigma(h)`.
    The split into V^+ (+) V^- needs p odd.
  - `W cap M_<D = (W^+ cap M_<D) (+) (W^- cap M_<D)`, and with columns ordered by weight descending,
    each summand is spanned by the rows whose pivot has weight < D.
  - Block-2 shifts follow from `u_{2,j} f = sigma(u_{1,j} sigma f)`.
  - My no-split mode, with all 2k shifts, gives the same numbers at k = 6, 7.
- **`d_5 <= 1 + 2k lambda` (claim 6).** Correct. `W' = span(g) + sum_v u_v L` lies in V_5, and
  `W' cap M_<5` lies in L, so W' is closed and contains g. One caveat: at the measured depths the bound
  is vacuous, since 1 + 2k lambda = 63,261 at k = 10 and 75,703 at k = 11, both above |M_<=5|
  (21,700 and 35,443). It only matters asymptotically, and only for upper growth.
- **`|Z| <= 2^(k+1)` (claim 7).** Correct. Fix x1. In S_3(x1, x2, X), the x2^2 coefficient is
  (x1 - X)^2.
  - If x1 != X, the polynomial is a quadratic in x2, so it has at most 2 roots.
  - If x1 = X, the x2 coefficient is -4(X^3 + aX + b) = -4 y_R^2. If y_R != 0 there is at most 1 root.
    If y_R = 0, the constant term is (3X^2 + a)^2, which is nonzero for a nonsingular curve.
  - Modulo the chains, I is the vanishing ideal of Z, because the evaluation map F_p[G] -> F_p^{H x H}
    is an isomorphism. So dim(I cap M_<=5) >= |M_<=5| - |Z|, and at k = 11 this is at least
    35,443 - 4,096 = 31,347 > 25,033.
- **LFD >= 6 from a deficient closed V_5.** If h is in I, has degree <= 5 and is not in V_5, then for
  every d <= 5 we have V_{max(d, deg h)} contained in V_5, so h is not in it. In the ordinary ring this
  relies on the modulo-chains premise; I checked that premise numerically at k = 4, 5 only.

## Label review

- Every value the producer marks exact (claims 1-4 and 8) reproduced with my independent engine and its
  deterministic full sweep. In each case the label is backed by a closure that was verified closed. No
  stall is used as a bound.
- For d_4, the identities as functions of k are correctly labelled empirical.
- Claims 9 and 10 are correctly labelled empirical. The arithmetic checks out:
  - C(16,5) - 4368 = 0 and C(18,5) - 8568 = 0 at k = 8 and 9.
  - C(20,5) - 15406 = 98 at k = 10.
  - C(22,5) - 21592 = 4742 at k = 11.
  - Equal values across instances also extend to equal part ranks and equal per-part profiles.
- Claim 11 has a measured part and a projected part.
  - The measured work ratios (0.80-0.83), the peak storage of 156.6M against 157.0M, and the throughput
    figures agree with the run records.
  - The projection arithmetic reproduces: the hour figures use the full per-phase bound divided by the
    measured throughput, and lambda is extrapolated linearly (3719, 3997, 4275).
  - The projections themselves are model outputs and untested. The README calls them "model"; the
    "empirical" label in the claim list is generous for the projected part.
  - My own engine needed about 6.3e12 multiply-adds for the k = 11 closure phase (its own accounting)
    and about 19.5 min per instance.
- Claim 12 (conjecture and methodology) is labelled appropriately.

## Not checked

- The heavy runs, at k = 12, 13, 14 and degree 6 at k = 10, are beyond the limits.
- So d_5(12) and d_5(13) are not verified.
- The 4-thread speedups are not verified.

Practical note: the producer's README commands write into the committed experiment directory by default
(validate, d4check, analyze, probe). Running them as written rewrites committed files. Use `--out` /
`--out-dir`, or a shadow root.

## Commands

```sh
cd /home/user/ml-cryptanalysis
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 scripts/ic_bounds_d5_verify.py small          # 18 closures k = 6, 7 (about 1.5 min)
python3 scripts/ic_bounds_d5_verify.py ordinary       # chain vs ordinary-ring closure, k = 4, 5 (2 min)
python3 scripts/ic_bounds_d5_verify.py depth8         # 22 saved systems (11 min)
python3 scripts/ic_bounds_d5_verify.py d4             # dim V_4, k = 6..13 (2.5 min)
python3 scripts/ic_bounds_d5_verify.py instance --k 11 --a 4021 --b 9917 --seed 12 --want   # etc.
python3 scripts/ic_bounds_d5_verify.py closure --k 10 --D 5 --full-sweep --minutes 22        # 4.6 min
python3 scripts/ic_bounds_d5_verify.py closure --k 10 --D 5 --prime 13313 --seed 11 --want --full-sweep --minutes 23
python3 scripts/ic_bounds_d5_verify.py closure --k 10 --D 5 --prime 15361 --a 4021 --b 9917 --seed 12 --want --full-sweep --minutes 23
python3 scripts/ic_bounds_d5_verify.py closure --k 11 --D 5 --full-sweep --minutes 22        # 19.2 min, 1.8 GB
python3 scripts/ic_bounds_d5_verify.py closure --k 11 --D 5 --a 4021 --b 9917 --seed 12 --want --full-sweep --minutes 23.5  # 19.8 min
```

Outputs are in `experiments/ic-bounds-d5/verify/`. Each JSON carries the sha256 of the two verifier
sources and, for closures, `basis_sha256` per part, a hash of the canonical RREF in my own column order.
It is not comparable to the producer's hash.
