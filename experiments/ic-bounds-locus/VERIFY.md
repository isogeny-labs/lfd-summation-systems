# Independent verification of ITEM 2: the exceptional locus

Verifier: a separate agent, 2026-11-08. This report checks the producer's claims in `README.md`.
It does not edit any producer file. The verifier's code and outputs are new files:
`scripts/ic_bounds_locus_verify*.py` and `experiments/ic-bounds-locus/verify/`.

## Independence and method

* **Own closure engine** (`scripts/ic_bounds_locus_verify.py`). It computes V_D modulo the chains and
  imports nothing from the producer's modules or scripts. It also avoids `ic_fall_certificates`,
  `ic_chain_shift_closure`, `ic_chain_block_closure` and `ic_curve_sweep`.
  * Generator: S_3 expanded by hand.
  * Zero set: brute force over H x H.
  * Linear algebra: only the `Ech` echelon class of `experiments/ic-curve-validation/code/lfd_engine.py`.
  * Dense columns for M_<=D: 794 or 1586 columns at k = 6 (D = 4 or 5), 1471 or 3473 at k = 7.
  * The loop stops only when every variable shift of a basis of V_D ∩ M_<=D-1 reduces to 0. Every
    result is then re-checked explicitly (`closed: true` on every row). So no value below is a stall.
* **LFD certificate.** dim V_D0 = |M_<=D0| - rank(eval of M_<=D0 on Z), and eval of M_<D0 has rank |Z|
  (local-to-global). For a unit ideal the certificate is 1 in V_D0. The exact LFD is the largest D with
  V_{D-1} != I_{D-1}, cross-checked against the last strict fall.
* **Ordinary-ring cross-check** (`ic_bounds_locus_verify_ordinary.py`). `lfd_engine.closure` runs in
  the full 12-variable polynomial ring with the chain equations.
* **Stratified subsample.** Seeds `random.Random(777+p)` for the targets, 4242 for the plane, 99 for the
  cross-depth points.
  * Every deviating target was recomputed: 44 at k = 6, 52 at k = 7.
  * Generic targets: per curve, one with t in H, one with t not in H, and 2 random. That is 52 at k = 6
    and 60 at k = 7.
* **Resources.** Every run used `ulimit -v 3000000` and one BLAS thread, with at most 2 processes at a
  time. Each run took under 25 minutes. Peak RSS was about 0.3 GB. The output is 0.8 MB.

## Commands (run from the repository root, with ulimit -v 3000000; OMP/OPENBLAS/MKL_NUM_THREADS=1)

```
python3 scripts/ic_bounds_locus_verify_sample.py         # verify/rowchecks.json, verify/sample.json (all 3185 rows)
python3 scripts/ic_bounds_locus_verify_plane.py          # verify/plane-checks.json, verify/plane-sample.json
python3 scripts/ic_bounds_locus_verify_proofs.py         # verify/proof-checks.json
python3 scripts/ic_bounds_locus_verify_predict.py        # verify/predict-checks.json
L=experiments/ic-bounds-locus/verify/lists; O=experiments/ic-bounds-locus/verify
python3 scripts/ic_bounds_locus_verify.py --prime 193 --k 6 --dmax 5 --targets "$(cat $L/p193-k6-A.txt)" --out $O/p193-k6-A.jsonl   # ~12 min
python3 scripts/ic_bounds_locus_verify.py --prime 193 --k 6 --dmax 5 --targets "$(cat $L/p193-k6-B.txt)" --out $O/p193-k6-B.jsonl   # ~12 min
python3 scripts/ic_bounds_locus_verify.py --prime 257 --k 7 --dmax 4 --targets "$(cat $L/p257-k7.txt)"    --out $O/p257-k7.jsonl      # ~10 min
python3 scripts/ic_bounds_locus_verify.py --prime 193 --k 6 --dmax 4 --targets "$(cat $L/plane-1.txt)"    --out $O/plane-1.jsonl      # ~10 min
python3 scripts/ic_bounds_locus_verify.py --prime 193 --k 6 --dmax 4 --targets "$(cat $L/plane-2.txt)"    --out $O/plane-2.jsonl
python3 scripts/ic_bounds_locus_verify.py --prime 193 --k 6 --dmax 4 --targets "$(cat $L/inv-p193-k6.txt)" --out $O/plane-inv6.jsonl
python3 scripts/ic_bounds_locus_verify.py --prime 257 --k 7 --dmax 4 --targets "$(cat $L/inv-p257-k7.txt)" --out $O/plane-inv7.jsonl
python3 scripts/ic_bounds_locus_verify.py --prime 257 --k 6 --dmax 4 --targets "$(cat $L/crossk.txt)"     --out $O/plane-crossk6.jsonl
python3 scripts/ic_bounds_locus_verify.py --prime 257 --k 7 --dmax 4 --targets "$(cat $L/crossk.txt)"     --out $O/plane-crossk7.jsonl
python3 scripts/ic_bounds_locus_verify_ordinary.py --prime 193 --k 6 --degrees 3,4 \
   --targets '22,23,1;1,0,40;1,0,8;0,1,0;0,1,85;111,126,28;5,0,0;1,0,0' --out $O/ordinary.jsonl
python3 scripts/ic_bounds_locus_verify.py --prime 257 --k 7 --dmax 5 --targets 72,135,228 --out $O/p257-k7-d5-probe.jsonl  # 97 s
python3 scripts/ic_bounds_locus_verify_compare.py        # verify/compare.json
```

The `--out` files are appended to, so delete them before a rerun. The list files were derived from
`sample.json` and `plane-sample.json`: deviating targets, then generic ones, split in half at k = 6. The
normalized points are (a/t^2, b/t^3) of the t != 0 sample targets that the plane sample does not
already contain.

## Results

### Row-level checks on all 3185 producer rows (`verify/rowchecks.json`)

* **Target lists are complete.** Every t with t^3+at+b a square, including 0, is present on every curve
  and nothing else is. Counts: 890 + 386 at k = 6, and 1153 + 756 at k = 7, as claimed.
* **Generator and |Z|.** The stored generator equals the verifier's NF(S_3) on all 3185 rows, and the
  stored |Z| equals a brute-force count on all 3185.
* **Internal consistency.** Stored falls, certificates and exact LFD agree with the stored dimensions on
  every row.

### Recomputation of targets (`verify/compare.json`)

* **208 targets** (96 at k = 6 to D = 5, 112 at k = 7 to D = 4): every dim V_D, low_D, "1 in V_D", |Z|,
  fall list and exact LFD agrees with the producer. That is **0 mismatches**, and every closure was
  verified closed.
* **k = 6 certificates.** 95 targets are certified by the verifier's own local-to-global certificate and
  1 by a unit ideal, (5,0,0) with 1 in V_4. The LFD agrees on all 96: 12 targets have LFD 4 (the 3
  with t = 0, ab = 0, and 9 C2), and 84 have LFD 5. At k = 7, the unit ideal (3,0,0) also has LFD 4 in
  both computations.
* **Class values all reproduce.** k = 6 | k = 7:
  * generic: 13/605/111 | 15/829/115
  * order two: 607 | 835
  * tangent: V_3 = 34 | 40, V_4 = 603 | 828
  * C1: 602 | 826, and 607 at the special point (32,103)
  * t = 0 with ab != 0: 612 | 849
  * t = 0 with a = 0: 44/773 | 52/1052
  * t = 0 with b = 0, a in H: 730 | 1343
  * t = 0 with b = 0, a not in H: unit | unit
  * C2: saturated | 857
* **Ordinary ring** (`verify/ordinary.jsonl`, lfd_engine, 8 targets x D = 3, 4, one per class plus
  generic). Both identities hold in all 16 cases: dim V_D(ordinary) = C(12+D,D) - |M_<=D| + dim V_D(chain),
  and the same for the low part.
* **Plane sample, p = 193, k = 6.** I recomputed 795 points:
  * every C1 point (186) and every C2 point (190);
  * the 3 special points, 36 points on 1+A+B = 0 and 30 on A = -3;
  * 312 generic points (300 random plus the normalized points of the generic sample targets).

  **0 dimension mismatches.** The one flagged row is definitional. At the cusp (-3,2), V_3 = I_3 and
  V_4 = I_4, but evaluation is not onto Z (rank 113 < |Z| = 127). The producer's flag
  `saturated4 = (codim = |Z|)` therefore reads false there.

### Point-set claims (`verify/plane-checks.json`)

* **Plane coverage.** The plane files hold 37249 distinct (A,B), all at t = 1. Class counts are exactly as
  claimed: 36486 / 193 / 191 / 190 / 186 / 1 / 1 / 1.
* **C1 cubic.** Over F_193 the cubic has exactly 190 zeros: the 186 class-602 points plus (-3,2),
  (-3,170), (162,58) and (32,103). No generic point lies on it. A refit through the 186 points has rank
  9 out of 10, and the unique cubic equals the claimed one.
  * F_257: the class-602 points in the line files plus `p257-crossk.json` give 19 points, rank 9, and the
    same coefficients as `c1-lift.json`. No other computed point lies on it.
  * F_449: 10 points, rank 9, a unique cubic.
* **C2.** The evaluation matrix has rank 190 out of 190 monomials of degree <= 18 (and 171 out of 171 at
  degree 17). Per line A = const the counts are 0:77, 1:67, 2:32, 3:11, 4:5, 6:1, and V_4 + |Z| = 794 at
  every point.
* **Lines.** 1+A+B = 0 has 192 points of class 607 plus the cusp. A = -3 has 191 tangent points plus the
  cusp and (-3,170).

### Proof claims (`verify/proof-checks.json`)

* **Identities hold.** I checked each with sympy:
  * S_3 with b = -t^3-at is exactly (x1x2 - t(x1+x2) - a - 2t^2)^2.
  * S_3 at t = 0, b = 0 is (x1x2 - a)^2.
  * S_3 is weighted homogeneous of weight 4.
  * The weight-2 matrix is [[1,-2t],[-2t,-2(a+t^2)]] with determinant -2(a+3t^2).
* **Tangent construction.** I instantiated it explicitly on 18 instances: (p,k) = (97,3), (97,5),
  (193,5), (193,6), (257,7), (257,8), each with 3t^2+a = 0. Every time h1 and h2 have weight 2, u12u21
  and u11u22 appear only where the proof says, and rank{g,h1,h2} = 3. Off the line, h1 has weight 3.
* **I read the written argument and it is correct.** With l_i = u_i1 - 2t u_i0 we have
  (u11+2tu10) l1 = u12 - 4t^2 u11, which needs k >= 3. And V_2 = span{g} because g has weight 2.
* **H-scaling.** On 10 random cases, NF S_3(x; mu t, mu^2 a, mu^3 b) = mu^4 phi(g) with
  phi(x^e) = mu^-e x^e, where mu is in H. phi is a diagonal, weight-preserving map that commutes with
  the shifts, so the proof holds.

### Invariance and depth

* **(A,B) invariance, tested on both sides with verifier values: 188 of 188 agree.** At k = 6 that is
  38 targets with t in H and 50 with t not in H. At k = 7, which the producer did not claim, it is 49 and
  51.
* **Cross-depth sample at p = 257.** 15 of 15 sampled points (4 C1, 4 order two, 4 C2, 3 tangent)
  reproduce at both k = 6 and k = 7.
* **New exact datum (k = 7, D = 5).** The C2 target t = 228 on (72,135), |Z| = 84, has V_5 = 3389 =
  3473 - 84 = I_5, with M_<=4 evaluating onto Z. That certifies **exact LFD 5**, with falls (4,5) and
  V_4 = 857 deficient. So at k = 7, C2 does not keep LFD 4, at least on this target. This answers one of
  the producer's open questions, for one target only.

### Held-out predictions (`verify/predict-checks.json`)

I re-applied the frozen rules with my own code, and the predicted sets are identical to the producer's
files.
* **k = 6:** 15 of 15 deviations predicted, 0 false alarms, 1 value miss (t = 89 on (63,189), 607 instead
  of 602).
* **k = 7:** 15 of 17 predicted, 0 false alarms, 0 value misses. The 2 misses are V_4 = 857 at t = 228 on
  (72,135) and t = 66 on (84,186).
* **Freezing.** Every prediction's `made_at_unix` is earlier than the last-write mtime of its held-out
  file. That is consistent with "frozen before", but an mtime is not proof of ordering.
* **What the k = 6 test measures.** The k = 6 C2 "rule" is a lookup in the computed plane. So stage 2 is
  a test of (A,B) invariance on 4 fresh curves, not a closed-form description of C2.

### README commands

These ran and reproduced the stored files:
* `ic_bounds_locus_curve.py`: identical output.
* `ic_bounds_locus_predict.py check` on (63,189): identical output.
* `ic_bounds_locus_plane.py --prime 193 --k 5 --fix B --values 0`: all 193 rows identical, apart from a
  new `certified` field.
* `ic_bounds_locus_lift.py`: same coefficients and reconstruction. The stored `c1-lift.json` is stale on
  point counts: it records 75 points at F_193, where a rerun gives 186, and 14 at F_257, where the claim
  says 19.

Not run:
* `ic_bounds_locus_report.py`, because it rewrites `summary.json` and the figures in place.
* The batch and plane full sweeps, because they are too long.

## Verdicts

| # | claim (short) | label | verdict |
|---|---|---|---|
| 1 | k=6: 1276 exact LFD in {4,5}; LFD 4 iff t=0,ab=0 or C2 | exact | confirmed (96-target stratified recompute incl. every deviating target; 1180 generic rows checked only for generator/Z/consistency) |
| 2 | generic values; 1 in V_4 only at (5,0,0), (3,0,0) | exact | confirmed (at k=7 "non-exceptional" includes the 8 V_4=857 targets no k=7 rule covers; counts are non-deviating by definition) |
| 3 | class values table k=6/k=7 | exact | confirmed |
| 4 | tangent: low_3 >= 3 proof | proof | confirmed (argument correct; V_3=6k-2 exact only at the computed k) |
| 5 | order-two square, t=0 b=0 zero set | proof | confirmed |
| 6 | H-orbit invariance | proof | confirmed |
| 7 | full (A,B) invariance 1268/1268 | empirical | confirmed (188/188 own, also at k=7) |
| 8 | complete plane classification p=193 | exact | confirmed on a sample of 795 points (the producer's dense loop stops only at a closed fixed point; unsampled generic points not recomputed) |
| 9 | C1 = cubic zero set minus 4 points over F_193; unique cubic at F_257, F_449 | exact | confirmed |
| 10 | C1 lifts to Q / a number field | conjecture | unverifiable (label correct; the reconstruction failure reproduces) |
| 11 | C2: 190 points, no degree <= 18 curve | exact | confirmed (all 190 points recomputed) |
| 12 | C2 is an absolutely irreducible curve of degree >= 19 | conjecture | unverifiable. The rank and the line counts are what 190 random points would also give, so they are no evidence that C2 is a curve at all. |
| 13 | "classes do not depend on the depth" | empirical | overstated: shown only for k=6 vs 7 (58 points at p=257, 15/15 recomputed). The k=5 data show no C1/C2 deficiency, C2 changes from saturated (LFD 4) at k=6 to LFD 5 at k=7 (new exact datum), and all values move with k. |
| 14 | k=6 prediction test 15/15, 0 FA, 14/15 values | exact | confirmed (the C2 part is a plane lookup, i.e. a test of claim 7; freeze order is supported only by mtimes) |
| 15 | k=7 prediction test 15/17, 0 FA | exact | confirmed |
| 16 | ordinary-ring re-check | exact | confirmed (own 16/16 cases) |
| 17 | simple features not exceptional | empirical | confirmed for what was checked: det = psi_4/psi_2 identity (sympy), and det=0, disc=0, psi_3(1)=0 loci mostly generic (190/192, 188/193, 190/193). The bounded feature search was not rerun. |

## Heavy runs (not run here)

The exact LFD at k = 7 for all 1909 targets with the verifier's dense engine would take about 97 s per
target at D = 5 (3473 columns, about 0.3 GB). That is about 51 CPU-hours single-threaded, which is why it
was not run. The command is the k = 7 line above with `--dmax 5` over `data/` and `heldout/`.
