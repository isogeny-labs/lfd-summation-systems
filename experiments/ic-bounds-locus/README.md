# Bounds program, item 2: the exceptional locus of the two-summand prime chain closure

System: PKM two-summand chain presentation, `{g_t} + chain`, with
`g_t = NF(S_3(x1, x2, t))` (`csd.ic_dreg.pkm_generator`), H the subgroup of
order 2^k of F_p^*, 2k chain variables. All dimensions below are **modulo the
chains** (squarefree chain monomials; equal to the ordinary closure for D >= 2,
re-checked independently in the ordinary ring, see "Validation").
`low_D = dim(V_D ∩ M_<=D-1)`; `|Z|` = scalar chain solutions.

Generic values (random targets): `dim V_3 = 2k+1`, `dim V_4 = 4k^2+172k-571`,
`low_4 = 4k+87` (k=6: 13, 605, 111; k=7: 15, 829, 115).

## Commands (every run used `ulimit -v 1500000`, `OMP/OPENBLAS/MKL_NUM_THREADS=1`)

```
# main sweeps: every rational target abscissa of every curve (one process per curve, < 5 min, ~0.2 GB)
EXTEND=6 scripts/ic_bounds_locus_batch.sh 193 6 5 experiments/ic-bounds-locus/data 0,1 1,0 22,23 111,126 46,29 68,122 61,121 0,5 5,0
scripts/ic_bounds_locus_batch.sh 257 7 4 experiments/ic-bounds-locus/data 0,1 1,0 0,3 3,0 88,32 194,158 200,59 87,256 148,227
# classification / feature search
python3 scripts/ic_bounds_locus_classify.py --prime 193 --k 6 --out classify-k6.json
python3 scripts/ic_bounds_locus_classify.py --prime 257 --k 7 --out classify-k7.json
# normalized plane (A, B) = (a/t^2, b/t^3) at t = 1
python3 scripts/ic_bounds_locus_plane.py --prime P --k K --fix B --values 0 --out plane/pP-kK-B0.jsonl          # lines B = 0
python3 scripts/ic_bounds_locus_plane.py --prime 193 --k 6 --fix A --values <A ≡ c mod 5> --no-certify --out plane/p193-k6-plane-c<c>.jsonl   # full plane, 5 chunks
python3 scripts/ic_bounds_locus_plane.py --prime 257 --k 6 --fix A --values 1,...,20 --no-certify --no-zeros --out plane/p257-k6-lines-A1-20.jsonl
python3 scripts/ic_bounds_locus_plane.py --prime 449|577 --k 6 --fix A --values 1,...,6 --no-certify --no-zeros --out plane/p449|577-k6-lines-A1-6.jsonl
python3 scripts/ic_bounds_locus_crossk.py --out plane/p257-crossk.json            # k = 6 vs k = 7 on the same (A, B)
python3 scripts/ic_bounds_locus_curve.py --plane 'plane/p193-k6-plane-c*.jsonl' --out c-fit-p193.json
python3 scripts/ic_bounds_locus_lift.py --out c1-lift.json                        # C1 cubic over F_193, F_257, F_449; lift to Q
# independent ordinary-ring re-check (lfd_engine, imports no repository closure engine)
python3 scripts/ic_bounds_locus_validate.py --data data/p193-k6-<curve>.jsonl --targets ... --degrees 3,4 --out validation/p193-k6-ordinary.json
# predictions on fresh curves (frozen before the curves were computed), then the check
python3 scripts/ic_bounds_locus_predict.py predict --rules rules-k7.json --prime 257 --k 7 --curve a,b --out predictions/...json
scripts/ic_bounds_locus_batch.sh 257 7 4 experiments/ic-bounds-locus/heldout 83,21 72,135 204,153 84,186 115,198 251,5
python3 scripts/ic_bounds_locus_predict.py check --pred predictions/...json --data heldout/...jsonl --out predictions/...-check.json
python3 scripts/ic_bounds_locus_report.py                                          # summary.json, plane-map.svg/png, verification.json
```

Curves (p = 193, k = 6, H = cubes of F_193^*, generator g0 = 125): (0,1), (0,5) [j = 0, b = 5 a non-cube],
(1,0), (5,0) [j = 1728, 5 not in H], random (22,23), (111,126), (46,29), (68,122), (61,121)
(`random.Random(20261107+p)`, first five draws). Held out: (63,189), (147,5) (draws 6, 7).
Curves (p = 257, k = 7, H = squares, g0 = 9): (0,1), (0,3), (1,0), (3,0), random (88,32), (194,158), (200,59),
(87,256), (148,227). Held out: (83,21), (72,135) (draws 6, 7) and (204,153), (84,186), (115,198), (251,5)
(`random.Random(20261108+p)`, first four draws). Every target t with x^3+ax+b a square (including 0) was used:
890 targets at k = 6 and 1153 at k = 7 (training); held out 386 targets (4 curves, k = 6) and 756 (6 curves, k = 7).
The extra fresh curves at k = 6 are (180,150), (172,33) (`random.Random(20261108+193)`, first two draws).

## Files

* `data/p<P>-k<K>-a<a>-b<b>.jsonl`: one row per target: generator (normal form, sha prefix), Z, per-degree
  `dim`, `low`, `one`, `codim`, evaluation rank, engine and exactness label, strict falls, status, exact LFD.
* `heldout/`: the same for the held-out curves; `predictions/`: frozen predictions and checks.
* `plane/`: normalized-plane sweeps; `c1-lift.json`, `c-fit-p193.json`: fitted exceptional curves.
* `validation/p193-k6-ordinary.json`: ordinary-ring re-check. `summary.json`, `plane-map.svg|png`,
  `verification.json` (sha256 of every script, engine and data file).


## Results

Evidence labels: **exact** = certified exact computation (certificate named), **proof** = argument given
here in full, **empirical** = pattern in data, **conjecture**. Nothing below is kernel-checked.

### 1. Per-target closure data (exact)

* k = 6, p = 193: 1276 targets (890 training, 386 held out). Every V_D for D <= 4 is a dense closure with
  `check_closed_exclusion_certificate` (or 1 derived); every V_5 is a shift closure that **reached the local
  target** |M_<=5| - |Z| with evaluation of M_<5 onto Z, so V_5 = I ∩ M_<=5 for all 1276 targets (no
  deficiency at D = 5 in any class), and the **exact LFD is certified for all 1276 targets** (1275
  local-to-global certificates at D = 5, 1 unit ideal with 1 in V_4). 39 targets re-run through
  `prime_local_profile` (D <= 4) agree.
* k = 7, p = 257: 1909 targets (1153 + 756), V_D through D = 4 (dense + exclusion certificate); 45
  re-runs agree. Exact LFD only for the unit ideal (3,0,0) (LFD 4); every other k = 7 target is
  undecided at D = 4 (V_4 deficient, as expected from growth.json). Exact LFD at k = 7 needs D = 5: see
  heavy runs.
* Independent ordinary-ring re-check (`lfd_engine.closure`, full polynomial ring, 12 variables, no
  repository closure engine): 13 targets, one per class plus generic, D = 3 and 4: all 26 ordinary
  dimensions and low dimensions equal C(2k+D, D) - |M_<=D| + quotient (`validation/p193-k6-ordinary.json`).

### 2. The classes (k = 6, p = 193 | k = 7, p = 257)

Normalized coordinates for t != 0: **A = a/t^2, B = b/t^3** (S_3 is weighted homogeneous: x, t weight 1,
a weight 2, b weight 3 in units of the scaling, S_3(mu x1, mu x2, mu t; mu^2 a, mu^3 b) = mu^4 S_3).

| class | condition | k = 6: V_3, V_4, low_4 | exact LFD (k = 6) | k = 7: V_3, V_4, low_4 | evidence |
|---|---|---|---|---|---|
| generic | none of the below | 13, 605, 111 | 5 | 15, 829, 115 | exact per target |
| order two | rhs(t) = t^3+at+b = 0, t != 0  (1+A+B = 0) | 13, 607, 113 | 5 | 15, 835, 117 | exact per target |
| tangent | 3t^2+a = 0, t != 0  (A = -3) | 34, 603, 111 (low_3 = 3) | 5 (falls 3,4,5) | 40, 828, 117 | exact; low_3 >= 3 proof below |
| C1 | cubic C1(A, B) = 0 (below) | 13, 602, 111 | 5 | 15, 826, 115 | exact per target; cubic exact at p = 193 |
| C2 | (A, B) in the 190-point set C2 (p = 193) | 13, = |M_<=4| - |Z| (saturated: V_4 = I ∩ M_<=4), low_4 = 299 - |Z| | **4** (single fall at 4) | 15, 857, 117 | exact per target; no curve of degree <= 18 |
| t = 0, ab != 0 | t = 0 | 13, 612, 117 | 5 | 15, 849, 123 | exact |
| t = 0, a = 0 (j = 0) | t = 0, a = 0 (also tangent) | 44, 773, 278 (low_3 = 4) | 4 (falls 3, 4) | 52, 1052, 150 | exact |
| t = 0, b = 0, a in H (j = 1728) | (0,0) target, a^(2^k) = 1 | 13, 730, 235, |Z| = 64 | 4 | 15, 1343, 350, |Z| = 128 | exact |
| t = 0, b = 0, a not in H | (0,0) target, a^(2^k) != 1 | unit ideal, 1 in V_4 | 4 | unit ideal, 1 in V_4 | exact (matches the kernel-checked j = 1728 family) |
| special points | (A,B) = (-3, 2) cusp (1+A+B = 0 = A+3, singular cubic); C1 ∩ {A = -3} = (-3, 170); (162, 58) on C1; (32, 103) on C1 | (216, 681, 216); (34, 600, 111); (13, 599, 111); (13, 607, 113) | — | — | exact (p = 193 plane) |

**Exact LFD at k = 6, p = 193 (all 1276 targets, 13 curves, exact):** LFD = 4 exactly when t = 0 and ab = 0,
or the normalized point (A, B) lies in C2; LFD = 5 for every other target. 12 of 1276 targets have LFD 4
(t = 0 with ab = 0: 3; C2: 9). Strict falls: (4,5) generic; (3,4,5) tangent; (3,4) for t = 0, a = 0; (4) for LFD 4.

Counts in the training data, k = 6: generic 861, order two 10, tangent 2, C1 3, C2 8, t = 0 (ab != 0) 3,
t = 0 j = 0 1, t = 0 j = 1728 2 (29 deviating of 890). k = 7: generic 1118, order two 12, tangent 4, C1 5,
C2 6, t = 0 (ab != 0) 5, t = 0 j = 0 1, t = 0 j = 1728 2 (35 deviating of 1153). The (0,5) j = 0 twist has no
deviating target at all (no rational t = 0, no rational 2-torsion, no rational solution of 3t^2 = 0 other than t = 0).

Conditions that were tested and are **not** exceptional by themselves (empirical, all p = 193 targets and the
full plane): t in H, -t in H, t^2 in H, 1/t in H, a or b in H, membership of every coefficient ratio of g in H,
t - e in H for a 2-torsion abscissa e, x(2R) in H, the j-invariant (j = 0 and j = 1728 curves behave
generically away from the classes above), the discriminant 4A^3+27B^2 = 0 away from the cusp (188 generic
points), R of order 3 (psi_3 = 0), R of order 4 (the determinant of the 3x3 coefficient matrix of g equals
-2(A^3+5A^2+4AB-5A+8B^2-20B-1) = psi_4/psi_2 up to a constant; its points are generic), and the order of the
point with abscissa 1 on the normalized curve or its twist (no small-order pattern on C1 or C2).
The bounded feature search (`classify-k6.json`, `classify-k7.json`: unions of <= 4 conjunctions of <= 2
literals over 40 polynomial and H-membership atoms) describes the V_3 deviations exactly (3t^2+a = 0) but
finds no description of the full deviating set; C1 and C2 are needed.

### 3. The normalized plane (A, B)

* **H-scaling invariance (proof).** For mu in H, x -> x/mu is an F_p-algebra automorphism of
  F_p[x1, x2]/(x1^(2^k)-1, x2^(2^k)-1) that multiplies each chain monomial X^(e1,e2) (0 <= e_i < 2^k) by
  mu^-(e1+e2), so it preserves Hamming weights, maps V_D onto V_D and preserves every dimension. By weighted
  homogeneity, NF(S_3(x1, x2, mu t; mu^2 a, mu^3 b)) = mu^4 NF(S_3(x1/mu, x2/mu, t; a, b)). Hence (a, b, t)
  and (mu^2 a, mu^3 b, mu t) have identical closure dimensions for mu in H.
* **Full invariance (empirical).** For mu not in H the substitution is not an algebra map (it rescales the
  wrap x^(2^k) = 1), yet every k = 6 target with t != 0 has the same class as the t = 1 point
  (A, B) = (a/t^2, b/t^3): V_3, low_3, V_4, low_4 equal, or both saturated (V_4 + |Z| = 794): **1268 of 1268**
  agree (434 with t in H, 834 with t not in H; `summary.json` "AB_invariance_k6"); 12 further spot checks on
  three sheets t in {1, 2, original} (scratch). Saturated values differ between sheets only through |Z|.
* **Complete plane, p = 193, k = 6 (exact; dense closures, the redundant exclusion re-check skipped, see
  below):** all 37249 points (A, B) at t = 1: generic 36486; 1+A+B = 0: 192 (+ cusp); A = -3: 191 (+ cusp
  + 1); C1: 186; C2: 190; special points 3 (+ (32,103)). Figure: `plane-map.svg`/`.png`.
* **C1 (exact at p = 193).** The 602-class is exactly the zero set over F_193 of the cubic
  `A^3 + 30A^2B + 72AB^2 + 50B^3 + 96A^2 + 131AB + 8B^2 + 166A + 162B + 15` minus its 4 special points
  (190 zeros: 186 of class 602, (-3,2), (-3,170), (162,58), (32,103)); no generic point lies on it.
  The same class is a cubic over F_257 (19 points, unique cubic, coefficients in `c1-lift.json`) and over
  F_449 (10 points, unique). Rational reconstruction from 193·257·449 = 22,270,849 gives no small-height
  rational cubic: the lift to Q (or to a number field) is open.
* **C2 (exact at p = 193).** 190 points, V_4 + |Z| = 794 on every one; no nonzero polynomial of degree
  <= 18 vanishes on all of them (the evaluation matrix of the 190 monomials of degree <= 18 has full rank
  190). The number of C2 points per line A = const is 0 on 77 lines, 1 on 67, 2 on 32, 3 on 11, 4 on 5,
  6 on 1 (roughly Poisson(1), as for the rational points of a high-degree curve). Conjecture: C2 is an
  absolutely irreducible plane curve of degree >= 19.
* **k-independence (empirical).** At p = 257, all 58 normalized points that are exceptional at k = 6 or
  k = 7 were computed at both depths: order two (24 points) 607/835, C1 (16) 602/826, C2 (15)
  saturated/857, tangent (3) 603/828, with no point changing class (`plane/p257-crossk.json`). On the line
  B = 0 the same lines and the class-C1/C2 points move with p (p = 193: C2 at A = 20, 162; p = 257: C2 at
  14; p = 449: C2 at 184, 221, C1 at 342, 412, 418), while 1+A+B = 0 and A = -3 persist.
* **k = 5 contrast (exact, B = 0 lines, p = 193 and 257).** At k = 5 every point of the line is either a
  unit ideal with 1 in V_4 or saturated (V_4 = I ∩ M_<=4); only V_3 deviates (28 = 6k-2 on A = -3).
  This is the depth-6 step seen from the locus side: below k = 6 there is no deficiency for C1/C2 to
  change.

### 4. Proofs for the polynomial classes

* **Tangent (proof that low_3 >= 3 when 3t^2+a = 0, t != 0, k >= 3).** The weight-2 part of g is
  u11u21 - 2t(u11u20 + u10u21) - 2(t^2+a)u10u20 = (u11, u10) M (u21, u20)^T with
  M = [[1, -2t], [-2t, -2(t^2+a)]], det M = -2(3t^2+a). So it factors as l1 l2 with
  l_i = u_{i,1} - 2t u_{i,0} iff 3t^2+a = 0. Then h1 = (u11 + 2t u10) g = u11 g + 2t u10 g is in V_3
  (two products of g by a variable) and, using u10^2 = u11, u11^2 = u12,
  (u11 + 2t u10) l1 = u12 - 4t^2 u11 has weight 1, so h1 = (u12 - 4t^2 u11) l2 + (u11 + 2t u10)(g - l1 l2)
  has weight <= 2 (g - l1 l2 has weight <= 1). The monomial u12 u21 has coefficient 1 in h1 and does not
  occur in g or in h2 (the block-swapped element, which contains u22 u11 with coefficient 1, absent from g
  and h1); hence g, h1, h2 are linearly independent and dim(V_3 ∩ M_<=2) >= 3 > 1 = dim V_2: a strict fall
  at 3. (The exact values 3 and V_3 = 6k-2 at k = 5, 6, 7 are exact computations, not proved in general.)
* **Order two (proof of the structure).** If t^3+at+b = 0, then S_3(x1, x2, t) = (x1x2 - t(x1+x2) - a - 2t^2)^2
  identically (substitute b = -t^3 - at; checked symbolically), so g = NF(q^2) with
  q = u10u20 - t(u10+u20) - (a+2t^2). The chain algebra is reduced (x^(2^k)-1 is separable over F_p), so
  q is in the ideal although it is not in span(g); the closure picks up 2 extra dimensions at D = 4
  (exact), not at D = 3.
* **t = 0.** g = u11u21 - 2a u10u20 - 4b(u10+u20) + a^2: the coefficients of u11, u21, u11u20, u10u21 vanish.
  With a = 0 the top form u11u21 is rank one (the tangent case at t = 0). With b = 0 and a in H the solution
  set is the anti-diagonal x1 x2 = a (|Z| = 2^k); with a not in H it is empty (unit ideal).

### 5. Predictions on fresh curves (frozen before the curves were computed)

* **k = 7** (`rules-k7.json`, sha256 ac439458...; polynomial rules plus the cubic C1 over F_257 fitted on
  training data): 6 fresh curves, 756 targets, 17 deviating. Predicted 15, all correct with the predicted
  values (V_3 = 40, V_4 = 835, 849, 826), **0 false alarms**; 2 misses, both class C2 (V_4 = 857, t = 228 on
  (72,135) and t = 66 on (84,186)), for which no k = 7 rule existed. 4 of 6 curves exact.
* **k = 6, stage 1** (`rules-k6.json`, polynomial rules only, frozen before the plane sweep):
  (63,189): predicted {0, 65, 81, 128, 132, 173}, (147,5): {12, 68}; all predicted targets deviate with the
  predicted values; one miss (t = 89 on (63,189), a C1 point).
* **k = 6, stage 2** (`rules-k6-plane.json`, sha256 9c524feb...: polynomial rules + cubic C1 over F_193 +
  C2 by lookup of (a/t^2, b/t^3) in the complete plane, frozen after the plane and before the curves):
  4 fresh curves, 386 targets, 15 deviating: **all 15 predicted, 0 false alarms, 4 of 4 curves exact**,
  including the C2 target t = 90 on (180,150) with predicted and certified exact LFD 4. 14 of 15 predicted
  values matched; the exception is t = 89 on (63,189), whose normalized point is the special point
  (32,103) of C1 (V_4 = 607, as the plane itself records; the rule's single value 602 was wrong there).

## What is not established

* No statement here is kernel-checked. The general-k formulas (V_3 = 6k-2 on the tangent line, the
  increments of the other classes) are data at k = 5, 6, 7 only.
* C1 over Q and C2 as a curve are not determined; C2 is only known as a point set at p = 193 (and by
  sampled lines at p = 257, 449, 577).
* Exact LFD at k = 7 is not computed (needs D = 5 at 3473 columns; heavy run below). The k = 6 classes
  that have LFD 4 (C2, t = 0 with ab = 0) are the only k = 6 targets below the generic 5.
* The full (A, B) invariance is empirical (1268/1268 targets); only the H-orbit invariance is proved.
* The plane sweep skipped the redundant exclusion re-check (`--no-certify`): its dimensions are the dense
  engine's fixed points; the B = 0 lines and all 3185 per-target rows used the certificates, and every
  per-target class agrees with the plane.

## Run notes

The first launch of the k = 6 sweep (dense D = 5, ~0.6 GB per process) coincided with a cgroup OOM kill of
another process in the shared 14.3 GB memory cgroup (an unrelated long computation at 13.5 GB RSS);
the sweep was restarted with the shift engine at D = 5 (peak RSS ~0.18 GB per process) and `ulimit -v 1500000`.

