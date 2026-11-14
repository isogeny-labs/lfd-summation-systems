# Towards tight LFD bounds for ECDLP summation systems: the six-item programme

Date: 2026-11-08. This document synthesises six item reports and their independent
verifications:

| item | directory | producer report | verifier report |
|---|---|---|---|
| 1. Three or more summands (m = 3) | `experiments/ic-bounds-m3/` | `README.md` | `VERIFY.md` |
| 2. Exceptional locus (m = 2, prime) | `experiments/ic-bounds-locus/` | `README.md` | `VERIFY.md` |
| 3. Degree-5 closure growth d_5(k) | `experiments/ic-bounds-d5/` | `README.md` | `VERIFY.md` |
| 4. Exact value at depth 10 | `experiments/ic-bounds-exact10/` | `README.md` | `VERIFY.md` |
| 5. Fall-space structure | `experiments/ic-bounds-fall/` | `README.md` | `VERIFY.md` |
| 6. Binary and Koblitz at larger size | `experiments/ic-bounds-binary/` | `README.md` | `VERIFY.md` |

Notation follows `formal/ICFallDegrees/Closure.lean` and `paper/curve-lfd-bounds.tex`:
V_D is the degree-D closure, LFD the ordinary iterated last fall degree, Z the scalar
solution set, M_<=D the monomials (chain masks) of degree at most D, I_D = I ∩ M_<=D.
All prime computations work modulo the chains, which equals the ordinary closure for
D >= 2 (repository premise; item 3's verifier re-checked it numerically at k = 4, 5 and
items 1, 2, 4 cross-checked selected cases in the full ordinary ring with
`experiments/ic-curve-validation/code/lfd_engine.py`).

Evidence labels: **kernel** (Lean-checked), **proof** (written argument, checked by a
second agent where stated), **exact** (certified computation; certificate named),
**empirical** (pattern in data), **conjecture**. A closure stall is never used as a bound.
Nothing produced by the six items is kernel-checked.

Kernel-checked results that the programme builds on (unchanged by this work):

* K1. LFD <= mk+1 for prime chain systems (chains plus reduced inputs, every m).
* K2. LFD <= r+1 for binary and Koblitz systems.
* K3. LFD >= 4 for every unit-ideal two-summand prime system with 3t^2+a != 0, every depth.
* K4. Exact LFD 4 for the (0,0)-target j = 1728 family when a^(2^k) != 1.

---

## 1. Every claim, with its label after verification

"Verifier" names the check in the item's `VERIFY.md`. Where the label or wording changed,
the change is in **bold**. Counts of independent recomputation are given where the
verifier sampled.

### Item 1: three summands (m = 3), prime subgroup chains, k = 1..4

| # | claim | producer | after verification |
|---|---|---|---|
| 1.1 | Direct encoding {chains, NF(S_4)}: LFD histogram k=1 {4:100, 3:5, 0:6}; k=2 {7:105, 6:11}; k=3 {7:117}; k=4 {7:117}. Certificates: unit (1 in V_D) or local-to-global. | exact | exact. Verifier recomputed all cases at k = 1..3 and 26/117 at k = 4: 0 mismatches. |
| 1.2 | At k = 2 the direct LFD equals the kernel bound mk+1 = 7 in all 105 cases with deg NF(S_4) = 6, so K1 is attained for m = 3, k = 2. | exact | exact (105/105 recomputed). |
| 1.3 | Carry lemma (Prop. 1): chains + {g}, g of degree d >= 2, some chain variable u with deg NF(ug) <= d, g not supported in one level set => LFD >= d+1. | proof | proof (checked step by step). |
| 1.4 | Corollary 2: direct m = 3, every k >= 3, every p with 2^k \| p-1, every curve, every target with 5t^3+3at+6b != 0: LFD >= 7. Uses c333 = -8(5t^3+3at+6b) and the x3^4 coefficient S_3(x1,x2,t)^2. | proof | proof (identities re-derived independently; 27 off-grid instances at k = 3, 4 all give exact LFD 7). |
| 1.5 | On 5t^3+3at+6b = 0 (deg NF(S_4) = 5, includes the (0,0) target on j = 1728): LFD 6 at k=2 (11), 7 at k=3 (10), 7 at k=4 (11). | exact | exact (all 32 recomputed). |
| 1.6 | Direct LFD = deg NF(S_4) + 1 on every deg-6 case at k = 2..4 and every deg-3 case at k = 1. | empirical | empirical. **The count is 318 deg-6 cases, not 316.** The remark that the degree-7 certificate "must stop working at some depth" is a heuristic, not a result. |
| 1.7 | Split encoding with a free auxiliary z (I_free): k=1 {5:83, 4:28}; k=2 {5:116}; k=3 {6:111, 5:6} (the six are y_R = 0 targets); k=4 {6:117}. | exact | exact. Recomputed k = 1, 2 in full, 29/117 at k = 3, **5/117 at k = 4**. |
| 1.8 | Field-equation split I_field = I_free + (z^p - z): closures equal for D < p; LFD equals the free value when dim R/I_free = dim R/I_field, else LFD >= p. Counts 100/87/80/69 and 11/29/37/48 at k = 1..4. | proof | proof (Remark 3 checked); counts exact. |
| 1.9 | m = 2 reference: LFD 3 at k = 1, 4 at k = 2..4. Gaps direct − m2 = +3 at k = 2..4; split − m2 = +1 at k = 2 and +2 at k = 3, 4. | exact | **overstated, corrected:** split − m2 at k = 3 is +2 in 111 cases and **+1 in 6** (the y_R = 0 targets). Everything else exact. |
| 1.10 | Ordinary-ring validator agrees: direct 140/140, split 122/122 (k = 1, 2; split only where I_free = I_field). | exact | exact (rerun byte-identical; verifier's own inputs 38/38). |
| 1.11 | S_4 = Res_X(S_3(x1,x2,X), S_3(x3,t,X)) is correct; every rational decomposition with x_i in H lies in Z; batched engine matches dense engines. | exact | exact (independent resultant; 8331 brute-force decomposition triples all in Z). |

### Item 2: exceptional locus, two summands, k = 6 (p = 193) and k = 7 (p = 257)

| # | claim | producer | after verification |
|---|---|---|---|
| 2.1 | k = 6: all 1276 rational targets on 13 curves have exact LFD in {4, 5}; LFD 4 exactly for t = 0 with ab = 0 (3) or (A,B) = (a/t^2, b/t^3) in C2 (9). V_5 = I_5 for all. | exact | exact. Verifier recomputed 96 targets (all 44 deviating + 52 generic) with its own certificate; the other 1180 generic rows were checked for generator, \|Z\| and internal consistency only. |
| 2.2 | Generic values k=6: V_3 = 13, V_4 = 605, low_4 = 111; k=7: 15, 829, 115. 1 in V_4 only for t = 0, b = 0, a not in H. | exact | exact. At k = 7 "non-exceptional" means "non-deviating"; the 8 V_4 = 857 targets are covered by no k = 7 rule. |
| 2.3 | Class values at k = 6 \| k = 7 (order two 607\|835; tangent V_3 = 6k−2, V_4 603\|828; t = 0 subcases; C1 602\|826; C2 saturated\|857). | exact | exact. The k = 7 value 1343 for t = 0, b = 0, a in H is deficient (eval rank 120 < \|Z\| = 128), so it is not I_4. |
| 2.4 | Tangent: 3t^2+a = 0, t != 0, k >= 3 => dim(V_3 ∩ M_<=2) >= 3 > dim V_2, a strict fall at 3. | proof | proof. V_3 = 6k−2 is exact only at k = 5, 6, 7. |
| 2.5 | Order two (t^3+at+b = 0): S_3(x1,x2,t) = (x1x2 − t(x1+x2) − a − 2t^2)^2; for t = 0, b = 0 the zero set is x1x2 = a. | proof | proof. |
| 2.6 | H-orbit invariance: (a,b,t) and (mu^2 a, mu^3 b, mu t) give identical closure dimensions for mu in H. | proof | proof. |
| 2.7 | Full invariance: the class depends only on (A,B) for t != 0 (1268/1268 at k = 6). | empirical | empirical (188/188 verifier-side, also at k = 7). Only the H-orbit part is proved. |
| 2.8 | Complete normalized plane at p = 193, k = 6 (37249 points) with the class counts. | exact | exact on the 795 recomputed points; the producer's loop stops only at closed fixed points, but the 36k unsampled generic points were not recomputed. |
| 2.9 | Over F_193 the 602-class is the zero set of an explicit cubic minus 4 special points; unique cubic over F_257 (19 points) and F_449 (10 points). | exact | exact. `c1-lift.json` has stale point counts (75, 14); a rerun gives the same coefficients. |
| 2.10 | C1 is the reduction of one curve over Q or a number field. | conjecture | conjecture (small-height CRT lift over 193·257·449 failed). |
| 2.11 | C2 at p = 193: 190 points, V_4 + \|Z\| = 794, no nonzero polynomial of degree <= 18 vanishes on them. | exact | exact (all 190 recomputed). |
| 2.12 | C2 is an absolutely irreducible plane curve of degree >= 19. | conjecture | conjecture, **essentially unsupported**: full rank at degree 18 and Poisson-like line counts are what 190 random points give. There is no evidence that C2 is a curve. |
| 2.13 | The classes do not depend on the depth (58 points at p = 257 keep their class k = 6 -> 7). | empirical | **overstated; restricted to "the classes at the 58 points coincide at k = 6 and k = 7".** k = 5 shows no C1/C2 deficiency, all values move with k, and the LFD is not depth-independent (row V2.1). |
| V2.1 | **New (verifier):** the C2 target t = 228 on (72,135), p = 257, k = 7, \|Z\| = 84 has V_5 = 3389 = I_5 with M_<=4 onto Z: exact LFD 5, whereas C2 has LFD 4 at k = 6. | — | exact (local-to-global), one target. |
| 2.14 | k = 6 held-out prediction: 15/15 deviations predicted, 0 false alarms, 14/15 values. | exact | exact. The C2 part is a lookup in the computed plane (a test of 2.7, not a closed-form rule); "frozen before" rests on file mtimes only. |
| 2.15 | k = 7 held-out prediction: 15/17 predicted with correct values, 0 false alarms; the 2 misses are C2. | exact | exact. |
| 2.16 | Ordinary-ring check at D = 3, 4 matches the chain closure (26 cases). | exact | exact (verifier 16/16). |
| 2.17 | Simple features (t, a, b, ratios in H; j; discriminant; torsion 3, 4) are not exceptional by themselves; the 3x3 coefficient determinant is psi_4/psi_2. | empirical | empirical (feature search not rerun; the determinant identity checked). |

### Item 3: degree-5 closure growth

| # | claim | producer | after verification |
|---|---|---|---|
| 3.1 | Split engine `csd/ic_bounds_d5.py` agrees with the dense engine on 18 closures at k = 6, 7 and with the shift engine on 12. | exact | exact. |
| 3.2 | All 20 depth-8 systems reproduce dim V_4 = 1061 and the saved D = 5 certificates. | exact | exact (22 systems incl. 2 at k = 7). |
| 3.3 | k = 10, three instances: dim V_5 = 18,569, low 3163, profile [0,0,2,138,3163,18569], deficient-closed => LFD >= 6. | exact | exact (full deterministic sweep; resume hash reproduced by the verifier because the producer's resume record is missing). |
| 3.4 | d_5(11) = 25,033 on two instances (p = 12289; unit seed 10; proper seed 12 with \|Z\| = 348), deficient-closed => LFD >= 6. | exact | exact (independent engine, both instances). |
| 3.5 | Split-closure certificate (sigma-split, block-1 shifts, final sweep) certifies W = V_D. | proof | proof. |
| 3.6 | d_5 <= 1 + 2k·lambda, lambda = dim(V_5 ∩ M_<5). | proof | proof. **Vacuous at the measured depths** (63,261 and 75,703 exceed \|M_<=5\| = 21,700 and 35,443). |
| 3.7 | \|Z\| <= 2^(k+1) for every two-summand subgroup system; hence any k = 11 instance with dim V_5 = 25,033 has LFD >= 6. | proof | proof. |
| 3.8 | dim V_4 = 4k^2+172k−571 and dim(V_4 ∩ M_<=3) = 4k+87 at k = 6..13 (two instances per depth). | exact values; formula empirical | same. For k >= 8 the formula is now explained per instance by item 5 (row 5.3). |
| 3.9 | d_5 and the profiles are identical on every instance at k = 10 (3) and k = 11 (2). | empirical | empirical. |
| 3.10 | Missing top-weight leading terms 0, 0, 98, 4742 at k = 8..11; k = 10 may be a transition point. | empirical | empirical. |
| 3.11 | Cost model: measured work ratios, storage and throughput; projections for k = 12, 13, 14. | empirical | measured part empirical; **projections downgraded to "model output, untested"**. |
| 3.12 | Four values d_5(10..13) only fix a cubic; d_5 decides LFD >= 6 per instance, nothing about upper bounds, genericity or 6 -> 7; a step every four depths. | conjecture | methodology correct; the step-every-four-depths is a conjecture. |

### Item 4: exact value at depth 10 (preparation)

| # | claim | producer | after verification |
|---|---|---|---|
| 4.1 | Block engine `csd/ic_chain_block_closure.py` agrees with dense on "1 in V_D" in 72/72 rows (k = 5..7, D = 4..6). | exact | exact (versions 865032b8 and de5e3970). |
| 4.2 | In the 45 rows with 1 not in V_D, the block engine's stall rank equals dense dim V_D. | empirical | empirical. A block-engine stall is not a certificate. |
| 4.3 | New split engine `csd/ic_bounds_exact10.py` matches dense in 72/72 rows (scratch and from V_{D−1}); its "reached" and "closed-deficient" statuses are certificates. | exact | exact. |
| 4.4 | Negative control rejected in 56/56 rows. | exact | exact, **weak control**: the control space is far smaller than V_D, so sensitivity to a small deficit is untested. |
| 4.5 | split_profile reproduces exact LFD and falls on 24/24 (instance, k) pairs. | exact | exact. |
| 4.6 | Mid-scale k = 8 (3 systems) and k = 9 (1 system): exact LFD 5. | exact | exact. |
| 4.7 | Runner at k = 8, D = 6 (p = 257, x = 190): dim V_6 = 14641 = \|M_<=6\| − \|Z\|. | exact | exact (value; interrupt-and-resume only read from logs). |
| 4.8 | Soundness of the split certificates. | proof | proof. |
| 4.9 | LFD >= 7 at depth k iff f_6(k) < delta_5(k) (given ev(M_<5) onto Z); delta_5(10) = 3131 on the unit instance. | proof | proof. |
| 4.10 | The original depth-10 degree-6 runner needed > 15 GB at set-up; four OOM kills near 22:34. | exact | shape arithmetic exact; attribution of the kills unverifiable. **README stale:** pid 11870 was itself OOM-killed at about 23:31 (anon-rss 8.85 GB). |
| 4.11 | The three k = 10 instances (x = 9362, 0 points; x = 7412, 98; x = 7776, 60). | exact | exact. |
| 4.12 | Dense split engine at D = 6: 6.6/17/45/113/520 GB for k = 10..14; time brackets. | empirical | model; memory arithmetic reproduces; **the cols^2 time bracket does not reproduce** (least squares gives 2.2/7.5/22/61/153 h, not 1.8/6/18/49/123 h); all times are extrapolated 4–34x beyond calibration. |
| 4.13 | k_7 (first depth with LFD >= 7, m = 2) lies in 12..17. | conjecture | conjecture. |
| 4.14 | (main session, after this table) Exact LFD 6 at k = 10 on the unit instance: V_5 closed without 1 (bit-identical rerun), 1 ∈ V_6 by the block engine (both swap parts full rank). So k_7 >= 11. | exact | **independently confirmed (B7, 2026-11-08):** `csd/ic_bounds_d5_verify.py` from the generator, 453 rounds, V_6 = M_<=6 (30318 + 30142), 1 ∈ V_6, full deterministic sweep of 21700 low rows passed (`experiments/ic-bounds-d5/verify/k10-p12289-s10-D6.json`). |

### Item 5: fall-space structure (degree 4, two summands)

| # | claim | producer | after verification |
|---|---|---|---|
| 5.1 | p = 40961, y^2 = x^3+21223x+9887, target (4748,524): dim V_4(k) = 605, 829, 1061, 1301, 1549, 1805, 2069 at k = 6..12; low_3 = 4k+87, low_2 = 2; 1 not in V_4; every closure verified closed. | exact | exact (independent closure; hashes match). |
| 5.2 | V_4(k) = W_k := ⊕_h tau_h(V_4(7) ∩ M_<=4−wt h) (h on levels >= 7, wt h <= 2) for k = 7..12; W_k passes the exclusion certificate built without V_4(k); V_4(k) ∩ M_{levels<=6} = V_4(7). | exact | exact. |
| 5.3 | Reduction theorem: W_k contains the generator rows, 1 in W_k iff 1 in V_4(7), and W_8 closed => W_k closed for every k >= 8. Hence 1 not in V_4(k) and dim V_4(k) <= 4k^2+172k−571 for all k >= 8 when the depth-8 certificate holds. | proof | proof. "Instance-free" applies to the reduction step only; the dimension statement uses this instance's base profile (2, 115, 829). |
| 5.4 | The base cannot be depth 6 (42 wrap-around dimensions; Free(V_4(6)) has the right dimension but is not closed). | exact | exact (k = 7..10). |
| 5.5 | F_4(k) = V_4 ∩ M_<=3 is a direct sum of level-set pieces with the stated counts; F_4(k) = F_4(7) ⊕ free translates of g and w for k = 8..10. | exact | exact. |
| 5.6 | V_4 ∩ M_<=2 = span(g, w), w = Q_w g with Q_w in closed form over Q(a,b,t), denominator a+3t^2. | proof | proof (generic over Q(a,b,t); the specialization checked on one instance). |
| 5.7 | Two p = 12289 instances (random curve; j = 1728 probe target x = 9362) have the same dims at k = 6..9, the same V_4(7) pivot set, and pass the free-extension certificate at k = 8, 9. | exact | exact. |
| 5.8 | Multiplier filtration on unit instance p = 1073153 (translation weight up to 10; kappa_4 matches for w <= 7). | exact | exact. ("10 = 2k−2" holds only at k = 6.) |
| 5.9 | Degree 5 on the primary instance: V_5 = I_5 at k = 6, 7, 8 (LFD <= 5). | exact | exact. **Round counts are 8, 10, 12** under the convention used for degree 4, not 9, 11, 13. |
| 5.10 | Degree-4 closure needs 9 rounds at every k = 6..12; birth histograms identical except two free families. | exact | exact. The README one-liner to reproduce it fails (JSON without trailing newline). |
| 5.11 | V_5(k) is the free extension of a wrap-free base V_5(k_0), k_0 >= 11, giving d_5(k) as a quartic. | conjecture | conjecture. Its premise "1 in V_5(k) for k <= 9" refers to the unit-ideal probe family, not to the primary instance. |
| 5.12 | No closed form found for the 829-dimensional base V_4(7); its pivot set is the same on three instances. | empirical | empirical. (README Section 2 calls 40817 prime; it is 7^4·17. Harmless.) |

### Item 6: binary and Koblitz

| # | claim | producer | after verification |
|---|---|---|---|
| 6.1 | All 1871 cases (1726 distinct systems) have an exact LFD: 1355 proper (local-to-global), 516 unit; none capped. | exact | exact (verifier recertified all 1871 with independent code). |
| 6.2 | m = 2, n = 11..14: max LFD 4 at k = 8, 4 at k = 9, 5 at k = 10; below the written min(k,l)+2 = 10, 11, 12 and K2 (17, 19, 21). | exact | exact. |
| 6.3 | First certified binary m = 2 systems with LFD 5 in the repository (30 systems at k = 10). | exact | **refuted as a priority claim**: `experiments/ic-binary-inverse` already has a certified full-field n = k = 8 Koblitz system with LFD 5. Restated: first certified LFD 5 for subspace bases with k < n. **28 distinct systems (30 cases).** |
| 6.4 | At k = 10, LFD 5 for \|Z\| in [472, 1056], LFD <= 4 for \|Z\| in [46, 352]; no separation at k = 8, 9. | empirical | empirical; on distinct systems the counts are **28 vs 127** (30 vs 137 are cases). Grid at n = 11, k = 10 is half done. |
| 6.5 | Mean \|Z\| tracks 2^(2k−n); unordered decompositions ≈ 0.45\|Z\|. | empirical | empirical. |
| 6.6 | Direct m = 3 S_4 descent: max LFD 6 at k = 2, 7 at k = 3, 4, 5 (695 cases); exceeds the m = 2 formula k+2. | exact | values exact; **"exceeds k+2" overstated: at k = 5, 7 = k+2.** It exceeds k+2 at k = 2, 3, 4. |
| 6.7 | The direct m = 3 maximum stays 7 from k = 3 to k = 5 (k = 5: 8 systems, one curve, one basis). | empirical | empirical. |
| 6.8 | Split m = 3: max LFD 4, 5, 5 at k = 2, 3, 4 (608 cases); paired order-two 54 >, 103 =, 4 < of 161. | exact | exact; k = 4 coverage partial. |
| 6.9 | P1: D* >= min{D : sum_{i<D} C(r,i) >= \|Z\|}. | proof | proof. |
| 6.10 | P2: Z is read off V_D* with r multiplication matrices of size \|Z\|. | proof | proof (the basis can be chosen without knowing Z). |
| 6.11 | P3: binary S_3, S_4 depend only on b; Koblitz (0,1), (1,1) give identical systems; order-two Z = {x1x2 = sqrt b}. | proof | proof; the explanation of the 145 duplicate cases is incomplete (other duplicates exist). |
| 6.12 | P4: invariant subspaces of F_{2^11}, F_{2^13}. | proof | proof + exact. |
| 6.13 | P5: D* <= d_full (static Macaulay degree); d_full = 12 = r at k = 4 for m3-direct. | exact | exact (23-row recheck). |
| 6.14 | D* − LFD is 0 in 1817 rows and 1 in 54; largest certificate matrices 21,700 and 16,384 columns. | exact | exact. |
| 6.15 | HKYY degree max(LFD, \|Z\|) >= r in 1041/1355 proper systems. | exact | exact (cases, not distinct systems). |
| 6.16 | Independent validation: rehash 1871, engine A 20/20, MutantXL 17/17. | exact | exact (3 of the 11 "m3" MutantXL rows are m3-split). |

---

## 2. Per item: established, failed, still needed

### Item 1, three summands

**Established.** Exact LFD for every case of a 461-case grid (5 primes, 6 curves, k = 1..4)
in three encodings. The direct encoding has LFD 7 at k = 3, 4 on every case, and K1's
bound mk+1 = 7 is attained at k = 2. A general lower bound is proved: LFD >= 7 for the
direct m = 3 system at every k >= 3 off the locus 5t^3+3at+6b = 0 (Corollary 2, via the
carry lemma). The free-auxiliary split gives LFD 6 at k = 3, 4 (5 for y_R = 0 targets);
adding z^p − z gives either the free value or LFD >= p.

**Failed or corrected.** Two summary statements were imprecise (318 not 316 cases;
split − m2 is +1 for six k = 3 targets). The degree-7 heuristic ("must stop working")
is not a result.

**Still needed.** An upper bound LFD <= 7 for the direct system at general k (or the depth
where it rises); an explanation of the exceptional locus' step 6 -> 7 at k = 3 (the carry
lemma does not apply there); a lower bound for the split; k = 5 data (heavy runs H9, H10);
a cost comparison of direct and split by Macaulay matrix size, not LFD.

### Item 2, exceptional locus

**Established.** At k = 6, p = 193, the exact LFD of every rational target on 13 curves is
4 or 5, and the deviations from the generic closure profile are fully classified at that
(p, k): order two, tangent (3t^2+a = 0), four t = 0 subcases, C1 (a cubic over F_193 minus
4 points), C2 (190 points), and 2 special points. Three facts are proved: the tangent fall at
3, the square identity for order-two targets, and invariance under H-scaling. Held-out
prediction tests at k = 6 and k = 7 had 0 false alarms.

**Failed or corrected.** "Classes do not depend on the depth" is overstated; the verifier
found exact LFD 5 for a C2 target at k = 7, where C2 had LFD 4 at k = 6. C1 did not lift to Q
at small height. The claim that C2 is a curve has no support. Full (A,B) invariance is
empirical only.

**Still needed.** A definition of C1 over Q (more primes, H6) and a decision whether C2 is a
curve at all (second plane, H5); an all-target exact sweep at k = 7 (H4); a proof of (A,B)
invariance for mu not in H; and the genericity conditions in a form a theorem can use
(Section 4, item N6).

### Item 3, degree-5 growth

**Established.** A split engine about 30x faster and 10x leaner than the shift engine at
k = 10, validated against two engines. d_5(10) = 18,569 (three instances) and
d_5(11) = 25,033 (two instances), each deficient-closed, so LFD >= 6 on those five
instances. |Z| <= 2^(k+1) is proved, so the k = 11 conclusion holds for any instance
attaining d_5(11) = 25,033. dim V_4 = 4k^2+172k−571 exactly at k = 6..13.

**Failed or corrected.** The bound d_5 <= 1 + 2k·lambda is vacuous at the measured depths.
The cost projections are untested model outputs.

**Still needed.** d_5(12), d_5(13) (H2, H3) and lambda at those depths; a genericity
argument (equal values on 2–3 instances are evidence only); d_6 for the step to 7.

### Item 4, exact value at depth 10

**Established (updated after the main session's run finished, 2026-11-08).** Exact LFD 6 at
depth 10 on the unit instance (p = 12289, y^2 = x^3 + x, target abscissa 9362), **exact**:
`experiments/ic-chain-depth10/degree6.json`. Lower side: V_5 verified closed with 1 not in V_5
(`probe.json`), reproduced bit for bit by a from-scratch rerun (`rerun.json`, basis sha256
`d4fe1334…`). Upper side: the block engine `csd/ic_chain_block_closure.py` from that V_5, eight
rounds of random combinations of variable shifts of low rows (every row lies in V_6, so the
certificate is one-sided and sound whatever the randomness); both swap parts reach full rank,
30318 + 30142 = 60460 = |M_<=6|, so V_6 = M_<=6 and 1 ∈ V_6. In the criterion of row 4.9,
f_6(10) = delta_5(10) = 3131. Round 7 shows the first falls: rank 59371 exceeds the
18569 + C(20,6) = 57329 a fall-free closure can reach, and round 8 shifts 18569 + 2042 = 20611
low rows. Also: a validated split engine with a deterministic closure certificate for unit and
proper ideals; a runner tested at k = 8, D = 6 (V_6 = I_6, exact LFD 5 there); and the memory
model that ends the dense method at k = 13.

**Failed or corrected.** The first depth-10 degree-6 runs were OOM-killed, including the
patched pid 11870; the run that finished embeds V_5 in checkpointed 1000-row batches. The
cols^2 time bracket does not reproduce. The negative control is weak.

**Still needed.** The two proper instances at D = 6 (H8); an independent second-engine
cross-check of the unit verdict (B7); for k >= 11, a new engine or a much bigger machine
(H14–H16).

### Item 5, fall-space structure

**Established.** The k-dependence of the degree-4 closure is explained. For k >= 7 the
closure is the free high-level extension of the depth-7 closure. That gives
dim V_4(k) = 4k^2+172k−571 and the fall-space dimension 4k+87, with an explicit closed
containing space W_k that excludes 1. A written reduction theorem turns one finite
certificate at depth 8 into the statement for every k >= 8 on that instance. The weight-2
fall element w has a closed form over Q(a,b,t) whose denominator a+3t^2 is the K3
nondegeneracy quantity. The fall space decomposes into explicit level-set pieces.

**Failed.** No closed form in (a,b,t) for the 829-dimensional base V_4(7); depth 6 cannot
serve as the base.

**Still needed.** Genericity of the base: which (p,a,b,t) share the pivot set and pass the
depth-8 certificate (H1b); the degree-5 analogue (base k_0 >= 11, H1a and H13).

### Item 6, binary and Koblitz

**Established.** 1871 exact certified values. Two-summand LFD stays far below the written
and kernel bounds (max 5 at k = 10, against 12 and 21). The direct m = 3 maximum is 7 at
k = 3, 4, 5. The split encoding is at most 5. Five written proofs (P1–P5) describe the
certificate degree D* and the cost of reading Z off it, and the data needed for a cost
statement is recorded per row.

**Failed or corrected.** The priority claim for binary LFD 5 is refuted (an earlier
full-field n = 8 system exists). Case and system counts were mixed. "Direct m = 3 exceeds
k+2" fails at k = 5.

**Still needed.** Grid completion (H11, H12), m = 2 at k = 11 (H17), a proved bound on D*,
a stopping rule that does not use |Z|, an operation count of the iterated closure, and
relation yield.

---

## 3. Heavy runs, in the order the main session should run them

Machine state measured at 2026-11-08 04:31 UTC: 15 GB RAM (about 14.3 GB usable in the
cgroup), 4 cores, **2.7 GB free disk** (less than the 5 GB assumed when the items were
planned), and pid 16280 (`scripts/run_ic_chain_depth10_degree6.py --combos 3000`) running at
7.0 GB RSS. Item 4 models that runner's late-round peak at about 10 GB. Two OOM incidents
have already happened on this host. Prefix every run with
`cd /home/user/ml-cryptanalysis; ulimit -v <cap>`, keep the BLAS thread variables at the
values shown, delete scratch checkpoints after each run, and do not run item scripts
whose default output path is a committed file without `--out` (the item 3 and item 4
verifiers flagged this).

The memory and time figures are the producers' estimates. Item 3's projections and all of
item 4's times are untested model outputs (rows 3.11, 4.12).

### Phase A: alongside pid 16280 (keep the total of everything else under about 1.5 GB)

| order | id | item | purpose | command | RAM | disk | time |
|---|---|---|---|---|---|---|---|
| A1 | H5 | 2 | Second plane sweep for C2: is it a curve? (and C1 points over F_257) | `for c in 0..9: python3 scripts/ic_bounds_locus_plane.py --prime 257 --k 6 --fix A --values <A≡c mod 10> --no-certify --out experiments/ic-bounds-locus/plane/p257-k6-plane-c<c>.jsonl` (each chunk about 25 min) | 0.2 GB | 0.02 GB | 4 h |
| A2 | H6 | 2 | Lift C1 to Q: about 10 more primes p ≡ 1 mod 64, about 15 lines each, then CRT | `python3 scripts/ic_bounds_locus_plane.py --prime $p --k 6 --fix A --values 1..15 --no-certify --no-zeros --out plane/p$p-k6-lines.jsonl` for p in 577, 641, 769, 1153, 1217, 1409, 1601, … (primes ≡ 1 mod 64 only; the producer's list also names 1667, which is ≡ 3 mod 64, and 1729 = 7·13·19, so check each p); then `python3 scripts/ic_bounds_locus_lift.py` | 0.3 GB | 0.01 GB | 5 h |
| A3 | H4 | 2 | Exact LFD at k = 7 for all 1909 targets (does C2 have LFD 4 or 5 at k = 7?) | `ulimit -v 4000000; for c in 0,1 1,0 0,3 3,0 88,32 194,158 200,59 87,256 148,227; do python3 scripts/ic_bounds_locus.py run --prime 257 --k 7 --curve $c --dmax 5 --extend-to 6 --out experiments/ic-bounds-locus/data-k7d5/p257-k7-$c.jsonl; done` (the held-out curves likewise) | 1.2 GB | 0.05 GB | 6 h |

A1 and A2 can run together. A3 runs alone beside pid 16280, or waits for Phase B.

### Phase B: after pid 16280 ends (this machine, one large job at a time)

| order | id | item | purpose | command | RAM | disk | time |
|---|---|---|---|---|---|---|---|
| B0 | — | 4 | **Done: reached, 1 ∈ V_6, exact LFD 6 (`experiments/ic-chain-depth10/degree6.json`); B2 is not needed.** Read pid 16280's result. "reached" (1 in V_6) gives exact LFD 6 at k = 10 for the unit instance; "closed-deficient" after a full verification pass gives LFD >= 7. A stall without the pass is neither. If it died or stalled unverified, run H7 instead. | — | — | — | — |
| B1 | H1a | 5 | Low profile of the depth-10 V_5 and the test V_5(10) ∩ M_{levels<=8} against the free extension (locates the degree-5 base). Only the owner of `~/.cache/ic-chain-depth10/` may run it, before that checkpoint is removed. | load the 770 MB basis (hash `d4fe1334…`) with `csd.ic_bounds_fall`; `F.low_profile`, `intersect_columns` with `level_mask_columns(grid, 2**9-1)` (exact snippet in `experiments/ic-bounds-fall/README.md`) | 8 GB | 1 GB | 0.5 h |
| B2 | H7 | 4 | Only if B0 gave no certified answer: exact k = 10 LFD, unit instance, with a deterministic verification pass. **Needs 4 GB disk: free disk first.** | `ulimit -v 12000000; export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4; python3 scripts/ic_bounds_exact10_degree6.py --k 10 --degree 6 --prime 12289 --a 1 --b 0 --target-x 9362 --v5 ~/.cache/ic-chain-depth10/rerun --v5-sha256 d4fe133411ff396d4b607165035615cde9e796f112f4a26bf23db1ef2612ee6b --tag k10-D6-p12289-x9362 --ckpt ~/.cache/ic-bounds-exact10/k10-unit` | 7 GB | 4 GB | 10 h |
| B3 | H1b | 5 | Genericity of the depth-7 base: pivot set and depth-8 certificate over many primes and random (a,b,t), recording which instances deviate. This is the input to the all-depth theorem (N5, N6). Reduce the stored output to pivot sets and certificate verdicts to stay within disk. | `for p in 12289 40961 61441 1073153; do for s in $(seq 1 200); do python3 scripts/ic_bounds_fall_analyze.py --prime $p --seed $s --kmax 8 --skip-moebius --unit-prime 0 --extra-instances '' --analysis-name analysis-k8-s$s.json; done; done` | 2 GB | 2 GB (reduce) | 3 h |
| B4 | H2 | 3 | d_5(12), two instances (unit seed 10, \|Z\| = 1352; random seed 12, \|Z\| = 1399). Sequential, or both at once with 2 threads each. | `ulimit -v 6000000; export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4; python3 scripts/ic_bounds_d5_probe.py --k 12 --save-every 1200` and `… --k 12 --a 4021 --b 9917 --seed 12 --want-solutions --save-every 1200` | 2.3 GB each | 0.8 GB each | 1 h each |
| B5 | H3 | 3 | d_5(13), p = 40961, two instances. **Strictly sequential, deleting the checkpoint between them** (2 × 1.8 GB exceeds the free disk). | `ulimit -v 9000000; export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4; python3 scripts/ic_bounds_d5_probe.py --k 13 --prime 40961 --save-every 1800`, then `… --k 13 --prime 40961 --a 4021 --b 9917 --seed 12 --want-solutions --save-every 1800`; then `python3 scripts/ic_bounds_d5_analyze.py` (writes into the committed directory) | 4.2 GB | 1.8 GB | 2.5 h each |
| B6 | H8 | 4 | Exact k = 10 LFD on the two proper instances (p = 13313, x = 7412, 98 points; p = 15361, a = 4021, b = 9917, x = 7776, 60 points). One at a time; **each needs 4 GB disk**. | `ulimit -v 12000000; python3 scripts/ic_bounds_exact10_degree6.py --k 10 --degree 6 --prime 13313 --a 1 --b 0 --target-x 7412 --tag k10-D6-p13313-x7412 --ckpt ~/.cache/ic-bounds-exact10/k10-p13313`; same with `--prime 15361 --a 4021 --b 9917 --target-x 7776 --tag k10-D6-p15361-x7776 --ckpt ~/.cache/ic-bounds-exact10/k10-p15361` | 7 GB | 4 GB | 15 h each |
| B7 | H13' | 3 | Independent cross-check of the k = 10 degree-6 verdict with item 3's engine (validation gate: a second engine). | `ulimit -v 6000000; export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4; python3 scripts/ic_bounds_d5_probe.py --k 10 --degree 6 --prev split --save-every 1800` | 2.5 GB | 0.95 GB | 5.3 h |
| B8 | H11 | 6 | Complete m = 2 n = 11, k = 10 (24 cases), together with H12 below. | `ulimit -v 3000000; export OMP_NUM_THREADS=1; python3 scripts/ic_bounds_binary_run.py --family m2 --n 11 --k 10 --workers 2 --submit-seconds 960 --out experiments/ic-bounds-binary/results-m2.jsonl` (repeat until "0 cases to run") | 2.6 GB | 0.01 GB | 2.5 h |
| B8 | H12 | 6 | Full m3-direct k = 5 grid (152 cases): does the plateau at 7 hold? | `ulimit -v 3000000; python3 scripts/ic_bounds_binary_run.py --family m3-direct --n 6 7 8 9 --k 5 --max-degree 14 --workers 2 --submit-seconds 960 --out experiments/ic-bounds-binary/results-m3-direct.jsonl` | 2.6 GB | 0.01 GB | 10 h |
| B9 | H9 | 1 | Direct m = 3 at k = 5: does 7 persist? (16384 monomials of degree <= 7) | `export OMP_NUM_THREADS=1; python3 -c "from csd import ic_dreg as I, ic_bounds_m3 as M; print(M.exact_profile('direct', I.PrimeCurve(97,0,1), 92, 5, engine='batched', max_columns=17000)['results'])"`, repeated for p in 97, 193, 257, 7681 and several targets | 6 GB | 0.01 GB | 3 h per instance |
| B10 | H10 | 1 | Split m = 3 (free auxiliary) at k = 5 (17548 mixed monomials of degree <= 6) | `… M.exact_profile('split', I.PrimeCurve(97,0,1), 92, 5, engine='batched', max_columns=18000) …` | 6.5 GB | 0.01 GB | 4 h per instance |
| B11 | H11b | 6 | Complete m3-split k = 4, n = 8, 9 (79 cases) | `ulimit -v 6000000; python3 scripts/ic_bounds_binary_run.py --family m3-split --n 8 9 --k 4 --max-degree 14 --max-bytes 2500000000 --workers 1 --out experiments/ic-bounds-binary/results-m3-split.jsonl` | 3 GB | 0.01 GB | 8 h |

B9 and B10 must not run at the same time (12.5 GB together). B8's two runs together use
about 5 GB and four cores. B4 must come before B5 so that the k = 12 result decides whether
the k = 13 runs are worth their disk.

### Phase C: needs a bigger machine (or more disk than this host has)

| id | item | purpose | command | RAM | disk | time |
|---|---|---|---|---|---|---|
| H13 | 5 | V_5(11) on the probe instance and the random-curve instance: is V_5(11) the free extension of V_5(10), or is a depth-11 base needed? (needs a checkpointing variant of the export script) | `ulimit -v 40000000; python3 scripts/ic_bounds_fall_export.py --prime 12289 --seed <instance> --k --k5 11 --product-cap 12000` | 32 GB | 10 GB | 24 h |
| H14 (done, section 7) | 3 | Fifth value d_5(14), p = 65537, uint32 storage (one check of a cubic) | `ulimit -v 20000000; export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4; python3 scripts/ic_bounds_d5_probe.py --k 14 --prime 65537 --save-every 1800` | 13.5 GB (needs the whole host) | 7.5 GB | 6 h (34 h worst case) |
| H15 (done, section 7) | 4 | k = 11, D = 6 on the unit instance: "reached" gives exact LFD 6, "closed-deficient" gives k_7 <= 11. Step 1 (verified V_5 at k = 11, `--degree 5`, 4 GB, est. 8 h) can be replaced by item 3's 15-minute split-engine closure if its basis is saved in the format the item 4 runner reads; that compatibility is not established. | `python3 scripts/ic_bounds_exact10_degree6.py --k 11 --degree 6 --prime 12289 --a 1 --b 0 --seed 10 --v5 <k=11 V_5 dir> --tag k11-D6-p12289-unit --ckpt $CKPT/k11-D6` | 18 GB | 13 GB | 60 h |
| H16 (done, section 7) | 4 | k = 12 and k = 13 (p = 40961), D = 6. k = 13 is the last depth the dense split engine reaches; k = 14 needs at least 0.5 TB or a new engine. | `python3 scripts/ic_bounds_exact10_degree6.py --k 12 --degree 6 --prime 12289 --a 1 --b 0 --seed 10 --tag k12-D6-p12289-unit --ckpt $CKPT/k12-D6`; `… --k 13 --prime 40961 … --tag k13-D6-p40961-unit --ckpt $CKPT/k13-D6` | 45 GB / 115 GB | 40 GB / 100 GB | about 300 h / 1400 h (extrapolated) |
| H17 | 6 | m = 2 at k = 11 (n = 13, 14, 15): does binary LFD reach 6? D = 6 has 110,056 columns, about 30 GB dense, so this host can certify only through D = 5. | `ulimit -v 8000000; python3 scripts/ic_bounds_binary_run.py --family m2 --n 13 14 15 --k 11 --max-columns 120000 --max-bytes 6000000000 --workers 1 --out experiments/ic-bounds-binary/results-m2-k11.jsonl` | 8 GB here (about 30 GB for D = 6) | 0.02 GB | 30 h |

---

## 4. The next proofs

Each is stated with its hypotheses, the containing space or certificate it rests on, and
the genericity conditions found so far. None is Lean-checked. N1, N2, N4 and N7 are short
and sit directly on existing kernel material; N5 and N6 are the route to an all-depth
generic statement.

**N1. Carry lemma (item 1, Prop. 1).** Let F = C ∪ {g}, where C are the chain equations of
any m blocks and depth k, and g is a chain normal form of degree d >= 2. Suppose a chain
variable u = u_{i,j} satisfies deg NF(u·g) <= d, and g, as a function on H^m, is not
supported in one level set {x_i^(2^j) = c}. Then LFD(F) >= d+1. Containing space:
W = (C ∩ M_<=d) + F_p·g is closed and contains V_d, because the chain Gröbner basis puts
the top forms of C in (u^2) and they cannot cancel g's squarefree top. NF(u g) is in I with
degree <= d and is not in W (C is radical). This needs only the chain Gröbner basis lemma
already used for K1.

**N2. Direct m = 3 lower bound (item 1, Cor. 2).** For every k >= 3, prime p with
2^k | p−1, curve y^2 = x^3+ax+b and target abscissa t with 5t^3+3at+6b != 0 mod p, the
system {chains, NF(S_4(x1,x2,x3,t))} has LFD >= 7. Inputs to formalize: the two
identities over Z[a,b,t], (i) [x1^3 x2^3 x3^3] S_4 = −8(5t^3+3at+6b) and
(ii) [x3^4] S_4 = S_3(x1,x2,t)^2; the consequence that (3,3,3) is the unique degree-6 mask;
and the slice count 6·2^k < 4^k for k >= 3. S_4 has 540 terms over Z[a,b,x1..x4]; the
identities are polynomial ring equalities (a `ring_nf` or certified-coefficient check).
Genericity condition: 5t^3+3at+6b != 0. The exceptional locus has LFD 7 at k = 3, 4
by exact computation, but no proof.

**N3. Two-summand tangent and order-two facts (item 2).** (a) If 3t^2+a = 0, t != 0 and
k >= 3, then dim(V_3 ∩ M_<=2) >= 3 > dim V_2 = 1. Explicit elements:
h1 = u11·g + 2t·u10·g and h2 symmetrically, independent of g by the monomials u12u21
and u22u11. This complements K3, whose hypothesis is 3t^2+a != 0. (b) If t^3+at+b = 0,
then S_3(x1,x2,t) = q^2 with q = x1x2 − t(x1+x2) − a − 2t^2, so q is in I (the chain algebra
is reduced). (c) For mu in H, (a,b,t) and (mu^2 a, mu^3 b, mu t) have equal closure
dimensions (weighted homogeneity; x -> x/mu is a weight-preserving automorphism commuting
with the shifts).

**N4. Two counting facts (item 3).** (a) |Z| <= 2^(k+1) for every two-summand subgroup
system (for fixed x1, S_3(x1, x2, X) is a nonzero polynomial of degree <= 2 in x2). Hence
dim(I ∩ M_<=D) >= |M_<=D| − 2^(k+1). Any verified closure with
dim V_D < |M_<=D| − 2^(k+1) gives LFD >= D+1 for every instance with that closure
dimension, whatever its solution set. (b) d_5 <= 1 + 2k·lambda. (c) (item 4) Given ev(M_<5)
onto Z: LFD >= 7 iff f_6(k) < delta_5(k), where f_6 = dim(V_6 ∩ M_<=5) − dim V_5 and
delta_5 = |M_<=5| − |Z| − dim V_5.

**N5. Reduction theorem for the degree-4 closure (item 5).** Work in
A_k = F_p[(Z/2^k)^2] with weight = popcount. Let S = V_4(7) for a fixed (p, a, b, c, t), and
W_k = ⊕_{h} tau_h(S ∩ M_<=4−wt h), with h ranging over exponent pairs supported on levels
>= 7 of weight <= 2. Then (i) W_k contains every admissible generator row, (ii) 1 is in W_k
iff 1 is in S, and (iii) if W_8 is closed under admissible variable products in A_8, then
W_k is closed in A_k for every k >= 8. The proof is carry combinatorics: a high variable
adds one bit to h; a low variable acts inside A_8 with at most one carry into level 7. It is
instance-free and is the natural Lean target. The base facts (V_4(7) of dimension 829,
1 not in V_4(7), W_8 closed: certificates of size 829 × 3473 and 1061 × 2517) are finite
exclusion certificates. They are larger than what the kernel has replayed without
`native_decide`. Consequence on any instance passing the certificate: 1 is not in V_4(k)
and dim V_4(k) <= 4k^2+172k−571 for every k >= 8. So a unit-ideal instance has LFD >= 5,
and a proper one has LFD >= 5 whenever |Z| < |M_<=4| − (4k^2+172k−571).

A point to settle before citing N5 across depths: for k >= 2, g = NF(S_3(c u1, c u2, t)) has
coordinates in A_k that do not depend on k (exponents <= 2 < 2^k), so the closure in A_k is
defined for every k. It is a summation system over a field F_q only when 2^k | q−1, and
closure dimensions are invariant under scalar extension. Written out, this would let N5
speak about one fixed (a, b, c, t) over the extensions F_{p^e}. That is still not the
cryptographic family, where p changes with k. This remark is a proof sketch written here
and has not been checked by a second agent.

**N6. Genericity of the base (items 2 and 5 together).** The statement to prove: there is a
proper Zariski-closed set E in (a,b,t)-space (or, if (A,B) invariance is proved, in the
normalized plane) such that for (a,b,t) outside E, the depth-7 closure has the pivot set
recorded in `experiments/ic-bounds-fall/p40961-s7/k7-D4.npz` and the depth-8 certificate
passes. The route is a parametric (pivot-polynomial) closure over Q(a,b,t), as in manuscript
Section 7, item 1: E is contained in the union of the vanishing loci of the pivots. Known
components of E:

* a+3t^2 = 0 (the tangent line; w degenerates; the denominator of Q_w);
* t^3+at+b = 0 (order-two targets);
* t = 0, with subcases a = 0, b = 0 (a in H or not), and ab != 0;
* C1: at p = 193, a cubic in (A,B), lift to Q unknown;
* C2: 190 points at p = 193, of unknown nature; and the special points (−3,2), (−3,170),
  (162,58), (32,103).

Three facts constrain the statement. Full (A,B) invariance is only empirical, so the plane
version needs a proof for mu not in H. C2 changes its LFD between k = 6 and k = 7 (row
V2.1), so E must be allowed to depend on the depth used for the base. If C2 is not the
reduction of a curve, E may not be the reduction of a fixed variety, and the genericity
statement must be per prime. Heavy runs H1b, H5 and H6 decide which form is possible.
Combined with N5, a proof of N6 would give "LFD >= 5 at every depth k >= 8 for generic
unit-ideal two-summand systems". That is one step above K3, and it would be the first
all-depth statement obtained from computed bases.

**N7. Certificate soundness, so that exact values can be replayed.** (a) The local-to-global
criterion in the general form of item 1 (neither radicality nor rational points required):
if dim V_D = |M_<=D| − dim R/I and M_<D maps onto R/I, then V_D = I_D, V_e = I_e for
e >= D, and the exact LFD is the last strict fall. (b) The split-closure certificate (items
3, 4): sigma-invariance of V_D, W = W^+ ⊕ W^- for p odd, and a clean sweep of the k block-1
shifts implies W = V_D. (c) Item 6's P1 (D* >= least D with sum_{i<D} C(r,i) >= |Z|) and P5
(D* <= d_full). These turn the computations into checkable certificates; they are not
bounds by themselves.

**N8. The upper-bound side (open, not yet formalizable).** No item produced an upper bound
below K1/K2 at general k. The candidates are LFD <= 7 for direct m = 3 at every k >= 3 off
5t^3+3at+6b = 0, and the m = 2 conjectural formula of manuscript Section 7 (LFD = the
least D with d_D(k) = |M_<=D| − |Z| and M_<D onto Z). Each needs a closed form for d_D(k)
plus a proof that the closure meets I_D at that degree. N5 supplies the k-dependence of d_4
only. An upper bound needs the analogous structure at the certificate degree (5, then 6),
which is conjecture 5.11 and heavy run H13.

---

## 5. What the data still cannot establish

* **Any bound at cryptographic sizes.** All exact values are at k <= 13 for d_4, k <= 11
  for d_5, k <= 9 for exact two-summand LFD (prime), k <= 4 for m = 3 prime, k <= 10 for
  binary m = 2 and k <= 5 for binary m = 3, over primes below 2^21 or fields up to F_{2^14}.
  A cryptographic system has k and p far beyond these. The only statements that reach every
  depth are K1–K4, Corollary 2 (N2), the tangent and order-two facts (N3), |Z| <= 2^(k+1)
  (N4a), and the per-instance reduction theorem (N5). Bounds at cryptographic sizes need
  proofs; no amount of data supplies them.
* **Upper bounds.** (Section 7 updates the depths: exact 6 up to k = 13 and LFD >= 7 at k = 14 on one
  unit instance; the statements below were written on 2026-11-08.) Every new lower bound here (LFD >= 6 at k = 10, 11; LFD >= 7 for direct
  m = 3) is exact or proved. No new upper bound below K1/K2 holds beyond the computed depths.
  The plateaus (direct m = 3 at 7; binary m = 2 at 4–5) are empirical.
* **Genericity.** Equal values on 2–3 instances at k = 10, 11 (item 3), and identical pivot
  sets on three instances (item 5), are evidence. Item 2 shows that the exceptional set is
  not simple (C1, C2, special points) and that it moves with depth.
* **Growth laws.** (Superseded in part by section 7: d_5 is now known at k = 10..14.) dim V_4 = 4k^2+172k−571 is explained for k >= 8 on certified instances.
  d_5(k) has two verified values; four values would fix only a cubic. The depth k_7 where the
  two-summand LFD reaches 7 is a conjecture (12..17). The binary threshold at |Z| ≈ 470 is
  a pattern in one grid.
* **Cost.** LFD is not the cost of index calculus. Item 6 records D*, |Z| and Macaulay
  sizes, but a cost statement still needs a proved bound on D*, the operation count of the
  iterated closure, a stopping rule that does not use |Z|, and relation yield. No item makes
  an ECDLP speed claim. None of these results is an improvement over the index-calculus
  baseline in the sense of the repository's success criterion, and none is a leaderboard
  row; their gap to the matched Pollard rho is not measured here.
* **The ordinary ring.** The modulo-chains premise (equal to the ordinary closure for
  D >= 2) underlies every prime-field number here. It was cross-checked in the full ordinary
  ring only at small k (k <= 5, and D <= 4 at k = 6). Its general proof should be cited, or
  formalized, alongside N7.

---

## 5. Status after the proof track (2026-11-08)

* **N1, N2: kernel-checked.** `CarryLower.chain_carry_lower`, `M3CarryLower.m3_lfd_ge_seven`:
  LFD >= 7 for the direct three-summand system at every k >= 3 when 5t^3+3at+6b != 0, over any
  field with 2 != 0 (the slice count is not needed: for k >= 3 the normal form is injective
  on the supports of S_4 and x1·S_4).
* **Chain bridge: kernel-checked.** `ChainClosureBridge.chain_mem_closure_iff`: for D >= 2 the
  ordinary closure is (chain ideal ∩ M_<=D) ⊕ (closure modulo the chains). The premise used by
  every computation is now a theorem.
* **N5: kernel-checked, conditional on the base.** `DepthReduction.W_closed` (instance-free,
  any k0, D, m) and `ChainModelMap.n5_lfd_ge_five`: a depth-7 base S containing the generator,
  excluding 1, with W_8(S) closed gives LFD >= 5 at every k >= 8 for unit ideals. Correction
  to the statement in §4: W_k sums over high patterns of weight <= 4, not <= 2.
  Written proof: `docs/proofs/N5-reduction.md`.
* **N6: analysed; the generic statement is out of reach.** `docs/proofs/N6-genericity.md`:
  closure dimension is not semicontinuous in either direction; the generic certificate is
  a parametric basis with denominator of degree >= 480 in the normalized plane (~10^10
  coefficients). Proved (written): Theorem 6.1, one certified point gives LFD >= 5 for all
  k >= 8 on its whole G_m-orbit (p−1 scalings); full (A,B)-invariance of the closure dimension
  is false, invariance of the certificate is true. Recommended heavy run: the certificate at
  all 37,249 plane points at p = 193 (~52 CPU-hours), giving the exact exceptional set E_193.
* Lean audit: 478 declarations, no project axioms, no placeholders.
* **Base certificate: kernel-checked for two instances (2026-11-08).** `N5CertInstanceP40961`
  (p = 40961, (a,b,t) = (21223, 9887, 4748), S = V_4(7), dim 829) and `N5CertInstanceP11`
  (p = 11, (2,5,7)): (B1), (B2), (G) discharged by 43 `decide +kernel` chunks each (~14 min,
  ~2·10^7 kernel ops), so 1 ∉ V_4(k) for every k ≥ 8 unconditionally. LFD ≥ 5 still assumes
  the unit ideal. Audit: 487 declarations.
* **B7 done:** the k = 10 exact value is confirmed by the independent verifier engine.
* **B4, B5 done (2026-11-09/10):** d_5(12) = 32057 and d_5(13) = 39657, each on two instances
  (unit-curve family and the random curve), identical profiles, closed by full deterministic
  sweeps (`experiments/ic-bounds-d5/verify/k1{2,3}-*-D5.json`): LFD >= 6 at k = 12, 13.
  Differences 6464, 7024, 7600; second differences 560, 576 (not a quadratic).

---

## 6. Status after the bounds track (2026-11-10)

* **Upper bound (proof, reviewed):** LFD ≤ k+2 for two-summand prime chains, p ≥ 5,
  nonsingular curve, subgroup and coset bases (`docs/proofs/U1-upper-bounds.md`, T1;
  red-team `docs/reviews-U1-upper-bound-redteam-20261110.md`: every step confirmed, 149,122
  ordinary-ring systems at 2 ≤ k ≤ 5 with a negative control). Attained at k = 2 only; not
  formalized; kernel bound stays 2k+1. No slope below 1 is proved.
* **Kernel (audit 513 declarations):** unconditional LFD ≥ 5 at every k ≥ 8 for two concrete
  instances (`N5CertWitnessInstance*`); depth restriction (`DepthRestrictionM2`); growth
  schema (`GrowthSchema`); exact LFD = 5 at k = 6 for p = 1153 (`ExactK6InstanceP1153`) and
  1 ∉ V_4(k) for every k ≥ 6 there (`ExactK6Deeper`).
* **Growth law (conjecture):** LFD(k) = ⌊(k+14)/4⌋ for generic two-summand systems; for unit
  ideals equivalent to 1 ∈ V_D(4D−11) and 1 ∉ V_D(4D−10) per degree; holds for D = 4, 5; first
  open case D = 6 (depths 13, 14) **(now measured on one unit instance: holds for D = 6, see
  section 7)**. 24 certified sweep systems on the exceptional loci
  (3t²+a = 0 at k = 2; y² = x³+x with target (0,0) at k = 6, 7) and one planted target
  (p = 449, k = 6) deviate.
* **Three summands (exact):** LFD 7 at k = 2..5 on twelve instances; V_7 = I ∩ M_≤7; no degree
  ≥ 7 falls behind at k ≤ 5. k = 6, D = 7 **has run on p7681-r: LFD ≥ 8 at k = 6 (section 7)**.
* **Binary/Koblitz:** two preregistered laws refuted; LFD 5 from k = 8; no bound in |Z| alone.

---

## 7. Status after the GPU runs and the frozen predictions (2026-11-10)

Sources: `experiments/ic-bounds-k14/` (`RESULTS.md`, `README.md`, `runs/*.json`, `POD.md`),
`experiments/ic-bounds-freeext/` (`predictions-2026-11-10.json`, `outcomes-2026-11-10.json`,
`records/D5-*.json`), `experiments/ic-bounds-m3-depth/` (`closures.jsonl`, `table.md`,
`summary.json`), and the manuscript `paper/curve-lfd-bounds.tex` (Section 7 of the paper, paragraphs
"Depths eleven to fourteen" and "Three summands"). New engine: `csd/ic_bounds_big.py`, a third
implementation, validated by equality of canonical reduced-echelon basis hashes with
`csd/ic_bounds_d5_verify.py` on all local cases and by a self-test on the card (CPU, GPU, reference).
Every run named below ended with a deterministic full sweep (every shift of every final low row
reduces to zero) and zero repairs. Label of every statement not marked otherwise: **exact**. Nothing here
is kernel-checked, and none of it is a cost or speed statement.

| item | before | now |
|---|---|---|
| Unit instance (p = 65537, y^2 = x^3 + x, R = (4271, 11720)) | not available | **proof (written, not formalized):** the solution set is empty at every depth k <= 15 (2-descent homomorphism x(P) mod squares for the 2-torsion point (0,0), a Frobenius case analysis, H inside the squares for k <= 15); enumeration confirms it at k = 10..14. Paper: Lemma "A unit instance at depths up to fifteen" |
| Upper side at k = 11, 12, 13 | LFD >= 6 only (exact value 6 only at k = 10) | exact value 6 at k = 10..13 on the unit instance: V_6 = M_<=6 and 1 in V_6 (dimensions 60,460 / 110,056 / 190,051 / 313,912), 1 not in V_5 (dimensions 18,569 / 25,033 / 32,057 / 39,657). Prediction `d6-full-through-13` held |
| Test at k = 14 (1 in V_6(13) against 1 not in V_6(14)) | beyond the engines | 1 in V_6(13); at k = 14 V_6 is closed, dimension 448,922 of 499,178 (codimension 50,256), 1 not in V_6(14), low profile (0, 0, 2, 154, 4451, 73587): **LFD >= 7 at k = 14**; with the kernel depth restriction and the lemma also at k = 15. Prediction `d6-falls-at-14` held. The exact value at k = 14 is open |
| d_5(14) | model projection | 47,849, profile (0, 0, 2, 154, 4323, 47849), closed, on the unit instance and on a random curve (p = 65537, seed 12, 4068 solutions). Prediction `d5-14` held |
| Profile recursion p_j(k+1) = p_j(k) + 2 p_{j-1}(k) + p_{j-2}(k) | D = 5 for k = 10..13 | D = 5 for k = 10..14 (checked on the five profiles). **Correction:** it is an identity of dimensions from k = 10; as subspaces the free extension of V_5(10) is not closed and differs from V_5(11); the base depth is 11 (W_12(V_5(11)) closed, equal to V_5(12) on the two instances where V_5(12) was recomputed; `records/D5-*.json`) |
| Growth law LFD(k) = floor((k+14)/4), unit ideals: 1 in V_D(4D-11), 1 not in V_D(4D-10) | holds D = 4, 5; D = 6 open | holds for D = 4, 5, 6 on the unit instance; still a **conjecture** for generic systems. First deficient depth 4D-10 with codimensions 10, 189, 3131, 50,256 for D = 3..6 (no formula for the codimension) |
| Random-curve replication (p = 65537, y^2 = x^3 + 4021 x + 9917) | not available | D = 5 at k = 13, 14: same dimensions and profiles; D = 6 at k = 13: V_6 = ideal part (312,884 = 313,912 - 1028), evaluation rank 1028, so exact value 6 there. D = 6 at k = 14: V_6 closed with dimension 448,922 (the same dimension, profile and part ranks as on the unit instance) against 495,110 for the ideal part, so LFD >= 7 there too |
| Three summands at k = 6 (instance p7681-r: p = 7681, a = 1876, b = 2678, t = 6223, 48 solutions) | k = 6, D = 7 queued | V_7 closed (full sweep), dimension 19,958 against 62,956 for I cap M_<=7, 1 not in V_7: **LFD >= 8 at k = 6** on this instance (one instance, one engine `scripts/ic_bounds_m3_depth.py`). V_7 = ideal part at k = 2..5 (exact value 7 on twelve instances) is unchanged |
| Heavy runs H14, H15, H16 | planned | done on a rented GPU on the p = 65537 unit instance (H14 as planned; H15 and H16 on that instance instead of p = 12289 and 40961) |

Frozen predictions (`predictions-2026-11-10.json`, committed before the runs): `d5-14`,
`d6-full-through-13`, `d6-falls-at-14` held; `d6-base-15` untested; `d5-closed-form` and `lfd-table`
untested beyond k = 14 (checked parts in `outcomes-2026-11-10.json`).

### Open items

1. **Exact value at k = 14.** Needs V_7(14): 1,683,218 columns, out of reach of the present engines.
   Whether the value is 7 as conjectured is untested; so is every depth beyond 14.
2. **Replication status.** One unit instance (j = 1728), one engine and one run at k = 14, D = 6; the
   k = 14 basis (about 45 GB) was not kept, only its hashes. The random-curve replication covers D = 5
   and D = 6 at k = 13 and 14 (same engine, one run each; bases not kept). A replay by a second engine or on other hardware at k >= 11 (D = 6) has not been done.
3. **Untested predictions.** `d6-base-15` (V_6 free-extension structure from k = 14), the closed form
   of `d5-closed-form` for k >= 15, and the `lfd-table` values from k = 14 on.
4. **Three summands.** Second instance and second engine at k = 6; exact value at k = 6 (degree-8
   closure, 106,762 columns); k = 7 (degree-7 closure, 198,440 columns). No prediction for depths beyond 6.
5. **Uniformity over curves, primes and targets.** The generic statement (N6) is unchanged: the
   exceptional loci move with depth and the evidence at D = 6 is one curve with j = 1728 plus one random curve.
6. **Upper bounds.** Unchanged: proved slope 1 in k (written), measured slope about 1/4; no mechanism relating
   V_D(k) to V_{D+1}(k+4) is proved.
7. **Optional: formalize the unit-instance lemma** (descent homomorphism and the Frobenius case analysis;
   or a kernel check of the enumeration at one depth).
8. **Housekeeping.** `experiments/ic-bounds-m3-depth/README.md` still describes the k = 6, D = 7 run as
   queued and not run; its `closures.jsonl`, `table.md` and `summary.json` already hold the result.

## 8. Status after the second batch and the relaunch (2026-11-11)

This section supersedes open items 1, 2 and 8 of section 7. Records and tables:
`experiments/ic-bounds-k14/RESULTS.md` (sections "Later runs on the second pod" and "Relaunch"), methods in
`LIFT.md` and `CERTCHECK.md`, outcomes in `experiments/ic-bounds-freeext/outcomes-2026-11-10.json`.

| Item | Status | Label |
|---|---|---|
| Exact value at k = 14 | **7** on three instances: unit ideals over p = 65537 and p = 163841 (1 in V_7, V_7 = M_<=7) and a random curve over p = 114689 (V_7 = the degree-7 part of the ideal, 1,680,848 = 1,683,218 - 2,370) | exact, by the degree lift (new method: written proof, 170 validation runs against direct closures, ideal-bound controls at k = 10 and k = 14; no independent certificate of V_7) |
| LFD >= 7 at k = 14 | four instances over three primes (65537 unit and random curve, 163841 unit, 114689 random curve), each dim V_6 = 448,922 | exact; second run with equal hashes on one instance; independent certificate check on three (65537 unit, 163841 unit, 114689 random curve); second run with equal hashes also on the 163841 unit instance |
| LFD = 6 at k = 13 | four instances over three primes | exact |
| Second engine | reference engine equals the GPU engine by basis hash at k = 14, D = 5 (rank 47,849) and confirms full rank at k = 10 and 11, D = 6 | exact |
| d_5 at k = 15, 16 | 56,649 and 66,073 with the predicted profiles (prediction `d5-closed-form` through 16) | exact |

Open after this:
1. Done on 2026-11-12: the lift on the unit instance over p = 163841 (second run of the closure with equal
   hashes, independent check, lift: 1 in V_7). See `RESULTS.md`, section "Job B".
2. Depth 15 and beyond: V_6(15) needs about 147 GB per swap part on the card (the 141 GB H200 is too small;
   180 GB and 288 GB cards were out of stock), then the lift would decide whether the value stays 7. The
   prediction `d6-base-15` is untested.
3. An independent certificate for a lift result (for example a second implementation of the lift, or a
   sampled derivation check), and a kernel-checkable form of any of the k >= 10 results.
4. Items 3 to 7 of section 7 are unchanged.

Operational lesson recorded in `experiments/ic-bounds-k14/RELAUNCH.md`: results on a rented machine must be
copied off it by a loop that does not depend on the agent session, and the lifetime cap must fit the funded
balance; one chain of results was lost on 2026-11-11 when the provider deleted a pod on an empty balance.

## 9. Status after the formalization pass (2026-11-12)

| Item | Status | Label |
|---|---|---|
| LFD ≤ k+2, two summands | Every depth k ≥ 2, every field of characteristic other than two, nonsingular curve, every target, subgroup base, equal depths: `SelectorExtension.lfd_le_depth_add_two_of_two_ne` (with roots of unity in the field: `SelectorBound.lfd_le_depth_add_two`). The hypotheses p ≥ 5 and 2^k \| p − 1 of the written statement are not needed. | kernel |
| Same bound for coset bases, unequal depths, raw generator (k+4) | `docs/proofs/U1-upper-bounds.md` | written |
| Exact LFD = 4 at k = 2 for unit ideals with 3t² + a ≠ 0 | `SelectorExtension.exact_lfd_four_at_depth_two` | kernel |
| Unit-instance lemma | `UnitInstance.unit_of_primitiveRoot`, `unit_zmod`; instances p = 65537 (k ≤ 15), p = 163841 (k ≤ 14) | kernel |
| Depth ≥ 14 corollaries | `DepthFourteen.p65537_summary`, `p163841_lower`: from the four recorded closures as hypotheses to 1 ∉ V_5(k) for k ≥ 10, 1 ∉ V_6(k) for k ≥ 14, LFD = 6 at k = 10..13, LFD = 7 at k = 14, no bound ≤ 6 at k = 15 | kernel, relative to recorded computations |

Audit: 547 declarations, no project axioms, no placeholders (`experiments/ic-lean/verification.json`).

Open after this: the singular case of the selector bound; a kernel bound with slope below 1 in k; any
kernel-checkable form of the recorded closures at k ≥ 10; the growth conjecture LFD(k) = ⌊(k+14)/4⌋.

## 10. Growth in the depth (2026-11-12)

| Item | Status | Label |
|---|---|---|
| 1 ∉ V_D(k) for every k ≥ k_D = D + Λ(D−1) ≤ 14D + 10 (k_D = 27, 39, 52, 64, 77 for D = 2..6), characteristic 0 or p > 2^(k+1), nonsingular curve, t ≠ 0 | `docs/proofs/N8-growth.md`, Theorem 5.1 | written proof; adversarial review `docs/reviews-N8-growth-redteam-20261112.md` (no fatal finding; the review predates the revised text); not formalized |
| ⌊(k−10)/14⌋ + 1 ≤ LFD(k) ≤ k + 2 for k ≥ 38 on unit-ideal instances with t ≠ 0; such instances exist at every depth (p a prime factor of 2^(2^k) + 1, y² = x³ + (n − X²)x, X and n non-squares) | Corollary 5.2; upper side and unit ideal are kernel (`SelectorExtension`, `UnitInstance.unit_zmod`) | lower side written |
| Analytic input: spread bound for minimal vanishing sums of monomials on the curve | Brownawell–Masser Theorem B in characteristic 0; characteristic p with p-th power coefficients attributed to Voloch (1985), original not retrieved; self-contained proof in the file | written |
| Splitting property FEI_k(w) on data | fails below 4w − 4, holds at k = 4..9 (w = 2), k = 8, 9 (w = 3), k = 12 (w = 4) | exact |

Open: the sharp statements 1 ∉ V_D(4D − 10) and 1 ∈ V_D(4D − 11) for general D; FEI_k(w) from 4w − 4 on for
every w (would give LFD ≥ ⌊(k+12)/4⌋ on unit ideals); the case t = 0; a kernel form (a conditional
formalization with the analytic input as a hypothesis is in progress, `formal/ICFallDegrees/Growth*.lean`).

## 11. Binary systems at k = 11 to 16 (2026-11-12)

Records and scoring: `experiments/ic-bounds-binary-deep/RESULTS.md`, `outcomes-2026-11-12.json` (57 of 57 systems
certified; predictions frozen before any closure at k >= 11, sha256 `7feee5db…e33e1`).

| Item | Result | Label |
|---|---|---|
| Bases with few solutions, n = 2k - 1, k = 11..16 | LFD 4 on 11 systems | exact |
| (14, 11), (17, 12), (16, 12): k + e = 19, 20 | LFD 4 on 24 affine systems (law B predicted 5) | exact |
| (18, 13): k + e = 21 | LFD 5 on 3 systems (post-hoc, A, A_Z predicted 4) | exact |
| Trace-zero hyperplane, n = k + 1 | LFD 5 at k = 11, 6 at k = 12 and 13 (every law predicted 5) | exact |
| Scaled subfield, order-two target | LFD 6 at k = 10 (two curves), 7 at k = 12: the inverse-family value k/2 + 1 | exact |
| A law for the affine class | none of the four frozen laws holds on all cells (best 45 of 52) | refuted as stated |

Replication (`REPLICATION.md`, frozen before the run): a second GPU execution with other batch and tile sizes reproduces all 57 records at every degree including the basis hashes (186 of 186 on proper subspaces); the CPU device reproduces the 46 smaller ones. Limits: one engine; no second implementation at k >= 11; k = 15, 16 have one system each.

## 12. Growth theorem in Lean, conditional (2026-11-12)

Audit: 582 declarations, no project axioms, no placeholders. Design note: `docs/proofs/N8-lean-plan.md`.

| Item | Lean | Label |
|---|---|---|
| Reduction: S_D(k) closed iff the splitting property; V_D ⊂ S_D; low-part criterion | `GrowthReduction.not_one_of_fei`, `GrowthCriterion.closed_iff_fei`, `not_one_of_crit` | kernel |
| Counting: Λ, window, no straddling; Λ(1..5) = 25, 36, 48, 59, 71; Λ(w) ≤ 13w + 23 | `GrowthCounting.*`, `GrowthMain.Lam_one` … `Lam_five` | kernel |
| Splitting from the spread bound | `GrowthSpread.fei_of_spread` | kernel, given `SpreadBound` |
| 1 ∉ V_D(k) for k ≥ D + Λ(D−1), and for k ≥ 14D + 10 | `GrowthMain.growth_theorem`, `growth_theorem_linear` | kernel, given `CurveModel` |
| Unit instances at every depth (p a prime factor of 2^(2^k) + 1) and ⌊(k−10)/14⌋ < d ≤ k + 2 for every last fall degree d, k ≥ 24 | `GrowthExistence.unit_instance` (unconditional), `growth_instance_window` | kernel, given `CurveModel` |
| The curve model: function field, Semaev relation for ξ₁ = x(P), ξ₂ = x(R − P), and the vanishing-sum bound (Theorem 3.3) | not formalized | written (published theorem plus the file's proof) |

The Lean statements do not mention "nonsingular" or "t ≠ 0": those are what the written proof uses to produce a
curve model. Not covered in Lean: Corollary 5.3 (free extension of the ideal spaces) and Proposition 7.2.

Label note (2026-11-12, after an owner comment on PR #101): Theorem `thm:growth` is labelled **written proof**. The
adversarial review covers an earlier blob (`6c98619e`) and does not carry over to the revised file (git blob `0e58d23c`; `a775f94e…` is its sha256, not a blob id). A second
review pinned to commit `09dd9760` (`docs/reviews-N8-growth-redteam2-20261112.md`) found no fatal error; the label
stays written proof. The Lean development is conditional on `GrowthMain.CurveModel`, which is
not kernel-proved.

## 13. Block-of-four upper bound: conjectural work (2026-11-12)

`docs/proofs/N9-block-four.md` (status: conjecture; no bound with slope below 1 is proved). Scripts
`scripts/ic_bounds_block4_*.py`, records `experiments/ic-bounds-block4/`.

| Item | Status | Label |
|---|---|---|
| Reduction: the block step BF(D) (1 ∈ V_D(4D−11) ⇒ 1 ∈ V_{D+1}(4D−7)) for all D ≥ 4 gives LFD ≤ ⌊(k+14)/4⌋ on unit ideals | written in N9 | proof |
| BF(3), BF(4), BF(5) | from the k14 records | exact |
| BF(6) (1 ∈ V_7(17)) | not computed | open |
| Unequal depths: 1 ∈ V_D(k₁,k₂) iff min(k₁,k₂) ≤ 4D−11 | D = 3 (84 pairs), D = 4 (89 pairs), D = 5 (5 pairs), two instances | exact per pair; empirical as a law |
| dim V_D(k₁,k₂) depends only on k₁ + k₂ above the threshold; frozen prediction for V_5(10,11) confirmed | D = 3, 4, one D = 5 case | exact per pair; empirical as a law |
| Doubling map w ↦ Q·w(x^{2^r}) as an all-depth lemma | refuted (no weight-≤1 multiplier at D = 3, 4) | exact |
| Top-block selectors | cost r degrees per r bits | refuted as a slope-1/4 route |
| Lift charging one degree per four bits of one coordinate | no proof | conjecture |

## 14. Three, four and five summands (2026-11-12)

`experiments/ic-bounds-msummands/` (README, predictions frozen before the closures, sha256 `deb4b304…aad7`). Generators of
S_(m+1) built by evaluation on the grid, validated against the existing m = 2 and m = 3 generators.

| Cell | Result | Label |
|---|---|---|
| Prime m = 3, k = 2, 3, 4 | LFD 7 (18 instances) | exact |
| Prime m = 4, k = 2 / 3 / 4 | LFD 9 / 13 / 13 (6 instances each); at k = 4 deg g = 12 < mk = 16 | exact |
| Prime m = 5, k = 2 / 3 | LFD 11 / 16 (6 instances each) | exact |
| Prime m = 3, k = 6 | LFD ≥ 8 (recorded); D = 8 closure stopped at the 2-hour limit | open |
| Prediction "LFD = deg g + 1, one strict fall" | 48 of 48 prime instances | exact per instance; empirical as a law |
| Prediction "the bound mk + 1 is attained" | fails on all 6 instances at m = 4, k = 4 | refuted |
| Binary m = 3, 4, 5, r ≤ 21 | 57 exact values (48 local, 9 pod), all ≤ min(r, d_in + 2); point rule min(r, d_in + 1): 29 of 32 local (scored where d_in = r), 8 of 9 on the pod (d_in a bound at r = 21); 2 not certified | exact; empirical as a law |
