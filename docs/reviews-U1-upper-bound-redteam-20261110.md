# Red-team review of U1: T1 "LFD ≤ k + 2" for two-summand prime chain systems

Date: 2026-11-10. Role: Red Team. This file challenges the argument and reports tests. It
changes no status and edits no evidence; the Coordinator and the Validator decide what is
recorded.

Verdict words: **confirmed**, **gap** (missing argument stated), **refuted** (counterexample
given). Defect words: **wrong**, **unsupported**, **unpriced**. Evidence labels as in the
reviewed document: kernel, proof, exact, empirical, conjecture.

## 0. What was reviewed

- `docs/proofs/U1-upper-bounds.md`, sha256
  `676c1c50dbec22726f4c063b7f123a1d462a0c139e9cd2ebf484da89edb55113`, committed in `185ae42a`
  and unchanged at `085c208f` (branch `claude/lfd-proofs-20261108`; both commits are titled
  "WIP: bounds-track agent snapshot (in progress; not imported, not audited)").
- `experiments/ic-bounds-upper/`: scripts and `results-*` files of the same commit. The `*.log`
  files its README names are excluded by `.gitignore` (`experiments/**/*.log`). They are not
  frozen evidence and I did not use them.
- Definitions: `formal/ICFallDegrees/Closure.lean` (`InClosure`, `HasLFDBound`),
  `LocalReductionBound.lean` (`lfd_bound_of_local`), `ChainClosureBridge.lean`
  (`chain_mem_closure_iff`, `quot_mul_mem`), the statement of
  `ChainTranslation.normalForm_X_mul_monomial`, `PrimeM2General.lean` (the generator).
- Companion proof: Theorem `thm:uniform` of `paper/binary-lfd.tex` (bound min(k,ℓ)+2).
- Not done: no Lean build. The kernel status of the bridge is taken from
  `formal/ICFallDegrees.lean` and `Audit.lean`, which import it and print its axioms.
  P3 (kernel) was not reviewed, as instructed.
- T1 is a degree bound for an algebraic system. It makes no cost, ratio or rho claim, so no arm
  is priced and nothing here is evidence about ECDLP cost.

## 1. Verdict on T1

I could not break T1. On the scale the dispatch asked for, **T1 is confirmed** in this scope:
p ≥ 5 prime, k ≥ 1 with 2^k ∣ p−1, a, b, t ∈ F_p with 4a³+27b² ≠ 0, chains of depth k in both
blocks, generator NF(S_3(x₁,x₂,t)), subgroup base or coset base x_i = s_i·u_{i,0}:
**LFD ≤ k + 2**.

I read every step of Sections 1 to 2.5 and found no error and no gap that affects the
statement. The degree accounting closes at budget k₁+2 with nothing hidden: the selector is
multiplied once, into the generator (k₁ + 2), and every later product is one variable times a
row of weight at most the budget minus one (k₁+1 in the fibre derivations, max(r−1, k₁+1) in the
induction). The passage to the ordinary polynomial ring is the kernel theorem
`chain_mem_closure_iff`, whose reduction step costs no extra degree
(`normalForm_reduction : InClosure G p.totalDegree (p − NF p)`).

I then tested the statement in the ordinary ring with the independent validator's engine, on
instances chosen where the bound is tight or the proof branches. Result: **236,267 systems
(227,160 normalized, 9,107 raw), no exact LFD above the bound, none left uncertified.** The
bound is attained on 202,040 of the normalized ones, and a negative control with generators
outside the hypothesis does exceed it, so the test could have failed. No certified value in the
repository exceeds k+2 either (1,294 two-summand prime systems, k = 2..9).

Four qualifications, detailed in Section 2:

1. T1 is a written proof on top of a kernel-checked bridge. The theorem itself is not
   formalized.
2. T1 needs a nonsingular curve. The kernel bound 2k+1 does not.
3. For unequal depths the proof is the same, but the bridge is kernel-checked only for equal
   depths.
4. For equal depths the bound is attained only at k = 1 and k = 2. For k ≥ 3 no computed
   instance attains it.

## 2. Findings, most severe first

No finding refutes T1 or any other labelled claim of the document. The findings concern scope,
labels and wording.

**F1. T1 cannot "replace" the kernel bound 2k+1; it sits beside it (scope; unsupported if
worded as a replacement).**
`prime_m2_lfd_bound` (paper Corollary `cor:m2`) is kernel-checked for every field of
characteristic ≠ 2, every a, b, t and k ≥ 2, with no condition 2^k ∣ p−1. T1 is a written
proof, needs 2^k ∣ p−1 as stated, and needs 4a³+27b² ≠ 0: for a singular cubic the fibre
polynomial at α = t vanishes identically when t is the double root and lies in the base, and
the proof stops there. (The condition 2^k ∣ p−1 looks removable by scalar extension, since LFD
is invariant under it; the document does not claim this and I did not check it.) The paper currently
says the bound holds "at every depth, prime, curve and target"; with k+2 that sentence needs the
hypothesis. Measured on singular parameters (outside T1): 2,570 systems (every singular (a,b)
and every t at (p,k) = (5,1), (13,1), (5,2), (13,2), (17,2), (29,2), (17,3); thinned at
(41,3), (17,4)), 52 of them with a vanishing fibre: exact LFD ≤ k+2 in all, and LFD ≤ 3 on the 52.
So k+2 is not falsified there, and it is not proved there.
Cheapest settlement: keep 2k+1 as the kernel entry and add k+2 as a written entry with the
nonsingularity hypothesis (Section 7).

**F2. The two "exact" supports the document gives for T1 are not independent tests of it
(evidence; unsupported as independent confirmation).**
- Section 2.6 (`selector_bound_check.py`) replays the proof in the group algebra, the same
  model in which the proof is written, and relies on the same bridge. It is a correct and
  complete instance-level check there (it builds m − π(m) for every monomial, all degrees). It
  has no negative control and no coset instance.
- Section 2.7 (`lfd_structure.py`) takes the generation condition for D' ≥ k+2 from Theorem 1
  (its docstring says so), scans only to k+2 under an `assert`, and writes
  `bound_k_plus_2_violations=0` as a literal constant. It establishes V_D' = I_D' for
  LFD ≤ D' ≤ k+2. For k ≥ 3 its "exact LFD" is conditional on the theorem under review for
  the degrees k+3..2k (for k ≤ 2 the kernel bound covers that tail). The sentence "No instance
  exceeds k + 2" therefore could not have come out otherwise above k+2.
Neither script is wrong. Before this review the only test of T1 outside the proof's own model
was the comparison with the certified sweep, where the bound is slack by at least 1 for every
k ≥ 3 and is reached only at k = 2 (217 of 220 systems). Section 5 below supplies the missing
test: ordinary ring, independent engine, a stabilisation certificate that does not use
Theorem 1, tight cases, a negative control.

**F3. "The constant 2 is optimal" is true only in a narrow sense (wording).**
T1′ is confirmed as stated: LFD = k+2 occurs at k = 1 and k = 2, so no bound k + c with c < 2
holds for all k. It does not say the bound is sharp. At k = 1 it equals the dimension bound.
For k ≥ 3 the largest exact value is 4, 4, 4, 5, 5, 5, 5, 6 at k = 3..10 against 5..12, and
nothing excludes, for instance, LFD ≤ k+1 for all k ≥ 3. In the paper's range k ≥ 2 the bound
is attained at exactly one depth. The paper must not call the bound sharp or tight without that
restriction. Two exact facts from Section 5 can be cited instead:
at k = 2, LFD = 4 for every instance with 3t²+a ≠ 0 and LFD = 3 for every instance with
3t²+a = 0 (132,303 systems, unit and proper ideals, no exception); and for unequal depths the
bound min(k₁,k₂)+2 is attained for min(k₁,k₂) ∈ {1, 2} at every longer depth tested.

**F4. Unequal depths are outside the kernel-checked bridge, and k = 1 is outside the Lean
presentation (gap in formal support, not in the mathematics).**
`chain_mem_closure_iff` is stated for `chainInputs K m k`, m blocks of one depth k, and
`PrimeM2General.inputs` requires 2 ≤ k. The document says so for unequal depths. The generic
theorem `mem_closure_iff` needs only `ReducedIdealTrivial G`, so the missing kernel pieces are
that property and the translation lemma for chains of two depths. Exact ordinary-ring support
from this review: 58,720 systems with k₁ ≠ k₂, no violation (Section 5). The repository
defines no unequal-depth presentation, so the paper should state this clause as a remark.

**F5. Propositions 4 and 5 are stated one degree weaker than they are used (statement/use
mismatch; the stronger statement is true).**
The proof of Proposition 4 shows V_{D'+1} ⊂ I_{D'} + Σ_v v·I_{D'} for every D' ≥ δ. Hence
LFD ≤ D gives generation for every D' ≥ max(D−1, δ), not only D' ≥ D. The text uses the stronger
form twice: "a failure of generation at D' forces LFD ≥ D'+2", and "LFD ≥ 4 follows from
Proposition 5 alone (435 > 17·2)" for k = 8. From Proposition 5 as written (D' ≥ D) the k = 8
numbers give only LFD ≥ 3. With D' ≥ max(D−1, 2) they give LFD ≥ 4. The claim is right; the
displayed statements should be the strong ones.

**F6. Two parenthetical values have no certificate (unsupported).**
- Section 6.2: "(The true value is 5.)" for (p; a, b; t) = (257; 1, 3; 190), k = 8. This
  instance is not in the certified sweeps (their p = 257, k = 8 subgroup systems are (1,0;190),
  (1,0;242), (108,166;38), (108,166;78)). The producer's own data give LFD ≥ 5 exactly
  (`results-depth-restriction.json`: dim V_4(8) = 1061 < 2255 = dim I_4). No upper certificate
  exists for this instance. Cheapest test: V_5(8) for it (window 6,885) and generation at
  D' = 5..9, or one run of the certified pipeline.
- Section 2.8: "The exact values (4, 5, 6 for k ≤ 13)". Exact values exist for k ≤ 10. At
  k = 11, 12, 13 only LFD ≥ 6 is certified (`docs/IC_BOUNDS_PROGRAM.md`, items H15, H16).

**F7. E2 is worded the wrong way round (wrong as worded; the arithmetic is right).**
"4D−10 is the first depth at which f_D fits below |M_≤D|" is false: f_5(7) = 2377 < 3473 and
f_4(3) = −19 < 57. What holds is that 4D−11 is the **last** depth with f_D(k) > |M_≤D(2k)|
(f_D exceeds only at k = 4, 5 for D = 4 and k = 8, 9 for D = 5). `slope_counts.py` computes
exactly this ("last_depth_free_impossible"), and Section 6.1 concedes it in a subordinate
clause; the table row E2 and the phrase "from the first depth at which the count allows it"
should be changed.

**F8. The raw-generator bound k+4 is correct and loose (not an error).**
For D ≥ 4 the raw and the normalized inputs have the same closure
(`PrimeM2General.raw_normalization_identity`: the difference is a chain combination of degree 4),
so the raw presentation satisfies LFD ≤ max(4, min(k₁,k₂)+2), the form used in the binary
paper. Measured: exact LFD = 4 on all 9,107 raw systems tested.

**F9. C1 should name its exceptional loci (wording).**
24 certified sweep systems deviate from ⌊(k+14)/4⌋: 3 at k = 2 with LFD 3, all with 3t²+a = 0,
and 21 at k = 6, 7 with LFD 4, all on y² = x³+x with target (0,0), 9 of them proper ideals.
Both loci are known in the paper. "Generic" in C1 must exclude them, and the reduction to U_D,
L_D is for unit ideals only, as the document says.

**F10. E1 "proof above the base" is conditional (label).**
It assumes the free-extension decomposition, which is exact only at the certified depths. The
label should read "proof, given the free extension". The empirical part is consistent with a
spot check here (24 pairs in the ordinary ring, image inside and onto).

**F11. Formalization and evidence hygiene (process).**
Section 7 is an outline, not started. Under the owner rule "Mathematical results should include
Lean proofs", T1 stays "written proof, formalization open" with a proof-map entry. The reviewed
files sit in commits marked "not imported, not audited". The README of
`experiments/ic-bounds-upper/` names two log files as outputs (`selector.log`,
`depth-restriction.log`) that git ignores.

## 3. Proof audit, step by step

| step | content | verdict |
|---|---|---|
| §1, (1.1) | wt(a·b) ≤ wt(a)+wt(b) in F_p[Z/N₁ × Z/N₂] | confirmed (a carry or a reduction modulo 2^k never raises the Hamming weight) |
| §1, bridge and (1.2) | ordinary closure = chain part ⊕ V_D; LFD ≤ D ⇔ I_D' ⊂ V_D' for D' ≥ D | confirmed for equal depths: `chain_mem_closure_iff` (D ≥ 2, any further inputs) and `normalForm_X_mul_monomial` (the variable acts by e_i ↦ e_i + 2^j mod 2^k). `InQuotClosure` is the document's V_D. Unequal depths: same argument, not formalized (F4) |
| Lemma 1.1 | budgeted products | confirmed; kernel counterpart `quot_mul_mem` |
| Theorem 1, statement | g of x₂-degree ≤ 2, no fibre polynomial identically zero ⇒ LFD ≤ k₁ + max(2, δ) | confirmed. Remark: for k₂ = 1 the decomposition g₀+g₁x₂+g₂x₂² is not unique in A (x₂² = 1), so the hypothesis is about a chosen lift; for S_3 the lift is S_3. No effect on the result |
| Lemma 2.1 | e_α: weight k₁, eigen-property, partition of unity, A = ⊕ B·e_α | confirmed (needs p ∤ N₁ and μ_{N₁} ⊂ F_p; both hold) |
| Lemma 2.2 | S_3(α, X, t) ≠ 0 for every α, t ∈ F_p when 4a³+27b² ≠ 0 | confirmed. Identities checked symbolically over Z: X²-coefficient (α−t)²; at α = t the X-coefficient is −4f(t) and the constant is f′(t)² − 8t·f(t). Brute force: 6,669,000 fibres, none zero (Section 5) |
| Lemma 2.3, deg q = 0 | 1 ∈ U(q) | confirmed |
| Lemma 2.3, deg q = 1 | w_j = u_j − β^{2^j}; (u_j + β_j)·w_j = w_{j+1}; end 1 − β^{N₂} | confirmed; every product has an operand of weight ≤ 1 |
| Lemma 2.3, deg q = 2 | recursion w_j, w′_j; X·w_j + b_j·Q₁ = w′_j; u_j·w_j + a_j·w_j + b_j·w′_j = w_{j+1} | confirmed by expansion, including j = 0 (w₀ = 0, w′₀ = Q₁), j = k₂, and k₂ = 1. w′_j has weight 2 and is never multiplied |
| Lemma 2.3, cases r = 0, 1, 2 | span{w_{k₂}, w′_{k₂}} = (d)/(q) inside span{1, X} | confirmed. Double root in H₂ gives r = 1 and is covered. Roots outside H₂ and irreducible q give r = 0. q = c(X²−1) at k₂ = 1 gives q̄ = 0, r = 2, all corrections zero |
| Lemma 2.4 | U(q_α)·e_α ⊂ V_D, D = k₁ + max(2, δ) | confirmed: wt(e_α)+wt(g) ≤ k₁+2; rows w·e_α with wt(w) ≤ 1 have weight ≤ k₁+1 ≤ D−1 |
| §2.5, decomposition | I = ⊕ q̄_α B e_α; A = I ⊕ S; wt(S) ≤ k₁+1 | confirmed (R_α ⊂ {1, X}) |
| §2.5, (C1), (C2) | 1 − π(1) and v·ρe_α − π(v·ρe_α) in V_D | confirmed |
| §2.5, induction | m − π(m) ∈ V_{max(D, wt m)} | confirmed: the operand has weight ≤ max(r−1, k₁+1) ≤ max(D, r) − 1 |
| §2.5, conclusion | I_D' ⊂ V_D' for D' ≥ D, then (1.2) | confirmed |
| Corollary 1, cosets | fibre S_3(s₁α, s₂X, t) ≠ 0 | confirmed. With s₁ = s₂ the system is the subgroup system of a rescaled curve (S_3 is weighted homogeneous); s₁ ≠ s₂ is new and was tested here |
| Corollary 1, unequal depths | min(k₁,k₂)+2 by symmetry of S_3 | proof confirmed in A; bridge not formalized (F4) |
| Corollary 1, raw input | min(k₁,k₂)+4 | confirmed and loose (F8) |

## 4. The seven attack points of the dispatch

**(a) Degree accounting under the ordinary rule.** Confirmed. The proof lives in V_D, whose
only rule is one variable times an element of weight ≤ D−1. Products used: e_α·g
(k₁ + 2, by Lemma 1.1); u_{2,j} times rows w·e_α with wt(w) ≤ 1 (operand weight ≤ k₁+1); one
variable times m′ − π(m′) (operand weight ≤ max(r−1, k₁+1)). In the ordinary ring each becomes a
multiplication of total cost ≤ max(D, r) followed by a chain reduction at the same degree. The
tight cases confirm there is no hidden +1: at k = 2 the exact LFD is never 5 on 132,303 systems
(it is 4 on 127,599 of them); for depths (1, k₂) it is 3 and never 4 on 42,592 systems.

**(b) The idempotent e_α.** Confirmed. As a polynomial it is Π_j (1 + α^{−2^j}u_{1,j})/2, squarefree
of total degree k₁. It is used as a multiplier exactly once, on the input g of degree 2. No
derived row is ever multiplied by e_α. Afterwards e_α only rides along inside rows.

**(c) The univariate lemma.** Confirmed, all degenerate cases (see the table). Tested as a
statement in the ordinary ring: the chain system of q̄ has exact LFD ≤ 2 for all 31,646
polynomials q ≠ 0 of degree ≤ 2 tried (every q at p = 5, 13, 17 and k₂ ≤ 4; all split, linear,
constant and a sample of irreducible q at (41, 3), (97, 5)).

**(d) The fibre condition.** Confirmed for every t, including t = 0, f(t) = 0, 3t²+a = 0, and
for j = 0, 1728. The argument is two polynomial identities and "f(t) = f′(t) = 0 means a double
root". Useful by-product: disc_X S_3(α, X, t) = 16·f(α)·f(t), so order-two targets make every
fibre polynomial a square; these were included in the tests.

**(e) `lfd_bound_of_local`.** Its hypotheses are met with A′ = functions on Z, s = the reduced
representative in S, σ = k₁+1, D = k₁+2: `section_id` is A = I ⊕ S; `bounded` is wt(S) ≤ k₁+1;
`room` is σ+1 ≤ D; `vanish` is clear; `local_reductions` quantifies over all polynomials p of
degree ≤ D, reduced or not, and follows from the induction of §2.5 for wt ≤ D together with the
bridge, because NF(p − s(q p)) = NF(p) − π(NF p). The written proof does not in fact call this
theorem: §2.5 runs the induction in A for all degrees and ends with (1.2). Both routes are valid.

**(f) Cosets and unequal depths.** See Corollary 1 rows above, F4, and Section 5.

**(g) Characteristic.** p is odd because 2^k ∣ p−1. The document fixes p ≥ 5 in Section 1, while
row T1 of its table says "every prime p". For p = 3 (k = 1 only) the argument works verbatim for
curves y² = x³+ax+b with a ≠ 0, but those are not all curves in characteristic 3, and at k = 1
the bound equals the dimension bound. The paper should say p ≥ 5.

## 5. Independent tests

### 5.1 Method

Own builder (`experiments/ic-bounds-upper-redteam/rt_lib.py`): chains and generator written
from the S_3 formula, as polynomials in the k₁+k₂ chain variables; coset scales; unequal depths;
raw option; zero set by enumeration. Engine: `experiments/ic-curve-validation/code/lfd_engine.py`
(sha256 `75b34ce1e31d7b677409a3bc153691de66c97cf4125dfc646d8b7b7b127510ea`), function
`exact_lfd`: the ordinary-ring closure V_{F,D} on all monomials of degree ≤ D, I_D from
evaluation ranks, and the certificate "V_D = I_D and the monomials of
degree < D separate Z", which implies V_D' = I_D' for all D' ≥ D without Theorem 1. Nothing is
imported from `experiments/ic-bounds-upper/`.

Target, fixed in `rt_run.py` before the runs: a **violation** is an exact LFD above
min(k₁,k₂)+2; an instance is **uncertified** if the scan (to bound+1, or to the largest degree
with at most 5,200 monomials) ends without a certificate. Either would have been reported.

Checks of the test itself:

- S_3 as coded is the summation polynomial: 326,280 triples of abscissae on all nonsingular
  curves over F_5, F_7, F_11, F_13, F_17 against the group law.
- Engine closure equals the literal Lean rule (`closure_naive`) on 66 (instance, D) pairs.
- Builder and engine reproduce the certified sweep: exact LFD and |Z| agree on 690 of 690
  distinct two-summand prime systems with k ≤ 4 (231 subgroup, 459 random coset).
- The histograms of the document's table in 2.7 are reproduced exactly in the ordinary ring:
  k = 1 (p = 5, 7, 13): 252 with LFD ≤ 2, 2,170 with LFD 3; k = 2 (p = 5, 13, 17): 416 and
  6,336; k = 3, p = 17: 4,624 with LFD 4.
- Negative control (`rt_control.py`): same chains and pipeline, generator a random element of
  chain degree 2 (so fibre polynomials of degree 3). Exact LFD exceeds min(k₁,k₂)+2 on 280 of
  300 systems at depths (1,3), 300 of 300 at (1,4), 58 of 60 at (2,4). With S_3 at the same
  depths it never does (18,496 systems). At equal depth 2 random generators also stay ≤ 4, so
  equal depths alone would not discriminate; the unequal depths do.

### 5.2 Results

All values are exact LFDs in the ordinary ring. "Ratio" is the coset ratio s₂/s₁ (x₁ ∈ H,
x₂ ∈ cH); with s₁ = s₂ a coset system is a subgroup system of a rescaled curve, so only the
ratio matters. "Special" means: targets t = 0, f(t) = 0, 3t²+a = 0, t in the base, t in the base
with f(t) = 0 (a constant fibre polynomial), on j = 0, j = 1728 and random curves, plus random
(a, b, t).

| depths (k₁,k₂) | instances | systems | bound | exact LFD: count | at the bound |
|---|---|---|---|---|---|
| (1,1) | every nonsingular (a,b), every t; p = 3, 5, 7, 11, 13; every ratio | 19,318 | 3 | 0: 180; 2: 1,470; 3: 17,668 | 17,668 |
| (2,2) | every nonsingular (a,b), every t; p = 5; 13 (3 ratios); 17 (4 ratios); 29 (2 ratios); 37 | 121,060 | 4 | 3: 4,336; 4: 116,724 | 116,724 |
| (2,2) | special; p = 41, 53, 61, 73, 97, 101, 109, 113, 12289, 40961; 3 ratios | 11,243 | 4 | 3: 368; 4: 10,875 | 10,875 |
| (3,3) | every nonsingular (a,b), every t; p = 17; both ratios | 9,248 | 5 | 4: 9,248 | 0 |
| (3,3) | special; p = 41, 73, 89, 97, 113, 12289; 2 ratios | 5,847 | 5 | 4: 5,847 | 0 |
| (4,4) | special; p = 17, 97, 113, 193, 241, 257, 12289 | 1,580 | 6 | 4: 1,580 | 0 |
| (5,5) | special; p = 97, 193, 257 | 144 | 7 | 4: 144 | 0 |
| (1,2), (2,1) | every nonsingular (a,b), every t; p = 5, 13, 17 (4 ratios at 17) | 25,248 | 3 | 3: 25,248 | all |
| (1,3), (1,4) | every nonsingular (a,b), every t; p = 17 | 13,872 | 3 | 3: 13,872 | all |
| (1,5) … (1,12) | special; p = 97, 193, 257, 769, 12289, 40961 | 3,472 | 3 | 3: 3,472 | all |
| (2,3), (2,4) | every nonsingular (a,b), every t; p = 17 | 13,872 | 4 | 3: 768; 4: 13,104 | 13,104 |
| (2,5) … (2,10) | special; same six primes | 1,152 | 4 | 3: 75; 4: 1,077 | 1,077 |
| (3,4) … (3,8) | special; p = 17 … 40961 | 976 | 5 | 4: 976 | 0 |
| (4,5), (4,6), (4,7) | special; p = 97, 193, 257 | 112 | 6 | 4: 112 | 0 |
| (5,6) | special; p = 193, 257 | 16 | 7 | 4: 16 | 0 |
| **all normalized** | | **227,160** | | **violations: 0; uncertified: 0** | **202,040** |
| raw generator, (1,1), (2,2), (1,2), (3,3) | planes at p = 5, 13, 17; special at (3,3) | 9,107 | k+4 | 4: 9,107 | 0 |

Of the normalized systems, 58,720 have k₁ ≠ k₂ and 98,309 have a coset ratio other than 1.

Readings:

- **No violation and no uncertified system** in 236,267 runs (7 jobs, 8,942 s of one core).
- **Smaller depth 1** (61,910 systems): LFD = 3 whenever the depths differ (42,592 of 42,592).
- **Smaller depth 2** (147,327 systems, equal and unequal depths, unit and proper ideals):
  LFD = 3 exactly when 3t²+a = 0 and LFD = 4 otherwise, with no exception. So the bound is the
  exact value off the tangent locus, for every longer depth tested (up to 10).
- **Smaller depth ≥ 3** (17,923 systems): LFD = 4 always; the bound is never attained.

These three patterns are empirical; each value is exact for its system.

### 5.3 Certified data of the repository

`sweep_compare.py`: all two-summand prime rows of `experiments/ic-curve-sweep` (1,262 systems,
422 subgroup and 840 random coset, k = 2..7), `ic-chain-depth8` (22) and `ic-chain-depth9` (10),
joined to their results. Every row is exact (lower = upper). Maximum 4 at k = 2..5 and 5 at
k = 6..9. No exact value and no lower bound exceeds k+2. The bound is attained only at k = 2.

## 6. The other claims

| claim | verdict |
|---|---|
| T1′ (constant 2; attained at k = 1, 2) | **confirmed** as stated; see F3 for what it does not say. Witness p = 13, y² = x³+x+2, t = 3, k = 2 reproduced in the ordinary ring: dim V_D = 5, 25, 70 against dim I_D = 15, 35, 70 for D = 2, 3, 4, which is 1, 5, 16 against 11, 15, 16 modulo the chains; so 1 ∉ V_3 and LFD = 4 |
| P4 (filtered generation criterion) | **confirmed**; stated one degree weaker than used (F5) |
| P5 (counting) | **confirmed**; same remark (F5). Example k = 8, p = 257 recomputed independently (`rt_p5.py`): \|Z\| = 262 and 238; dim I_2 = 2; dim I_3 = 435 and 459; dim(I_2 + Σ v·I_2) = 34 = 17·2, generation fails at D′ = 2; dim(I_3 + Σ v·I_3) = 2255 and 2279 = dim I_4, generation holds at D′ = 3. D*(k) = 3, 5, 7, 8, 12, 16, 23, 30 for k = 8..128 confirmed; 2·H⁻¹(1/2) = 0.22006; k = 16 numbers 3205, 4951, 5372.4 confirmed. The 0.22·k lower bound is conditional, as the document says |
| E1 | empirical part consistent with a spot check (24 pairs, ordinary ring); "proof above the base" is conditional (F10) |
| E2 | arithmetic **confirmed** (f_4, f_5, Q_4, Q_5, all table entries, the four D = 5 profiles); wording **wrong** (F7) |
| C1 | conjecture; the equivalence with U_D, L_D for unit ideals is correct given P3 and the kernel fact 1 ∉ V_3 for 3t²+a ≠ 0; exceptional loci must be named (F9) |
| Corollary 3.2 | **confirmed** by reading (ε = indicator of H_k², weight 2) |
| §2.6 exact check | valid in the group algebra; not an independent test (F2) |
| §2.7 histograms | reproduced exactly for k = 1, 2 and for k = 3 at p = 17; for k ≥ 3 the script's values are conditional on T1 (F2) |
| §2.8 "slope exactly 1 for the selector method" | discussion, not a theorem; consistent with the tests: for a short block of depth 1 or 2 the selector bound is the exact value |
| §6.2 "(The true value is 5.)", §2.8 "(4, 5, 6 for k ≤ 13)" | **unsupported** as written (F6) |

## 7. Wording the paper may use

Safe statement:

> **Theorem (selector bound; written proof).** Let p ≥ 5 be prime, k ≥ 2 with 2^k ∣ p−1, and
> a, b, t ∈ F_p with 4a³+27b² ≠ 0. The normalized two-summand chain presentation
> {chain equations} ∪ {NF(S_3(x₁,x₂,t))} over F_p, for the subgroup base and for every coset
> base x_i = c·u_{i,0}, satisfies LFD ≤ k+2.

Safe description of its status:

> The proof is the fibre-selector argument of the binary case, carried out in the chain
> quotient F_p[(Z/2^k)²] with Hamming-weight degree. The passage to the ordinary polynomial ring
> is kernel-checked (`ChainClosureBridge.chain_mem_closure_iff`); the theorem itself is not
> formalized. It was confirmed by exact computation in the ordinary ring on 149,122 systems
> with 2 ≤ k ≤ 5, including every nonsingular (a, b, t) over F_5, F_13, F_17, F_29, F_37 at
> k = 2 and over F_17 at k = 3. The bound is
> attained at k = 2: there LFD = 4 whenever 3t²+a ≠ 0. For k ≥ 3 no computed system attains it
> (the largest exact values are 4 for k ≤ 5, 5 for k = 6..9 and 6 at k = 10), so it is not
> claimed to be sharp. The kernel-checked bound remains 2k+1, which needs no hypothesis on
> a, b, t.

In the two-sided table: keep "2k+1; mk+1 for all m (Lean)" under "Proved upper" and put
"≤ k+2, nonsingular curve (written)" in the "Written" column of the prime row, as the binary
row does with min(k,ℓ)+2.

Optional remark, labelled written proof with an unformalized bridge: for chains of depths
k₁ ≠ k₂ the same argument gives LFD ≤ min(k₁,k₂)+2, and exact computation finds equality
whenever min(k₁,k₂) ≤ 2 and 3t²+a ≠ 0.

Do not write: "replaces the dimension bound"; "kernel-checked" or "Lean" next to k+2; "for
every curve and target" or "for all a, b, t" without nonsingularity; "every prime" (say p ≥ 5);
"sharp", "tight" or "optimal" without "at k = 2"; "exact values 4, 5, 6 for k ≤ 13".

## 8. Cheapest falsification not ruled out

For T1 itself I have none to offer. Each instance of T1 reduces to one finite fact,
V_{k+2} = I_{k+2}, and every instance I could reach is certified. The universal statement rests
on the written proof, which has no step that depends on the size of k or p. What remains
untested:

1. The build and axiom status of the bridge modules at the reviewed commit. Test:
   `lake build ICFallDegrees.ChainClosureBridge ICFallDegrees.ChainTranslation` under the lock
   and the `#print axioms` lines 470 to 486 of `Audit.lean`. Minutes.
2. Unequal depths with both blocks deep were only sampled: exact LFD 4 on all 128 systems at
   depths (4,5), (4,6), (4,7), (5,6) and on all 8 at (3,8). Deeper pairs were not run.
3. Equal depths k ≥ 6 were not run here in the ordinary ring (one certificate at k = 6 needs the
   closure on 6,188 monomials; time). They are covered by the certified sweeps, whose values
   stay ≤ 6 through k = 10.

Outside T1: whether k+2 holds for singular (a, b) is open (F1); the cheapest decisive test would
be a singular instance with a vanishing fibre at larger k, where the measured LFD has been 3.

## 9. Reproduction

Scripts, summaries and console outputs: `experiments/ic-bounds-upper-redteam/` (README there).

    cd experiments/ic-bounds-upper-redteam
    ulimit -v 1500000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
    python3 rt_sanity.py; python3 rt_sweep_cross.py; python3 rt_control.py
    for j in eq12 eq3 eq45 uneq17 uneqdeep uneqmid raw; do python3 rt_run.py $j; done
    python3 rt_univariate.py; python3 rt_singular.py; python3 rt_p5.py; python3 rt_e1.py
    python3 -I sweep_compare.py

Environment: Python 3.13.18, numpy 2.5.3, one thread, under 1.5 GB. Per-instance rows
(`rt-*.jsonl`, about 45 MB) were not kept; the summaries list every family with its histogram.

## 10. Not reviewed

- P3 and its Lean modules (kernel; excluded by the dispatch). `rt_e1.py` only spot-checks the
  statement on 24 pairs.
- The numbers of Section 5.2 of the document other than E1 (tests T2′, T3, T4, T5, T6), the
  files `lift_of_one.py`, `sweep_monotone.py`, `depth_restriction.py`, and the attribution
  "Known: U_5, L_5 (exact, k = 9, 10)" in Section 5.1.
- The interpretive parts of Sections 5.3, 6.1 and 6.3, and Remark 2.9 (m ≥ 3), which the
  document itself marks as not proved.
- The producer's scripts were read, not rerun.
