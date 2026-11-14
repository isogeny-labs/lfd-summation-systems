# Red-team review 2 of the growth theorem (N8), commit-pinned (2026-11-12)

**Pinned artifacts (hashes checked before the review; working tree identical to the commit for every file below).**

* Commit `09dd976093937c4e5a1ed8aa10c81f87d06ad2d1`.
* `docs/proofs/N8-growth.md`: blob `0e58d23c75993973bffa1a0545bb10488790f281`.
* `docs/proofs/N8-lean-plan.md`: blob `a0326459a36bece22a31d730857727a34599faf2`.
* `paper/curve-lfd-bounds.tex`: blob `4538941b48a0c03bb97d0a7e2baf729f8eb38266`. Read: `thm:growth` and the
  paragraphs after it (l. 957–1041), plus the abstract (l. 527–532) and the introduction (l. 620–627), because
  point 3 asks whether the paper implies more anywhere.
* Lean files at the commit: `GrowthCounting` `5ff34744`, `GrowthReduction` `e33fcc3f`, `GrowthSpread` `7dddb90d`,
  `GrowthMain` `4e1af514`, `GrowthCriterion` `33ccb4da`, `GrowthExistence` `bbdaf4e7`. Background, read only:
  `Closure` `0d6a5444`, `PrimeM2General` `1dd75165`, `ChainQuotient` `d9d9653b`, `ChainClosureBridge` `49203fa6`,
  `UnitInstance` `632176ed`, `SelectorExtension` `ba384e50`.

**Reviewer and scope.** I am a breakthrough-tier red-team session, separate from the producer, the coordinator and
the first reviewer. My brief was to find errors, not to confirm the proof. Labels used: **kernel**, **proof**,
**exact**, **empirical**, **conjecture**. Every label below is a recommendation; I do not mark anything as
supported, closed or rejected. I edited no existing file, made no commit, and changed no status. My new files are
this review and `experiments/ic-bounds-growth/review2/` (four scripts with their outputs). Lean was used only
through the lake lock: two `--no-build` checks and two scratch `#print axioms` files kept outside the repository.
I read the earlier review (`docs/reviews-N8-growth-redteam-20261112.md`, blob `6c98619e`) only after forming my
own view (§13).

---

## 1. Verdict

**No fatal error to T1 and no fatal error to a secondary claim.** I re-derived the characteristic-p Wronskian
argument, Lemma 3.2, Theorem 3.3, Lemma 3.4, Lemmas 4.1–4.4, Theorem 5.1, Propositions 2.1/2.2/7.1/7.2 and
Corollaries 5.2/5.3 by hand. Every step holds. Numerical checks with my own code found no counterexample to any
finite statement or to Theorem 3.3 over F^p.

The findings are twelve presentation problems and one minor repairable gap (§2). They have three sources:
stale text that predates the Lean development at the same commit; wording in the paper's abstract, introduction
and proof sketch; and the scope of one **exact** label.

| claim | label in N8 (blob 0e58d23c) | recommended label |
|---|---|---|
| T1 (Thm 5.1) | proof | **proof** (unchanged). Also: the reduction and counting are **kernel, conditional on `GrowthMain.CurveModel`** (`GrowthMain.growth_theorem`). `CurveModel` itself is **proof** (N8 §1, §3) and is not in Lean. |
| C1 (Cor 5.2) | proof given T1; bridge, upper bound and unit lemma kernel | (a), (b): **proof** given T1, and **kernel conditional on `CurveModel`** (`growth_lfd_gt`, `growth_lfd_gt_linear`, `growth_lfd_window`). (c), existence of unit instances at every k ≥ 2 over a prime factor of 2^(2^k)+1: **kernel, unconditional** (`GrowthExistence.fermat_factor_dvd`, `unit_instance`). "Unbounded": **proof** given T1. |
| T2 (Thm 3.1) | proof | **proof**. No novelty is claimed. The characteristic-p attribution to Voloch is unverified; the proof written in the file is correct. |
| T3 (Props 2.1, 2.2) | proof | **kernel**, unconditional for any prime 𝔭 ∋ S_3 with x_i, x_i^N − 1 ∉ 𝔭 (`GrowthCriterion.closed_iff_fei`, `GrowthReduction.not_one_of_fei`, `closure_normalForm_mem_idealSpace`). The file understates it. |
| T4 (Cor 5.3) | proof | **proof** (checked; not in Lean) |
| Prop 7.1 | proof | **kernel** (`GrowthCriterion.not_one_of_crit`, `crit_of_fei`). The file understates it. |
| Prop 7.2 | proof | **proof** (checked) |
| E1, E3, E4, R1 | exact / empirical | unchanged (not re-run) |
| E2 | exact | **exact on the computed range only**: D = 3, k ≤ 7; D = 4, k ≤ 8 (finding R2-9) |

**Is the artifact set sufficient at this tier?** Yes for review of T1 at the **proof** level, and of the Lean
development at the **kernel-conditional** level. Every step is written out. The Lean modules are present,
imported in `formal/ICFallDegrees.lean` (l. 111–116), listed in `Audit.lean` (l. 556–590), up to date and
dependent only on standard axioms. It is **not** a kernel closure of T1: the analytic hypothesis `CurveModel` (the
function field of E and Theorem 3.3) is not proved in Lean. That is an explicit, named gap, not a hidden one.
T1 is a mathematical lower bound, not a cryptanalytic break, so problem size, ICV1 curve, model of computation
and cost accounting do not apply.

---

## 2. Findings

Severity scale: fatal to T1 / fatal to a secondary claim / repairable gap / presentation.

| id | location | severity | reason | repair | cheapest test that would overturn this finding |
|---|---|---|---|---|---|
| R2-1 | N8 l. 10–11 ("Nothing in this file is kernel-checked"); §8 l. 700 ("No Lean file was written; nothing is to be added to `ICFallDegrees.lean` or the audit"); Cor 5.2(a) proof l. 353–354 and §6.3 l. 574–578 ("The translation of S_D(k) into such a W … has not been formalized for this W") | presentation (stale) | At the same commit, `GrowthReduction` instantiates the bridge with W = S_D(k) (§10). The six Growth modules are imported and audited. T3, Prop 7.1, Lemma 3.4, Lemma 4.1 (except Λ(w) ≥ 4w, which Thm 5.1 does not use), Lemmas 4.2–4.4 and Thm 5.1-from-Thm 3.3 are kernel-checked. The proof file says the opposite. | Replace §8 with a pointer to `N8-lean-plan.md`. Relabel T3 and Prop 7.1 **kernel**. Add "kernel conditional on `CurveModel`" to T1, C1(a)(b). Delete the two "not formalized" sentences. | `grep -n Growth formal/ICFallDegrees.lean formal/ICFallDegrees/Audit.lean` |
| R2-2 | `N8-lean-plan.md` l. 6–8 | presentation (stale) | The note says the modules are "not yet imported … or listed in `Audit.lean`". They are both, at this commit. | Update the status line. | same grep |
| R2-3 | `docs/IC_BOUNDS_PROGRAM.md` l. 665 (changed in this commit) | presentation (pinning) | It names the revised proof file as blob `a775f94e…`. No such object exists in the repository. The revised blob is `0e58d23c`. A commit-pinned review cannot follow a pin that resolves to nothing. | Correct the hash. | `git cat-file -t a775f94e` (fails); `git rev-parse 09dd9760:docs/proofs/N8-growth.md` |
| R2-4 | paper abstract l. 527–529; introduction l. 623–626 | presentation (evidence status) | "kernel-checked, with a published bound on vanishing sums in function fields as an explicit hypothesis". The Lean hypothesis is not a published bound. It is `CurveModel`/`SpreadBound`, which is N8's Theorem 3.3: Theorem 3.1 plus the curve-specific Lemma 3.2 and step (3.1), a written proof. The published input exists only as Theorem 3.1's type: Brownawell–Masser in characteristic 0 (quoted second-hand) and "of the type of Voloch" in characteristic p (statement not retrieved). Every existence instance is in characteristic p. The theorem environment and l. 1003–1010 state the status correctly. | Write "with the vanishing-sum bound for this curve (Theorem 3.3 of the proof file; a written proof of Brownawell–Masser/Voloch type) as an explicit hypothesis". | Read `GrowthSpread.SpreadBound`: it quantifies over circuits of the monomials ξ₁^a ξ₂^b of one model. No published theorem has that form. |
| R2-5 | paper abstract l. 530–532 | presentation (scope) | The abstract states the bound for "unit-ideal instances with nonzero target abscissa". It drops the hypothesis char 0 or p > 2^(k+1) of `thm:growth`, and nonsingularity is left implicit. Under 2^k ∣ p − 1 the hypothesis fails only for p = 2^k + 1. For k ≥ 38 that case needs a prime Fermat number F_m with m ≥ 33. To my knowledge F_33 is the first Fermat number of unknown status, so the case is not formally excluded at k = 2^33. Over general fields (the paper's other theorems allow any characteristic ≠ 2) the omission is real. | Add "(characteristic 0 or > 2^(k+1))" to the abstract. | Exhibit a unit instance with t ≠ 0 and 2^k < p < 2^(k+1) where the abstract's bound is claimed. T1 does not cover it (N8 Remark (i): open). |
| R2-6 | paper sketch l. 986–987 ("t≠0 keeps the two coordinate functions from sharing a pole"); N8 §0 l. 29–30 | presentation (wrong reason) | The poles of ξ₁, ξ₂ are at O and R. These differ for every t, because R ≠ O. What t ≠ 0 does is make ξ₂ a unit at O and ξ₁ a unit at R. At t = 0, ξ₂ vanishes at O and ξ₁ at R, to order 1, or order 2 when b = 0. I verified this numerically (§8). N8 §0's first half-sentence ("poles at different points") is automatic; its second half is the real content. | "t ≠ 0 makes each coordinate function a unit at the pole of the other". | §8 table, rows t = 0. |
| R2-7 | paper sketch l. 989–991 ("with differences taken modulo p") | presentation (imprecise) | What is proved is ‖2Δa‖_p + ‖2Δb‖_p ≤ 6(m−1)(m−2), i.e. *twice* the differences modulo p. The bound ‖Δ‖_p ≤ 3(m−1)(m−2) does not follow: for Δ = (p−1)/2, ‖2Δ‖_p = 1 but ‖Δ‖_p = (p−1)/2. | "twice the exponent differences, modulo p". | Read N8 Thm 3.3 and its proof: only ord ≡ −2Δ (mod p) is used. |
| R2-8 | paper l. 1019–1021 | presentation (nit) | "nonsingular and t≠0 … do not appear in the Lean statements". `nonsingular` is a hypothesis of `GrowthMain.growth_lfd_window` and is discharged in `GrowthExistence.growth_instance_window` (upper bound, `instance_nonsingular`). Also, the lower half of "⌊(k−10)/14⌋ < d … for k ≥ 24" is trivial for 24 ≤ k < 38 (it says d > 1). | "do not appear in the Lean lower-bound statements"; mention that the bound is non-trivial from k = 38. | `#check @GrowthMain.growth_lfd_window` |
| R2-9 | N8 E2 row l. 23; §7.1 l. 603 | presentation (scope of an **exact** label) | "Crit_D(k) holds exactly for k ≥ 4D − 10" / "Crit_D(k) holds if and only if k ≥ 4D − 10", labelled exact. The computed range is D = 3, k ≤ 7 and D = 4, k ≤ 8 (Section 6.2 tables). Read literally, as a statement for every k, it is not established: nothing covers D = 4, 9 ≤ k < 52. The read-off table in §6.2 states the range correctly. | Add the range to E2 and to §7.1. | none needed (it is a wording fix) |
| R2-10 | N8 T1 row l. 17, second sentence | repairable gap (minor, secondary sentence) | "every element of V_D(k) … vanishes on the curve S_3(x₁, x₂, t) = 0". The proof gives f̃ ∈ 𝔭, i.e. vanishing on the image of P ↦ (x(P), x(R−P)). Equating that image with the zero set of S_3(·,·,t) needs the converse half of Semaev's theorem: every zero comes from points with P₁ ± P₂ ± R = O. The file cites only the forward half (l. 75–77) and says "no irreducibility of S_3 is needed". The statement is true by the converse. The Lean counterpart, `growth_closure_vanishes`, states vanishing at (ξ₁, ξ₂), which is exactly f̃ ∈ 𝔭. | Cite the converse direction, or write "lies in 𝔭 (vanishes on the image curve of E)". | A zero of S_3(·,·,t) not of the form (x(P), x(R−P)). There is none, by Semaev; this only checks the citation. |
| R2-11 | N8 Remark 1 l. 405–407 (unequal depths) | presentation (nit) | Lemma 4.3 bounds windows in both coordinates. With k₁ ≠ k₂ it needs p > 2^(max(k₁,k₂)+1), which the remark does not state. It is an unlabelled remark. | Add the condition. | — |
| R2-12 | N8 Thm 3.1 proof l. 151–153, 171 | presentation (nit) | Two steps are left to the reader. (i) "ker δ … a field between F^p and F, different from F, and [F : F^p] = p" leaves "hence ker δ = F^p" implicit. (ii) "δ_t is the formal derivative" on k̄((t)) is asserted without reference. Both are standard and correct. | One clause and one citation (e.g. Stichtenoth, differentials). | — |
| R2-13 | N8 l. 7–9, 10–11 vs. the commit message of 09dd9760 | presentation (process) | The file says the earlier review "applies to the earlier blob … would need a re-review". This review is that re-review, for blob 0e58d23c. Its verdict should not be carried to later blobs either. | Record this review's pin next to the labels. | — |

None of these changes Λ, k_D, the thresholds, or any statement of T1, C1 (a)–(c), T4, Prop 7.1 or Prop 7.2.

---

## 3. Point 1: the characteristic-p Wronskian argument (Thm 3.1/3.3) and Lemma 3.4

Checked line by line by hand. Numerical support is in §9 and in `fp_circuits.py`.

* **Hasse vs. ordinary derivatives.** The file uses ordinary derivatives δ^i, δ = d/dz, together with independence
  over ker δ. In characteristic p, ker δ = F^p: it contains F^p, it misses z, and [F : F^p] = p is prime. This
  pairing is the correct one. A family independent over F^p has at most p members, so r = m − 1 ≤ p and only
  δ^0..δ^(p−1) occur; δ^p = 0 on F is never reached. Hasse derivatives would be needed only for independence over
  the constants k̄. That statement is false here, and the file uses the Frobenius control to say so (l. 234–238).
* **(a) Wronskian criterion.** The induction holds in any characteristic. A zero determinant gives a column
  dependence with first coefficient 1. Differentiating row i and subtracting row i+1 gives an (r−1)-system with
  determinant W(w₂..w_r) ≠ 0, so δ(c_j) = 0 and the c_j lie in ker δ.
* **Separating variable and change of derivation.** z exists because k̄ is perfect. A uniformizer t_P has order 1,
  not divisible by p, so it is separating. Der_k̄(F) is one-dimensional over F, so δ = δ(t)·δ_t, and
  W_δ = c^(r(r−1)/2)·W_{δ_t} with ord_P c = −ord_P(dz/dt_P).
* **(c) Local bounds.** (3.2) is the termwise estimate over permutations, with loss Σπ(j) = r(r−1)/2. (3.3) holds
  because δ_t(t^(pe)) = 0, so δ_t^l(t^(pe) y) = t^(pe) δ_t^l(y). This works for any e ∈ Z, including the zeros and
  poles of non-constant F^p coefficients (§9).
* **(d) Orders at every pole and zero.** W = ±W(v_j : j ≠ i) for every i, by multilinearity and Σv = 0. At every P
  one chooses i_P of minimal order. At P ∈ 𝒮 the loss is ≤ n₀. At P ∉ 𝒮 there is no loss: hypothesis 3 gives
  orders ≡ 0 mod p there, including any zeros and poles of the φ_e. The −n₀·ord_P(dz) term covers the ramification
  of z.
* **Genus and divisor terms.** Σ_P ord_P(W) = 0, Σ_P ord_P(dz) = 2γ − 2 = 0 for E (z separating), Σ_P A_P = 0, and
  −Σμ_P = H(v). Hence H ≤ n₀(2γ − 2 + |𝒮|).
* **Constant 6(m−1)(m−2).** n₀|𝒮| ≤ 6·(m−1)(m−2)/2 = 3(m−1)(m−2). (3.1) gives 2H ≥ |ord_O(v_e/v_e′)| + |ord_R(v_e/v_e′)|.
  That needs O ≠ R, which always holds. The orders are ≡ −2Δa and −2Δb (mod p), since ord(φ_e) ≡ 0. So
  ‖2Δa‖_p + ‖2Δb‖_p ≤ 6(m−1)(m−2). In characteristic 0 this is Brownawell–Masser's pairwise bound with k = m − 1,
  g = 1, |S| ≤ 6 (the derivation is re-done in N8 l. 136–141; I agree with it).
* **Use of t ≠ 0.** It is used exactly for ord_O(ξ₂) = 0 and ord_R(ξ₁) = 0 (Lemma 3.2). §8 shows that it is the
  only exclusion needed.
* **Lemma 3.4 (F^p-circuits vs. K-coefficients).** The lemma needs only F₀-linear algebra inside F and K ⊂ k̄ ⊂ F₀
  (k̄ is perfect, so k̄ ⊂ F^p). It is applied to **one** polynomial f = Σc_e x^e with c_e ∈ K. The value
  y = Σ_{M₁} c_e ξ^e lies in two F₀-spans. If y ≠ 0, a circuit would meet both parts. So y = 0. The M₁-part is a
  polynomial with coefficients in K whose image in F is 0, so it lies in 𝔭 = ker(K[x] → F). No F^p-linear
  structure on 𝔭 is used, and 𝔭 is never extended to F₀[x]. The Lean version (`split_of_no_crossing`, `cf`,
  `linearCombination_cf`, `cf_topPart` in `fei_of_spread`) does exactly this through `IsScalarTower K F₀ F`. It is
  kernel-checked with standard axioms.

**Verdict, point 1: correct.** No gap beyond the nits R2-12 and R2-6/R2-7, which concern the paper's sketch.

## 4. Point 2: exact k_D and the characteristic quantifiers

* **Λ and k_D** were recomputed exactly (§6): Λ(1..5) = 25, 36, 48, 59, 71, and k_2..k_6 = 27, 39, 52, 64, 77. The
  whole N8 table (D ≤ 10, 12, 16, 20) matches. k_D ≤ 14D + 10 holds for D ≤ 121.
* **Every use of "char 0 or p > 2^(k+1)" (p ≥ 2N + 1):**
  - T1, C1, Thm 5.1, Cor 5.2: as hypotheses.
  - Lemmas 4.2–4.4: p ≥ 2N + 1, the same condition.
  - Thm 3.1/3.3 and Lemma 3.4 need no bound on p beyond p > 3.
  - Cor 5.3 needs p > 2^(k+2), stated (Lemma 4.3 at depth k + 1).
  - Remark (i), 2^k < p < 2^(k+1), is correctly called open.
  - Only Remark 1 lacks its condition (R2-11).
* **Existence instances.** If p ∣ 2^(2^k) + 1, then 2^(2^k) ≡ −1, so 2 has order exactly 2^(k+1), so
  2^(k+1) ∣ p − 1 and p ≥ 2^(k+1) + 1 > 2^(k+1). This holds for every k ≥ 2, and is kernel-checked in
  `GrowthExistence.fermat_factor_dvd`.
* **Hypotheses of the instance** (X, n non-squares, a = n − X², b = 0, t = X):
  - p is odd;
  - t ≠ 0, since 0 is a square;
  - a ≠ 0, so 4a³ + 27b² = 4a³ ≠ 0;
  - t³ + at = nX is a product of two non-squares, hence a nonzero square, so Y ≠ 0 exists;
  - t is not a square;
  - 2^(k+1) ∣ p − 1.

  These are exactly the hypotheses of `UnitInstance.unit_zmod` (printed: p prime, 2 ≤ k, 2^(k+1) ∣ p − 1,
  Y² = t³ + at, Y ≠ 0, ¬IsSquare t) and of `GrowthExistence.unit_instance`. They also satisfy T1: R = (X, Y) is
  not 2-torsion and t ≠ 0. b = 0 is allowed; §8 checks the divisors for b = 0.
* The Lean thresholds use p = 0 ∨ 2·2^k + 1 ≤ p, which is the same condition.

**Verdict, point 2: correct.** The abstract omits the quantifier (R2-5).

## 5. Point 3: FEI and the Lean hypothesis `GrowthMain.CurveModel`

The statements below were printed with `pp.coercions` from the built modules.

* **What the written proof supplies.**
  - A tower K ⊂ k̄ ⊂ F₀ ⊂ F, with F = k̄(E) and F₀ = k̄ (char 0) or F^p (char p). It satisfies `IsScalarTower`.
  - ξ = (x(P), x(R−P)). Both entries are non-constant, hence ≠ 0 (`ne_zero`).
  - S3at a b t ξ = 0 by Semaev's identity (`on_curve`). The Lean `S3at` is the N8 formula term for term.
  - `SpreadBound F₀ ξ p` with p = char K.
* **Quantifiers of `SpreadBound`.** It ranges over every `J : Finset (Fin 2 →₀ ℕ)`, every circuit, and all e, e′ ∈ J,
  and asks for ∃ z₁ z₂ ∈ ℤ with |2(e₀ − e′₀) − p z₁| + |2(e₁ − e′₁) − p z₂| ≤ 6(|J|−1)(|J|−2), all in ℤ.
  - This is N8 Thm 3.3 for all finite J ⊂ Z²_{≥0}, with ‖n‖_p = min_z |n − pz|; the two minima are independent.
    For p = 0 it reduces to |·|.
  - `Circuit` is "not `LinearIndepOn` on J, and `LinearIndepOn` on every J′ ⊂ J" for the family `mono ξ` indexed by
    exponents. So two exponents with proportional monomials form a 2-circuit, as in N8.
  - |J| ≥ 2 is automatic, since mono ≠ 0. `mono ξ e = ξ 0 ^ e 0 * ξ 1 ^ e 1`, with ξ 0 = x(P), whose pole is at O,
    matching the first coordinate in Thm 3.3's proof.
* **Is p tied to the characteristic?** Not in `growth_theorem`, `growth_not_one` or `fei_of_spread`. This is not a
  gap: the Lean proof is combinatorial and needs only `SpreadBound … p` and p = 0 ∨ p ≥ 2N + 1 for the *same* p.
  The written proof supplies both with p = char K. The tie is made by construction in `growth_lfd_gt_zmod` and
  `growth_instance` (ZMod p and `CurveModel … p`).
  - The hypothesis cannot be satisfied trivially. With F₀ ⊇ F^p in characteristic p, {(0,0),(p,0)} is a circuit
    with bound 0. So `SpreadBound` fails for every p′ ∤ 2p, and p′ ∈ {1, 2} is excluded by hp.
* **Status in the paper.**
  - Correct: the theorem label "[Growth; written proof]", and l. 1003–1010 ("stated here unconditionally, as a
    written proof; in Lean it is conditional … explicit hypothesis … not an axiom, and not proved there").
  - Overstated: the abstract and the introduction ("published bound … as an explicit hypothesis", R2-4), and the
    abstract's missing characteristic hypothesis (R2-5).
  - Minor: l. 1019–1021 (R2-8).
  - Elsewhere (l. 1756) the paper says only "written proof".

**Verdict, point 3: the translation is faithful.** The paper's precise statement is right; two summary sentences
imply more than is shown.

## 6. Point 4: the numerical Λ bound (own checker)

`experiments/ic-bounds-growth/review2/lambda_check.py` uses its own code and imports nothing from the repository.

* C(n, ≤w) comes from Pascal rows, cross-checked against `math.comb`.
* Λ(w) is the largest ℓ ≥ 1 with 2^(ℓ−1) < 7776·C(2ℓ, ≤w)², computed with exact integers for ℓ < L₀(w).
* The tail ℓ ≥ L₀ is handled rigorously by an independent argument: C(n, ≤w) ≤ (n+1)^w, and
  7776(2ℓ+1)^(2w)/2^(ℓ−1) is non-increasing from L₀ on.

Results for w ≤ 120 (`lambda_check_w120.json`):

* Λ(1..5) = 25, 36, 48, 59, 71. Λ(w) for w = 1..60: 25, 36, 48, 59, 71, 83, 95, 106, 118, 130, 142, …, 714, 726.
* The full N8 table (Λ(D−1) and k_D for D = 2–10, 12, 16, 20) matches.
* 4w + 21 ≤ Λ(w) ≤ 13w + 12 for 1 ≤ w ≤ 120. So Λ(w) ≤ 13w + 23 holds with slack 11, and k_D ≤ 14D + 10 holds for
  D ≤ 121. Every ℓ ≤ Λ(w) satisfies the inequality (no gaps).
* ρ is strictly increasing on [4w, L₀]. The satisfying set in [4w, ∞) is an initial segment.
* The exact check 7776·C(2ℓ, ≤w)² ≤ 2^(ℓ−1) at ℓ = 13w + 24 holds for every w checked.
* The linear form ⌊(k−10)/14⌋ + 1 never exceeds 1 + max{D : k_D ≤ k} for 38 ≤ k ≤ 1700. LFD ≥ 8 first holds at
  k = 90.
* Constants in the proof of Lemma 4.1:
  - differences for w = 1..7: 7.77, 7.73, 8.06, 8.55, 9.11, 9.71, 10.34; N8 quotes truncations, which is
    conservative;
  - 2·log₂(5.437·16) = 12.886 < 12.89;
  - (8/7)^4 = 1.706 < 2, and 2(7/8)^4 = 1.172 > 1.17;
  - c = 11.9925, the root of 4H₂(1/(2c)) = 1, so the slope is 1/12.99;
  - k_D/D ∈ [12.898, 12.914] for D = 55..61.
* Lean's numeric steps hold: 17^52·16^2 ≤ 2^13·16^52, 17^4 ≤ 2·16^4, and the base case
  15552·17^96 ≤ 2^24·16^96. The Chernoff estimate C(n, ≤w)·16^n ≤ 17^n·16^w holds for every n, w in range.
* **Lean definitions** (statements printed):
  - `Csum n w = Σ_{i<w+1} n.choose i`.
  - `LambdaBound w lam := ∀ ℓ, lam < ℓ → 7776·Csum(2ℓ) w² ≤ 2^(ℓ−1)`, i.e. "the inequality fails above lam".
  - `Lam w := Nat.findGreatest (2^(ℓ−1) < 7776·Csum(2ℓ) w²) (13w+23)`. With `lambdaBound_linear`
    (`LambdaBound w (13w+23)`) and `lambdaBound_Lam`/`Lam_greatest` this is the text's Λ(w).
  - `Lam_one … Lam_five : Lam w = 25, 36, 48, 59, 71` via `Lam_eq_of`: the inequality holds at n and fails on
    (n, 13w+23], evaluated by `decide +kernel` on a kernel-evaluable `csumFast` proved equal to `Csum`.
  - All of these have standard axioms only.

**Verdict, point 4: confirmed.**

## 7. Point 5: the no-wrap lift identity (1.1) at the top-run boundary

* **By hand, on a monomial x^a.** If a_i + 2^j < N, both sides equal x^(a + 2^j e_i). If a_i + 2^j ≥ N, which
  includes the boundary a_i = N − 2^j, then (v x^a)~ = x^(a + (2^j − N)e_i). The right side is
  x_i^(2^j)x^a − x_i^(2^j)x^a + x_i^(2^j−N)x^a, which is the same. At the boundary the i-exponent becomes 0. The
  inequality defining f̃_hi is "≥", which is correct.
* a_i + 2^j ≥ 2N is impossible, since a_i ≤ N − 1 and 2^j ≤ N/2.
* "Top run" (bits j..k−1 all one) ⟺ a_i ≥ N − 2^j ⟺ a_i + 2^j ≥ N, for 0 ≤ a_i < N.
* **Brute force** (`lift_check.py`, own code, coefficients in GF(1000003)):
  - every monomial in [0, N)², both i and every j < k, for k = 2..6 (61,888 cases);
  - 8,040 random and boundary-heavy combinations, with a_i ∈ {N−2^j−1, N−2^j, N−2^j+1, N−1, 0} for every
    second exponent;
  - both the written form (Laurent) and the Lean form x_i^(N−2^j)·(vf)~ = x_i^N(f̃ − f̃_hi) + f̃_hi, which has
    non-negative exponents;
  - the divisibility x_i^(N−2^j) ∣ f̃_hi: **0 failures**;
  - the top-run equivalence, no double wrap, and wt((a + 2^j) mod N) ≤ wt(a) + 1, for every a < 2^k with
    k ≤ 14: **0 failures**.
* **Further brute force (own code).**
  - Lemma 4.2: k ≤ 7, p = 0 and the four smallest primes ≥ 2N+1, all σ ≤ N+1, all a₀. At most 2 blocks occur
    (the bound is 6).
  - Lemma 4.4: k ≤ 9, every λ, σ < 2^(λ−1), j ≥ λ−1, every straddling pair, with the smallest admissible p (the
    weakest case). 0 failures.
  - The count of Lemma 4.3 (pairs of weight ≤ w in one pair of aligned blocks ≤ C(2ℓ, ≤w)): k ≤ 6, 0 failures.
* **Lean `GrowthReduction.lift_identity`** (printed) is (1.1) multiplied by x_i^(N−2^j). It holds for every reduced
  chain polynomial and every j : Fin k. `ChainTranslation.normalForm_X_mul_monomial` shows that the reduced product
  is the cyclic shift (b_i + 2^j) mod 2^k, i.e. multiplication in A_k.

**Verdict, point 5: correct.**

---

## 8. Addendum (a): Lemma 3.2, the divisors of ξ₁ and ξ₂

**By hand.**
* div(ξ₁) = div(x) = (P₀) + (−P₀) − 2(O), with x(±P₀) = 0. When b = 0 this is 2(P₀) − 2(O), since P₀ = (0,0) is
  2-torsion.
* τ: P ↦ R − P is an involutive automorphism of the curve, so div(ξ₂) = τ(div x) = (R − P₀) + (R + P₀) − 2(R), with
  orders preserved. This holds for every R ≠ O, including R 2-torsion (t a root of x³ + ax + b).
* Coincidences:
  - R ∓ P₀ = O ⟺ R = ±P₀ ⟺ t = 0;
  - R ± P₀ = R and ±P₀ = O are impossible;
  - R ± P₀ = ∓P₀ ⟺ R = ∓2P₀ (t = a²/(4b) ≠ 0): a zero of ξ₂ coincides with a zero of ξ₁;
  - R + P₀ = R − P₀ ⟺ b = 0;
  - R = ±R is irrelevant.
* ord_O(ξ₂) = 0 ⟺ O ∉ {R, R ± P₀} ⟺ t ≠ 0. ord_R(ξ₁) = 0 ⟺ R ∉ {O, ±P₀} ⟺ t ≠ 0.
* **So t ≠ 0 excludes exactly R = ±P₀.** No other case needs to be excluded. The other coincidences (b = 0,
  R 2-torsion, R = ±2P₀, R 2-torsion with b = 0) only shrink 𝒮. They do not touch the orders at O and R used by
  Thm 3.3, and |𝒮| ≤ 6 still holds.

**Numerically.** `divisor_check.py lemma32` is my own Laurent-series engine over F_{p²}:
* uniformizer x − x₀ at points with y₀ ≠ 0, y at points with y₀ = 0, and z = −x/y at O (Silverman's w(z));
* ξ₂ is computed from the chord formula in the function field;
* the full divisors are obtained by scanning all rational places. Completeness: all predicted points are
  F_{p²}-rational, and the zeros total 2 and the poles total 2.

Two curves per configuration (p = 13 and p = 17), all 14 runs matching the prediction exactly:

| configuration | (p, a, b, t) | div(ξ₂) found | ord_O ξ₂ | ord_R ξ₁ | |𝒮| | Lemma 3.2 conclusion (ord_O = −2a, ord_R = −2b on 4 monomials) |
|---|---|---|---|---|---|---|
| generic | (13,1,1,1), (17,1,1,1) | (R∓P₀) − 2(R) | 0 | 0 | 6 | holds |
| b = 0 | (13,1,0,1), (17,1,0,1) | 2(R+P₀) − 2(R) | 0 | 0 | 4 | holds |
| R 2-torsion | (13,0,1,4), (17,0,1,16) | (R∓P₀) − 2(R) | 0 | 0 | 6 | holds |
| R 2-torsion and b = 0 | (13,1,0,5), (17,1,0,4) | 2(R+P₀) − 2(R) | 0 | 0 | 4 | holds |
| R = ±2P₀ | (13,1,1,10), (17,1,1,13) | (−P₀ or P₀) + (other) − 2(R) | 0 | 0 | 5 | holds |
| t = 0 (excluded) | (13,0,1,0), (17,0,1,0) | (O) + (·) − 2(R) | **1** | **1** | 3 | **fails** |
| t = 0, b = 0 (excluded) | (13,1,0,0), (17,1,0,0) | 2(O) − 2(R) | **2** | **2** | 2 | **fails** |

The conclusion holds in exactly the runs with t ≠ 0 (`lemma32.json`).

## 9. Addendum (b): non-constant F^p coefficients in Theorem 3.1

**Argument.** In Thm 3.3 the coefficients φ_e ∈ F^p can be non-constant, with zeros and poles anywhere. Their
orders are p·ord(ψ_e), where φ_e = ψ_e^p.
* At places outside 𝒮, hypothesis 3 is exactly "orders ≡ 0 mod p", not "units". (3.3) shows no loss there, because
  t^(pe) is a δ_t-constant for any e ∈ Z.
* At places of 𝒮, (3.2) holds for arbitrary orders.
* The φ-orders enter H(v) and the A_P, and they cancel consistently in the summation.
* At O and R, ord(φ_e/φ_e′) ≡ 0 mod p. This is why the conclusion is about ‖2Δ‖_p and not |2Δ|.

The argument treats the non-constant case correctly.

**Test with a relation whose coefficients are non-constant p-th powers** (`divisor_check.py twisted`). The test
case:
* curves (p, a, b, t) = (11, 1, 3, 2) and (13, 1, 0, 2); the second is an existence-type b = 0 curve;
* monomials ξ^(m_e), m_e = e + p·u_e, e ∈ {0,1,2}², with twists u_e up to 9;
* the relation Σ_e φ_e ξ^(m_e) = 0 with φ_e = c_e·h^p·ξ^(−p·u_e), where c_e are the S_3 coefficients;
* h = (x₁ − α)/(x₂ − β) has zeros and poles outside 𝒮, e.g. zeros at (1,4), (1,7) and a double pole at (8, 5w)
  for p = 11.

Results, all exact over F_{p²}:
* The relation vanishes identically as a Laurent series at every rational place.
* {ξ^(m_e)} is F^p-dependent: the Wronskian with ordinary derivatives vanishes at 10 points. Every 8 of the 9 are
  independent. The rank over the constants is 9, so the relation **needs** non-constant coefficients.
* ord(φ_e) ≡ 0 (mod p) at every place, including the zeros and poles of h outside 𝒮.
* H(v) = 8 ≤ n₀|𝒮| (168, resp. 112). |ord_O| + |ord_R| ≤ 2H for every pair. ord_O(v_e/v_e′) ≡ −2Δ (mod p).
* Theorem 3.3's mod-p spread is at most 8 ≤ 336. The integer spread is 398 and 470, above 336. So the
  characteristic-0 form of the bound is **false** for F^p-circuits and the mod-p form is needed. It holds.
* **The local Wronskian inequality of step (d)**, ord_P W ≥ −n₀·ord_P(dx) + (A_P − μ_P) − n₀[P ∈ 𝒮], was computed
  from Laurent series by DVR elimination of the 8×8 Wronskian matrix, with δ = d/dx₁. It holds at all places
  tested: 13 for p = 11 and 11 for p = 13.
  - At the zeros and poles of h there is no loss outside 𝒮, as (3.3) says. Equality holds at 7 of the 8 such
    places (ord W = 8p, or 8·(−2p) at the double pole of h); at (11,4) for p = 13, ord W = −103 against the
    bound −104.
  - At the ramification points of x it holds with equality, ord W = −28 = −n₀·ord(dx).
  - At the places of 𝒮 the actual loss is 9 or 18, below n₀ = 28. For example, at O for p = 11: ord W = −124
    against the bound −134.

**Further F^p tests** (`fp_circuits.py`, own local-expansion Wronskian, residues mod p):
* curves p = 61 (two curves, one with b = 0), p = 101 (two curves, one with b = 0) and p = 131;
* m = 2: no residue e ≠ 0 has ξ^e ∈ F^p;
* m = 3, exhaustive over all pairs of residues (p = 61 twice, p = 101): no three-term F^p-circuit at all; every
  pair is eliminated by a point with W ≠ 0;
* m = 4: 60 random bases {0, e₁, e₂} with pairwise spread > 36, each scanned against every residue (up to 1.03
  million extensions per curve): no dependent extension;
* control: the nine S_3 monomials are dependent at every point, and every eight are independent.

## 10. Addendum (c): does Lean instantiate the bridge with W = S_D?

**Confirmed.**
* `GrowthReduction.not_one_of_fei` is, verbatim, `ChainClosureBridge.prime_m2_not_one k depth a b t two_ne D
  (idealSpace 𝔭 k D) (normalized_mem_idealSpace …) (fun v w hw hdeg => idealSpace_closed …)
  (one_not_mem_idealSpace k 𝔭 D)`.
* `closure_normalForm_mem_idealSpace` proves that the chain normal form of every element of the ordinary
  degree-D closure lies in `idealSpace 𝔭 k D`, i.e. V_D ⊂ S_D. It does this via
  `closure_normalForm_mem` + `quot_closure_minimal`. The chain inputs have normal form 0
  (`normalForm_chain_eq_zero`).

**`idealSpace` is exactly S_D(k) of the text.** By definition,
`idealSpace 𝔭 k D = {w ∈ reducedSpace(chainInputs K 2 k) | totalDegree w ≤ D ∧ expand K 2 k w ∈ 𝔭}`.
* The reduced space is the span of the squarefree chain monomials: `standard_le_one` and `digits_reduced`.
  These monomials correspond to the basis x^b, b ∈ [0, 2^k)², of A_k.
* Their degree is wt(b₁) + wt(b₂) (`digits_degree`). So "totalDegree ≤ D" is "weight ≤ D" (`boxed_expand`,
  `lower_totalDegree_le`).
* `expand` (u_{i,j} ↦ x_i^(2^j)) sends the monomial of b to x^b with the same coefficient (`expand_digits`). On
  reduced polynomials it is a bijection onto polynomials with exponents in [0, 2^k)², with inverse `lower`
  (`expand_lower`). This is the text's lift.
* 𝔭 is an arbitrary prime here. `GrowthMain` instantiates it as ker(aeval ξ), the text's 𝔭.

**Hypotheses discharged, against Props 2.1 and 2.2:**
* `gen`: the normal form of the generator is the generator (`normalForm_normalized`). It is reduced, of degree 2,
  with lift S_3 (`expand_normalized`) ∈ 𝔭. This is the text's "g ∈ S_D(k) (k ≥ 2, g̃ = S_3)".
* `products`: `idealSpace_closed`. For w ∈ S_D with totalDegree w < D and every chain variable v,
  normalForm(X v·w) ∈ S_D, under FEI_k(D−1) with 𝔭 prime and x_i ∉ 𝔭. This is Prop 2.1, "FEI ⇒ closed". The proof
  is the cleared lift identity plus primality. The direction used does not need x_i^N − 1 ∉ 𝔭, matching the
  text.
* `one`: 1 ∉ S_D, because 𝔭 is proper (`one_not_mem_idealSpace`).
* Prop 2.2 is `not_one_of_fei` together with `closure_normalForm_mem_idealSpace`. The converse direction of
  Prop 2.1 is `GrowthCriterion.closed_iff_fei` (with x_i^N − 1 ∉ 𝔭).
* All have standard axioms only (printed).

So the N8 sentences saying that this translation "has not been formalized" are stale (R2-1). The first review's
remark to the same effect (its §6 and §7) predates the Lean files.

---

## 11. What I verified myself, and how, versus what I only read

**Verified myself:**
* All pins (`git rev-parse`, `git hash-object`).
* Lean, under the lake lock:
  - `lake build --no-build` for `GrowthExistence` and `GrowthCriterion` (up to date against the sources, which
    covers all six Growth modules and their dependencies);
  - three scratch files with `#print axioms` (38 declarations, including all 35 Growth entries of `Audit.lean`
    l. 556–590; every one depends only on `[propext, Classical.choice, Quot.sound]`);
  - the elaborated statements of `CurveModel`, `SpreadBound`, `Circuit`, `mono`, `LambdaBound`, `Lam`, `Csum`,
    `FEI`, `Boxed`, `lift_identity`, `growth_theorem`, `growth_closure_vanishes`, `unit_instance`,
    `growth_instance_window`, `unit_zmod`, `idealSpace`, `not_one_of_fei`, `closure_normalForm_mem_idealSpace`,
    `idealSpace_closed`, `normalized_mem_idealSpace`, `one_not_mem_idealSpace`, `expand_normalized`,
    `boxed_expand`, `expand_digits`, `expand_lower`, `lower_totalDegree_le`, `prime_m2_not_one`, `digits_degree`,
    `normalForm_X_mul_monomial`.
* Hand derivations of every step listed in §1.
* The four scripts in `experiments/ic-bounds-growth/review2/` (§§6–9).
* Arithmetic spot checks of N8 §6.2:
  - |M_≤3| and |M_≤4| headers; the D = 5 monomial counts and 55,455;
  - dim S_3 = 4k + 98 against the tables; 3855 + 294 + 302 = 4451 and its match with p₄ of V_6(14) in
    `experiments/ic-bounds-k14/RESULTS.md` l. 21;
  - 286 = 2·142 + 2.
* The instances of Section 6 are nondegenerate: t ≠ 0, R not 2-torsion, 4a³ + 27b² ≠ 0, for (12289, 4021, 9917,
  4408), (40961, 21223, 9887, 4748), (12289, 1, 0, 9362) and (65537, 1, 0, 4271).

**Only read:**
* the internals of the audited background modules (`ChainClosureBridge`, `ChainQuotient`, `UnitInstance`,
  `SelectorExtension`, `PrimeM2General`) beyond their statements;
* the exact data of N8 §§6–7 (closures, FEI verdicts, the D = 5, k = 12 certificate, E3);
* the attributions to Brownawell–Masser and Voloch.

## 12. What I could not verify

* The original statements of Brownawell–Masser (1986) Theorem B and of Voloch (1985). I did not retrieve them. The
  characteristic-p proof written in N8 is self-contained and I checked it, so correctness does not depend on them;
  only the attribution does.
* The exact closure and kernel data of N8 §6.2, §7.2 and E1–E4 (not re-run; no heavy lock taken).
* The F^p tests are searches for counterexamples. They are exhaustive only for m ≤ 3 on the tested curves, so they
  are evidence, not proof, for Thm 3.3 at m ≥ 4.
* That `CurveModel` is satisfiable in Lean. The function field of E and Wronskians are absent from the pinned
  Mathlib. This is the named gap of the Lean development.

## 13. Where I differ from the first review (read after forming the above)

The first review applies to blob `6c98619e`. I agree with its verdicts (T1, C1, T4, Prop 7.1, Prop 7.2: proof; no
fatal error) and with its re-derivations. Its F1–F12 are resolved in blob `0e58d23c`; I checked the F1 repair. The
differences:

1. Its §6/§7 say that nobody has formalized W = S_D and that Prop 7.1 would need "a small Lean file". At commit
   09dd9760 both exist (`GrowthReduction`, `GrowthCriterion`). I therefore recommend **kernel** for T3 and Prop 7.1,
   and the conditional kernel label for T1/C1. This is new information, not an error in that review.
2. It did not see the paper's abstract and introduction wording (R2-4, R2-5) or the sketch (R2-6, R2-7), which were
   written later.
3. It gave E2 no scope caveat; I add R2-9.
4. Minor: its §3 says Thm 3.1 (d) "holds with Hasse or ordinary derivatives alike". That is true only because r ≤ p
   is forced by F^p-independence. With Hasse derivatives of orders ≥ p the criterion is over the constants and the
   statement changes. Not a disagreement about T1.

## 14. Commands run (all under `ulimit -v 1500000`, single thread, except the Lean commands)

```
git rev-parse 09dd9760:docs/proofs/N8-growth.md; git hash-object docs/proofs/N8-growth.md   (and all files above)
cd formal && flock /tmp/lfd-locks/lake.lock sh -c 'ulimit -v 7000000; ~/.elan/bin/lake build --no-build ICFallDegrees.GrowthExistence'
cd formal && flock /tmp/lfd-locks/lake.lock sh -c 'ulimit -v 7000000; ~/.elan/bin/lake build --no-build ICFallDegrees.GrowthCriterion'
cd formal && flock /tmp/lfd-locks/lake.lock sh -c 'ulimit -v 7000000; ~/.elan/bin/lake env lean <scratchpad>/redteam2/Axioms.lean'
cd formal && flock /tmp/lfd-locks/lake.lock sh -c 'ulimit -v 7000000; ~/.elan/bin/lake env lean <scratchpad>/redteam2/Bridge.lean'
cd formal && flock /tmp/lfd-locks/lake.lock sh -c 'ulimit -v 7000000; ~/.elan/bin/lake env lean <scratchpad>/redteam2/Axioms2.lean'
cd experiments/ic-bounds-growth/review2
python3 -I lambda_check.py 60 lambda_check.json            # 4 s, ALL CHECKS PASS
python3 -I lambda_check.py 120 lambda_check_w120.json      # 52 s, ALL CHECKS PASS
python3 -I lift_check.py 6 14 lift_check.json              # 10 s, ALL CHECKS PASS
python3 -I fp_circuits.py 61 5 7 3 12 30 1 m3              # NO VIOLATION FOUND (also: 61 19 0 2; 101 5 7 3 m3; 131 7 11 5 nom3; 101 30 0 2 nom3)
python3 -I divisor_check.py lemma32 lemma32.json           # 3 s, LEMMA 3.2 CHECKS PASS
python3 -I divisor_check.py twisted 11 1 3 2 twisted_11_1_3_2.json   # TWISTED-RELATION CHECKS PASS
python3 -I divisor_check.py twisted 13 1 0 2 twisted_13_1_0_2.json   # TWISTED-RELATION CHECKS PASS
```

`<scratchpad>` is `/tmp/claude-0/-home-user/a6e79f2a-a5ec-5904-b377-343a14909fff/scratchpad`. The scratch Lean files
are not part of the repository.
