# SD1 in Lean: plan for the static solving degree sd^stat = δ + maxwt(supp g^{-1})

Date: 2026-11-15. Source statement: `docs/proofs/SD1-solving-degree.md` Sections 1–4 (Lemma 2.1 (a), (b);
Theorem 3.1; Corollaries 4.1–4.4). Sections 5–7 of SD1 (proper ideals, binary case) are out of scope.

Status: **plan only. Nothing here has been compiled.** No Lean toolchain is installed in this container
(Section 7). Every Lean statement below is a sketch, "as precise as possible without checking it
compiles". Declaration names marked *existing* were verified by reading the file named; names marked
*new* are proposals. SD1 itself is **[pcc:unverified]** (no pcc claim id), so every Lean theorem planned
here would formalize an unreviewed written proof. Assertions about the mathematics that are not in SD1
carry [pcc:unverified].

## 1. Target and route

**Target (conditional on nothing but the unit hypothesis).** For a field `K` with `2 ≠ 0`, `m k : ℕ`,
and a generator `ĝ : P (m*k) K` that is reduced modulo the chains (as the normalized S_3 is,
`ChainModelMap.normalized_reduced`) and non-zero, with `1 ∈ Ideal.span (chainInputs K m k ∪ {ĝ})`:
there is a unique reduced `h` with `ĝ * h - 1 ∈ Ideal.span (chainInputs K m k)`, and

    ∀ D, (1 : P (m*k) K) ∈ StaticSpan (chainInputs K m k ∪ {ĝ}) D ↔ ĝ.totalDegree + h.totalDegree ≤ D

so `sd^stat = δ + h.totalDegree`.

**Route: stay in the reduced polynomial ring; do not build the group algebra.** This departs from the
suggestion to define `A = K[(Z/N)^m]` through `MonoidAlgebra`/`Finsupp`. SD1's quotient model
`A = K[G]`, `wt`, `M_{≤d}` and `ι` are all already present in the repository in polynomial form:

| SD1 object | Lean object (existing) |
|---|---|
| `NF : R → A` | `ChainClosureBridge.normalForm (chainInputs K m k) (chainInputs_nonzero m k)` (linear: `ChainClosureBridge.nfMap`) |
| basis `x^e`, `e ∈ G`, and the reduced lift `ι` | reduced monomials `monomial (ChainTranslation.digits k b) 1`, `b i < 2^k`; `ChainModelMap.digits_injOn`, `ChainQuotient.bitExp_digits`/`ChainModelMap.bitExp_digits`, `ChainQuotient.bitExp_injective` |
| `wt(e) = Σ popcount(e_i)` | `BinaryWeight.wt`; `ChainModelMap.digits_degree : (digits k b).degree = ∑ i, wt (b i)` (needs `∀ i, b i < 2^k`) |
| `M_{≤d}` | the submodule `reducedSpace (chainInputs K m k) ⊓ {p | p.totalDegree ≤ d}` (new `MLe`) |
| `wt(a)` for `a ∈ A` | `(ι a).totalDegree` = `totalDegree` of the reduced representative (new `maxwt` is simply `totalDegree` on `reducedSpace`) |
| `g = NF ĝ` | `normalForm … ĝ` (equals `ĝ` when `ĝ` is reduced) |
| product in `A` | `normalForm (x * y)`; compatibility lemma `ChainClosureBridge.normalForm_mul_normalForm` |
| `g^{-1}` | the unique reduced `h` with `ĝ*h - 1 ∈ Ideal.span chains` (Step 5) |

Consequences, both better than the written proof:

* The weight inequality (1.1) and the "reachable exponent" lemma `wt(e) ≤ d ⇔ ∃ monomial μ, deg μ ≤ d,
  NF μ = x^e` are **not needed** in Lean. Both directions of Lemma 2.1 (b)'s first claim reduce to
  `ChainClosureBridge.normalForm_totalDegree_le` (NF does not raise degree) and
  `ChainClosureBridge.normalForm_of_reduced` (a reduced polynomial is its own NF). The sketch in
  Section 3 Step 3 shows the details. If the exponent form is still wanted for a statement in terms
  of `wt`, it is one corollary (Step 3b) from `digits_degree`.
* The group-algebra inverse (`Finsupp`/`MonoidAlgebra` inverse and its uniqueness) is replaced by
  `reduced_unique_of_reduced_zero` + `reduced_chain_ideal_trivial` (uniqueness) and the extraction
  of `h` from `1 ∈ Ideal.span (… ∪ {ĝ})` (existence, Step 5). A `MonoidAlgebra` model would only be
  needed for Corollary 4.2 (the Fourier/character-sum formula for `c_top`), which is outside this plan
  (Section 5).

Not present anywhere in `formal/` (grep, 2026-11-15, `grep -ril "macaulay\|static" formal/`): there is
no definition of a static or Macaulay span. The two hits are comments only:
`QuadraticNoLow.lean:5` ("refutations in the full iterated closure, not merely a static Macaulay
span") and `GrowthConcreteClosureStep.lean:12` ("one-shot Macaulay matrix"). `FilteredSpan.lean`
(`rowSpace`, `closure_le_rowSpace`) is an RREF row space for lower bounds, not a static span; it is not
reused. `InClosure` (`Closure.lean:28`) is the **iterated** closure: its `mul` constructor multiplies
already derived elements, so it is V_{F,D}, not Mac_D. This is the central modelling point: Lemma 2.1 (a)
is a statement about the one-layer span, and the existing `normalForm_reduction`
(`ChainClosureBridge.lean:56`, an `InClosure G` statement) does **not** give it (Step 2).

## 2. New definitions (file `StaticSpan.lean`, namespace `ICFallDegrees.StaticSpan`)

```lean
import ICFallDegrees.ChainClosureBridge
import ICFallDegrees.ChainModelMap      -- reduced_of_le_one, digits_degree

open MvPolynomial FiniteQuotient ChainQuotient ChainClosureBridge

variable {r : ℕ} {K : Type*} [Field K]

/-- Mac_D: span of the one-layer multiples `μ * f`, `f ∈ F`, `μ` a monomial, `deg (μ f) ≤ D`. -/
def StaticSpan (F : Set (P r K)) (D : ℕ) : Submodule K (P r K) :=
  Submodule.span K {p | ∃ f ∈ F, ∃ e : Fin r →₀ ℕ,
    p = monomial e 1 * f ∧ (monomial e 1 * f).totalDegree ≤ D}

/-- sd^stat for a unit ideal, as a predicate (least D is `Nat.find` when it exists). -/
def StaticOne (F : Set (P r K)) (D : ℕ) : Prop := (1 : P r K) ∈ StaticSpan F D

/-- M_{≤d}: reduced polynomials (squarefree in the chain variables) of total degree ≤ d. -/
def MLe (m k d : ℕ) : Submodule K (P (m*k) K) :=
  reducedSpace (chainInputs K m k) ⊓ {p | p.totalDegree ≤ d}   -- needs a `Submodule` for the degree filter
```

For `MLe` use `Submodule.span K {monomial e 1 | …}` or the structure literal with
`totalDegree_add` / `totalDegree_smul_le` (`Closure.lean:25`, `Mathlib` `MvPolynomial.totalDegree_add`);
the second is easier because `closureSpace` (`Closure.lean:92`) already shows the shape of the
literal. Provide `mem_MLe : p ∈ MLe m k d ↔ p ∈ reducedSpace _ ∧ p.totalDegree ≤ d` and `MLe_mono`.
Difficulty: trivial. The one-layer span is stated with the `monomial e 1` form of SD1 ("μ a monomial");
`mul_mem_staticSpan` (Step 1) shows it equals the `a * f` form with `a.totalDegree + f.totalDegree ≤ D`
used by `InClosure.mul`.

## 3. Step by step

Axiom expectation for every step: `[propext, Classical.choice, Quot.sound]`, the profile of the whole
bridge (`N8-lean-plan.md` §1, U1 plan §5). `normalForm` is defined by `Classical.choose`
(`ChainClosureBridge.lean:50`), so `Classical.choice` is unavoidable and expected. No `native_decide`,
no `sorry`, no new `axiom`.

### Step 1. Basic facts on `StaticSpan` and comparison with the iterated closure (Corollary 4.1, first inequality)

```lean
theorem mul_mem_staticSpan {F : Set (P r K)} {D : ℕ} {f : P r K} (hf : f ∈ F)
    (a : P r K) (h : a.totalDegree + f.totalDegree ≤ D) : a * f ∈ StaticSpan F D
theorem StaticSpan_mono {F G} (h : F ⊆ G) {D E} (hDE : D ≤ E) : StaticSpan F D ≤ StaticSpan G E
theorem StaticSpan_le_ideal (F : Set (P r K)) (D : ℕ) : (StaticSpan F D : Set _) ⊆ Ideal.span F
theorem staticSpan_le_closureSpace (F : Set (P r K)) (D : ℕ) : StaticSpan F D ≤ closureSpace F D
theorem hasLFDBound_of_staticOne {F : Set (P r K)} {D} (h : StaticOne F D) : HasLFDBound F D
theorem staticOne_gt_of_not_lfd {F : Set (P r K)} {d D} (h : StaticOne F D)
    (hd : ∀ e ≤ d, ¬ HasLFDBound F e) : d < D
```

* Reuse: `Closure.lean` `InClosure.input`, `InClosure.mul`, `closureSpace`, `closure_minimal`,
  `HasLFDBound`, `hasLFDBound_of_one` (`Closure.lean:122`), `InClosure.sound`; Mathlib
  `MvPolynomial.totalDegree_mul_of_isDomain` (used in sibling files, e.g.
  `crypto/research/guarded_fall_degree_20261028/formal/FilteredClosure.lean:127`, signature
  `(ha : a ≠ 0) (hb : b ≠ 0) : (a*b).totalDegree = a.totalDegree + b.totalDegree`),
  `MvPolynomial.le_totalDegree`, `MvPolynomial.as_sum`.
* New: all six lemmas. `mul_mem_staticSpan`: write `a = Σ_{e ∈ a.support} monomial e (coeff e a)`
  (`as_sum`), each `monomial e c = c • monomial e 1`, and `(monomial e 1 * f).totalDegree ≤ deg e + f.totalDegree ≤ a.totalDegree + f.totalDegree`
  via `le_totalDegree` and `totalDegree_mul`. `staticSpan_le_closureSpace`: `Submodule.span_le`, then
  `InClosure.mul (monomial e 1) (InClosure.input hf _) _` with the degree budget from the domain
  lemma (case `f = 0` is `InClosure.zero`).
* This gives the **cheap unconditional corollary** (Corollary 4.1, "d_lf ≤ sd^stat"):
  `hasLFDBound_of_staticOne`, and `staticOne_gt_of_not_lfd` makes every existing lower bound in the
  repository a lower bound on sd^stat: e.g. `CarryLower.chain_carry_lower` (`CarryLower.lean:76`,
  `∀ e ≤ d, ¬ HasLFDBound (chainInputs K m k ∪ {g}) e`) and `GrowthMain.growth_lfd_gt`
  (`N8-lean-plan.md` §1) give `d < D` whenever `StaticOne … D`.
* Difficulty: easy (about 80 lines). The equality "sd^CG = d_lf" of Corollary 4.1 is **not** Lean-able
  here: `sd^CG` is not defined in `formal/` and the claim is Caminata–Gorla Theorem 3.1 (cited in SD1,
  quotes "checked by the coordinator, not by lane D"). Label it **cited, not formalized**.

### Step 2. Lemma 2.1 (a): ker NF ∩ R_{≤D} ⊆ one-layer chain span

```lean
theorem reduce_static (m k D : ℕ) (e : Fin (m*k) →₀ ℕ) (he : e.degree ≤ D) :
    ∃ r : P (m*k) K, r ∈ reducedSpace (chainInputs K m k) ∧ r.totalDegree ≤ e.degree ∧
      monomial e 1 - r ∈ StaticSpan (chainInputs K m k) D
theorem reduce_static_poly (p : P (m*k) K) (D : ℕ) (hp : p.totalDegree ≤ D) :
    ∃ r, r ∈ reducedSpace (chainInputs K m k) ∧ r.totalDegree ≤ p.totalDegree ∧
      p - r ∈ StaticSpan (chainInputs K m k) D
theorem kernel_mem_static (two_ne : (2 : K) ≠ 0) (p : P (m*k) K) (D : ℕ)
    (hp : p.totalDegree ≤ D) (hnf : normalForm (chainInputs K m k) (chainInputs_nonzero m k) p = 0) :
    p ∈ StaticSpan (chainInputs K m k) D
theorem static_chain_le_kernel (two_ne : (2 : K) ≠ 0) {p} (hp : p ∈ StaticSpan (chainInputs K m k) D) :
    normalForm (chainInputs K m k) (chainInputs_nonzero m k) p = 0
```

* **What is new and why.** The written proof divides by a Gröbner basis and checks that every division [pcc:unverified]
  step is a one-layer multiple. The existing division (`FiniteQuotient.exists_reduced_remainder`,
  `FiniteQuotient.lean:56`, surfaced as `normalForm_reduction : InClosure G p.totalDegree (p - NF p)`)
  yields an **iterated** closure statement, which is weaker than Mac_D. Do not try to extract
  the static statement from it. Instead prove (a) by a direct strong induction on `e.degree`, which is
  shorter than reusing the order machinery [pcc:unverified: this is an alternative to SD1's
  Buchberger-criterion proof, not in SD1; the mathematics is the same one-step rewrite]:
  if `e` is squarefree (`∀ v, e v ≤ 1`) take `r = monomial e 1` (reduced by
  `ChainModelMap.reduced_of_le_one`, `ChainModelMap.lean:350`). Otherwise some `e v ≥ 2`;
  write `e = 2·single v + e'`. With `chainPolynomial K m k v = X v ^ 2 - chainTail K m k v`
  (`ChainQuotient.lean:79`), `monomial e 1 - monomial e'' 1 = monomial e' 1 * chainPolynomial v`
  where `e'' = e' + single (next v) 1` or `e'` (tail `1`, end of chain). Degree of the multiple is
  `e.degree` (domain lemma), and `e''.degree < e.degree`, so the induction hypothesis applies with the
  same `D`. `StaticSpan_mono` is not even needed because `D` is fixed.
  Add `Submodule.add_mem` to chain the two differences.
* Reuse: `ChainQuotient.chainPolynomial`, `chainTail`, `chainTail_totalDegree_le` (`:93`),
  `ChainModelMap.reduced_of_le_one` (`:350`), `ChainQuotient.standard_le_one` (`:131`, the converse:
  reduced ⇒ squarefree), `ChainClosureBridge.normalForm_eq_of` (`:71`, uniqueness, to identify `r`
  with `NF p`), `ChainQuotient.reduced_chain_ideal_trivial` (`:244`, the `zero` hypothesis, needs
  `2 ≠ 0`), `Ideal.span` soundness for `StaticSpan_le_ideal`.
* `kernel_mem_static`: `p - r ∈ StaticSpan ⊆ Ideal.span chains`, so `normalForm_eq_of` gives `NF p = r`;
  with `NF p = 0` get `r = 0` and `p = p - r`. This also covers `D < 2` with no extra case split (the
  written proof treats it separately); `D` enters only through `hp`.
* `static_chain_le_kernel` (the "⊇" of (a)): every generator `monomial e 1 * c_v` lies in `Ideal.span`
  chains, so `ChainClosureBridge.normalForm_of_mem` (`:82`) gives `0`; `nfMap` is linear
  (`Submodule.span_induction`).
* Difficulty: **medium** (about 150 lines). The only fiddly part is the dependent index `next v`
  (`Fin (m*k)` with `row`/`depth`, `ChainQuotient.var`, `depth_val`, `var_val`); the same index
  gymnastics occur in `ChainTranslation.digits_update_of_bit_zero` (`:171`) and
  `ChainClosureBridge.inputs_castSucc_mem_chainInputs` (`:424`), and can be copied.

### Step 3. Lemma 2.1 (b): NF(Mac_D) = g · M_{≤D−δ}

```lean
theorem nf_mem_MLe (p : P (m*k) K) :
    normalForm (chainInputs K m k) (chainInputs_nonzero m k) p ∈ MLe m k p.totalDegree
theorem MLe_eq_nf_image (two_ne : (2 : K) ≠ 0) (d : ℕ) :
    (MLe m k d : Set _) = normalForm … '' {p | p.totalDegree ≤ d}
theorem nf_mul_eq (two_ne : (2 : K) ≠ 0) (a b : P (m*k) K) :
    normalForm … (a * b) = normalForm … (normalForm … a * normalForm … b)
theorem nf_staticSpan (two_ne : (2 : K) ≠ 0) {ĝ : P (m*k) K} (D : ℕ) :
    Submodule.map (nfMap …) (StaticSpan (chainInputs K m k ∪ {ĝ}) D) =
      Submodule.map ((nfMap …).comp (LinearMap.mulLeft K ĝ)) (MLe m k (D - ĝ.totalDegree))
```

with the convention `MLe m k d = ⊥` for the truncated subtraction case `D < δ` (use
`if D < ĝ.totalDegree then ⊥ else MLe m k (D - ĝ.totalDegree)`; truncated `Nat.sub` would give
`MLe … 0`, which is wrong. Alternatively state it with `∀ p, p ∈ StaticSpan … D → ∃ h ∈ MLe m k d,
δ + d ≤ D ∧ NF p = NF (ĝ*h)`.) Several generators: replace `{ĝ}` by a finite family and sum; only the
`m = 3` (one generator `chainS4`) and `m = 2` (one generator) cases are used, so the several-generator
form is **not planned**.

* Reuse: `normalForm_totalDegree_le` (`:63`), `normalForm_reduced` (`:53`), `normalForm_of_reduced`
  (`:78`), `normalForm_mul_normalForm` (`:116`, gives `NF(a*p) = NF(a*NF p)`; `nf_mul_eq` is two
  applications plus `mul_comm`), `normalForm_of_mem` (`:82`), `nfMap`/`nfMap_apply` (`:131`).
* New content. First claim of (b) in polynomial form: `nf_mem_MLe` is `normalForm_reduced` +
  `normalForm_totalDegree_le`; `MLe_eq_nf_image` "⊆" uses `normalForm_of_reduced` (a reduced `p` with
  degree ≤ d is its own NF), "⊇" uses `nf_mem_MLe`. So `span NF(R_{≤d}) = M_{≤d}` follows without
  the weight inequality (1.1) and without the monomial-by-monomial argument of SD1.
  Second claim (`nf_staticSpan`): `Submodule.map_span` turns it into the images of the generators.
  For `f ∈ chains`, `NF(μ f) = 0` (`normalForm_of_mem`, since `μ f ∈ Ideal.span chains`). For `f = ĝ`:
  `NF(μ ĝ) = NF(NF μ * ĝ)` by `nf_mul_eq`/`normalForm_mul_normalForm`, with `NF μ ∈ MLe (deg μ)` and
  `deg μ = deg(μ ĝ) − δ ≤ D − δ` (domain lemma, `ĝ ≠ 0`). The reverse inclusion takes `h ∈ MLe` to the
  row `h * ĝ` via `mul_mem_staticSpan` (Step 1), using `h.totalDegree + δ ≤ D`.
* Step 3b (exponent form, optional): for `b : Fin m → ℕ` with `∀ i, b i < 2^k`,
  `(monomial (digits k b) 1) ∈ MLe m k d ↔ ∑ i, wt (b i) ≤ d`, from `ChainModelMap.digits_degree` (`:87`),
  `ChainTranslation.digits_reduced` (`:147`). Only needed to restate the main theorem in `wt` language.
* Difficulty: **medium-easy** (about 120 lines). Remark (ii) of SD1 (δ is the *ordinary* degree of `ĝ`,
  not `wt g`) is automatic: the statement quantifies over `ĝ.totalDegree`, and the normalized generators
  are reduced anyway (`ChainModelMap.normalized_reduced`, `ChainModelMap.lean:396`). If one wants raw S_3
  (`PrimeM2General.rawS3`, δ = 4, `wt g = 2`), apply the theorem with `ĝ = rawS3` and `g = NF rawS3 =
  normalized` (`PrimeM2General.raw_normalization_identity`, `PrimeM2General.lean:41`); nothing extra.

### Step 4. Reduced inverse: existence and uniqueness (group-algebra inverse replacement)

```lean
theorem exists_reduced_inverse (ĝ : P (m*k) K)
    (unit : (1 : P (m*k) K) ∈ Ideal.span (chainInputs K m k ∪ {ĝ})) :
    ∃ h, h ∈ reducedSpace (chainInputs K m k) ∧ ĝ * h - 1 ∈ Ideal.span (chainInputs K m k)
theorem reduced_inverse_unique (two_ne : (2 : K) ≠ 0) {ĝ h h'}
    (hh : h ∈ reducedSpace _ ∧ ĝ * h - 1 ∈ Ideal.span (chainInputs K m k))
    (hh' : h' ∈ reducedSpace _ ∧ ĝ * h' - 1 ∈ Ideal.span (chainInputs K m k)) : h = h'
def reducedInverse (...) : P (m*k) K := Classical.choose (exists_reduced_inverse ĝ unit)
```

* Existence (new): `Ideal.span_union`/`Ideal.mem_span_insert`-style decomposition gives
  `1 = q * ĝ + c` with `c ∈ Ideal.span chains`; take `h = normalForm q`; then `ĝ * h - 1 ≡ ĝ * q - 1 ≡ 0`
  using `sub_normalForm_mem` (`:60`). Mathlib: `Ideal.mem_span_insert`, `Ideal.span_union`,
  `Ideal.mem_span_singleton'`.
* Uniqueness: `ĝ (h − h') ∈ Ideal.span chains` does **not** by itself give `h = h'`; one needs the
  inverse property. Do it as in SD1 ("inverses are unique") in the quotient: `h − h' ≡ h·(ĝ h' − 1) −
  h'(ĝ h − 1) ≡ 0` (polynomial identity `h − h' = h*(ĝ*h' − 1) − h'*(ĝ*h − 1)` by `ring`), so
  `h − h' ∈ Ideal.span chains` and is reduced (`(reducedSpace _).sub_mem`), hence `0` by
  `reduced_chain_ideal_trivial` (the existing `ReducedIdealTrivial`, `FiniteQuotient.lean:47`) —
  this is `FiniteQuotient.reduced_unique_of_reduced_zero` (`:99`) or `ChainClosureBridge.normalForm_eq_of`.
* Difficulty: **easy** (about 60 lines).
* The group algebra route (rejected, for reference): a `MonoidAlgebra K (Multiplicative (Fin m → ZMod (2^k)))`
  model with `A^×` and a `Units`-based `g⁻¹` would need an isomorphism with the reduced ring and a
  transport of `totalDegree` to the support weight; Step 4 avoids all of it.

### Step 5. Theorem 3.1 (the main theorem)

```lean
theorem staticOne_iff (two_ne : (2 : K) ≠ 0) {ĝ h : P (m*k) K} (hg : ĝ ≠ 0)
    (hinv : h ∈ reducedSpace (chainInputs K m k) ∧ ĝ * h - 1 ∈ Ideal.span (chainInputs K m k)) (D : ℕ) :
    StaticOne (chainInputs K m k ∪ {ĝ}) D ↔ ĝ.totalDegree + h.totalDegree ≤ D
theorem staticSolvingDegree (two_ne : (2 : K) ≠ 0) {ĝ} (hg : ĝ ≠ 0)
    (unit : (1 : P (m*k) K) ∈ Ideal.span (chainInputs K m k ∪ {ĝ})) :
    ∃ h, reducedInverse-spec h ∧ IsLeast {D | StaticOne (chainInputs K m k ∪ {ĝ}) D}
      (ĝ.totalDegree + h.totalDegree)
```

Proof (mirrors SD1 Theorem 3.1):

* "⇒": `1 ∈ StaticSpan … D`; apply `nfMap` and `nf_staticSpan` (Step 3): `1 = NF 1 = NF(ĝ h₁)` with
  `h₁ ∈ MLe (D − δ)` and `δ ≤ D`. Then `ĝ h₁ − 1 ∈ Ideal.span chains` and `h₁` reduced, so
  `h₁ = h` by `reduced_inverse_unique` (Step 4); thus `h.totalDegree ≤ D − δ`.
  Needs `normalForm_one` (`ChainClosureBridge.lean:126`, with `ChainClosureBridge.chain_one_reduced`,
  `:368`) for `NF 1 = 1`.
* "⇐": put `P := h * ĝ`. Then `P ∈ StaticSpan … D` by `mul_mem_staticSpan` (`h.totalDegree + δ ≤ D`).
  `P − 1 ∈ Ideal.span chains` by hypothesis `hinv`, `(P − 1).totalDegree ≤ D`
  (`totalDegree_sub_le`, `Closure.lean:25`; `totalDegree_mul`), so `kernel_mem_static` (Step 2) with
  `NF (P − 1) = 0` (`normalForm_of_mem`) gives `P − 1 ∈ StaticSpan chains D ⊆ StaticSpan (chains ∪ {ĝ}) D`
  (`StaticSpan_mono`). Hence `1 = P − (P − 1)`.
* `staticSolvingDegree`: existence of `h` is Step 4; `IsLeast` from `staticOne_iff` and monotonicity
  of `StaticSpan` in `D`.
* The SD1 equivalences (i)⇔(ii)⇔(iii) of Theorem 3.1: (i)⇔(ii) is the `h` of Step 4 (`1 ∈ Ideal.span F`
  ⇔ `∃ h` reduced inverse; "⇐" is `ĝ h − (ĝ h − 1) = 1`). (iii) (no zero on the chain variety) is the
  weak Nullstellensatz over `K̄` and is **not planned**; the repository already has the grid form:
  `UnitInstance.one_mem_of_no_grid_zero` (`UnitInstance.lean:137`, combinatorial Nullstellensatz) and
  `GridQuotient`. Instances of the unit hypothesis are available unconditionally:
  `UnitInstance.unit_of_primitiveRoot` (`:265`) and `UnitInstance.unit_zmod` (`:352`), both with
  `b = 0`, `t` a non-square, `Y² = t³ + a t`, `Y ≠ 0`, `2^(k+1) ∣ p − 1`, concluding
  `(1 : P k K) ∈ Ideal.span (Set.range (inputs k depth a 0 t))`; `UnitInstance.unit_p65537` (`:372`) and
  `unit_p163841` (`:390`) are fixed-parameter instances (`2 ≤ k ≤ 15` and `2 ≤ k ≤ 14`).
  To use them, rewrite the range as `chainInputs K 2 k ∪ {normalized k depth a b t}` (one inclusion is
  `ChainClosureBridge.inputs_castSucc_mem_chainInputs`, `:424`, the other `ChainQuotient.chainInputs_two_subset`,
  `:361`; **new** lemma `range_inputs_eq : Set.range (inputs k depth a b t) = chainInputs K 2 k ∪ {normalized k depth a b t}`,
  `Fin.lastCases` split, easy).
* Difficulty: **medium** (about 200 lines including the Step 5 glue). This is the step that is also
  cheap to cross-check against the exact data of SD1 §3 (68/68 unit instances), but only after a
  toolchain exists (Section 7).

### Step 6. Corollary 4.1 (upper bound and the top-coefficient characterization)

```lean
theorem staticSolvingDegree_le (two_ne : (2 : K) ≠ 0) … :
    ĝ.totalDegree + h.totalDegree ≤ ĝ.totalDegree + m * k
theorem totalDegree_reduced_eq_iff (h : P (m*k) K) (hred : h ∈ reducedSpace (chainInputs K m k)) :
    h.totalDegree = m * k ↔ coeff (∑ v, Finsupp.single v 1) h ≠ 0
```

* Reuse: `ChainQuotient.reduced_totalDegree_le` (`ChainQuotient.lean:333`, reduced ⇒ `totalDegree ≤ m*k`),
  `standard_degree_bound` (`:137`), `standard_le_one` (`:131`).
* New: `totalDegree_reduced_eq_iff`. "⇐": the all-ones exponent `e* := ∑ v, single v 1` has degree `m*k`
  (`Finsupp.degree`, `Finsupp.sum_fintype`) and lies in the support. "⇒": some `e ∈ support` has
  `e.degree = m*k` with `∀ v, e v ≤ 1`; `Finset.sum_eq_card_nsmul`-style: sum of `m*k` numbers `≤ 1`
  equal to `m*k` forces all equal to 1 (`Finset.sum_lt_sum` or `Finset.card_eq_sum_ones`), so `e = e*`.
* Meaning of `c_top` in `wt` language: `e*` is `digits k (fun _ => 2^k − 1)`
  (`ChainTranslation.digits`, `bit` lemmas `bit_le_one`); `wt (2^k − 1) = k`
  (`BinaryWeight.wt_two_pow_add`, `:95`, or induction). Optional, only to match SD1's `x^{(N−1,…,N−1)}`.
* Difficulty: **easy-medium** (about 80 lines).
* The "iff" form of SD1 ("sd^stat = δ + mk ⇔ c_top ≠ 0") then follows from Steps 5 and 6 immediately.

### Step 7. Corollary 4.4 (generic gap) and the m = 3 version

```lean
theorem static_gap_m2 {K} [Field K] (two_ne : (2:K) ≠ 0) (k) (depth : 2 ≤ k) {ζ} (root : IsPrimitiveRoot ζ (2^k))
    (a b t : K) (nonsingular : 4*a^3 + 27*b^2 ≠ 0)
    (unit : (1 : FiniteQuotient.P (2*k) K) ∈ Ideal.span (Set.range (PrimeM2General.inputs k depth a b t)))
    (hδ : (PrimeM2General.normalized k depth a b t).totalDegree = 2) {h} (hinv : …)
    (ctop : coeff (∑ v, Finsupp.single v 1) h ≠ 0) :
    ∃ D, IsLeast {D | StaticOne (Set.range (PrimeM2General.inputs k depth a b t)) D} (2*k + 2)
       ∧ HasLFDBound (Set.range (inputs …)) (k + 2) ∧ k ≤ (2*k + 2) - (k + 2)
```

* Reuse (all existing): `SelectorBound.lfd_le_depth_add_two` (the U1 theorem, `HasLFDBound … (k+2)`, per [pcc:unverified]
  U1-lean-plan.md §5; for the `2 ≠ 0` form `SelectorExtension.lfd_le_depth_add_two_of_two_ne`, same
  plan §5 table), `ChainModelMap.normalized_reduced`, Steps 5–6. The gap
  `sd^stat − d_lf ≥ k` is then `omega`: `2k + 2 − (k + 2) = k`. The hypotheses are exactly SD1's: "under U1".
  `sd^stat ≥ d_lf` is `hasLFDBound_of_staticOne` only as an upper-bound statement; the Lean statement
  above records the two facts side by side (`HasLFDBound … (k+2)` and `StaticOne … (2k+2)` least), which
  is what "gap ≥ k" means without defining `d_lf` as a minimum (compare `IsLastFallDegree`, `Closure.lean:119`).
  A statement with the exact `d_lf` would need `IsLastFallDegree`, not available in general.
* New: `hδ`. `PrimeM2General.normalized = top (active) (−2(t²+a)) (−2t) + tail` (`PrimeM2General.lean:32`);
  `totalDegree = 2` holds when the leading quadratic part is non-zero; the read of `top`
  (`GeneralQuadraticTop`) is needed to discharge `hδ` rather than assume it [pcc:unverified: that it is 2
  for every nonsingular curve and every `t`; SD1 states δ = 2 for "the normalized S_3" without a
  case analysis; the quadratic coefficient `(x₁ − x₂)²` of S_3 suggests degree 2, but the normalized
  generator's top part has coefficient `−2(t² + a)`, which can vanish when `t² = −a`]. Keep `hδ` as a
  hypothesis unless that case is analysed. This is a real gap in the written source: SD1 Corollary 4.1
  with `δ = 2` silently assumes it.
* `c_top ≠ 0` is a hypothesis, as in SD1 (open problem 1: no proof for all nonsingular curves). No
  existing Lean result supplies it; a concrete instance would need a verified computation of `h`
  (for `k = 2, p = 17`, say, a `decide`/`norm_num` certificate on `h` with `hinv` checked by
  `Ideal.mem_span`-witnesses), which belongs to a separate certificate step in the style of
  `ExactK6*`/`N5Cert*`, not to this plan.
* `m = 3`: `M3CarryLower.chainS4_normal_form` (`:447`, reduced + congruent to the S_4 lift) and
  `M3CarryLower.chainS4_totalDegree` (`:457`, `totalDegree = 6` for `5t³ + 3at + 6b ≠ 0`, `k = n + 3`,
  `2 ≠ 0`) give SD1's δ₃ = 6 on the generic locus, so `staticOne_iff` instantiates to `6 + 3k`.
  Unit-ness for `m = 3` is not provided by `UnitInstance` (two summands); that stays a hypothesis.
  The D ≤ 3 + … lower bounds via `CarryLower` still apply through Step 1.
* Difficulty: **easy given Steps 1–6** (about 60 lines), except `hδ` (medium, depends on `top`).

### Step 8 (optional, not in SD1's claims). Depth monotonicity

`DepthRestriction.restrict_closure` (`:106`) has a static analogue: `restrict` maps `StaticSpan` at depth `k+1`
into `StaticSpan` at depth `k` (same row-by-row argument, `restrict_chainInput`, `:92`). SD1 does not state
it and it is not needed. Not planned. [pcc:unverified: whether sd^stat is monotone in the depth is not
asserted by SD1]

## 4. Dependency graph and order of work

    Step 1 (StaticSpan, comparison)   Step 4 (reduced inverse)
        |                                 |
    Step 2 (Lemma a: kernel)              |
        |                                 |
    Step 3 (Lemma b: NF image) -----------+
        \                                /
         +--------> Step 5 (Theorem) <--+
                       |
                  Step 6 (top coefficient, upper bound)
                       |
                  Step 7 (gap, m = 3)

Suggested modules: `StaticSpan.lean` (Steps 1, 2, 3), `StaticInverse.lean` (Step 4),
`StaticSolvingDegree.lean` (Steps 5, 6, 7). Estimated total: 700–800 lines, comparable to
`SelectorBound.lean` (469 lines) plus `SelectorFibre.lean` (414). Hardest steps: 2 (index handling for the
chain rewrite) and 5 (glue); everything else is short.

## 5. Not covered (and why)

* **Corollary 4.2 (character sum `c_top = N^{-m} Σ_α (α₁⋯α_m)/g(α)`).** Needs the discrete Fourier
  transform on `(Z/N)^m`, `H = μ_N ⊂ K` and `p ∤ N`: a `MonoidAlgebra`/`AddMonoidAlgebra` model and
  `IsPrimitiveRoot` orthogonality (`IsPrimitiveRoot.geom_sum_eq_zero`). This is where the group-algebra
  inverse of the original request would be needed. Existing relatives: `GridQuotient`, `SelectorGrid`
  (`eval_expand`), `SelectorIdempotent`/`SelectorFibre.E`. Medium-hard; separate plan.
* **Corollary 4.3 (|E| ≤ N² − 1, root counting for `P(t)`).** Polynomial in `t` over `ZMod p` with
  orbit grouping; Mathlib `Polynomial.card_roots'` suffices; needs 4.2 first.
* **Corollary 4.4's measured law** `⌊(k+14)/4⌋`: conjecture in SD1, nothing to formalize.
* Weak Nullstellensatz clause (iii) of Theorem 3.1 over `K̄`; `sd^CG = d_lf` (Caminata–Gorla Theorem 3.1,
  cited); Sections 5 (proper ideals, static Gröbner degree), 6 (binary, Boolean ring), 7 (cost remarks).
* The `exact` evidence of SD1 (68/68 unit instances, brute force 18,828/18,828): computational, outside
  Lean. A Lean theorem would add a **kernel** label to the proof, not to those counts.

## 6. Labels after formalization

If Steps 1–6 land: Lemma 2.1 (a), (b) (single generator), Theorem 3.1 (the formula and (i)⇔(ii)), Corollary 4.1
(bounds, top-coefficient iff, `d_lf ≤ sd^stat` in the `HasLFDBound` form) become **kernel** for the stated
hypotheses. They remain **[pcc:unverified]** until a pcc claim is issued: a Lean proof does not by itself
issue one. `sd^CG = d_lf` stays "cited". Corollary 4.4 becomes **kernel, conditional on `hδ`, `c_top ≠ 0`
and the unit hypothesis**, with U1 already kernel (`U1-lean-plan.md` §5).

## 7. Disk and toolchain requirement, and the blocker in this container

* Toolchain: `formal/lean-toolchain` = `leanprover/lean4:v4.19.0`. [pcc:unverified]
* Mathlib: `formal/lakefile.toml` requires `mathlib` at `rev = "v4.19.0"`; `formal/lake-manifest.json` pins
  commit `c44e0c8ee63ca166450922a373c7409c5d26b00b` (`inputRev` `v4.19.0`), together with its
  dependencies (`plausible` and others, see the manifest).
* Build discipline of this repository (`N8-lean-plan.md` §6): each module is built alone through the
  lock, `flock /tmp/lfd-locks/lake.lock sh -c 'ulimit -v 7000000; ~/.elan/bin/lake build ICFallDegrees.<Module>'`,
  i.e. a 7 GB address-space limit; `Mathlib.Tactic` imports can exceed it (`GrowthCounting` avoids them).
  The new modules should import only `Mathlib.Algebra.MvPolynomial.*` through the existing bridge files, as
  `ChainClosureBridge` already does.
* Disk: the Mathlib checkout, the downloaded `.olean` cache (`lake exe cache get`) and the elan toolchain
  together need on the order of 8–12 GB [pcc:unverified: figure from general experience with Mathlib of
  this vintage, not measured in this repository]. Plus the repository's own `.lake/build` for
  `ICFallDegrees` (the modules imported by the new files only).
* **Exact blocker here (verified 2026-11-15):** `df -h /` reports `/dev/vda 252G 33G 6.2G 85% /`, i.e.
  **6.2 GB free**; `~/.elan` does not exist, `which lake lean elan` prints nothing, and `formal/` contains
  no `.lake` directory. So (1) no Lean/Lake binary is installed, and (2) the free space is below the
  expected size of Mathlib plus its cache. No build, `#print axioms` or `lake env lean` check was or can
  be run here. Everything above is unchecked plan.
* To unblock: a runner with ≥ 15 GB free (margin), `elan` installed, `lake exe cache get` run once in `formal/`,
  then the per-module build above and the `#print axioms` lines below. Per the repository's default route
  this is an execution job (`isolab-eks` in the `crypto` repository's instructions); it is not authorized
  by this document.

## 8. Lines for integration (to be added only after a successful build)

`formal/ICFallDegrees.lean`:

```lean
import ICFallDegrees.StaticSpan
import ICFallDegrees.StaticInverse
import ICFallDegrees.StaticSolvingDegree
```

`formal/ICFallDegrees/Audit.lean` (expect `[propext, Classical.choice, Quot.sound]` for each):

```lean
#print axioms ICFallDegrees.StaticSpan.staticSpan_le_closureSpace
#print axioms ICFallDegrees.StaticSpan.hasLFDBound_of_staticOne
#print axioms ICFallDegrees.StaticSpan.staticOne_gt_of_not_lfd
#print axioms ICFallDegrees.StaticSpan.reduce_static_poly
#print axioms ICFallDegrees.StaticSpan.kernel_mem_static
#print axioms ICFallDegrees.StaticSpan.static_chain_le_kernel
#print axioms ICFallDegrees.StaticSpan.nf_staticSpan
#print axioms ICFallDegrees.StaticSpan.exists_reduced_inverse
#print axioms ICFallDegrees.StaticSpan.reduced_inverse_unique
#print axioms ICFallDegrees.StaticSpan.staticOne_iff
#print axioms ICFallDegrees.StaticSpan.staticSolvingDegree
#print axioms ICFallDegrees.StaticSpan.totalDegree_reduced_eq_iff
#print axioms ICFallDegrees.StaticSpan.static_gap_m2
```

Audit expectations: no declaration mentions `sorryAx`, `Lean.ofReduceBool` or a user axiom. The `Classical.choice`
dependence comes from `normalForm` and `Classical.choose` in `reducedInverse`; if `reducedInverse` is
avoided by stating results with `∃ h`/`hinv`, the profile is unchanged. A failed audit on any of the
Step 2 or Step 5 declarations would indicate an `omega`/`decide` fallback to `native_decide` and must be fixed,
not recorded.

## 9. Sources checked

Files read in full or in the cited ranges, 2026-11-15: `docs/proofs/SD1-solving-degree.md` (all),
`docs/proofs/N8-lean-plan.md` (all), `docs/proofs/U1-lean-plan.md` (Sections 1–5),
`formal/ICFallDegrees/Closure.lean` (28–143), `ChainTranslation.lean` (declaration list, 295–313),
`ChainClosureBridge.lean` (44–160, 292–418, 455–489), `ChainQuotient.lean` (72–250, 328–360),
`ChainModelMap.lean` (declaration list, 86–112, 215–232, 346–402), `DepthRestriction.lean` (1–40, declaration list),
`UnitInstance.lean` (1–30, 262–275, 350–360), `CarryLower.lean` (36–87), `FiniteQuotient.lean` (declaration list,
30–50), `FilteredSpan.lean` (1–75), `PrimeM2General.lean` (13–42), `M3CarryLower.lean` (445–462),
`formal/lakefile.toml`, `formal/lake-manifest.json`, `formal/lean-toolchain`. research-kb was not queried
for this plan (a planning note for a formalization of an in-repository proof; no novelty claim is made).
