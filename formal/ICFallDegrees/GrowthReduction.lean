import ICFallDegrees.ChainModelMap
import Mathlib.RingTheory.MvPolynomial.Basic
import Mathlib.RingTheory.Ideal.Prime

noncomputable section

/-!
# Growth theorem, reduction to a splitting property of an ideal
(Section 2 of `docs/proofs/N8-growth.md`: Propositions 2.1 (the direction used) and 2.2)

`Q K = K[x₁, x₂]` is the ring of the two summand coordinates and `𝔭` a prime ideal of it with
`x_i ∉ 𝔭`.  The chain variables expand by `ChainQuotient.expand K 2 k : u_{i,j} ↦ x_i^(2^j)`;
on a reduced (squarefree) chain polynomial `w` this is the **lift** `f̃` of the written proof:
the polynomial with the same coefficients on the exponents in `[0, 2^k)²`.

* `topPart k i j f` is the top-run part `f̃_hi = Σ_{a_i ≥ 2^k - 2^j} c_a x^a`.
* `FEI 𝔭 k w` is the splitting property `FEI_k(w)` of the written proof.
* `idealSpace 𝔭 k D` is the no-wrap ideal space `S_D(k)`, as a subspace of the reduced chain
  polynomials: reduced, of degree (= binary weight) `≤ D`, with lift in `𝔭`.
* `lift_identity` is the identity (1.1), multiplied by `x_i^(2^k - 2^j)` to stay in the
  polynomial ring.
* `idealSpace_closed`: under `FEI 𝔭 k (D-1)`, `S_D(k)` is closed under the degree-`D` rule
  (Proposition 2.1, the direction needed).
* `not_one_of_fei` (Proposition 2.2): under `FEI 𝔭 k (D-1)` and `S_3 ∈ 𝔭`, the constant `1` is
  not in the ordinary degree-`D` closure of the audited system
  `PrimeM2General.inputs k depth a b t`; `closure_normalForm_mem_idealSpace` is the inclusion
  `V_D(k) ⊂ S_D(k)`.

The converse direction of Proposition 2.1 and Proposition 7.1 are in `GrowthCriterion`.
-/

namespace ICFallDegrees.GrowthReduction

open MvPolynomial FiniteQuotient ChainQuotient BinaryWeight ChainTranslation ChainModelMap
open PrimeM2General (inputs normalized active rawS3 raw_normalization_identity)

variable {K : Type*} [Field K]

/-- The polynomial ring of the two summand coordinates. -/
abbrev Q (K : Type*) [Field K] := MvPolynomial (Fin 2) K

/-! ## The top-run part and the splitting property -/

/-- The top-run part of a polynomial in coordinate `i` at level `j`, depth `k`: the monomials
whose `i`-th exponent is at least `2^k - 2^j` (`f̃_hi` in the written proof). -/
def topPart (k : ℕ) (i : Fin 2) (j : ℕ) : Q K →ₗ[K] Q K :=
  (basisMonomials (Fin 2) K).constr K
    (fun d => if 2 ^ k - 2 ^ j ≤ d i then monomial d (1 : K) else 0)

theorem topPart_monomial (k : ℕ) (i : Fin 2) (j : ℕ) (d : Fin 2 →₀ ℕ) (c : K) :
    topPart k i j (monomial d c) = if 2 ^ k - 2 ^ j ≤ d i then monomial d c else 0 := by
  have base : topPart (K := K) k i j (monomial d 1) =
      if 2 ^ k - 2 ^ j ≤ d i then monomial d (1 : K) else 0 := by
    have h := (basisMonomials (Fin 2) K).constr_basis K
      (fun d => if 2 ^ k - 2 ^ j ≤ d i then monomial d (1 : K) else 0) d
    rw [coe_basisMonomials] at h
    exact h
  have hc : monomial d c = c • monomial d (1 : K) := by
    rw [smul_monomial, smul_eq_mul, mul_one]
  rw [hc, map_smul, base]
  split_ifs
  · rfl
  · rw [smul_zero]

/-- The coefficients of the top-run part. -/
theorem coeff_topPart (k : ℕ) (i : Fin 2) (j : ℕ) (f : Q K) (d : Fin 2 →₀ ℕ) :
    coeff d (topPart k i j f) = if 2 ^ k - 2 ^ j ≤ d i then coeff d f else 0 := by
  classical
  conv_lhs => rw [f.as_sum]
  rw [map_sum, coeff_sum]
  simp only [topPart_monomial]
  by_cases hd : d ∈ f.support
  · rw [Finset.sum_eq_single d]
    · split_ifs
      · rw [coeff_monomial, if_pos rfl]
      · rw [coeff_zero]
    · intro e _ hne
      split_ifs
      · rw [coeff_monomial, if_neg hne]
      · rw [coeff_zero]
    · intro h
      exact absurd hd h
  · rw [Finset.sum_eq_zero]
    · rw [not_mem_support_iff.mp hd]
      split_ifs <;> rfl
    · intro e he
      have hne : e ≠ d := fun h => hd (h ▸ he)
      split_ifs
      · rw [coeff_monomial, if_neg hne]
      · rw [coeff_zero]

/-- Exponents in `[0, 2^k)²` and binary weight at most `w`. -/
def Boxed (k w : ℕ) (f : Q K) : Prop :=
  ∀ d ∈ f.support, d 0 < 2 ^ k ∧ d 1 < 2 ^ k ∧ wt (d 0) + wt (d 1) ≤ w

/-- **The splitting property `FEI_k(w)`** (Section 2 of the written proof): for every element of
`𝔭` with exponents in `[0, 2^k)²` and weight `≤ w`, every coordinate `i` and every level
`j < k`, the top-run part lies in `𝔭`. -/
def FEI (𝔭 : Ideal (Q K)) (k w : ℕ) : Prop :=
  ∀ f ∈ 𝔭, Boxed k w f → ∀ (i : Fin 2) (j : ℕ), j < k → topPart k i j f ∈ 𝔭

/-! ## Lifts of reduced chain polynomials -/

/-- A function on `Fin 2` as a finitely supported exponent. -/
def fs (b : Fin 2 → ℕ) : Fin 2 →₀ ℕ := Finsupp.equivFunOnFinite.symm b

@[simp] theorem fs_apply (b : Fin 2 → ℕ) (i : Fin 2) : fs b i = b i := rfl

theorem bitExp_digits_fs {k : ℕ} (b : Fin 2 → ℕ) (hb : ∀ i, b i < 2 ^ k) :
    bitExp (digits k b) = fs b := by
  ext i
  rw [fs_apply]
  exact congrFun (bitExp_digits b hb) i

theorem digits_bitExp {k : ℕ} (e : Fin (2 * k) →₀ ℕ) (he : ∀ v, e v ≤ 1) :
    digits k (fun i => bitExp e i) = e := by
  apply bitExp_injective _ _ (digits_le_one _) he
  ext i
  exact congrFun (bitExp_digits (fun i => bitExp e i) (fun i => bitExp_lt e he i)) i

theorem expand_digits {k : ℕ} (b : Fin 2 → ℕ) (hb : ∀ i, b i < 2 ^ k) (c : K) :
    expand K 2 k (monomial (digits k b) c) = monomial (fs b) c := by
  rw [expand_monomial, bitExp_digits_fs b hb]

theorem support_expand {k : ℕ} (w : P (2 * k) K) (d : Fin 2 →₀ ℕ)
    (hd : d ∈ (expand K 2 k w).support) : ∃ e ∈ w.support, d = bitExp e := by
  classical
  have expansion : expand K 2 k w = ∑ e ∈ w.support, monomial (bitExp e) (coeff e w) := by
    conv_lhs => rw [w.as_sum]
    rw [map_sum]
    exact Finset.sum_congr rfl (fun e _ => expand_monomial 2 k e _)
  rw [expansion] at hd
  obtain ⟨e, he, hd⟩ := Finset.mem_biUnion.mp (support_sum hd)
  have single := support_monomial_subset hd
  rw [Finset.mem_singleton] at single
  exact ⟨e, he, single⟩

/-- The lift of a reduced chain polynomial of degree `≤ w` has exponents in the box and binary
weight `≤ w`: the degree of a squarefree chain monomial is the binary weight of its exponent. -/
theorem boxed_expand {k : ℕ} (w : P (2 * k) K) (hw : w ∈ reducedSpace (chainInputs K 2 k))
    (n : ℕ) (hdeg : w.totalDegree ≤ n) : Boxed k n (expand K 2 k w) := by
  intro d hd
  obtain ⟨e, he, rfl⟩ := support_expand w d hd
  have sq : ∀ v, e v ≤ 1 := fun v => standard_le_one 2 k e (hw e he) v
  have lt := bitExp_lt e sq
  refine ⟨lt 0, lt 1, ?_⟩
  have hdig : e.degree = wt (bitExp e 0) + wt (bitExp e 1) := by
    have h := digits_degree (k := k) (fun i => bitExp e i) lt
    rw [digits_bitExp e sq, Fin.sum_univ_two] at h
    exact h
  have hle := le_totalDegree he
  have hdeg' : e.degree = e.sum (fun _ n => n) := rfl
  omega

theorem X_pow_mul_monomial (i : Fin 2) (n : ℕ) (d : Fin 2 →₀ ℕ) (c : K) :
    (X i : Q K) ^ n * monomial d c = monomial (Finsupp.single i n + d) c := by
  rw [X_pow_eq_monomial, monomial_mul, one_mul]

/-- **The lift identity (1.1)**, cleared of the negative power: for a reduced chain polynomial
`w` with lift `f̃ = expand w`, the lift `T` of the reduced product `u_{i,j} · w` satisfies
`x_i^(N - 2^j) · T = x_i^N · (f̃ - f̃_hi) + f̃_hi`, `N = 2^k`. -/
theorem lift_identity (two_ne : (2 : K) ≠ 0) {k : ℕ} (w : P (2 * k) K)
    (hw : w ∈ reducedSpace (chainInputs K 2 k)) (i : Fin 2) (j : Fin k) :
    (X i : Q K) ^ (2 ^ k - 2 ^ (j : ℕ)) *
        expand K 2 k (ChainClosureBridge.normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k)
          (X (var 2 k i j) * w)) =
      X i ^ 2 ^ k * (expand K 2 k w - topPart k i j (expand K 2 k w)) +
        topPart k i j (expand K 2 k w) := by
  classical
  have zero := reduced_chain_ideal_trivial (K := K) 2 k two_ne
  let L₁ : P (2 * k) K →ₗ[K] Q K :=
    (LinearMap.mulLeft K ((X i : Q K) ^ (2 ^ k - 2 ^ (j : ℕ)))).comp
      ((expand K 2 k).toLinearMap.comp
        ((ChainClosureBridge.nfMap (chainInputs K 2 k) (chainInputs_nonzero 2 k) zero).comp
          (LinearMap.mulLeft K (X (var 2 k i j)))))
  let L₂ : P (2 * k) K →ₗ[K] Q K :=
    ((LinearMap.mulLeft K ((X i : Q K) ^ 2 ^ k)).comp (LinearMap.id - topPart k i j) +
      topPart k i j).comp (expand K 2 k).toLinearMap
  change L₁ w = L₂ w
  have s₁ : L₁ w = ∑ e ∈ w.support, L₁ (monomial e (coeff e w)) := by
    conv_lhs => rw [w.as_sum]
    rw [map_sum]
  have s₂ : L₂ w = ∑ e ∈ w.support, L₂ (monomial e (coeff e w)) := by
    conv_lhs => rw [w.as_sum]
    rw [map_sum]
  rw [s₁, s₂]
  apply Finset.sum_congr rfl
  intro e he
  have sq : ∀ v, e v ≤ 1 := fun v => standard_le_one 2 k e (hw e he) v
  set b : Fin 2 → ℕ := fun i => bitExp e i with hb
  have lt : ∀ i', b i' < 2 ^ k := bitExp_lt e sq
  have hdig : digits k b = e := digits_bitExp e sq
  set c := coeff e w
  have hjk : 2 ^ (j : ℕ) < 2 ^ k := Nat.pow_lt_pow_right (by norm_num) j.isLt
  have hmon : monomial e c = c • monomial (digits k b) (1 : K) := by
    rw [hdig, smul_monomial, smul_eq_mul, mul_one]
  have left : L₁ (monomial e c) =
      monomial (Finsupp.single i (2 ^ k - 2 ^ (j : ℕ)) +
        fs (Function.update b i ((b i + 2 ^ (j : ℕ)) % 2 ^ k))) c := by
    have hb' : ∀ i', Function.update b i ((b i + 2 ^ (j : ℕ)) % 2 ^ k) i' < 2 ^ k := by
      intro i'
      rcases eq_or_ne i' i with rfl | hne
      · rw [Function.update_self]
        exact Nat.mod_lt _ (Nat.two_pow_pos k)
      · rw [Function.update_of_ne hne]
        exact lt i'
    simp only [L₁, LinearMap.comp_apply, LinearMap.mulLeft_apply, ChainClosureBridge.nfMap_apply,
      AlgHom.toLinearMap_apply]
    rw [hmon, mul_smul_comm, ChainClosureBridge.normalForm_smul _ _ zero,
      normalForm_X_mul_monomial two_ne i j b lt, smul_monomial, smul_eq_mul, mul_one,
      expand_digits _ hb', X_pow_mul_monomial]
  have right : L₂ (monomial e c) =
      X i ^ 2 ^ k * (monomial (fs b) c - if 2 ^ k - 2 ^ (j : ℕ) ≤ b i then monomial (fs b) c else 0) +
        if 2 ^ k - 2 ^ (j : ℕ) ≤ b i then monomial (fs b) c else 0 := by
    simp only [L₂, LinearMap.comp_apply, LinearMap.add_apply, LinearMap.sub_apply,
      LinearMap.id_apply, LinearMap.mulLeft_apply, AlgHom.toLinearMap_apply]
    rw [← hdig, expand_digits b lt, topPart_monomial, fs_apply]
  rw [left, right]
  by_cases top : 2 ^ k - 2 ^ (j : ℕ) ≤ b i
  · rw [if_pos top, sub_self, mul_zero, zero_add]
    refine congrArg (fun d => monomial d c) (Finsupp.ext fun i' => ?_)
    rw [Finsupp.add_apply, fs_apply, fs_apply, Finsupp.single_apply]
    rcases eq_or_ne i i' with rfl | hne
    · rw [if_pos rfl, Function.update_self]
      have hsum : b i + 2 ^ (j : ℕ) = (b i + 2 ^ (j : ℕ) - 2 ^ k) + 2 ^ k := by omega
      rw [hsum, Nat.add_mod_right, Nat.mod_eq_of_lt (by have := lt i; omega)]
      omega
    · rw [if_neg hne, Function.update_of_ne (Ne.symm hne), zero_add]
  · rw [if_neg top, sub_zero, add_zero, X_pow_mul_monomial]
    refine congrArg (fun d => monomial d c) (Finsupp.ext fun i' => ?_)
    rw [Finsupp.add_apply, Finsupp.add_apply, fs_apply, fs_apply, Finsupp.single_apply,
      Finsupp.single_apply]
    rcases eq_or_ne i i' with rfl | hne
    · rw [if_pos rfl, if_pos rfl, Function.update_self, Nat.mod_eq_of_lt (by omega)]
      omega
    · rw [if_neg hne, if_neg hne, Function.update_of_ne (Ne.symm hne)]

/-! ## The no-wrap ideal space and its closedness -/

/-- **The no-wrap ideal space `S_D(k)`**: reduced chain polynomials of degree (binary weight)
at most `D` whose lift lies in `𝔭`. -/
def idealSpace (𝔭 : Ideal (Q K)) (k D : ℕ) : Submodule K (P (2 * k) K) where
  carrier := {w | w ∈ reducedSpace (chainInputs K 2 k) ∧ w.totalDegree ≤ D ∧ expand K 2 k w ∈ 𝔭}
  zero_mem' := ⟨(reducedSpace _).zero_mem, by simp, by simp⟩
  add_mem' := by
    rintro x y ⟨hx1, hx2, hx3⟩ ⟨hy1, hy2, hy3⟩
    refine ⟨(reducedSpace _).add_mem hx1 hy1, (totalDegree_add _ _).trans (max_le hx2 hy2), ?_⟩
    rw [map_add]
    exact 𝔭.add_mem hx3 hy3
  smul_mem' := by
    rintro c x ⟨hx1, hx2, hx3⟩
    refine ⟨(reducedSpace _).smul_mem c hx1, (totalDegree_smul_le c x).trans hx2, ?_⟩
    rw [map_smul, smul_eq_C_mul]
    exact 𝔭.mul_mem_left _ hx3

theorem mem_idealSpace {𝔭 : Ideal (Q K)} {k D : ℕ} {w : P (2 * k) K} :
    w ∈ idealSpace 𝔭 k D ↔
      w ∈ reducedSpace (chainInputs K 2 k) ∧ w.totalDegree ≤ D ∧ expand K 2 k w ∈ 𝔭 :=
  Iff.rfl

/-- **Proposition 2.1 (the direction used).** If `FEI_k(D-1)` holds for a prime `𝔭` with
`x_i ∉ 𝔭`, then `S_D(k)` is closed under the degree-`D` rule: reduced multiplication by a chain
variable of an element of degree `< D`. -/
theorem idealSpace_closed (two_ne : (2 : K) ≠ 0) (𝔭 : Ideal (Q K)) [𝔭.IsPrime]
    (hX : ∀ i : Fin 2, (X i : Q K) ∉ 𝔭) (k D : ℕ) (fei : FEI 𝔭 k (D - 1))
    (v : Fin (2 * k)) (w : P (2 * k) K) (hw : w ∈ idealSpace 𝔭 k D) (hdeg : w.totalDegree < D) :
    ChainClosureBridge.normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k) (X v * w) ∈
      idealSpace 𝔭 k D := by
  obtain ⟨hred, _, hmem⟩ := hw
  refine ⟨ChainClosureBridge.normalForm_reduced _ _ _, ?_, ?_⟩
  · apply (ChainClosureBridge.normalForm_totalDegree_le _ _ _).trans
    apply (totalDegree_mul _ _).trans
    rw [totalDegree_X]
    omega
  · have hv := var_row_depth v
    have id := lift_identity two_ne w hred (row v) (depth v)
    rw [hv] at id
    have boxed := boxed_expand w hred (D - 1) (by omega)
    have top := fei _ hmem boxed (row v) (depth v) (depth v).isLt
    have rhs : X (row v) ^ 2 ^ k * (expand K 2 k w - topPart k (row v) (depth v) (expand K 2 k w)) +
        topPart k (row v) (depth v) (expand K 2 k w) ∈ 𝔭 :=
      𝔭.add_mem (𝔭.mul_mem_left _ (𝔭.sub_mem hmem top)) top
    rw [← id] at rhs
    rcases Ideal.IsPrime.mem_or_mem inferInstance rhs with h | h
    · exact absurd (Ideal.IsPrime.mem_of_pow_mem inferInstance _ h) (hX (row v))
    · exact h

/-! ## The generator -/

/-- The third summation polynomial `S_3(x₁, x₂, t)` of `y² = x³ + a x + b`, in the two summand
coordinates (the formula of Section 1 of the written proof, and of `PrimeM2General.rawS3`). -/
def S3poly (a b t : K) : Q K :=
  (X 0 - X 1) ^ 2 * C (t ^ 2) - 2 * ((X 0 + X 1) * (X 0 * X 1 + C a) + C (2 * b)) * C t +
    (X 0 * X 1 - C a) ^ 2 - C (4 * b) * (X 0 + X 1)

theorem expand_active (k : ℕ) (depth : 2 ≤ k) :
    expand K 2 k (X (active k depth 0)) = X 0 ∧
    expand K 2 k (X (active k depth 1)) = X 0 ^ 2 ∧
    expand K 2 k (X (active k depth 2)) = X 1 ∧
    expand K 2 k (X (active k depth 3)) = X 1 ^ 2 := by
  have v00 : active k depth 0 = var 2 k 0 ⟨0, by omega⟩ := Fin.ext (by simp [var_val, active])
  have v01 : active k depth 1 = var 2 k 0 ⟨1, by omega⟩ := Fin.ext (by simp [var_val, active])
  have v10 : active k depth 2 = var 2 k 1 ⟨0, by omega⟩ := Fin.ext (by simp [var_val, active])
  have v11 : active k depth 3 = var 2 k 1 ⟨1, by omega⟩ :=
    Fin.ext (by simp [var_val, active]; omega)
  refine ⟨?_, ?_, ?_, ?_⟩
  · rw [v00, expand_X, row_var, depth_var]; simp
  · rw [v01, expand_X, row_var, depth_var]; simp
  · rw [v10, expand_X, row_var, depth_var]; simp
  · rw [v11, expand_X, row_var, depth_var]; simp

/-- The lift of the normalized generator is `S_3` (for every depth `k ≥ 2`). -/
theorem expand_normalized (k : ℕ) (depth : 2 ≤ k) (a b t : K) :
    expand K 2 k (normalized k depth a b t) = S3poly a b t := by
  obtain ⟨h0, h1, h2, h3⟩ := expand_active (K := K) k depth
  have id := congrArg (expand K 2 k) (raw_normalization_identity k depth a b t)
  simp only [map_sub, map_add, map_mul, map_pow, h0, h1, h2, h3, sub_self, mul_zero,
    add_zero] at id
  rw [← sub_eq_zero.mp id]
  simp only [rawS3, S3poly, map_sub, map_add, map_mul, map_pow, h0, h2, algHom_C, map_ofNat,
    MvPolynomial.algebraMap_eq]

/-! ## Proposition 2.2 -/

theorem normalized_mem_idealSpace (two_ne : (2 : K) ≠ 0) (k : ℕ) (depth : 2 ≤ k) (a b t : K)
    (𝔭 : Ideal (Q K)) (gen : S3poly a b t ∈ 𝔭) (D : ℕ)
    (hD : (normalized k depth a b t).totalDegree ≤ D) :
    ChainClosureBridge.normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k)
      (normalized k depth a b t) ∈ idealSpace 𝔭 k D := by
  rw [normalForm_normalized two_ne k depth a b t]
  refine ⟨normalized_reduced k depth a b t, hD, ?_⟩
  rw [expand_normalized]
  exact gen

theorem one_not_mem_idealSpace (k : ℕ) (𝔭 : Ideal (Q K)) [𝔭.IsPrime] (D : ℕ) :
    (1 : P (2 * k) K) ∉ idealSpace 𝔭 k D := by
  rintro ⟨_, _, h⟩
  rw [map_one] at h
  exact (Ideal.IsPrime.ne_top inferInstance) ((Ideal.eq_top_iff_one 𝔭).mpr h)

/-- **Proposition 2.2, `1 ∉ V_D(k)`.** If `𝔭` is a prime ideal of `K[x₁, x₂]` containing `S_3`
and neither coordinate, and `FEI_k(D-1)` holds, then `1` is not in the ordinary degree-`D`
closure of the two-summand chain presentation. -/
theorem not_one_of_fei (two_ne : (2 : K) ≠ 0) (k : ℕ) (depth : 2 ≤ k) (a b t : K)
    (𝔭 : Ideal (Q K)) [𝔭.IsPrime] (hX : ∀ i : Fin 2, (X i : Q K) ∉ 𝔭)
    (gen : S3poly a b t ∈ 𝔭) (D : ℕ) (fei : FEI 𝔭 k (D - 1)) :
    ¬ InClosure (Set.range (inputs k depth a b t)) D (1 : P (2 * k) K) :=
  ChainClosureBridge.prime_m2_not_one k depth a b t two_ne D (idealSpace 𝔭 k D)
    (normalized_mem_idealSpace two_ne k depth a b t 𝔭 gen D)
    (fun v w hw hdeg => idealSpace_closed two_ne 𝔭 hX k D fei v w hw hdeg)
    (one_not_mem_idealSpace k 𝔭 D)

/-- **Proposition 2.2, `V_D(k) ⊂ S_D(k)`.** Under the same hypotheses the chain normal form of
every element of the ordinary degree-`D` closure lies in the no-wrap ideal space: its lift
vanishes modulo `𝔭`. -/
theorem closure_normalForm_mem_idealSpace (two_ne : (2 : K) ≠ 0) (k : ℕ) (depth : 2 ≤ k)
    (a b t : K) (𝔭 : Ideal (Q K)) [𝔭.IsPrime] (hX : ∀ i : Fin 2, (X i : Q K) ∉ 𝔭)
    (gen : S3poly a b t ∈ 𝔭) (D : ℕ) (fei : FEI 𝔭 k (D - 1)) {p : P (2 * k) K}
    (h : InClosure (Set.range (inputs k depth a b t)) D p) :
    ChainClosureBridge.normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k) p ∈
      idealSpace 𝔭 k D := by
  have zero := reduced_chain_ideal_trivial (K := K) 2 k two_ne
  have hq := ChainClosureBridge.closure_normalForm_mem (nonzero := chainInputs_nonzero 2 k) zero h
  refine ChainClosureBridge.quot_closure_minimal _ D (idealSpace 𝔭 k D) ?_
    (fun v w hw hdeg => idealSpace_closed two_ne 𝔭 hX k D fei v w hw hdeg) hq
  rintro f ⟨i, rfl⟩ hd
  revert hd
  refine Fin.lastCases ?_ (fun j => ?_) i <;> intro hd
  · simp only [inputs, Fin.lastCases_last] at hd ⊢
    exact normalized_mem_idealSpace two_ne k depth a b t 𝔭 gen D hd
  · rw [ChainClosureBridge.normalForm_chain_eq_zero 2 k
      (ChainClosureBridge.inputs_castSucc_mem_chainInputs k depth a b t j) two_ne]
    exact (idealSpace 𝔭 k D).zero_mem

/-- Unit-ideal consequence: under `FEI_k(D-1)` no `e ≤ D` is a last fall degree bound. -/
theorem lfd_gt_of_fei (two_ne : (2 : K) ≠ 0) (k : ℕ) (depth : 2 ≤ k) (a b t : K)
    (𝔭 : Ideal (Q K)) [𝔭.IsPrime] (hX : ∀ i : Fin 2, (X i : Q K) ∉ 𝔭)
    (gen : S3poly a b t ∈ 𝔭) (D : ℕ) (fei : FEI 𝔭 k (D - 1))
    (unit : (1 : P (2 * k) K) ∈ Ideal.span (Set.range (inputs k depth a b t))) :
    ∀ e ≤ D, ¬ HasLFDBound (Set.range (inputs k depth a b t)) e := by
  intro e he bound
  have derived := one_of_hasLFDBound unit bound
  exact not_one_of_fei two_ne k depth a b t 𝔭 hX gen D fei (derived.mono he)

end ICFallDegrees.GrowthReduction
