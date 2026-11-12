import ICFallDegrees.ChainQuotient
import ICFallDegrees.QuotientLowerBound
import ICFallDegrees.VariableClosure

noncomputable section

/-!
# The ordinary/quotient closure bridge

Let `G` be a reduced-ideal-trivial set of nonzero polynomials (for example the
chain equations of `ChainQuotient`).  Graded-lex division gives every
polynomial `p` a normal form `normalForm G p` in the reduced space of `G`, with
`p - normalForm G p` derivable from `G` inside the degree budget `p.totalDegree`
(`FiniteQuotient.exists_reduced_remainder`), and reduced-ideal triviality makes
that normal form unique and `K`-linear.

The *quotient closure* `InQuotClosure G F D` is the closure computed in the
reduced space: it contains the normal forms of the inputs of degree `≤ D` and
is closed under `w ↦ normalForm G (X v * w)` for its elements of degree `< D`.
This is the "closure modulo the chains" that the computational pipeline
(`csd/ic_bounds_fall.py`, `check_containing_space`) works with, once squarefree
chain monomials are read as exponent pairs.

The bridge proved here, for `G ⊆ F` with all of `G` of degree `≤ D`:

* `mem_closure_iff`: `InClosure F D p ↔ p.totalDegree ≤ D ∧ InQuotClosure G F D (normalForm G p)`,
  so the ordinary closure is the sum of the chain part `{c ∈ ⟨G⟩ : deg c ≤ D}` and the
  quotient closure, and the two have only `0` in common;
* `one_mem_iff`: `1` is in the ordinary closure iff it is in the quotient closure;
* `not_one_of_witness`: a submodule `W` that contains the normal forms of the inputs,
  is closed under `normalForm G (X v * ·)` on its elements of degree `< D`, and excludes
  `1`, refutes `InClosure F D 1` in the ordinary ring (this half needs neither `G ⊆ F`
  nor a degree bound on `G`), hence `LFD > D` for a unit ideal.

Everything is then specialised to the chain equations and to the normalized
two-summand presentation `PrimeM2General.inputs`.
-/

namespace ICFallDegrees.ChainClosureBridge

open MvPolynomial FiniteQuotient

variable {r : ℕ} {K : Type*} [Field K]

section NormalForm

variable (G : Set (P r K)) (nonzero : (0 : P r K) ∉ G)

/-- The normal form modulo `G` supplied by graded-lex division. -/
def normalForm (p : P r K) : P r K :=
  Classical.choose (exists_reduced_remainder G nonzero p)

theorem normalForm_reduced (p : P r K) : normalForm G nonzero p ∈ reducedSpace G :=
  (Classical.choose_spec (exists_reduced_remainder G nonzero p)).1

theorem normalForm_reduction (p : P r K) :
    InClosure G p.totalDegree (p - normalForm G nonzero p) :=
  (Classical.choose_spec (exists_reduced_remainder G nonzero p)).2

theorem sub_normalForm_mem (p : P r K) : p - normalForm G nonzero p ∈ Ideal.span G :=
  (normalForm_reduction G nonzero p).sound

theorem normalForm_totalDegree_le (p : P r K) :
    (normalForm G nonzero p).totalDegree ≤ p.totalDegree :=
  remainder_degree_le (normalForm_reduction G nonzero p)

variable (zero : ReducedIdealTrivial G)
include zero

/-- Uniqueness: any reduced polynomial congruent to `p` is its normal form. -/
theorem normalForm_eq_of {p z : P r K} (hz : z ∈ reducedSpace G) (hpz : p - z ∈ Ideal.span G) :
    normalForm G nonzero p = z := by
  apply reduced_unique_of_reduced_zero G zero (normalForm_reduced G nonzero p) hz
  have identity : normalForm G nonzero p - z = (p - z) - (p - normalForm G nonzero p) := by ring
  rw [identity]
  exact (Ideal.span G).sub_mem hpz (sub_normalForm_mem G nonzero p)

theorem normalForm_of_reduced {p : P r K} (hp : p ∈ reducedSpace G) :
    normalForm G nonzero p = p :=
  normalForm_eq_of G nonzero zero hp (by simp)

theorem normalForm_of_mem {p : P r K} (hp : p ∈ Ideal.span G) : normalForm G nonzero p = 0 :=
  normalForm_eq_of G nonzero zero (reducedSpace G).zero_mem (by simpa using hp)

theorem normalForm_zero : normalForm G nonzero (0 : P r K) = 0 :=
  normalForm_of_mem G nonzero zero (Ideal.span G).zero_mem

theorem normalForm_add (p q : P r K) :
    normalForm G nonzero (p + q) = normalForm G nonzero p + normalForm G nonzero q := by
  apply normalForm_eq_of G nonzero zero
  · exact (reducedSpace G).add_mem (normalForm_reduced G nonzero p) (normalForm_reduced G nonzero q)
  · have identity : p + q - (normalForm G nonzero p + normalForm G nonzero q) =
        (p - normalForm G nonzero p) + (q - normalForm G nonzero q) := by ring
    rw [identity]
    exact (Ideal.span G).add_mem (sub_normalForm_mem G nonzero p) (sub_normalForm_mem G nonzero q)

theorem normalForm_smul (c : K) (p : P r K) :
    normalForm G nonzero (c • p) = c • normalForm G nonzero p := by
  apply normalForm_eq_of G nonzero zero
  · exact (reducedSpace G).smul_mem c (normalForm_reduced G nonzero p)
  · have identity : c • p - c • normalForm G nonzero p = C c * (p - normalForm G nonzero p) := by
      simp only [smul_eq_C_mul]
      ring
    rw [identity]
    exact (Ideal.span G).mul_mem_left _ (sub_normalForm_mem G nonzero p)

theorem normalForm_sub (p q : P r K) :
    normalForm G nonzero (p - q) = normalForm G nonzero p - normalForm G nonzero q := by
  have h := normalForm_add G nonzero zero (p - q) q
  rw [sub_add_cancel] at h
  rw [h]
  ring

/-- Reduction commutes with multiplication: `a * p` and `a * normalForm p` have the same
normal form. -/
theorem normalForm_mul_normalForm (a p : P r K) :
    normalForm G nonzero (a * p) = normalForm G nonzero (a * normalForm G nonzero p) := by
  apply normalForm_eq_of G nonzero zero (normalForm_reduced G nonzero _)
  have identity : a * p - normalForm G nonzero (a * normalForm G nonzero p) =
      a * (p - normalForm G nonzero p) +
        (a * normalForm G nonzero p - normalForm G nonzero (a * normalForm G nonzero p)) := by ring
  rw [identity]
  exact (Ideal.span G).add_mem ((Ideal.span G).mul_mem_left a (sub_normalForm_mem G nonzero p))
    (sub_normalForm_mem G nonzero _)

theorem normalForm_one (one_reduced : (1 : P r K) ∈ reducedSpace G) :
    normalForm G nonzero (1 : P r K) = 1 :=
  normalForm_of_reduced G nonzero zero one_reduced

/-- The normal form as a `K`-linear map. -/
def nfMap : P r K →ₗ[K] P r K where
  toFun := normalForm G nonzero
  map_add' := normalForm_add G nonzero zero
  map_smul' := normalForm_smul G nonzero zero

@[simp] theorem nfMap_apply (p : P r K) : nfMap G nonzero zero p = normalForm G nonzero p := rfl

end NormalForm

section QuotientClosure

variable (G : Set (P r K)) (nonzero : (0 : P r K) ∉ G)

/-- The closure computed in the reduced space of `G`: normal forms of the inputs,
linear combinations, and `normalForm (X v * w)` for elements `w` of degree below `D`. -/
inductive InQuotClosure (F : Set (P r K)) (D : ℕ) : P r K → Prop
  | input {f} (member : f ∈ F) (degree : f.totalDegree ≤ D) :
      InQuotClosure F D (normalForm G nonzero f)
  | zero : InQuotClosure F D 0
  | add {f g} : InQuotClosure F D f → InQuotClosure F D g → InQuotClosure F D (f + g)
  | smul (c : K) {f} : InQuotClosure F D f → InQuotClosure F D (c • f)
  | mul (v : Fin r) {w} : InQuotClosure F D w → w.totalDegree < D →
      InQuotClosure F D (normalForm G nonzero (X v * w))

variable {G nonzero}

namespace InQuotClosure

variable {F : Set (P r K)} {D : ℕ}

theorem reduced {w : P r K} (h : InQuotClosure G nonzero F D w) : w ∈ reducedSpace G := by
  induction h with
  | input _ _ => exact normalForm_reduced G nonzero _
  | zero => exact (reducedSpace G).zero_mem
  | add _ _ ih₁ ih₂ => exact (reducedSpace G).add_mem ih₁ ih₂
  | smul c _ ih => exact (reducedSpace G).smul_mem c ih
  | mul _ _ _ _ => exact normalForm_reduced G nonzero _

theorem degree_le {w : P r K} (h : InQuotClosure G nonzero F D w) : w.totalDegree ≤ D := by
  induction h with
  | input _ hd => exact (normalForm_totalDegree_le G nonzero _).trans hd
  | zero => simp
  | add _ _ ih₁ ih₂ => exact (totalDegree_add _ _).trans (max_le ih₁ ih₂)
  | smul c _ ih => exact (totalDegree_smul_le c _).trans ih
  | mul v _ low _ =>
      apply (normalForm_totalDegree_le G nonzero _).trans
      apply (totalDegree_mul _ _).trans
      rw [totalDegree_X]
      omega

/-- Quotient-closure elements are ordinary-closure elements, once the reducing
system is among the inputs within the budget. -/
theorem lift (included : G ⊆ F) (gdeg : ∀ g ∈ G, g.totalDegree ≤ D) {w : P r K}
    (h : InQuotClosure G nonzero F D w) : InClosure F D w := by
  have transport : ∀ p : P r K, p.totalDegree ≤ D →
      InClosure F D (p - normalForm G nonzero p) := by
    intro p hp
    exact InClosure.transport hp (fun g hg _ => .input (included hg) (gdeg g hg))
      (normalForm_reduction G nonzero p)
  induction h with
  | input member hd =>
      have h₁ : InClosure F D _ := .input member hd
      simpa using h₁.sub (transport _ hd)
  | zero => exact .zero
  | add _ _ ih₁ ih₂ => exact ih₁.add ih₂
  | smul c _ ih => exact ih.smul c
  | @mul v w _ low ih =>
      have budget : (X v : P r K).totalDegree + w.totalDegree ≤ D := by rw [totalDegree_X]; omega
      have h₁ : InClosure F D (X v * w) := ih.mul (X v) budget
      have hd : (X v * w : P r K).totalDegree ≤ D := h₁.degree_le
      simpa using h₁.sub (transport _ hd)

end InQuotClosure

/-- The quotient closure as a submodule. -/
def quotClosureSpace (F : Set (P r K)) (D : ℕ) : Submodule K (P r K) where
  carrier := InQuotClosure G nonzero F D
  zero_mem' := .zero
  add_mem' := InQuotClosure.add
  smul_mem' := InQuotClosure.smul

theorem mem_quotClosureSpace {F : Set (P r K)} {D : ℕ} {w : P r K} :
    w ∈ quotClosureSpace (G := G) (nonzero := nonzero) F D ↔ InQuotClosure G nonzero F D w :=
  Iff.rfl

/-- Least-space property of the quotient closure. -/
theorem quot_closure_minimal (F : Set (P r K)) (D : ℕ) (W : Submodule K (P r K))
    (inputs : ∀ f ∈ F, f.totalDegree ≤ D → normalForm G nonzero f ∈ W)
    (products : ∀ (v : Fin r) (w : P r K), w ∈ W → w.totalDegree < D →
      normalForm G nonzero (X v * w) ∈ W)
    {w : P r K} (h : InQuotClosure G nonzero F D w) : w ∈ W := by
  induction h with
  | input member degree => exact inputs _ member degree
  | zero => exact W.zero_mem
  | add _ _ ih₁ ih₂ => exact W.add_mem ih₁ ih₂
  | smul c _ ih => exact W.smul_mem c ih
  | mul v _ low ih => exact products v _ ih low

variable (zero : ReducedIdealTrivial G)
include zero

/-- Multiplying a quotient-closure element by any polynomial within the budget and
reducing stays in the quotient closure: the variable rule generates the multiplier rule. -/
theorem quot_mul_mem {F : Set (P r K)} {D : ℕ} (a : P r K) {w : P r K}
    (hw : InQuotClosure G nonzero F D w) (budget : a.totalDegree + w.totalDegree ≤ D) :
    InQuotClosure G nonzero F D (normalForm G nonzero (a * w)) := by
  classical
  have monomials : ∀ e : Fin r →₀ ℕ, e.degree + w.totalDegree ≤ D →
      InQuotClosure G nonzero F D (normalForm G nonzero (monomial e 1 * w)) := by
    intro e
    induction e using Finsupp.induction with
    | zero =>
        intro _
        rw [show (monomial (0 : Fin r →₀ ℕ) (1 : K)) = 1 from rfl, one_mul,
          normalForm_of_reduced G nonzero zero hw.reduced]
        exact hw
    | single_add i n e _ _ ih =>
        intro bound
        have degrees : (Finsupp.single i n + e).degree = n + e.degree := by simp
        rw [degrees] at bound
        have base := ih (by omega)
        have powers : ∀ j : ℕ, j + e.degree + w.totalDegree ≤ D →
            InQuotClosure G nonzero F D (normalForm G nonzero (X i ^ j * (monomial e 1 * w))) := by
          intro j
          induction j with
          | zero => intro _; simpa using base
          | succ j hj =>
              intro hbound
              have before := hj (by omega)
              have degree : (normalForm G nonzero (X i ^ j * (monomial e 1 * w))).totalDegree <
                  D := by
                apply lt_of_le_of_lt (normalForm_totalDegree_le G nonzero _)
                apply lt_of_le_of_lt (totalDegree_mul _ _)
                have dp := totalDegree_mul (monomial e (1 : K)) w
                have de : (monomial e (1 : K)).totalDegree ≤ e.degree :=
                  totalDegree_monomial_le e (1 : K)
                have dx : (X i ^ j : P r K).totalDegree = j := by simp
                omega
              have derived := InQuotClosure.mul (G := G) (nonzero := nonzero) i before degree
              rw [← normalForm_mul_normalForm G nonzero zero] at derived
              simpa [pow_succ, mul_assoc, mul_left_comm] using derived
        have derived := powers n bound
        simpa only [X_pow_eq_monomial, ← mul_assoc, monomial_mul, mul_one] using derived
  have expansion : normalForm G nonzero (a * w) =
      ∑ e ∈ a.support, coeff e a • normalForm G nonzero (monomial e 1 * w) := by
    conv_lhs => rw [a.as_sum, Finset.sum_mul]
    rw [← nfMap_apply G nonzero zero, map_sum]
    apply Finset.sum_congr rfl
    intro e _
    rw [nfMap_apply, ← normalForm_smul G nonzero zero]
    congr 1
    rw [smul_eq_C_mul, ← mul_assoc, C_mul_monomial, mul_one]
  rw [expansion]
  apply (quotClosureSpace (G := G) (nonzero := nonzero) F D).sum_mem
  intro e he
  apply (quotClosureSpace (G := G) (nonzero := nonzero) F D).smul_mem
  have edeg : e.degree ≤ a.totalDegree := le_totalDegree he
  exact monomials e (by omega)

/-- Every ordinary-closure element reduces into the quotient closure. No relation between
`G` and `F` is needed for this direction. -/
theorem closure_normalForm_mem {F : Set (P r K)} {D : ℕ} {p : P r K} (h : InClosure F D p) :
    InQuotClosure G nonzero F D (normalForm G nonzero p) := by
  induction h with
  | input member degree => exact .input member degree
  | zero => rw [normalForm_zero G nonzero zero]; exact .zero
  | add _ _ ih₁ ih₂ => rw [normalForm_add G nonzero zero]; exact ih₁.add ih₂
  | smul c _ ih => rw [normalForm_smul G nonzero zero]; exact ih.smul c
  | @mul a f _ budget ih =>
      rw [normalForm_mul_normalForm G nonzero zero]
      apply quot_mul_mem zero a ih
      have := normalForm_totalDegree_le G nonzero f
      omega

/-- The bridge: ordinary closure membership is a degree bound on the chain part plus
quotient closure membership of the normal form. -/
theorem mem_closure_iff (F : Set (P r K)) (D : ℕ) (included : G ⊆ F)
    (gdeg : ∀ g ∈ G, g.totalDegree ≤ D) (p : P r K) :
    InClosure F D p ↔ p.totalDegree ≤ D ∧ InQuotClosure G nonzero F D (normalForm G nonzero p) := by
  constructor
  · intro h
    exact ⟨h.degree_le, closure_normalForm_mem zero h⟩
  · rintro ⟨hd, hq⟩
    have reduced := hq.lift included gdeg
    have chain : InClosure F D (p - normalForm G nonzero p) :=
      InClosure.transport hd (fun g hg _ => .input (included hg) (gdeg g hg))
        (normalForm_reduction G nonzero p)
    simpa using chain.add reduced

/-- For reduced polynomials the two closures coincide. -/
theorem mem_closure_iff_of_reduced (F : Set (P r K)) (D : ℕ) (included : G ⊆ F)
    (gdeg : ∀ g ∈ G, g.totalDegree ≤ D) {p : P r K} (hp : p ∈ reducedSpace G) :
    InClosure F D p ↔ InQuotClosure G nonzero F D p := by
  rw [mem_closure_iff zero F D included gdeg, normalForm_of_reduced G nonzero zero hp]
  constructor
  · exact fun h => h.2
  · exact fun h => ⟨h.degree_le, h⟩

theorem one_mem_iff (F : Set (P r K)) (D : ℕ) (included : G ⊆ F)
    (gdeg : ∀ g ∈ G, g.totalDegree ≤ D) (one_reduced : (1 : P r K) ∈ reducedSpace G) :
    InClosure F D (1 : P r K) ↔ InQuotClosure G nonzero F D 1 :=
  mem_closure_iff_of_reduced zero F D included gdeg one_reduced

/-- A quotient-side containing space excluding `1` refutes `1` in the ordinary closure. -/
theorem not_one_of_witness (F : Set (P r K)) (D : ℕ) (one_reduced : (1 : P r K) ∈ reducedSpace G)
    (W : Submodule K (P r K))
    (inputs : ∀ f ∈ F, f.totalDegree ≤ D → normalForm G nonzero f ∈ W)
    (products : ∀ (v : Fin r) (w : P r K), w ∈ W → w.totalDegree < D →
      normalForm G nonzero (X v * w) ∈ W)
    (one : (1 : P r K) ∉ W) : ¬ InClosure F D (1 : P r K) := by
  intro h
  have hq := closure_normalForm_mem (nonzero := nonzero) zero h
  rw [normalForm_one G nonzero zero one_reduced] at hq
  exact one (quot_closure_minimal F D W inputs products hq)

/-- Unit ideal plus a quotient-side witness at degree `D`: no `e ≤ D` is an LFD bound. -/
theorem lfd_gt_of_witness (F : Set (P r K)) (D : ℕ) (one_reduced : (1 : P r K) ∈ reducedSpace G)
    (W : Submodule K (P r K))
    (inputs : ∀ f ∈ F, f.totalDegree ≤ D → normalForm G nonzero f ∈ W)
    (products : ∀ (v : Fin r) (w : P r K), w ∈ W → w.totalDegree < D →
      normalForm G nonzero (X v * w) ∈ W)
    (one : (1 : P r K) ∉ W) (unit : (1 : P r K) ∈ Ideal.span F) :
    ∀ e ≤ D, ¬ HasLFDBound F e := by
  intro e he bound
  have derived := one_of_hasLFDBound unit bound
  exact not_one_of_witness zero F D one_reduced W inputs products one (derived.mono he)

end QuotientClosure

/-! ## Specialisation to the chain equations -/

section Chains

open ChainQuotient

variable (m k : ℕ)

theorem chain_one_reduced : (1 : P (m*k) K) ∈ reducedSpace (chainInputs K m k) := by
  classical
  intro e he g hg
  obtain ⟨v, rfl⟩ := hg
  rw [leading_chainPolynomial, Finsupp.single_le_iff]
  have : e = 0 := by
    rw [show (1 : P (m*k) K) = monomial 0 1 from C_1.symm.trans C_apply] at he
    rw [support_monomial, if_neg one_ne_zero, Finset.mem_singleton] at he
    exact he
  rw [this]
  simp

theorem normalForm_chain_eq_zero {f : P (m*k) K} (hf : f ∈ chainInputs K m k)
    (two_ne : (2 : K) ≠ 0) :
    normalForm (chainInputs K m k) (chainInputs_nonzero m k) f = 0 :=
  normalForm_of_mem _ _ (reduced_chain_ideal_trivial m k two_ne) (Ideal.subset_span hf)

/-- The witness theorem for any input list containing the chain equations of `m` summands at
depth `k`: only the non-chain inputs need their normal forms in `W`. -/
theorem not_one_of_chain_witness (two_ne : (2 : K) ≠ 0) (F : Set (P (m*k) K)) (D : ℕ)
    (W : Submodule K (P (m*k) K))
    (inputs : ∀ f ∈ F, f ∉ chainInputs K m k → f.totalDegree ≤ D →
      normalForm (chainInputs K m k) (chainInputs_nonzero m k) f ∈ W)
    (products : ∀ (v : Fin (m*k)) (w : P (m*k) K), w ∈ W → w.totalDegree < D →
      normalForm (chainInputs K m k) (chainInputs_nonzero m k) (X v * w) ∈ W)
    (one : (1 : P (m*k) K) ∉ W) : ¬ InClosure F D (1 : P (m*k) K) := by
  apply not_one_of_witness (reduced_chain_ideal_trivial m k two_ne) F D (chain_one_reduced m k) W
    _ products one
  intro f hf hd
  by_cases chain : f ∈ chainInputs K m k
  · rw [normalForm_chain_eq_zero m k chain two_ne]
    exact W.zero_mem
  · exact inputs f hf chain hd

/-- The exact bridge for the chain presentation: for `D ≥ 2` and any further inputs `R`,
ordinary closure membership is "chain part of degree `≤ D`" plus quotient closure membership. -/
theorem chain_mem_closure_iff (two_ne : (2 : K) ≠ 0) (R : Set (P (m*k) K)) (D : ℕ) (hD : 2 ≤ D)
    (p : P (m*k) K) :
    InClosure (chainInputs K m k ∪ R) D p ↔
      p.totalDegree ≤ D ∧ InQuotClosure (chainInputs K m k) (chainInputs_nonzero m k)
        (chainInputs K m k ∪ R) D (normalForm (chainInputs K m k) (chainInputs_nonzero m k) p) := by
  apply mem_closure_iff (reduced_chain_ideal_trivial m k two_ne) _ D Set.subset_union_left
  rintro g ⟨v, rfl⟩
  exact (chainPolynomial_totalDegree_le m k v).trans hD

end Chains

/-! ## The normalized two-summand presentation -/

section PrimeM2

open ChainQuotient PrimeM2General

variable (k : ℕ) (depth : 2 ≤ k) (a b t : K)

/-- Every chain input of the two-summand presentation is a chain polynomial of `ChainQuotient`. -/
theorem inputs_castSucc_mem_chainInputs (j : Fin (2*k)) :
    inputs k depth a b t (Fin.castSucc j) ∈ chainInputs K 2 k := by
  refine ⟨j, ?_⟩
  simp only [inputs, Fin.lastCases_castSucc]
  unfold chainPolynomial ChainQuotient.chainTail PrimeM2General.chainTail
  have hj := j.isLt
  have decomposition := Nat.mod_add_div (j : ℕ) k
  have remainder := Nat.mod_lt (j : ℕ) (by omega : 0 < k)
  have quotient : (j : ℕ) / k < 2 := (Nat.div_lt_iff_lt_mul (by omega)).mpr hj
  rw [depth_val]
  have cases : (j : ℕ) / k = 0 ∨ (j : ℕ) / k = 1 := by
    rcases (Nat.lt_succ_iff.mp quotient).lt_or_eq with q | q
    · exact Or.inl (Nat.lt_one_iff.mp q)
    · exact Or.inr q
  split_ifs with h₁ h₂ h₂
  · exfalso
    rcases cases with q | q <;> rw [q] at decomposition <;>
      simp only [Nat.mul_zero, Nat.mul_one, Nat.add_zero] at decomposition <;> omega
  · rw [sub_eq_add_neg]
    congr 2
    congr 1
    apply Fin.ext
    rw [var_val, row_val]
    simp only [Fin.val_mk]
    linarith [Nat.mod_add_div (j : ℕ) k]
  · rw [sub_eq_add_neg]
  · exfalso
    rcases cases with q | q <;> rw [q] at decomposition <;>
      simp only [Nat.mul_zero, Nat.mul_one, Nat.add_zero] at decomposition <;> omega

/-- Quotient-side witness for the normalized two-summand presentation at every depth: a
submodule containing the normal form of the generator, closed under reduced variable products
on its elements of degree `< D`, and excluding `1`. -/
theorem prime_m2_not_one (two_ne : (2 : K) ≠ 0) (D : ℕ)
    (W : Submodule K (FiniteQuotient.P (2*k) K))
    (gen : (normalized k depth a b t).totalDegree ≤ D →
      normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k) (normalized k depth a b t) ∈ W)
    (products : ∀ (v : Fin (2*k)) (w : FiniteQuotient.P (2*k) K), w ∈ W → w.totalDegree < D →
      normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k) (X v * w) ∈ W)
    (one : (1 : FiniteQuotient.P (2*k) K) ∉ W) :
    ¬ InClosure (Set.range (inputs k depth a b t)) D (1 : FiniteQuotient.P (2*k) K) := by
  apply not_one_of_chain_witness 2 k two_ne _ D W _ products one
  rintro f ⟨i, rfl⟩ notchain hd
  revert notchain hd
  refine Fin.lastCases ?_ (fun j => ?_) i <;> intro notchain hd
  · simp only [inputs, Fin.lastCases_last] at hd ⊢
    exact gen hd
  · exact absurd (inputs_castSucc_mem_chainInputs k depth a b t j) notchain

/-- Unit-ideal consequence: no `e ≤ D` is an LFD bound. -/
theorem prime_m2_lfd_gt (two_ne : (2 : K) ≠ 0) (D : ℕ)
    (W : Submodule K (FiniteQuotient.P (2*k) K))
    (gen : (normalized k depth a b t).totalDegree ≤ D →
      normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k) (normalized k depth a b t) ∈ W)
    (products : ∀ (v : Fin (2*k)) (w : FiniteQuotient.P (2*k) K), w ∈ W → w.totalDegree < D →
      normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k) (X v * w) ∈ W)
    (one : (1 : FiniteQuotient.P (2*k) K) ∉ W)
    (unit : (1 : FiniteQuotient.P (2*k) K) ∈ Ideal.span (Set.range (inputs k depth a b t))) :
    ∀ e ≤ D, ¬ HasLFDBound (Set.range (inputs k depth a b t)) e := by
  intro e he bound
  have derived := one_of_hasLFDBound unit bound
  exact prime_m2_not_one k depth a b t two_ne D W gen products one (derived.mono he)

end PrimeM2

end ICFallDegrees.ChainClosureBridge
