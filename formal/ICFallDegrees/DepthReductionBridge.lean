import ICFallDegrees.ChainClosureBridge
import ICFallDegrees.DepthReduction

noncomputable section

/-!
# N5 assembled: base certificates ⇒ `1 ∉ V_D` at every depth `k ≥ k0+1`

This joins the two halves:

* `ChainClosureBridge.prime_m2_not_one`: a submodule of reduced polynomials that contains the
  normal form of the generator, is closed under reduced variable products on its elements of
  degree `< D`, and excludes `1`, refutes `InClosure (inputs k …) D 1` in the ordinary ring;
* `DepthReduction.W_closed` and `single_zero_mem_W_iff`: the free high extension `W k0 D S k`
  of a base `S` is closed at every depth `k ≥ k0+1` once it is closed at depth `k0+1`, and
  contains `1` iff `S` does.

The join is a linear map `Φ` from the exponent model at depth `k` into the polynomial ring,
with the three properties collected in `ModelMap`: it reflects `1`, it reflects the degree
filtration, and it intertwines the model translation `translate k0 k i j` with
`normalForm (u_{i,j} * ·)`.  Mathematically `Φ` is the monomial identification
`δ_(l,q) ↦ ∏ u_{i,j}^{bit_j(l_i + 2^k0 q_i)}`, and `ModelMap.step` is the statement that the
chain normal form of `u_{i,j}` times a squarefree monomial is the squarefree monomial of the
translated exponent.  That identification is *not* proved here: `ModelMap` is a hypothesis of
every theorem in this file, and discharging it (by constructing `Φ` and proving the carry
identity `u_{i,j} · m_e ≡ m_{e + 2^j e_i} (mod chains)`) is the remaining formal step.

With `k0 = 7`, `D = 4` this is N5 for the two-summand prime presentation: the base facts
`1 ∉ V_4(7)` and "`W_8` is closed" (the `829 × 3473` and `1061 × 2517` certificates) are the
hypotheses `one` and `base`, and the conclusion is `1 ∉ V_4(k)` for every `k ≥ 8`, hence
`LFD ≥ 5` for a unit ideal.
-/

namespace ICFallDegrees.DepthReductionBridge

open MvPolynomial FiniteQuotient ChainClosureBridge ChainQuotient PrimeM2General DepthReduction

variable {K : Type*} [Field K]

/-- The interface between the depth-`k` exponent model and the reduced polynomials of the
two-summand chain presentation: `Φ` reflects `1`, reflects the degree filtration on boxed
elements, and intertwines model translations with reduced variable products. -/
structure ModelMap (k0 k D : ℕ) (Φ : Al 2 K →ₗ[K] FiniteQuotient.P (2*k) K) : Prop where
  one : ∀ w, Boxed k0 k w → Φ w = 1 → w = Finsupp.single ((0 : Fin 2 → ℕ), (0 : Fin 2 → ℕ)) 1
  degree : ∀ w, Boxed k0 k w → (Φ w).totalDegree < D → WtLT D w
  step : ∀ w, Boxed k0 k w → ∀ (i : Fin 2) (j : Fin k),
    normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k) (X (var 2 k i j) * Φ w) =
      Φ (translate k0 k i (j : ℕ) w)

/-- The generator of the presentation lies in the base: some `v ∈ S` in the low box and within
the weight budget maps under `Φ` to the chain normal form of the normalized generator. -/
def GeneratorInBase (k0 D : ℕ) (S : Submodule K (Lo 2 K)) (k : ℕ) (hk2 : 2 ≤ k) (a b t : K)
    (Φ : Al 2 K →ₗ[K] FiniteQuotient.P (2*k) K) : Prop :=
  (normalized k hk2 a b t).totalDegree ≤ D →
    ∃ v ∈ S, (∀ l ∈ v.support, LowBox k0 l) ∧ (∀ l ∈ v.support, lowWt l ≤ D) ∧
      Φ (embed 0 v) = normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k) (normalized k hk2 a b t)

/-- **N5, assembled.** Base at depth `k0`: `1 ∉ S` and the free extension closed at depth
`k0+1`.  Then for every depth `k ≥ k0+1` and every faithful model map `Φ`, `1` is not in the
degree-`D` closure of the two-summand presentation. -/
theorem not_one_of_base (k0 D : ℕ) (S : Submodule K (Lo 2 K))
    (base : Closed k0 (k0 + 1) D (W k0 D S (k0 + 1)))
    (one : Finsupp.single (0 : Fin 2 → ℕ) (1 : K) ∉ S)
    (k : ℕ) (hk : k0 + 1 ≤ k) (hk2 : 2 ≤ k) (two_ne : (2 : K) ≠ 0) (a b t : K)
    (Φ : Al 2 K →ₗ[K] FiniteQuotient.P (2*k) K) (hΦ : ModelMap k0 k D Φ)
    (gen : GeneratorInBase k0 D S k hk2 a b t Φ) :
    ¬ InClosure (Set.range (inputs k hk2 a b t)) D (1 : FiniteQuotient.P (2*k) K) := by
  apply prime_m2_not_one k hk2 a b t two_ne D ((W k0 D S k).map Φ)
  · intro hd
    obtain ⟨v, hv, hvbox, hvwt, hΦv⟩ := gen hd
    exact ⟨embed 0 v, embed_zero_mem_W hv hvbox hvwt, hΦv⟩
  · rintro v' w' ⟨w, hw, rfl⟩ hdeg
    have hclosed := W_closed base k hk w hw (hΦ.degree w hw.1 hdeg) (ChainQuotient.row v')
      (ChainQuotient.depth v') (ChainQuotient.depth v').isLt
    refine ⟨translate k0 k (ChainQuotient.row v') (ChainQuotient.depth v') w, hclosed, ?_⟩
    have := hΦ.step w hw.1 (ChainQuotient.row v') (ChainQuotient.depth v')
    rw [var_row_depth] at this
    exact this.symm
  · rintro ⟨w, hw, hΦw⟩
    have := hΦ.one w hw.1 hΦw
    subst this
    exact one (single_zero_mem_W_iff.mp hw)

/-- Unit-ideal consequence: no `e ≤ D` is an LFD bound at any depth `k ≥ k0+1`. -/
theorem lfd_gt_of_base (k0 D : ℕ) (S : Submodule K (Lo 2 K))
    (base : Closed k0 (k0 + 1) D (W k0 D S (k0 + 1)))
    (one : Finsupp.single (0 : Fin 2 → ℕ) (1 : K) ∉ S)
    (k : ℕ) (hk : k0 + 1 ≤ k) (hk2 : 2 ≤ k) (two_ne : (2 : K) ≠ 0) (a b t : K)
    (Φ : Al 2 K →ₗ[K] FiniteQuotient.P (2*k) K) (hΦ : ModelMap k0 k D Φ)
    (gen : GeneratorInBase k0 D S k hk2 a b t Φ)
    (unit : (1 : FiniteQuotient.P (2*k) K) ∈ Ideal.span (Set.range (inputs k hk2 a b t))) :
    ∀ e ≤ D, ¬ HasLFDBound (Set.range (inputs k hk2 a b t)) e := by
  intro e he bound
  have derived := one_of_hasLFDBound unit bound
  exact not_one_of_base k0 D S base one k hk hk2 two_ne a b t Φ hΦ gen (derived.mono he)

/-- N5 with the depth-7 base and degree 4: a unit-ideal two-summand instance whose base
excludes `1` and whose depth-8 free extension is closed has `LFD ≥ 5` at every depth `k ≥ 8`. -/
theorem n5_lfd_ge_five (S : Submodule K (Lo 2 K))
    (base : Closed 7 8 4 (W 7 4 S 8))
    (one : Finsupp.single (0 : Fin 2 → ℕ) (1 : K) ∉ S)
    (k : ℕ) (hk : 8 ≤ k) (hk2 : 2 ≤ k) (two_ne : (2 : K) ≠ 0) (a b t : K)
    (Φ : Al 2 K →ₗ[K] FiniteQuotient.P (2*k) K) (hΦ : ModelMap 7 k 4 Φ)
    (gen : GeneratorInBase 7 4 S k hk2 a b t Φ)
    (unit : (1 : FiniteQuotient.P (2*k) K) ∈ Ideal.span (Set.range (inputs k hk2 a b t))) :
    ∀ e < 5, ¬ HasLFDBound (Set.range (inputs k hk2 a b t)) e := by
  intro e he
  exact lfd_gt_of_base 7 4 S base one k hk hk2 two_ne a b t Φ hΦ gen unit e (by omega)

/-- Prime-field entry point of the assembled statement. -/
theorem prime_field_n5_lfd_ge_five (p : ℕ) [Fact p.Prime] (large : 5 ≤ p)
    (S : Submodule (ZMod p) (Lo 2 (ZMod p)))
    (base : Closed 7 8 4 (W 7 4 S 8))
    (one : Finsupp.single (0 : Fin 2 → ℕ) (1 : ZMod p) ∉ S)
    (k : ℕ) (hk : 8 ≤ k) (hk2 : 2 ≤ k) (a b t : ZMod p)
    (Φ : Al 2 (ZMod p) →ₗ[ZMod p] FiniteQuotient.P (2*k) (ZMod p)) (hΦ : ModelMap 7 k 4 Φ)
    (gen : GeneratorInBase 7 4 S k hk2 a b t Φ)
    (unit : (1 : FiniteQuotient.P (2*k) (ZMod p)) ∈
      Ideal.span (Set.range (inputs k hk2 a b t))) :
    ∀ e < 5, ¬ HasLFDBound (Set.range (inputs k hk2 a b t)) e :=
  n5_lfd_ge_five S base one k hk hk2 (PrimeFieldConclusions.two_ne_zero p large) a b t Φ hΦ gen
    unit

end ICFallDegrees.DepthReductionBridge
