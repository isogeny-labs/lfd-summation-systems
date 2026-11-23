import ICFallDegrees.SummandGrowthTransfer
import ICFallDegrees.GrowthReduction

/-!
Growth for every number of summands through a two-summand slice.

`not_one_of_fei_inputs` is the exclusion of `1` from the growth reduction with the
generator replaced by an arbitrary input set: only membership of the non-chain inputs
of degree at most `D` in the ideal space is used. `slice_not_one` composes it with the
degree-controlled transfer of `SummandGrowthTransfer` (preprocessing cost `0`).

Not asserted here: that a particular substitution maps the m-summand direct presentation
into such an input set. For Semaev polynomials that needs the defining zero property of
`S_{m+1}`, which is not formalized; the paper's Corollary `cor:growthm` supplies it in
writing. The splitting hypothesis `FEI` is the same explicit hypothesis as in
`GrowthReduction.not_one_of_fei`.
-/
noncomputable section
namespace ICFallDegrees.GrowthSlice

open MvPolynomial FiniteQuotient ChainQuotient GrowthReduction

variable {K : Type*} [Field K]

/-- Generator independence of the growth reduction. -/
theorem not_one_of_fei_inputs (two_ne : (2 : K) ≠ 0) (k : ℕ) (𝔭 : Ideal (Q K)) [𝔭.IsPrime]
    (hX : ∀ i : Fin 2, (X i : Q K) ∉ 𝔭) (D : ℕ) (fei : FEI 𝔭 k (D - 1))
    (F : Set (P (2 * k) K))
    (inputs : ∀ f ∈ F, f ∉ chainInputs K 2 k → f.totalDegree ≤ D →
      ChainClosureBridge.normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k) f ∈
        idealSpace 𝔭 k D) :
    ¬ InClosure F D (1 : P (2 * k) K) :=
  ChainClosureBridge.not_one_of_chain_witness 2 k two_ne F D (idealSpace 𝔭 k D) inputs
    (fun v w hw hdeg => idealSpace_closed two_ne 𝔭 hX k D fei v w hw hdeg)
    (one_not_mem_idealSpace k 𝔭 D)

/-- The slice: a substitution that does not raise degrees and sends every source input of
degree at most `D` to a chain input of two blocks or to an element whose normal form lies in
the ideal space excludes `1` from the degree-`D` closure of the source. -/
theorem slice_not_one {σ : Type*} (two_ne : (2 : K) ≠ 0) (k : ℕ) (𝔭 : Ideal (Q K))
    [𝔭.IsPrime] (hX : ∀ i : Fin 2, (X i : Q K) ∉ 𝔭) (D : ℕ) (fei : FEI 𝔭 k (D - 1))
    (φ : MvPolynomial σ K →ₐ[K] P (2 * k) K) (linear : ∀ i, (φ (X i)).totalDegree ≤ 1)
    (F : Set (MvPolynomial σ K))
    (slice : ∀ f ∈ F, f.totalDegree ≤ D → φ f ∉ chainInputs K 2 k →
      ChainClosureBridge.normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k) (φ f) ∈
        idealSpace 𝔭 k D) :
    ¬ InClosure F D (1 : MvPolynomial σ K) := by
  intro h
  let G : Set (P (2 * k) K) := {g | ∃ f ∈ F, f.totalDegree ≤ D ∧ φ f = g}
  have hG : InClosure G D (φ 1) :=
    SummandGrowthTransfer.map_closure_budget φ linear F G D D le_rfl
      (fun f hf hd => InClosure.input ⟨f, hf, hd, rfl⟩
        ((AffineTransport.degree_le φ linear f).trans hd)) h
  rw [map_one] at hG
  refine not_one_of_fei_inputs two_ne k 𝔭 hX D fei G ?_ hG
  rintro g ⟨f, hf, hd, rfl⟩ notchain _
  exact slice f hf hd notchain

/-- Unit-ideal consequence for the source system: no `e ≤ D` is a last fall degree bound. -/
theorem slice_lfd_gt {σ : Type*} (two_ne : (2 : K) ≠ 0) (k : ℕ) (𝔭 : Ideal (Q K))
    [𝔭.IsPrime] (hX : ∀ i : Fin 2, (X i : Q K) ∉ 𝔭) (D : ℕ) (fei : FEI 𝔭 k (D - 1))
    (φ : MvPolynomial σ K →ₐ[K] P (2 * k) K) (linear : ∀ i, (φ (X i)).totalDegree ≤ 1)
    (F : Set (MvPolynomial σ K))
    (slice : ∀ f ∈ F, f.totalDegree ≤ D → φ f ∉ chainInputs K 2 k →
      ChainClosureBridge.normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k) (φ f) ∈
        idealSpace 𝔭 k D)
    (unit : (1 : MvPolynomial σ K) ∈ Ideal.span F) :
    ¬ HasLFDBound F D := fun bound =>
  slice_not_one two_ne k 𝔭 hX D fei φ linear F slice (one_of_hasLFDBound unit bound)

#print axioms not_one_of_fei_inputs
#print axioms slice_not_one
#print axioms slice_lfd_gt

end ICFallDegrees.GrowthSlice
