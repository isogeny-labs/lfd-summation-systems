import Mathlib.Algebra.Module.LinearMap.Basic

/-!
Elementary identities for a proposed proper-ideal growth argument.

These lemmas do NOT establish a growth theorem. The chain Euler-derivative
commutator and the three vanishing hypotheses must still be instantiated
from the concrete growth splitting theorem. The first lemma records the
cofactor identity supplied by Frobenius in a finite chain quotient.
-/
namespace ICFallDegrees.GrowthJetPrototype

/-- If Frobenius fixes `h`, the square `h^2` generates `h` as an ideal element. -/
theorem frobenius_cofactor {R : Type*} [Monoid R] (h : R) (p : ℕ)
    (hp : 2 ≤ p) (hf : h ^ p = h) : h = h ^ 2 * h ^ (p - 2) := by
  calc
    h = h ^ p := hf.symm
    _ = h ^ (2 + (p - 2)) := by
      rw [Nat.add_comm 2 (p - 2), Nat.sub_add_cancel hp]
    _ = h ^ 2 * h ^ (p - 2) := pow_add h 2 (p - 2)

/-- A commutator reduces first-derivative invariance to three scalar vanishings.
No analytic splitting property or chain specialization is assumed proved here. -/
theorem shift_jet_zero {K V W : Type*} [Ring K]
    [AddCommGroup V] [Module K V] [AddCommGroup W] [Module K W]
    (φ : V →ₗ[K] W) (δ T P : V →ₗ[K] V) (f : V) (c n : K)
    (hcomm : δ (T f) = T (δ f) + c • T f - n • T (P f))
    (hplain : φ (T f) = 0) (hderiv : φ (T (δ f)) = 0)
    (hwrap : φ (T (P f)) = 0) :
    φ (T f) = 0 ∧ φ (δ (T f)) = 0 := by
  constructor
  · exact hplain
  · rw [hcomm]
    simp only [map_sub, map_add, map_smul, hplain, hderiv, hwrap,
      smul_zero, zero_add, sub_zero]

#print axioms frobenius_cofactor
#print axioms shift_jet_zero
end ICFallDegrees.GrowthJetPrototype
