import ICFallDegrees.Closure
import Mathlib.Data.Finsupp.Weight
import Mathlib.Algebra.MvPolynomial.Eval
import Mathlib.Algebra.BigOperators.Ring.Finset

/-! Degree-safe transport under affine polynomial substitutions. Input
membership is checked at the original budget; ideal equality alone is not used. -/
noncomputable section
namespace ICFallDegrees.AffineTransport
open MvPolynomial
open scoped BigOperators
variable {σ τ K : Type*} [CommRing K]

private theorem monomial_degree (f : MvPolynomial σ K →ₐ[K] MvPolynomial τ K)
    (linear : ∀ i, (f (X i)).totalDegree≤1) (e : σ →₀ ℕ) (c : K) :
    (f (monomial e c)).totalDegree≤e.degree := by
  have constants : f (C c)=C c := f.commutes c
  rw [monomial_eq,map_mul,constants]
  have product : (f (e.prod fun i j => (X i : MvPolynomial σ K)^j)).totalDegree≤e.degree := by
    simp only [Finsupp.prod,map_prod,map_pow]
    apply (totalDegree_finset_prod _ _).trans
    apply Finset.sum_le_sum
    intro i _
    exact (totalDegree_pow _ _).trans (by simpa using Nat.mul_le_mul_left (e i) (linear i))
  exact (totalDegree_mul _ _).trans (by simpa only [totalDegree_C,zero_add] using product)

theorem degree_le (f : MvPolynomial σ K →ₐ[K] MvPolynomial τ K)
    (linear : ∀ i, (f (X i)).totalDegree≤1) (p : MvPolynomial σ K) :
    (f p).totalDegree≤p.totalDegree := by
  classical
  have identity : f p=∑ e∈p.support, f (monomial e (coeff e p)) := by
    conv_lhs => rw [p.as_sum]
    rw [map_sum]
  rw [identity]
  apply totalDegree_finsetSum_le
  intro e he
  exact (monomial_degree f linear e (coeff e p)).trans (le_totalDegree he)

theorem map_closure (f : MvPolynomial σ K →ₐ[K] MvPolynomial τ K)
    (linear : ∀ i, (f (X i)).totalDegree≤1)
    (F : Set (MvPolynomial σ K)) (G : Set (MvPolynomial τ K)) (D : ℕ)
    (inputs : ∀ p∈F, p.totalDegree≤D → InClosure G D (f p))
    {p} (member : InClosure F D p) : InClosure G D (f p) := by
  induction member with
  | input belongs degree => exact inputs _ belongs degree
  | zero => simpa using (InClosure.zero : InClosure G D 0)
  | add _ _ ih ih' => simpa only [map_add] using ih.add ih'
  | smul c _ ih => simpa only [map_smul] using ih.smul c
  | mul a hp budget ih =>
      simpa only [map_mul] using ih.mul (f a)
        ((Nat.add_le_add (degree_le f linear a) (degree_le f linear _)).trans budget)

end ICFallDegrees.AffineTransport
