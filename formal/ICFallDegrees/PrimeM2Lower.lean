import ICFallDegrees.PrimeM2Sharp
import ICFallDegrees.PrimeFieldConclusions

/-! An all-depth lower bound for every curve and target.

For the normalized two-summand chain presentation `PrimeM2General.inputs k depth a b t` at every
depth `k ≥ 2`, over every field of characteristic other than two, the constant `1` is not in the
degree-3 closure as soon as `3*t^2 + a ≠ 0`.  Hence every such system whose ideal is the unit
ideal (no scalar solution in the subgroup grid) has `LFD ≥ 4`.

This generalizes `PrimeM2Sharp.no_one_at_three`, which treats the single family `b = t = 0`
(where the condition reads `a ≠ 0`), to all curves `y^2 = x^3 + a x + b` and all target
abscissae `t`.  The proof instantiates the quadratic-top theorem `QuadraticNoLow.no_one_at_three`
with `α = -2(t^2+a)` and `β = -2t`; its degree-3 syzygy hypothesis is the nondegeneracy
`α - β^2 = -2(3t^2 + a) ≠ 0`, the same condition that separates ordinary FFD 3 from 4 in
`PrimeM2General.exact_ordinary_ffd`. -/
noncomputable section
open Classical
namespace ICFallDegrees.PrimeM2Lower
open MvPolynomial HomogeneousFirstFall GeneralQuadraticTop GeneralQuadraticClassification
open PrimeM2General
variable {K : Type*} [Field K]
abbrev P (k : ℕ) (K : Type*) [Field K] := PrimeM2General.P k K

theorem nondegenerate_identity (a t : K) : -2*(t^2+a)-(-2*t)^2=(-2)*(3*t^2+a) := by ring

/-- `1` is not derivable at degree 3, for every depth, curve and target with `3t^2 + a ≠ 0`. -/
theorem no_one_at_three (two_ne : (2 : K)≠0) (k : ℕ) (depth : 2≤k) (a b t : K)
    (nondegenerate : 3*t^2+a≠0) :
    ¬InClosure (Set.range (inputs k depth a b t)) 3 1 := by
  let tails : Fin (2*k+1) → P k K := Fin.lastCases (tail k depth a b t) (chainTail k)
  have presentation : inputs k depth a b t=
      (fun i => q (active k depth) (-2*(t^2+a)) (-2*t) i+tails i) := by
    funext i
    refine Fin.lastCases ?_ (fun i => ?_) i
    · simp [inputs,normalized,tails,q]
    · simp [inputs,q,tails]
  rw [presentation]
  apply QuadraticNoLow.no_one_at_three
  · exact homogeneous_q _ _ _
  · intro i
    refine Fin.lastCases ?_ (fun i => ?_) i
    · have bound := tail_small k depth a b t
      simpa only [tails,Fin.lastCases_last] using (show (tail k depth a b t).totalDegree≤1 by omega)
    · have bound := chain_tail_small (K:=K) k i
      simpa only [tails,Fin.lastCases_castSucc] using (show (chainTail k i : P k K).totalDegree≤1 by omega)
  · intro h relation
    exact low_syzygy_zero _ (pivot (active k depth))
      (coefficient_pivots _ (active_injective k depth) (-2*(t^2+a)) (-2*t)) 2 2 le_rfl relation
  · intro h relation
    apply degree_three_syzygy_zero _ (active_injective k depth) (-2*(t^2+a)) (-2*t) _ relation
    rw [nondegenerate_identity]
    exact mul_ne_zero (neg_ne_zero.mpr two_ne) nondegenerate

/-- Every unit-ideal two-summand chain system with `3t^2 + a ≠ 0` has `LFD ≥ 4`, at every depth:
no `e < 4` is an LFD bound. -/
theorem lfd_ge_four_of_unit (two_ne : (2 : K)≠0) (k : ℕ) (depth : 2≤k) (a b t : K)
    (nondegenerate : 3*t^2+a≠0)
    (unit : (1 : P k K)∈Ideal.span (Set.range (inputs k depth a b t))) :
    ∀ e<4, ¬HasLFDBound (Set.range (inputs k depth a b t)) e := by
  intro e small bound
  have derived := one_of_hasLFDBound unit bound
  exact no_one_at_three two_ne k depth a b t nondegenerate (derived.mono (by omega))

/-- Prime-field entry point: every prime `p ≥ 5`, depth `k ≥ 2`, curve and target. -/
theorem prime_field_lfd_ge_four_of_unit (p : ℕ) [Fact p.Prime] (large : 5≤p) (k : ℕ) (depth : 2≤k)
    (a b t : ZMod p) (nondegenerate : 3*t^2+a≠0)
    (unit : (1 : P k (ZMod p))∈Ideal.span (Set.range (inputs k depth a b t))) :
    ∀ e<4, ¬HasLFDBound (Set.range (inputs k depth a b t)) e :=
  lfd_ge_four_of_unit (PrimeFieldConclusions.two_ne_zero p large) k depth a b t nondegenerate unit

end ICFallDegrees.PrimeM2Lower
