import ICFallDegrees.GrowthPrimeE6BoxReplacementGap
import ICFallDegrees.GrowthCircuitReplacementFin
import ICFallDegrees.GrowthCircuitErasedBasis
import ICFallDegrees.GrowthCircuitInterface
import Mathlib.Tactic

/-!
# N8: actual bounded E6 circuit pair separation in positive characteristic

Reindex the independent complement of a minimal E6 monomial circuit,
extract a genuine nonzero-coefficient replacement, and apply the actual
positive-characteristic E6 cofactor/pole bound. This replaces the
char-zero-only circuit pair bound for exponent sets in one finite box.
No circuit-transfer premise from Q is required.
-/

noncomputable section
namespace ICFallDegrees.GrowthPrimeE6BoxCircuitPair

open ICFallDegrees.GrowthSpread
open ICFallDegrees.GrowthRatFuncModel
open ICFallDegrees.GrowthRatFuncWronskianNormalization
open ICFallDegrees.GrowthCircuitInterface
open ICFallDegrees.GrowthCircuitErasedBasis
open ICFallDegrees.GrowthCircuitReplacementFin
open ICFallDegrees.GrowthPrimeE6BoxReplacementGap

variable {K : Type*} [Field K]

/-- In the *actual E6 model* every pair in any minimal monomial
dependence inside the N-by-N box satisfies the sharp cofactor bound
provided the characteristic exceeds the cleared polynomial degrees. -/
theorem e6_bounded_circuit_pair_le_four_budget
    (N : ℕ)
    (J : Finset (Fin 2 →₀ ℕ))
    (hbox : ∀ d ∈ J, d 0 < N ∧ d 1 < N)
    (hcirc : Circuit K (mono (xi (6 : K))) J)
    (x y : Fin 2 →₀ ℕ)
    (hx : x ∈ J) (hy : y ∈ J) (hxy : y ≠ x)
    (hchar : ringChar K = 0 ∨ 2 * N ≤ ringChar K)
    (h9 : (8 : K) + 1 ≠ 0) :
    |2 * ((x 0 : ℤ) - y 0)| +
      |2 * ((x 1 : ℤ) - y 1)| ≤
        4 * ((∑ j : Fin (J.erase x).card, 2 * j.val : ℕ) : ℤ) := by
  classical
  let er : Fin (J.erase x).card → (Fin 2 →₀ ℕ) :=
    fun i => (((Finset.equivFin (J.erase x)).symm i).1)
  let e : Fin (J.erase x).card → ℕ × ℕ :=
    fun i => (er i 0, er i 1)
  let enew : ℕ × ℕ := (x 0, x 1)
  have h8 : (6 : K) + 2 = 8 := by ring
  have hind0 :=
    erase_reindex_linearIndependent (mono (xi (6 : K))) J hcirc x hx
  have hind : LinearIndependent K
      (fun j : Fin (J.erase x).card => monomial (8 : K) (e j)) := by
    simpa only [e, er, mono_xi_eq_monomial, h8] using hind0
  have hboxe : ∀ j : Fin (J.erase x).card, (e j).1 < N ∧ (e j).2 < N := by
    intro j
    have hj : er j ∈ J.erase x :=
      (((Finset.equivFin (J.erase x)).symm j).2)
    exact hbox (er j) (Finset.mem_of_mem_erase hj)
  obtain ⟨coeff, hcoeff, hrel⟩ :=
    exists_replacement_fin (mono (xi (6 : K))) J hcirc x hx
  have hrel' : monomial (8 : K) enew =
      ∑ j : Fin (J.erase x).card, coeff j • monomial (8 : K) (e j) := by
    simpa only [enew, e, er, mono_xi_eq_monomial, h8] using hrel
  obtain ⟨i, hi⟩ := erase_reindex_contains J x y hx hy hxy
  have hbound := e6_box_replacement_l1_le_four_budget
    (J.erase x).card N e enew coeff i hboxe hind
    (hcoeff i) hrel' hchar h9
  simpa only [enew, e, er, hi] using hbound

#print axioms e6_bounded_circuit_pair_le_four_budget

end ICFallDegrees.GrowthPrimeE6BoxCircuitPair
