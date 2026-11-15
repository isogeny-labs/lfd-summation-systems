import ICFallDegrees.GrowthCircuitReplacementFin
import ICFallDegrees.GrowthCircuitErasedBasis
import ICFallDegrees.GrowthActualReplacementSpread
import ICFallDegrees.GrowthCircuitInterface
import Mathlib.Tactic

/-!
# N8: genuine pairwise exponent control on arbitrary minimal circuits

Every pair of distinct indices in a GrowthSpread.Circuit for the actual
rational section comes from a nonzero coefficient in an independent
reindexed replacement relation. The sharp 4B cofactor theorem therefore
applies. This theorem does not postulate a replacement relation, minor,
pole identity, or SpreadBound.
-/

noncomputable section
namespace ICFallDegrees.GrowthCircuitPairBound

open ICFallDegrees.GrowthSpread
open ICFallDegrees.GrowthCircuitInterface
open ICFallDegrees.GrowthCircuitErasedBasis
open ICFallDegrees.GrowthCircuitReplacementFin
open ICFallDegrees.GrowthActualReplacementSpread
open ICFallDegrees.GrowthRatFuncModel
open ICFallDegrees.GrowthRatFuncWronskianNormalization

variable {K : Type*} [Field K] [CharZero K]

theorem distinct_circuit_pair_le_four_budget
    (a : K) (ha : a + 3 ≠ 0)
    (J : Finset (Fin 2 →₀ ℕ))
    (hcirc : Circuit K (mono (xi a)) J)
    (x y : Fin 2 →₀ ℕ)
    (hx : x ∈ J) (hy : y ∈ J) (hxy : y ≠ x) :
    |2 * ((x 0 : ℤ) - y 0)| +
      |2 * ((x 1 : ℤ) - y 1)| ≤
        4 * ((∑ j : Fin (J.erase x).card, 2 * j.val : ℕ) : ℤ) := by
  classical
  let er : Fin (J.erase x).card → (Fin 2 →₀ ℕ) :=
    fun i => (((Finset.equivFin (J.erase x)).symm i).1)
  let e : Fin (J.erase x).card → ℕ × ℕ :=
    fun i => (er i 0, er i 1)
  let enew : ℕ × ℕ := (x 0, x 1)
  have hind0 :=
    erase_reindex_linearIndependent (mono (xi a)) J hcirc x hx
  have hind : LinearIndependent K
      (fun j : Fin (J.erase x).card => monomial (a + 2) (e j)) := by
    simpa only [e, er, mono_xi_eq_monomial] using hind0
  obtain ⟨coeff, hcoeff, hrel⟩ :=
    exists_replacement_fin (mono (xi a)) J hcirc x hx
  have hrel' : monomial (a + 2) enew =
      ∑ j : Fin (J.erase x).card, coeff j • monomial (a + 2) (e j) := by
    simpa only [enew, e, er, mono_xi_eq_monomial] using hrel
  obtain ⟨i, hi⟩ := erase_reindex_contains J x y hx hy hxy
  have ha' : (a + 2) + 1 ≠ 0 := by
    convert ha using 1 <;> ring
  have hbound := replacement_l1_le_four_cofactor_budget
    (J.erase x).card (a + 2) ha' e enew coeff i hind
    (hcoeff i) hrel'
  simpa only [enew, e, er, hi] using hbound

end ICFallDegrees.GrowthCircuitPairBound
