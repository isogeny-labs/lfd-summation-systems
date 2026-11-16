import ICFallDegrees.GrowthNonTorsionTargetLocal

noncomputable section

/-!
# The target-point valuation on the actual function field

The constructed DVR's maximal ideal defines a normalized adic valuation,
extended to the fraction field of the original coordinate ring. Its
uniformizer x-t has multiplicative value ofAdd(-1); elements with nonzero
point evaluation have value one.
-/

namespace ICFallDegrees.GrowthNonTorsionTargetValuation

open Polynomial WeierstrassCurve.Affine
open GrowthNonTorsionInfinity GrowthNonTorsionInfinitySecant GrowthNonTorsionTargetLocal
open scoped Multiplicative

private def maximalPlace (B : Type*) [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] : IsDedekindDomain.HeightOneSpectrum B :=
  ⟨IsLocalRing.maximalIdeal B, inferInstance, IsDiscreteValuationRing.not_a_field B⟩

variable {K : Type*} [Field K] (W : WeierstrassCurve.Affine K)

def targetValuation (t Y : K) (hR : W.Equation t Y)
    (hns : W.Nonsingular t Y) (hY : Y ≠ W.negY t Y) : Valuation W.FunctionField ℤₘ₀ := by
  letI := pointLocal_dvr W t Y hR hns hY
  exact (maximalPlace (PointLocalRing W t Y hR)).valuation W.FunctionField

theorem targetValuation_parameter (t Y : K) (hR : W.Equation t Y)
    (hns : W.Nonsingular t Y) (hY : Y ≠ W.negY t Y) :
    targetValuation W t Y hR hns hY (fieldX W - constants W t) =
      ↑(Multiplicative.ofAdd (-1:ℤ)) := by
  letI := pointLocal_dvr W t Y hR hns hY
  let p := maximalPlace (PointLocalRing W t Y hR)
  have he : fieldX W - constants W t =
      algebraMap (PointLocalRing W t Y hR) W.FunctionField (localParameter W t Y hR) := by
    rw [localParameter,← IsScalarTower.algebraMap_apply]
    simp only [CoordinateRing.XClass,map_sub,AdjoinRoot.mk_C,fieldX,constants,
      RingHom.comp_apply,coordinateX,Polynomial.C_sub]
    rfl
  change p.valuation W.FunctionField _ = _
  rw [he,p.valuation_of_algebraMap]
  exact p.intValuation_singleton (localParameter_ne_zero W t Y hR)
    (pointLocal_maximal_principal W t Y hR hns hY)

theorem targetValuation_unit (t Y : K) (hR : W.Equation t Y)
    (hns : W.Nonsingular t Y) (hY : Y ≠ W.negY t Y)
    (u : W.CoordinateRing) (hu : pointEval W t Y hR u ≠ 0) :
    targetValuation W t Y hR hns hY (algebraMap W.CoordinateRing W.FunctionField u) = 1 := by
  letI := pointLocal_dvr W t Y hR hns hY
  let p := maximalPlace (PointLocalRing W t Y hR)
  let f := algebraMap W.CoordinateRing (PointLocalRing W t Y hR)
  have hu' : IsUnit (f u) :=
    (IsLocalization.AtPrime.isUnit_to_map_iff _ (targetIdeal W t Y hR) u).2 hu
  change p.valuation W.FunctionField _ = _
  rw [IsScalarTower.algebraMap_apply W.CoordinateRing (PointLocalRing W t Y hR),
    p.valuation_of_algebraMap]
  apply le_antisymm (p.intValuation_le_one _)
  apply le_of_not_gt
  intro hlt
  have hm : f u ∈ IsLocalRing.maximalIdeal (PointLocalRing W t Y hR) :=
    Ideal.dvd_span_singleton.mp ((p.intValuation_lt_one_iff_dvd _).mp hlt)
  exact (IsLocalRing.not_mem_maximalIdeal.mpr hu') hm

theorem targetValuation_X (t Y : K) (hR : W.Equation t Y)
    (hns : W.Nonsingular t Y) (hY : Y ≠ W.negY t Y) (ht : t ≠ 0) :
    targetValuation W t Y hR hns hY (fieldX W) = 1 :=
  targetValuation_unit W t Y hR hns hY (coordinateX W) (by simpa using ht)

theorem targetValuation_constant (t Y : K) (hR : W.Equation t Y)
    (hns : W.Nonsingular t Y) (hY : Y ≠ W.negY t Y) (c : K) (hc : c ≠ 0) :
    targetValuation W t Y hR hns hY (constants W c) = 1 := by
  simpa only [constants,RingHom.comp_apply,← IsScalarTower.algebraMap_apply] using
    targetValuation_unit W t Y hR hns hY (algebraMap K W.CoordinateRing c) (by simpa using hc)

end ICFallDegrees.GrowthNonTorsionTargetValuation
