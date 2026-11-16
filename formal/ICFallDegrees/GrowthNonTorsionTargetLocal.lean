import ICFallDegrees.GrowthNonTorsionInfinitySecant
import Mathlib.RingTheory.DiscreteValuationRing.TFAE
import Mathlib.RingTheory.Localization.LocalizationLocalization
import Mathlib.RingTheory.DedekindDomain.AdicValuation

noncomputable section

/-!
# The actual local ring at a non-two-torsion point

Evaluation on the coordinate ring identifies the point ideal. Localizing
Mathlib's product of the point and opposite-point ideals makes the latter
the unit ideal, so x-t generates the local maximal ideal. Noetherianity,
nonzeroness of x-t and the local PID criterion then give an actual DVR.
-/

namespace ICFallDegrees.GrowthNonTorsionTargetLocal

open Polynomial WeierstrassCurve.Affine
open GrowthNonTorsionInfinity GrowthNonTorsionInfinitySecant

variable {K : Type*} [Field K] (W : WeierstrassCurve.Affine K)

def pointEval (t Y : K) (hR : W.Equation t Y) : W.CoordinateRing →+* K :=
  AdjoinRoot.lift (Polynomial.evalRingHom t) Y (by
    simpa only [Polynomial.eval₂_evalRingHom] using hR)

@[simp] theorem pointEval_constant (t Y c : K) (hR : W.Equation t Y) :
    pointEval W t Y hR (algebraMap K W.CoordinateRing c) = c := by
  change AdjoinRoot.lift _ _ _ (AdjoinRoot.of _ (Polynomial.C c)) = c
  simp [AdjoinRoot.lift_of]

@[simp] theorem pointEval_X (t Y : K) (hR : W.Equation t Y) :
    pointEval W t Y hR (coordinateX W) = t := by
  change AdjoinRoot.lift _ _ _ (AdjoinRoot.of _ Polynomial.X) = t
  simp [AdjoinRoot.lift_of]

@[simp] theorem pointEval_Y (t Y : K) (hR : W.Equation t Y) :
    pointEval W t Y hR (coordinateY W) = Y := by
  exact AdjoinRoot.lift_root _

theorem pointEval_surjective (t Y : K) (hR : W.Equation t Y) :
    Function.Surjective (pointEval W t Y hR) := by
  intro c
  exact ⟨algebraMap K W.CoordinateRing c, pointEval_constant W t Y c hR⟩

@[simp] theorem pointEval_YClass (t Y : K) (hR : W.Equation t Y) (p : K[X]) :
    pointEval W t Y hR (CoordinateRing.YClass W p) = Y-p.eval t := by
  simp [CoordinateRing.YClass,pointEval,AdjoinRoot.lift_mk]

@[simp] theorem pointEval_XClass (t Y c : K) (hR : W.Equation t Y) :
    pointEval W t Y hR (CoordinateRing.XClass W c) = t-c := by
  simp [CoordinateRing.XClass,pointEval,AdjoinRoot.lift_mk]

def pointIdeal (t Y : K) : Ideal W.CoordinateRing :=
  CoordinateRing.XYIdeal W t (Polynomial.C Y)

theorem pointIdeal_maximal (t Y : K) (hR : W.Equation t Y) :
    (pointIdeal W t Y).IsMaximal := by
  exact Ideal.Quotient.maximal_of_isField _
    ((CoordinateRing.quotientXYIdealEquiv W hR).toRingEquiv.isField K (Field.toIsField K))

theorem pointIdeal_eq_ker (t Y : K) (hR : W.Equation t Y) :
    pointIdeal W t Y = RingHom.ker (pointEval W t Y hR) := by
  apply (pointIdeal_maximal W t Y hR).eq_of_le
    (RingHom.ker_isMaximal_of_surjective _ (pointEval_surjective W t Y hR)).ne_top
  rw [pointIdeal,CoordinateRing.XYIdeal,Ideal.span_le]
  intro u hu
  rcases hu with rfl | rfl
  · change pointEval W t Y hR (CoordinateRing.XClass W t) = 0
    change pointEval W t Y hR (AdjoinRoot.mk _ (Polynomial.C (X-C t))) = 0
    simp [pointEval,AdjoinRoot.lift_mk]
  · change pointEval W t Y hR (CoordinateRing.YClass W (C Y)) = 0
    change pointEval W t Y hR (AdjoinRoot.mk _ (X-Polynomial.C (C Y))) = 0
    simp [pointEval,AdjoinRoot.lift_mk]

instance coordinateRing_noetherian : IsNoetherianRing W.CoordinateRing := by
  letI := Module.Finite.of_basis (CoordinateRing.basis W)
  exact isNoetherianRing_iff.mpr (isNoetherian_of_tower K[X] (inferInstance))

def targetIdeal (t Y : K) (hR : W.Equation t Y) : Ideal W.CoordinateRing :=
  RingHom.ker (pointEval W t Y hR)

instance targetIdeal_prime (t Y : K) (hR : W.Equation t Y) :
    (targetIdeal W t Y hR).IsPrime := RingHom.ker_isPrime _

abbrev PointLocalRing (t Y : K) (hR : W.Equation t Y) :=
  Localization.AtPrime (targetIdeal W t Y hR)

instance pointLocal_noetherian (t Y : K) (hR : W.Equation t Y) :
    IsNoetherianRing (PointLocalRing W t Y hR) :=
  IsLocalization.isNoetherianRing (targetIdeal W t Y hR).primeCompl _ (inferInstance)

theorem pointLocal_maximal_principal (t Y : K) (hR : W.Equation t Y)
    (hns : W.Nonsingular t Y) (hY : Y ≠ W.negY t Y) :
    IsLocalRing.maximalIdeal (PointLocalRing W t Y hR) =
      Ideal.span {algebraMap W.CoordinateRing (PointLocalRing W t Y hR)
        (CoordinateRing.XClass W t)} := by
  let f := algebraMap W.CoordinateRing (PointLocalRing W t Y hR)
  have hunit : Ideal.map f (CoordinateRing.XYIdeal W t (C (W.negY t Y))) = ⊤ := by
    apply Ideal.eq_top_of_isUnit_mem _
      (Ideal.mem_map_of_mem f (Ideal.subset_span (Set.mem_insert_of_mem _ (Set.mem_singleton _))))
    apply (IsLocalization.AtPrime.isUnit_to_map_iff _ (targetIdeal W t Y hR) _).2
    change pointEval W t Y hR (CoordinateRing.YClass W (C (W.negY t Y))) ≠ 0
    simpa using sub_ne_zero.mpr hY
  have hm := congrArg (Ideal.map f) (CoordinateRing.XYIdeal_neg_mul hns)
  rw [Ideal.map_mul,hunit,Ideal.top_mul] at hm
  calc
    _ = Ideal.map f (targetIdeal W t Y hR) := Localization.AtPrime.map_eq_maximalIdeal.symm
    _ = Ideal.map f (pointIdeal W t Y) :=
      congrArg (Ideal.map f) (pointIdeal_eq_ker W t Y hR).symm
    _ = Ideal.map f (CoordinateRing.XIdeal W t) := hm
    _ = _ := by simp only [CoordinateRing.XIdeal,Ideal.map_span,Set.image_singleton,f]

def localParameter (t Y : K) (hR : W.Equation t Y) : PointLocalRing W t Y hR :=
  algebraMap W.CoordinateRing (PointLocalRing W t Y hR) (CoordinateRing.XClass W t)

theorem localParameter_ne_zero (t Y : K) (hR : W.Equation t Y) :
    localParameter W t Y hR ≠ 0 := by
  simpa only [localParameter,map_zero] using
    (IsLocalization.injective (PointLocalRing W t Y hR)
      (targetIdeal W t Y hR).primeCompl_le_nonZeroDivisors).ne
        (CoordinateRing.XClass_ne_zero W t)

theorem pointLocal_not_field (t Y : K) (hR : W.Equation t Y)
    (hns : W.Nonsingular t Y) (hY : Y ≠ W.negY t Y) :
    ¬ IsField (PointLocalRing W t Y hR) := by
  intro hf
  have hm := IsLocalRing.isField_iff_maximalIdeal_eq.mp hf
  rw [pointLocal_maximal_principal W t Y hR hns hY] at hm
  exact localParameter_ne_zero W t Y hR (Ideal.span_singleton_eq_bot.mp hm)

theorem pointLocal_dvr (t Y : K) (hR : W.Equation t Y)
    (hns : W.Nonsingular t Y) (hY : Y ≠ W.negY t Y) :
    IsDiscreteValuationRing (PointLocalRing W t Y hR) := by
  apply ((IsDiscreteValuationRing.TFAE _ (pointLocal_not_field W t Y hR hns hY)).out 4 0).mp
  rw [pointLocal_maximal_principal W t Y hR hns hY]
  exact ⟨_,rfl⟩

end ICFallDegrees.GrowthNonTorsionTargetLocal
