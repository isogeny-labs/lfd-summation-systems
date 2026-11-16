import Mathlib.AlgebraicGeometry.EllipticCurve.Group
import Mathlib.RingTheory.Valuation.ExtendToLocalization
import Mathlib.RingTheory.Norm.Basic
import Mathlib.Tactic.Abel
import Mathlib.Tactic.GCongr
import Mathlib.Tactic.NormNum

noncomputable section

/-!
# An actual infinity valuation on the elliptic function field

The norm degree on the coordinate ring is a weighted degree with weights
2 and 3. We prove its addition inequality, package a valuation, and extend
it to the fraction field. No genus-one height or circuit-spread bound is
assumed here.
-/

namespace ICFallDegrees.GrowthNonTorsionInfinity

open Polynomial WeierstrassCurve.Affine
open scoped nonZeroDivisors Multiplicative

variable {K : Type*} [Field K] (W : WeierstrassCurve.Affine K)

private theorem twice_degree_add_le (p q : Polynomial K) :
    2 • (p + q).degree ≤ max (2 • p.degree) (2 • q.degree) := by
  rcases le_total p.degree q.degree with h | h
  · have hd := (Polynomial.degree_add_le p q).trans (max_eq_right h).le
    exact (show 2 • (p+q).degree ≤ 2 • q.degree by gcongr).trans (le_max_right _ _)
  · have hd := (Polynomial.degree_add_le p q).trans (max_eq_left h).le
    exact (show 2 • (p+q).degree ≤ 2 • p.degree by gcongr).trans (le_max_left _ _)

theorem norm_degree_add_le (u v : W.CoordinateRing) :
    (Algebra.norm K[X] (u + v)).degree ≤
      max (Algebra.norm K[X] u).degree (Algebra.norm K[X] v).degree := by
  obtain ⟨p,q,rfl⟩ := CoordinateRing.exists_smul_basis_eq u
  obtain ⟨r,s,rfl⟩ := CoordinateRing.exists_smul_basis_eq v
  have he : (p • (1:W.CoordinateRing) + q • CoordinateRing.mk W Polynomial.X) +
      (r • (1:W.CoordinateRing) + s • CoordinateRing.mk W Polynomial.X) =
      (p+r) • (1:W.CoordinateRing) + (q+s) • CoordinateRing.mk W Polynomial.X := by
    simp only [add_smul]
    abel
  rw [he, CoordinateRing.degree_norm_smul_basis, CoordinateRing.degree_norm_smul_basis,
    CoordinateRing.degree_norm_smul_basis]
  apply max_le
  · apply (twice_degree_add_le p r).trans
    exact max_le ((le_max_left _ _).trans (le_max_left _ _))
      ((le_max_left _ _).trans (le_max_right _ _))
  · have h := add_le_add_right (twice_degree_add_le q s) (3:WithBot ℕ)
    rw [← max_add_add_right] at h
    apply h.trans
    exact max_le ((le_max_right _ _).trans (le_max_left _ _))
      ((le_max_right _ _).trans (le_max_right _ _))

def coordinateValue (u : W.CoordinateRing) : ℤₘ₀ := by
  classical
  exact if u = 0 then 0 else ↑(Multiplicative.ofAdd ((Algebra.norm K[X] u).natDegree : ℤ))

private theorem norm_ne_zero {u : W.CoordinateRing} (hu : u ≠ 0) :
    Algebra.norm K[X] u ≠ 0 :=
  (Algebra.norm_ne_zero_iff_of_basis (CoordinateRing.basis W)).2 hu

theorem coordinateValue_mul (u v : W.CoordinateRing) :
    coordinateValue W (u*v) = coordinateValue W u * coordinateValue W v := by
  classical
  by_cases hu : u = 0
  · simp [coordinateValue,hu]
  by_cases hv : v = 0
  · simp [coordinateValue,hv]
  rw [coordinateValue, coordinateValue, coordinateValue,
    if_neg hu, if_neg hv, if_neg (mul_ne_zero hu hv),
    ← WithZero.coe_mul, WithZero.coe_inj, ← ofAdd_add,
    map_mul, Polynomial.natDegree_mul (norm_ne_zero W hu) (norm_ne_zero W hv)]
  simp

theorem coordinateValue_add_le (u v : W.CoordinateRing) :
    coordinateValue W (u+v) ≤ max (coordinateValue W u) (coordinateValue W v) := by
  classical
  by_cases hu : u = 0
  · simp [hu,coordinateValue]
  by_cases hv : v = 0
  · simp [hv,coordinateValue]
  by_cases huv : u+v = 0
  · simp [huv,coordinateValue]
  have hd := norm_degree_add_le W u v
  rw [Polynomial.degree_eq_natDegree (norm_ne_zero W huv),
    Polynomial.degree_eq_natDegree (norm_ne_zero W hu),
    Polynomial.degree_eq_natDegree (norm_ne_zero W hv)] at hd
  have hn : (Algebra.norm K[X] (u+v)).natDegree ≤
      max (Algebra.norm K[X] u).natDegree (Algebra.norm K[X] v).natDegree := by
    exact WithBot.coe_le_coe.mp (by simpa only [WithBot.coe_max] using hd)
  rw [coordinateValue,coordinateValue,coordinateValue,if_neg hu,if_neg hv,if_neg huv,
    le_max_iff,WithZero.coe_le_coe,Multiplicative.ofAdd_le,
    WithZero.coe_le_coe,Multiplicative.ofAdd_le,← le_max_iff]
  exact_mod_cast hn

def coordinateValuation : Valuation W.CoordinateRing ℤₘ₀ where
  toFun := coordinateValue W
  map_zero' := by simp [coordinateValue]
  map_one' := by simp [coordinateValue]
  map_mul' := coordinateValue_mul W
  map_add_le_max' := coordinateValue_add_le W

theorem coordinateValue_ne_zero {u : W.CoordinateRing} (hu : u ≠ 0) :
    coordinateValue W u ≠ 0 := by
  classical
  simp [coordinateValue,hu]

def infinityValuation : Valuation W.FunctionField ℤₘ₀ :=
  (coordinateValuation W).extendToLocalization
    (S := W.CoordinateRing⁰) (by
      intro u hu
      change coordinateValue W u ≠ 0
      exact coordinateValue_ne_zero W (mem_nonZeroDivisors_iff_ne_zero.mp hu))
    W.FunctionField

theorem infinityValuation_map (u : W.CoordinateRing) :
    infinityValuation W (algebraMap W.CoordinateRing W.FunctionField u) =
      coordinateValue W u :=
  Valuation.extendToLocalization_apply_map_apply _ _ _ _

def coordinateX : W.CoordinateRing := algebraMap K[X] W.CoordinateRing X

def coordinateY : W.CoordinateRing := CoordinateRing.mk W Polynomial.X

def fieldX : W.FunctionField := algebraMap W.CoordinateRing W.FunctionField (coordinateX W)

def fieldY : W.FunctionField := algebraMap W.CoordinateRing W.FunctionField (coordinateY W)

def constants : K →+* W.FunctionField :=
  (algebraMap W.CoordinateRing W.FunctionField).comp
    ((algebraMap K[X] W.CoordinateRing).comp Polynomial.C)

theorem norm_coordinateX : Algebra.norm K[X] (coordinateX W) = (X:K[X]) ^ 2 := by
  simpa [coordinateX] using
    Algebra.norm_algebraMap_of_basis (CoordinateRing.basis W) (X:K[X])

theorem norm_coordinateY : Algebra.norm K[X] (coordinateY W) =
    -(X ^ 3 + C W.a₂ * X ^ 2 + C W.a₄ * X + C W.a₆) := by
  simpa [coordinateY] using CoordinateRing.norm_smul_basis W (0:K[X]) 1

theorem norm_coordinateY_natDegree : (Algebra.norm K[X] (coordinateY W)).natDegree = 3 := by
  rw [norm_coordinateY,Polynomial.natDegree_neg]
  simpa using (Polynomial.natDegree_cubic (a := (1:K)) (b := W.a₂)
    (c := W.a₄) (d := W.a₆) one_ne_zero)

theorem infinityValuation_X : infinityValuation W (fieldX W) =
    ↑(Multiplicative.ofAdd (2:ℤ)) := by
  have hx : coordinateX W ≠ 0 := by
    apply (Algebra.norm_ne_zero_iff_of_basis (CoordinateRing.basis W)).1
    rw [norm_coordinateX]
    exact pow_ne_zero 2 Polynomial.X_ne_zero
  rw [fieldX,infinityValuation_map,coordinateValue,if_neg hx,norm_coordinateX,
    Polynomial.natDegree_X_pow]
  norm_num

theorem infinityValuation_Y : infinityValuation W (fieldY W) =
    ↑(Multiplicative.ofAdd (3:ℤ)) := by
  have hy : coordinateY W ≠ 0 := by
    apply (Algebra.norm_ne_zero_iff_of_basis (CoordinateRing.basis W)).1
    intro hz
    have hd := norm_coordinateY_natDegree W
    rw [hz,Polynomial.natDegree_zero] at hd
    omega
  rw [fieldY,infinityValuation_map,coordinateValue,if_neg hy,norm_coordinateY_natDegree]
  norm_num

theorem infinityValuation_constant (c : K) (hc : c ≠ 0) :
    infinityValuation W (constants W c) = 1 := by
  have hn : Algebra.norm K[X] (algebraMap K[X] W.CoordinateRing (C c)) = (C c)^2 := by
    simpa using Algebra.norm_algebraMap_of_basis (CoordinateRing.basis W) (C c)
  have hz : algebraMap K[X] W.CoordinateRing (C c) ≠ 0 := by
    apply (Algebra.norm_ne_zero_iff_of_basis (CoordinateRing.basis W)).1
    rw [hn]
    exact pow_ne_zero 2 (Polynomial.C_ne_zero.mpr hc)
  change infinityValuation W (algebraMap W.CoordinateRing W.FunctionField
    (algebraMap K[X] W.CoordinateRing (C c))) = 1
  rw [infinityValuation_map,coordinateValue,if_neg hz,hn,Polynomial.natDegree_pow,
    Polynomial.natDegree_C]
  rfl

theorem infinityValuation_constant_le_one (c : K) :
    infinityValuation W (constants W c) ≤ 1 := by
  by_cases hc : c = 0
  · simp [hc]
  · rw [infinityValuation_constant W c hc]

theorem infinityValuation_uniformizer :
    infinityValuation W (fieldX W / fieldY W) = ↑(Multiplicative.ofAdd (-1:ℤ)) := by
  rw [Valuation.map_div,infinityValuation_X,infinityValuation_Y,
    ← WithZero.coe_div,← ofAdd_sub]
  norm_num



end ICFallDegrees.GrowthNonTorsionInfinity
