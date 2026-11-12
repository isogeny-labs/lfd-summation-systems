import ICFallDegrees.PrimeM2Truncated
import ICFallDegrees.PrimeM2Sharp
import Mathlib.Data.ZMod.Basic
import Mathlib.Algebra.Field.ZMod

/-! Prime-field entry points and explicit admissible finite instances. -/
noncomputable section
open Classical
namespace ICFallDegrees.PrimeFieldConclusions
open MvPolynomial HomogeneousFirstFall PrimeM2General

theorem two_ne_zero (p : ℕ) (large : 5≤p) : (2 : ZMod p)≠0 := by
  intro zero
  have divides := (ZMod.natCast_zmod_eq_zero_iff_dvd 2 p).mp zero
  have bound := Nat.le_of_dvd (by decide : 0<2) divides
  omega

theorem ordinary_ffd (p : ℕ) [Fact p.Prime] (large : 5≤p) (k : ℕ) (depth : 2≤k)
    (a b t : ZMod p) :
    IsKoszulFirstFall (inputs k depth a b t) (if 3*t^2+a=0 then 3 else 4) := by
  have result := exact_ordinary_ffd (two_ne_zero p large) k depth a b t
  by_cases condition : 3*t^2+a=0
  · simpa only [if_pos condition] using result
  · simpa only [if_neg condition] using result

theorem truncated_ffd (p : ℕ) [Fact p.Prime] (large : 5≤p) (k : ℕ) (depth : 2≤k)
    (a b t : ZMod p) :
    TruncatedFirstFall.IsFirstFall p (inputs k depth a b t) (if 3*t^2+a=0 then 3 else 4) := by
  have result := PrimeM2Truncated.exact_truncated_ffd p large (two_ne_zero p large) k depth a b t
  by_cases condition : 3*t^2+a=0
  · simpa only [if_pos condition] using result
  · simpa only [if_neg condition] using result

/-- The target is on the nonsingular short Weierstrass curve for every
nonzero a in an admissible prime field. -/
theorem sharp_curve_data (p : ℕ) [Fact p.Prime] (large : 5≤p) (a : ZMod p) (nonzero : a≠0) :
    (0 : ZMod p)^2=0^3+a*0+0 ∧ 4*a^3+27*0^2≠0 := by
  refine ⟨by simp,?_⟩
  have four_ne : (4 : ZMod p)≠0 := by
    have identity : (4 : ZMod p)=2*2 := by ring
    rw [identity]
    exact mul_ne_zero (two_ne_zero p large) (two_ne_zero p large)
  simpa only [zero_pow (by decide : 2≠0),mul_zero,add_zero] using mul_ne_zero four_ne (pow_ne_zero 3 nonzero)

/-- Exact LFD on a genuine subgroup-chain instance. Geometric nonsingularity
and target membership follow from sharp_curve_data; subgroup size is explicit. -/
theorem sharp_lfd (p : ℕ) [Fact p.Prime] (large : 5≤p) (k : ℕ) (depth : 2≤k)
    (_subgroup : 2^k∣p-1) (a : ZMod p) (nonzero : a≠0) (outside : a^(2^k)≠1) :
    IsLastFallDegree (Set.range (inputs k depth a 0 0)) 4 :=
  PrimeM2Sharp.exact_lfd_four (two_ne_zero p large) k depth a nonzero outside

instance : Fact (Nat.Prime 17) := ⟨by decide⟩
instance : Fact (Nat.Prime 97) := ⟨by decide⟩

theorem prime17_depth2_lfd : IsLastFallDegree (Set.range (inputs 2 (by decide) (2 : ZMod 17) 0 0)) 4 :=
  sharp_lfd 17 (by decide) 2 (by decide) (by decide) 2 (by decide) (by decide)

theorem prime17_depth3_lfd : IsLastFallDegree (Set.range (inputs 3 (by decide) (3 : ZMod 17) 0 0)) 4 :=
  sharp_lfd 17 (by decide) 3 (by decide) (by decide) 3 (by decide) (by decide)

theorem prime97_depth4_lfd : IsLastFallDegree (Set.range (inputs 4 (by decide) (2 : ZMod 97) 0 0)) 4 :=
  sharp_lfd 97 (by decide) 4 (by decide) (by decide) 2 (by decide) (by decide)

end ICFallDegrees.PrimeFieldConclusions
