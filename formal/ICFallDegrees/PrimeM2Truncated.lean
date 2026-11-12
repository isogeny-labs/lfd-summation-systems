import ICFallDegrees.PrimeM2General
import ICFallDegrees.TruncatedFirstFall

noncomputable section
open Classical
namespace ICFallDegrees.PrimeM2Truncated
open MvPolynomial HomogeneousFirstFall GeneralQuadraticTop PrimeM2General
variable {K : Type*} [Field K]

/-- Exact FFD in K[z]/(z_i^p), with Koszul and Frobenius rows both declared
trivial. For prime fields p>3, p≥5 and 2≠0 hold automatically. -/
theorem exact_truncated_ffd (p : ℕ) [Fact p.Prime] [CharP K p]
    (large : 5≤p) (two_ne : (2 : K)≠0) (k : ℕ) (depth : 2≤k) (a b t : K) :
    TruncatedFirstFall.IsFirstFall p (inputs k depth a b t) (if 3*t^2+a=0 then 3 else 4) := by
  have tops : topForms (inputs k depth a b t)=q (active k depth) (-2*(t^2+a)) (-2*t) := by
    funext i; exact (degree_and_top k depth a b t i).2
  have degrees : (fun i => (inputs k depth a b t i).totalDegree)=(fun _ => 2) := by
    funext i; exact (degree_and_top k depth a b t i).1
  have determinant : -2*(t^2+a)-(-2*t)^2= -2*(3*t^2+a) := by ring
  unfold TruncatedFirstFall.IsFirstFall
  rw [tops,degrees]
  by_cases zero : 3*t^2+a=0
  · rw [if_pos zero]
    have coefficient : -2*(t^2+a)=(-2*t)^2 := by
      apply sub_eq_zero.mp
      rw [determinant,zero,mul_zero]
    rw [coefficient]
    exact TruncatedFirstFall.first_fall_three p large _ (active_injective k depth) (-2*t)
  · rw [if_neg zero]
    apply TruncatedFirstFall.first_fall_four p large _ (active_injective k depth)
    rw [determinant]
    exact mul_ne_zero (neg_ne_zero.mpr two_ne) zero

/-- All extra quotient rows are actual cycles for the tops of this system. -/
theorem frobenius_rows_are_syzygies (p : ℕ) [Fact p.Prime] [CharP K p]
    (k : ℕ) (depth : 2≤k) (a t : K) (i : Fin (2*k+1)) :
    ∑ j, TruncatedFirstFall.frobeniusRow p (q (active k depth) (-2*(t^2+a)) (-2*t)) i j*
      TruncatedFirstFall.image p (q (active k depth) (-2*(t^2+a)) (-2*t) j)=0 := by
  apply TruncatedFirstFall.frobenius_row_syzygy
  intro j
  refine Fin.lastCases ?_ (fun j => ?_) j
  · simp [q,top]
  · simp [q]

end ICFallDegrees.PrimeM2Truncated
