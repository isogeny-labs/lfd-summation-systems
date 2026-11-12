import ICFallDegrees.Closure
import Mathlib.Tactic.Ring

/-! Uniform degree-four refutations along a square chain. All reductions
are charged in the ordinary polynomial ring, including terminal equations. -/
noncomputable section
namespace ICFallDegrees.SquareChainRefutation
open MvPolynomial
variable {σ K : Type*} [Field K]

def row (z : ℕ → MvPolynomial σ K) (a : K) (j : ℕ) : MvPolynomial σ K :=
  z (j+1)-C (2*a^(2^j))*z j+C ((a^(2^j))^2)

def multiplier (z : ℕ → MvPolynomial σ K) (a : K) (j : ℕ) : MvPolynomial σ K :=
  z (j+1)+C (2*a^(2^j))*z j+C ((a^(2^j))^2)

theorem row_degree (z : ℕ → MvPolynomial σ K) (bounded : ∀ j, (z j).totalDegree≤2)
    (a : K) (j : ℕ) : (row z a j).totalDegree≤2 := by
  unfold row
  apply (totalDegree_add _ _).trans
  apply max_le _ (by rw [totalDegree_C]; omega)
  apply (totalDegree_sub_le _ _).trans
  exact max_le (bounded _) (by
      have h := totalDegree_mul (C (2*a^(2^j))) (z j)
      rw [totalDegree_C,zero_add] at h
      exact h.trans (bounded j))

theorem multiplier_degree (z : ℕ → MvPolynomial σ K) (bounded : ∀ j, (z j).totalDegree≤2)
    (a : K) (j : ℕ) : (multiplier z a j).totalDegree≤2 := by
  unfold multiplier
  apply (totalDegree_add _ _).trans
  apply max_le _ (by rw [totalDegree_C]; omega)
  exact (totalDegree_add _ _).trans
    (max_le (bounded _) (by
      have h := totalDegree_mul (C (2*a^(2^j))) (z j)
      rw [totalDegree_C,zero_add] at h
      exact h.trans (bounded j)))

theorem propagation_identity (z : ℕ → MvPolynomial σ K) (a : K) (j : ℕ) :
    row z a (j+1)=multiplier z a j*row z a j+
      C (4*(a^(2^j))^2)*(z j^2-z (j+1))-(z (j+1)^2-z (j+2)) := by
  have power : a^(2^(j+1))=(a^(2^j))^2 := by rw [pow_succ,pow_mul]
  simp only [row,multiplier,power,map_mul,map_pow,map_ofNat]
  ring

/-- The degree bound is independent of chain depth. The assumptions are
literal closure certificates for the square-chain defects and initial row. -/
theorem one_at_four (F : Set (MvPolynomial σ K)) (z : ℕ → MvPolynomial σ K)
    (bounded : ∀ j, (z j).totalDegree≤2) (a : K) (k : ℕ)
    (defects : ∀ j≤k, InClosure F 4 (z j^2-z (j+1)))
    (initial : InClosure F 4 (row z a 0)) (terminal : z k=1) (terminal_next : z (k+1)=1)
    (outside : a^(2^k)≠1) : InClosure F 4 1 := by
  have rows : ∀ j≤k, InClosure F 4 (row z a j) := by
    intro j
    induction j with
    | zero => intro _; exact initial
    | succ j ih =>
        intro bound
        have prev := ih (by omega)
        have multiplied := prev.mul (multiplier z a j) (by
          have := multiplier_degree z bounded a j
          have := row_degree z bounded a j
          omega)
        have correction := (defects j (by omega)).smul (4*(a^(2^j))^2)
        rw [smul_eq_C_mul] at correction
        rw [propagation_identity]
        exact (multiplied.add correction).sub (defects (j+1) bound)
  have final := rows k le_rfl
  have value : row z a k=C ((1-a^(2^k))^2) := by
    simp only [row,terminal,terminal_next,map_mul,map_pow,map_sub,map_one,map_ofNat]
    ring
  rw [value] at final
  have nonzero : (1-a^(2^k))^2≠0 := pow_ne_zero _ (sub_ne_zero.mpr (Ne.symm outside))
  have scaled := final.smul (((1-a^(2^k))^2)⁻¹)
  simpa only [smul_eq_C_mul,←map_mul,inv_mul_cancel₀ nonzero,map_one] using scaled

end ICFallDegrees.SquareChainRefutation
