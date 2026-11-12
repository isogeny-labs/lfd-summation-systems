import ICFallDegrees.Closure
import Mathlib.Tactic.Ring

/-! Ordinary degree-three inverse/carry rewrites. Newly derived quadratic
rows may be multiplied at their actual degree; the initial derivation
budget is retained separately. No characteristic assumption is needed. -/
noncomputable section
namespace ICFallDegrees.CyclicInverseRewrite
open MvPolynomial
variable {σ K : Type*} [CommRing K]

private theorem degree_pair (u v w : MvPolynomial σ K)
    (hu : u.totalDegree≤1) (hv : v.totalDegree≤1) (hw : w.totalDegree≤1) :
    (u*v-w).totalDegree≤2 := by
  exact (totalDegree_sub_le _ _).trans (max_le
    ((totalDegree_mul _ _).trans (by omega)) (by omega))

private theorem degree_square (u v : MvPolynomial σ K)
    (hu : u.totalDegree≤1) (hv : v.totalDegree≤1) : (u^2-v).totalDegree≤2 := by
  have square := totalDegree_pow u 2
  exact (totalDegree_sub_le _ _).trans (max_le (by omega) (by omega))

theorem carry (F : Set (MvPolynomial σ K)) (u v next : MvPolynomial σ K)
    (hu : u.totalDegree≤1) (hv : v.totalDegree≤1) (hn : next.totalDegree≤1)
    (square : InClosure F 2 (u^2-next)) (unit : InClosure F 2 (u*v-1)) :
    InClosure F 3 (next*v-u) := by
  have sd := degree_square u next hu hn
  have ud : (u*v-1).totalDegree≤2 := by
    exact degree_pair u v 1 hu hv (by simp)
  have left := (square.mono (by decide : 2≤3)).mul v (by omega)
  have right := (unit.mono (by decide : 2≤3)).mul u (by omega)
  have identity : next*v-u=u*(u*v-1)-v*(u^2-next) := by ring
  rw [identity]
  exact right.sub left

theorem adjacent (F : Set (MvPolynomial σ K)) (u v next nextnext : MvPolynomial σ K)
    (hu : u.totalDegree≤1) (hv : v.totalDegree≤1)
    (hn : next.totalDegree≤1) (hnn : nextnext.totalDegree≤1)
    (known : InClosure F 3 (next*v-u))
    (square : InClosure F 2 (next^2-nextnext)) :
    InClosure F 3 (u*next-v*nextnext) := by
  have kd := degree_pair next v u hn hv hu
  have sd := degree_square next nextnext hn hnn
  have left := known.mul next (by omega)
  have right := (square.mono (by decide : 2≤3)).mul v (by omega)
  have identity : u*next-v*nextnext=v*(next^2-nextnext)-next*(next*v-u) := by ring
  rw [identity]
  exact right.sub left

theorem multiply_quadratic (F : Set (MvPolynomial σ K))
    (p a : MvPolynomial σ K) (D : ℕ) (room : 3≤D)
    (known : InClosure F 3 p) (quadratic : p.totalDegree≤2)
    (budget : a.totalDegree+2≤D) : InClosure F D (a*p) := by
  exact (known.mono room).mul a (by omega)

end ICFallDegrees.CyclicInverseRewrite
