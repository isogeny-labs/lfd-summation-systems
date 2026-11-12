import Mathlib.Algebra.CharP.Lemmas
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Ring

/-! The scalar binary S3 fiber is never a zero polynomial and has at most
two roots. Boolean coordinate-fiber interpolation is a separate bridge. -/
namespace ICFallDegrees.BinaryFiber
variable {K : Type*} [Field K]

def quadratic (A B C y : K) : K := A*y^2+B*y+C

theorem quadratic_roots (A B C u : K) (nonzero : A≠0)
    (root : quadratic A B C u=0) :
    ∀ v, quadratic A B C v=0 → v=u ∨ v= -B/A-u := by
  intro v hv
  have factor : (v-u)*(A*(v+u)+B)=0 := by
    calc
      _ = quadratic A B C v-quadratic A B C u := by unfold quadratic; ring
      _ = 0 := by rw [hv,root]; ring
  rcases mul_eq_zero.mp factor with equal | sum
  · exact Or.inl (sub_eq_zero.mp equal)
  · right
    have value : v+u= -B/A := by
      apply (eq_div_iff nonzero).mpr
      linear_combination sum
    calc
      v = (v+u)-u := by ring
      _ = -B/A-u := by rw [value]

theorem coefficients_nonzero (x t b : K) (hb : b≠0) :
    (x+t)^2≠0 ∨ x*t≠0 ∨ (x*t)^2+b≠0 := by
  by_cases hA : (x+t)^2=0
  · right
    by_cases hB : x*t=0
    · right
      simpa [hB] using hb
    · exact Or.inl hB
  · exact Or.inl hA

/-- Once a root exists, a second root exhausts the scalar fiber. The
constant case is excluded using b!=0, including x=t=0. -/
theorem s3_fiber_two (x t b u : K) (hb : b≠0)
    (root : quadratic ((x+t)^2) (x*t) ((x*t)^2+b) u=0) :
    ∃ other : K, ∀ v,
      quadratic ((x+t)^2) (x*t) ((x*t)^2+b) v=0 → v=u ∨ v=other := by
  by_cases hA : (x+t)^2=0
  · have hB : x*t≠0 := by
      intro zero
      have C := coefficients_nonzero x t b hb
      simp only [hA,ne_eq,not_true_eq_false,false_or,zero] at C
      have bad := root
      simp [quadratic,hA,zero] at bad
      exact hb bad
    refine ⟨u,?_⟩
    intro v hv
    left
    have equal : (x*t)*(v-u)=0 := by
      have first := root
      simp only [quadratic,hA,zero_mul,zero_add] at first hv
      linear_combination hv-first
    exact sub_eq_zero.mp ((mul_eq_zero.mp equal).resolve_left hB)
  · exact ⟨-(x*t)/(x+t)^2-u,quadratic_roots _ _ _ u hA root⟩

theorem s3_identity [CharP K 2] (x y t b : K) :
    (x*y)^2+(x*t)^2+(y*t)^2+x*y*t+b=
      quadratic ((x+t)^2) (x*t) ((x*t)^2+b) y := by
  unfold quadratic
  rw [add_pow_char]
  ring

end ICFallDegrees.BinaryFiber
