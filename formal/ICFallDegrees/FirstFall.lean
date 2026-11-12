import ICFallDegrees.Closure

namespace ICFallDegrees

open MvPolynomial

variable {σ K : Type*} [CommRing K]

/-- Affine first *strict* fall, using the full iterated V-space.
This is distinct from homogeneous syzygies modulo Koszul syzygies. -/
def HasStrictFall (F : Set (MvPolynomial σ K)) (D : ℕ) : Prop :=
  ∃ f, InClosure F D f ∧ f.totalDegree < D ∧ ¬ InClosure F (D-1) f

def IsFirstStrictFall (F : Set (MvPolynomial σ K)) (D : ℕ) : Prop :=
  HasStrictFall F D ∧ ∀ e < D, ¬ HasStrictFall F e

theorem strictFall_le_lfd_bound {F : Set (MvPolynomial σ K)} {D d : ℕ}
    (fall : HasStrictFall F D) (bound : HasLFDBound F d) : D ≤ d := by
  obtain ⟨f,hf,degree,absent⟩ := fall
  by_contra bigger
  have hd : d < D := by omega
  have derived := bound f hf.sound
  have budget : max d f.totalDegree ≤ D-1 := by omega
  exact absent (derived.mono budget)

end ICFallDegrees
