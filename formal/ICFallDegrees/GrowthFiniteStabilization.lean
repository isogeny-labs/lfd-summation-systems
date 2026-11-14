import ICFallDegrees.Closure
import Mathlib.LinearAlgebra.FiniteDimensional.Lemmas
import Mathlib.Tactic

/-!
# N8: finite-stage stabilization of bounded-degree row spaces

A degree-D polynomial space on finitely many variables is finite dimensional.
Any increasing iterative construction of subspaces in such an ambient space
has a stationary adjacent pair after at most `finrank` iterations. If the
iteration comes from a deterministic step operator, the stationary pair
implies that all later stages agree.

These are general unconditional finite-dimensional linear algebra results,
not a claim yet that a particular iteration computes `InClosure`.
The latter equivalence is the next required interface.
-/

namespace ICFallDegrees.GrowthFiniteStabilization

variable {K V : Type*} [Field K] [AddCommGroup V] [Module K V]
  [FiniteDimensional K V]

/-- Every increasing finite-dimensional chain has a stationary adjacent
pair by stage `finrank K V`. No fixed-point assumption is hidden here. -/
theorem exists_stationary_step (W : ℕ → Submodule K V)
    (grow : ∀ n, W n ≤ W (n + 1)) :
    ∃ n ≤ Module.finrank K V, W n = W (n + 1) := by
  classical
  let d := Module.finrank K V
  by_contra h
  have hg : ∀ n, n ≤ d →
      Module.finrank K ↥(W n) < Module.finrank K ↥(W (n + 1)) := by
    intro n hn
    have hne : W n ≠ W (n + 1) := by
      intro heq
      exact h ⟨n, hn, heq⟩
    exact Submodule.finrank_lt_finrank_of_lt
      (lt_of_le_of_ne (grow n) hne)
  have hbound : ∀ n : ℕ, n ≤ d + 1 →
      n ≤ Module.finrank K ↥(W n) := by
    intro n
    induction n with
    | zero =>
        intro _
        omega
    | succ n ih =>
        intro hn
        have hlo := ih (by omega)
        have hstep := hg n (by omega)
        omega
  have hlarge := hbound (d + 1) (le_refl _)
  have hsmall := Submodule.finrank_le (W (d + 1))
  omega

/-- A deterministic row-extension operator stays fixed forever once
successive spaces agree. There is a certificate stage bounded by
the ambient finite dimension, not by the characteristic. -/
theorem eventual_stabilization (T : Submodule K V → Submodule K V)
    (W : ℕ → Submodule K V)
    (grow : ∀ n, W n ≤ W (n + 1))
    (step : ∀ n, W (n + 1) = T (W n)) :
    ∃ n ≤ Module.finrank K V,
      T (W n) = W n ∧ ∀ m, W (n + m) = W n := by
  obtain ⟨n, hn, heq⟩ := exists_stationary_step W grow
  have hfixed : T (W n) = W n := by
    calc
      T (W n) = W (n + 1) := (step n).symm
      _ = W n := heq.symm
  refine ⟨n, hn, hfixed, ?_⟩
  intro m
  induction m with
  | zero => simp
  | succ m ih =>
      calc
        W (n + (m + 1)) = W (n + m + 1) := congrArg W (Nat.add_succ n m)
        _ = T (W (n + m)) := step (n + m)
        _ = T (W n) := congrArg T ih
        _ = W n := hfixed

#print axioms exists_stationary_step
#print axioms eventual_stabilization

end ICFallDegrees.GrowthFiniteStabilization
