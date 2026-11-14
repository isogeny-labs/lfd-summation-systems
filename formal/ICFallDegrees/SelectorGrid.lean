import ICFallDegrees.ChainClosureBridge
import Mathlib.RingTheory.RootsOfUnity.PrimitiveRoots

noncomputable section

/-!
# The grid lemma for chain presentations

Over a field that contains a primitive `2^k`-th root of unity, a polynomial in the chain
variables that vanishes at every chain point (`u_{i,j} = x_i^(2^j)` with `x_i^(2^k) = 1`)
is derivable from the chain equations alone, within its own total degree.

This is the only place where "modulo the chains" enters the selector bound: every identity
of the group algebra is afterwards checked by evaluating at chain points.
-/

namespace ICFallDegrees.SelectorGrid

open MvPolynomial FiniteQuotient ChainQuotient ChainClosureBridge

variable {K : Type*} [Field K]

/-- The chain point of a tuple `x` of summand values. -/
def pt {m k : ℕ} (x : Fin m → K) : Fin (m*k) → K := fun v => x (row v)^2^(depth v : ℕ)

@[simp] theorem pt_var {m k : ℕ} (x : Fin m → K) (i : Fin m) (j : Fin k) :
    pt x (var m k i j) = x i^2^(j : ℕ) := by
  simp [pt]

theorem eval_expand (m k : ℕ) (x : Fin m → K) (q : P (m*k) K) :
    eval x (expand K m k q) = eval (pt x) q := by
  induction q using MvPolynomial.induction_on with
  | C c => simp [expand]
  | add p q hp hq => rw [map_add, map_add, map_add, hp, hq]
  | mul_X p i hp =>
      rw [map_mul, map_mul, map_mul, hp, expand_X, eval_X, map_pow, eval_X]
      rfl

/-- Chain equations vanish at chain points of roots of unity. -/
theorem eval_chainPolynomial (m k : ℕ) (x : Fin m → K) (hx : ∀ i, x i^2^k = 1) (v : Fin (m*k)) :
    eval (pt x) (chainPolynomial K m k v) = 0 := by
  rw [← eval_expand, expand_chainPolynomial]
  split_ifs
  · simp
  · simp [hx]

/-- Every element of the chain ideal vanishes at chain points of roots of unity. -/
theorem eval_of_mem_span (m k : ℕ) (x : Fin m → K) (hx : ∀ i, x i^2^k = 1) {p : P (m*k) K}
    (hp : p ∈ Ideal.span (chainInputs K m k)) : eval (pt x) p = 0 := by
  have le : Ideal.span (chainInputs K m k) ≤ RingHom.ker (eval (pt x)) := by
    apply Ideal.span_le.mpr
    rintro g ⟨v, rfl⟩
    exact eval_chainPolynomial m k x hx v
  exact le hp

/-- A squarefree polynomial vanishing on the grid of `2^k`-th roots of unity is zero. -/
theorem reduced_eq_zero_of_grid (m k : ℕ) {ζ : K} (root : IsPrimitiveRoot ζ (2^k))
    (p : P (m*k) K) (reduced : p ∈ reducedSpace (chainInputs K m k))
    (vanish : ∀ x : Fin m → K, (∀ i, x i^2^k = 1) → eval (pt x) p = 0) : p = 0 := by
  classical
  have multilinear : ∀ e ∈ p.support, ∀ v, e v ≤ 1 :=
    fun e he v => standard_le_one m k e (reduced e he) v
  have expansion : expand K m k p = ∑ e ∈ p.support, monomial (bitExp e) (coeff e p) := by
    conv_lhs => rw [p.as_sum]
    rw [map_sum]
    apply Finset.sum_congr rfl
    intro e _
    rw [expand_monomial]
  let S : Finset K := Polynomial.nthRootsFinset (2^k) (1 : K)
  have cardS : S.card = 2^k := root.card_nthRootsFinset
  have rootS : ∀ x ∈ S, x^2^k = 1 := by
    intro x hx
    exact (Polynomial.mem_nthRootsFinset (by positivity) (1 : K)).mp hx
  have vanish' : ∀ x : Fin m → K, (∀ i, x i ∈ S) → eval x (expand K m k p) = 0 := by
    intro x hx
    rw [eval_expand]
    exact vanish x (fun i => rootS _ (hx i))
  have degree : ∀ i, (expand K m k p).degreeOf i < S.card := by
    intro i
    rw [cardS, degreeOf_lt_iff (by positivity)]
    intro d hd
    rw [expansion] at hd
    obtain ⟨e, he, hd⟩ := Finset.mem_biUnion.mp (support_sum hd)
    have single := support_monomial_subset hd
    rw [Finset.mem_singleton] at single
    rw [single]
    exact bitExp_lt e (multilinear e he) i
  have zero : expand K m k p = 0 :=
    eq_zero_of_eval_zero_at_prod_finset (expand K m k p) (fun _ => S) degree vanish'
  by_contra nonzero
  have nonempty : p.support.Nonempty :=
    Finset.nonempty_iff_ne_empty.mpr (fun h => nonzero (support_eq_empty.mp h))
  obtain ⟨e₀, he₀⟩ := nonempty
  have coefficient := congrArg (coeff (bitExp e₀)) zero
  rw [expansion, coeff_sum, coeff_zero] at coefficient
  simp only [coeff_monomial] at coefficient
  rw [Finset.sum_eq_single e₀] at coefficient
  · rw [if_pos rfl] at coefficient
    exact (mem_support_iff.mp he₀) coefficient
  · intro e he hne
    rw [if_neg]
    intro heq
    exact hne (bitExp_injective e e₀ (multilinear e he) (multilinear e₀ he₀) heq)
  · intro h
    exact absurd he₀ h

/-- **Grid lemma.** A polynomial that vanishes at every chain point of roots of unity is in
the closure of the chain equations at its own total degree. -/
theorem grid_mem_closure (m k : ℕ) {ζ : K} (root : IsPrimitiveRoot ζ (2^k)) (p : P (m*k) K)
    (vanish : ∀ x : Fin m → K, (∀ i, x i^2^k = 1) → eval (pt x) p = 0) :
    InClosure (chainInputs K m k) p.totalDegree p := by
  have reduction := normalForm_reduction (chainInputs K m k) (chainInputs_nonzero m k) p
  have zero : normalForm (chainInputs K m k) (chainInputs_nonzero m k) p = 0 := by
    apply reduced_eq_zero_of_grid m k root _ (normalForm_reduced _ _ p)
    intro x hx
    have difference := eval_of_mem_span m k x hx reduction.sound
    rw [map_sub, vanish x hx, zero_sub, neg_eq_zero] at difference
    exact difference
  rwa [zero, sub_zero] at reduction

/-- Enlarging the input set and the degree keeps derivations. -/
theorem closure_of_subset {σ : Type*} {G F : Set (MvPolynomial σ K)} {d D : ℕ}
    {p : MvPolynomial σ K} (included : G ⊆ F) (hdD : d ≤ D) (h : InClosure G d p) :
    InClosure F D p :=
  InClosure.transport hdD (fun _ member degree => .input (included member) (degree.trans hdD)) h

/-- The grid lemma inside any input list that contains the chain equations. -/
theorem grid_mem (m k : ℕ) {ζ : K} (root : IsPrimitiveRoot ζ (2^k)) (F : Set (P (m*k) K))
    (included : chainInputs K m k ⊆ F) (D : ℕ) (p : P (m*k) K) (degree : p.totalDegree ≤ D)
    (vanish : ∀ x : Fin m → K, (∀ i, x i^2^k = 1) → eval (pt x) p = 0) : InClosure F D p :=
  closure_of_subset included degree (grid_mem_closure m k root p vanish)

end ICFallDegrees.SelectorGrid
