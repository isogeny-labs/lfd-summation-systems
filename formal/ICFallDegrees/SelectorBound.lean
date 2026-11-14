import ICFallDegrees.SelectorFibre
import ICFallDegrees.SelectorAssembly
import Mathlib.LinearAlgebra.Lagrange
import Mathlib.FieldTheory.Finite.Basic

noncomputable section

/-!
# The selector bound: `LFD ≤ k + 2` for the two-summand prime chain presentation

Theorem 1 of `docs/proofs/U1-upper-bounds.md` (paper: Theorem `thm:selector`), for the
normalized presentation `PrimeM2General.inputs` at every depth `k ≥ 2`, over every field that
contains a primitive `2^k`-th root of unity, for every nonsingular cubic `x^3 + a x + b`
(`4 a^3 + 27 b^2 ≠ 0`) and every target `t`.

The projector is `N p = ∑_α E α * lin (I_α p)`, where `α` runs over the `2^k`-th roots of
unity, `E α` is the selector of `α` in the first summand, and `I_α p` is the Lagrange
interpolant of `β ↦ p(α, β)` on the fibre of `α` (the roots of `S_3(α, ·, t)` among the roots
of unity; at most two because the fibre polynomial is not zero). `N` has degree at most
`k + 1`, kills the ideal, and its corrections are derivable at degree `k + 2` by the fibre
lemma; `SelectorAssembly.lfd_bound_of_projector` concludes.

The bound is not claimed to be sharp for `k ≥ 3`, and it does not replace the bound `2k+1`
of `ChainQuotient.prime_m2_lfd_bound`, which needs neither roots of unity nor
nonsingularity.
-/

open Classical

namespace ICFallDegrees.SelectorBound

open MvPolynomial FiniteQuotient ChainQuotient SelectorGrid SelectorFibre
open PrimeM2General (inputs normalized)

variable {K : Type*} [Field K]

section Projector

variable (k : ℕ) (a b t : K)

/-- The `2^k`-th roots of unity. -/
def H (k : ℕ) : Finset K := Polynomial.nthRootsFinset (2^k) (1 : K)

theorem mem_H {α : K} : α ∈ H k ↔ α^2^k = 1 := by
  unfold H
  exact Polynomial.mem_nthRootsFinset (by positivity) (1 : K)

/-- The fibre of `α`: the roots of `S_3(α, ·, t)` among the roots of unity. -/
def fib (α : K) : Finset K :=
  (H k).filter (fun β => c₂ t α*β^2 + c₁ a b t α*β + c₀ a b t α = 0)

theorem mem_fib {α β : K} :
    β ∈ fib k a b t α ↔ β^2^k = 1 ∧ c₂ t α*β^2 + c₁ a b t α*β + c₀ a b t α = 0 := by
  unfold fib
  rw [Finset.mem_filter, mem_H]

theorem fib_card_le (two_ne : (2 : K) ≠ 0) (nonsingular : 4*a^3 + 27*b^2 ≠ 0) (α : K) :
    (fib k a b t α).card ≤ 2 := by
  obtain ⟨q, hq⟩ : ∃ q : Polynomial K, q = Polynomial.C (c₂ t α)*Polynomial.X^2 +
      Polynomial.C (c₁ a b t α)*Polynomial.X + Polynomial.C (c₀ a b t α) := ⟨_, rfl⟩
  have ev : ∀ β : K, q.eval β = c₂ t α*β^2 + c₁ a b t α*β + c₀ a b t α := by
    intro β
    rw [hq]
    simp
  have qne : q ≠ 0 := by
    intro h
    have e0 := ev 0
    have e1 := ev 1
    have em := ev (-1)
    rw [h, Polynomial.eval_zero] at e0 e1 em
    apply fibre_nonzero two_ne a b t nonsingular α
    have h0 : c₀ a b t α = 0 := by linear_combination -e0
    have h1 : c₁ a b t α = 0 := by
      have double : (2 : K)*c₁ a b t α = 0 := by linear_combination em - e1
      exact (mul_eq_zero.mp double).resolve_left two_ne
    have h2 : c₂ t α = 0 := by linear_combination -e1 - h1 - h0
    exact ⟨h2, h1, h0⟩
  have sub : (fib k a b t α).val ⊆ q.roots := by
    intro β hβ
    rw [Polynomial.mem_roots qne, Polynomial.IsRoot.def, ev]
    exact ((mem_fib k a b t).mp (Finset.mem_val.mp hβ)).2
  have deg : q.natDegree ≤ 2 := by
    rw [hq]
    exact Polynomial.natDegree_quadratic_le
  exact (Polynomial.card_le_degree_of_subset_roots sub).trans deg

theorem injOn (α : K) : Set.InjOn (id : K → K) (↑(fib k a b t α) : Set K) :=
  Function.injective_id.injOn

theorem interpolate_degree (two_ne : (2 : K) ≠ 0) (nonsingular : 4*a^3 + 27*b^2 ≠ 0) (α : K)
    (r : K → K) : (Lagrange.interpolate (fib k a b t α) id r).degree < 2 := by
  have h := Lagrange.degree_interpolate_lt r (injOn k a b t α)
  have card : ((fib k a b t α).card : WithBot ℕ) ≤ 2 := by
    exact_mod_cast fib_card_le k a b t two_ne nonsingular α
  exact lt_of_lt_of_le h card

theorem eval_of_degree_lt_two (f : Polynomial K) (hf : f.degree < 2) (β : K) :
    f.eval β = f.coeff 0 + f.coeff 1*β := by
  have hnat : f.natDegree < 2 := by
    by_cases h0 : f = 0
    · rw [h0, Polynomial.natDegree_zero]
      norm_num
    · exact (Polynomial.natDegree_lt_iff_degree_lt h0).mpr (by exact_mod_cast hf)
  rw [Polynomial.eval_eq_sum_range' hnat]
  simp [Finset.sum_range_succ]

/-- The interpolant of `β ↦ p(α, β)` on the fibre of `α`. -/
def interp (α : K) (p : P (2*k) K) : Polynomial K :=
  Lagrange.interpolate (fib k a b t α) id (fun β => eval (pt ![α, β]) p)

theorem interp_degree (two_ne : (2 : K) ≠ 0) (nonsingular : 4*a^3 + 27*b^2 ≠ 0) (α : K)
    (p : P (2*k) K) : (interp k a b t α p).degree < 2 :=
  interpolate_degree k a b t two_ne nonsingular α _

theorem eval_interp {α β : K} (hβ : β ∈ fib k a b t α) (p : P (2*k) K) :
    (interp k a b t α p).eval β = eval (pt ![α, β]) p :=
  Lagrange.eval_interpolate_at_node (fun β => eval (pt ![α, β]) p) (injOn k a b t α) hβ

/-- A polynomial `c₀ + c₁ X` as a polynomial in the first chain variable of summand 1. -/
def lin (k : ℕ) : Polynomial K →ₗ[K] P (2*k) K where
  toFun f := C (f.coeff 0) + C (f.coeff 1)*w k 0
  map_add' f g := by
    simp only [Polynomial.coeff_add, map_add]
    ring
  map_smul' c f := by
    simp only [Polynomial.coeff_smul, smul_eq_mul, map_mul, RingHom.id_apply, smul_eq_C_mul]
    ring

theorem lin_apply (f : Polynomial K) : lin k f = C (f.coeff 0) + C (f.coeff 1)*w k 0 := rfl

theorem lin_totalDegree_le (f : Polynomial K) : (lin k f).totalDegree ≤ 1 := by
  rw [lin_apply]
  refine (totalDegree_add _ _).trans (max_le (by simp) ?_)
  refine (totalDegree_mul _ _).trans ?_
  rw [totalDegree_C, zero_add]
  exact w_totalDegree_le k 0

theorem eval_lin (hk : 0 < k) (f : Polynomial K) (hf : f.degree < 2) (x : Fin 2 → K)
    (hx : x 1^2^k = 1) : eval (pt x) (lin k f) = f.eval (x 1) := by
  rw [lin_apply, map_add, map_mul, eval_C, eval_C, eval_w k x hx hk.le, pow_zero, pow_one,
    eval_of_degree_lt_two f hf]

/-- The projector onto the section `⊕_α span(R_α) e_α`. -/
def N : P (2*k) K →ₗ[K] P (2*k) K where
  toFun p := ∑ α ∈ H k, E 2 k 0 α*lin k (interp k a b t α p)
  map_add' p q := by
    show (∑ α ∈ H k, _) = (∑ α ∈ H k, _) + ∑ α ∈ H k, _
    rw [← Finset.sum_add_distrib]
    apply Finset.sum_congr rfl
    intro α _
    have e : (fun β => eval (pt ![α, β]) (p+q)) =
        (fun β => eval (pt ![α, β]) p) + (fun β => eval (pt ![α, β]) q) := by
      funext β
      simp
    have h : interp k a b t α (p+q) = interp k a b t α p + interp k a b t α q := by
      unfold interp
      rw [e, map_add]
    rw [h, map_add, mul_add]
  map_smul' c p := by
    show (∑ α ∈ H k, _) = c • ∑ α ∈ H k, _
    rw [Finset.smul_sum]
    apply Finset.sum_congr rfl
    intro α _
    have e : (fun β => eval (pt ![α, β]) (c • p)) = c • (fun β => eval (pt ![α, β]) p) := by
      funext β
      simp [smul_eval]
    have h : interp k a b t α (c • p) = c • interp k a b t α p := by
      unfold interp
      rw [e, map_smul]
    rw [h, map_smul, mul_smul_comm]

theorem N_apply (p : P (2*k) K) :
    N k a b t p = ∑ α ∈ H k, E 2 k 0 α*lin k (interp k a b t α p) := rfl

theorem pt_pair_zero (α β : K) (j : Fin k) : pt ![α, β] (var 2 k 0 j) = α^2^(j : ℕ) :=
  (pt_var _ _ _).trans rfl

theorem pt_pair_one (α β : K) (j : Fin k) : pt ![α, β] (var 2 k 1 j) = β^2^(j : ℕ) :=
  (pt_var _ _ _).trans rfl

theorem pair_root {α β : K} (hα : α^2^k = 1) (hβ : β^2^k = 1) :
    ∀ i : Fin 2, (![α, β] : Fin 2 → K) i^2^k = 1 := by
  rw [Fin.forall_fin_two]
  exact ⟨hα, hβ⟩

variable {k a b t}

theorem eval_N (hk : 0 < k) {ζ : K} (root : IsPrimitiveRoot ζ (2^k))
    (nonsingular : 4*a^3 + 27*b^2 ≠ 0) (p : P (2*k) K) (x : Fin 2 → K)
    (hx : ∀ i, x i^2^k = 1) :
    eval (pt x) (N k a b t p) = (interp k a b t (x 0) p).eval (x 1) := by
  have two_ne := two_ne_zero_of_root hk root
  rw [N_apply, map_sum, Finset.sum_eq_single (x 0)]
  · rw [map_mul, eval_E 2 k root 0 (x 0) (hx 0) x (hx 0), if_pos rfl, one_mul,
      eval_lin k hk _ (interp_degree k a b t two_ne nonsingular _ _) x (hx 1)]
  · intro α hα hne
    rw [map_mul, eval_E 2 k root 0 α ((mem_H k).mp hα) x (hx 0), if_neg (Ne.symm hne), zero_mul]
  · intro h
    exact absurd ((mem_H k).mpr (hx 0)) h

theorem eval_N_pair (hk : 0 < k) {ζ : K} (root : IsPrimitiveRoot ζ (2^k))
    (nonsingular : 4*a^3 + 27*b^2 ≠ 0) (p : P (2*k) K) {α β : K} (hα : α^2^k = 1)
    (hβ : β^2^k = 1) :
    eval (pt ![α, β]) (N k a b t p) = (interp k a b t α p).eval β :=
  (eval_N hk root nonsingular p ![α, β] (pair_root k hα hβ)).trans rfl

theorem N_totalDegree_le (p : P (2*k) K) : (N k a b t p).totalDegree ≤ k+1 := by
  rw [N_apply]
  refine (totalDegree_finset_sum _ _).trans (Finset.sup_le fun α _ => ?_)
  refine (totalDegree_mul _ _).trans ?_
  have h1 := E_totalDegree_le 2 k (0 : Fin 2) α
  have h2 := lin_totalDegree_le k (interp k a b t α p)
  omega

theorem N_absorb (hk : 0 < k) {ζ : K} (root : IsPrimitiveRoot ζ (2^k))
    (nonsingular : 4*a^3 + 27*b^2 ≠ 0) (i : Fin (2*k)) (p : P (2*k) K) :
    N k a b t (X i*N k a b t p) = N k a b t (X i*p) := by
  rw [N_apply k a b t (X i*N k a b t p), N_apply k a b t (X i*p)]
  apply Finset.sum_congr rfl
  intro α hα
  have h : interp k a b t α (X i*N k a b t p) = interp k a b t α (X i*p) := by
    unfold interp
    apply Lagrange.interpolate_eq_of_values_eq_on
    intro β hβ
    show eval (pt ![α, β]) (X i*N k a b t p) = eval (pt ![α, β]) (X i*p)
    rw [map_mul, map_mul, eval_N_pair hk root nonsingular p ((mem_H k).mp hα)
      ((mem_fib k a b t).mp hβ).1, eval_interp k a b t hβ p]
  rw [h]

theorem eval_ideal (depth : 2 ≤ k) {α β : K} (hα : α^2^k = 1) (hβ : β ∈ fib k a b t α)
    {f : P (2*k) K} (hf : f ∈ Ideal.span (Set.range (inputs k depth a b t))) :
    eval (pt ![α, β]) f = 0 := by
  obtain ⟨hβ1, hβ2⟩ := (mem_fib k a b t).mp hβ
  have hx := pair_root k hα hβ1
  have le : Ideal.span (Set.range (inputs k depth a b t)) ≤ RingHom.ker (eval (pt ![α, β])) := by
    apply Ideal.span_le.mpr
    rintro _ ⟨i, rfl⟩
    rw [SetLike.mem_coe, RingHom.mem_ker]
    refine Fin.lastCases ?_ (fun j => ?_) i
    · have last : inputs k depth a b t (Fin.last (2*k)) = normalized k depth a b t := by
        simp only [inputs, Fin.lastCases_last]
      rw [last, eval_normalized]
      exact hβ2
    · obtain ⟨v, hv⟩ := ChainClosureBridge.inputs_castSucc_mem_chainInputs k depth a b t j
      rw [← hv]
      exact eval_chainPolynomial 2 k _ hx v
  exact le hf

theorem N_vanish (depth : 2 ≤ k) (f : P (2*k) K)
    (hf : f ∈ Ideal.span (Set.range (inputs k depth a b t))) : N k a b t f = 0 := by
  rw [N_apply]
  apply Finset.sum_eq_zero
  intro α hα
  have h : interp k a b t α f = 0 := by
    unfold interp
    have values : ∀ β ∈ fib k a b t α,
        (fun β => eval (pt ![α, β]) f) β = (0 : K → K) β := by
      intro β hβ
      show eval (pt ![α, β]) f = 0
      exact eval_ideal depth ((mem_H k).mp hα) hβ hf
    rw [Lagrange.interpolate_eq_of_values_eq_on _ _ values, map_zero]
  rw [h, map_zero, mul_zero]

theorem N_one (depth : 2 ≤ k) {ζ : K} (root : IsPrimitiveRoot ζ (2^k))
    (nonsingular : 4*a^3 + 27*b^2 ≠ 0) :
    InClosure (Set.range (inputs k depth a b t)) (k+2) (1 - N k a b t 1) := by
  have hk : 0 < k := by omega
  have two_ne := two_ne_zero_of_root hk root
  have term : ∀ α ∈ H k, InClosure (Set.range (inputs k depth a b t)) (k+2)
      (E 2 k 0 α*(w k k - C ((interp k a b t α 1).coeff 0) -
        C ((interp k a b t α 1).coeff 1)*w k 0)) := by
    intro α hα
    refine ((fibre (a := a) (b := b) (t := t) depth root nonsingular α ((mem_H k).mp hα))
      k le_rfl _ _).1 ?_
    intro β hβ hq
    have mem : β ∈ fib k a b t α := (mem_fib k a b t).mpr ⟨hβ, hq⟩
    rw [← eval_of_degree_lt_two _ (interp_degree k a b t two_ne nonsingular α 1) β,
      eval_interp k a b t mem, map_one, hβ]
  have unity : InClosure (Set.range (inputs k depth a b t)) (k+2)
      (1 - ∑ α ∈ H k, E 2 k 0 α) := by
    apply grid_mem 2 k root _ (chainInputs_two_subset k depth a b t) (k+2)
    · refine (totalDegree_sub_le _ _).trans (max_le (by simp) ?_)
      refine (totalDegree_finset_sum _ _).trans (Finset.sup_le fun α _ => ?_)
      exact (E_totalDegree_le 2 k 0 α).trans (by omega)
    · intro x hx
      rw [map_sub, map_one, map_sum]
      have each : ∀ α ∈ H k, eval (pt x) (E 2 k 0 α) = if x 0 = α then 1 else 0 :=
        fun α hα => eval_E 2 k root 0 α ((mem_H k).mp hα) x (hx 0)
      rw [Finset.sum_congr rfl each, Finset.sum_ite_eq, if_pos ((mem_H k).mpr (hx 0)), sub_self]
  have ident : 1 - N k a b t 1 = (1 - ∑ α ∈ H k, E 2 k 0 α) +
      ∑ α ∈ H k, E 2 k 0 α*(w k k - C ((interp k a b t α 1).coeff 0) -
        C ((interp k a b t α 1).coeff 1)*w k 0) := by
    rw [N_apply, w_top]
    simp only [lin_apply, mul_sub, mul_add, mul_one, Finset.sum_sub_distrib,
      Finset.sum_add_distrib]
    ring
  rw [ident]
  exact (closureSpace (Set.range (inputs k depth a b t)) (k+2)).add_mem unity
    ((closureSpace (Set.range (inputs k depth a b t)) (k+2)).sum_mem term)

/-- The correction for a chain variable of the selected summand: it acts on `E α` by the
scalar `α^(2^j)`, so the correction vanishes on the grid. -/
theorem step_zero (depth : 2 ≤ k) {ζ : K} (root : IsPrimitiveRoot ζ (2^k))
    (nonsingular : 4*a^3 + 27*b^2 ≠ 0) (α : K) (hα : α ∈ H k) (p : P (2*k) K) (j : Fin k) :
    InClosure (Set.range (inputs k depth a b t)) (k+2)
      (X (var 2 k 0 j)*(E 2 k 0 α*lin k (interp k a b t α p)) -
        E 2 k 0 α*lin k (interp k a b t α (X (var 2 k 0 j)*N k a b t p))) := by
  have hk : 0 < k := by omega
  have hα' := (mem_H k).mp hα
  have fdeg : (interp k a b t α p).degree < (fib k a b t α).card :=
    Lagrange.degree_interpolate_lt _ (injOn k a b t α)
  have hg : interp k a b t α (X (var 2 k 0 j)*N k a b t p) =
      (α^2^(j : ℕ)) • interp k a b t α p := by
    refine (Lagrange.eq_interpolate_of_eval_eq (s := fib k a b t α) (v := id)
      (r := fun β => eval (pt ![α, β]) (X (var 2 k 0 j)*N k a b t p)) (injOn k a b t α)
      (lt_of_le_of_lt (Polynomial.degree_smul_le _ _) fdeg) ?_).symm
    intro β hβ
    show Polynomial.eval β ((α^2^(j : ℕ)) • interp k a b t α p) =
      eval (pt ![α, β]) (X (var 2 k 0 j)*N k a b t p)
    rw [map_mul, eval_X, pt_pair_zero, eval_N_pair hk root nonsingular p hα'
      ((mem_fib k a b t).mp hβ).1, Polynomial.eval_smul, smul_eq_mul]
  rw [hg, map_smul]
  apply grid_mem 2 k root _ (chainInputs_two_subset k depth a b t) (k+2)
  · have h1 := E_totalDegree_le 2 k (0 : Fin 2) α
    have h2 := lin_totalDegree_le k (interp k a b t α p)
    refine (totalDegree_sub_le _ _).trans (max_le ?_ ?_)
    · refine (totalDegree_mul _ _).trans ?_
      have h3 := totalDegree_mul (E 2 k 0 α) (lin k (interp k a b t α p))
      rw [totalDegree_X]
      omega
    · refine (totalDegree_mul _ _).trans ?_
      have h3 := totalDegree_smul_le (α^2^(j : ℕ)) (lin k (interp k a b t α p))
      omega
  · intro x hx
    rw [map_sub, map_mul, map_mul, map_mul, eval_X, pt_var, eval_E 2 k root 0 α hα' x (hx 0),
      smul_eval]
    split_ifs with h
    · rw [h]
      ring
    · ring

/-- The correction for a chain variable of the other summand: the fibre lemma. -/
theorem step_one (depth : 2 ≤ k) {ζ : K} (root : IsPrimitiveRoot ζ (2^k))
    (nonsingular : 4*a^3 + 27*b^2 ≠ 0) (α : K) (hα : α ∈ H k) (p : P (2*k) K) (j : Fin k) :
    InClosure (Set.range (inputs k depth a b t)) (k+2)
      (X (var 2 k 1 j)*(E 2 k 0 α*lin k (interp k a b t α p)) -
        E 2 k 0 α*lin k (interp k a b t α (X (var 2 k 1 j)*N k a b t p))) := by
  have hk : 0 < k := by omega
  have two_ne := two_ne_zero_of_root hk root
  have hα' := (mem_H k).mp hα
  obtain ⟨f, hf⟩ : ∃ f : Polynomial K, f = interp k a b t α p := ⟨_, rfl⟩
  have fdeg : f.degree < 2 := by
    rw [hf]
    exact interp_degree k a b t two_ne nonsingular α p
  obtain ⟨I₁, hI₁⟩ : ∃ I₁ : Polynomial K,
      I₁ = Lagrange.interpolate (fib k a b t α) id (fun β => β^2^(j : ℕ)) := ⟨_, rfl⟩
  obtain ⟨I₂, hI₂⟩ : ∃ I₂ : Polynomial K,
      I₂ = Lagrange.interpolate (fib k a b t α) id (fun β => β^(2^(j : ℕ)+1)) := ⟨_, rfl⟩
  have deg₁ : I₁.degree < 2 := by
    rw [hI₁]
    exact interpolate_degree k a b t two_ne nonsingular α _
  have deg₂ : I₂.degree < 2 := by
    rw [hI₂]
    exact interpolate_degree k a b t two_ne nonsingular α _
  have hg : interp k a b t α (X (var 2 k 1 j)*N k a b t p) = f.coeff 0 • I₁ + f.coeff 1 • I₂ := by
    have values : ∀ β ∈ fib k a b t α,
        (fun β => eval (pt ![α, β]) (X (var 2 k 1 j)*N k a b t p)) β =
          (f.coeff 0 • (fun β : K => β^2^(j : ℕ)) +
            f.coeff 1 • (fun β : K => β^(2^(j : ℕ)+1))) β := by
      intro β hβ
      show eval (pt ![α, β]) (X (var 2 k 1 j)*N k a b t p) =
        f.coeff 0*β^2^(j : ℕ) + f.coeff 1*β^(2^(j : ℕ)+1)
      rw [map_mul, eval_X, pt_pair_one, eval_N_pair hk root nonsingular p hα'
        ((mem_fib k a b t).mp hβ).1, ← hf, eval_of_degree_lt_two f fdeg β, pow_succ]
      ring
    unfold interp
    rw [Lagrange.interpolate_eq_of_values_eq_on _ _ values, map_add, map_smul, map_smul,
      ← hI₁, ← hI₂]
  have row₁ : InClosure (Set.range (inputs k depth a b t)) (k+2)
      (E 2 k 0 α*(w k j - C (I₁.coeff 0) - C (I₁.coeff 1)*w k 0)) := by
    refine ((fibre (a := a) (b := b) (t := t) depth root nonsingular α hα')
      j j.isLt.le _ _).1 ?_
    intro β hβ hq
    have mem : β ∈ fib k a b t α := (mem_fib k a b t).mpr ⟨hβ, hq⟩
    rw [← eval_of_degree_lt_two I₁ deg₁ β, hI₁]
    exact Lagrange.eval_interpolate_at_node _ (injOn k a b t α) mem
  have row₂ : InClosure (Set.range (inputs k depth a b t)) (k+2)
      (E 2 k 0 α*(w k 0*w k j - C (I₂.coeff 0) - C (I₂.coeff 1)*w k 0)) := by
    refine ((fibre (a := a) (b := b) (t := t) depth root nonsingular α hα')
      j j.isLt.le _ _).2 ?_
    intro β hβ hq
    have mem : β ∈ fib k a b t α := (mem_fib k a b t).mpr ⟨hβ, hq⟩
    rw [← eval_of_degree_lt_two I₂ deg₂ β, hI₂]
    exact Lagrange.eval_interpolate_at_node _ (injOn k a b t α) mem
  have hX : (X (var 2 k 1 j) : P (2*k) K) = w k j := (w_lt k j.isLt).symm
  have ident : X (var 2 k 1 j)*(E 2 k 0 α*lin k f) -
      E 2 k 0 α*lin k (f.coeff 0 • I₁ + f.coeff 1 • I₂) =
      f.coeff 0 • (E 2 k 0 α*(w k j - C (I₁.coeff 0) - C (I₁.coeff 1)*w k 0)) +
        f.coeff 1 • (E 2 k 0 α*(w k 0*w k j - C (I₂.coeff 0) - C (I₂.coeff 1)*w k 0)) := by
    rw [hX]
    simp only [lin_apply, Polynomial.coeff_add, Polynomial.coeff_smul, smul_eq_mul, map_add,
      map_mul, smul_eq_C_mul]
    ring
  rw [hg, ← hf, ident]
  exact (row₁.smul _).add (row₂.smul _)

theorem N_step (depth : 2 ≤ k) {ζ : K} (root : IsPrimitiveRoot ζ (2^k))
    (nonsingular : 4*a^3 + 27*b^2 ≠ 0) (i : Fin (2*k)) (p : P (2*k) K) :
    InClosure (Set.range (inputs k depth a b t)) (k+2)
      (X i*N k a b t p - N k a b t (X i*N k a b t p)) := by
  have e₁ := N_apply k a b t p
  have e₂ := N_apply k a b t (X i*N k a b t p)
  have ident : X i*N k a b t p - N k a b t (X i*N k a b t p) =
      ∑ α ∈ H k, (X i*(E 2 k 0 α*lin k (interp k a b t α p)) -
        E 2 k 0 α*lin k (interp k a b t α (X i*N k a b t p))) := by
    rw [e₂, Finset.sum_sub_distrib, ← Finset.mul_sum, ← e₁]
  rw [ident]
  apply (closureSpace (Set.range (inputs k depth a b t)) (k+2)).sum_mem
  intro α hα
  have both : ∀ (r : Fin 2) (j : Fin k),
      InClosure (Set.range (inputs k depth a b t)) (k+2)
        (X (var 2 k r j)*(E 2 k 0 α*lin k (interp k a b t α p)) -
          E 2 k 0 α*lin k (interp k a b t α (X (var 2 k r j)*N k a b t p))) := by
    rw [Fin.forall_fin_two]
    exact ⟨step_zero depth root nonsingular α hα p, step_one depth root nonsingular α hα p⟩
  have h := both (ChainQuotient.row i) (ChainQuotient.depth i)
  rw [var_row_depth] at h
  exact h

end Projector

/-- **Selector bound.** For every depth `k ≥ 2`, every field containing a primitive `2^k`-th
root of unity, every nonsingular cubic `x^3 + a x + b` and every target `t`, the normalized
two-summand chain presentation has last fall degree at most `k + 2`. -/
theorem lfd_le_depth_add_two (k : ℕ) (depth : 2 ≤ k) {ζ : K} (root : IsPrimitiveRoot ζ (2^k))
    (a b t : K) (nonsingular : 4*a^3 + 27*b^2 ≠ 0) :
    HasLFDBound (Set.range (inputs k depth a b t)) (k+2) :=
  SelectorAssembly.lfd_bound_of_projector _ (N k a b t) (k+1) (k+2) le_rfl
    N_totalDegree_le
    (N_absorb (by omega) root nonsingular)
    (N_vanish depth)
    (N_one depth root nonsingular)
    (N_step depth root nonsingular)

/-- A prime field with `2^k ∣ p - 1` contains a primitive `2^k`-th root of unity. -/
theorem exists_primitiveRoot_zmod (p : ℕ) [Fact p.Prime] (k : ℕ) (dvd : 2^k ∣ p - 1) :
    ∃ ζ : ZMod p, IsPrimitiveRoot ζ (2^k) := by
  obtain ⟨g, hg⟩ := IsCyclic.exists_generator (α := (ZMod p)ˣ)
  have order : orderOf g = p - 1 := by
    rw [orderOf_eq_card_of_forall_mem_zpowers hg, Nat.card_eq_fintype_card, ZMod.card_units]
  have prim : IsPrimitiveRoot g (p - 1) := by
    rw [← order]
    exact IsPrimitiveRoot.orderOf g
  obtain ⟨c, hc⟩ := dvd
  have pos : 0 < p - 1 := by
    have := (Fact.out : p.Prime).two_le
    omega
  have power := IsPrimitiveRoot.pow pos prim (by rw [hc, mul_comm] : p - 1 = c*2^k)
  exact ⟨((g^c : (ZMod p)ˣ) : ZMod p), IsPrimitiveRoot.coe_units_iff.mpr power⟩

/-- Prime-field entry point: `2^k ∣ p - 1` (so `p ≥ 5` for `k ≥ 2`), nonsingular curve,
every target. -/
theorem prime_field_lfd_le_depth_add_two (p : ℕ) [Fact p.Prime] (k : ℕ) (depth : 2 ≤ k)
    (root : 2^k ∣ p - 1) (a b t : ZMod p) (nonsingular : 4*a^3 + 27*b^2 ≠ 0) :
    HasLFDBound (Set.range (inputs k depth a b t)) (k+2) := by
  obtain ⟨ζ, hζ⟩ := exists_primitiveRoot_zmod p k root
  exact lfd_le_depth_add_two k depth hζ a b t nonsingular

end ICFallDegrees.SelectorBound
