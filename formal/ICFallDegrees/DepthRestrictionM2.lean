import ICFallDegrees.DepthRestriction
import ICFallDegrees.ChainClosureBridge

/-!
# Depth restriction for the normalized two-summand presentation

Instantiation of `DepthRestriction` for `PrimeM2General.inputs`: the substitution `u_{i,k} ↦ 1`
sends the normalized generator of depth `k+1` to the normalized generator of depth `k`.  Hence,
for every degree `D ≥ 2`, every field, curve and target:

* `one_of_deeper`: `1 ∈ V_D` at depth `k+1` implies `1 ∈ V_D` at depth `k`;
* `not_one_deeper`, `not_one_of_le`: an exclusion `1 ∉ V_D` at depth `k` holds at every depth
  `k' ≥ k`;
* `lfd_gt_of_le`: a unit-ideal instance at depth `k' ≥ k` has no last-fall-degree bound `≤ D`
  when `1 ∉ V_D` at depth `k`.

See `docs/proofs/U1-upper-bounds.md`, Proposition 3 and Corollary 3.1.
-/

noncomputable section

namespace ICFallDegrees.DepthRestrictionM2

open MvPolynomial FiniteQuotient ChainQuotient DepthRestriction

variable {K : Type*} [Field K]

theorem restrict_C (m k : ℕ) (c : K) : restrict K m k (C c) = C c := by
  simp [restrict, MvPolynomial.algebraMap_eq]

variable (k : ℕ) (depth : 2 ≤ k) (a b t : K)

/-- The four active variables (levels 0 and 1 of both blocks) keep their place. -/
theorem restrict_active_aux (n : Fin 4) (i : Fin 2) (j : ℕ) (hj : j < k)
    (source : (PrimeM2General.active (k+1) (by omega) n : ℕ) = j + (k+1) * (i : ℕ))
    (target : (PrimeM2General.active k depth n : ℕ) = j + k * (i : ℕ)) :
    restrict K 2 k (X (PrimeM2General.active (k+1) (by omega) n)) =
      X (PrimeM2General.active k depth n) := by
  have hv : PrimeM2General.active (k+1) (by omega) n = var 2 (k+1) i ⟨j, by omega⟩ := by
    apply Fin.ext
    rw [var_val]
    exact source
  have ht : PrimeM2General.active k depth n = var 2 k i ⟨j, hj⟩ := by
    apply Fin.ext
    rw [var_val]
    exact target
  rw [restrict_X, hv, restrictVar_var, dif_pos hj, ht]

theorem restrict_active0 :
    restrict K 2 k (X (PrimeM2General.active (k+1) (by omega) 0)) =
      X (PrimeM2General.active k depth 0) :=
  restrict_active_aux k depth 0 0 0 (by omega) (by simp [PrimeM2General.active])
    (by simp [PrimeM2General.active])

theorem restrict_active1 :
    restrict K 2 k (X (PrimeM2General.active (k+1) (by omega) 1)) =
      X (PrimeM2General.active k depth 1) :=
  restrict_active_aux k depth 1 0 1 (by omega) (by simp [PrimeM2General.active])
    (by simp [PrimeM2General.active])

theorem restrict_active2 :
    restrict K 2 k (X (PrimeM2General.active (k+1) (by omega) 2)) =
      X (PrimeM2General.active k depth 2) :=
  restrict_active_aux k depth 2 1 0 (by omega) (by simp [PrimeM2General.active])
    (by simp [PrimeM2General.active])

theorem restrict_active3 :
    restrict K 2 k (X (PrimeM2General.active (k+1) (by omega) 3)) =
      X (PrimeM2General.active k depth 3) :=
  restrict_active_aux k depth 3 1 1 (by omega) (by simp [PrimeM2General.active]; omega)
    (by simp [PrimeM2General.active]; omega)

/-- The normalized generator of depth `k+1` restricts to the normalized generator of depth `k`. -/
theorem restrict_normalized :
    restrict K 2 k (PrimeM2General.normalized (k+1) (by omega) a b t) =
      PrimeM2General.normalized k depth a b t := by
  unfold PrimeM2General.normalized PrimeM2General.tail GeneralQuadraticTop.top
  simp only [map_add, map_mul, restrict_C, restrict_active0 k depth, restrict_active1 k depth,
    restrict_active2 k depth, restrict_active3 k depth]

variable (D : ℕ)

/-- Degree-`D` derivations of the two-summand presentation restrict from depth `k+1` to `k`. -/
theorem restrict_closure (hD : 2 ≤ D) {p : P (2*(k+1)) K}
    (hp : InClosure (Set.range (PrimeM2General.inputs (k+1) (by omega) a b t)) D p) :
    InClosure (Set.range (PrimeM2General.inputs k depth a b t)) D (restrict K 2 k p) := by
  refine AffineTransport.map_closure _ restrict_linear _ _ D ?_ hp
  rintro f ⟨i, rfl⟩
  refine Fin.lastCases ?_ (fun j => ?_) i
  · intro hd
    simp only [PrimeM2General.inputs, Fin.lastCases_last] at hd ⊢
    have hdeg := (restrict_degree_le (K := K) (m := 2) (k := k)
      (PrimeM2General.normalized (k+1) (by omega) a b t)).trans hd
    rw [restrict_normalized k depth a b t] at hdeg ⊢
    exact .input ⟨Fin.last _, by simp [PrimeM2General.inputs]⟩ hdeg
  · intro _
    have hc := ChainClosureBridge.inputs_castSucc_mem_chainInputs (k+1) (by omega) a b t j
    rcases restrict_chainInput hc with hmem | hz
    · have hsub := chainInputs_two_subset k depth a b t hmem
      obtain ⟨v, hv⟩ := hmem
      exact .input hsub ((hv ▸ chainPolynomial_totalDegree_le 2 k v).trans hD)
    · rw [hz]
      exact .zero

/-- If `1` is derived at degree `D` at depth `k+1`, it is derived at degree `D` at depth `k`. -/
theorem one_of_deeper (hD : 2 ≤ D)
    (one : InClosure (Set.range (PrimeM2General.inputs (k+1) (by omega) a b t)) D 1) :
    InClosure (Set.range (PrimeM2General.inputs k depth a b t)) D 1 := by
  simpa using restrict_closure k depth a b t D hD one

/-- Exclusion of `1` at depth `k` persists at depth `k+1`. -/
theorem not_one_deeper (hD : 2 ≤ D)
    (excluded : ¬ InClosure (Set.range (PrimeM2General.inputs k depth a b t)) D 1) :
    ¬ InClosure (Set.range (PrimeM2General.inputs (k+1) (by omega) a b t)) D 1 :=
  fun one => excluded (one_of_deeper k depth a b t D hD one)

/-- Exclusion of `1` at depth `k` persists at every larger depth. -/
theorem not_one_of_le (hD : 2 ≤ D)
    (excluded : ¬ InClosure (Set.range (PrimeM2General.inputs k depth a b t)) D 1)
    (k' : ℕ) (hk : k ≤ k') (depth' : 2 ≤ k') :
    ¬ InClosure (Set.range (PrimeM2General.inputs k' depth' a b t)) D 1 := by
  induction k', hk using Nat.le_induction with
  | base => exact excluded
  | succ n hn ih => exact not_one_deeper n (depth.trans hn) a b t D hD (ih (depth.trans hn))

/-- A unit-ideal instance at depth `k' ≥ k` has no last-fall-degree bound `e ≤ D` when `1` is
excluded from the degree-`D` closure at depth `k`. -/
theorem lfd_gt_of_le (hD : 2 ≤ D)
    (excluded : ¬ InClosure (Set.range (PrimeM2General.inputs k depth a b t)) D 1)
    (k' : ℕ) (hk : k ≤ k') (depth' : 2 ≤ k')
    (unit : (1 : P (2*k') K) ∈ Ideal.span (Set.range (PrimeM2General.inputs k' depth' a b t))) :
    ∀ e ≤ D, ¬ HasLFDBound (Set.range (PrimeM2General.inputs k' depth' a b t)) e := by
  intro e he bound
  exact not_one_of_le k depth a b t D hD excluded k' hk depth'
    ((one_of_hasLFDBound unit bound).mono he)

end ICFallDegrees.DepthRestrictionM2
