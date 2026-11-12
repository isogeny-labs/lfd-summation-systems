import ICFallDegrees.ChainQuotient
import ICFallDegrees.AffineTransport

/-!
# Restriction of a chain presentation to the next smaller depth

For `m` summands the chain presentation of depth `k+1` maps onto the one of depth `k` by the
substitution `u_{i,k} ↦ 1` (restriction of functions from the subgroup of order `2^(k+1)` to the
subgroup of order `2^k`).  The substitution sends every chain equation of depth `k+1` to a chain
equation of depth `k` or to zero, and it does not raise degrees.  Hence a degree-`D` derivation
at depth `k+1` restricts to a degree-`D` derivation at depth `k`:

* `restrict_closure`: closure membership is transported, for arbitrary further inputs;
* `one_of_deeper`: if `1` is derived at degree `D` at depth `k+1`, it is at depth `k`;
* `not_one_deeper`: a depth-`k` system with `1 ∉ V_D` keeps `1 ∉ V_D` at depth `k+1`;
* `lfd_bound_of_deeper_unit`: a last-fall-degree bound for a unit ideal at depth `k+1` holds at
  depth `k`.

Nothing here uses the group-algebra model or the characteristic.  See
`docs/proofs/U1-upper-bounds.md`, Proposition 3.
-/

noncomputable section

namespace ICFallDegrees.DepthRestriction

open MvPolynomial FiniteQuotient ChainQuotient

variable (K : Type*) [Field K] (m k : ℕ)

/-- Image of a depth-`k+1` chain variable: the same variable at depth `k`, or `1` at the top. -/
def restrictVar (v : Fin (m*(k+1))) : P (m*k) K :=
  if h : (depth v).val < k then X (var m k (row v) ⟨(depth v).val, h⟩) else 1

/-- The substitution `u_{i,k} ↦ 1`. -/
def restrict : P (m*(k+1)) K →ₐ[K] P (m*k) K := aeval (restrictVar K m k)

variable {K m k}

@[simp] theorem restrict_X (v : Fin (m*(k+1))) : restrict K m k (X v) = restrictVar K m k v := by
  simp [restrict]

theorem restrict_linear (v : Fin (m*(k+1))) : (restrict K m k (X v)).totalDegree ≤ 1 := by
  rw [restrict_X]
  unfold restrictVar
  split_ifs <;> simp

theorem restrict_degree_le (p : P (m*(k+1)) K) :
    (restrict K m k p).totalDegree ≤ p.totalDegree :=
  AffineTransport.degree_le _ restrict_linear p

theorem restrictVar_var (i : Fin m) (j : Fin (k+1)) :
    restrictVar K m k (var m (k+1) i j) =
      if h : j.val < k then X (var m k i ⟨j.val, h⟩) else 1 := by
  unfold restrictVar
  by_cases h : j.val < k
  · have h' : (depth (var m (k+1) i j)).val < k := by simpa using h
    rw [dif_pos h', dif_pos h]
    congr 2
    · simp
    · apply Fin.ext; simp
  · have h' : ¬ (depth (var m (k+1) i j)).val < k := by simpa using h
    rw [dif_neg h', dif_neg h]

theorem chainTail_var (n : ℕ) (i : Fin m) (j : Fin n) :
    chainTail K m n (var m n i j) =
      if h : j.val+1 < n then X (var m n i ⟨j.val+1, h⟩) else 1 := by
  unfold chainTail
  by_cases h : j.val+1 < n
  · have h' : (depth (var m n i j)).val+1 < n := by simpa using h
    rw [dif_pos h', dif_pos h]
    congr 2
    · simp
    · apply Fin.ext; simp
  · have h' : ¬ (depth (var m n i j)).val+1 < n := by simpa using h
    rw [dif_neg h', dif_neg h]

/-- A chain equation of depth `k+1` restricts to the chain equation of the same variable at
depth `k`, and the top equation `u_{i,k}^2 - 1` restricts to zero. -/
theorem restrict_chainPolynomial (i : Fin m) (j : Fin (k+1)) :
    restrict K m k (chainPolynomial K m (k+1) (var m (k+1) i j)) =
      if h : j.val < k then chainPolynomial K m k (var m k i ⟨j.val, h⟩) else 0 := by
  unfold chainPolynomial
  rw [map_sub, map_pow, restrict_X, restrictVar_var, chainTail_var]
  by_cases h : j.val < k
  · have h1 : j.val+1 < k+1 := by omega
    rw [dif_pos h, dif_pos h1, dif_pos h, restrict_X, restrictVar_var, chainTail_var]
  · have h1 : ¬ j.val+1 < k+1 := by omega
    rw [dif_neg h, dif_neg h1, dif_neg h]
    simp

theorem restrict_chainInput {f : P (m*(k+1)) K} (hf : f ∈ chainInputs K m (k+1)) :
    restrict K m k f ∈ chainInputs K m k ∨ restrict K m k f = 0 := by
  obtain ⟨v, rfl⟩ := hf
  rw [← var_row_depth v, restrict_chainPolynomial]
  by_cases h : (depth v).val < k
  · rw [dif_pos h]
    exact Or.inl ⟨_, rfl⟩
  · rw [dif_neg h]
    exact Or.inr rfl

variable (R : Set (P (m*(k+1)) K)) (R' : Set (P (m*k) K)) (D : ℕ)

/-- Degree-`D` derivations at depth `k+1` restrict to degree-`D` derivations at depth `k`,
provided the restriction of each further input is derivable at depth `k`. -/
theorem restrict_closure (hD : 2 ≤ D)
    (inputs : ∀ r ∈ R, r.totalDegree ≤ D →
      InClosure (chainInputs K m k ∪ R') D (restrict K m k r))
    {p : P (m*(k+1)) K} (hp : InClosure (chainInputs K m (k+1) ∪ R) D p) :
    InClosure (chainInputs K m k ∪ R') D (restrict K m k p) := by
  refine AffineTransport.map_closure _ restrict_linear _ _ D ?_ hp
  intro f hf hd
  rcases hf with hf | hf
  · rcases restrict_chainInput hf with hc | hz
    · obtain ⟨v, hv⟩ := hc
      rw [← hv]
      exact .input (Or.inl ⟨v, rfl⟩) ((chainPolynomial_totalDegree_le m k v).trans hD)
    · rw [hz]
      exact .zero
  · exact inputs f hf hd

/-- If `1` is derived at degree `D` at depth `k+1`, it is derived at degree `D` at depth `k`. -/
theorem one_of_deeper (hD : 2 ≤ D)
    (inputs : ∀ r ∈ R, r.totalDegree ≤ D →
      InClosure (chainInputs K m k ∪ R') D (restrict K m k r))
    (one : InClosure (chainInputs K m (k+1) ∪ R) D 1) :
    InClosure (chainInputs K m k ∪ R') D 1 := by
  simpa using restrict_closure R R' D hD inputs one

/-- Contrapositive: exclusion of `1` at depth `k` persists at depth `k+1`. -/
theorem not_one_deeper (hD : 2 ≤ D)
    (inputs : ∀ r ∈ R, r.totalDegree ≤ D →
      InClosure (chainInputs K m k ∪ R') D (restrict K m k r))
    (excluded : ¬ InClosure (chainInputs K m k ∪ R') D 1) :
    ¬ InClosure (chainInputs K m (k+1) ∪ R) D 1 :=
  fun one => excluded (one_of_deeper R R' D hD inputs one)

/-- A last-fall-degree bound for a unit ideal at depth `k+1` is a bound at depth `k`. -/
theorem lfd_bound_of_deeper_unit (hD : 2 ≤ D)
    (inputs : ∀ r ∈ R, r.totalDegree ≤ D →
      InClosure (chainInputs K m k ∪ R') D (restrict K m k r))
    (unit : (1 : P (m*(k+1)) K) ∈ Ideal.span (chainInputs K m (k+1) ∪ R))
    (bound : HasLFDBound (chainInputs K m (k+1) ∪ R) D) :
    HasLFDBound (chainInputs K m k ∪ R') D :=
  hasLFDBound_of_one (one_of_deeper R R' D hD inputs (one_of_hasLFDBound unit bound))

end ICFallDegrees.DepthRestriction
