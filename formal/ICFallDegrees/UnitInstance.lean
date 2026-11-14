import ICFallDegrees.DepthRestrictionM2
import ICFallDegrees.GridQuotient
import Mathlib.RingTheory.Polynomial.Cyclotomic.Basic
import Mathlib.NumberTheory.LegendreSymbol.JacobiSymbol
import Mathlib.Tactic.NormNum.LegendreSymbol
import Mathlib.Tactic.NormNum.Prime

/-!
# The unit instance: a two-summand chain system without solutions

For a curve `y² = x³ + a x` (`b = 0`) over a field `K` of characteristic other than two, and a
target abscissa `t` which is not a square in `K` but for which `t³ + a t` is a nonzero square, the
third summation polynomial `S₃(x₁, x₂, t)` has no zero with `x₁`, `x₂` nonzero squares of `K`
(`S3_ne_zero`).  The proof is elementary (no elliptic-curve group law):

* `S₃` is symmetric, and as a quadratic in any one argument its discriminant is `16 c c'`, where
  `c(x) = x³ + a x` is evaluated at the two other arguments (`S3_disc`).  Hence a zero forces
  `c(x₁)` and `c(x₂)` to be squares `y₁²`, `y₂²` (`cubic_isSquare`);
* off the diagonal, `(x₁-x₂)² t = (x₁+x₂)(x₁x₂+a) ± 2 y₁ y₂`, and then
  `x₁ x₂ t (x₁-x₂)² = (x₂ y₁ ± x₁ y₂)²`, so `t` is a square;
* on the diagonal `S₃ = (x²-a)² - 4 c(x) t`, so `t` is a square or `x = 0`.

If moreover `K` contains a primitive `2^(k+1)`-th root of unity, every root of `x^(2^k) = 1` is
a nonzero square and `x^(2^k) - 1` splits with distinct roots; the chain relations put every
variable in that group, so the system `PrimeM2General.inputs k _ a 0 t` has no zero on the grid
and, by the combinatorial Nullstellensatz, generates the unit ideal (`unit_of_primitiveRoot`,
`unit_zmod`).  The concrete instances are `p = 65537`, `(a, t) = (1, 4271)`, `2 ≤ k ≤ 15`
(`unit_p65537`) and `p = 163841`, `(a, t) = (3, 73730)`, `2 ≤ k ≤ 14` (`unit_p163841`).

Everything in this file is kernel-checked; no recorded computation enters.
-/

noncomputable section

namespace ICFallDegrees.UnitInstance

/-! ## Part 1: no zero of `S₃` at nonzero squares -/

section Scalar

variable {K : Type*} [Field K]

/-- The third summation polynomial of `y² = x³ + a x`, as a function of three abscissas. -/
def S3 (a x₁ x₂ t : K) : K :=
  (x₁-x₂)^2*t^2-2*((x₁+x₂)*(x₁*x₂+a))*t+(x₁*x₂-a)^2

theorem S3_symm_left (a x₁ x₂ t : K) : S3 a x₁ x₂ t = S3 a x₂ x₁ t := by
  unfold S3; ring

theorem S3_symm_right (a x₁ x₂ t : K) : S3 a x₁ x₂ t = S3 a x₁ t x₂ := by
  unfold S3; ring

/-- The discriminant of `S₃` as a quadratic in its last argument is `16 c(x₁) c(x₂)`. -/
theorem S3_disc (a x₁ x₂ t : K) :
    ((x₁-x₂)^2*t-(x₁+x₂)*(x₁*x₂+a))^2 =
      4*(x₁^3+a*x₁)*(x₂^3+a*x₂)+(x₁-x₂)^2*S3 a x₁ x₂ t := by
  unfold S3; ring

/-- On the diagonal, `S₃ = (x²-a)² - 4 c(x) t`. -/
theorem S3_diag (a x t : K) : S3 a x x t = (x^2-a)^2-4*(x^3+a*x)*t := by
  unfold S3; ring

/-- At a zero of `S₃` whose last argument has a nonzero square cubic value, the cubic value of
the second argument is a square.  No distinctness is needed. -/
theorem cubic_isSquare (two_ne : (2 : K) ≠ 0) {a x₁ x₂ t Y : K}
    (zero : S3 a x₁ x₂ t = 0) (curve : Y^2 = t^3+a*t) (hY : Y ≠ 0) :
    ∃ y : K, y^2 = x₂^3+a*x₂ := by
  have disc := S3_disc a x₂ t x₁
  rw [← S3_symm_right, S3_symm_left, zero, mul_zero, add_zero] at disc
  refine ⟨((x₂-t)^2*x₁-(x₂+t)*(x₂*t+a))/(2*Y), ?_⟩
  rw [div_pow, disc, ← curve]
  field_simp
  ring

/-- **The no-solution lemma.**  Over a field with `2 ≠ 0`, let `Y² = t³ + a t` with `Y ≠ 0`
and `t` not a square.  Then `S₃(x₁, x₂, t) ≠ 0` for all nonzero squares `x₁`, `x₂`. -/
theorem S3_ne_zero (two_ne : (2 : K) ≠ 0) {a t Y : K} (curve : Y^2 = t^3+a*t) (hY : Y ≠ 0)
    (nonsquare : ¬ IsSquare t) {x₁ x₂ : K} (sq₁ : IsSquare x₁) (sq₂ : IsSquare x₂)
    (ne₁ : x₁ ≠ 0) (ne₂ : x₂ ≠ 0) : S3 a x₁ x₂ t ≠ 0 := by
  intro zero
  obtain ⟨y₂, hy₂⟩ := cubic_isSquare two_ne zero curve hY
  obtain ⟨y₁, hy₁⟩ := cubic_isSquare two_ne ((S3_symm_left a x₂ x₁ t).trans zero) curve hY
  apply nonsquare
  by_cases diagonal : x₁ = x₂
  · subst diagonal
    rw [S3_diag, ← hy₁] at zero
    by_cases hy : y₁ = 0
    · exfalso
      rw [hy] at zero hy₁
      have h1 : x₁^2-a = 0 := by
        have : (x₁^2-a)^2 = 0 := by linear_combination zero
        exact pow_eq_zero_iff (two_ne_zero) |>.mp this
      have h2 : x₁^2+a = 0 := by
        have : x₁*(x₁^2+a) = 0 := by linear_combination -hy₁
        exact (mul_eq_zero.mp this).resolve_left ne₁
      have h3 : (2 : K)*x₁^2 = 0 := by linear_combination h1+h2
      exact ne₁ (pow_eq_zero_iff (two_ne_zero) |>.mp ((mul_eq_zero.mp h3).resolve_left two_ne))
    · refine ⟨(x₁^2-a)/(2*y₁), ?_⟩
      field_simp
      linear_combination -zero
  · obtain ⟨s₁, hs₁⟩ := sq₁
    obtain ⟨s₂, hs₂⟩ := sq₂
    have hd : x₁-x₂ ≠ 0 := sub_ne_zero.mpr diagonal
    have hs₁0 : s₁ ≠ 0 := by rintro rfl; exact ne₁ (by simpa using hs₁)
    have hs₂0 : s₂ ≠ 0 := by rintro rfl; exact ne₂ (by simpa using hs₂)
    have disc := S3_disc a x₁ x₂ t
    rw [zero, mul_zero, add_zero] at disc
    have square : ((x₁-x₂)^2*t-(x₁+x₂)*(x₁*x₂+a))^2 = (2*y₁*y₂)^2 := by
      rw [disc, ← hy₁, ← hy₂]; ring
    rcases sq_eq_sq_iff_eq_or_eq_neg.mp square with root | root
    · have key : x₁*x₂*t*(x₁-x₂)^2 = (x₂*y₁+x₁*y₂)^2 := by
        linear_combination x₁*x₂*root-x₂^2*hy₁-x₁^2*hy₂
      refine ⟨(x₂*y₁+x₁*y₂)/(s₁*s₂*(x₁-x₂)), ?_⟩
      field_simp
      rw [hs₁, hs₂] at key ⊢
      linear_combination key
    · have key : x₁*x₂*t*(x₁-x₂)^2 = (x₂*y₁-x₁*y₂)^2 := by
        linear_combination x₁*x₂*root-x₂^2*hy₁-x₁^2*hy₂
      refine ⟨(x₂*y₁-x₁*y₂)/(s₁*s₂*(x₁-x₂)), ?_⟩
      field_simp
      rw [hs₁, hs₂] at key ⊢
      linear_combination key

end Scalar

/-! ## Part 2: an ideal without zeros on a finite grid is the unit ideal -/

section Grid

open MvPolynomial

variable {K : Type*} [Field K]

/-- Finite-grid Nullstellensatz in the form used here: an ideal which contains the grid
polynomials `∏_{s ∈ S i} (X i - s)` and has, at every grid point, a member that does not vanish
there, is the unit ideal. -/
theorem one_mem_of_no_grid_zero {σ : Type*} [Fintype σ] [DecidableEq σ] (S : σ → Finset K)
    (I : Ideal (MvPolynomial σ K))
    (grid : ∀ i, ∏ s ∈ S i, (X i - C s) ∈ I)
    (nozero : ∀ x : σ → K, (∀ i, x i ∈ S i) → ∃ f ∈ I, eval x f ≠ 0) :
    (1 : MvPolynomial σ K) ∈ I := by
  classical
  by_cases nonempty : ∀ i, (S i).Nonempty
  swap
  · push_neg at nonempty
    obtain ⟨i, hi⟩ := nonempty
    rw [Finset.not_nonempty_iff_eq_empty] at hi
    simpa [hi] using grid i
  choose! f hfI hf using nozero
  let e : (σ → K) → MvPolynomial σ K :=
    fun p => ∏ i, ∏ s ∈ (S i).erase (p i), (X i - C s)
  let g : MvPolynomial σ K :=
    ∑ p ∈ Fintype.piFinset S, C ((eval p (f p * e p))⁻¹) * (f p * e p)
  have gI : g ∈ I := Ideal.sum_mem _ fun p hp =>
    I.mul_mem_left _ (I.mul_mem_right _ (hfI p (Fintype.mem_piFinset.mp hp)))
  have vanishing : ∀ x : σ → K, (∀ i, x i ∈ S i) → eval x (1 - g) = 0 := by
    intro x hx
    have he : eval x (e x) ≠ 0 := by
      simp only [e, map_prod, map_sub, eval_X, eval_C]
      refine Finset.prod_ne_zero_iff.mpr fun i _ => Finset.prod_ne_zero_iff.mpr fun s hs => ?_
      exact sub_ne_zero.mpr (Finset.ne_of_mem_erase hs).symm
    have hval : eval x (f x * e x) ≠ 0 := by
      rw [map_mul]; exact mul_ne_zero (hf x hx) he
    have sum : eval x g = 1 := by
      simp only [g, map_sum]
      rw [Finset.sum_eq_single x]
      · rw [map_mul, eval_C, inv_mul_cancel₀ hval]
      · intro p _ hpx
        obtain ⟨i, hi⟩ := Function.ne_iff.mp hpx
        have : eval x (e p) = 0 := by
          simp only [e, map_prod, map_sub, eval_X, eval_C]
          refine Finset.prod_eq_zero (Finset.mem_univ i) ?_
          exact Finset.prod_eq_zero (Finset.mem_erase.mpr ⟨fun h => hi h.symm, hx i⟩)
            (sub_self _)
        simp only [map_mul, this, mul_zero]
      · intro hnot
        exact absurd (Fintype.mem_piFinset.mpr hx) hnot
    rw [map_sub, map_one, sum, sub_self]
  obtain ⟨h, -, combination⟩ :=
    combinatorial_nullstellensatz_exists_linearCombination S nonempty (1 - g) vanishing
  have diff : 1 - g ∈ I := by
    rw [combination, Finsupp.linearCombination_apply]
    refine Ideal.sum_mem _ fun i _ => ?_
    exact I.mul_mem_left (h i) (grid i)
  simpa using I.add_mem diff gI

end Grid

/-! ## The chain system -/

section Chain

open MvPolynomial PrimeM2General

variable {K : Type*} [Field K] (k : ℕ) (depth : 2 ≤ k) (a b t : K)

/-- The chain relation `x_i² - x_{i+1}` (inside a block) is in the ideal. -/
theorem chain_step_mem (i : Fin (2*k)) (h : ¬ (i.val+1 = k ∨ i.val+1 = 2*k)) :
    (X i^2 - X ⟨i.val+1, by have := i.isLt; omega⟩ : P k K) ∈
      Ideal.span (Set.range (inputs k depth a b t)) := by
  have member : inputs k depth a b t (Fin.castSucc i) ∈
      Ideal.span (Set.range (inputs k depth a b t)) := Ideal.subset_span ⟨_, rfl⟩
  simp only [inputs, Fin.lastCases_castSucc, PrimeM2General.chainTail, dif_neg h] at member
  rwa [← sub_eq_add_neg] at member

/-- The chain relation `x_i² - 1` (at the end of a block) is in the ideal. -/
theorem chain_end_mem (i : Fin (2*k)) (h : i.val+1 = k ∨ i.val+1 = 2*k) :
    (X i^2 - 1 : P k K) ∈ Ideal.span (Set.range (inputs k depth a b t)) := by
  have member : inputs k depth a b t (Fin.castSucc i) ∈
      Ideal.span (Set.range (inputs k depth a b t)) := Ideal.subset_span ⟨_, rfl⟩
  simp only [inputs, Fin.lastCases_castSucc, PrimeM2General.chainTail, dif_pos h] at member
  rwa [← sub_eq_add_neg] at member

/-- A variable at distance `d` from the end of its block satisfies `x^(2^d) = 1` modulo the
ideal. -/
theorem chain_pow_mem (n : ℕ) : ∀ i : Fin (2*k),
    (i.val+(n+1) = k ∨ (k ≤ i.val ∧ i.val+(n+1) = 2*k)) →
    (X i^2^(n+1) - 1 : P k K) ∈ Ideal.span (Set.range (inputs k depth a b t)) := by
  induction n with
  | zero =>
      intro i h
      simpa using chain_end_mem k depth a b t i (by omega)
  | succ n ih =>
      intro i h
      have hi := i.isLt
      have notend : ¬ (i.val+1 = k ∨ i.val+1 = 2*k) := by omega
      have step := chain_step_mem k depth a b t i notend
      have next := ih ⟨i.val+1, by omega⟩ (by
        show (i.val+1)+(n+1) = k ∨ (k ≤ i.val+1 ∧ (i.val+1)+(n+1) = 2*k)
        omega)
      obtain ⟨c, hc⟩ := sub_dvd_pow_sub_pow (X i^2 : P k K) (X ⟨i.val+1, by omega⟩) (2^(n+1))
      have split : (X i^2^(n+1+1) - 1 : P k K) =
          ((X i^2)^2^(n+1) - X ⟨i.val+1, by omega⟩^2^(n+1)) +
            (X ⟨i.val+1, by omega⟩^2^(n+1) - 1) := by
        rw [← pow_mul, ← pow_succ']; ring
      rw [split, hc]
      exact Ideal.add_mem _ (Ideal.mul_mem_right _ _ step) next

/-- Every variable satisfies `x^(2^k) = 1` modulo the ideal. -/
theorem var_pow_mem (i : Fin (2*k)) :
    (X i^2^k - 1 : P k K) ∈ Ideal.span (Set.range (inputs k depth a b t)) := by
  have hi := i.isLt
  obtain ⟨n, hn, hnk⟩ : ∃ n, (i.val+(n+1) = k ∨ (k ≤ i.val ∧ i.val+(n+1) = 2*k)) ∧ n+1 ≤ k := by
    by_cases small : i.val < k
    · exact ⟨k-i.val-1, Or.inl (by omega), by omega⟩
    · exact ⟨2*k-i.val-1, Or.inr ⟨by omega, by omega⟩, by omega⟩
  have base := chain_pow_mem k depth a b t n i hn
  have exponent : 2^k = 2^(n+1)*2^(k-(n+1)) := by
    rw [← pow_add]; congr 1; omega
  obtain ⟨c, hc⟩ := sub_dvd_pow_sub_pow (X i^2^(n+1) : P k K) 1 (2^(k-(n+1)))
  rw [exponent, pow_mul]
  rw [one_pow] at hc
  rw [hc]
  exact Ideal.mul_mem_right _ _ base

/-- Evaluation of the raw `S₃` generator for `b = 0`. -/
theorem eval_rawS3 (x : Fin (2*k) → K) :
    eval x (rawS3 k depth a 0 t) = S3 a (x (active k depth 0)) (x (active k depth 2)) t := by
  simp only [rawS3, S3, map_sub, map_add, map_mul, map_pow, eval_X, eval_C, map_ofNat]
  ring

/-- **The unit ideal.**  If `K` contains a primitive `2^(k+1)`-th root of unity, `t` is not a
square and `t³ + a t` is a nonzero square, the two-summand chain system of `y² = x³ + a x` with
target abscissa `t` generates the unit ideal at depth `k`. -/
theorem unit_of_primitiveRoot {ζ : K} (root : IsPrimitiveRoot ζ (2^(k+1))) {Y : K}
    (curve : Y^2 = t^3+a*t) (hY : Y ≠ 0) (nonsquare : ¬ IsSquare t) :
    (1 : P k K) ∈ Ideal.span (Set.range (inputs k depth a 0 t)) := by
  classical
  have pos : 0 < 2^k := by positivity
  have root' : IsPrimitiveRoot (ζ^2) (2^k) := by
    have quotient : 2^(k+1)/2 = 2^k := by
      rw [pow_succ]; exact Nat.mul_div_cancel _ two_pos
    have := root.pow_of_dvd two_ne_zero (dvd_pow_self 2 (Nat.succ_ne_zero k))
    rwa [quotient] at this
  have two_ne : (2 : K) ≠ 0 := by
    have minus : IsPrimitiveRoot (ζ^2^k) 2 := by
      have quotient : 2^(k+1)/2^k = 2 := by
        rw [pow_succ]; exact Nat.mul_div_cancel_left _ pos
      have := root.pow_of_dvd (pow_ne_zero k two_ne_zero) (pow_dvd_pow 2 (Nat.le_succ k))
      rwa [quotient] at this
    intro two
    have h1 : ζ^2^k = 1 := by
      rw [minus.eq_neg_one_of_two_right]; linear_combination -two
    exact minus.ne_one one_lt_two h1
  have square_of_root : ∀ x : K, x^2^k = 1 → IsSquare x ∧ x ≠ 0 := by
    intro x hx
    haveI : NeZero (2^k) := ⟨pos.ne'⟩
    obtain ⟨i, -, rfl⟩ := root'.eq_pow_of_pow_eq_one hx
    exact ⟨⟨ζ^i, by ring⟩, pow_ne_zero _ (pow_ne_zero _ (root.ne_zero (by positivity)))⟩
  refine one_mem_of_no_grid_zero (fun _ => Polynomial.nthRootsFinset (2^k) (1 : K)) _ ?_ ?_
  · intro i
    have identity := congrArg (Polynomial.aeval (X i : P k K))
      (Polynomial.X_pow_sub_one_eq_prod pos root')
    simp only [map_sub, map_pow, Polynomial.aeval_X, map_one, map_prod, Polynomial.aeval_C,
      MvPolynomial.algebraMap_eq] at identity
    rw [← identity]
    exact var_pow_mem k depth a 0 t i
  · intro x hx
    have roots : ∀ i, x i^2^k = 1 := fun i => (Polynomial.mem_nthRootsFinset pos 1).mp (hx i)
    by_cases normal : eval x (normalized k depth a 0 t) = 0
    swap
    · exact ⟨_, Ideal.subset_span ⟨Fin.last _, by simp [inputs]⟩, normal⟩
    have first := chain_step_mem k depth a 0 t (active k depth 0)
      (by simp [active]; omega)
    have second := chain_step_mem k depth a 0 t (active k depth 2)
      (by simp [active]; omega)
    have first_eq : (⟨(active k depth 0).val+1, by simp [active]; omega⟩ : Fin (2*k)) =
        active k depth 1 := by apply Fin.ext; simp [active]
    have second_eq : (⟨(active k depth 2).val+1, by simp [active]; omega⟩ : Fin (2*k)) =
        active k depth 3 := by apply Fin.ext; simp [active]
    rw [first_eq] at first
    rw [second_eq] at second
    by_cases c₁ : eval x (X (active k depth 0)^2 - X (active k depth 1) : P k K) = 0
    swap
    · exact ⟨_, first, c₁⟩
    by_cases c₂ : eval x (X (active k depth 2)^2 - X (active k depth 3) : P k K) = 0
    swap
    · exact ⟨_, second, c₂⟩
    exfalso
    have identity := congrArg (eval x) (raw_normalization_identity k depth a 0 t)
    rw [map_sub, map_add, map_mul, map_mul, c₁, c₂, normal, mul_zero, mul_zero, add_zero,
      sub_zero, eval_rawS3] at identity
    obtain ⟨sq₁, ne₁⟩ := square_of_root _ (roots (active k depth 0))
    obtain ⟨sq₂, ne₂⟩ := square_of_root _ (roots (active k depth 2))
    exact S3_ne_zero two_ne curve hY nonsquare sq₁ sq₂ ne₁ ne₂ identity

end Chain

/-! ## Prime fields -/

section PrimeField

open PrimeM2General

/-- A prime field with `n ∣ p - 1` contains a primitive `n`-th root of unity. -/
theorem exists_primitiveRoot_zmod (p : ℕ) [Fact p.Prime] (n : ℕ) (hn : n ∣ p-1) :
    ∃ ζ : ZMod p, IsPrimitiveRoot ζ n := by
  have prime := (Fact.out : p.Prime)
  obtain ⟨g, hg⟩ := IsCyclic.exists_ofOrder_eq_natCard (α := (ZMod p)ˣ)
  rw [Nat.card_eq_fintype_card, ZMod.card_units_eq_totient, Nat.totient_prime prime] at hg
  obtain ⟨m, hm⟩ := hn
  have hp : p-1 ≠ 0 := by have := prime.two_le; omega
  have hm0 : m ≠ 0 := by rintro rfl; exact hp (by simpa using hm)
  have order : orderOf (g^m) = n := by
    rw [orderOf_pow_of_dvd hm0 (by rw [hg, hm]; exact dvd_mul_left m n), hg, hm,
      Nat.mul_div_cancel _ (Nat.pos_of_ne_zero hm0)]
  refine ⟨((g^m : (ZMod p)ˣ) : ZMod p), ?_⟩
  rw [IsPrimitiveRoot.coe_units_iff, ← order]
  exact IsPrimitiveRoot.orderOf _

/-- **The unit ideal over a prime field** with `2^(k+1) ∣ p - 1`. -/
theorem unit_zmod (p : ℕ) [Fact p.Prime] (k : ℕ) (depth : 2 ≤ k) (hp : 2^(k+1) ∣ p-1)
    {a t Y : ZMod p} (curve : Y^2 = t^3+a*t) (hY : Y ≠ 0) (nonsquare : ¬ IsSquare t) :
    (1 : P k (ZMod p)) ∈ Ideal.span (Set.range (inputs k depth a 0 t)) := by
  obtain ⟨ζ, root⟩ := exists_primitiveRoot_zmod p (2^(k+1)) hp
  exact unit_of_primitiveRoot k depth a t root curve hY nonsquare

instance fact_prime_65537 : Fact (Nat.Prime 65537) := ⟨by norm_num⟩

instance fact_prime_163841 : Fact (Nat.Prime 163841) := ⟨by norm_num⟩

theorem nonsquare_4271 : ¬ IsSquare (4271 : ZMod 65537) := by
  have symbol : jacobiSym 4271 65537 = -1 := by norm_num
  simpa using ZMod.nonsquare_of_jacobiSym_eq_neg_one symbol

theorem curve_4271 : (11720 : ZMod 65537)^2 = 4271^3+1*4271 := by decide

theorem nonzero_11720 : (11720 : ZMod 65537) ≠ 0 := by decide

/-- The paper's unit instance: `y² = x³ + x` over `𝔽_65537`, target `(4271, 11720)`, generates
the unit ideal at every depth `2 ≤ k ≤ 15`. -/
theorem unit_p65537 (k : ℕ) (depth : 2 ≤ k) (hk : k ≤ 15) :
    (1 : P k (ZMod 65537)) ∈
      Ideal.span (Set.range (inputs k depth (1 : ZMod 65537) 0 4271)) := by
  refine unit_zmod 65537 k depth ?_ curve_4271 nonzero_11720 nonsquare_4271
  have : 65537-1 = 2^16 := by norm_num
  rw [this]
  exact pow_dvd_pow 2 (by omega)

theorem nonsquare_73730 : ¬ IsSquare (73730 : ZMod 163841) := by
  have symbol : jacobiSym 73730 163841 = -1 := by norm_num
  simpa using ZMod.nonsquare_of_jacobiSym_eq_neg_one symbol

theorem curve_73730 : (28291 : ZMod 163841)^2 = 73730^3+3*73730 := by decide

theorem nonzero_28291 : (28291 : ZMod 163841) ≠ 0 := by decide

/-- The second unit instance: `y² = x³ + 3x` over `𝔽_163841`, target `(73730, 28291)`, generates
the unit ideal at every depth `2 ≤ k ≤ 14`. -/
theorem unit_p163841 (k : ℕ) (depth : 2 ≤ k) (hk : k ≤ 14) :
    (1 : P k (ZMod 163841)) ∈
      Ideal.span (Set.range (inputs k depth (3 : ZMod 163841) 0 73730)) := by
  refine unit_zmod 163841 k depth ?_ curve_73730 nonzero_28291 nonsquare_73730
  have : 163841-1 = 2^15*5 := by norm_num
  rw [this]
  exact Dvd.dvd.mul_right (pow_dvd_pow 2 (by omega)) 5

end PrimeField

end ICFallDegrees.UnitInstance
