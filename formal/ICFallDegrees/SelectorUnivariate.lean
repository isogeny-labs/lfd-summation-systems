import Mathlib.RingTheory.RootsOfUnity.PrimitiveRoots
import Mathlib.Algebra.Polynomial.Splits
import Mathlib.Algebra.Polynomial.Degree.SmallDegree
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.FieldSimp

/-!
# The univariate lemma of the selector bound

Lemma 2.3 of `docs/proofs/U1-upper-bounds.md`, in an abstract commutative `K`-algebra.
The data are a "selector" `E`, a chain `w 0, …, w k` with `w k = 1`, a row space `Cl`
(the degree-bounded closure) and a space `Low` of rows that may still be multiplied by a
chain variable. From the row `E * (w 1 - s * w 0 + n)` of a monic quadratic
`Q = X^2 - s X + n` the recursion produces, for every `j ≤ k`, the rows

  `E * (w j - a_j - b_j * w 0)`  and  `E * (w 0 * w j - …)`,

where `a_j + b_j X` is the remainder of `X^(2^j)` modulo `Q`. The two rows at `j = k` span the
image of `1 - X^(2^k)` in `K[X]/(Q)`; when the field contains the `2^k`-th roots of unity this
gives every `E * (c₀ + c₁ * w 0)` with `c₀ + c₁ X` vanishing on the roots of `Q` among them.
No gcd is computed: the three cases are a nonzero determinant, a root `(1-a)/b` of `Q` that
is a root of unity, and `Q ∣ X^(2^k) - 1`.
-/

namespace ICFallDegrees.SelectorUnivariate

variable {K : Type*} [Field K]

/-- The remainder `a_j + b_j X` of `X^(2^j)` modulo `X^2 - s X + n`, as the pair `(a_j, b_j)`. -/
def rem (s n : K) : ℕ → K × K
  | 0 => (0, 1)
  | j+1 => ((rem s n j).1^2 - n*(rem s n j).2^2,
      2*(rem s n j).1*(rem s n j).2 + s*(rem s n j).2^2)

theorem rem_root (s n θ : K) (hθ : θ^2 - s*θ + n = 0) (j : ℕ) :
    (rem s n j).1 + (rem s n j).2*θ = θ^2^j := by
  induction j with
  | zero => simp [rem]
  | succ j ih =>
      have sq : θ^2^(j+1) = (θ^2^j)^2 := by rw [pow_succ, pow_mul]
      rw [sq, ← ih]
      simp only [rem]
      linear_combination (-(rem s n j).2^2)*hθ

section Polynomial
open Polynomial

theorem rem_dvd (s n : K) (j : ℕ) :
    (X^2 - C s*X + C n : K[X]) ∣ X^2^j - C (rem s n j).1 - C (rem s n j).2*X := by
  induction j with
  | zero => simp [rem]
  | succ j ih =>
      obtain ⟨u, hu⟩ := ih
      refine ⟨u*(X^2^j + C (rem s n j).1 + C (rem s n j).2*X) + C ((rem s n j).2^2), ?_⟩
      have sq : (X : K[X])^2^(j+1) = (X^2^j)^2 := by rw [pow_succ, pow_mul]
      rw [sq]
      simp only [rem, map_sub, map_add, map_mul, map_pow, map_ofNat]
      linear_combination (X^2^j + C (rem s n j).1 + C (rem s n j).2*X)*hu

/-- A monic quadratic dividing `X^N - 1` has two distinct roots among the `N`-th roots of
unity, when the field contains a primitive one. -/
theorem two_roots {N : ℕ} {ζ : K} (root : IsPrimitiveRoot ζ N) (hN : 0 < N) (s n : K)
    (dvd : (X^2 - C s*X + C n : K[X]) ∣ X^N - C 1) :
    ∃ β₁ β₂ : K, β₁ ≠ β₂ ∧ β₁^N = 1 ∧ β₂^N = 1 ∧
      β₁^2 - s*β₁ + n = 0 ∧ β₂^2 - s*β₂ + n = 0 := by
  classical
  obtain ⟨Q, hQ⟩ : ∃ Q : K[X], Q = X^2 - C s*X + C n := ⟨_, rfl⟩
  rw [← hQ] at dvd
  have form : Q = C 1*X^2 + C (-s)*X + C n := by
    rw [hQ, C_1, C_neg]; ring
  have degQ : Q.natDegree = 2 := by
    rw [form]; exact natDegree_quadratic one_ne_zero
  have Qne : Q ≠ 0 := by
    intro h
    rw [h, natDegree_zero] at degQ
    exact absurd degQ (by norm_num)
  have big_ne : (X^N - C 1 : K[X]) ≠ 0 := X_pow_sub_C_ne_zero hN (1 : K)
  have big_splits : Splits (RingHom.id K) (X^N - C 1 : K[X]) := by
    rw [splits_iff_card_roots, natDegree_X_pow_sub_C]
    exact root.card_nthRoots_one
  have Q_splits : Splits (RingHom.id K) Q := splits_of_splits_of_dvd _ big_ne big_splits dvd
  have cardQ : Multiset.card Q.roots = 2 := by
    rw [splits_iff_card_roots.mp Q_splits, degQ]
  have le : Q.roots ≤ nthRoots N (1 : K) := roots.le_of_dvd big_ne dvd
  have nodup : Q.roots.Nodup := Multiset.nodup_of_le le root.nthRoots_one_nodup
  obtain ⟨β₁, β₂, h12⟩ := Multiset.card_eq_two.mp cardQ
  have m1 : β₁ ∈ Q.roots := by rw [h12]; simp
  have m2 : β₂ ∈ Q.roots := by rw [h12]; simp
  have ne : β₁ ≠ β₂ := by
    intro h
    rw [h12, h] at nodup
    simp at nodup
  have ev : ∀ β ∈ Q.roots, β^N = 1 ∧ β^2 - s*β + n = 0 := by
    intro β hβ
    refine ⟨(mem_nthRoots hN).mp (Multiset.mem_of_le le hβ), ?_⟩
    have h := (mem_roots Qne).mp hβ
    rw [hQ] at h
    simpa [IsRoot.def] using h
  exact ⟨β₁, β₂, ne, (ev _ m1).1, (ev _ m2).1, (ev _ m1).2, (ev _ m2).2⟩

end Polynomial

section Rows

variable {R : Type*} [CommRing R] [Algebra K R]

local notation "ι" => algebraMap K R

/-- The recursion of Lemma 2.3: all remainder rows up to level `k`. -/
theorem rows (Cl Low : Submodule K R) (E : R) (w : ℕ → R) (k : ℕ) (s n : K) (hk : 0 < k)
    (lowE : E ∈ Low) (low : ∀ j, E*w j ∈ Low)
    (mul : ∀ j < k, ∀ y ∈ Cl, y ∈ Low → w j*y ∈ Cl)
    (sq : ∀ j < k, E*(w j^2 - w (j+1)) ∈ Cl)
    (quad : E*(w 1 - ι s*w 0 + ι n) ∈ Cl) :
    ∀ j ≤ k, E*(w j - ι (rem s n j).1 - ι (rem s n j).2*w 0) ∈ Cl ∧
      E*(w 0*w j - ι (-n*(rem s n j).2) - ι ((rem s n j).1 + s*(rem s n j).2)*w 0) ∈ Cl := by
  have scale : ∀ (c : K) (y : R), y ∈ Cl → ι c*y ∈ Cl := by
    intro c y hy
    rw [← Algebra.smul_def]
    exact Cl.smul_mem c hy
  have scaleLow : ∀ (c : K) (y : R), y ∈ Low → ι c*y ∈ Low := by
    intro c y hy
    rw [← Algebra.smul_def]
    exact Low.smul_mem c hy
  have lowrow : ∀ (j : ℕ) (a b : K), E*(w j - ι a - ι b*w 0) ∈ Low := by
    intro j a b
    have h := Low.sub_mem (Low.sub_mem (low j) (scaleLow a E lowE)) (scaleLow b _ (low 0))
    have e : E*(w j - ι a - ι b*w 0) = E*w j - ι a*E - ι b*(E*w 0) := by ring
    rw [e]
    exact h
  have quad' : E*(w 0^2 - ι s*w 0 + ι n) ∈ Cl := by
    have h := Cl.add_mem (sq 0 hk) quad
    have e : E*(w 0^2 - ι s*w 0 + ι n) = E*(w 0^2 - w (0+1)) + E*(w 1 - ι s*w 0 + ι n) := by
      simp only [zero_add]
      ring
    rw [e]
    exact h
  have second : ∀ (j : ℕ) (a b : K), E*(w j - ι a - ι b*w 0) ∈ Cl →
      E*(w 0*w j - ι (-n*b) - ι (a + s*b)*w 0) ∈ Cl := by
    intro j a b hr
    have h := Cl.add_mem (mul 0 hk _ hr (lowrow j a b)) (scale b _ quad')
    have e : E*(w 0*w j - ι (-n*b) - ι (a + s*b)*w 0) =
        w 0*(E*(w j - ι a - ι b*w 0)) + ι b*(E*(w 0^2 - ι s*w 0 + ι n)) := by
      simp only [map_neg, map_mul, map_add]
      ring
    rw [e]
    exact h
  intro j
  induction j with
  | zero =>
      intro _
      have first : E*(w 0 - ι (rem s n 0).1 - ι (rem s n 0).2*w 0) ∈ Cl := by
        have e : E*(w 0 - ι (rem s n 0).1 - ι (rem s n 0).2*w 0) = 0 := by simp [rem]
        rw [e]
        exact Cl.zero_mem
      exact ⟨first, second 0 _ _ first⟩
  | succ j ih =>
      intro hj
      obtain ⟨hr, hr'⟩ := ih (by omega)
      have h := Cl.sub_mem (Cl.add_mem (Cl.add_mem (mul j (by omega) _ hr (lowrow j _ _))
        (scale (rem s n j).1 _ hr)) (scale (rem s n j).2 _ hr')) (sq j (by omega))
      have e : E*(w (j+1) - ι (rem s n (j+1)).1 - ι (rem s n (j+1)).2*w 0) =
          w j*(E*(w j - ι (rem s n j).1 - ι (rem s n j).2*w 0)) +
            ι (rem s n j).1*(E*(w j - ι (rem s n j).1 - ι (rem s n j).2*w 0)) +
            ι (rem s n j).2*(E*(w 0*w j - ι (-n*(rem s n j).2) -
              ι ((rem s n j).1 + s*(rem s n j).2)*w 0)) - E*(w j^2 - w (j+1)) := by
        simp only [rem, map_neg, map_mul, map_add, map_sub, map_pow, map_ofNat]
        ring
      have first : E*(w (j+1) - ι (rem s n (j+1)).1 - ι (rem s n (j+1)).2*w 0) ∈ Cl := by
        rw [e]
        exact h
      exact ⟨first, second (j+1) _ _ first⟩

/-- The end of Lemma 2.3: the two top rows give every linear row vanishing on the roots of
`Q` among the `N`-th roots of unity. -/
theorem linear_mem (Cl : Submodule K R) (E x : R) {N : ℕ} (s n a b : K)
    (root_val : ∀ θ : K, θ^2 - s*θ + n = 0 → a + b*θ = θ^N)
    (unit_case : a = 1 → b = 0 → ∃ β₁ β₂ : K, β₁ ≠ β₂ ∧ β₁^N = 1 ∧ β₂^N = 1 ∧
      β₁^2 - s*β₁ + n = 0 ∧ β₂^2 - s*β₂ + n = 0)
    (row₁ : E*(1 - ι a - ι b*x) ∈ Cl)
    (row₂ : E*(x - ι (-n*b) - ι (a + s*b)*x) ∈ Cl)
    (c₀ c₁ : K) (vanish : ∀ β : K, β^N = 1 → β^2 - s*β + n = 0 → c₀ + c₁*β = 0) :
    E*(ι c₀ + ι c₁*x) ∈ Cl := by
  have scale : ∀ (c : K) (y : R), y ∈ Cl → ι c*y ∈ Cl := by
    intro c y hy
    rw [← Algebra.smul_def]
    exact Cl.smul_mem c hy
  by_cases hdet : (1-a)*(1-a-s*b) + n*b^2 = 0
  · by_cases hb : b = 0
    · have ha : a = 1 := by
        rw [hb] at hdet
        have square : (1-a)^2 = 0 := by linear_combination hdet
        have h := (pow_eq_zero_iff (by norm_num : (2 : ℕ) ≠ 0)).mp square
        linear_combination -h
      obtain ⟨β₁, β₂, ne, h1, h2, q1, q2⟩ := unit_case ha hb
      have e1 := vanish β₁ h1 q1
      have e2 := vanish β₂ h2 q2
      have product : c₁*(β₁ - β₂) = 0 := by linear_combination e1 - e2
      have hc₁ : c₁ = 0 := by
        rcases mul_eq_zero.mp product with h | h
        · exact h
        · exact absurd (sub_eq_zero.mp h) ne
      have hc₀ : c₀ = 0 := by
        rw [hc₁] at e1
        linear_combination e1
      rw [hc₀, hc₁]
      simp
    · obtain ⟨θ, hθ⟩ : ∃ θ : K, b*θ = 1-a := ⟨(1-a)/b, mul_div_cancel₀ _ hb⟩
      have scaled : b^2*(θ^2 - s*θ + n) = 0 := by
        linear_combination hdet + (b*θ + (1-a) - s*b)*hθ
      have hQ : θ^2 - s*θ + n = 0 := by
        rcases mul_eq_zero.mp scaled with h | h
        · exact absurd h (pow_ne_zero 2 hb)
        · exact h
      have hθN : θ^N = 1 := by
        rw [← root_val θ hQ]
        linear_combination hθ
      have hv := vanish θ hθN hQ
      obtain ⟨d, hd⟩ : ∃ d : K, d*b = -c₁ := ⟨-c₁/b, div_mul_cancel₀ _ hb⟩
      have k1 : c₁ = d*(-b) := by linear_combination hd
      have k0 : c₀ = d*(1-a) := by
        linear_combination hv - θ*hd + d*hθ
      have e : E*(ι (d*(1-a)) + ι (d*(-b))*x) = ι d*(E*(1 - ι a - ι b*x)) := by
        simp only [map_mul, map_sub, map_neg, map_one]
        ring
      rw [k0, k1, e]
      exact scale d _ row₁
  · obtain ⟨d₀, hd₀⟩ : ∃ d₀ : K, c₀ = d₀*((1-a)*(1-a-s*b) + n*b^2) :=
      ⟨c₀/((1-a)*(1-a-s*b) + n*b^2), (div_mul_cancel₀ _ hdet).symm⟩
    obtain ⟨d₁, hd₁⟩ : ∃ d₁ : K, c₁ = d₁*((1-a)*(1-a-s*b) + n*b^2) :=
      ⟨c₁/((1-a)*(1-a-s*b) + n*b^2), (div_mul_cancel₀ _ hdet).symm⟩
    have h := Cl.add_mem
      (scale d₀ _ (Cl.add_mem (scale (1-a-s*b) _ row₁) (scale b _ row₂)))
      (scale d₁ _ (Cl.add_mem (scale (-n*b) _ row₁) (scale (1-a) _ row₂)))
    have e : E*(ι (d₀*((1-a)*(1-a-s*b) + n*b^2)) + ι (d₁*((1-a)*(1-a-s*b) + n*b^2))*x) =
        ι d₀*(ι (1-a-s*b)*(E*(1 - ι a - ι b*x)) + ι b*(E*(x - ι (-n*b) - ι (a + s*b)*x))) +
        ι d₁*(ι (-n*b)*(E*(1 - ι a - ι b*x)) + ι (1-a)*(E*(x - ι (-n*b) - ι (a + s*b)*x))) := by
      simp only [map_mul, map_add, map_sub, map_neg, map_pow, map_one]
      ring
    rw [hd₀, hd₁, e]
    exact h

/-- **Univariate lemma.** With `w k = 1` and the `2^k`-th roots of unity in the field: for
every `j ≤ k`, the rows `E * (w j - c₀ - c₁ * w 0)` and `E * (w 0 * w j - c₀ - c₁ * w 0)` are
in `Cl` as soon as `c₀ + c₁ X` takes the right values on the roots of `Q` among the roots of
unity. -/
theorem fibre_rows (Cl Low : Submodule K R) (E : R) (w : ℕ → R) (k : ℕ) (s n : K) (hk : 0 < k)
    {ζ : K} (root : IsPrimitiveRoot ζ (2^k))
    (lowE : E ∈ Low) (low : ∀ j, E*w j ∈ Low)
    (mul : ∀ j < k, ∀ y ∈ Cl, y ∈ Low → w j*y ∈ Cl)
    (sq : ∀ j < k, E*(w j^2 - w (j+1)) ∈ Cl)
    (quad : E*(w 1 - ι s*w 0 + ι n) ∈ Cl) (top : w k = 1) :
    ∀ j ≤ k, ∀ c₀ c₁ : K,
      ((∀ β : K, β^2^k = 1 → β^2 - s*β + n = 0 → c₀ + c₁*β = β^2^j) →
        E*(w j - ι c₀ - ι c₁*w 0) ∈ Cl) ∧
      ((∀ β : K, β^2^k = 1 → β^2 - s*β + n = 0 → c₀ + c₁*β = β^(2^j+1)) →
        E*(w 0*w j - ι c₀ - ι c₁*w 0) ∈ Cl) := by
  have all := rows Cl Low E w k s n hk lowE low mul sq quad
  obtain ⟨top₁, top₂⟩ := all k le_rfl
  rw [top] at top₁ top₂
  rw [mul_one] at top₂
  have linear : ∀ c₀ c₁ : K, (∀ β : K, β^2^k = 1 → β^2 - s*β + n = 0 → c₀ + c₁*β = 0) →
      E*(ι c₀ + ι c₁*w 0) ∈ Cl := by
    intro c₀ c₁ vanish
    refine linear_mem Cl E (w 0) s n (rem s n k).1 (rem s n k).2
      (fun θ hθ => rem_root s n θ hθ k) ?_ top₁ top₂ c₀ c₁ vanish
    intro ha hb
    have dvd := rem_dvd s n k
    rw [ha, hb] at dvd
    have dvd' : (Polynomial.X^2 - Polynomial.C s*Polynomial.X + Polynomial.C n : Polynomial K) ∣
        Polynomial.X^2^k - Polynomial.C 1 := by
      simpa using dvd
    exact two_roots root (by positivity) s n dvd'
  intro j hj c₀ c₁
  obtain ⟨row₁, row₂⟩ := all j hj
  constructor
  · intro values
    have correction := linear ((rem s n j).1 - c₀) ((rem s n j).2 - c₁) (by
      intro β hβ hQ
      have h1 := rem_root s n β hQ j
      have h2 := values β hβ hQ
      linear_combination h1 - h2)
    have e : E*(w j - ι c₀ - ι c₁*w 0) =
        E*(w j - ι (rem s n j).1 - ι (rem s n j).2*w 0) +
          E*(ι ((rem s n j).1 - c₀) + ι ((rem s n j).2 - c₁)*w 0) := by
      simp only [map_sub]
      ring
    rw [e]
    exact Cl.add_mem row₁ correction
  · intro values
    have correction := linear (-n*(rem s n j).2 - c₀) ((rem s n j).1 + s*(rem s n j).2 - c₁) (by
      intro β hβ hQ
      have h1 := rem_root s n β hQ j
      have h2 := values β hβ hQ
      rw [pow_succ] at h2
      linear_combination β*h1 - h2 - (rem s n j).2*hQ)
    have e : E*(w 0*w j - ι c₀ - ι c₁*w 0) =
        E*(w 0*w j - ι (-n*(rem s n j).2) - ι ((rem s n j).1 + s*(rem s n j).2)*w 0) +
          E*(ι (-n*(rem s n j).2 - c₀) + ι ((rem s n j).1 + s*(rem s n j).2 - c₁)*w 0) := by
      simp only [map_sub]
      ring
    rw [e]
    exact Cl.add_mem row₂ correction

end Rows

end ICFallDegrees.SelectorUnivariate
