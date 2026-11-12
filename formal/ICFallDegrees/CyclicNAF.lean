import Mathlib.Data.Int.Basic
import Lean.Elab.Tactic.Omega

/-! Rigidity of bounded cyclic carries between nonadjacent signed words.
The existence of these carries from modular equality is a separate written
bridge; no modular-arithmetic or polynomial hypothesis is hidden here. -/
namespace ICFallDegrees.CyclicNAF

theorem zero_forward (a b c : ℕ → ℤ)
    (ha : ∀ i, -1≤a i ∧ a i≤1) (hb : ∀ i, -1≤b i ∧ b i≤1)
    (hc : ∀ i, -1≤c i ∧ c i≤1)
    (na : ∀ i, a i≠0 → a (i+1)=0) (nb : ∀ i, b i≠0 → b (i+1)=0)
    (equation : ∀ i, 2*c (i+1)=c i+a i-b i) (i : ℕ) (zero : c i=0) :
    c (i+1)=0 := by
  by_contra nonzero
  have ai := ha i
  have bi := hb i
  have ci := hc (i+1)
  have first := equation i
  have an : a i≠0 := by omega
  have bn : b i≠0 := by omega
  have az := na i an
  have bz := nb i bn
  have second := equation (i+1)
  omega

theorem carry_constant (a b c : ℕ → ℤ) (n : ℕ) (positive : 0<n)
    (period : ∀ i, c (i+n)=c i)
    (ha : ∀ i, -1≤a i ∧ a i≤1) (hb : ∀ i, -1≤b i ∧ b i≤1)
    (hc : ∀ i, -1≤c i ∧ c i≤1)
    (na : ∀ i, a i≠0 → a (i+1)=0) (nb : ∀ i, b i≠0 → b (i+1)=0)
    (equation : ∀ i, 2*c (i+1)=c i+a i-b i) : ∀ i, c i=c 0 := by
  have iterate : ∀ i, c i=0 → ∀ k, c (i+k)=0 := by
    intro i hi k
    induction k with
    | zero => simpa using hi
    | succ k ih =>
        simpa [Nat.add_assoc] using zero_forward a b c ha hb hc na nb equation (i+k) ih
  have backwards : ∀ i, c (i+1)=0 → c i=0 := by
    intro i hi
    have h := iterate (i+1) hi (n-1)
    have index : i+1+(n-1)=i+n := by omega
    rw [index,period] at h
    exact h
  have step : ∀ i, c (i+1)=c i := by
    intro i
    by_cases zero : c i=0
    · simpa [zero] using zero_forward a b c ha hb hc na nb equation i zero
    · have next : c (i+1)≠0 := fun h => zero (backwards i h)
      have ai := ha i
      have bi := hb i
      have ci := hc i
      have cj := hc (i+1)
      have first := equation i
      omega
  intro i
  induction i with
  | zero => rfl
  | succ i ih => exact (step i).trans ih

theorem digit_rigidity (a b c : ℕ → ℤ) (n : ℕ) (positive : 0<n)
    (period : ∀ i, c (i+n)=c i)
    (ha : ∀ i, -1≤a i ∧ a i≤1) (hb : ∀ i, -1≤b i ∧ b i≤1)
    (hc : ∀ i, -1≤c i ∧ c i≤1)
    (na : ∀ i, a i≠0 → a (i+1)=0) (nb : ∀ i, b i≠0 → b (i+1)=0)
    (equation : ∀ i, 2*c (i+1)=c i+a i-b i) :
    (∀ i, a i=b i) ∨
    ((c 0=1 ∨ c 0= -1) ∧ ∀ i, a i-b i=c 0 ∧ a (i+1)=c 0-a i) := by
  have constant := carry_constant a b c n positive period ha hb hc na nb equation
  have difference : ∀ i, a i-b i=c 0 := by
    intro i
    have eq := equation i
    rw [constant i,constant (i+1)] at eq
    omega
  by_cases zero : c 0=0
  · left
    intro i
    have h := difference i
    omega
  · right
    have c0 := hc 0
    have signed : c 0=1 ∨ c 0= -1 := by omega
    refine ⟨signed,?_⟩
    intro i
    refine ⟨difference i,?_⟩
    have d := difference i
    have dn := difference (i+1)
    have ai := ha i
    have bi := hb i
    have aj := ha (i+1)
    have bj := hb (i+1)
    by_cases az : a i=0
    · have bn : b i≠0 := by omega
      have bz := nb i bn
      omega
    · have an := na i az
      omega

theorem odd_uniqueness (a b c : ℕ → ℤ) (k : ℕ)
    (period_a : ∀ i, a (i+(2*k+1))=a i)
    (period_c : ∀ i, c (i+(2*k+1))=c i)
    (ha : ∀ i, -1≤a i ∧ a i≤1) (hb : ∀ i, -1≤b i ∧ b i≤1)
    (hc : ∀ i, -1≤c i ∧ c i≤1)
    (na : ∀ i, a i≠0 → a (i+1)=0) (nb : ∀ i, b i≠0 → b (i+1)=0)
    (equation : ∀ i, 2*c (i+1)=c i+a i-b i) : ∀ i, a i=b i := by
  rcases digit_rigidity a b c (2*k+1) (by omega) period_c ha hb hc na nb equation with equal | exceptional
  · exact equal
  · obtain ⟨signed,flip⟩ := exceptional
    have twice : ∀ i, a (i+2)=a i := by
      intro i
      have h := (flip i).2
      have j : a (i+2)=c 0-a (i+1) := by
        simpa [Nat.add_assoc] using (flip (i+1)).2
      omega
    have even : ∀ j, a (2*j)=a 0 := by
      intro j
      induction j with
      | zero => rfl
      | succ j ih =>
          have h := twice (2*j)
          have index : 2*(j+1)=2*j+2 := by omega
          rw [index,h,ih]
    have last := (flip (2*k)).2
    have closed := period_a 0
    have same := even k
    simp only [Nat.zero_add] at closed
    omega

end ICFallDegrees.CyclicNAF
