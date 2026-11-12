import ICFallDegrees.CyclicNAF
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Data.Int.ModEq
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.RingTheory.Coprime.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith

/-! Construct bounded periodic carries from modular equality. Unlike the
previous rigidity theorem, the public conclusions assume no carry witness. -/
namespace ICFallDegrees.CyclicCarryWindow
open scoped BigOperators

def window (a : ℕ → ℤ) (i n : ℕ) : ℤ :=
  ∑ j ∈ Finset.range n, (2 : ℤ)^j*a (i+j)

private theorem powers (n : ℕ) : (∑ j ∈ Finset.range n, (2 : ℤ)^j)=2^n-1 := by
  induction n with
  | zero => simp
  | succ n ih => rw [Finset.sum_range_succ,ih,pow_succ]; ring

private theorem window_sub (a b : ℕ → ℤ) (i n : ℕ) :
    window (fun j => a j-b j) i n=window a i n-window b i n := by
  simp [window,mul_sub,Finset.sum_sub_distrib]

private theorem window_period (a : ℕ → ℤ) (n : ℕ)
    (period : ∀ i, a (i+n)=a i) (i : ℕ) : window a (i+n) n=window a i n := by
  apply Finset.sum_congr rfl
  intro j _
  have index : i+n+j=i+j+n := by omega
  rw [index,period]

private theorem window_shift (a : ℕ → ℤ) (i n : ℕ)
    (period : a (i+n)=a i) :
    2*window a (i+1) n=window a i n+(2^n-1)*a i := by
  have first : window a i (n+1)=window a i n+2^n*a (i+n) := by
    exact Finset.sum_range_succ _ _
  have second : window a i (n+1)=2*window a (i+1) n+a i := by
    unfold window
    rw [Finset.sum_range_succ']
    simp only [pow_zero,mul_one,one_mul,Nat.add_zero]
    rw [Finset.mul_sum]
    congr 1
    apply Finset.sum_congr rfl
    intro j _
    have index : i+(j+1)=i+1+j := by omega
    rw [index,pow_succ]
    ring
  rw [period] at first
  linarith

private theorem has_zero (a : ℕ → ℤ) (n : ℕ)
    (period : ∀ i, a (i+(n+1))=a i)
    (nonadjacent : ∀ i, a i≠0 → a (i+1)=0) (i : ℕ) :
    ∃ j<n+1, a (i+j)=0 := by
  by_cases zero : a i=0
  · exact ⟨0,by omega,by simpa using zero⟩
  · have next := nonadjacent i zero
    by_cases hn : n=0
    · subst n
      have same := period i
      simp only [Nat.zero_add] at same
      exact False.elim (zero (same.symm.trans next))
    · exact ⟨1,by omega,next⟩

private theorem window_upper (a : ℕ → ℤ) (n i : ℕ)
    (bound : ∀ j, a j≤1) (zero : ∃ j<n, a (i+j)=0) :
    window a i n<2^n-1 := by
  have compare : ∀ j∈Finset.range n, (2 : ℤ)^j*a (i+j)≤2^j := by
    intro j _
    simpa using mul_le_mul_of_nonneg_left (bound (i+j)) (le_of_lt (pow_pos (by decide : (0 : ℤ)<2) j))
  have strict : ∃ j∈Finset.range n, (2 : ℤ)^j*a (i+j)<2^j := by
    obtain ⟨j,hj,hz⟩ := zero
    exact ⟨j,Finset.mem_range.mpr hj,by simp [hz,pow_pos (by decide : (0 : ℤ)<2) j]⟩
  exact (Finset.sum_lt_sum compare strict).trans_eq (powers n)

private theorem window_bounds (a : ℕ → ℤ) (n i : ℕ)
    (period : ∀ i, a (i+(n+1))=a i)
    (bound : ∀ j, -1≤a j ∧ a j≤1)
    (nonadjacent : ∀ i, a i≠0 → a (i+1)=0) :
    -(2^(n+1)-1)<window a i (n+1) ∧ window a i (n+1)<2^(n+1)-1 := by
  have zero := has_zero a n period nonadjacent i
  have upper := window_upper a (n+1) i (fun j => (bound j).2) zero
  have negative : ∀ j, -a j≤1 := by intro j; have h := (bound j).1; omega
  have negzero : ∃ j<n+1, -a (i+j)=0 := by
    obtain ⟨j,hj,hz⟩ := zero
    exact ⟨j,hj,by rw [hz,neg_zero]⟩
  have lower := window_upper (fun j => -a j) (n+1) i negative negzero
  have negwindow : window (fun j => -a j) i (n+1)= -window a i (n+1) := by
    simp [window,mul_neg,Finset.sum_neg_distrib]
  rw [negwindow] at lower
  exact ⟨by linarith,upper⟩

/-- Modular equality supplies a periodic integer carry sequence in [-1,1]. -/
theorem exists_bounded_carries (a b : ℕ → ℤ) (n : ℕ)
    (pa : ∀ i, a (i+(n+1))=a i) (pb : ∀ i, b (i+(n+1))=b i)
    (ha : ∀ i, -1≤a i ∧ a i≤1) (hb : ∀ i, -1≤b i ∧ b i≤1)
    (na : ∀ i, a i≠0 → a (i+1)=0) (nb : ∀ i, b i≠0 → b (i+1)=0)
    (same : (2^(n+1)-1 : ℤ) ∣ window a 0 (n+1)-window b 0 (n+1)) :
    ∃ c : ℕ → ℤ, (∀ i, c (i+(n+1))=c i) ∧
      (∀ i, -1≤c i ∧ c i≤1) ∧ (∀ i, 2*c (i+1)=c i+a i-b i) := by
  let M : ℤ := 2^(n+1)-1
  let d : ℕ → ℤ := fun i => a i-b i
  have positive : 0<M := by
    dsimp [M]
    have h : (1 : ℤ)≤2^n := one_le_pow₀ (by decide)
    rw [pow_succ]
    omega
  have coprime : IsCoprime M (2 : ℤ) := by
    refine ⟨-1,2^n,?_⟩
    dsimp [M]
    rw [pow_succ]
    ring
  have pd : ∀ i, d (i+(n+1))=d i := by intro i; dsimp [d]; rw [pa,pb]
  have divisible : ∀ i, M∣window d i (n+1) := by
    intro i
    induction i with
    | zero => simpa [M,window_sub,d] using same
    | succ i ih =>
        apply coprime.dvd_of_dvd_mul_left
        rw [window_shift d i (n+1) (pd i)]
        exact dvd_add ih (dvd_mul_right M (d i))
  let c : ℕ → ℤ := fun i => window d i (n+1)/M
  have image : ∀ i, M*c i=window d i (n+1) := fun i => Int.mul_ediv_cancel' (divisible i)
  refine ⟨c,?_,?_,?_⟩
  · intro i
    dsimp [c]
    rw [window_period d (n+1) pd i]
  · intro i
    have ab := window_bounds a n i pa ha na
    have bb := window_bounds b n i pb hb nb
    have value := image i
    rw [window_sub] at value
    change M*c i=window a i (n+1)-window b i (n+1) at value
    have low : (-2 : ℤ)<c i := (mul_lt_mul_left positive).mp (by dsimp [M] at *; linarith)
    have high : c i<2 := (mul_lt_mul_left positive).mp (by dsimp [M] at *; linarith)
    constructor <;> omega
  · intro i
    apply mul_left_cancel₀ (ne_of_gt positive)
    calc
      M*(2*c (i+1))=2*(M*c (i+1)) := by ring
      _=2*window d (i+1) (n+1) := by rw [image]
      _=window d i (n+1)+M*d i := window_shift d i (n+1) (pd i)
      _=M*(c i+a i-b i) := by rw [←image i]; dsimp [d]; ring

/-- The only distinct congruent cyclic nonadjacent words have alternating
opposite signs. The carry witness is constructed, not assumed. -/
theorem rigidity_of_modEq (a b : ℕ → ℤ) (n : ℕ)
    (pa : ∀ i, a (i+(n+1))=a i) (pb : ∀ i, b (i+(n+1))=b i)
    (ha : ∀ i, -1≤a i ∧ a i≤1) (hb : ∀ i, -1≤b i ∧ b i≤1)
    (na : ∀ i, a i≠0 → a (i+1)=0) (nb : ∀ i, b i≠0 → b (i+1)=0)
    (same : window a 0 (n+1) ≡ window b 0 (n+1) [ZMOD (2^(n+1)-1)]) :
    (∀ i, a i=b i) ∨ ∃ epsilon : ℤ, (epsilon=1 ∨ epsilon= -1) ∧
      ∀ i, a i-b i=epsilon ∧ a (i+1)=epsilon-a i := by
  have divides : (2^(n+1)-1 : ℤ)∣window a 0 (n+1)-window b 0 (n+1) := by
    simpa only [neg_sub] using (dvd_neg.mpr same.dvd)
  obtain ⟨c,period,bounds,equation⟩ := exists_bounded_carries a b n pa pb ha hb na nb divides
  rcases CyclicNAF.digit_rigidity a b c (n+1) (by omega) period ha hb bounds na nb equation with equal | other
  · exact Or.inl equal
  · exact Or.inr ⟨c 0,other⟩

/-- Actual odd cyclic NAF uniqueness from equality modulo 2^(2k+1)-1;
no carry witness or uniqueness hypothesis is assumed. -/
theorem odd_unique_of_modEq (a b : ℕ → ℤ) (k : ℕ)
    (pa : ∀ i, a (i+(2*k+1))=a i) (pb : ∀ i, b (i+(2*k+1))=b i)
    (ha : ∀ i, -1≤a i ∧ a i≤1) (hb : ∀ i, -1≤b i ∧ b i≤1)
    (na : ∀ i, a i≠0 → a (i+1)=0) (nb : ∀ i, b i≠0 → b (i+1)=0)
    (same : window a 0 (2*k+1) ≡ window b 0 (2*k+1) [ZMOD (2^(2*k+1)-1)]) :
    ∀ i, a i=b i := by
  have divides : (2^(2*k+1)-1 : ℤ)∣window a 0 (2*k+1)-window b 0 (2*k+1) := by
    simpa only [neg_sub] using (dvd_neg.mpr same.dvd)
  obtain ⟨c,period,bounds,equation⟩ := exists_bounded_carries a b (2*k) pa pb ha hb na nb divides
  exact CyclicNAF.odd_uniqueness a b c k pa period ha hb bounds na nb equation

end ICFallDegrees.CyclicCarryWindow
