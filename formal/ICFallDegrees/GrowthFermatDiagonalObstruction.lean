import ICFallDegrees.GrowthFixedE6Existence
import ICFallDegrees.GrowthExistence
import Mathlib.Tactic

/-!
# N8: the fixed E6 system is not the unit ideal over Fermat-factor primes

The prime-field half of the fixed-curve growth theorem
(`GrowthFixedE6PrimeTheta`, `GrowthFixedE6PrimeSequence`) is stated for
primes `p ≡ 1 (mod 2^(k+1))` beyond a bound that depends on the depth `k`
and is not computed. This file shows that the unit-ideal hypothesis
is not automatic for admissible primes: for every depth `k ≥ 2`, every prime factor `p` of the Fermat
number `2^(2^k) + 1` satisfies `p ≡ 1 (mod 2^(k+1))`, and over `F_p` the
system of `y² = x³ + 6x - 7` with target abscissa `1` has the common zero
`x₁ = x₂ = 4`, so it is a proper ideal.

`e6_diagonal_grid_root` is the statement of ml-cryptanalysis `0712c0d5`:
the point `(4, 4)` lies on the grid and on the eliminant. The remaining
statements, added in the companion repository on 2026-11-14, are about the
polynomial system itself (`PrimeM2General.inputs`).
-/

noncomputable section
namespace ICFallDegrees.GrowthFermatDiagonalObstruction

open MvPolynomial
open ICFallDegrees.PrimeM2General
open ICFallDegrees.GeneralQuadraticTop (top)

/-- If 2^(2^k) = -1, then 4 is a 2^k-th root of unity and the
fixed E6 eliminant h(x,z)=(x-1)z-x-8 vanishes at (4,4).
The equality h(4,4)=0 holds already over the integers. -/
theorem e6_diagonal_grid_root (R : Type*) [CommRing R] (k : ℕ)
    (hF : (2 : R) ^ (2 ^ k) + 1 = 0) :
    (4 : R) ^ (2 ^ k) = 1 ∧
      (((4 : R) - 1) * 4 - 4 - 8 = 0) := by
  have hpow : (2 : R) ^ (2 ^ k) = -1 := by
    linear_combination hF
  constructor
  · calc
      (4 : R) ^ (2 ^ k) = ((2 : R) ^ 2) ^ (2 ^ k) := by norm_num
      _ = (2 : R) ^ (2 * (2 ^ k)) := by rw [← pow_mul]
      _ = (2 : R) ^ ((2 ^ k) * 2) := by rw [Nat.mul_comm]
      _ = ((2 : R) ^ (2 ^ k)) ^ 2 := by rw [pow_mul]
      _ = 1 := by rw [hpow]; ring
  · ring

variable (K : Type*) [Field K]

/-- The diagonal point of the system at depth `k`: the chain variable of level
`j` of either summand takes the value `4^(2^j)`. -/
def point (k : ℕ) (i : Fin (2 * k)) : K :=
  if i.val < k then (4 : K) ^ (2 ^ i.val) else (4 : K) ^ (2 ^ (i.val - k))

/-- Every chain equation vanishes at the diagonal point. -/
theorem eval_chain (k : ℕ) (h4 : (4 : K) ^ (2 ^ k) = 1) (i : Fin (2 * k)) :
    eval (point K k) (X i ^ 2 + chainTail k i : P k K) = 0 := by
  have hi := i.isLt
  unfold chainTail
  split_ifs with hend
  · have hsq : point K k i ^ 2 = 1 := by
      unfold point
      rcases hend with h | h
      · have hlt : i.val < k := by omega
        have hexp : 2 ^ i.val * 2 = 2 ^ k := by
          rw [← pow_succ, h]
        rw [if_pos hlt, ← pow_mul, hexp, h4]
      · have hge : ¬ i.val < k := by omega
        have hexp : 2 ^ (i.val - k) * 2 = 2 ^ k := by
          rw [← pow_succ]
          congr 1
          omega
        rw [if_neg hge, ← pow_mul, hexp, h4]
    simp only [map_add, map_pow, eval_X, map_neg, map_one, hsq]
    ring
  · have hnext : point K k i ^ 2 = point K k ⟨i.val + 1, by omega⟩ := by
      unfold point
      by_cases hlt : i.val < k
      · have hlt' : i.val + 1 < k := by omega
        rw [if_pos hlt, if_pos hlt', ← pow_mul, ← pow_succ]
      · have hge' : ¬ i.val + 1 < k := by omega
        have hexp : i.val + 1 - k = i.val - k + 1 := by omega
        rw [if_neg hlt, if_neg hge', ← pow_mul, ← pow_succ, hexp]
    simp only [map_add, map_pow, eval_X, map_neg, hnext]
    ring

/-- The normalized summation polynomial of `y² = x³ + 6x - 7` with target
abscissa `1` vanishes at the diagonal point. This needs no hypothesis on the
field: it is the integer identity `(4 - 1)·4 - 4 - 8 = 0`, squared. -/
theorem eval_normalized (k : ℕ) (depth : 2 ≤ k) :
    eval (point K k) (normalized k depth (6 : K) (-7 : K) 1) = 0 := by
  have h0 : point K k (active k depth 0) = 4 := by
    have hk : 0 < k := by omega
    simp [point, active, hk]
  have h1 : point K k (active k depth 1) = 16 := by
    have hk : 1 < k := by omega
    simp [point, active, hk]
    norm_num
  have h2 : point K k (active k depth 2) = 4 := by
    simp [point, active]
  have h3 : point K k (active k depth 3) = 16 := by
    simp [point, active]
    norm_num
  simp only [normalized, top, tail, map_add, map_mul, eval_C, eval_X, h0, h1, h2, h3]
  norm_num

/-- **The system is a proper ideal whenever `2^(2^k) = -1` in the field.** -/
theorem e6_not_unit_of_fermat (k : ℕ) (depth : 2 ≤ k)
    (hF : (2 : K) ^ (2 ^ k) + 1 = 0) :
    (1 : P k K) ∉
      Ideal.span (Set.range (inputs k depth (6 : K) (-7 : K) 1)) := by
  have h4 := (e6_diagonal_grid_root K k hF).1
  intro hone
  have hle : Ideal.span (Set.range (inputs k depth (6 : K) (-7 : K) 1)) ≤
      RingHom.ker (eval (point K k)) := by
    rw [Ideal.span_le]
    rintro _ ⟨j, rfl⟩
    rw [SetLike.mem_coe, RingHom.mem_ker]
    induction j using Fin.lastCases with
    | last => simpa [inputs] using eval_normalized K k depth
    | cast i => simpa [inputs] using eval_chain K k h4 i
  have hzero := hle hone
  rw [RingHom.mem_ker, map_one] at hzero
  exact one_ne_zero hzero

/-- Prime factors of the Fermat number `2^(2^k) + 1`: the prime is
`≡ 1 (mod 2^(k+1))`, exceeds `2^(k+1)`, and the fixed E6 system at depth `k`
over `F_p` is a proper ideal. -/
theorem e6_not_unit_of_fermat_factor (p : ℕ) [Fact p.Prime] (k : ℕ) (depth : 2 ≤ k)
    (h : p ∣ 2 ^ 2 ^ k + 1) :
    p ≡ 1 [MOD 2 ^ (k + 1)] ∧ 2 ^ (k + 1) < p ∧
      (1 : MvPolynomial (Fin (2 * k)) (ZMod p)) ∉
        Ideal.span (Set.range (inputs k depth (6 : ZMod p) (-7 : ZMod p) 1)) := by
  have prime : p.Prime := Fact.out
  have hdvd := GrowthExistence.fermat_factor_dvd p k h
  have hp1 : 1 ≤ p := prime.one_lt.le
  have hpos : 0 < p - 1 := by have := prime.two_le; omega
  have hz : ((2 ^ 2 ^ k + 1 : ℕ) : ZMod p) = 0 :=
    (ZMod.natCast_zmod_eq_zero_iff_dvd _ p).mpr h
  have hF : (2 : ZMod p) ^ (2 ^ k) + 1 = 0 := by
    push_cast at hz
    exact hz
  refine ⟨?_, ?_, e6_not_unit_of_fermat (ZMod p) k depth hF⟩
  · exact ((Nat.modEq_iff_dvd' hp1).mpr hdvd).symm
  · have := Nat.le_of_dvd hpos hdvd
    omega

/-- **Bad admissible primes exist at every depth.** For every `k ≥ 2` there is
a prime `p ≡ 1 (mod 2^(k+1))`, `p > 2^(k+1)`, over which the fixed E6 system
at depth `k` is a proper ideal. At these primes the unit-ideal argument
cannot establish the numerical LFD lower bound. Whether that numerical
bound holds for these proper ideals is not decided by this theorem. -/
theorem exists_bad_admissible_prime (k : ℕ) (depth : 2 ≤ k) :
    ∃ p : ℕ, Nat.Prime p ∧ p ≡ 1 [MOD 2 ^ (k + 1)] ∧ 2 ^ (k + 1) < p ∧
      ∀ [Fact p.Prime],
        (1 : MvPolynomial (Fin (2 * k)) (ZMod p)) ∉
          Ideal.span (Set.range (inputs k depth (6 : ZMod p) (-7 : ZMod p) 1)) := by
  have hne : 2 ^ 2 ^ k + 1 ≠ 1 := by
    have : 0 < 2 ^ 2 ^ k := by positivity
    omega
  obtain ⟨p, hp, hdiv⟩ := Nat.exists_prime_and_dvd hne
  haveI : Fact p.Prime := ⟨hp⟩
  obtain ⟨hmod, hgt, hnot⟩ := e6_not_unit_of_fermat_factor p k depth hdiv
  exact ⟨p, hp, hmod, hgt, fun {_} => hnot⟩

#print axioms e6_diagonal_grid_root
#print axioms e6_not_unit_of_fermat
#print axioms e6_not_unit_of_fermat_factor
#print axioms exists_bad_admissible_prime

end ICFallDegrees.GrowthFermatDiagonalObstruction
