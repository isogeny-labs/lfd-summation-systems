import Mathlib.Algebra.Module.LinearMap.Basic
import Mathlib.Tactic

/-! Concrete cyclic coefficient operators. These identities do not prove
the analytic wrapping-part vanishing or a proper-ideal growth theorem. -/
namespace ICFallDegrees.GrowthJetClosure

variable {K A : Type*} [CommRing K]

abbrev Coefficients (A K : Type*) (N : ℕ) := A → Fin N → K

/-- Euler derivative in the distinguished exponent block. -/
def euler (N : ℕ) : Coefficients A K N →ₗ[K] Coefficients A K N where
  toFun f a b := (b.val : K) * f a b
  map_add' f g := by funext a b; simp [mul_add]
  map_smul' c f := by funext a b; simp [mul_left_comm]

/-- Cyclic monomial multiplication, expressed on coefficient arrays. -/
def cyclicShift (N s : ℕ) (hN : 0 < N) :
    Coefficients A K N →ₗ[K] Coefficients A K N where
  toFun f a b := f a ⟨(b.val + N - s) % N, Nat.mod_lt _ hN⟩
  map_add' _ _ := rfl
  map_smul' _ _ := rfl

/-- Source coefficients whose distinguished exponent wraps on shifting. -/
def wrappingPart (N s : ℕ) : Coefficients A K N →ₗ[K] Coefficients A K N where
  toFun f a b := if N - s ≤ b.val then f a b else 0
  map_add' f g := by
    funext a b
    change (if N - s ≤ b.val then (f + g) a b else 0) = _
    by_cases h : N - s ≤ b.val <;> simp only [h, if_true, if_false, Pi.add_apply, add_zero]
  map_smul' c f := by
    funext a b
    change (if N - s ≤ b.val then (c • f) a b else 0) = _
    by_cases h : N - s ≤ b.val <;> simp only [h, if_true, if_false, Pi.smul_apply, smul_zero, RingHom.id_apply]

theorem inverse_index (N s b : ℕ) (hN : 0 < N) (hs : s < N) (hb : b < N) :
    let q := (b + N - s) % N
    (b < s → q + s = b + N ∧ N - s ≤ q) ∧
    (¬ b < s → q + s = b ∧ ¬ N - s ≤ q) := by
  dsimp
  constructor
  · intro h
    have hlt : b + N - s < N := by omega
    rw [Nat.mod_eq_of_lt hlt]
    constructor <;> omega
  · intro h
    have heq : b + N - s = (b - s) + N := by omega
    have hm : (b + N - s) % N = b - s := by
      rw [heq, Nat.add_mod_right, Nat.mod_eq_of_lt (show b - s < N by omega)]
    rw [hm]
    constructor <;> omega

/-- Arbitrary-depth truncated-shift/Euler commutator. For chains take
N=2^k and s=2^j, j<k. No numerical-testing hypothesis is used. -/
theorem euler_cyclicShift (N s : ℕ) (hN : 0 < N) (hs : s < N)
    (f : Coefficients A K N) :
    euler N (cyclicShift N s hN f) =
      cyclicShift N s hN (euler N f) + s • cyclicShift N s hN f -
        N • cyclicShift N s hN (wrappingPart N s f) := by
  funext a b
  let q : Fin N := ⟨(b.val + N - s) % N, Nat.mod_lt _ hN⟩
  obtain ⟨hwrap, hnowrap⟩ := inverse_index N s b.val hN hs b.isLt
  by_cases hb : b.val < s
  · obtain ⟨heq, hw⟩ := hwrap hb
    have heq' : q.val + s = b.val + N := heq
    have hc : (q.val : K) + (s : K) = (b.val : K) + (N : K) := by
      simpa only [Nat.cast_add] using congrArg (Nat.cast (R := K)) heq'
    change (b.val : K) * f a q =
      (q.val : K) * f a q + s • f a q - N • (if N - s ≤ q.val then f a q else 0)
    rw [if_pos hw]
    simp only [nsmul_eq_mul]
    linear_combination -(congrArg (fun t : K => t * f a q) hc)
  · obtain ⟨heq, hw⟩ := hnowrap hb
    have heq' : q.val + s = b.val := heq
    have hc : (q.val : K) + (s : K) = (b.val : K) := by simpa only [Nat.cast_add] using congrArg (Nat.cast (R := K)) heq'
    change (b.val : K) * f a q =
      (q.val : K) * f a q + s • f a q - N • (if N - s ≤ q.val then f a q else 0)
    rw [if_neg hw]
    simp only [nsmul_eq_mul, mul_zero, sub_zero]
    rw [← add_mul, hc]

/-- A support budget, independent of the coefficient field. Taking the
weight to be the sum of two binary popcounts gives the chain degree. -/
def Supported (N : ℕ) (weight : A → Fin N → ℕ) (D : ℕ)
    (f : Coefficients A K N) : Prop :=
  ∀ a b, f a b ≠ 0 → weight a b ≤ D

/-- The concrete two-block chain specialization, with j<k and no
characteristic restriction. A is the first exponent block. -/
theorem chain_euler_shift (k j : ℕ) (hj : j < k)
    (f : Coefficients (Fin (2 ^ k)) K (2 ^ k)) :
    euler (2 ^ k) (cyclicShift (2 ^ k) (2 ^ j) (by positivity) f) =
      cyclicShift (2 ^ k) (2 ^ j) (by positivity) (euler (2 ^ k) f) +
      (2 ^ j) • cyclicShift (2 ^ k) (2 ^ j) (by positivity) f -
      (2 ^ k) • cyclicShift (2 ^ k) (2 ^ j) (by positivity)
        (wrappingPart (2 ^ k) (2 ^ j) f) :=
  euler_cyclicShift _ _ (by positivity) (Nat.pow_lt_pow_right (by decide) hj) f

theorem euler_supported (N : ℕ) (weight : A → Fin N → ℕ) (D : ℕ)
    (f : Coefficients A K N) (hf : Supported N weight D f) :
    Supported N weight D (euler N f) := by
  intro a b h
  apply hf a b
  intro hz
  apply h
  change (b.val : K) * f a b = 0
  rw [hz, mul_zero]

theorem wrappingPart_supported (N s : ℕ) (weight : A → Fin N → ℕ) (D : ℕ)
    (f : Coefficients A K N) (hf : Supported N weight D f) :
    Supported N weight D (wrappingPart N s f) := by
  intro a b h
  apply hf a b
  intro hz
  apply h
  change (if N - s ≤ b.val then f a b else 0) = 0
  rw [hz]
  split_ifs <;> rfl

/-- The first-jet space: both a row and its Euler derivative lie in S. -/
def jetSpace (N : ℕ) (S : Submodule K (Coefficients A K N)) :
    Submodule K (Coefficients A K N) := S ⊓ S.comap (euler N)

/-- A shift in another coordinate commutes with the Euler operator. -/
theorem euler_other_coordinate (N : ℕ) (pullback : A → A)
    (f : Coefficients A K N) :
    euler N (fun a b => f (pullback a) b) =
      fun a b => euler N f (pullback a) b := rfl

/-- Concrete commutator plus ordinary and wrapping closure imply first-jet
closure. The hypotheses are only imposed on low-weight rows. In the
intended application they still need to be derived from FEI. -/
theorem jetSpace_shift (N s : ℕ) (hN : 0 < N) (hs : s < N)
    (weight : A → Fin N → ℕ) (D : ℕ)
    (S : Submodule K (Coefficients A K N))
    (plain : ∀ f ∈ S, Supported N weight D f → cyclicShift N s hN f ∈ S)
    (wrap : ∀ f ∈ S, Supported N weight D f →
      cyclicShift N s hN (wrappingPart N s f) ∈ S)
    (f : Coefficients A K N) (hf : f ∈ jetSpace N S)
    (hdegree : Supported N weight D f) :
    cyclicShift N s hN f ∈ jetSpace N S := by
  obtain ⟨hplain, hderiv⟩ := hf
  constructor
  · exact plain f hplain hdegree
  · change euler N (cyclicShift N s hN f) ∈ S
    rw [euler_cyclicShift N s hN hs f]
    exact S.sub_mem
      (S.add_mem (plain _ hderiv (euler_supported N weight D f hdegree))
        (S.nsmul_mem (plain f hplain hdegree) s))
      (S.nsmul_mem (wrap f hplain hdegree) N)

/-- Any element with derivative outside S is excluded from the first-jet
space, even if the element itself lies in S. -/
theorem not_mem_jetSpace_of_derivative (N : ℕ)
    (S : Submodule K (Coefficients A K N)) (f : Coefficients A K N)
    (h : euler N f ∉ S) : f ∉ jetSpace N S := fun hf => h hf.2

#print axioms inverse_index
#print axioms euler_cyclicShift
#print axioms chain_euler_shift
#print axioms euler_supported
#print axioms wrappingPart_supported
#print axioms euler_other_coordinate
#print axioms jetSpace_shift
#print axioms not_mem_jetSpace_of_derivative
end ICFallDegrees.GrowthJetClosure
