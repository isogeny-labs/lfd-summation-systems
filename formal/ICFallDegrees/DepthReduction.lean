import ICFallDegrees.BinaryWeight
import Mathlib.LinearAlgebra.Finsupp.Defs
import Mathlib.Algebra.Module.Submodule.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Piecewise
import Mathlib.Data.Fintype.BigOperators

noncomputable section

/-!
# Depth reduction for the degree-`D` closure (the carry-combinatorics core of N5)

Modulo the chain equations, the `m`-summand presentation at depth `k` lives in the group
algebra `K[(Z/2^k)^m]`; a squarefree chain monomial is an exponent vector, its degree is the
binary weight, and the chain variable `u_{i,j}` is the translation by `2^j` in block `i`.
Here an exponent `a_i ∈ [0, 2^k)` is written as `a_i = l_i + 2^k0 q_i` with a *low part*
`l_i < 2^k0` and a *high pattern* `q_i < 2^(k-k0)`, so the depth-`k` algebra is the finitely
supported functions on `(Fin m → ℕ) × (Fin m → ℕ)` whose support lies in the box
`InBox k0 k`.  The weight is `wt l + wt q` (`BinaryWeight.wt_add_mul_two_pow`).  The
translation `translate k0 k i j` adds `2^j` to block `i` with carries, carrying out of the
low window into the high pattern and wrapping around at the top level `k-1`
(`step`, inverted on the box by `unstep`).

Given a subspace `S` of low elements (the depth-`k0` base, in the application `V_4(7)`), the
free high extension is

  `W k0 D S k = { w in the box : for every pattern q, the q-component of w lies in S
                 and every monomial l of it has wt l + wt q ≤ D }`,

that is `⊕_q τ_q (S ∩ M_{≤ D - wt q})` written through its components.  The theorems:

* `single_zero_mem_W_iff`: `1 ∈ W k0 D S k ↔ 1 ∈ S` (part (ii) of N5);
* `W_closed`: if `W k0 D S (k0+1)` is closed under all admissible variable translations at
  depth `k0+1`, then `W k0 D S k` is closed at every depth `k ≥ k0+1` (part (iii) of N5).

Nothing about `S` is assumed: not that it is a closure, not its dimension, not that it is
closed itself.  The instance-specific facts (`1 ∉ V_4(7)`, closedness of `W_8`) enter only
as hypotheses of the assembled statements in `DepthReductionBridge`.
-/

namespace ICFallDegrees.DepthReduction

open BinaryWeight

variable {m : ℕ} {K : Type*} [Field K]

/-! ## Exponents, weights, boxes -/

/-- Low elements: finitely supported functions on low exponent vectors. -/
abbrev Lo (m : ℕ) (K : Type*) [Field K] := (Fin m → ℕ) →₀ K

/-- Exponents: a low part and a high pattern per block. -/
abbrev Ex (m : ℕ) := (Fin m → ℕ) × (Fin m → ℕ)

/-- The depth-uniform ambient algebra (as a `K`-module). -/
abbrev Al (m : ℕ) (K : Type*) [Field K] := Ex m →₀ K

def lowWt (l : Fin m → ℕ) : ℕ := ∑ i, wt (l i)

def pwt (a : Ex m) : ℕ := lowWt a.1 + lowWt a.2

@[simp] theorem lowWt_zero : lowWt (0 : Fin m → ℕ) = 0 := by simp [lowWt]

theorem lowWt_update (l : Fin m → ℕ) (i : Fin m) (y : ℕ) :
    lowWt (Function.update l i y) + wt (l i) = lowWt l + wt y := by
  classical
  unfold lowWt
  have h1 : ∑ j, wt (Function.update l i y j) =
      ∑ j, Function.update (fun j => wt (l j)) i (wt y) j := by
    apply Finset.sum_congr rfl
    intro j _
    exact Function.apply_update (fun _ n => wt n) l i y j
  rw [h1, Finset.sum_update_of_mem (Finset.mem_univ i),
    Finset.sum_eq_add_sum_diff_singleton (Finset.mem_univ i) (fun j => wt (l j))]
  ring

def LowBox (k0 : ℕ) (l : Fin m → ℕ) : Prop := ∀ i, l i < 2 ^ k0

def InBox (k0 k : ℕ) (a : Ex m) : Prop := LowBox k0 a.1 ∧ LowBox (k - k0) a.2

def Boxed (k0 k : ℕ) (w : Al m K) : Prop := ∀ a ∈ w.support, InBox k0 k a

def WtLT (D : ℕ) (w : Al m K) : Prop := ∀ a ∈ w.support, pwt a < D

theorem lowBox_update {k0 : ℕ} {l : Fin m → ℕ} (hl : LowBox k0 l) (i : Fin m) {y : ℕ}
    (hy : y < 2 ^ k0) : LowBox k0 (Function.update l i y) := by
  intro i'
  rcases eq_or_ne i' i with rfl | h
  · simpa using hy
  · rw [Function.update_of_ne h]
    exact hl i'

theorem lowBox_zero (k0 : ℕ) : LowBox k0 (0 : Fin m → ℕ) := fun _ => Nat.two_pow_pos _

/-- The unit vector `e_i`. -/
def unit (i : Fin m) : Fin m → ℕ := Function.update 0 i 1

theorem unit_apply_same (i : Fin m) : unit i i = 1 := by simp [unit]

theorem unit_ne_zero (i : Fin m) : unit i ≠ 0 := by
  intro h
  have := congrFun h i
  simp [unit] at this

theorem lowBox_unit {n : ℕ} (hn : 0 < n) (i : Fin m) : LowBox n (unit i) := by
  intro i'
  have : 1 < 2 ^ n := Nat.one_lt_two_pow (by omega)
  rcases eq_or_ne i' i with rfl | h
  · simpa [unit] using this
  · simp only [unit, Function.update_of_ne h, Pi.zero_apply]
    exact Nat.two_pow_pos _

/-! ## Translations with carries -/

variable (k0 k : ℕ) (i : Fin m) (j : ℕ)

/-- A low variable (`j < k0`): add `2^j` to the low part, carrying into the high pattern. -/
def stepLow (a : Ex m) : Ex m :=
  if a.1 i + 2 ^ j < 2 ^ k0 then (Function.update a.1 i (a.1 i + 2 ^ j), a.2)
  else (Function.update a.1 i (a.1 i + 2 ^ j - 2 ^ k0),
    Function.update a.2 i ((a.2 i + 1) % 2 ^ (k - k0)))

/-- A high variable (`j ≥ k0`): add `2^(j-k0)` to the high pattern, cyclically. -/
def stepHigh (a : Ex m) : Ex m :=
  (a.1, Function.update a.2 i ((a.2 i + 2 ^ (j - k0)) % 2 ^ (k - k0)))

def step (a : Ex m) : Ex m := if j < k0 then stepLow k0 k i j a else stepHigh k0 k i j a

def unstepLow (a : Ex m) : Ex m :=
  if 2 ^ j ≤ a.1 i then (Function.update a.1 i (a.1 i - 2 ^ j), a.2)
  else (Function.update a.1 i (a.1 i + 2 ^ k0 - 2 ^ j),
    Function.update a.2 i ((a.2 i + 2 ^ (k - k0) - 1) % 2 ^ (k - k0)))

def unstepHigh (a : Ex m) : Ex m :=
  (a.1, Function.update a.2 i ((a.2 i + 2 ^ (k - k0) - 2 ^ (j - k0)) % 2 ^ (k - k0)))

def unstep (a : Ex m) : Ex m :=
  if j < k0 then unstepLow k0 k i j a else unstepHigh k0 k i j a

variable {k0 k i j}

theorem step_inBox {a : Ex m} (ha : InBox k0 k a) :
    InBox k0 k (step k0 k i j a) := by
  have hN : 0 < 2 ^ (k - k0) := Nat.two_pow_pos _
  have hpj : 0 < 2 ^ j := Nat.two_pow_pos _
  have hpk : 0 < 2 ^ k0 := Nat.two_pow_pos _
  unfold step
  split_ifs with hjl
  · unfold stepLow
    have hjk : 2 ^ j < 2 ^ k0 := Nat.pow_lt_pow_right (by norm_num) hjl
    split_ifs with hc
    · exact ⟨lowBox_update ha.1 i hc, ha.2⟩
    · refine ⟨lowBox_update ha.1 i (by have := ha.1 i; omega), lowBox_update ha.2 i ?_⟩
      exact Nat.mod_lt _ hN
  · exact ⟨ha.1, lowBox_update ha.2 i (Nat.mod_lt _ hN)⟩

theorem unstep_inBox {a : Ex m} (ha : InBox k0 k a) :
    InBox k0 k (unstep k0 k i j a) := by
  have hN : 0 < 2 ^ (k - k0) := Nat.two_pow_pos _
  have hpj : 0 < 2 ^ j := Nat.two_pow_pos _
  have hpk : 0 < 2 ^ k0 := Nat.two_pow_pos _
  unfold unstep
  split_ifs with hjl
  · unfold unstepLow
    have hjk : 2 ^ j < 2 ^ k0 := Nat.pow_lt_pow_right (by norm_num) hjl
    split_ifs with hc
    · exact ⟨lowBox_update ha.1 i (by have := ha.1 i; omega), ha.2⟩
    · refine ⟨lowBox_update ha.1 i (by omega), lowBox_update ha.2 i (Nat.mod_lt _ hN)⟩
  · exact ⟨ha.1, lowBox_update ha.2 i (Nat.mod_lt _ hN)⟩

theorem step_unstep (hk : k0 < k) (hj : j < k) {a : Ex m} (ha : InBox k0 k a) :
    step k0 k i j (unstep k0 k i j a) = a := by
  have hN : 1 ≤ 2 ^ (k - k0) := Nat.one_le_two_pow
  have hpj : 0 < 2 ^ j := Nat.two_pow_pos _
  have hpk : 0 < 2 ^ k0 := Nat.two_pow_pos _
  obtain ⟨l, q⟩ := a
  have hl : l i < 2 ^ k0 := ha.1 i
  have hq : q i < 2 ^ (k - k0) := ha.2 i
  unfold step unstep
  split_ifs with hjl
  · unfold unstepLow stepLow
    have hjk : 2 ^ j < 2 ^ k0 := Nat.pow_lt_pow_right (by norm_num) hjl
    by_cases hc : 2 ^ j ≤ l i
    · rw [if_pos hc]
      dsimp only
      rw [Function.update_self, if_pos (by omega), Function.update_idem, Nat.sub_add_cancel hc,
        Function.update_eq_self]
    · rw [if_neg hc]
      dsimp only
      rw [Function.update_self, Function.update_self, if_neg (by omega), Function.update_idem,
        Function.update_idem]
      have h1 : l i + 2 ^ k0 - 2 ^ j + 2 ^ j - 2 ^ k0 = l i := by omega
      rw [h1, Function.update_eq_self, mod_add_sub_cancel _ _ _ hq hN, Function.update_eq_self]
  · unfold unstepHigh stepHigh
    dsimp only
    have hd : 2 ^ (j - k0) ≤ 2 ^ (k - k0) := Nat.pow_le_pow_right (by norm_num) (by omega)
    rw [Function.update_self, Function.update_idem, mod_add_sub_cancel _ _ _ (ha.2 i) hd,
      Function.update_eq_self]

theorem unstep_step (hk : k0 < k) (hj : j < k) {a : Ex m} (ha : InBox k0 k a) :
    unstep k0 k i j (step k0 k i j a) = a := by
  have hN : 1 ≤ 2 ^ (k - k0) := Nat.one_le_two_pow
  have hpj : 0 < 2 ^ j := Nat.two_pow_pos _
  have hpk : 0 < 2 ^ k0 := Nat.two_pow_pos _
  obtain ⟨l, q⟩ := a
  have hl : l i < 2 ^ k0 := ha.1 i
  have hq : q i < 2 ^ (k - k0) := ha.2 i
  unfold step unstep
  split_ifs with hjl
  · unfold unstepLow stepLow
    have hjk : 2 ^ j < 2 ^ k0 := Nat.pow_lt_pow_right (by norm_num) hjl
    by_cases hc : l i + 2 ^ j < 2 ^ k0
    · rw [if_pos hc]
      dsimp only
      rw [Function.update_self, if_pos (by omega), Function.update_idem, Nat.add_sub_cancel,
        Function.update_eq_self]
    · rw [if_neg hc]
      dsimp only
      rw [Function.update_self, Function.update_self, if_neg (by omega), Function.update_idem,
        Function.update_idem]
      have h1 : l i + 2 ^ j - 2 ^ k0 + 2 ^ k0 - 2 ^ j = l i := by omega
      rw [h1, Function.update_eq_self, mod_sub_add_cancel _ _ _ hq hN, Function.update_eq_self]
  · unfold unstepHigh stepHigh
    dsimp only
    have hd : 2 ^ (j - k0) ≤ 2 ^ (k - k0) := Nat.pow_le_pow_right (by norm_num) (by omega)
    rw [Function.update_self, Function.update_idem, mod_sub_add_cancel _ _ _ (ha.2 i) hd,
      Function.update_eq_self]

theorem step_injOn (hk : k0 < k) (hj : j < k) :
    Set.InjOn (step k0 k i j) {a : Ex m | InBox k0 k a} := by
  intro a ha b hb hab
  rw [← unstep_step hk hj ha, ← unstep_step hk hj hb, hab]

/-- One translation raises the weight by at most one, whatever the carries and wraps. -/
theorem pwt_step_le (hk : k0 < k) (hj : j < k) {a : Ex m} (ha : InBox k0 k a) :
    pwt (step k0 k i j a) ≤ pwt a + 1 := by
  have hpj : 0 < 2 ^ j := Nat.two_pow_pos _
  have hpk : 0 < 2 ^ k0 := Nat.two_pow_pos _
  obtain ⟨l, q⟩ := a
  have hl : l i < 2 ^ k0 := ha.1 i
  have hq : q i < 2 ^ (k - k0) := ha.2 i
  unfold step pwt
  split_ifs with hjl
  · unfold stepLow
    split_ifs with hc
    · dsimp only
      have h1 := lowWt_update l i (l i + 2 ^ j)
      have h2 := wt_add_two_pow_le j (l i)
      omega
    · dsimp only
      have h1 := lowWt_update l i (l i + 2 ^ j - 2 ^ k0)
      have h2 := carry_wt k0 j (l i) hjl hl (not_lt.mp hc)
      have h3 := lowWt_update q i ((q i + 1) % 2 ^ (k - k0))
      have h4 := wt_mod_two_pow_add_le (k - k0) 0 (q i) (by omega) hq
      rw [pow_zero] at h4
      omega
  · unfold stepHigh
    dsimp only
    have h1 := lowWt_update q i ((q i + 2 ^ (j - k0)) % 2 ^ (k - k0))
    have h2 := wt_mod_two_pow_add_le (k - k0) (j - k0) (q i) (by omega) hq
    omega

/-! ## The translation operators, components and embeddings -/

variable (k0 k i j)

/-- Multiplication by the chain variable `u_{i,j}` at depth `k`, as a linear map. -/
def translate : Al m K →ₗ[K] Al m K := Finsupp.lmapDomain K K (step k0 k i j)

variable {k0 k i j}

theorem translate_support (w : Al m K) :
    (translate (K := K) k0 k i j w).support ⊆ w.support.image (step k0 k i j) := by
  classical
  exact Finsupp.mapDomain_support

theorem translate_boxed {w : Al m K} (hw : Boxed k0 k w) :
    Boxed k0 k (translate (K := K) k0 k i j w) := by
  classical
  intro b hb
  obtain ⟨a, ha, rfl⟩ := Finset.mem_image.mp (translate_support w hb)
  exact step_inBox (i := i) (j := j) (hw a ha)

theorem translate_apply (hk : k0 < k) (hj : j < k) {w : Al m K} (hw : Boxed k0 k w) {b : Ex m}
    (hb : InBox k0 k b) : translate (K := K) k0 k i j w b = w (unstep k0 k i j b) := by
  have hS : (w.support : Set (Ex m)) ⊆ {a | InBox k0 k a} := fun a ha => hw a ha
  have key := Finsupp.mapDomain_apply' {a : Ex m | InBox k0 k a} (f := step k0 k i j) w hS
    (step_injOn (i := i) hk hj) (unstep_inBox (i := i) (j := j) hb)
  rw [step_unstep hk hj hb] at key
  exact key

theorem translate_apply_of_not_inBox {w : Al m K} (hw : Boxed k0 k w)
    {b : Ex m} (hb : ¬ InBox k0 k b) : translate (K := K) k0 k i j w b = 0 :=
  Finsupp.not_mem_support_iff.mp (fun h => hb (translate_boxed hw b h))

/-- The `q`-component of an element, read as a low element. -/
def piece (q : Fin m → ℕ) : Al m K →ₗ[K] Lo m K :=
  Finsupp.lcomapDomain (fun l : Fin m → ℕ => (l, q)) (fun _ _ h => (Prod.ext_iff.mp h).1)

@[simp] theorem piece_apply (q : Fin m → ℕ) (w : Al m K) (l : Fin m → ℕ) :
    piece q w l = w (l, q) := rfl

/-- A low element placed at high pattern `q`. -/
def embed (q : Fin m → ℕ) : Lo m K →ₗ[K] Al m K :=
  Finsupp.lmapDomain K K (fun l : Fin m → ℕ => (l, q))

theorem embed_apply_same (q : Fin m → ℕ) (v : Lo m K) (l : Fin m → ℕ) :
    embed q v (l, q) = v l :=
  Finsupp.mapDomain_apply (fun _ _ h => (Prod.ext_iff.mp h).1) v l

theorem embed_apply_ne (q : Fin m → ℕ) (v : Lo m K) (l : Fin m → ℕ) {q' : Fin m → ℕ}
    (hq : q' ≠ q) : embed q v (l, q') = 0 := by
  apply Finsupp.mapDomain_notin_range
  rintro ⟨l', hl'⟩
  exact hq (Prod.ext_iff.mp hl').2.symm

theorem piece_embed (q : Fin m → ℕ) (v : Lo m K) : piece q (embed q v) = v := by
  ext l
  rw [piece_apply, embed_apply_same]

theorem piece_embed_ne (q : Fin m → ℕ) (v : Lo m K) {q' : Fin m → ℕ} (hq : q' ≠ q) :
    piece q' (embed q v) = 0 := by
  ext l
  rw [piece_apply, embed_apply_ne q v l hq, Finsupp.zero_apply]

theorem embed_boxed {k0 k : ℕ} {q : Fin m → ℕ} (hq : LowBox (k - k0) q) {v : Lo m K}
    (hv : ∀ l ∈ v.support, LowBox k0 l) : Boxed k0 k (embed q v) := by
  classical
  intro a ha
  obtain ⟨l, hl, rfl⟩ := Finset.mem_image.mp (Finsupp.mapDomain_support ha)
  exact ⟨hv l hl, hq⟩

theorem piece_support_lowBox {k0 k : ℕ} {w : Al m K} (hw : Boxed k0 k w) (q : Fin m → ℕ) :
    ∀ l ∈ (piece q w).support, LowBox k0 l := by
  intro l hl
  rw [Finsupp.mem_support_iff, piece_apply] at hl
  exact (hw (l, q) (Finsupp.mem_support_iff.mpr hl)).1

theorem piece_single (l0 q0 : Fin m → ℕ) (c : K) (q : Fin m → ℕ) :
    piece q (Finsupp.single (l0, q0) c) = if q = q0 then Finsupp.single l0 c else 0 := by
  classical
  by_cases hq : q = q0
  · subst hq
    rw [if_pos rfl]
    ext l
    rw [piece_apply, Finsupp.single_apply, Finsupp.single_apply]
    simp [Prod.ext_iff]
  · rw [if_neg hq]
    ext l
    rw [piece_apply, Finsupp.single_apply, Finsupp.zero_apply, if_neg]
    intro h
    exact hq (Prod.ext_iff.mp h).2.symm

/-! ## The free high extension `W` -/

variable (k0 D : ℕ) (S : Submodule K (Lo m K)) (k : ℕ)

/-- `W k0 D S k`: elements of the depth-`k` box whose every high component lies in `S`, with
the weight budget `wt l + wt q ≤ D` on each monomial. -/
def W : Submodule K (Al m K) where
  carrier := {w | Boxed k0 k w ∧ ∀ q, piece q w ∈ S ∧
    ∀ l ∈ (piece q w).support, lowWt l + lowWt q ≤ D}
  zero_mem' := by
    refine ⟨fun a ha => by simp at ha, fun q => ⟨by simp, fun l hl => by simp at hl⟩⟩
  add_mem' := by
    classical
    rintro w₁ w₂ ⟨hb₁, hp₁⟩ ⟨hb₂, hp₂⟩
    refine ⟨fun a ha => ?_, fun q => ⟨?_, fun l hl => ?_⟩⟩
    · rcases Finset.mem_union.mp (Finsupp.support_add ha) with h | h
      · exact hb₁ a h
      · exact hb₂ a h
    · rw [map_add]
      exact S.add_mem (hp₁ q).1 (hp₂ q).1
    · rw [map_add] at hl
      rcases Finset.mem_union.mp (Finsupp.support_add hl) with h | h
      · exact (hp₁ q).2 l h
      · exact (hp₂ q).2 l h
  smul_mem' := by
    rintro c w ⟨hb, hp⟩
    refine ⟨fun a ha => hb a (Finsupp.support_smul ha), fun q => ⟨?_, fun l hl => ?_⟩⟩
    · rw [map_smul]
      exact S.smul_mem c (hp q).1
    · rw [map_smul] at hl
      exact (hp q).2 l (Finsupp.support_smul hl)

variable {k0 D S k}

theorem mem_W {w : Al m K} : w ∈ W k0 D S k ↔ Boxed k0 k w ∧ ∀ q, piece q w ∈ S ∧
    ∀ l ∈ (piece q w).support, lowWt l + lowWt q ≤ D := Iff.rfl

/-- Part (ii) of N5: the constant `1 = δ_0` lies in `W` iff it lies in the base. -/
theorem single_zero_mem_W_iff :
    Finsupp.single ((0 : Fin m → ℕ), (0 : Fin m → ℕ)) (1 : K) ∈ W k0 D S k ↔
      Finsupp.single (0 : Fin m → ℕ) (1 : K) ∈ S := by
  classical
  constructor
  · intro h
    have := (h.2 0).1
    rwa [piece_single, if_pos rfl] at this
  · intro hS
    refine ⟨fun a ha => ?_, fun q => ?_⟩
    · have : a = (0, 0) := Finset.mem_singleton.mp (Finsupp.support_single_subset ha)
      subst this
      exact ⟨lowBox_zero k0, lowBox_zero (k - k0)⟩
    · rw [piece_single]
      split_ifs with hq
      · subst hq
        refine ⟨hS, fun l hl => ?_⟩
        have : l = 0 := Finset.mem_singleton.mp (Finsupp.support_single_subset hl)
        subst this
        simp
      · exact ⟨S.zero_mem, fun l hl => by simp at hl⟩

/-- A low element of `S` within the budget, placed at pattern `0`, lies in `W`. -/
theorem embed_zero_mem_W {v : Lo m K} (hv : v ∈ S) (hvbox : ∀ l ∈ v.support, LowBox k0 l)
    (hvwt : ∀ l ∈ v.support, lowWt l ≤ D) : embed 0 v ∈ W k0 D S k := by
  refine ⟨embed_boxed (lowBox_zero _) hvbox, fun q => ?_⟩
  by_cases hq : q = 0
  · subst hq
    rw [piece_embed]
    exact ⟨hv, fun l hl => by simpa using hvwt l hl⟩
  · rw [piece_embed_ne 0 v hq]
    exact ⟨S.zero_mem, fun l hl => by simp at hl⟩

/-- Closure under every admissible variable translation at depth `k`. -/
def Closed (k0 k D : ℕ) (W : Submodule K (Al m K)) : Prop :=
  ∀ w ∈ W, WtLT D w → ∀ (i : Fin m) (j : ℕ), j < k → translate k0 k i j w ∈ W

/-! ## Components of a translate -/

/-- High variable: the components are permuted. -/
theorem piece_translate_high (hk : k0 < k) (hj : j < k) (hjh : k0 ≤ j) {w : Al m K}
    (hw : Boxed k0 k w) {q : Fin m → ℕ} (hq : LowBox (k - k0) q) :
    piece q (translate (K := K) k0 k i j w) =
      piece (Function.update q i ((q i + 2 ^ (k - k0) - 2 ^ (j - k0)) % 2 ^ (k - k0))) w := by
  ext l
  rw [piece_apply, piece_apply]
  by_cases hl : LowBox k0 l
  · rw [translate_apply hk hj hw ⟨hl, hq⟩]
    simp only [unstep, if_neg (not_lt.mpr hjh), unstepHigh]
  · rw [translate_apply_of_not_inBox hw (fun h => hl h.1)]
    symm
    exact Finsupp.not_mem_support_iff.mp (fun h => hl (hw _ h).1)

/-- The predecessor pattern of `q` in block `i` (cyclic). -/
def qPrev (k0 k : ℕ) (i : Fin m) (q : Fin m → ℕ) : Fin m → ℕ :=
  Function.update q i ((q i + 2 ^ (k - k0) - 1) % 2 ^ (k - k0))

theorem lowBox_qPrev {q : Fin m → ℕ} (hq : LowBox (k - k0) q) : LowBox (k - k0) (qPrev k0 k i q) :=
  lowBox_update hq i (Nat.mod_lt _ (Nat.two_pow_pos _))

theorem two_pow_succ_sub (k0 : ℕ) : 2 ^ (k0 + 1 - k0) = 2 := by
  rw [Nat.add_sub_cancel_left, pow_one]

/-- Low variable: the `q`-component of the depth-`k` translate is the no-carry part of the
depth-`(k0+1)` translate of the `q`-component plus the carried part of the depth-`(k0+1)`
translate of the predecessor component. -/
theorem piece_translate_low (hk : k0 < k) (hj : j < k) (hjl : j < k0) {w : Al m K}
    (hw : Boxed k0 k w) {q : Fin m → ℕ} (hq : LowBox (k - k0) q) :
    piece q (translate (K := K) k0 k i j w) =
      piece 0 (translate k0 (k0 + 1) i j (embed 0 (piece q w))) +
        piece (unit i) (translate k0 (k0 + 1) i j (embed 0 (piece (qPrev k0 k i q) w))) := by
  classical
  have h2 := two_pow_succ_sub k0
  have hk1 : k0 < k0 + 1 := Nat.lt_succ_self k0
  have hj1 : j < k0 + 1 := by omega
  have hb0 : LowBox (k0 + 1 - k0) (0 : Fin m → ℕ) := lowBox_zero _
  have hbu : LowBox (k0 + 1 - k0) (unit i) := lowBox_unit (by omega) i
  ext l
  rw [Finsupp.add_apply, piece_apply, piece_apply, piece_apply]
  by_cases hl : LowBox k0 l
  · rw [translate_apply hk hj hw ⟨hl, hq⟩,
      translate_apply hk1 hj1 (embed_boxed hb0 (piece_support_lowBox hw q)) ⟨hl, hb0⟩,
      translate_apply hk1 hj1 (embed_boxed hb0 (piece_support_lowBox hw _)) ⟨hl, hbu⟩]
    simp only [unstep, if_pos hjl, unstepLow]
    try dsimp only
    by_cases hc : 2 ^ j ≤ l i
    · rw [if_pos hc, if_pos hc, if_pos hc, embed_apply_same, embed_apply_ne _ _ _ (unit_ne_zero i),
        piece_apply, add_zero]
    · rw [if_neg hc, if_neg hc, if_neg hc]
      have e1 : Function.update (0 : Fin m → ℕ) i
          (((0 : Fin m → ℕ) i + 2 ^ (k0 + 1 - k0) - 1) % 2 ^ (k0 + 1 - k0)) = unit i := by
        rw [h2]
        simp [unit]
      have e2 : Function.update (unit i) i
          ((unit i i + 2 ^ (k0 + 1 - k0) - 1) % 2 ^ (k0 + 1 - k0)) = 0 := by
        rw [h2, unit_apply_same]
        simp [unit]
      rw [e1, e2, embed_apply_ne _ _ _ (unit_ne_zero i), embed_apply_same, piece_apply, zero_add]
      rfl
  · rw [translate_apply_of_not_inBox hw (fun h => hl h.1),
      translate_apply_of_not_inBox (embed_boxed hb0 (piece_support_lowBox hw q))
        (fun h => hl h.1),
      translate_apply_of_not_inBox (embed_boxed hb0 (piece_support_lowBox hw _))
        (fun h => hl h.1), add_zero]

/-! ## Part (iii): closedness propagates from depth `k0+1` to every depth -/

/-- Components of an element of `W` of weight below `D`, placed at pattern `0`, are admissible
inputs for the depth-`(k0+1)` closedness hypothesis. -/
theorem embed_piece_mem_W {w : Al m K} (hw : w ∈ W k0 D S k) (q : Fin m → ℕ) :
    embed 0 (piece q w) ∈ W k0 D S (k0 + 1) :=
  embed_zero_mem_W (hw.2 q).1 (piece_support_lowBox hw.1 q)
    (fun l hl => by have := (hw.2 q).2 l hl; omega)

theorem embed_piece_wtLT {w : Al m K} (hlow : WtLT D w) (q : Fin m → ℕ) :
    WtLT D (embed (K := K) 0 (piece q w)) := by
  classical
  intro a ha
  obtain ⟨l, hl, rfl⟩ := Finset.mem_image.mp (Finsupp.mapDomain_support ha)
  rw [Finsupp.mem_support_iff, piece_apply] at hl
  have := hlow (l, q) (Finsupp.mem_support_iff.mpr hl)
  simp only [pwt, lowWt_zero] at this ⊢
  omega

/-- **Depth reduction (N5 (iii)).** If the free extension of `S` is closed under all admissible
variable translations at depth `k0+1`, it is closed at every depth `k ≥ k0+1`. -/
theorem W_closed (base : Closed k0 (k0 + 1) D (W k0 D S (k0 + 1))) (k : ℕ) (hk : k0 + 1 ≤ k) :
    Closed k0 k D (W k0 D S k) := by
  classical
  intro w hw hlow i j hj
  have hk' : k0 < k := by omega
  have hbox := hw.1
  refine ⟨translate_boxed hbox, fun q => ⟨?_, fun l hl => ?_⟩⟩
  · by_cases hq : LowBox (k - k0) q
    · by_cases hjl : j < k0
      · rw [piece_translate_low hk' hj hjl hbox hq]
        apply S.add_mem
        · exact ((base _ (embed_piece_mem_W hw q) (embed_piece_wtLT hlow q) i j (by omega)).2
            0).1
        · exact ((base _ (embed_piece_mem_W hw (qPrev k0 k i q))
            (embed_piece_wtLT hlow (qPrev k0 k i q)) i j (by omega)).2 (unit i)).1
      · rw [piece_translate_high hk' hj (not_lt.mp hjl) hbox hq]
        exact (hw.2 _).1
    · have : piece q (translate (K := K) k0 k i j w) = 0 := by
        ext l
        rw [piece_apply, translate_apply_of_not_inBox hbox (fun h => hq h.2),
          Finsupp.zero_apply]
      rw [this]
      exact S.zero_mem
  · rw [Finsupp.mem_support_iff, piece_apply] at hl
    obtain ⟨a, ha, hstep⟩ :=
      Finset.mem_image.mp (translate_support w (Finsupp.mem_support_iff.mpr hl))
    have h1 := pwt_step_le (i := i) hk' hj (hbox a ha)
    have h2 := hlow a ha
    rw [hstep] at h1
    simp only [pwt] at h1 h2
    omega

end ICFallDegrees.DepthReduction
