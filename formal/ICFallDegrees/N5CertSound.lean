import ICFallDegrees.N5Cert
import ICFallDegrees.ChainModelMap

/-!
# Soundness of the N5 base certificate checker

`N5Cert` is a computable checker on the certificate data `d : Data` (a reduced echelon basis of a
depth-7 subspace, packed into naturals).  This file gives the data its meaning and proves that
the checks imply the three per-instance hypotheses of `ChainModelMap.not_one_of_base` /
`n5_lfd_ge_five` with `k0 = 7`, `D = 4`:

* `S d` is the span of the row vectors `rowVec d i`, `i < d.nRows`, each `δ_{pivot} + (free part
  read from the packed row)`;
* `member d u = true → toLo u ∈ S d` (`toLo_mem_of_member`): the echelon membership test;
* `one_not_mem`: `δ₀ ∉ S d` (B1);
* `closed`: `Closed 7 8 4 (W 7 4 (S d) 8)` (B2), through the finite shift checks on the rows of
  pivot weight `≤ 3` and the component formulas `piece_translate_low/high`;
* `Φ_embed_eCode` and the generator lemmas, for (G).

Everything is relative to `Valid d`, the conjunction of the Bool checks, which the per-instance
file discharges by `decide +kernel`.  Only `propext`, `Classical.choice` and `Quot.sound` are used.
-/

noncomputable section

namespace ICFallDegrees.N5CertSound

open ICFallDegrees.N5Cert DepthReduction BinaryWeight

/-! ## Nat lemmas on packed slots -/

section Nat

theorem get16_lt (enc i : ℕ) : get16 enc i < 65536 := by
  unfold get16
  have : (65535 : ℕ) = 2 ^ 16 - 1 := by norm_num
  rw [this, Nat.and_two_pow_sub_one_eq_mod]
  exact Nat.mod_lt _ (by norm_num)

/-- `∑ s_f B^f`. -/
def packSlots (B : ℕ) : List ℕ → ℕ
  | [] => 0
  | s :: rest => s + B * packSlots B rest

theorem unpack_length (B n R : ℕ) : (unpack B n R).length = n := by
  induction n generalizing R with
  | zero => rfl
  | succ n ih => simp [unpack, ih]

theorem unpack_packSlots {B : ℕ} (hB : 0 < B) :
    ∀ sl : List ℕ, (∀ s ∈ sl, s < B) → unpack B sl.length (packSlots B sl) = sl := by
  intro sl
  induction sl with
  | nil => intro _; rfl
  | cons s rest ih =>
    intro h
    have hs : s < B := h s (List.mem_cons_self)
    simp only [List.length_cons, unpack, packSlots]
    rw [Nat.add_mul_mod_self_left, Nat.mod_eq_of_lt hs, Nat.add_mul_div_left _ _ hB,
      Nat.div_eq_of_lt hs, zero_add, ih (fun x hx => h x (List.mem_cons_of_mem _ hx))]

theorem packSlots_unpack {B : ℕ} (hB : 0 < B) :
    ∀ (n R : ℕ), R < B ^ n → packSlots B (unpack B n R) = R := by
  intro n
  induction n with
  | zero =>
    intro R hR
    simp at hR
    simp [unpack, packSlots, hR]
  | succ n ih =>
    intro R hR
    simp only [unpack, packSlots]
    rw [ih (R / B) (by rw [Nat.div_lt_iff_lt_mul hB, ← pow_succ]; exact hR), Nat.mod_add_div]

theorem packSlots_zipWith_add (B : ℕ) :
    ∀ (a b : List ℕ), a.length = b.length →
      packSlots B (List.zipWith (· + ·) a b) = packSlots B a + packSlots B b := by
  intro a
  induction a with
  | nil => intro b hb; cases b with | nil => rfl | cons => simp at hb
  | cons x a ih =>
    intro b hb
    cases b with
    | nil => simp at hb
    | cons y b =>
      simp only [List.zipWith_cons_cons, packSlots, List.length_cons] at hb ⊢
      rw [ih b (by omega)]
      ring

theorem packSlots_map_mul (B x : ℕ) :
    ∀ s : List ℕ, packSlots B (s.map (x * ·)) = x * packSlots B s := by
  intro s
  induction s with
  | nil => simp [packSlots]
  | cons y s ih => simp only [List.map_cons, packSlots, ih]; ring

/-- The unit slot vector of length `n` with value `c` at position `f`. -/
def unitSl : ℕ → ℕ → ℕ → List ℕ
  | 0, _, _ => []
  | n+1, 0, c => c :: List.replicate n 0
  | n+1, f+1, c => 0 :: unitSl n f c

theorem unitSl_length : ∀ n f c, (unitSl n f c).length = n
  | 0, _, _ => rfl
  | n+1, 0, _ => by simp [unitSl]
  | n+1, f+1, c => by simp [unitSl, unitSl_length n f c]

theorem packSlots_replicate_zero (B : ℕ) : ∀ n, packSlots B (List.replicate n 0) = 0
  | 0 => rfl
  | n+1 => by simp [List.replicate_succ, packSlots, packSlots_replicate_zero B n]

theorem packSlots_unitSl (B : ℕ) :
    ∀ n f c, f < n → packSlots B (unitSl n f c) = c * B ^ f
  | 0, _, _, h => by omega
  | n+1, 0, c, _ => by simp [unitSl, packSlots, packSlots_replicate_zero]
  | n+1, f+1, c, h => by
    simp only [unitSl, packSlots, packSlots_unitSl B n f c (by omega)]
    ring

theorem mem_unitSl_le : ∀ n f c, ∀ s ∈ unitSl n f c, s ≤ c
  | 0, _, _, s, h => by simp [unitSl] at h
  | n+1, 0, c, s, h => by
    simp only [unitSl, List.mem_cons, List.mem_replicate] at h
    rcases h with rfl | ⟨_, rfl⟩ <;> omega
  | n+1, f+1, c, s, h => by
    simp only [unitSl, List.mem_cons] at h
    rcases h with rfl | h
    · exact Nat.zero_le _
    · exact mem_unitSl_le n f c s h

theorem mem_replicate_zero_le {n s : ℕ} (h : s ∈ List.replicate n 0) : s ≤ 0 := by
  rw [List.mem_replicate] at h; omega

theorem mem_zipWith_add_le {A C : ℕ} :
    ∀ (a b : List ℕ), (∀ x ∈ a, x ≤ A) → (∀ y ∈ b, y ≤ C) →
      ∀ s ∈ List.zipWith (· + ·) a b, s ≤ A + C
  | [], _, _, _, s, h => by simp at h
  | _ :: _, [], _, _, s, h => by simp at h
  | x :: a, y :: b, ha, hb, s, h => by
    simp only [List.zipWith_cons_cons, List.mem_cons] at h
    rcases h with rfl | h
    · exact Nat.add_le_add (ha x (List.mem_cons_self)) (hb y (List.mem_cons_self))
    · exact mem_zipWith_add_le a b (fun x hx => ha x (List.mem_cons_of_mem _ hx))
        (fun y hy => hb y (List.mem_cons_of_mem _ hy)) s h

theorem slotsZero_sound (B p : ℕ) :
    ∀ n R, slotsZero B p n R = true → ∀ s ∈ unpack B n R, p ∣ s := by
  intro n
  induction n with
  | zero => intro R _ s hs; simp [unpack] at hs
  | succ n ih =>
    intro R h s hs
    simp only [slotsZero, Bool.and_eq_true, beq_iff_eq] at h
    simp only [unpack, List.mem_cons] at hs
    rcases hs with rfl | hs
    · exact Nat.dvd_of_mod_eq_zero h.1
    · exact ih _ h.2 s hs

theorem slotsLt_sound {B : ℕ} (hB : 0 < B) (bound : ℕ) :
    ∀ n R, slotsLt B bound n R = true → (∀ s ∈ unpack B n R, s < bound) ∧ R < B ^ n := by
  intro n
  induction n with
  | zero =>
    intro R h
    simp only [slotsLt, beq_iff_eq] at h
    subst h
    exact ⟨fun s hs => by simp [unpack] at hs, by simp⟩
  | succ n ih =>
    intro R h
    simp only [slotsLt, Bool.and_eq_true, decide_eq_true_eq] at h
    obtain ⟨h1, h2⟩ := ih _ h.2
    refine ⟨fun s hs => ?_, ?_⟩
    · simp only [unpack, List.mem_cons] at hs
      rcases hs with rfl | hs
      · exact h.1
      · exact h1 s hs
    · rw [pow_succ, ← Nat.div_lt_iff_lt_mul hB]
      exact h2

theorem allBelow_sound (P : ℕ → Bool) : ∀ n, allBelow P n = true → ∀ i < n, P i = true := by
  intro n
  induction n with
  | zero => intro _ i hi; omega
  | succ n ih =>
    intro h i hi
    simp only [allBelow, Bool.and_eq_true] at h
    rcases Nat.lt_succ_iff_lt_or_eq.mp hi with hi | rfl
    · exact ih h.2 i hi
    · exact h.1

theorem rangeAll_sound (P : ℕ → Bool) (lo : ℕ) :
    ∀ n, rangeAll P lo n = true → ∀ i, lo ≤ i → i < lo + n → P i = true := by
  intro n
  induction n with
  | zero => intro _ i _ hi; omega
  | succ n ih =>
    intro h i hlo hi
    simp only [rangeAll, Bool.and_eq_true] at h
    rcases Nat.lt_or_ge i (lo + n) with h' | h'
    · exact ih h.2 i hlo h'
    · have : i = lo + n := by omega
      subst this
      exact h.1

theorem popcount_eq_wt : ∀ fuel n, n < 2 ^ fuel → popcount fuel n = wt n := by
  intro fuel
  induction fuel with
  | zero =>
    intro n hn
    simp at hn
    subst hn
    simp [popcount]
  | succ f ih =>
    intro n hn
    rw [popcount, wt_eq n, ih (n / 2) (by rw [pow_succ] at hn; omega)]

end Nat

/-! ## Meaning of the data -/

variable {K : Type*} [Field K]

/-- The low exponent vector of a code `m = l₀ + 128 l₁`. -/
def col (m : ℕ) : Fin 2 → ℕ := ![m % 128, m / 128]

theorem col_zero (m : ℕ) : col m 0 = m % 128 := rfl
theorem col_one (m : ℕ) : col m 1 = m / 128 := rfl

theorem col_injOn {m m' : ℕ} (h : col m = col m') : m = m' := by
  have h0 := congrFun h 0
  have h1 := congrFun h 1
  simp only [col_zero, col_one] at h0 h1
  omega

theorem col_eq_zero_iff (m : ℕ) : col m = 0 ↔ m = 0 := by
  constructor
  · intro h
    have h0 := congrFun h 0
    have h1 := congrFun h 1
    simp only [col_zero, col_one, Pi.zero_apply] at h0 h1
    omega
  · rintro rfl
    funext i
    fin_cases i <;> rfl

theorem lowBox_col {m : ℕ} (hm : m < 16384) : LowBox 7 (col m) := by
  intro i
  fin_cases i
  · show m % 128 < 2 ^ 7
    omega
  · show m / 128 < 2 ^ 7
    omega

theorem lowWt_col {m : ℕ} (hm : m < 65536) : lowWt (col m) = popcount 16 m := by
  rw [popcount_eq_wt 16 m hm]
  unfold lowWt
  rw [Fin.sum_univ_two, col_zero, col_one]
  have := wt_add_mul_two_pow 7 (m % 128) (m / 128) (Nat.mod_lt _ (by norm_num))
  rw [show (2 : ℕ) ^ 7 = 128 by norm_num] at this
  rw [← this, Nat.mod_add_div]

variable (K) in
/-- The basis vector of a code. -/
def eCode (m : ℕ) : Lo 2 K := Finsupp.single (col m) 1

variable (K) in
/-- The low element of a sparse list of entries. -/
def toLo (L : List Entry) : Lo 2 K := (L.map fun e => (e.2 : K) • eCode K e.1).sum

theorem toLo_nil : toLo K ([] : List Entry) = 0 := rfl

theorem toLo_cons (e : Entry) (L : List Entry) :
    toLo K (e :: L) = (e.2 : K) • eCode K e.1 + toLo K L := by
  simp [toLo]

theorem toLo_support_subset (L : List Entry) :
    (toLo K L).support ⊆ (L.map fun e => col e.1).toFinset := by
  induction L with
  | nil => simp [toLo_nil]
  | cons e L ih =>
    rw [toLo_cons]
    intro l hl
    rcases Finset.mem_union.mp (Finsupp.support_add hl) with h | h
    · have := Finsupp.support_smul h
      rw [eCode] at this
      have := Finsupp.support_single_subset this
      simp only [Finset.mem_singleton] at this
      subst this
      simp
    · have := ih h
      simp only [List.map_cons, List.toFinset_cons, Finset.mem_insert]
      exact Or.inr this

/-- The meaning of a slot vector, slots starting at free column `f`. -/
def ofSlots (d : Data) : List ℕ → ℕ → Lo 2 K
  | [], _ => 0
  | s :: rest, f => (s : K) • eCode K (d.freeCode f) + ofSlots d rest (f + 1)

theorem ofSlots_zipWith_add (d : Data) :
    ∀ (a b : List ℕ) (f : ℕ), a.length = b.length →
      ofSlots (K := K) d (List.zipWith (· + ·) a b) f = ofSlots d a f + ofSlots d b f
  | [], b, f, hb => by cases b with | nil => simp [ofSlots] | cons => simp at hb
  | _ :: _, [], _, hb => by simp at hb
  | x :: a, y :: b, f, hb => by
    simp only [List.zipWith_cons_cons, ofSlots, Nat.cast_add, add_smul]
    rw [ofSlots_zipWith_add d a b (f + 1) (by simpa using hb)]
    abel

theorem ofSlots_map_mul (d : Data) (x : ℕ) :
    ∀ (s : List ℕ) (f : ℕ), ofSlots (K := K) d (s.map (x * ·)) f = (x : K) • ofSlots d s f
  | [], f => by simp [ofSlots]
  | y :: s, f => by
    simp only [List.map_cons, ofSlots, Nat.cast_mul, mul_smul, smul_add]
    rw [ofSlots_map_mul d x s (f + 1)]

theorem ofSlots_replicate_zero (d : Data) :
    ∀ (n f : ℕ), ofSlots (K := K) d (List.replicate n 0) f = 0
  | 0, _ => rfl
  | n+1, f => by simp [List.replicate_succ, ofSlots, ofSlots_replicate_zero d n (f + 1)]

theorem ofSlots_unitSl (d : Data) :
    ∀ (n g c f : ℕ), g < n →
      ofSlots (K := K) d (unitSl n g c) f = (c : K) • eCode K (d.freeCode (f + g))
  | 0, _, _, _, h => by omega
  | n+1, 0, c, f, _ => by simp [unitSl, ofSlots, ofSlots_replicate_zero]
  | n+1, g+1, c, f, h => by
    simp only [unitSl, ofSlots, Nat.cast_zero, zero_smul, zero_add]
    rw [ofSlots_unitSl d n g c (f + 1) (by omega)]
    congr 3
    omega

theorem ofSlots_eq_zero_of_dvd (d : Data) {p : ℕ} (hp : (p : K) = 0) :
    ∀ (s : List ℕ) (f : ℕ), (∀ x ∈ s, p ∣ x) → ofSlots (K := K) d s f = 0
  | [], _, _ => rfl
  | x :: s, f, h => by
    simp only [ofSlots]
    obtain ⟨q, hq⟩ := h x (List.mem_cons_self)
    rw [ofSlots_eq_zero_of_dvd d hp s (f + 1) (fun y hy => h y (List.mem_cons_of_mem _ hy)), hq,
      Nat.cast_mul, hp, zero_mul, zero_smul, add_zero]

/-- Evaluating a slot vector at a low exponent that is none of its columns gives `0`. -/
theorem ofSlots_apply_eq_zero (d : Data) :
    ∀ (s : List ℕ) (f : ℕ) (l : Fin 2 → ℕ), (∀ g, g < s.length → col (d.freeCode (f + g)) ≠ l) →
      ofSlots (K := K) d s f l = 0
  | [], _, _, _ => rfl
  | x :: s, f, l, h => by
    simp only [ofSlots, Finsupp.add_apply, Finsupp.smul_apply, eCode, Finsupp.single_apply]
    rw [ofSlots_apply_eq_zero d s (f + 1) l (fun g hg => by
      rw [show f + 1 + g = f + (g + 1) by omega]
      exact h (g + 1) (by simpa using hg))]
    rw [if_neg (by simpa using h 0 (by simp)), smul_zero, add_zero]

/-- The free part of row `i`, read from its packed form. -/
def rowFree (d : Data) (i : ℕ) : Lo 2 K := ofSlots d (unpack d.B d.nFree (d.row i)) 0

/-- Row `i`: the pivot with coefficient `1` plus its free part. -/
def rowVec (d : Data) (i : ℕ) : Lo 2 K := eCode K (d.piv i) + rowFree d i

variable (K) in
/-- The subspace of the certificate: the span of all rows. -/
def S (d : Data) : Submodule K (Lo 2 K) :=
  Submodule.span K (Set.range fun i : Fin d.nRows => rowVec d i)

variable (K) in
/-- The span of the rows with pivot weight `≤ 3`. -/
def S3 (d : Data) : Submodule K (Lo 2 K) :=
  Submodule.span K (Set.range fun i : Fin d.low.length => rowVec d (d.low.get i).1)

theorem rowVec_mem_S (d : Data) {i : ℕ} (hi : i < d.nRows) : rowVec d i ∈ S K d :=
  Submodule.subset_span ⟨⟨i, hi⟩, rfl⟩

/-! ## Validity of the data -/

/-- The conjunction of the kernel checks. -/
structure Valid (d : Data) : Prop where
  global : globalOK d = true
  rows : ∀ i < d.nRows, rowOK d i = true
  free : ∀ f < d.nFree, freeOK d f = true
  lowIdx : lowIdxOK d = true
  low : ∀ r ∈ d.low, lowRowOK d r = true
  gen : member d d.gen = true
  nRows_lt : d.nRows < 1024

section Valid

variable {d : Data}

theorem B_pos (d : Data) : 0 < d.B := Nat.two_pow_pos _

theorem bound_lt_B (hd : Valid d) : 512 * d.p * d.p < d.B := by
  have := hd.global
  simp only [globalOK, Bool.and_eq_true, decide_eq_true_eq] at this
  exact this.1.1

theorem tbl_zero (hd : Valid d) : 1024 ≤ d.tblGet 0 := by
  have := hd.global
  simp only [globalOK, Bool.and_eq_true, decide_eq_true_eq] at this
  exact this.1.2

theorem p_ge_three (hd : Valid d) : 3 ≤ d.p := by
  have := hd.global
  simp only [globalOK, Bool.and_eq_true, decide_eq_true_eq] at this
  exact this.2

theorem tbl_piv (hd : Valid d) {i : ℕ} (hi : i < d.nRows) : d.tblGet (d.piv i) = i + 1 := by
  have := hd.rows i hi
  simp only [rowOK, Bool.and_eq_true, beq_iff_eq, decide_eq_true_eq] at this
  exact this.1.1

theorem piv_lt (hd : Valid d) {i : ℕ} (hi : i < d.nRows) : d.piv i < 16384 := by
  have := hd.rows i hi
  simp only [rowOK, Bool.and_eq_true, beq_iff_eq, decide_eq_true_eq] at this
  exact this.1.2

theorem row_slots_lt (hd : Valid d) {i : ℕ} (hi : i < d.nRows) :
    ∀ s ∈ unpack d.B d.nFree (d.row i), s < d.p := by
  have := hd.rows i hi
  simp only [rowOK, Bool.and_eq_true, beq_iff_eq, decide_eq_true_eq] at this
  exact (slotsLt_sound (B_pos d) d.p _ _ this.2).1

theorem row_lt (hd : Valid d) {i : ℕ} (hi : i < d.nRows) : d.row i < d.B ^ d.nFree := by
  have := hd.rows i hi
  simp only [rowOK, Bool.and_eq_true, beq_iff_eq, decide_eq_true_eq] at this
  exact (slotsLt_sound (B_pos d) d.p _ _ this.2).2

theorem tbl_free (hd : Valid d) {f : ℕ} (hf : f < d.nFree) : 1024 ≤ d.tblGet (d.freeCode f) := by
  have := hd.free f hf
  simp only [freeOK, Bool.and_eq_true, decide_eq_true_eq] at this
  exact this.1

theorem freeCode_lt (hd : Valid d) {f : ℕ} (hf : f < d.nFree) : d.freeCode f < 16384 := by
  have := hd.free f hf
  simp only [freeOK, Bool.and_eq_true, decide_eq_true_eq] at this
  exact this.2

theorem piv_injective (hd : Valid d) {i j : ℕ} (hi : i < d.nRows) (hj : j < d.nRows) (h : d.piv i = d.piv j) :
    i = j := by
  have h1 := tbl_piv hd hi
  have h2 := tbl_piv hd hj
  rw [h] at h1
  omega

theorem freeCode_ne_piv (hd : Valid d) {f j : ℕ} (hf : f < d.nFree) (hj : j < d.nRows) : d.freeCode f ≠ d.piv j := by
  intro h
  have h1 := tbl_free hd hf
  have h2 := tbl_piv hd hj
  rw [h] at h1
  have := hd.nRows_lt
  omega

theorem piv_ne_zero (hd : Valid d) {j : ℕ} (hj : j < d.nRows) : d.piv j ≠ 0 := by
  intro h
  have h1 := tbl_piv hd hj
  have h2 := tbl_zero hd
  rw [h] at h1
  have := hd.nRows_lt
  omega

/-- The echelon property: the rows are the identity on the pivot coordinates. -/
theorem rowVec_apply_piv (hd : Valid d) {i j : ℕ} (hi : i < d.nRows) (hj : j < d.nRows) :
    rowVec (K := K) d i (col (d.piv j)) = if i = j then 1 else 0 := by
  unfold rowVec rowFree
  rw [Finsupp.add_apply, ofSlots_apply_eq_zero d _ 0 _ (fun g hg h => ?_), add_zero, eCode,
    Finsupp.single_apply]
  · by_cases hij : i = j
    · subst hij; simp
    · rw [if_neg hij, if_neg]
      intro h
      exact hij (piv_injective hd hi hj (col_injOn h))
  · rw [unpack_length] at hg
    rw [zero_add] at h
    exact freeCode_ne_piv hd hg hj (col_injOn h)

/-- Every element of `S` is the combination of the rows with its pivot coordinates. -/
theorem repr (hd : Valid d) {v : Lo 2 K} (hv : v ∈ S K d) :
    v = ∑ j : Fin d.nRows, v (col (d.piv j)) • rowVec d j := by
  obtain ⟨c, hc⟩ := (Submodule.mem_span_range_iff_exists_fun K).mp hv
  have key : ∀ j : Fin d.nRows, v (col (d.piv j)) = c j := by
    intro j
    rw [← hc, Finsupp.finset_sum_apply]
    simp only [Finsupp.smul_apply, smul_eq_mul]
    rw [Finset.sum_eq_single j]
    · rw [rowVec_apply_piv hd j.isLt j.isLt, if_pos rfl, mul_one]
    · intro i _ hij
      rw [rowVec_apply_piv hd i.isLt j.isLt, if_neg (fun h => hij (Fin.ext h)), mul_zero]
    · intro h; exact absurd (Finset.mem_univ j) h
  calc v = ∑ j : Fin d.nRows, c j • rowVec d j := hc.symm
    _ = ∑ j : Fin d.nRows, v (col (d.piv j)) • rowVec d j :=
      Finset.sum_congr rfl (fun j _ => by rw [key j])

/-- (B1): the constant is not in `S`. -/
theorem one_not_mem (hd : Valid d) : Finsupp.single (0 : Fin 2 → ℕ) (1 : K) ∉ S K d := by
  intro h
  have hr := repr hd h
  have hz : ∀ j : Fin d.nRows, (Finsupp.single (0 : Fin 2 → ℕ) (1 : K)) (col (d.piv j)) = 0 := by
    intro j
    rw [Finsupp.single_apply, if_neg]
    intro h0
    exact piv_ne_zero hd j.isLt ((col_eq_zero_iff _).mp h0.symm)
  simp only [hz, zero_smul, Finset.sum_const_zero] at hr
  exact one_ne_zero (Finsupp.single_eq_zero.mp hr)

/-- Elements of `S` supported on monomials of weight `≤ 3` lie in the span of the low rows. -/
theorem mem_S3_of_lowWt (hd : Valid d) {v : Lo 2 K} (hv : v ∈ S K d) (hw : ∀ l ∈ v.support, lowWt l ≤ 3) :
    v ∈ S3 K d := by
  rw [repr hd hv]
  apply Submodule.sum_mem
  intro j _
  by_cases hp : popcount 16 (d.piv j) ≤ 3
  · apply Submodule.smul_mem
    have hmem : (j : ℕ) ∈ d.low.map Prod.fst := by
      have := hd.lowIdx
      simp only [lowIdxOK, beq_iff_eq] at this
      rw [← this, List.mem_filter, List.mem_range]
      exact ⟨j.isLt, by simpa using hp⟩
    obtain ⟨r, hr, hrj⟩ := List.mem_map.mp hmem
    obtain ⟨n, hn⟩ := List.mem_iff_get.mp hr
    apply Submodule.subset_span
    refine ⟨n, ?_⟩
    simp only [hn, hrj]
  · have h0 : v (col (d.piv j)) = 0 := by
      rw [← Finsupp.not_mem_support_iff]
      intro hs
      have := hw _ hs
      rw [lowWt_col (by have := piv_lt hd j.isLt; omega)] at this
      exact hp this
    rw [h0, zero_smul]
    exact Submodule.zero_mem _

/-! ## Membership -/

/-- The slot vector a free entry contributes. -/
def contribFree (d : Data) (e : Entry) : List ℕ := unitSl d.nFree (d.tblGet e.1 - 1024) e.2

def svFree (d : Data) : List Entry → List ℕ
  | [] => List.replicate d.nFree 0
  | e :: L => List.zipWith (· + ·) (contribFree d e) (svFree d L)

theorem svFree_length (d : Data) : ∀ L, (svFree d L).length = d.nFree
  | [] => by simp [svFree]
  | e :: L => by rw [svFree, List.length_zipWith, contribFree, unitSl_length, svFree_length d L, Nat.min_self]

/-- Unfolding `freeEntryOK`. -/
theorem freeEntryOK_iff (d : Data) (e : Entry) :
    freeEntryOK d e = true ↔ 0 < e.2 ∧ e.2 < d.p ∧ e.1 < 16384 ∧ ∃ v, d.tblGet e.1 = v + 1 ∧
      1024 ≤ v ∧ v + 1 - 1024 < d.nFree ∧ d.freeCode (v + 1 - 1024) = e.1 := by
  unfold freeEntryOK
  rcases hv : d.tblGet e.1 with _ | v
  · simp only [Bool.and_eq_true, decide_eq_true_eq, Bool.false_eq_true, and_false, false_iff,
      not_and, not_exists]
    intro _ _ _ x hx
    omega
  · simp only [Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq]
    constructor
    · rintro ⟨⟨⟨h1, h2⟩, h3⟩, ⟨h4, h5⟩, h6⟩
      exact ⟨h1, h2, h3, v, rfl, h4, h5, h6⟩
    · rintro ⟨h1, h2, h3, v', hv', h4, h5, h6⟩
      obtain rfl : v' = v := by omega
      exact ⟨⟨⟨h1, h2⟩, h3⟩, ⟨h4, h5⟩, h6⟩

theorem packFree_eq (d : Data) : ∀ L : List Entry, (∀ e ∈ L, freeEntryOK d e = true) →
    packFree d L = packSlots d.B (svFree d L)
  | [], _ => by simp [packFree, svFree, packSlots_replicate_zero]
  | e :: L, h => by
    obtain ⟨_, _, _, v, hv, hv1, hv2, _⟩ := (freeEntryOK_iff d e).mp (h e (List.mem_cons_self))
    simp only [packFree, svFree, hv]
    rw [packSlots_zipWith_add _ _ _ (by rw [svFree_length, contribFree, unitSl_length]),
      packFree_eq d L (fun x hx => h x (List.mem_cons_of_mem _ hx)), contribFree, hv,
      packSlots_unitSl _ _ _ _ (by omega), Nat.shiftLeft_eq, Data.B, pow_mul,
      show v + 1 - 1024 = v + 1 - 1024 from rfl]

theorem toLo_eq_ofSlots_svFree (d : Data) : ∀ L : List Entry, (∀ e ∈ L, freeEntryOK d e = true) →
    toLo K L = ofSlots d (svFree d L) 0
  | [], _ => by simp [toLo_nil, svFree, ofSlots_replicate_zero]
  | e :: L, h => by
    obtain ⟨_, _, _, v, hv, hv1, hv2, hv3⟩ := (freeEntryOK_iff d e).mp (h e (List.mem_cons_self))
    rw [toLo_cons, svFree, ofSlots_zipWith_add _ _ _ _ (by rw [svFree_length, contribFree, unitSl_length]),
      toLo_eq_ofSlots_svFree d L (fun x hx => h x (List.mem_cons_of_mem _ hx)), contribFree, hv,
      ofSlots_unitSl _ _ _ _ _ (by omega), zero_add, hv3]

theorem svFree_le (d : Data) : ∀ L : List Entry, (∀ e ∈ L, freeEntryOK d e = true) →
    ∀ s ∈ svFree d L, s ≤ L.length * d.p
  | [], _, s, hs => by simp [svFree] at hs; omega
  | e :: L, h, s, hs => by
    obtain ⟨_, hp, _, _⟩ := (freeEntryOK_iff d e).mp (h e (List.mem_cons_self))
    simp only [svFree] at hs
    have := mem_zipWith_add_le _ _ (mem_unitSl_le _ _ _)
      (svFree_le d L (fun x hx => h x (List.mem_cons_of_mem _ hx))) s hs
    simp only [List.length_cons]
    nlinarith

/-- The slot vector an entry contributes to the echelon accumulator. -/
def contribAcc (d : Data) (e : Entry) : List ℕ :=
  match d.tblGet e.1 with
  | 0 => List.replicate d.nFree 0
  | v+1 => if v < 1024 then (unpack d.B d.nFree (d.row v)).map (e.2 * ·)
      else unitSl d.nFree (v + 1 - 1024) (d.p - e.2)

def svAcc (d : Data) : List Entry → List ℕ
  | [] => List.replicate d.nFree 0
  | e :: L => List.zipWith (· + ·) (contribAcc d e) (svAcc d L)

theorem contribAcc_length (d : Data) (e : Entry) : (contribAcc d e).length = d.nFree := by
  unfold contribAcc
  split
  · simp
  · split <;> simp [unpack_length, unitSl_length]

theorem svAcc_length (d : Data) : ∀ L, (svAcc d L).length = d.nFree
  | [] => by simp [svAcc]
  | e :: L => by rw [svAcc, List.length_zipWith, contribAcc_length, svAcc_length d L, Nat.min_self]

/-- The pivot part of a sparse vector: `∑_{pivot entries} x • rowVec`. -/
def pivPart (d : Data) : List Entry → Lo 2 K
  | [] => 0
  | e :: L => (match d.tblGet e.1 with
      | 0 => 0
      | v+1 => if v < 1024 then (e.2 : K) • rowVec d v else 0) + pivPart d L

theorem entryOK_iff (d : Data) (e : Entry) :
    entryOK d e = true ↔ 0 < e.2 ∧ e.2 < d.p ∧ e.1 < 16384 ∧ ∃ v, d.tblGet e.1 = v + 1 ∧
      ((v < 1024 ∧ v < d.nRows ∧ d.piv v = e.1) ∨
        (1024 ≤ v ∧ v + 1 - 1024 < d.nFree ∧ d.freeCode (v + 1 - 1024) = e.1)) := by
  unfold entryOK
  rcases hv : d.tblGet e.1 with _ | v
  · simp only [Bool.and_eq_true, decide_eq_true_eq, Bool.false_eq_true, and_false, false_iff,
      not_and, not_exists]
    intro _ _ _ x hx
    omega
  · by_cases h : v < 1024
    · simp only [Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq, h, if_true]
      constructor
      · rintro ⟨⟨⟨h1, h2⟩, h3⟩, h4, h5⟩
        exact ⟨h1, h2, h3, v, rfl, Or.inl ⟨h, h4, h5⟩⟩
      · rintro ⟨h1, h2, h3, v', hv', h4 | h4⟩
        · obtain rfl : v' = v := by omega
          exact ⟨⟨⟨h1, h2⟩, h3⟩, h4.2⟩
        · omega
    · simp only [Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq, h, if_false]
      constructor
      · rintro ⟨⟨⟨h1, h2⟩, h3⟩, h4, h5⟩
        exact ⟨h1, h2, h3, v, rfl, Or.inr ⟨by omega, h4, h5⟩⟩
      · rintro ⟨h1, h2, h3, v', hv', h4 | h4⟩
        · omega
        · obtain rfl : v' = v := by omega
          exact ⟨⟨⟨h1, h2⟩, h3⟩, h4.2⟩

theorem accum_eq (hd : Valid d) : ∀ u : List Entry, (∀ e ∈ u, entryOK d e = true) →
    accum d u = packSlots d.B (svAcc d u)
  | [], _ => by simp [accum, svAcc, packSlots_replicate_zero]
  | e :: u, h => by
    obtain ⟨_, _, _, v, hv, hcase⟩ := (entryOK_iff d e).mp (h e List.mem_cons_self)
    rw [accum, svAcc, packSlots_zipWith_add _ _ _ (by rw [contribAcc_length, svAcc_length]),
      accum_eq hd u (fun x hx => h x (List.mem_cons_of_mem _ hx))]
    congr 1
    simp only [contribAcc, hv]
    rcases hcase with ⟨h1, h2, _⟩ | ⟨h1, h2, _⟩
    · rw [if_pos h1, if_pos h1, packSlots_map_mul, packSlots_unpack (B_pos d) _ _ (row_lt hd h2)]
    · rw [if_neg (by omega), if_neg (by omega), packSlots_unitSl _ _ _ _ h2, Nat.shiftLeft_eq, Data.B,
        pow_mul]

theorem svAcc_le (hd : Valid d) : ∀ u : List Entry, (∀ e ∈ u, entryOK d e = true) →
    ∀ s ∈ svAcc d u, s ≤ u.length * (d.p * d.p)
  | [], _, s, hs => by simp [svAcc] at hs; omega
  | e :: u, h, s, hs => by
    obtain ⟨_, hp, _, v, hv, hcase⟩ := (entryOK_iff d e).mp (h e (List.mem_cons_self))
    simp only [svAcc] at hs
    have hc : ∀ x ∈ contribAcc d e, x ≤ d.p * d.p := by
      intro x hx
      simp only [contribAcc, hv] at hx
      rcases hcase with ⟨h1, h2, _⟩ | ⟨h1, _, _⟩
      · rw [if_pos h1, List.mem_map] at hx
        obtain ⟨y, hy, rfl⟩ := hx
        have := row_slots_lt hd h2 y hy
        exact le_of_lt (Nat.mul_lt_mul'' hp this)
      · rw [if_neg (by omega)] at hx
        have := mem_unitSl_le _ _ _ x hx
        have h3 : d.p - e.2 ≤ d.p := Nat.sub_le _ _
        have h4 : d.p ≤ d.p * d.p := Nat.le_mul_self _
        omega
    have := mem_zipWith_add_le _ _ hc (svAcc_le hd u (fun x hx => h x (List.mem_cons_of_mem _ hx))) s hs
    simp only [List.length_cons]
    nlinarith

theorem toLo_add_ofSlots_svAcc (hd : Valid d) {p : ℕ} (hp : p = d.p) (hpK : (p : K) = 0) :
    ∀ u : List Entry, (∀ e ∈ u, entryOK d e = true) →
      toLo K u + ofSlots d (svAcc d u) 0 = pivPart d u
  | [], _ => by simp [toLo_nil, svAcc, ofSlots_replicate_zero, pivPart]
  | e :: u, h => by
    obtain ⟨_, hlt, _, v, hv, hcase⟩ := (entryOK_iff d e).mp (h e List.mem_cons_self)
    have ih := toLo_add_ofSlots_svAcc hd hp hpK u (fun x hx => h x (List.mem_cons_of_mem _ hx))
    rw [toLo_cons, svAcc, ofSlots_zipWith_add _ _ _ _ (by rw [contribAcc_length, svAcc_length]),
      pivPart, ← ih]
    simp only [contribAcc, hv]
    rcases hcase with ⟨h1, h2, h3⟩ | ⟨h1, h2, h3⟩
    · rw [if_pos h1, if_pos h1, ofSlots_map_mul]
      unfold rowVec rowFree
      rw [h3, smul_add]
      abel
    · rw [if_neg (by omega), if_neg (by omega), ofSlots_unitSl _ _ _ _ _ h2, zero_add, h3]
      subst hp
      rw [Nat.cast_sub hlt.le, hpK, zero_sub, neg_smul]
      abel

theorem pivPart_mem (hd : Valid d) : ∀ u : List Entry, (∀ e ∈ u, entryOK d e = true) → pivPart d u ∈ S K d
  | [], _ => Submodule.zero_mem _
  | e :: u, h => by
    obtain ⟨_, _, _, v, hv, hcase⟩ := (entryOK_iff d e).mp (h e (List.mem_cons_self))
    simp only [pivPart, hv]
    apply Submodule.add_mem _ _ (pivPart_mem hd u (fun x hx => h x (List.mem_cons_of_mem _ hx)))
    rcases hcase with ⟨h1, h2, _⟩ | ⟨h1, _, _⟩
    · rw [if_pos h1]
      exact Submodule.smul_mem _ _ (rowVec_mem_S d h2)
    · rw [if_neg (by omega)]
      exact Submodule.zero_mem _

/-- **Membership soundness.** -/
theorem toLo_mem_of_member (hd : Valid d) (hpK : ((d.p : ℕ) : K) = 0) {u : List Entry} (hu : member d u = true) :
    toLo K u ∈ S K d := by
  simp only [member, Bool.and_eq_true, decide_eq_true_eq, List.all_eq_true] at hu
  obtain ⟨⟨hlen, hok⟩, hz⟩ := hu
  have hacc := accum_eq hd u hok
  have hsl : ∀ s ∈ svAcc d u, s < d.B := by
    intro s hs
    have h1 := svAcc_le hd u hok s hs
    have h2 := bound_lt_B hd
    have : u.length * (d.p * d.p) ≤ 512 * d.p * d.p := by nlinarith
    omega
  have hun : unpack d.B d.nFree (accum d u) = svAcc d u := by
    rw [hacc, ← svAcc_length d u, unpack_packSlots (B_pos d) _ hsl]
  have hz' := slotsZero_sound d.B d.p _ _ hz
  rw [hun] at hz'
  have h0 : ofSlots (K := K) d (svAcc d u) 0 = 0 := ofSlots_eq_zero_of_dvd d hpK _ _ hz'
  have key := toLo_add_ofSlots_svAcc hd rfl hpK u hok
  rw [h0, add_zero] at key
  rw [key]
  exact pivPart_mem hd u hok

/-! ## The low rows -/

theorem lowRowOK_iff (d : Data) (r : ℕ × List Entry) :
    lowRowOK d r = true ↔ r.1 < d.nRows ∧ r.2.length ≤ 256 ∧ (∀ e ∈ r.2, freeEntryOK d e = true) ∧
      packFree d r.2 = d.row r.1 ∧ shiftsOK d r = true := by
  unfold lowRowOK
  rw [Bool.and_eq_true, Bool.and_eq_true, Bool.and_eq_true, Bool.and_eq_true, decide_eq_true_eq,
    decide_eq_true_eq, List.all_eq_true, beq_iff_eq]
  tauto

/-- The sparse form of a low row is its row vector. -/
theorem toLo_lowFull (hd : Valid d) {r : ℕ × List Entry} (hr : r ∈ d.low) : toLo K (lowFull d r) = rowVec d r.1 := by
  obtain ⟨hlt, hlen, hok, hpack, _⟩ := (lowRowOK_iff d r).mp (hd.low r hr)
  unfold lowFull rowVec rowFree
  rw [toLo_cons, Nat.cast_one, one_smul, toLo_eq_ofSlots_svFree d _ hok]
  congr 1
  have hsl : ∀ s ∈ svFree d r.2, s < d.B := by
    intro s hs
    have h1 := svFree_le d _ hok s hs
    have h2 := bound_lt_B hd
    have h3 := p_ge_three hd
    have : r.2.length * d.p ≤ 256 * d.p := Nat.mul_le_mul_right _ hlen
    nlinarith
  rw [← hpack, packFree_eq d _ hok, ← svFree_length d r.2, unpack_packSlots (B_pos d) _ hsl]

theorem shifts_mem (hd : Valid d) (hpK : ((d.p : ℕ) : K) = 0) {r : ℕ × List Entry} (hr : r ∈ d.low) {i j : ℕ}
    (hi : i < 2) (hj : j < 7) (c : Bool) :
    toLo K (shiftList i j c (lowFull d r)) ∈ S K d := by
  obtain ⟨_, _, _, _, hs⟩ := (lowRowOK_iff d r).mp (hd.low r hr)
  have := allBelow_sound _ _ hs j hj
  simp only [Bool.and_eq_true] at this
  apply toLo_mem_of_member hd hpK
  rcases Nat.lt_succ_iff_lt_or_eq.mp hi with hi | rfl
  · have : i = 0 := by omega
    subst this
    cases c
    · exact this.1.1.1
    · exact this.1.1.2
  · cases c
    · exact this.1.2
    · exact this.2

end Valid

/-! ## The translation operators on sparse lists -/

variable (K) in
/-- The no-carry part of `u_{i,j}` at depth `8` on a low element (pattern `0` component). -/
def T0 (i : Fin 2) (j : ℕ) : Lo 2 K →ₗ[K] Lo 2 K :=
  (piece 0).comp ((translate 7 8 i j).comp (embed 0))

variable (K) in
/-- The carried part (pattern `e_i` component). -/
def T1 (i : Fin 2) (j : ℕ) : Lo 2 K →ₗ[K] Lo 2 K :=
  (piece (unit i)).comp ((translate 7 8 i j).comp (embed 0))

theorem embed_zero_eCode (m : ℕ) :
    embed (K := K) 0 (eCode K m) = Finsupp.single (col m, (0 : Fin 2 → ℕ)) 1 := by
  unfold embed eCode
  exact Finsupp.mapDomain_single

theorem zero_ne_update (i : Fin 2) : (0 : Fin 2 → ℕ) ≠ Function.update 0 i 1 :=
  fun h => unit_ne_zero i h.symm
theorem col_update_zero_nocarry {m j : ℕ} (h : m % 128 + 2 ^ j < 128) :
    Function.update (col m) 0 (m % 128 + 2 ^ j) = col (m % 128 + 2 ^ j + 128 * (m / 128)) := by
  have := Nat.two_pow_pos j
  funext x
  fin_cases x
  · show m % 128 + 2 ^ j = (m % 128 + 2 ^ j + 128 * (m / 128)) % 128
    omega
  · show m / 128 = (m % 128 + 2 ^ j + 128 * (m / 128)) / 128
    omega
theorem two_pow_le_64 {j : ℕ} (hj : j < 7) : 2 ^ j ≤ 64 :=
  calc 2 ^ j ≤ 2 ^ 6 := Nat.pow_le_pow_right (by norm_num) (by omega)
    _ = 64 := by norm_num

theorem col_update_zero_carry {m j : ℕ} (hj : j < 7) (h : ¬ m % 128 + 2 ^ j < 128) :
    Function.update (col m) 0 (m % 128 + 2 ^ j - 128) =
      col (m % 128 + 2 ^ j - 128 + 128 * (m / 128)) := by
  obtain ⟨s, hs⟩ : ∃ s, m % 128 + 2 ^ j = s + 128 := ⟨m % 128 + 2 ^ j - 128, by omega⟩
  rw [hs, Nat.add_sub_cancel]
  have : m % 128 < 128 := Nat.mod_lt _ (by norm_num)
  have hs' : s < 128 := by have := two_pow_le_64 hj; omega
  funext x
  fin_cases x
  · show s = (s + 128 * (m / 128)) % 128
    omega
  · show m / 128 = (s + 128 * (m / 128)) / 128
    omega
theorem col_update_one_nocarry {m j : ℕ} (_h : m / 128 + 2 ^ j < 128) :
    Function.update (col m) 1 (m / 128 + 2 ^ j) = col (m % 128 + 128 * (m / 128 + 2 ^ j)) := by
  have := Nat.two_pow_pos j
  funext x
  fin_cases x
  · show m % 128 = (m % 128 + 128 * (m / 128 + 2 ^ j)) % 128
    omega
  · show m / 128 + 2 ^ j = (m % 128 + 128 * (m / 128 + 2 ^ j)) / 128
    omega
theorem col_update_one_carry {m j : ℕ} (h : ¬ m / 128 + 2 ^ j < 128) :
    Function.update (col m) 1 (m / 128 + 2 ^ j - 128) =
      col (m % 128 + 128 * (m / 128 + 2 ^ j - 128)) := by
  obtain ⟨s, hs⟩ : ∃ s, m / 128 + 2 ^ j = s + 128 := ⟨m / 128 + 2 ^ j - 128, by omega⟩
  rw [hs, Nat.add_sub_cancel]
  funext x
  fin_cases x
  · show m % 128 = (m % 128 + 128 * s) % 128
    omega
  · show s = (m % 128 + 128 * s) / 128
    omega
theorem toLo_option (o : Option Entry) :
    toLo K o.toList = match o with | none => 0 | some e => (e.2 : K) • eCode K e.1 := by
  cases o <;> simp [toLo_nil, toLo_cons]

theorem T_eCode_zero {j : ℕ} (hj : j < 7) (m c : ℕ) :
    T0 K 0 j ((c : K) • eCode K m) = toLo K (shiftEntry 0 j false (m, c)).toList ∧
    T1 K 0 j ((c : K) • eCode K m) = toLo K (shiftEntry 0 j true (m, c)).toList := by
  have hstep : ∀ a : Ex 2, translate (K := K) 7 8 (0 : Fin 2) j (Finsupp.single a (1 : K)) =
      Finsupp.single (step 7 8 (0 : Fin 2) j a) 1 :=
    fun a => ICFallDegrees.ChainModelMap.translate_single (0 : Fin 2) j a 1
  have h7 : (2 : ℕ) ^ 7 = 128 := by norm_num
  have h1 : (2 : ℕ) ^ (8 - 7) = 2 := by norm_num
  have hne : ((0 : Fin 2 → ℕ) = Function.update (0 : Fin 2 → ℕ) (0 : Fin 2) (1 : ℕ)) = False :=
    eq_false (zero_ne_update 0)
  have hne' : (Function.update (0 : Fin 2 → ℕ) (0 : Fin 2) (1 : ℕ) = (0 : Fin 2 → ℕ)) = False :=
    eq_false (zero_ne_update 0).symm
  simp only [T0, T1, LinearMap.comp_apply, map_smul, embed_zero_eCode, hstep, step, if_pos hj,
    stepLow, toLo_option, shiftEntry, unit, Pi.zero_apply, zero_add, h7, h1, Nat.one_mod, col_zero,
    if_true]
  by_cases hc : m % 128 + 2 ^ j < 128
  · simp only [if_pos hc, piece_single, hne, hne', if_true, if_false, Bool.false_eq_true, eCode, smul_zero,
      col_update_zero_nocarry hc]
    exact ⟨trivial, trivial⟩
  · simp only [if_neg hc, piece_single, hne, hne', if_true, if_false, Bool.false_eq_true, eCode, smul_zero,
      col_update_zero_carry hj hc]
    exact ⟨trivial, trivial⟩

theorem T_eCode_one {j : ℕ} (hj : j < 7) (m c : ℕ) :
    T0 K 1 j ((c : K) • eCode K m) = toLo K (shiftEntry 1 j false (m, c)).toList ∧
    T1 K 1 j ((c : K) • eCode K m) = toLo K (shiftEntry 1 j true (m, c)).toList := by
  have hstep : ∀ a : Ex 2, translate (K := K) 7 8 (1 : Fin 2) j (Finsupp.single a (1 : K)) =
      Finsupp.single (step 7 8 (1 : Fin 2) j a) 1 :=
    fun a => ICFallDegrees.ChainModelMap.translate_single (1 : Fin 2) j a 1
  have h7 : (2 : ℕ) ^ 7 = 128 := by norm_num
  have h1 : (2 : ℕ) ^ (8 - 7) = 2 := by norm_num
  have hne : ((0 : Fin 2 → ℕ) = Function.update (0 : Fin 2 → ℕ) (1 : Fin 2) (1 : ℕ)) = False :=
    eq_false (zero_ne_update 1)
  have hne' : (Function.update (0 : Fin 2 → ℕ) (1 : Fin 2) (1 : ℕ) = (0 : Fin 2 → ℕ)) = False :=
    eq_false (zero_ne_update 1).symm
  have h10 : ((1 : ℕ) = 0) = False := eq_false one_ne_zero
  simp only [T0, T1, LinearMap.comp_apply, map_smul, embed_zero_eCode, hstep, step, if_pos hj,
    stepLow, toLo_option, shiftEntry, unit, Pi.zero_apply, zero_add, h7, h1, Nat.one_mod, col_one,
    h10, if_false]
  by_cases hc : m / 128 + 2 ^ j < 128
  · simp only [if_pos hc, piece_single, hne, hne', if_true, if_false, Bool.false_eq_true, eCode, smul_zero,
      col_update_one_nocarry hc]
    exact ⟨trivial, trivial⟩
  · simp only [if_neg hc, piece_single, hne, hne', if_true, if_false, Bool.false_eq_true, eCode, smul_zero,
      col_update_one_carry hc]
    exact ⟨trivial, trivial⟩

theorem T_eCode {i : Fin 2} {j : ℕ} (hj : j < 7) (m c : ℕ) :
    T0 K i j ((c : K) • eCode K m) = toLo K (shiftEntry i j false (m, c)).toList ∧
    T1 K i j ((c : K) • eCode K m) = toLo K (shiftEntry i j true (m, c)).toList := by
  fin_cases i
  · exact T_eCode_zero hj m c
  · exact T_eCode_one hj m c

theorem T_toLo {i : Fin 2} {j : ℕ} (hj : j < 7) :
    ∀ L : List Entry,
      T0 K i j (toLo K L) = toLo K (shiftList i j false L) ∧
      T1 K i j (toLo K L) = toLo K (shiftList i j true L)
  | [] => by simp [toLo_nil, shiftList]
  | e :: L => by
    obtain ⟨ih0, ih1⟩ := T_toLo (i := i) hj L
    obtain ⟨h0, h1⟩ := T_eCode (K := K) (i := i) hj e.1 e.2
    rw [toLo_cons, map_add, map_add, ih0, ih1, h0, h1]
    unfold shiftList
    rw [List.filterMap_cons, List.filterMap_cons]
    constructor
    · cases shiftEntry (i : ℕ) j false (e.1, e.2) <;> simp [toLo_nil, toLo_cons]
    · cases shiftEntry (i : ℕ) j true (e.1, e.2) <;> simp [toLo_nil, toLo_cons]

/-! ## (B2): closedness from the finite checks -/

/-- Closedness of the free extension at depth `8` from: `S ∩ {weight ≤ 3} ⊆ S₃`, and `T0`, `T1`
map `S₃` into `S`.  The proof follows `DepthReduction.W_closed`. -/
theorem closed_of (S S3 : Submodule K (Lo 2 K))
    (hS3 : ∀ v ∈ S, (∀ l ∈ v.support, lowWt l ≤ 3) → v ∈ S3)
    (hT : ∀ (i : Fin 2) (j : ℕ), j < 7 → ∀ v ∈ S3, T0 K i j v ∈ S ∧ T1 K i j v ∈ S) :
    Closed 7 8 4 (W 7 4 S 8) := by
  classical
  intro w hw hlow i j hj
  have hk' : 7 < 8 := by norm_num
  have hbox := hw.1
  refine ⟨translate_boxed hbox, fun q => ⟨?_, fun l hl => ?_⟩⟩
  · by_cases hq : LowBox (8 - 7) q
    · have hpiece : ∀ q', piece q' w ∈ S3 := by
        intro q'
        apply hS3 _ (hw.2 q').1
        intro l hl
        rw [Finsupp.mem_support_iff, piece_apply] at hl
        have := hlow (l, q') (Finsupp.mem_support_iff.mpr hl)
        simp only [pwt] at this
        omega
      by_cases hjl : j < 7
      · rw [piece_translate_low hk' hj hjl hbox hq]
        apply S.add_mem
        · exact (hT i j hjl _ (hpiece q)).1
        · exact (hT i j hjl _ (hpiece (qPrev 7 8 i q))).2
      · rw [piece_translate_high hk' hj (not_lt.mp hjl) hbox hq]
        exact (hw.2 _).1
    · have : piece q (translate (K := K) 7 8 i j w) = 0 := by
        ext l
        rw [piece_apply, translate_apply_of_not_inBox hbox (fun h => hq h.2), Finsupp.zero_apply]
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

section Valid2

variable {d : Data}

theorem T_S3_le (hd : Valid d) (hpK : ((d.p : ℕ) : K) = 0) (i : Fin 2) {j : ℕ} (hj : j < 7) :
    ∀ v ∈ S3 K d, T0 K i j v ∈ S K d ∧ T1 K i j v ∈ S K d := by
  have key : ∀ (T : Lo 2 K →ₗ[K] Lo 2 K), (∀ r ∈ d.low, T (rowVec d r.1) ∈ S K d) →
      ∀ v ∈ S3 K d, T v ∈ S K d := by
    intro T hT v hv
    have : (S3 K d).map T ≤ S K d := by
      rw [S3, Submodule.map_span_le]
      rintro _ ⟨n, rfl⟩
      exact hT _ (List.get_mem _ _)
    exact this ⟨v, hv, rfl⟩
  intro v hv
  constructor
  · refine key (T0 K i j) (fun r hr => ?_) v hv
    rw [← toLo_lowFull hd hr, (T_toLo (i := i) hj _).1]
    exact shifts_mem hd hpK hr i.isLt hj false
  · refine key (T1 K i j) (fun r hr => ?_) v hv
    rw [← toLo_lowFull hd hr, (T_toLo (i := i) hj _).2]
    exact shifts_mem hd hpK hr i.isLt hj true

/-- (B2). -/
theorem closed (hd : Valid d) (hpK : ((d.p : ℕ) : K) = 0) : Closed 7 8 4 (W 7 4 (S K d) 8) :=
  closed_of (S K d) (S3 K d) (fun _ hv hw => mem_S3_of_lowWt hd hv hw)
    (fun i _ hj => T_S3_le hd hpK i hj)

/-- (G), first half: the generator list is in `S`. -/
theorem gen_mem (hd : Valid d) (hpK : ((d.p : ℕ) : K) = 0) : toLo K d.gen ∈ S K d :=
  toLo_mem_of_member hd hpK hd.gen

end Valid2

/-! ## (G): the model map on sparse generator lists -/

open MvPolynomial ChainTranslation ChainQuotient ChainModelMap in
theorem Φ_embed_eCode (k : ℕ) (m : ℕ) :
    Φ (K := K) 7 k (embed 0 (eCode K m)) = monomial (digits k (col m)) 1 := by
  rw [embed_zero_eCode, Φ_single]
  have : encode 7 (col m, (0 : Fin 2 → ℕ)) = col m := by
    funext i
    simp [encode]
  rw [this]

open ChainTranslation ChainQuotient in
theorem bit_of_lt {n j : ℕ} (hn : n < 4) (hj : 2 ≤ j) : bit n j = 0 := by
  unfold bit
  rw [Nat.div_eq_of_lt (lt_of_lt_of_le hn (Nat.pow_le_pow_right (by norm_num : 0 < 2) hj))]

open ChainTranslation ChainQuotient in
/-- The chain exponent of a two-bit exponent vector, as four singles. -/
theorem digits_two_bits {k : ℕ} (hk : 2 ≤ k) (b : Fin 2 → ℕ) (h0 : b 0 < 4) (h1 : b 1 < 4) :
    digits k b =
      Finsupp.single (var 2 k 0 ⟨0, by omega⟩) (bit (b 0) 0) +
      Finsupp.single (var 2 k 0 ⟨1, by omega⟩) (bit (b 0) 1) +
      Finsupp.single (var 2 k 1 ⟨0, by omega⟩) (bit (b 1) 0) +
      Finsupp.single (var 2 k 1 ⟨1, by omega⟩) (bit (b 1) 1) := by
  classical
  ext v
  simp only [digits_apply, Finsupp.add_apply, Finsupp.single_apply, Fin.ext_iff, var_val,
    Fin.val_zero, Fin.val_one, mul_zero, mul_one, add_zero, zero_add]
  have hv := v.isLt
  have hkpos : 0 < k := by omega
  have hrow : (row (m := 2) (k := k) v : ℕ) = v / k := row_val v
  have hdep : (depth (m := 2) (k := k) v : ℕ) = v % k := depth_val v
  have hb : ∀ i : Fin 2, b i < 4 := by intro i; fin_cases i <;> assumption
  have hk0 : k ≠ 0 := by omega
  have hk1 : 1 + k ≠ 0 := by omega
  have hk1' : k ≠ 1 := by omega
  have hk2 : 1 + k ≠ 1 := by omega
  have hkk : 1 + k ≠ k := by omega
  have hkk' : k ≠ 1 + k := by omega
  by_cases c0 : (v : ℕ) = 0
  · have hr : row (m := 2) (k := k) v = 0 := Fin.ext (by rw [hrow, Fin.val_zero, c0, Nat.zero_div])
    have hd : (depth (m := 2) (k := k) v : ℕ) = 0 := by rw [hdep, c0, Nat.zero_mod]
    rw [hr, hd]
    simp [c0, hk0, hk0.symm, hk1, hk1.symm, hk1', hk1'.symm, hk2, hk2.symm, hkk, hkk']
  by_cases c1 : (v : ℕ) = 1
  · have hr : row (m := 2) (k := k) v = 0 :=
      Fin.ext (by rw [hrow, Fin.val_zero, c1, Nat.div_eq_of_lt (by omega)])
    have hd : (depth (m := 2) (k := k) v : ℕ) = 1 := by rw [hdep, c1, Nat.mod_eq_of_lt (by omega)]
    rw [hr, hd]
    simp [c1, hk0, hk0.symm, hk1, hk1.symm, hk1', hk1'.symm, hk2, hk2.symm, hkk, hkk']
  by_cases ck : (v : ℕ) = k
  · have hr : row (m := 2) (k := k) v = 1 := Fin.ext (by rw [hrow, Fin.val_one, ck, Nat.div_self hkpos])
    have hd : (depth (m := 2) (k := k) v : ℕ) = 0 := by rw [hdep, ck, Nat.mod_self]
    rw [hr, hd]
    simp [ck, hk0, hk0.symm, hk1, hk1.symm, hk1', hk1'.symm, hk2, hk2.symm, hkk, hkk']
  by_cases ck1 : (v : ℕ) = 1 + k
  · have hr : row (m := 2) (k := k) v = 1 :=
      Fin.ext (by rw [hrow, Fin.val_one, ck1, Nat.add_div_right _ hkpos, Nat.div_eq_of_lt (by omega)])
    have hd : (depth (m := 2) (k := k) v : ℕ) = 1 := by
      rw [hdep, ck1, Nat.add_mod_right, Nat.mod_eq_of_lt (by omega)]
    rw [hr, hd]
    simp [ck1, hk0, hk0.symm, hk1, hk1.symm, hk1', hk1'.symm, hk2, hk2.symm, hkk, hkk']
  · have hd : 2 ≤ (depth (m := 2) (k := k) v : ℕ) := by
      rw [hdep]
      have hq : (v : ℕ) / k < 2 := (Nat.div_lt_iff_lt_mul hkpos).mpr (by omega)
      have hdm := Nat.div_add_mod (v : ℕ) k
      by_contra hlt
      push_neg at hlt
      interval_cases hq' : (v : ℕ) / k <;> interval_cases hm : (v : ℕ) % k <;> omega
    rw [bit_of_lt (hb _) hd]
    simp [Ne.symm c0, Ne.symm c1, Ne.symm ck, Ne.symm ck1]

end ICFallDegrees.N5CertSound

end
