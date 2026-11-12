import ICFallDegrees.N5CertSound
import ICFallDegrees.ShiftWitness

/-!
# N7 certificate: a low ideal element of weight `≤ 4` outside the base `S`

Extension of the N5 base certificate (`N5Cert`, `N5CertSound`) by one more object: a cofactor
`h = ∑ c_e x₀^{e₀} x₁^{e₁}` (the list `terms`) such that the product `u = h · g` with the
generator `g` (computed without any wrap: every shifted monomial stays in the low box
`[0,128)²`) has all its monomials of weight `≤ 4` and is **not** in the row space `S` of the
certified echelon basis.

Checker (plain `Nat`/`List` code, structural recursion only; evaluated by `decide +kernel`):

* `prodOK d F terms q`: every term keeps the generator inside the low box (`termOK`), and the
  sparse product `∑_terms c · shift_e(gen)`, accumulated by sorted merges (`addP`, `mulAcc`) and
  reduced modulo `p` (`normalize`), is literally the list `q`; the cofactor is cut into chunks,
  one such check per chunk, so that each kernel evaluation stays small;
* `sumOK p F parts u`: `u` is the sum of the partial products modulo `p`;
* `outsideOK d u f`: every entry of `u` is an admissible column of weight `≤ 4`, the echelon
  accumulator of `u` cannot overflow a slot (`coefSum u * p < B`), the free column `f` occurs
  only once among the free columns, and slot `f` of the accumulator is not divisible by `p`:
  the echelon residual of `u` is nonzero in free column `f`.

Soundness: `witness_lowIdeal`, `sum_lowIdeal` (`toLo u` is a `ShiftWitness.LowIdeal` element of
the generator),
`support_of_outsideOK` (low box and weight `≤ 4`), `not_mem_of_outsideOK` (`toLo u ∉ S`).
-/

namespace ICFallDegrees.N5CertWitness

open ICFallDegrees.N5Cert

/-- A cofactor term `(e₀, e₁, c)`: the monomial `c · x₀^{e₀} x₁^{e₁}`. -/
abbrev Term := Nat × Nat × Nat

/-- `c · x₀^{e₀} x₁^{e₁} · L` on sparse lists (no reduction modulo `p`, no wrap). -/
def shiftScale (t : Term) (L : List Entry) : List Entry :=
  L.map fun e => (e.1 + t.1 + 128 * t.2.1, t.2.2 * e.2)

/-- Sum of two sparse lists by a merge on the codes (fuel `n`; any fuel is sound). -/
def addP : Nat → List Entry → List Entry → List Entry
  | 0, A, B => A ++ B
  | _+1, [], B => B
  | _+1, a :: A, [] => a :: A
  | n+1, a :: A, b :: B =>
    if a.1 < b.1 then a :: addP n A (b :: B)
    else if b.1 < a.1 then b :: addP n (a :: A) B
    else (a.1, a.2 + b.2) :: addP n A B

/-- The product `(∑ terms) · g`, accumulated term by term. -/
def mulAcc (F : Nat) (g : List Entry) : List Term → List Entry
  | [] => []
  | t :: ts => addP F (shiftScale t g) (mulAcc F g ts)

/-- Reduce the coefficients modulo `p` and drop the zero ones. -/
def normalize (p : Nat) : List Entry → List Entry
  | [] => []
  | e :: L => if e.2 % p = 0 then normalize p L else (e.1, e.2 % p) :: normalize p L

/-- The term keeps every monomial of `g` inside the low box `[0,128)²`. -/
def termOK (g : List Entry) (t : Term) : Bool :=
  g.all fun e => e.1 % 128 + t.1 < 128 && e.1 / 128 + t.2.1 < 128

def eqList : List Entry → List Entry → Bool
  | [], [] => true
  | a :: A, b :: B => a.1 == b.1 && a.2 == b.2 && eqList A B
  | _, _ => false

def coefSum : List Entry → Nat
  | [] => 0
  | e :: L => e.2 + coefSum L

/-- `u` is the product of the cofactor `terms` with the generator of `d`, modulo `p`. -/
def prodOK (d : Data) (F : Nat) (terms : List Term) (u : List Entry) : Bool :=
  terms.all (termOK d.gen) && eqList (normalize d.p (mulAcc F d.gen terms)) u

/-- The sum of a list of sparse lists. -/
def sumParts (F : Nat) : List (List Entry) → List Entry
  | [] => []
  | q :: qs => addP F q (sumParts F qs)

/-- `u` is the sum of the partial products, modulo `p`. -/
def sumOK (p F : Nat) (parts : List (List Entry)) (u : List Entry) : Bool :=
  eqList (normalize p (sumParts F parts)) u

/-- `u` has weight `≤ 4` and a nonzero echelon residual in free column `f`. -/
def outsideOK (d : Data) (u : List Entry) (f : Nat) : Bool :=
  u.all (fun e => entryOK d e && decide (popcount 16 e.1 ≤ 4)) &&
    decide (coefSum u * d.p < d.B) && decide (f < d.nFree) &&
    allBelow (fun f' => f' == f || !(d.freeCode f' == d.freeCode f)) d.nFree &&
    !((unpack d.B d.nFree (accum d u)).getD f 0 % d.p == 0)

end ICFallDegrees.N5CertWitness

/-! ## Soundness -/

noncomputable section

namespace ICFallDegrees.N5CertWitness

open ICFallDegrees.N5Cert ICFallDegrees.N5CertSound DepthReduction BinaryWeight
open ICFallDegrees.ShiftWitness

variable {K : Type*} [Field K]

theorem toLo_append (A B : List Entry) : toLo K (A ++ B) = toLo K A + toLo K B := by
  simp [toLo, List.map_append, List.sum_append]

theorem toLo_addP : ∀ (n : ℕ) (A B : List Entry), toLo K (addP n A B) = toLo K A + toLo K B
  | 0, A, B => by rw [addP, toLo_append]
  | n+1, [], B => by rw [addP, toLo_nil, zero_add]
  | n+1, a :: A, [] => by rw [addP, toLo_nil, add_zero]
  | n+1, a :: A, b :: B => by
    rw [addP]
    split_ifs with h1 h2
    · rw [toLo_cons, toLo_addP n A (b :: B), toLo_cons a A]
      abel
    · rw [toLo_cons, toLo_addP n (a :: A) B, toLo_cons b B]
      abel
    · have hab : b.1 = a.1 := by omega
      rw [toLo_cons, toLo_addP n A B, toLo_cons a A, toLo_cons b B, hab]
      simp only [Nat.cast_add, add_smul]
      abel

theorem cast_mod {p : ℕ} (hp : (p : K) = 0) (x : ℕ) : ((x % p : ℕ) : K) = (x : K) := by
  conv_rhs => rw [← Nat.div_add_mod x p]
  rw [Nat.cast_add, Nat.cast_mul, hp, zero_mul, zero_add]

theorem toLo_normalize {p : ℕ} (hp : (p : K) = 0) :
    ∀ L : List Entry, toLo K (normalize p L) = toLo K L
  | [] => rfl
  | e :: L => by
    rw [normalize]
    split_ifs with h
    · rw [toLo_normalize hp L, toLo_cons, ← cast_mod hp e.2, h, Nat.cast_zero, zero_smul, zero_add]
    · rw [toLo_cons, toLo_cons, toLo_normalize hp L, cast_mod hp]

theorem eqList_sound : ∀ A B : List Entry, eqList A B = true → A = B
  | [], [], _ => rfl
  | [], _ :: _, h => by simp [eqList] at h
  | _ :: _, [], h => by simp [eqList] at h
  | a :: A, b :: B, h => by
    simp only [eqList, Bool.and_eq_true, beq_iff_eq] at h
    rw [eqList_sound A B h.2, Prod.ext h.1.1 h.1.2]

theorem col_add {m e0 e1 : ℕ} (h : m % 128 + e0 < 128) :
    col (m + e0 + 128 * e1) = col m + ![e0, e1] := by
  funext i
  fin_cases i
  · show (m + e0 + 128 * e1) % 128 = m % 128 + e0
    omega
  · show (m + e0 + 128 * e1) / 128 = m / 128 + e1
    omega

theorem toLo_shiftScale (t : Term) :
    ∀ L : List Entry, (∀ e ∈ L, e.1 % 128 + t.1 < 128) →
      toLo K (shiftScale t L) = (t.2.2 : K) • lowShift ![t.1, t.2.1] (toLo K L)
  | [], _ => by simp [shiftScale, toLo_nil]
  | e :: L, h => by
    have ih := toLo_shiftScale t L (fun x hx => h x (List.mem_cons_of_mem _ hx))
    have he := h e List.mem_cons_self
    unfold shiftScale at ih ⊢
    rw [List.map_cons, toLo_cons, ih, toLo_cons, map_add, smul_add, map_smul]
    congr 1
    simp only [eCode]
    rw [lowShift_single, col_add he, Nat.cast_mul, mul_smul]

theorem vec_eq (e0 e1 : ℕ) :
    (![e0, e1] : Fin 2 → ℕ) = e0 • Pi.single 0 1 + e1 • Pi.single 1 1 := by
  funext i
  fin_cases i <;> simp

/-- Shifting by `(e₀, e₁)` inside the low box preserves low ideal elements. -/
theorem lowIdeal_shift2 {g v : Lo 2 K} (h : LowIdeal 7 g v) (e0 e1 : ℕ)
    (hs : ∀ l ∈ v.support, l 0 + e0 < 128 ∧ l 1 + e1 < 128) :
    LowIdeal 7 g (lowShift ![e0, e1] v) := by
  rw [vec_eq, lowShift_add]
  refine LowIdeal.shift_nsmul 1 e1 (LowIdeal.shift_nsmul 0 e0 h (fun l hl => ?_)) ?_
  · have := (hs l hl).1
    norm_num
    exact this
  · intro l' hl'
    obtain ⟨l, hl, rfl⟩ := lowShift_support _ _ hl'
    have := (hs l hl).2
    simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul]
    norm_num
    exact this

theorem termOK_bounds {g : List Entry} {t : Term} (h : termOK g t = true) :
    ∀ e ∈ g, e.1 % 128 + t.1 < 128 ∧ e.1 / 128 + t.2.1 < 128 := by
  intro e he
  simp only [termOK, List.all_eq_true, Bool.and_eq_true, decide_eq_true_eq] at h
  exact h e he

theorem toLo_mulAcc_lowIdeal (F : ℕ) (g : List Entry) :
    ∀ terms : List Term, (∀ t ∈ terms, termOK g t = true) →
      LowIdeal 7 (toLo K g) (toLo K (mulAcc F g terms))
  | [], _ => by rw [mulAcc, toLo_nil]; exact LowIdeal.zero
  | t :: ts, h => by
    have ht := termOK_bounds (h t List.mem_cons_self)
    rw [mulAcc, toLo_addP, toLo_shiftScale t g (fun e he => (ht e he).1)]
    refine LowIdeal.add (LowIdeal.smul _ (lowIdeal_shift2 LowIdeal.gen _ _ ?_))
      (toLo_mulAcc_lowIdeal F g ts (fun x hx => h x (List.mem_cons_of_mem _ hx)))
    intro l hl
    have hmem := toLo_support_subset g hl
    rw [List.mem_toFinset, List.mem_map] at hmem
    obtain ⟨e, he, rfl⟩ := hmem
    exact ht e he

variable {d : Data}

/-- **Product soundness.**  The checked list is a low ideal element of the generator. -/
theorem witness_lowIdeal (hpK : ((d.p : ℕ) : K) = 0) {F : ℕ} {terms : List Term} {u : List Entry}
    (h : prodOK d F terms u = true) : LowIdeal 7 (toLo K d.gen) (toLo K u) := by
  simp only [prodOK, Bool.and_eq_true, List.all_eq_true] at h
  rw [← eqList_sound _ _ h.2, toLo_normalize hpK]
  exact toLo_mulAcc_lowIdeal F d.gen terms h.1

theorem all_nil {α : Type*} {P : α → Prop} : ∀ q ∈ ([] : List α), P q :=
  fun _ h => absurd h List.not_mem_nil

theorem all_cons {α : Type*} {P : α → Prop} {a : α} {L : List α} (ha : P a)
    (hL : ∀ q ∈ L, P q) : ∀ q ∈ a :: L, P q := by
  intro q hq
  rcases List.mem_cons.mp hq with rfl | h
  · exact ha
  · exact hL q h

/-- **Sum soundness.**  A checked sum of low ideal elements is a low ideal element. -/
theorem sum_lowIdeal {p : ℕ} (hpK : (p : K) = 0) {F : ℕ} {g : Lo 2 K} {parts : List (List Entry)}
    {u : List Entry} (h : ∀ q ∈ parts, LowIdeal 7 g (toLo K q))
    (hs : sumOK p F parts u = true) : LowIdeal 7 g (toLo K u) := by
  unfold sumOK at hs
  rw [← eqList_sound _ _ hs, toLo_normalize hpK]
  clear hs
  induction parts with
  | nil => rw [sumParts, toLo_nil]; exact LowIdeal.zero
  | cons q qs ih =>
    rw [sumParts, toLo_addP]
    exact LowIdeal.add (h q List.mem_cons_self) (ih (fun x hx => h x (List.mem_cons_of_mem _ hx)))

theorem outsideOK_spec {u : List Entry} {f : ℕ} (h : outsideOK d u f = true) :
    (∀ e ∈ u, entryOK d e = true ∧ popcount 16 e.1 ≤ 4) ∧ coefSum u * d.p < d.B ∧
      f < d.nFree ∧ (∀ f' < d.nFree, f' ≠ f → d.freeCode f' ≠ d.freeCode f) ∧
      (unpack d.B d.nFree (accum d u)).getD f 0 % d.p ≠ 0 := by
  unfold outsideOK at h
  simp only [Bool.and_eq_true, List.all_eq_true, decide_eq_true_eq, Bool.not_eq_true',
    beq_eq_false_iff_ne, ne_eq] at h
  obtain ⟨⟨⟨⟨h1, h2⟩, h3⟩, h4⟩, h5⟩ := h
  refine ⟨h1, h2, h3, fun f' hf' hne heq => ?_, h5⟩
  have := allBelow_sound _ _ h4 f' hf'
  simp only [Bool.or_eq_true, beq_iff_eq, Bool.not_eq_true', beq_eq_false_iff_ne, ne_eq] at this
  rcases this with h | h
  · exact hne h
  · exact h heq

/-- Low box and weight of the checked list. -/
theorem support_of_outsideOK {u : List Entry} {f : ℕ} (h : outsideOK d u f = true) :
    ∀ l ∈ (toLo K u).support, LowBox 7 l ∧ lowWt l ≤ 4 := by
  obtain ⟨h1, -⟩ := outsideOK_spec h
  intro l hl
  have hmem := toLo_support_subset u hl
  rw [List.mem_toFinset, List.mem_map] at hmem
  obtain ⟨e, he, rfl⟩ := hmem
  obtain ⟨hok, hwt⟩ := h1 e he
  obtain ⟨_, _, hlt, _⟩ := (entryOK_iff d e).mp hok
  exact ⟨lowBox_col hlt, by rw [lowWt_col (by omega)]; exact hwt⟩

/-- A sharper slot bound for the echelon accumulator: `(∑ coefficients) · p`. -/
theorem svAcc_le' (hd : Valid d) : ∀ u : List Entry, (∀ e ∈ u, entryOK d e = true) →
    ∀ s ∈ svAcc d u, s ≤ coefSum u * d.p
  | [], _, s, hs => by simp [svAcc] at hs; omega
  | e :: u, h, s, hs => by
    obtain ⟨hpos, hp, _, v, hv, hcase⟩ := (entryOK_iff d e).mp (h e (List.mem_cons_self))
    simp only [svAcc] at hs
    have hc : ∀ x ∈ contribAcc d e, x ≤ e.2 * d.p := by
      intro x hx
      simp only [contribAcc, hv] at hx
      rcases hcase with ⟨h1, h2, _⟩ | ⟨h1, _, _⟩
      · rw [if_pos h1, List.mem_map] at hx
        obtain ⟨y, hy, rfl⟩ := hx
        have := row_slots_lt hd h2 y hy
        exact Nat.mul_le_mul_left _ this.le
      · rw [if_neg (by omega)] at hx
        have := mem_unitSl_le _ _ _ x hx
        have h3 : d.p - e.2 ≤ d.p := Nat.sub_le _ _
        have h4 : d.p ≤ e.2 * d.p := Nat.le_mul_of_pos_left _ hpos
        omega
    have := mem_zipWith_add_le _ _ hc
      (svAcc_le' hd u (fun x hx => h x (List.mem_cons_of_mem _ hx))) s hs
    simp only [coefSum]
    rw [Nat.add_mul]
    exact this

/-- Evaluating a slot vector at a column that occurs only once among its columns. -/
theorem ofSlots_apply_unique (d : Data) :
    ∀ (s : List ℕ) (f0 g : ℕ), g < s.length →
      (∀ g' < s.length, g' ≠ g → d.freeCode (f0 + g') ≠ d.freeCode (f0 + g)) →
      ofSlots (K := K) d s f0 (col (d.freeCode (f0 + g))) = ((s.getD g 0 : ℕ) : K)
  | [], _, _, hg, _ => by simp at hg
  | x :: s, f0, 0, _, hu => by
    simp only [ofSlots, Finsupp.add_apply, Finsupp.smul_apply, eCode, add_zero,
      Finsupp.single_eq_same, smul_eq_mul, mul_one, List.getD_cons_zero]
    rw [ofSlots_apply_eq_zero d s (f0 + 1) _ (fun g' hg' hcol => ?_), add_zero]
    have := hu (g' + 1) (by simpa using hg') (by omega)
    rw [add_zero, ← add_assoc, add_right_comm] at this
    exact this (col_injOn hcol)
  | x :: s, f0, g + 1, hg, hu => by
    simp only [ofSlots, Finsupp.add_apply, Finsupp.smul_apply, eCode, List.getD_cons_succ]
    have hne : col (d.freeCode f0) ≠ col (d.freeCode (f0 + (g + 1))) := by
      intro hcol
      have := hu 0 (by simp) (by omega)
      rw [add_zero] at this
      exact this (col_injOn hcol)
    rw [Finsupp.single_eq_of_ne hne, smul_zero, zero_add]
    have ih := ofSlots_apply_unique d s (f0 + 1) g (by simpa using hg) (fun g' hg' hne' => by
      have := hu (g' + 1) (by simpa using hg') (by omega)
      rw [show f0 + (g' + 1) = f0 + 1 + g' by omega, show f0 + (g + 1) = f0 + 1 + g by omega]
        at this
      exact this)
    rw [show f0 + (g + 1) = f0 + 1 + g by omega]
    exact ih

/-- **Non-membership soundness.**  A nonzero echelon residual in a free column: `toLo u ∉ S`.
`hchar` says that `K` has characteristic exactly `d.p`. -/
theorem not_mem_of_outsideOK (hd : Valid d) (hpK : ((d.p : ℕ) : K) = 0)
    (hchar : ∀ n : ℕ, (n : K) = 0 → d.p ∣ n) {u : List Entry} {f : ℕ}
    (h : outsideOK d u f = true) : toLo K u ∉ S K d := by
  obtain ⟨h1, hB, hf, huniq, hres⟩ := outsideOK_spec h
  have hok : ∀ e ∈ u, entryOK d e = true := fun e he => (h1 e he).1
  intro hmem
  have hacc := accum_eq hd u hok
  have hsl : ∀ s ∈ svAcc d u, s < d.B := fun s hs => lt_of_le_of_lt (svAcc_le' hd u hok s hs) hB
  have hun : unpack d.B d.nFree (accum d u) = svAcc d u := by
    rw [hacc, ← svAcc_length d u, unpack_packSlots (B_pos d) _ hsl]
  rw [hun] at hres
  have key := toLo_add_ofSlots_svAcc hd rfl hpK u hok
  have hr : ofSlots (K := K) d (svAcc d u) 0 ∈ S K d := by
    have : ofSlots (K := K) d (svAcc d u) 0 = pivPart d u - toLo K u := by
      rw [← key]; abel
    rw [this]
    exact Submodule.sub_mem _ (pivPart_mem hd u hok) hmem
  have hzero : ofSlots (K := K) d (svAcc d u) 0 = 0 := by
    have hrep := repr hd hr
    have hz : ∀ j : Fin d.nRows, ofSlots (K := K) d (svAcc d u) 0 (col (d.piv j)) = 0 := by
      intro j
      apply ofSlots_apply_eq_zero
      intro g hg hcol
      rw [svAcc_length] at hg
      rw [zero_add] at hcol
      exact freeCode_ne_piv hd hg j.isLt (col_injOn hcol)
    simp only [hz, zero_smul, Finset.sum_const_zero] at hrep
    exact hrep
  have hval := ofSlots_apply_unique (K := K) d (svAcc d u) 0 f (by rw [svAcc_length]; exact hf)
    (fun g' hg' hne => by
      rw [svAcc_length] at hg'
      rw [zero_add, zero_add]
      exact huniq g' hg' hne)
  rw [hzero, Finsupp.zero_apply] at hval
  exact hres (Nat.mod_eq_zero_of_dvd (hchar _ hval.symm))

end ICFallDegrees.N5CertWitness
