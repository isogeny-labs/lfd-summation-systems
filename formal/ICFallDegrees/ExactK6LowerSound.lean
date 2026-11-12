import ICFallDegrees.ExactK6Cert
import ICFallDegrees.N5CertSound

/-!
# Soundness of the depth-6 containing-space checker (`1 ∉ V_4`)

The data is an `N5Cert.Data` record (reduced echelon basis, packed rows) whose codes are the
complement codes `m = 4095 - (e₀ + 64 e₁)` of the depth-6 squarefree chain monomials
(`ExactK6Cert`).  With `d.low = []` the depth-7 specific checks of `N5CertSound.Valid` are
vacuous, so `Valid d` holds and the echelon membership test `N5CertSound.toLo_mem_of_member`,
the representation `N5CertSound.repr` and the row space `N5CertSound.S` are used unchanged.

New here:

* `Ψ`, the model map `δ_{col m} ↦` the squarefree chain monomial of `m`, injective on the box
  (`Ψ_coeff`), so that it reflects `1` (`Ψ_ne_one`) and the degree filtration
  (`piv_weight_of_degree`);
* `shift_single`: reduced multiplication by `u_{i,j}` at depth 6 is `shiftCode i j` on codes
  (`ChainTranslation.normalForm_X_mul_monomial`, exact wrap `2^6 ≡ 0`);
* `not_one`: the witness theorem `ChainClosureBridge.prime_m2_not_one` at depth 6 with the
  witness `Ψ (S d)`.

Only `propext`, `Classical.choice`, `Quot.sound` are used.
-/

noncomputable section

namespace ICFallDegrees.ExactK6Lower

open ICFallDegrees.N5Cert ICFallDegrees.N5CertSound
open MvPolynomial FiniteQuotient ChainQuotient ChainTranslation ChainClosureBridge BinaryWeight
open DepthReduction

/-- Weights of complement codes: `wt(e₀) + wt(e₁) + popcount m = 12`. -/
theorem weight_fact : allBelow (fun m => popcount 6 ((4095 - m) % 64) +
    popcount 6 ((4095 - m) / 64) + popcount 16 m == 12) 4096 = true := by decide +kernel

variable {K : Type*} [Field K]

/-- The block exponents of the complement code carried by a low index. -/
def dec (l : Fin 2 → ℕ) : Fin 2 → ℕ :=
  ![(4095 - (l 0 + 128 * l 1)) % 64, (4095 - (l 0 + 128 * l 1)) / 64]

theorem dec_col (m : ℕ) : dec (col m) = ![(4095 - m) % 64, (4095 - m) / 64] := by
  unfold dec
  rw [col_zero, col_one, Nat.mod_add_div]

theorem dec_col_lt (m : ℕ) (i : Fin 2) : dec (col m) i < 2 ^ 6 := by
  rw [dec_col]
  fin_cases i
  · show (4095 - m) % 64 < 2 ^ 6
    omega
  · show (4095 - m) / 64 < 2 ^ 6
    omega

/-- The low indices of the depth-6 box. -/
def Bx : Set (Fin 2 → ℕ) := {l | ∃ m, m < 4096 ∧ l = col m}

theorem exponent_injOn : Set.InjOn (fun l : Fin 2 → ℕ => digits 6 (dec l)) Bx := by
  rintro _ ⟨m, hm, rfl⟩ _ ⟨m', hm', rfl⟩ h
  have h' := ChainModelMap.digits_injOn (dec_col_lt m) (dec_col_lt m') h
  rw [dec_col, dec_col] at h'
  have h0 := congrFun h' 0
  have h1 := congrFun h' 1
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons] at h0 h1
  have : m = m' := by omega
  rw [this]

variable (K) in
/-- The model map: `δ_l ↦` the squarefree chain monomial with block exponents `dec l`. -/
def Ψ : Lo 2 K →ₗ[K] P (2 * 6) K := Finsupp.lmapDomain K K (fun l => digits 6 (dec l))

theorem Ψ_single (l : Fin 2 → ℕ) (c : K) :
    Ψ K (Finsupp.single l c) = monomial (digits 6 (dec l)) c := by
  show Finsupp.mapDomain _ (Finsupp.single l c) = _
  rw [Finsupp.mapDomain_single, single_eq_monomial]

theorem Ψ_eCode (m : ℕ) : Ψ K (eCode K m) = monomial (digits 6 (dec (col m))) 1 :=
  Ψ_single _ _

theorem Ψ_toLo : ∀ L : List Entry,
    Ψ K (toLo K L) = (L.map fun e => (e.2 : K) • monomial (digits 6 (dec (col e.1))) (1 : K)).sum
  | [] => by simp [toLo_nil]
  | e :: L => by
    rw [toLo_cons, map_add, map_smul, Ψ_eCode, Ψ_toLo L, List.map_cons, List.sum_cons]

theorem Ψ_coeff {s : Lo 2 K} (hs : (s.support : Set (Fin 2 → ℕ)) ⊆ Bx) {l : Fin 2 → ℕ}
    (hl : l ∈ Bx) : coeff (digits 6 (dec l)) (Ψ K s) = s l := by
  change Finsupp.mapDomain (fun l : Fin 2 → ℕ => digits 6 (dec l)) s (digits 6 (dec l)) = s l
  exact Finsupp.mapDomain_apply' Bx s hs exponent_injOn hl

/-! ## The row space lives in the box -/

section Box

variable {d : Data}

theorem box_piv (hbox : ExactK6Cert.boxOK d = true) {i : ℕ} (hi : i < d.nRows) :
    d.piv i < 4096 := by
  simp only [ExactK6Cert.boxOK, Bool.and_eq_true] at hbox
  simpa using allBelow_sound _ _ hbox.1 i hi

theorem box_free (hbox : ExactK6Cert.boxOK d = true) {f : ℕ} (hf : f < d.nFree) :
    d.freeCode f < 4096 := by
  simp only [ExactK6Cert.boxOK, Bool.and_eq_true] at hbox
  simpa using allBelow_sound _ _ hbox.2 f hf

theorem eCode_supported {m : ℕ} (hm : m < 4096) : eCode K m ∈ Finsupp.supported K K Bx :=
  Finsupp.single_mem_supported K _ ⟨m, hm, rfl⟩

theorem ofSlots_supported (hbox : ExactK6Cert.boxOK d = true) : ∀ (L : List ℕ) (f : ℕ),
    f + L.length ≤ d.nFree → ofSlots (K := K) d L f ∈ Finsupp.supported K K Bx
  | [], _, _ => Submodule.zero_mem _
  | s :: L, f, h => by
    simp only [List.length_cons] at h
    exact Submodule.add_mem _ (Submodule.smul_mem _ _ (eCode_supported (box_free hbox (by omega))))
      (ofSlots_supported hbox L (f + 1) (by omega))

theorem S_le_supported (hbox : ExactK6Cert.boxOK d = true) :
    S K d ≤ Finsupp.supported K K Bx := by
  apply Submodule.span_le.mpr
  rintro _ ⟨i, rfl⟩
  unfold rowVec rowFree
  exact Submodule.add_mem _ (eCode_supported (box_piv hbox i.isLt))
    (ofSlots_supported hbox _ 0 (by rw [unpack_length]; omega))

theorem support_subset (hbox : ExactK6Cert.boxOK d = true) {s : Lo 2 K} (hs : s ∈ S K d) :
    (s.support : Set (Fin 2 → ℕ)) ⊆ Bx :=
  (Finsupp.mem_supported K s).mp (S_le_supported hbox hs)

/-! ## `1` is not in the image -/

theorem digits_dec_eq_zero {m : ℕ} (hm : m < 4096) (h : digits 6 (dec (col m)) = 0) :
    m = 4095 := by
  have h' := (ChainModelMap.digits_zero_iff (dec (col m)) (dec_col_lt m)).mp h
  rw [dec_col] at h'
  have h0 := congrFun h' 0
  have h1 := congrFun h' 1
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons, Pi.zero_apply] at h0 h1
  omega

theorem Ψ_ne_one (hd : Valid d) (hbox : ExactK6Cert.boxOK d = true)
    (hconst : ExactK6Cert.constOK d = true) {s : Lo 2 K} (hs : s ∈ S K d) :
    Ψ K s ≠ (1 : P (2 * 6) K) := by
  classical
  intro h
  have hz : ∀ j : Fin d.nRows, s (col (d.piv j)) = 0 := by
    intro j
    rw [← Ψ_coeff (support_subset hbox hs) ⟨d.piv j, box_piv hbox j.isLt, rfl⟩, h, coeff_one,
      if_neg]
    intro h0
    have := digits_dec_eq_zero (box_piv hbox j.isLt) h0.symm
    have h1 := tbl_piv hd j.isLt
    simp only [ExactK6Cert.constOK, decide_eq_true_eq] at hconst
    rw [this] at h1
    have := hd.nRows_lt
    omega
  have hr := N5CertSound.repr hd hs
  simp only [hz, zero_smul, Finset.sum_const_zero] at hr
  rw [hr, map_zero] at h
  exact zero_ne_one h

/-! ## Low rows -/

/-- The span of the rows listed in `low`. -/
def S3 (K : Type*) [Field K] (d : Data) (low : List (ℕ × List Entry)) : Submodule K (Lo 2 K) :=
  Submodule.span K (Set.range fun i : Fin low.length => rowVec d (low.get i).1)

theorem mem_S3 (hd : Valid d) {low : List (ℕ × List Entry)}
    (hidx : ExactK6Cert.lowIdxOK d low = true) {v : Lo 2 K} (hv : v ∈ S K d)
    (hw : ∀ j : Fin d.nRows, v (col (d.piv j)) ≠ 0 → 9 ≤ popcount 16 (d.piv j)) :
    v ∈ S3 K d low := by
  rw [N5CertSound.repr hd hv]
  apply Submodule.sum_mem
  intro j _
  by_cases h0 : v (col (d.piv j)) = 0
  · rw [h0, zero_smul]
    exact Submodule.zero_mem _
  · apply Submodule.smul_mem
    have hmem : (j : ℕ) ∈ low.map Prod.fst := by
      simp only [ExactK6Cert.lowIdxOK, beq_iff_eq] at hidx
      rw [← hidx, List.mem_filter, List.mem_range]
      exact ⟨j.isLt, by simpa using hw j h0⟩
    obtain ⟨r, hr, hrj⟩ := List.mem_map.mp hmem
    obtain ⟨n, hn⟩ := List.mem_iff_get.mp hr
    apply Submodule.subset_span
    refine ⟨n, ?_⟩
    simp only [hn, hrj]

theorem lowRowOK_iff (d : Data) (r : ℕ × List Entry) :
    ExactK6Cert.lowRowOK d r = true ↔ r.1 < d.nRows ∧ r.2.length ≤ 256 ∧
      (∀ e ∈ r.2, freeEntryOK d e = true) ∧ packFree d r.2 = d.row r.1 ∧
      ExactK6Cert.shiftsOK d r = true := by
  unfold ExactK6Cert.lowRowOK
  rw [Bool.and_eq_true, Bool.and_eq_true, Bool.and_eq_true, Bool.and_eq_true, decide_eq_true_eq,
    decide_eq_true_eq, List.all_eq_true, beq_iff_eq]
  tauto

/-- The sparse form of a low row is its row vector. -/
theorem toLo_lowFull (hd : Valid d) {r : ℕ × List Entry} (hr : ExactK6Cert.lowRowOK d r = true) :
    toLo K (lowFull d r) = rowVec d r.1 := by
  obtain ⟨hlt, hlen, hok, hpack, _⟩ := (lowRowOK_iff d r).mp hr
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
  rw [← hpack, packFree_eq d _ hok, ← svFree_length d r.2, unpack_packSlots (N5CertSound.B_pos d) _ hsl]

/-- Every code of the sparse form of a low row is in the box. -/
theorem lowFull_codes (hbox : ExactK6Cert.boxOK d = true) {r : ℕ × List Entry}
    (hr : ExactK6Cert.lowRowOK d r = true) : ∀ e ∈ lowFull d r, e.1 < 4096 := by
  obtain ⟨hlt, _, hok, _, _⟩ := (lowRowOK_iff d r).mp hr
  intro e he
  unfold lowFull at he
  rcases List.mem_cons.mp he with rfl | he
  · exact box_piv hbox hlt
  · obtain ⟨_, _, _, v, _, _, hv2, hv3⟩ := (freeEntryOK_iff d e).mp (hok e he)
    rw [← hv3]
    exact box_free hbox hv2

end Box

/-! ## The translations -/

section Shift

variable (two_ne : (2 : K) ≠ 0)
include two_ne

theorem shift_single {m : ℕ} (hm : m < 4096) (i : Fin 2) (j : Fin 6) :
    normalForm (chainInputs K 2 6) (chainInputs_nonzero 2 6)
        (X (var 2 6 i j) * monomial (digits 6 (dec (col m))) (1 : K)) =
      monomial (digits 6 (dec (col (ExactK6Cert.shiftCode i j m)))) 1 := by
  rw [normalForm_X_mul_monomial two_ne i j _ (dec_col_lt m)]
  congr 2
  rw [dec_col, dec_col]
  have e : (2 : ℕ) ^ 6 = 64 := by norm_num
  rw [e]
  unfold ExactK6Cert.shiftCode
  generalize (2 : ℕ) ^ (j : ℕ) = q
  congr 1
  funext i'
  fin_cases i <;> fin_cases i' <;>
    simp only [Fin.zero_eta, Fin.mk_one, Fin.isValue, Function.update_self, Fin.val_zero,
      Fin.val_one, if_true, one_ne_zero, if_false, Matrix.cons_val_zero, Matrix.cons_val_one,
      Matrix.head_cons, ne_eq, zero_ne_one, not_false_eq_true, Function.update_of_ne,
      Fin.one_eq_zero_iff, OfNat.ofNat_ne_one, Nat.reduceAdd] <;> omega

theorem shift_toLo (i : Fin 2) (j : Fin 6) : ∀ L : List Entry, (∀ e ∈ L, e.1 < 4096) →
    normalForm (chainInputs K 2 6) (chainInputs_nonzero 2 6) (X (var 2 6 i j) * Ψ K (toLo K L)) =
      Ψ K (toLo K (ExactK6Cert.shiftList i j L))
  | [], _ => by
    have zero := reduced_chain_ideal_trivial (K := K) 2 6 two_ne
    simp [ExactK6Cert.shiftList, toLo_nil, normalForm_zero _ _ zero]
  | e :: L, h => by
    have zero := reduced_chain_ideal_trivial (K := K) 2 6 two_ne
    have ih := shift_toLo i j L (fun x hx => h x (List.mem_cons_of_mem _ hx))
    unfold ExactK6Cert.shiftList at ih ⊢
    rw [List.map_cons, toLo_cons, toLo_cons, map_add, map_add, map_smul, map_smul, Ψ_eCode, Ψ_eCode,
      mul_add, normalForm_add _ _ zero, ih, mul_smul_comm, normalForm_smul _ _ zero,
      shift_single two_ne (h e (List.mem_cons_self))]

end Shift

/-! ## The witness theorem at depth 6 -/

open PrimeM2General in
/-- **`1 ∉ V_4` at depth 6.**  A valid echelon certificate whose low rows are closed under the
twelve translations, whose constant column is free and which contains the generator refutes `1`
in the degree-4 closure of the two-summand chain presentation at depth 6. -/
theorem not_one (d : Data) (hd : Valid d) (hpK : ((d.p : ℕ) : K) = 0) (two_ne : (2 : K) ≠ 0)
    (low : List (ℕ × List Entry)) (hbox : ExactK6Cert.boxOK d = true)
    (hconst : ExactK6Cert.constOK d = true) (hidx : ExactK6Cert.lowIdxOK d low = true)
    (hlow : ∀ r ∈ low, ExactK6Cert.lowRowOK d r = true) (a b t : K)
    (hgen : Ψ K (toLo K d.gen) = normalized 6 (by norm_num) a b t) :
    ¬ InClosure (Set.range (inputs 6 (by norm_num) a b t)) 4 (1 : FiniteQuotient.P (2 * 6) K) := by
  classical
  have zero := reduced_chain_ideal_trivial (K := K) 2 6 two_ne
  apply prime_m2_not_one 6 (by norm_num) a b t two_ne 4 ((S K d).map (Ψ K))
  · intro _
    rw [ChainModelMap.normalForm_normalized two_ne]
    exact ⟨toLo K d.gen, gen_mem hd hpK, hgen⟩
  · rintro v w ⟨s, hs, rfl⟩ hdeg
    -- the pivot coordinates of `s` have weight `≤ 3`
    have hwt : ∀ j : Fin d.nRows, s (col (d.piv j)) ≠ 0 → 9 ≤ popcount 16 (d.piv j) := by
      intro j hne
      have hm := box_piv hbox j.isLt
      have hmem : digits 6 (dec (col (d.piv j))) ∈ (Ψ K s).support := by
        rw [mem_support_iff, Ψ_coeff (support_subset hbox hs) ⟨d.piv j, hm, rfl⟩]
        exact hne
      have hle := le_totalDegree hmem
      have hrfl : ((digits 6 (dec (col (d.piv j)))).sum fun _ e => e) =
          (digits 6 (dec (col (d.piv j)))).degree := rfl
      rw [hrfl, ChainModelMap.digits_degree _ (dec_col_lt _), Fin.sum_univ_two, dec_col] at hle
      simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons] at hle
      rw [← popcount_eq_wt 6 _ (by omega), ← popcount_eq_wt 6 _ (by omega)] at hle
      have hf := allBelow_sound _ _ weight_fact (d.piv j) hm
      simp only [beq_iff_eq] at hf
      omega
    have hs3 := mem_S3 hd hidx hs hwt
    -- the linear map `s ↦ normalForm (X v * Ψ s)`
    let T : Lo 2 K →ₗ[K] FiniteQuotient.P (2 * 6) K :=
      (nfMap (chainInputs K 2 6) (chainInputs_nonzero 2 6) zero).comp
        ((LinearMap.mulLeft K (X v)).comp (Ψ K))
    have hT : ∀ x, T x = normalForm (chainInputs K 2 6) (chainInputs_nonzero 2 6) (X v * Ψ K x) :=
      fun _ => rfl
    rw [← hT]
    have hle : S3 K d low ≤ ((S K d).map (Ψ K)).comap T := by
      apply Submodule.span_le.mpr
      rintro _ ⟨n, rfl⟩
      have hr := hlow _ (List.get_mem low n)
      obtain ⟨_, _, _, _, hsh⟩ := (lowRowOK_iff d _).mp hr
      show T (rowVec d (low.get n).1) ∈ Submodule.map (Ψ K) (S K d)
      rw [hT, ← toLo_lowFull hd hr, ← var_row_depth v,
        shift_toLo two_ne (row v) (depth v) _ (lowFull_codes hbox hr)]
      refine ⟨_, toLo_mem_of_member hd hpK ?_, rfl⟩
      have hj := allBelow_sound _ _ hsh (depth v : ℕ) (depth v).isLt
      simp only [Bool.and_eq_true] at hj
      have hrow : (row v : ℕ) = 0 ∨ (row v : ℕ) = 1 := by have := (row v).isLt; omega
      rcases hrow with h0 | h1
      · rw [h0]; exact hj.1
      · rw [h1]; exact hj.2
    exact hle hs3
  · rintro ⟨s, hs, h⟩
    exact Ψ_ne_one hd hbox hconst hs h

end ICFallDegrees.ExactK6Lower

end
