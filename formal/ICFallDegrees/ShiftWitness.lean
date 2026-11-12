import ICFallDegrees.ChainModelMap

noncomputable section

/-!
# N7: the lower bound without the unit-ideal hypothesis

`HasLFDBound F d` asks that *every* element `f` of the ideal lies in the closure of degree
`max d (deg f)`.  Hence `LFD > D` follows from a single ideal element of degree `≤ D` that is
not in the degree-`D` closure; the unit `1` is only one possible choice of such an element,
and it is not available when the ideal is proper.

This file removes the unit-ideal hypothesis from `ChainModelMap.lfd_gt_of_base`:

* `normalForm_mem_map_W`: under the base hypotheses of N5 (the free extension of `S` closed at
  depth `k0+1`, generator in `S`), the chain normal form of every element of the degree-`D`
  closure at depth `k ≥ k0+1` lies in the image `Φ (W k0 D S k)` of the free extension;
* `not_inClosure_of_not_mem_W`: so `Φ w` is not in the closure for every boxed `w ∉ W`;
* `lfd_gt_of_witness`, `lfd_gt_of_base_witness`: an ideal element `Φ w` of degree `≤ D` with
  `w ∉ W` (in particular `w = embed 0 v` with `v ∉ S` in the low box, weights `≤ D`) refutes
  every LFD bound `e ≤ D`;
* `LowIdeal k0 g`: the low elements obtained from the generator vector `g` by linear
  combinations and unit shifts that stay inside the low box (no carry out of the low window),
  a depth-independent description of ideal elements; `lowIdeal_mem_ideal` maps it into the
  ideal of the two-summand presentation at every depth `k > k0`;
* `lfd_gt_of_lowIdeal` and `growth_schema`: the schema "for every `D` and `k0`, a base `S` with
  `W_{k0+1}(S)` closed, the generator in `S`, and a low ideal element of weight `≤ D` outside `S`
  gives `LFD > D` at every depth `k ≥ k0+1`".  The schema is conditional on such a base for
  each `D`; no base is constructed here.

No hypothesis `1 ∉ S` and no unit-ideal hypothesis appears.
-/

namespace ICFallDegrees.ShiftWitness

open MvPolynomial FiniteQuotient ChainQuotient BinaryWeight DepthReduction ChainTranslation
open ChainModelMap ChainClosureBridge PrimeM2General

variable {m k k0 : ℕ} {K : Type*} [Field K]

/-! ## Three more properties of the model map -/

/-- Every value of the model map is a reduced (squarefree) polynomial. -/
theorem Φ_reduced (w : Al m K) : Φ k0 k w ∈ reducedSpace (chainInputs K m k) := by
  induction w using Finsupp.induction_linear with
  | zero => rw [map_zero]; exact (reducedSpace _).zero_mem
  | add f g hf hg => rw [map_add]; exact (reducedSpace _).add_mem hf hg
  | single a c => rw [Φ_single]; exact digits_reduced _ c

/-- The model map is injective on boxed elements. -/
theorem Φ_injOn (hk : k0 ≤ k) {w w' : Al m K} (hw : Boxed k0 k w) (hw' : Boxed k0 k w')
    (h : Φ k0 k w = Φ k0 k w') : w = w' := by
  ext a
  by_cases ha : InBox k0 k a
  · rw [← Φ_coeff hk hw ha, ← Φ_coeff hk hw' ha, h]
  · rw [Finsupp.not_mem_support_iff.mp (fun h' => ha (hw a h')),
      Finsupp.not_mem_support_iff.mp (fun h' => ha (hw' a h'))]

/-- The model map preserves the weight filtration on boxed elements (converse of `Φ_degree`). -/
theorem Φ_totalDegree_le (hk : k0 ≤ k) {w : Al m K} (hw : Boxed k0 k w) {D : ℕ}
    (hwt : ∀ a ∈ w.support, pwt a ≤ D) : (Φ k0 k w).totalDegree ≤ D := by
  classical
  unfold MvPolynomial.totalDegree
  apply Finset.sup_le
  intro s hs
  have hsub : (Φ k0 k w).support ⊆ w.support.image (fun a : Ex m => digits k (encode k0 a)) :=
    Finsupp.mapDomain_support
  obtain ⟨a, ha, rfl⟩ := Finset.mem_image.mp (hsub hs)
  have h1 := pwt_eq_degree hk (hw a ha)
  have h2 : (digits k (encode k0 a)).degree = (digits k (encode k0 a)).sum (fun _ e => e) := rfl
  have h3 := hwt a ha
  omega

/-! ## The containing space of the closure, without the constant -/

section Containing

variable (k0) (D : ℕ) (S : Submodule K (Lo 2 K))
  (base : Closed k0 (k0 + 1) D (W k0 D S (k0 + 1)))
  (k : ℕ) (hk : k0 + 1 ≤ k) (hk2 : 2 ≤ k) (two_ne : (2 : K) ≠ 0) (a b t : K)
  (gen : GeneratorInBase' k0 D S k hk2 a b t)

include base hk two_ne gen

/-- Under the N5 base hypotheses, the chain normal form of every element of the degree-`D`
closure of the two-summand presentation at depth `k ≥ k0+1` lies in `Φ (W k0 D S k)`. -/
theorem normalForm_mem_map_W {p : FiniteQuotient.P (2 * k) K}
    (h : InClosure (Set.range (inputs k hk2 a b t)) D p) :
    normalForm (chainInputs K 2 k) (chainInputs_nonzero 2 k) p ∈
      (W k0 D S k).map (Φ (K := K) k0 k) := by
  have hq := closure_normalForm_mem (nonzero := chainInputs_nonzero 2 k)
    (reduced_chain_ideal_trivial 2 k two_ne) h
  have hΦ := modelMap (K := K) (k0 := k0) (k := k) two_ne (by omega) D
  refine quot_closure_minimal _ D _ ?_ ?_ hq
  · rintro f ⟨i, rfl⟩ hd
    revert hd
    refine Fin.lastCases ?_ (fun j => ?_) i <;> intro hd
    · simp only [inputs, Fin.lastCases_last] at hd ⊢
      obtain ⟨v, hv, hvbox, hvwt, hΦv⟩ := gen hd
      exact ⟨embed 0 v, embed_zero_mem_W hv hvbox hvwt,
        hΦv.trans (normalForm_normalized two_ne k hk2 a b t).symm⟩
    · rw [normalForm_chain_eq_zero 2 k (inputs_castSucc_mem_chainInputs k hk2 a b t j) two_ne]
      exact Submodule.zero_mem _
  · rintro v' w' ⟨w, hw, rfl⟩ hdeg
    have hclosed := W_closed base k hk w hw (hΦ.degree w hw.1 hdeg) (ChainQuotient.row v')
      (ChainQuotient.depth v') (ChainQuotient.depth v').isLt
    refine ⟨translate k0 k (ChainQuotient.row v') (ChainQuotient.depth v') w, hclosed, ?_⟩
    have := hΦ.step w hw.1 (ChainQuotient.row v') (ChainQuotient.depth v')
    rw [var_row_depth] at this
    exact this.symm

/-- A boxed element outside the free extension is not in the degree-`D` closure. -/
theorem not_inClosure_of_not_mem_W {w : Al 2 K} (hwbox : Boxed k0 k w) (hw : w ∉ W k0 D S k) :
    ¬ InClosure (Set.range (inputs k hk2 a b t)) D (Φ k0 k w) := by
  intro h
  have hmem := normalForm_mem_map_W k0 D S base k hk hk2 two_ne a b t gen h
  rw [normalForm_of_reduced _ _ (reduced_chain_ideal_trivial 2 k two_ne) (Φ_reduced w)] at hmem
  obtain ⟨w', hw', hΦw⟩ := hmem
  have : w' = w := Φ_injOn (by omega) hw'.1 hwbox hΦw
  exact hw (this ▸ hw')

/-- **Lower bound from one ideal element.**  An element `Φ w` of the ideal, of degree `≤ D`,
with `w` boxed and outside the free extension, refutes every LFD bound `e ≤ D`.  No unit-ideal
hypothesis. -/
theorem lfd_gt_of_witness {w : Al 2 K} (hwbox : Boxed k0 k w) (hw : w ∉ W k0 D S k)
    (ideal : Φ k0 k w ∈ Ideal.span (Set.range (inputs k hk2 a b t)))
    (degree : (Φ k0 k w).totalDegree ≤ D) :
    ∀ e ≤ D, ¬ HasLFDBound (Set.range (inputs k hk2 a b t)) e := by
  intro e he bound
  exact not_inClosure_of_not_mem_W k0 D S base k hk hk2 two_ne a b t gen hwbox hw
    ((bound _ ideal).mono (max_le he degree))

/-- The same with a low witness: `v` in the low box, of weight `≤ D`, outside the base `S`,
whose image at pattern `0` lies in the ideal. -/
theorem lfd_gt_of_base_witness {v : Lo 2 K} (hvbox : ∀ l ∈ v.support, LowBox k0 l)
    (hvwt : ∀ l ∈ v.support, lowWt l ≤ D) (hvS : v ∉ S)
    (ideal : Φ k0 k (embed 0 v) ∈ Ideal.span (Set.range (inputs k hk2 a b t))) :
    ∀ e ≤ D, ¬ HasLFDBound (Set.range (inputs k hk2 a b t)) e := by
  classical
  have hbox : Boxed k0 k (embed (K := K) 0 v) := embed_boxed (lowBox_zero _) hvbox
  refine lfd_gt_of_witness k0 D S base k hk hk2 two_ne a b t gen hbox ?_ ideal ?_
  · intro hW
    have := (hW.2 0).1
    rw [piece_embed] at this
    exact hvS this
  · apply Φ_totalDegree_le (by omega) hbox
    intro a' ha'
    obtain ⟨l, hl, rfl⟩ := Finset.mem_image.mp (Finsupp.mapDomain_support ha')
    simpa [pwt] using hvwt l hl

end Containing

/-! ## Low ideal elements: shifts without carry out of the low window -/

/-- Translation of low exponents by `e`, as a linear map. -/
def lowShift (e : Fin m → ℕ) : Lo m K →ₗ[K] Lo m K := Finsupp.lmapDomain K K (fun l => l + e)

theorem lowShift_single (e l : Fin m → ℕ) (c : K) :
    lowShift e (Finsupp.single l c) = Finsupp.single (l + e) c := by
  show Finsupp.mapDomain _ (Finsupp.single l c) = _
  rw [Finsupp.mapDomain_single]

theorem lowShift_support (e : Fin m → ℕ) (v : Lo m K) {l' : Fin m → ℕ}
    (h : l' ∈ (lowShift e v).support) : ∃ l ∈ v.support, l' = l + e := by
  classical
  obtain ⟨l, hl, rfl⟩ := Finset.mem_image.mp (Finsupp.mapDomain_support h)
  exact ⟨l, hl, rfl⟩

theorem lowShift_zero (v : Lo m K) : lowShift (0 : Fin m → ℕ) v = v := by
  show Finsupp.mapDomain _ v = v
  have : (fun l : Fin m → ℕ => l + 0) = id := by funext l; simp
  rw [this, Finsupp.mapDomain_id]

theorem lowShift_add (e e' : Fin m → ℕ) (v : Lo m K) :
    lowShift (e + e') v = lowShift e' (lowShift e v) := by
  show Finsupp.mapDomain _ v = Finsupp.mapDomain _ (Finsupp.mapDomain _ v)
  rw [← Finsupp.mapDomain_comp]
  congr 1
  funext l
  simp [Function.comp, add_assoc]

/-- The low elements generated from `g` by linear combinations and unit shifts that keep the
support inside the low box `[0, 2^k0)^m`.  Since no carry leaves the low window, the same
element is an ideal element at every depth `k > k0` (`lowIdeal_mem_ideal`). -/
inductive LowIdeal (k0 : ℕ) (g : Lo m K) : Lo m K → Prop
  | gen : LowIdeal k0 g g
  | zero : LowIdeal k0 g 0
  | add {v v'} : LowIdeal k0 g v → LowIdeal k0 g v' → LowIdeal k0 g (v + v')
  | smul (c : K) {v} : LowIdeal k0 g v → LowIdeal k0 g (c • v)
  | shift (i : Fin m) {v} : LowIdeal k0 g v → (∀ l ∈ v.support, l i + 1 < 2 ^ k0) →
      LowIdeal k0 g (lowShift (Pi.single i 1) v)

namespace LowIdeal

variable {g : Lo m K}

/-- Low ideal elements stay in the low box. -/
theorem lowBox (hg : ∀ l ∈ g.support, LowBox k0 l) {v : Lo m K} (h : LowIdeal k0 g v) :
    ∀ l ∈ v.support, LowBox k0 l := by
  classical
  induction h with
  | gen => exact hg
  | zero => intro l hl; simp at hl
  | add _ _ ih₁ ih₂ =>
    intro l hl
    rcases Finset.mem_union.mp (Finsupp.support_add hl) with h | h
    · exact ih₁ l h
    · exact ih₂ l h
  | smul c _ ih => intro l hl; exact ih l (Finsupp.support_smul hl)
  | shift i _ hs ih =>
    intro l' hl' i'
    obtain ⟨l, hl, rfl⟩ := lowShift_support _ _ hl'
    rcases eq_or_ne i' i with rfl | hne
    · simpa using hs l hl
    · simpa [Pi.single_eq_of_ne hne] using ih l hl i'

/-- Shift by `n` in coordinate `i`, when the support stays inside the low box. -/
theorem shift_nsmul (i : Fin m) : ∀ (n : ℕ) {v : Lo m K}, LowIdeal k0 g v →
    (∀ l ∈ v.support, l i + n < 2 ^ k0) → LowIdeal k0 g (lowShift (n • Pi.single i 1) v)
  | 0, v, h, _ => by rw [zero_smul, lowShift_zero]; exact h
  | n + 1, v, h, hs => by
    rw [succ_nsmul, lowShift_add]
    refine LowIdeal.shift i (shift_nsmul i n h (fun l hl => by have := hs l hl; omega)) ?_
    intro l' hl'
    obtain ⟨l, hl, rfl⟩ := lowShift_support _ _ hl'
    have := hs l hl
    simp only [Pi.add_apply, Pi.smul_apply, Pi.single_eq_same, smul_eq_mul, mul_one]
    omega

end LowIdeal

/-- A unit shift without carry is the chain variable `u_{i,0}` at every depth. -/
theorem translate_embed_zero (i : Fin m) {v : Lo m K}
    (hs : ∀ l ∈ v.support, l i + 1 < 2 ^ k0) :
    translate k0 k i 0 (embed 0 v) = embed 0 (lowShift (Pi.single i 1) v) := by
  show Finsupp.mapDomain _ (Finsupp.mapDomain _ v) = Finsupp.mapDomain _ (Finsupp.mapDomain _ v)
  rw [← Finsupp.mapDomain_comp, ← Finsupp.mapDomain_comp]
  apply Finsupp.mapDomain_congr
  intro l hl
  have h1 := hs l hl
  have hk0 : 0 < k0 := by
    rcases Nat.eq_zero_or_pos k0 with h | h
    · rw [h] at h1; omega
    · exact h
  simp only [Function.comp, step, if_pos hk0, stepLow, pow_zero, if_pos h1]
  congr 1
  funext i'
  rcases eq_or_ne i' i with rfl | hne
  · simp
  · simp [Function.update_of_ne hne, Pi.single_eq_of_ne hne]

/-- Low ideal elements map into every ideal that contains the chain equations and the image
of the generator vector. -/
theorem lowIdeal_mem_ideal (two_ne : (2 : K) ≠ 0) (hk : k0 < k)
    (I : Ideal (FiniteQuotient.P (m * k) K)) (hI : chainInputs K m k ⊆ I) {g : Lo m K}
    (hgbox : ∀ l ∈ g.support, LowBox k0 l) (hg : Φ k0 k (embed 0 g) ∈ I) {v : Lo m K}
    (h : LowIdeal k0 g v) : Φ k0 k (embed 0 v) ∈ I := by
  induction h with
  | gen => exact hg
  | zero => simp
  | add _ _ ih₁ ih₂ => rw [map_add, map_add]; exact I.add_mem ih₁ ih₂
  | smul c _ ih => rw [map_smul, map_smul, smul_eq_C_mul]; exact I.mul_mem_left _ ih
  | @shift i v hv hs ih =>
    have hkpos : 0 < k := by omega
    have hbox : Boxed k0 k (embed (K := K) 0 v) :=
      embed_boxed (lowBox_zero _) (LowIdeal.lowBox hgbox hv)
    have hstep := Φ_step two_ne hk hbox i ⟨0, hkpos⟩
    simp only [] at hstep
    rw [← translate_embed_zero (k := k) i hs, ← hstep]
    have hsub := sub_normalForm_mem (chainInputs K m k) (chainInputs_nonzero m k)
      (X (var m k i ⟨0, hkpos⟩) * Φ k0 k (embed 0 v))
    have hin : X (var m k i ⟨0, hkpos⟩) * Φ k0 k (embed 0 v) ∈ I := I.mul_mem_left _ ih
    have := I.sub_mem hin (Ideal.span_le.mpr hI hsub)
    simpa using this

/-! ## The lower bound from a low ideal element -/

/-- **N7.**  Base at depth `k0` (`W_{k0+1}(S)` closed, generator vector `g ∈ S`), and a low
ideal element `v` of weight `≤ D` outside `S`: no `e ≤ D` is an LFD bound of the two-summand
presentation at depth `k ≥ k0+1`.  There is no unit-ideal hypothesis and no hypothesis
`1 ∉ S`. -/
theorem lfd_gt_of_lowIdeal (k0 D : ℕ) (S : Submodule K (Lo 2 K))
    (base : Closed k0 (k0 + 1) D (W k0 D S (k0 + 1)))
    (k : ℕ) (hk : k0 + 1 ≤ k) (hk2 : 2 ≤ k) (two_ne : (2 : K) ≠ 0) (a b t : K)
    (g : Lo 2 K) (hgS : g ∈ S) (hgbox : ∀ l ∈ g.support, LowBox k0 l)
    (hgwt : ∀ l ∈ g.support, lowWt l ≤ D)
    (hgΦ : Φ k0 k (embed 0 g) = normalized k hk2 a b t)
    (v : Lo 2 K) (hv : LowIdeal k0 g v) (hvwt : ∀ l ∈ v.support, lowWt l ≤ D) (hvS : v ∉ S) :
    ∀ e ≤ D, ¬ HasLFDBound (Set.range (inputs k hk2 a b t)) e := by
  have gen : GeneratorInBase' k0 D S k hk2 a b t := fun _ => ⟨g, hgS, hgbox, hgwt, hgΦ⟩
  refine lfd_gt_of_base_witness k0 D S base k hk hk2 two_ne a b t gen
    (LowIdeal.lowBox hgbox hv) hvwt hvS ?_
  refine lowIdeal_mem_ideal two_ne (by omega) _ ?_ hgbox ?_ hv
  · exact (chainInputs_two_subset k hk2 a b t).trans Ideal.subset_span
  · rw [hgΦ]
    apply Ideal.subset_span
    exact ⟨Fin.last _, by simp [inputs]⟩

end ICFallDegrees.ShiftWitness
