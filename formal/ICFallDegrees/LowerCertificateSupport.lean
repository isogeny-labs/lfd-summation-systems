import ICFallDegrees.SparseLowerCertificate

noncomputable section

namespace ICFallDegrees.LowerCertificateSupport

open MvPolynomial SparsePolynomial SparseLowerCertificate

variable {r n : ℕ} {K : Type*} [CommRing K] [Nontrivial K] [DecidableEq K]

theorem checked_rref (inputs : List (Sparse r K)) (basis : Fin n → Sparse r K)
    (pivot : Fin n → Exponent r) (D : ℕ) (checked : check inputs basis pivot D=true) :
    ∀ i j, coeff (toExponent (pivot i)) (meaning (basis j))=if i=j then 1 else 0 := by
  have rows : ∀ i, rowCheck basis pivot D i=true := by
    have pair := Bool.and_eq_true_iff.mp checked
    simpa [List.all_eq_true,List.forall_mem_ofFn_iff] using pair.1
  intro i j
  rw [meaning_coefficient]
  exact (rowCheck_spec basis pivot D i (rows i)).2.2.2.1 j

theorem checked_closure_containment (inputs : List (Sparse r K)) (basis : Fin n → Sparse r K)
    (pivot : Fin n → Exponent r) (D : ℕ) (checked : check inputs basis pivot D=true) :
    closureSpace {p | p∈inputs.map meaning} D≤FilteredSpan.rowSpace (fun i => meaning (basis i)) := by
  have checks : (∀ i, rowCheck basis pivot D i=true) ∧
      (∀ p∈inputs, contains basis pivot p=true) := by
    simpa [check,List.all_eq_true,List.forall_mem_ofFn_iff] using checked
  have rows := fun i => rowCheck_spec basis pivot D i (checks.1 i)
  apply FilteredSpan.closure_le_rowSpace {p | p∈inputs.map meaning}
    (fun i => meaning (basis i)) (fun i => toExponent (pivot i)) D
    (checked_rref inputs basis pivot D checked)
  · intro i degree v
    rw [weight_toExponent] at degree
    have member := contains_sound basis pivot (mul (variableData v) (basis i))
      ((rows i).2.2.2.2 degree v)
    simpa only [meaning_mul,meaning_variable] using member
  · intro p hp
    obtain ⟨data,member,rfl⟩ := List.mem_map.mp hp
    exact contains_sound basis pivot data (checks.2 data member)

theorem rowSpace_coefficient_zero (basis : Fin n → MvPolynomial (Fin r) K)
    (e : Fin r →₀ ℕ) (zero : ∀ i, coeff e (basis i)=0)
    {p} (member : p∈FilteredSpan.rowSpace basis) : coeff e p=0 := by
  classical
  obtain ⟨c,hc⟩ := member
  rw [←hc]
  simp [FilteredSpan.rowMap,Fintype.linearCombination_apply,coeff_sum,
    smul_eq_C_mul,coeff_C_mul,zero]

theorem rowSpace_no_low (basis : Fin n → MvPolynomial (Fin r) K)
    (pivot : Fin n → Fin r →₀ ℕ)
    (rref : ∀ i j, coeff (pivot i) (basis j)=if i=j then 1 else 0)
    (d : ℕ) (minimum : ∀ i, d≤(pivot i).degree)
    {p} (member : p∈FilteredSpan.rowSpace basis) (degree : p.totalDegree<d) : p=0 := by
  classical
  obtain ⟨c,hc⟩ := member
  have zero : ∀ i, c i=0 := by
    intro i
    rw [←FilteredSpan.coefficient_rowMap basis pivot rref c i,hc]
    apply coeff_eq_zero_of_totalDegree_lt
    exact lt_of_lt_of_le degree (minimum i)
  rw [←hc]
  simp [FilteredSpan.rowMap,Fintype.linearCombination_apply,zero]

end ICFallDegrees.LowerCertificateSupport
