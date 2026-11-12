import ICFallDegrees.Closure

/-!
# Degree-aware algebraic certificates

The certificate contains input indices, scalars, and polynomial multipliers;
it contains no proof fields and cannot assert an external rank result.
The semantic checker is noncomputable because mathlib's `MvPolynomial`
representation is noncomputable. Its soundness is kernel-checked. A practical
computable sparse front end is a separate obligation, not claimed here.
-/

namespace ICFallDegrees

open MvPolynomial

variable {σ K : Type*} [CommRing K]

inductive Certificate (σ K : Type*) [CommRing K]
  | zero
  | input (index : ℕ)
  | add (left right : Certificate σ K)
  | scale (scalar : K) (parent : Certificate σ K)
  | mul (multiplier : MvPolynomial σ K) (parent : Certificate σ K)

noncomputable def checkCertificate (inputs : List (MvPolynomial σ K)) (D : ℕ) :
    Certificate σ K → Option (MvPolynomial σ K)
  | .zero => some 0
  | .input index => do
      let f ← inputs[index]?
      if f.totalDegree ≤ D then some f else none
  | .add left right => do
      let f ← checkCertificate inputs D left
      let g ← checkCertificate inputs D right
      some (f+g)
  | .scale scalar parent => do
      let f ← checkCertificate inputs D parent
      some (scalar • f)
  | .mul multiplier parent => do
      let f ← checkCertificate inputs D parent
      if multiplier.totalDegree+f.totalDegree ≤ D then some (multiplier*f) else none

theorem checkCertificate_sound (inputs : List (MvPolynomial σ K)) (D : ℕ)
    (certificate : Certificate σ K) {p}
    (checked : checkCertificate inputs D certificate=some p) :
    InClosure {f | f ∈ inputs} D p := by
  induction certificate generalizing p with
  | zero =>
      simp only [checkCertificate,Option.some.injEq] at checked
      subst p
      exact .zero
  | input index =>
      change (inputs[index]?.bind (fun f => if f.totalDegree ≤ D then some f else none))=some p at checked
      obtain ⟨f,hf,result⟩ := Option.bind_eq_some.mp checked
      by_cases hd : f.totalDegree ≤ D
      · simp only [if_pos hd,Option.some.injEq] at result
        subst p
        exact .input (List.mem_of_getElem? hf) hd
      · simp [hd] at result
  | add left right ihleft ihright =>
      change ((checkCertificate inputs D left).bind (fun f =>
        (checkCertificate inputs D right).bind (fun g => some (f+g))))=some p at checked
      obtain ⟨f,hf,rest⟩ := Option.bind_eq_some.mp checked
      obtain ⟨g,hg,result⟩ := Option.bind_eq_some.mp rest
      cases Option.some.inj result
      exact .add (ihleft hf) (ihright hg)
  | scale scalar parent ih =>
      change ((checkCertificate inputs D parent).bind (fun f => some (scalar • f)))=some p at checked
      obtain ⟨f,hf,result⟩ := Option.bind_eq_some.mp checked
      cases Option.some.inj result
      exact .smul scalar (ih hf)
  | mul multiplier parent ih =>
      change ((checkCertificate inputs D parent).bind (fun f =>
        if multiplier.totalDegree+f.totalDegree ≤ D then some (multiplier*f) else none))=some p at checked
      obtain ⟨f,hf,result⟩ := Option.bind_eq_some.mp checked
      by_cases hd : multiplier.totalDegree+f.totalDegree ≤ D
      · simp only [if_pos hd,Option.some.injEq] at result
        subst p
        exact .mul multiplier (ih hf) hd
      · simp [hd] at result

theorem checkCertificate_ideal_sound (inputs : List (MvPolynomial σ K)) (D : ℕ)
    (certificate : Certificate σ K) {p}
    (checked : checkCertificate inputs D certificate=some p) :
    p ∈ Ideal.span {f | f ∈ inputs} := (checkCertificate_sound inputs D certificate checked).sound

end ICFallDegrees
