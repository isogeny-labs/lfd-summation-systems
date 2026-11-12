import ICFallDegrees.SparsePolynomial

/-!
# Executable degree-aware certificate checking

All checking and arithmetic use sparse data. Theorems connect successful
checks to the mathematical full V-space, including multiplication of newly
fallen rows. Certificates contain no proof fields or external rank claims.
-/

namespace ICFallDegrees.SparseCertificate

open SparsePolynomial

inductive Proof (r : ℕ) (K : Type*)
  | zero
  | input (index : ℕ)
  | add (left right : Proof r K)
  | scale (scalar : K) (parent : Proof r K)
  | mul (multiplier : Sparse r K) (parent : Proof r K)

variable {r : ℕ} {K : Type*} [CommRing K] [DecidableEq K]

def check (inputs : List (Sparse r K)) (D : ℕ) : Proof r K → Option (Sparse r K)
  | .zero => some []
  | .input index => do
      let p ← inputs[index]?
      let normalized := normalize p
      if degree normalized≤D then some normalized else none
  | .add left right => do
      let p ← check inputs D left
      let q ← check inputs D right
      some (add p q)
  | .scale scalar parent => do
      let p ← check inputs D parent
      some (scale scalar p)
  | .mul multiplier parent => do
      let p ← check inputs D parent
      let normalized := normalize multiplier
      if degree normalized+degree p≤D then some (mul normalized p) else none

def oneData : Sparse r K := [(fun _ => 0,1)]

omit [DecidableEq K] in
theorem meaning_oneData : meaning (oneData : Sparse r K)=1 := by
  have zero : toExponent (fun _ : Fin r => 0)=0 := by ext i; rfl
  simp [oneData,meaning,zero,MvPolynomial.monomial_zero']

theorem check_sound (inputs : List (Sparse r K)) (D : ℕ)
    (certificate : Proof r K) {out : Sparse r K}
    (checked : check inputs D certificate=some out) :
    InClosure {p | p∈inputs.map meaning} D (meaning out) := by
  induction certificate generalizing out with
  | zero =>
      simp only [check,Option.some.injEq] at checked
      subst out
      exact InClosure.zero
  | input index =>
      change (inputs[index]?.bind (fun p =>
        if degree (normalize p)≤D then some (normalize p) else none))=some out at checked
      obtain ⟨p,hp,result⟩ := Option.bind_eq_some.mp checked
      by_cases hd : degree (normalize p)≤D
      · simp only [if_pos hd,Option.some.injEq] at result
        subst out
        have budget : (meaning p).totalDegree≤D := by
          rw [←meaning_normalize p]
          exact (meaning_degree_le _).trans hd
        have member : meaning p∈inputs.map meaning :=
          List.mem_map_of_mem (List.mem_of_getElem? hp)
        simpa only [meaning_normalize] using InClosure.input member budget
      · simp [hd] at result
  | add left right ihleft ihright =>
      change ((check inputs D left).bind (fun p =>
        (check inputs D right).bind (fun q => some (add p q))))=some out at checked
      obtain ⟨p,hp,rest⟩ := Option.bind_eq_some.mp checked
      obtain ⟨q,hq,result⟩ := Option.bind_eq_some.mp rest
      cases Option.some.inj result
      rw [meaning_add]
      exact InClosure.add (ihleft hp) (ihright hq)
  | scale scalar parent ih =>
      change ((check inputs D parent).bind (fun p => some (scale scalar p)))=some out at checked
      obtain ⟨p,hp,result⟩ := Option.bind_eq_some.mp checked
      cases Option.some.inj result
      rw [meaning_scale]
      exact InClosure.smul scalar (ih hp)
  | mul multiplier parent ih =>
      change ((check inputs D parent).bind (fun p =>
        if degree (normalize multiplier)+degree p≤D then
          some (mul (normalize multiplier) p) else none))=some out at checked
      obtain ⟨p,hp,result⟩ := Option.bind_eq_some.mp checked
      by_cases hd : degree (normalize multiplier)+degree p≤D
      · simp only [if_pos hd,Option.some.injEq] at result
        subst out
        rw [meaning_mul]
        exact InClosure.mul (meaning (normalize multiplier)) (ih hp)
          ((add_le_add (meaning_degree_le _) (meaning_degree_le _)).trans hd)
      · simp [hd] at result

theorem check_ideal_sound (inputs : List (Sparse r K)) (D : ℕ)
    (certificate : Proof r K) {out : Sparse r K}
    (checked : check inputs D certificate=some out) :
    meaning out∈Ideal.span {p | p∈inputs.map meaning} := (check_sound inputs D certificate checked).sound

theorem checked_one_lfd_bound (inputs : List (Sparse r K)) (D : ℕ)
    (certificate : Proof r K)
    (checked : check inputs D certificate=some oneData) :
    HasLFDBound {p | p∈inputs.map meaning} D := by
  apply hasLFDBound_of_one
  simpa only [meaning_oneData] using check_sound inputs D certificate checked

end ICFallDegrees.SparseCertificate
