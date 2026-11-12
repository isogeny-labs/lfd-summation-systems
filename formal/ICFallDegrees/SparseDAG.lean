import ICFallDegrees.SparseCertificate

/-! Shared, indexed sparse derivations. Every reference must point to an
already checked row; linear combinations need no external rank result. -/

namespace ICFallDegrees.SparseDAG

open SparsePolynomial

inductive Step (r : ℕ) (K : Type*)
  | input (index : ℕ)
  | linear (weights : List (ℕ × K))
  | mul (multiplier : Sparse r K) (parent : ℕ)

variable {r : ℕ} {K : Type*} [CommRing K] [DecidableEq K]

def evalLinear (rows : List (Sparse r K)) : List (ℕ × K) → Option (Sparse r K)
  | [] => some []
  | (index,c)::rest => do
      let p ← rows[index]?
      let q ← evalLinear rows rest
      some (add (scale c p) q)

def evalStep (inputs : List (Sparse r K)) (D : ℕ) (rows : List (Sparse r K)) :
    Step r K → Option (Sparse r K)
  | .input index => do
      let p ← inputs[index]?
      if degree (normalize p)≤D then some (normalize p) else none
  | .linear weights => evalLinear rows weights
  | .mul multiplier parent => do
      let p ← rows[parent]?
      if degree (normalize multiplier)+degree p≤D then some (mul (normalize multiplier) p) else none

def runFrom (inputs : List (Sparse r K)) (D : ℕ) (rows : List (Sparse r K)) :
    List (Step r K) → Option (List (Sparse r K))
  | [] => some rows
  | instruction::rest => do
      let p ← evalStep inputs D rows instruction
      runFrom inputs D (rows++[p]) rest

def run (inputs : List (Sparse r K)) (D : ℕ) (program : List (Step r K)) :
    Option (List (Sparse r K)) := runFrom inputs D [] program

def checkUnit (inputs : List (Sparse r K)) (D : ℕ) (program : List (Step r K)) : Bool :=
  match run inputs D program with
  | none => false
  | some rows => match rows.getLast? with
      | none => false
      | some p => (add p (scale (-1) SparseCertificate.oneData)).isEmpty

theorem evalLinear_sound (F : Set (MvPolynomial (Fin r) K)) (D : ℕ)
    (rows : List (Sparse r K)) (known : ∀ p∈rows, InClosure F D (meaning p))
    (weights : List (ℕ × K)) {out : Sparse r K}
    (checked : evalLinear rows weights=some out) : InClosure F D (meaning out) := by
  induction weights generalizing out with
  | nil => simp only [evalLinear,Option.some.injEq] at checked; subst out; exact .zero
  | cons item rest ih =>
      obtain ⟨index,c⟩ := item
      change (rows[index]?.bind (fun p => (evalLinear rows rest).bind
        (fun q => some (add (scale c p) q))))=some out at checked
      obtain ⟨p,hp,tail⟩ := Option.bind_eq_some.mp checked
      obtain ⟨q,hq,result⟩ := Option.bind_eq_some.mp tail
      cases Option.some.inj result
      rw [meaning_add,meaning_scale]
      exact .add (.smul c (known p (List.mem_of_getElem? hp))) (ih hq)

theorem evalStep_sound (inputs : List (Sparse r K)) (D : ℕ)
    (rows : List (Sparse r K))
    (known : ∀ p∈rows, InClosure {q | q∈inputs.map meaning} D (meaning p))
    (instruction : Step r K) {out : Sparse r K}
    (checked : evalStep inputs D rows instruction=some out) :
    InClosure {q | q∈inputs.map meaning} D (meaning out) := by
  cases instruction with
  | input index =>
      change (inputs[index]?.bind (fun p =>
        if degree (normalize p)≤D then some (normalize p) else none))=some out at checked
      obtain ⟨p,hp,result⟩ := Option.bind_eq_some.mp checked
      by_cases hd : degree (normalize p)≤D
      · simp only [if_pos hd,Option.some.injEq] at result
        subst out
        have member : meaning p∈inputs.map meaning := List.mem_map_of_mem (List.mem_of_getElem? hp)
        have budget : (meaning p).totalDegree≤D := by
          rw [←meaning_normalize p]
          exact (meaning_degree_le _).trans hd
        simpa only [meaning_normalize] using InClosure.input member budget
      · simp [hd] at result
  | linear weights => exact evalLinear_sound _ D rows known weights checked
  | mul multiplier parent =>
      change (rows[parent]?.bind (fun p =>
        if degree (normalize multiplier)+degree p≤D then
          some (mul (normalize multiplier) p) else none))=some out at checked
      obtain ⟨p,hp,result⟩ := Option.bind_eq_some.mp checked
      by_cases hd : degree (normalize multiplier)+degree p≤D
      · simp only [if_pos hd,Option.some.injEq] at result
        subst out
        rw [meaning_mul]
        exact InClosure.mul (meaning (normalize multiplier)) (known p (List.mem_of_getElem? hp))
          ((add_le_add (meaning_degree_le _) (meaning_degree_le _)).trans hd)
      · simp [hd] at result

theorem runFrom_sound (inputs : List (Sparse r K)) (D : ℕ)
    (program : List (Step r K)) {rows out : List (Sparse r K)}
    (known : ∀ p∈rows, InClosure {q | q∈inputs.map meaning} D (meaning p))
    (checked : runFrom inputs D rows program=some out) :
    ∀ p∈out, InClosure {q | q∈inputs.map meaning} D (meaning p) := by
  induction program generalizing rows with
  | nil => simp only [runFrom,Option.some.injEq] at checked; subst out; exact known
  | cons instruction rest ih =>
      change ((evalStep inputs D rows instruction).bind (fun p =>
        runFrom inputs D (rows++[p]) rest))=some out at checked
      obtain ⟨p,hp,result⟩ := Option.bind_eq_some.mp checked
      apply ih ?_ result
      intro q hq
      rcases List.mem_append.mp hq with old | new
      · exact known q old
      · have eq : q=p := List.mem_singleton.mp new
        subst q
        exact evalStep_sound inputs D rows known instruction hp

theorem checkUnit_sound (inputs : List (Sparse r K)) (D : ℕ)
    (program : List (Step r K)) (checked : checkUnit inputs D program=true) :
    InClosure {p | p∈inputs.map meaning} D (1 : MvPolynomial (Fin r) K) := by
  cases executed : run inputs D program with
  | none => simp [checkUnit,executed] at checked
  | some rows =>
      cases last : rows.getLast? with
      | none => simp [checkUnit,executed,last] at checked
      | some p =>
          have empty : add p (scale (-1) SparseCertificate.oneData)=[] := by
            simpa [checkUnit,executed,last,List.isEmpty_iff] using checked
          have equal := congrArg meaning empty
          have unit : meaning p=1 := by
            apply sub_eq_zero.mp
            simpa only [meaning_add,meaning_scale,SparseCertificate.meaning_oneData,
              meaning,neg_one_smul,←sub_eq_add_neg] using equal
          have known := runFrom_sound inputs D program (rows:=[]) (by simp) executed
          have member : p∈rows := List.mem_of_getLast? last
          simpa [unit] using known p member

abbrev Entry (r : ℕ) (K : Type*) := Step r K × Sparse r K

def checkEntry (inputs : List (Sparse r K)) (D : ℕ) (rows : List (Sparse r K))
    (entry : Entry r K) : Bool :=
  match evalStep inputs D rows entry.1 with
  | none => false
  | some p => (add p (scale (-1) entry.2)).isEmpty

def checkEntriesFrom (inputs : List (Sparse r K)) (D : ℕ)
    (rows : List (Sparse r K)) : List (Entry r K) → Bool
  | [] => true
  | entry::rest => checkEntry inputs D rows entry &&
      checkEntriesFrom inputs D (rows++[entry.2]) rest

def checkClaimedUnit (inputs : List (Sparse r K)) (D : ℕ) (entries : List (Entry r K)) : Bool :=
  match entries.getLast? with
  | none => false
  | some entry => checkEntriesFrom inputs D [] entries &&
      (add entry.2 (scale (-1) SparseCertificate.oneData)).isEmpty

theorem checkEntriesFrom_of_each (inputs : List (Sparse r K)) (D : ℕ)
    (entries : List (Entry r K)) (rows : List (Sparse r K))
    (each : ∀ i : Fin entries.length,
      checkEntry inputs D (rows++(entries.take i.val).map Prod.snd) (entries.get i)=true) :
    checkEntriesFrom inputs D rows entries=true := by
  induction entries generalizing rows with
  | nil => rfl
  | cons entry rest ih =>
      have first := each ⟨0,by simp⟩
      have first' : checkEntry inputs D rows entry=true := by simpa using first
      have tail : checkEntriesFrom inputs D (rows++[entry.2]) rest=true := by
        apply ih
        intro i
        have next := each ⟨i.val+1,by simpa using Nat.succ_lt_succ i.isLt⟩
        simpa [List.take_succ_cons,List.append_assoc] using next
      simp only [checkEntriesFrom,first',tail,Bool.and_self]

theorem meaning_eq_of_empty_sub (p q : Sparse r K)
    (empty : add p (scale (-1) q)=[]) : meaning p=meaning q := by
  have equal := congrArg meaning empty
  apply sub_eq_zero.mp
  simpa only [meaning_add,meaning_scale,meaning,neg_one_smul,←sub_eq_add_neg] using equal

theorem checkEntry_sound (inputs : List (Sparse r K)) (D : ℕ)
    (rows : List (Sparse r K))
    (known : ∀ p∈rows, InClosure {q | q∈inputs.map meaning} D (meaning p))
    (entry : Entry r K) (checked : checkEntry inputs D rows entry=true) :
    InClosure {q | q∈inputs.map meaning} D (meaning entry.2) := by
  cases executed : evalStep inputs D rows entry.1 with
  | none => simp [checkEntry,executed] at checked
  | some p =>
      have empty : add p (scale (-1) entry.2)=[] := by
        simpa [checkEntry,executed,List.isEmpty_iff] using checked
      rw [←meaning_eq_of_empty_sub p entry.2 empty]
      exact evalStep_sound inputs D rows known entry.1 executed

theorem checkEntriesFrom_sound (inputs : List (Sparse r K)) (D : ℕ)
    (entries : List (Entry r K)) {rows : List (Sparse r K)}
    (known : ∀ p∈rows, InClosure {q | q∈inputs.map meaning} D (meaning p))
    (checked : checkEntriesFrom inputs D rows entries=true) :
    ∀ p∈rows++entries.map Prod.snd, InClosure {q | q∈inputs.map meaning} D (meaning p) := by
  induction entries generalizing rows with
  | nil => simpa using known
  | cons entry rest ih =>
      have checks : checkEntry inputs D rows entry=true ∧
          checkEntriesFrom inputs D (rows++[entry.2]) rest=true := by
        simpa only [checkEntriesFrom,Bool.and_eq_true_iff] using checked
      have member := checkEntry_sound inputs D rows known entry checks.1
      have extended : ∀ p∈rows++[entry.2], InClosure {q | q∈inputs.map meaning} D (meaning p) := by
        intro p hp
        rcases List.mem_append.mp hp with old | new
        · exact known p old
        · have equal : p=entry.2 := List.mem_singleton.mp new
          simpa [equal] using member
      have result := ih extended checks.2
      simpa [List.append_assoc] using result

theorem checkClaimedUnit_sound (inputs : List (Sparse r K)) (D : ℕ)
    (entries : List (Entry r K)) (checked : checkClaimedUnit inputs D entries=true) :
    InClosure {q | q∈inputs.map meaning} D (1 : MvPolynomial (Fin r) K) := by
  cases last : entries.getLast? with
  | none => simp [checkClaimedUnit,last] at checked
  | some entry =>
      have checks : checkEntriesFrom inputs D [] entries=true ∧
          (add entry.2 (scale (-1) SparseCertificate.oneData)).isEmpty=true := by
        simpa [checkClaimedUnit,last] using checked
      have empty : add entry.2 (scale (-1) SparseCertificate.oneData)=[] := List.isEmpty_iff.mp checks.2
      have unit : meaning entry.2=1 := by
        simpa only [SparseCertificate.meaning_oneData] using
          meaning_eq_of_empty_sub entry.2 SparseCertificate.oneData empty
      have known := checkEntriesFrom_sound inputs D entries (rows:=[]) (by simp) checks.1
      have member : entry.2∈entries.map Prod.snd := List.mem_map_of_mem (List.mem_of_getLast? last)
      simpa [unit] using known entry.2 (by simpa using member)

end ICFallDegrees.SparseDAG
