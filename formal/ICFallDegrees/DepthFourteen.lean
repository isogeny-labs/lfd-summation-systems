import ICFallDegrees.UnitInstance

/-!
# Consequences of the recorded closures at depths 10 to 14, at every depth

The closures of the unit instances at depths `k = 10, …, 14` and degrees `5, 6, 7` are
**recorded computations** (`experiments/ic-bounds-k14/RESULTS.md`): exact runs of a closure
engine, *not* checked by the Lean kernel.  They enter this file only as explicit hypotheses of
theorems, named

* `h10 : ¬ InClosure (inputs 10 …) 5 1`  (recorded: `1 ∉ V_5(10)`),
* `h13 : InClosure (inputs 13 …) 6 1`    (recorded: `1 ∈ V_6(13)`),
* `h14 : ¬ InClosure (inputs 14 …) 6 1`  (recorded: `1 ∉ V_6(14)`),
* `h7 : InClosure (inputs 14 …) 7 1`     (recorded: `1 ∈ V_7(14)`, obtained by a degree lift
  in one run, without an independent certificate).

Nothing here kernel-checks these four statements.  What the kernel checks is every step from
them to the stated conclusions: the depth restriction (`DepthRestrictionM2`), and, where a
conclusion needs the unit ideal, the no-solution lemma and the finite-grid Nullstellensatz of
`UnitInstance`.

The general statements hold for every field, curve `(a, b)` and target `t`; the instance
statements are for `y² = x³ + x` over `𝔽_65537` with target `(4271, 11720)`, and for
`y² = x³ + 3x` over `𝔽_163841` with target `(73730, 28291)`.
-/

noncomputable section

namespace ICFallDegrees.DepthFourteen

open PrimeM2General

/-! ## General field, curve and target -/

section General

variable {K : Type*} [Field K] (a b t : K)

/-- Membership of `1` in the degree-`D` closure descends to every smaller depth. -/
theorem one_of_le (D : ℕ) (hD : 2 ≤ D) (k : ℕ) (depth : 2 ≤ k)
    (one : InClosure (Set.range (inputs k depth a b t)) D 1) :
    ∀ (k' : ℕ) (depth' : 2 ≤ k'), k' ≤ k → InClosure (Set.range (inputs k' depth' a b t)) D 1 := by
  induction k, depth using Nat.le_induction with
  | base =>
      intro k' depth' hk
      obtain rfl : k' = 2 := by omega
      exact one
  | succ n hn ih =>
      intro k' depth' hk
      rcases Nat.lt_or_ge k' (n+1) with small | large
      · exact ih (DepthRestrictionM2.one_of_deeper n hn a b t D hD one) k' depth' (by omega)
      · obtain rfl : k' = n+1 := by omega
        exact one

/-- From the recorded `1 ∉ V_6(14)` (hypothesis `h14`, not kernel-checked): `1 ∉ V_6(k)` for
every `k ≥ 14`. -/
theorem not_one_six_of_ge
    (h14 : ¬ InClosure (Set.range (inputs 14 (by norm_num) a b t)) 6 1)
    (k : ℕ) (hk : 14 ≤ k) (depth : 2 ≤ k) :
    ¬ InClosure (Set.range (inputs k depth a b t)) 6 1 :=
  DepthRestrictionM2.not_one_of_le 14 (by norm_num) a b t 6 (by norm_num) h14 k hk depth

/-- From the recorded `1 ∉ V_6(14)` (hypothesis `h14`, not kernel-checked): at every depth
`k ≥ 14` at which the ideal is the unit ideal, no `e ≤ 6` is a last-fall-degree bound. -/
theorem lfd_gt_six_of_unit
    (h14 : ¬ InClosure (Set.range (inputs 14 (by norm_num) a b t)) 6 1)
    (k : ℕ) (hk : 14 ≤ k) (depth : 2 ≤ k)
    (unit : (1 : P k K) ∈ Ideal.span (Set.range (inputs k depth a b t))) :
    ∀ e ≤ 6, ¬ HasLFDBound (Set.range (inputs k depth a b t)) e :=
  DepthRestrictionM2.lfd_gt_of_le 14 (by norm_num) a b t 6 (by norm_num) h14 k hk depth unit

/-- From the recorded `1 ∉ V_5(10)` (hypothesis `h10`, not kernel-checked): `1 ∉ V_5(k)` for
every `k ≥ 10`. -/
theorem not_one_five_of_ge
    (h10 : ¬ InClosure (Set.range (inputs 10 (by norm_num) a b t)) 5 1)
    (k : ℕ) (hk : 10 ≤ k) (depth : 2 ≤ k) :
    ¬ InClosure (Set.range (inputs k depth a b t)) 5 1 :=
  DepthRestrictionM2.not_one_of_le 10 (by norm_num) a b t 5 (by norm_num) h10 k hk depth

/-- From the recorded `1 ∈ V_6(13)` (hypothesis `h13`, not kernel-checked): `1 ∈ V_6(k)` for
every `2 ≤ k ≤ 13`. -/
theorem one_six_of_le
    (h13 : InClosure (Set.range (inputs 13 (by norm_num) a b t)) 6 1)
    (k : ℕ) (depth : 2 ≤ k) (hk : k ≤ 13) :
    InClosure (Set.range (inputs k depth a b t)) 6 1 :=
  one_of_le a b t 6 (by norm_num) 13 (by norm_num) h13 k depth hk

/-- From the recorded `1 ∈ V_6(13)` and `1 ∉ V_5(10)` (hypotheses `h13`, `h10`, not
kernel-checked): the last fall degree is exactly `6` at every depth `10 ≤ k ≤ 13`.  No separate
unit-ideal hypothesis is needed: `h13` already places `1` in the ideal at these depths. -/
theorem lfd_eq_six
    (h13 : InClosure (Set.range (inputs 13 (by norm_num) a b t)) 6 1)
    (h10 : ¬ InClosure (Set.range (inputs 10 (by norm_num) a b t)) 5 1)
    (k : ℕ) (depth : 2 ≤ k) (lower : 10 ≤ k) (upper : k ≤ 13) :
    IsLastFallDegree (Set.range (inputs k depth a b t)) 6 :=
  exactLFD_of_first_one (by norm_num) (one_six_of_le a b t h13 k depth upper)
    (not_one_five_of_ge a b t h10 k lower depth)

/-- From the recorded `1 ∈ V_7(14)` and `1 ∉ V_6(14)` (hypotheses `h7`, `h14`, not
kernel-checked; `h7` comes from a degree lift without an independent certificate): the last
fall degree is exactly `7` at depth `14`. -/
theorem lfd_eq_seven
    (h7 : InClosure (Set.range (inputs 14 (by norm_num) a b t)) 7 1)
    (h14 : ¬ InClosure (Set.range (inputs 14 (by norm_num) a b t)) 6 1) :
    IsLastFallDegree (Set.range (inputs 14 (by norm_num) a b t)) 7 :=
  exactLFD_of_first_one (by norm_num) h7 h14

end General

/-! ## Prime fields with a kernel-checked unit ideal -/

section PrimeField

/-- From the recorded `1 ∉ V_6(14)` (hypothesis `h14`, not kernel-checked), for a prime field
instance `y² = x³ + a x` whose target abscissa `t` is a non-square with `t³ + a t = Y² ≠ 0`: at
every depth `k ≥ 14` with `2^(k+1) ∣ p - 1`, no `e ≤ 6` is a last-fall-degree bound.  The unit
ideal is kernel-checked (`UnitInstance.unit_zmod`). -/
theorem lfd_gt_six_zmod (p : ℕ) [Fact p.Prime] {a t Y : ZMod p}
    (curve : Y^2 = t^3+a*t) (hY : Y ≠ 0) (nonsquare : ¬ IsSquare t)
    (h14 : ¬ InClosure (Set.range (inputs 14 (by norm_num) a 0 t)) 6 1)
    (k : ℕ) (hk : 14 ≤ k) (depth : 2 ≤ k) (hp : 2^(k+1) ∣ p-1) :
    ∀ e ≤ 6, ¬ HasLFDBound (Set.range (inputs k depth a 0 t)) e :=
  lfd_gt_six_of_unit a 0 t h14 k hk depth
    (UnitInstance.unit_zmod p k depth hp curve hY nonsquare)

end PrimeField

/-! ## The instance `p = 65537`, `y² = x³ + x`, target `(4271, 11720)` -/

section P65537

open UnitInstance

/-- The chain system of the unit instance over `𝔽_65537` at depth `k`. -/
abbrev system (k : ℕ) (depth : 2 ≤ k) : Set (P k (ZMod 65537)) :=
  Set.range (inputs k depth (1 : ZMod 65537) 0 4271)

/-- Recorded `1 ∉ V_6(14)` (hypothesis, not kernel-checked) gives `1 ∉ V_6(k)` for every
`k ≥ 14`. -/
theorem p65537_not_one_six
    (h14 : ¬ InClosure (system 14 (by norm_num)) 6 1)
    (k : ℕ) (hk : 14 ≤ k) (depth : 2 ≤ k) : ¬ InClosure (system k depth) 6 1 :=
  not_one_six_of_ge 1 0 4271 h14 k hk depth

/-- Recorded `1 ∉ V_6(14)` (hypothesis, not kernel-checked) gives: no `e ≤ 6` is a
last-fall-degree bound at depths `14` and `15`.  The unit ideal at these depths is
kernel-checked (`UnitInstance.unit_p65537`). -/
theorem p65537_lfd_gt_six
    (h14 : ¬ InClosure (system 14 (by norm_num)) 6 1)
    (k : ℕ) (hk : 14 ≤ k) (hk' : k ≤ 15) (depth : 2 ≤ k) :
    ∀ e ≤ 6, ¬ HasLFDBound (system k depth) e :=
  lfd_gt_six_of_unit 1 0 4271 h14 k hk depth (unit_p65537 k depth hk')

/-- Recorded `1 ∉ V_5(10)` (hypothesis, not kernel-checked) gives `1 ∉ V_5(k)` for every
`k ≥ 10`. -/
theorem p65537_not_one_five
    (h10 : ¬ InClosure (system 10 (by norm_num)) 5 1)
    (k : ℕ) (hk : 10 ≤ k) (depth : 2 ≤ k) : ¬ InClosure (system k depth) 5 1 :=
  not_one_five_of_ge 1 0 4271 h10 k hk depth

/-- Recorded `1 ∈ V_6(13)` (hypothesis, not kernel-checked) gives `1 ∈ V_6(k)` for every
`2 ≤ k ≤ 13`. -/
theorem p65537_one_six
    (h13 : InClosure (system 13 (by norm_num)) 6 1)
    (k : ℕ) (depth : 2 ≤ k) (hk : k ≤ 13) : InClosure (system k depth) 6 1 :=
  one_six_of_le 1 0 4271 h13 k depth hk

/-- Recorded `1 ∈ V_6(13)` and `1 ∉ V_5(10)` (hypotheses, not kernel-checked) give the exact
last fall degree `6` at every depth `10 ≤ k ≤ 13`. -/
theorem p65537_lfd_eq_six
    (h13 : InClosure (system 13 (by norm_num)) 6 1)
    (h10 : ¬ InClosure (system 10 (by norm_num)) 5 1)
    (k : ℕ) (depth : 2 ≤ k) (lower : 10 ≤ k) (upper : k ≤ 13) :
    IsLastFallDegree (system k depth) 6 :=
  lfd_eq_six 1 0 4271 h13 h10 k depth lower upper

/-- Recorded `1 ∈ V_7(14)` (degree lift, no independent certificate) and `1 ∉ V_6(14)`
(hypotheses, not kernel-checked) give the exact last fall degree `7` at depth `14`. -/
theorem p65537_lfd_eq_seven
    (h7 : InClosure (system 14 (by norm_num)) 7 1)
    (h14 : ¬ InClosure (system 14 (by norm_num)) 6 1) :
    IsLastFallDegree (system 14 (by norm_num)) 7 :=
  lfd_eq_seven 1 0 4271 h7 h14

/-- All recorded computations of the unit instance over `𝔽_65537` together (hypotheses, not
kernel-checked): exclusion at degree 5 from depth 10 on, exclusion at degree 6 from depth 14
on, exact last fall degree 6 at depths 10 to 13 and 7 at depth 14, and no bound `≤ 6` at
depth 15. -/
theorem p65537_summary
    (h10 : ¬ InClosure (system 10 (by norm_num)) 5 1)
    (h13 : InClosure (system 13 (by norm_num)) 6 1)
    (h14 : ¬ InClosure (system 14 (by norm_num)) 6 1)
    (h7 : InClosure (system 14 (by norm_num)) 7 1) :
    (∀ (k : ℕ) (depth : 2 ≤ k), 10 ≤ k → ¬ InClosure (system k depth) 5 1) ∧
    (∀ (k : ℕ) (depth : 2 ≤ k), 14 ≤ k → ¬ InClosure (system k depth) 6 1) ∧
    (∀ (k : ℕ) (depth : 2 ≤ k), 10 ≤ k → k ≤ 13 → IsLastFallDegree (system k depth) 6) ∧
    IsLastFallDegree (system 14 (by norm_num)) 7 ∧
    (∀ e ≤ 6, ¬ HasLFDBound (system 15 (by norm_num)) e) :=
  ⟨fun k depth hk => p65537_not_one_five h10 k hk depth,
    fun k depth hk => p65537_not_one_six h14 k hk depth,
    fun k depth lower upper => p65537_lfd_eq_six h13 h10 k depth lower upper,
    p65537_lfd_eq_seven h7 h14,
    p65537_lfd_gt_six h14 15 (by norm_num) (by norm_num) (by norm_num)⟩

end P65537

/-! ## The instance `p = 163841`, `y² = x³ + 3x`, target `(73730, 28291)` -/

section P163841

open UnitInstance

/-- The chain system of the unit instance over `𝔽_163841` at depth `k`. -/
abbrev system' (k : ℕ) (depth : 2 ≤ k) : Set (P k (ZMod 163841)) :=
  Set.range (inputs k depth (3 : ZMod 163841) 0 73730)

/-- Recorded `1 ∉ V_6(14)` over `𝔽_163841` (hypothesis, not kernel-checked) gives: no `e ≤ 6`
is a last-fall-degree bound at depth `14`, where the unit ideal is kernel-checked
(`UnitInstance.unit_p163841`), and `1 ∉ V_6(k)` for every `k ≥ 14`. -/
theorem p163841_lower
    (h14 : ¬ InClosure (system' 14 (by norm_num)) 6 1) :
    (∀ (k : ℕ) (depth : 2 ≤ k), 14 ≤ k → ¬ InClosure (system' k depth) 6 1) ∧
    (∀ e ≤ 6, ¬ HasLFDBound (system' 14 (by norm_num)) e) :=
  ⟨fun k depth hk => not_one_six_of_ge 3 0 73730 h14 k hk depth,
    lfd_gt_six_of_unit 3 0 73730 h14 14 le_rfl (by norm_num)
      (unit_p163841 14 (by norm_num) le_rfl)⟩

/-- Recorded `1 ∈ V_6(13)` and `1 ∉ V_5(13)` over `𝔽_163841` (hypotheses, not kernel-checked)
give the exact last fall degree `6` at depth `13`. -/
theorem p163841_lfd_eq_six
    (h13 : InClosure (system' 13 (by norm_num)) 6 1)
    (h5 : ¬ InClosure (system' 13 (by norm_num)) 5 1) :
    IsLastFallDegree (system' 13 (by norm_num)) 6 :=
  exactLFD_of_first_one (by norm_num) h13 h5

end P163841

end ICFallDegrees.DepthFourteen
