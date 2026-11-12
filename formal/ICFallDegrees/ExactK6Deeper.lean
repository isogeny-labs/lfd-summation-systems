import ICFallDegrees.ExactK6InstanceP1153
import ICFallDegrees.DepthRestrictionM2

/-! The depth-6 exclusion `1 ∉ V_4` of the instance `p = 1153`, `(a, b, t) = (139, 1094, 185)`
persists at every deeper depth, by the depth restriction.  Together with the kernel-checked
exact value at depth 6 this gives, for one concrete system, `1 ∉ V_4(k)` for every `k ≥ 6`,
and `LFD ≥ 5` at every depth `k ≥ 6` at which the ideal is the unit ideal. -/

namespace ICFallDegrees.ExactK6Deeper
open ICFallDegrees ExactK6InstanceP1153

/-- `1` is not in the degree-4 closure at any depth `k ≥ 6`. -/
theorem not_one_in_closure_four_of_ge (k : ℕ) (hk6 : 6 ≤ k) (hk2 : 2 ≤ k) :
    ¬ InClosure (Set.range (PrimeM2General.inputs k hk2 (139 : ZMod 1153) 1094 185)) 4 1 :=
  DepthRestrictionM2.not_one_of_le 6 ExactK6InstanceP1153.hk 139 1094 185 4 (by norm_num)
    ExactK6InstanceP1153.not_one_in_closure_four k hk6 hk2

/-- At every depth `k ≥ 6` where the ideal is the unit ideal, no `e ≤ 4` is an LFD bound. -/
theorem lfd_ge_five_of_unit (k : ℕ) (hk6 : 6 ≤ k) (hk2 : 2 ≤ k)
    (unit : (1 : FiniteQuotient.P (2*k) (ZMod 1153)) ∈
      Ideal.span (Set.range (PrimeM2General.inputs k hk2 (139 : ZMod 1153) 1094 185))) :
    ∀ e ≤ 4, ¬ HasLFDBound (Set.range (PrimeM2General.inputs k hk2 (139 : ZMod 1153) 1094 185)) e :=
  DepthRestrictionM2.lfd_gt_of_le 6 ExactK6InstanceP1153.hk 139 1094 185 4 (by norm_num)
    ExactK6InstanceP1153.not_one_in_closure_four k hk6 hk2 unit

end ICFallDegrees.ExactK6Deeper
