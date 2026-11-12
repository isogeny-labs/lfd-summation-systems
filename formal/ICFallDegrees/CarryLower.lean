import ICFallDegrees.QuotientLowerBound
import ICFallDegrees.ChainQuotient

noncomputable section

/-!
# The carry lower bound (N1)

Let `G` be a reduced-ideal-trivial subsystem (for example the squaring chains
`ChainQuotient.chainInputs K m k`, any `m`, `k`, characteristic other than two)
and let `g` be a reduced polynomial (a chain normal form) of total degree `d`.
Suppose a multiple `a * g` of `g` has a reduced representative `h` of degree
at most `d` (for the chains: `h = NF(u * g)` for a chain variable `u` whose
carries do not raise the degree) which is not a scalar multiple of `g`.
Then no `e ≤ d` is a last-fall-degree bound of `G ∪ {g}`: `LFD ≥ d + 1`.

Proof (Proposition 1 of `experiments/ic-bounds-m3/README.md`). The degree-`d`
closure of `G ∪ {g}` lies in the quotient span `{p | ∃ c, p - c*g ∈ (G)}`
(`QuotientLowerBound.closure_le_quotient_span`: a row of degree `< d` in that
span already lies in the ideal of `G`, because its reduced remainder has
degree `< d` while a nonzero multiple of `g` has degree `d`). The element
`h` lies in the ideal of `G ∪ {g}` and has degree at most `d`, but it is not in
the quotient span: `h - c*g ∈ (G)` with both sides reduced forces `h = c*g`.

Compared with the written proof, the top-degree-form argument for step (ii)
is replaced by reduction with the literal degree budget
(`FiniteQuotient.exists_reduced_remainder`) and the uniqueness of reduced
representatives; no Groebner standard-representation statement is needed, and
no lower bound `2 ≤ d` is used.
-/

namespace ICFallDegrees.CarryLower

open MvPolynomial FiniteQuotient

variable {r : ℕ} {K : Type*} [Field K]

/-- General carry lower bound for a reduced-ideal-trivial subsystem `G`. -/
theorem not_hasLFDBound_of_carry
    (G : Set (P r K)) (nonzero : (0 : P r K)∉G) (zero : ReducedIdealTrivial G)
    (g h a : P r K) (d : ℕ)
    (g_reduced : g∈reducedSpace G) (g_degree : g.totalDegree=d)
    (h_reduced : h∈reducedSpace G) (h_degree : h.totalDegree≤d)
    (carry : h-a*g∈Ideal.span G)
    (not_multiple : ∀ c : K, h≠c • g) :
    ∀ e≤d, ¬HasLFDBound (G∪{g}) e := by
  intro e he bound
  -- The degree-`d` closure lies in the quotient span of `g`.
  have contained : closureSpace (G∪{g}) d≤quotientSpan (Ideal.span G) g := by
    apply closure_le_quotient_span (G∪{g}) G nonzero zero g g_reduced d g_degree
    rintro p (hp | hp) _
    · exact ⟨0,by simpa using Ideal.subset_span hp⟩
    · rw [Set.mem_singleton_iff.mp hp]
      exact ⟨1,by simp⟩
  -- `h` lies in the ideal of the system.
  have member : h∈Ideal.span (G∪{g}) := by
    have hg : g∈Ideal.span (G∪{g}) := Ideal.subset_span (Or.inr rfl)
    have mono : Ideal.span G≤Ideal.span (G∪{g}) := Ideal.span_mono Set.subset_union_left
    have identity : h=(h-a*g)+a*g := by ring
    rw [identity]
    exact (Ideal.span (G∪{g})).add_mem (mono carry) ((Ideal.span (G∪{g})).mul_mem_left a hg)
  have derived := bound h member
  have closed : InClosure (G∪{g}) d h := derived.mono (max_le he h_degree)
  obtain ⟨c,hc⟩ := contained (show h∈closureSpace (G∪{g}) d from closed)
  have scaled : C c*g∈reducedSpace G := by
    simpa only [smul_eq_C_mul] using (reducedSpace G).smul_mem c g_reduced
  have same := reduced_unique_of_reduced_zero G zero h_reduced scaled hc
  exact not_multiple c (by rw [same,smul_eq_C_mul])

open ChainQuotient in
/-- **Carry lemma (N1)** for squaring-chain factor bases. Let `C` be the chain
equations of `m` summands at depth `k` (`X_v^2 - X_next(v)`, last `X^2 - 1`) over a
field with `2 ≠ 0`, `g` a chain normal form of total degree `d`, and `h` a chain
normal form of degree at most `d` congruent to `X u * g` modulo the chains, not a
scalar multiple of `g`. Then `LFD(C ∪ {g}) ≥ d+1`: no `e ≤ d` is an LFD bound. -/
theorem chain_carry_lower (m k : ℕ) (two_ne : (2 : K)≠0)
    (g h : P (m*k) K) (u : Fin (m*k)) (d : ℕ)
    (g_reduced : g∈reducedSpace (chainInputs K m k)) (g_degree : g.totalDegree=d)
    (h_reduced : h∈reducedSpace (chainInputs K m k)) (h_degree : h.totalDegree≤d)
    (carry : h-X u*g∈Ideal.span (chainInputs K m k))
    (not_multiple : ∀ c : K, h≠c • g) :
    ∀ e≤d, ¬HasLFDBound (chainInputs K m k∪{g}) e := by
  apply not_hasLFDBound_of_carry (chainInputs K m k) (chainInputs_nonzero m k)
    (reduced_chain_ideal_trivial m k two_ne) g h (X u) d g_reduced g_degree
    h_reduced h_degree carry not_multiple

end ICFallDegrees.CarryLower
