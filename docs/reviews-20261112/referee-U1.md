# Referee report: U1 written proof of the selector bound (Theorem `thm:selector`)

Date: 2026-11-12. Reviewed: `docs/proofs/U1-upper-bounds.md` (sha256 `676c1c50…5113`),
`paper/curve-lfd-bounds.tex` Section 2 and lines 824–861 (`thm:selector`), lines 1300–1306
(`thm:bridge`), the Lean statements `Closure.InClosure` / `HasLFDBound`,
`ChainClosureBridge.InQuotClosure` / `chain_mem_closure_iff`, `PrimeM2General.normalized`,
`LocalReductionBound.lfd_bound_of_local` (statements read, **not built**), and afterwards
`docs/reviews-U1-upper-bound-redteam-20261110.md`. The repository was not modified. My scripts are
in `/tmp/claude-0/-home-user/c25f3f7b-b111-516d-a596-5c5b49ca2dae/scratchpad/`
(`ordlfd.py`, `sweep.py`, `lem.py`, `quot.py`, `quot2.py`, `control.py`, `qctl.py`, `ratio.py`,
`bridge.py`); none of them imports repository code.

## Verdict

**The mathematical statement of `thm:selector` (p ≥ 5, k ≥ 2, 2^k | p−1, 4a³+27b² ≠ 0, subgroup
or coset base, LFD ≤ k+2) is proved by U1 Sections 1–2.5, conditional on one thing I could not
check: that the Lean bridge `chain_mem_closure_iff` actually builds** (no `lake` in this
environment). I found no fatal gap and no fixable mathematical gap in Sections 1–2.5. Every step
was re-derived by hand; Lemma 2.2, Lemma 2.3, the bridge identity and the final bound were also
checked by independent computation, with no counterexample.

What is **not** in order is the surrounding material: the proof sketch printed in the paper does
not describe a valid argument as literally written (finding 1); the computational evidence the
paper and U1 cite for this theorem is not in the repository (finding 2); the paper's negative
control does not apply to the population it is attached to (finding 3); and U1 is byte-identical
to the text the prior red team reviewed, so none of that review's findings is resolved in U1
itself (finding 4).

## Findings

**1. Paper proof sketch is not a correct summary of the proof (gap fixable; paper,
`curve-lfd-bounds.tex:839–845`).**
(a) "so every element of the ideal supported on the fibre x₁ = α is derived within weight k+2" is
false as written: such elements e_α·b have weight up to 2k and cannot lie in V_{k+2} when
k+2 < 2k. What U1 proves is only that the *correction rows* (u_jρ − π_α(u_jρ))·e_α and the e_α
with empty fibre lie in V_{k+2}.
(b) "Summing over α gives I ∩ M_{≤D} ⊂ V_D" is not a valid step: for a ∈ I ∩ M_{≤D} the pieces
a·e_α have weight up to wt(a)+k > D, so the fibrewise statement cannot be summed inside the
budget. U1 does not do this; it runs an induction over monomials (§2.5) using the projector π
onto a complement S.
(c) The sketch never states the fact that makes the local-to-global step work: the complement
S = ⊕ span(R_α)e_α has weight ≤ k+1 = D−1 (hypotheses `bounded` and `room` of
`lfd_bound_of_local`). Without it "I ∩ M_{≤D} ⊂ V_D at D = k+2" does not imply LFD ≤ k+2.
Repair: replace the three sentences by "the reduced representatives ρ·e_α (ρ ∈ {1, x₂}) have
weight ≤ k+1; the corrections v·ρe_α − π(v·ρe_α) and 1 − π(1) lie in V_{k+2}; induction on
monomials gives m − π(m) ∈ V_{max(k+2, wt m)}".

**2. Cited evidence is absent from the repository (gap fixable; U1 lines 5, 210–260, 668–671;
paper lines 849–861 and 1623–1626).**
`experiments/ic-bounds-upper/` (U1 §2.6, §2.7, §8) and `experiments/ic-bounds-upper-redteam/`
(source of the paper's 149,122 / 132,303 / 2,570 counts) do not exist in this checkout; a
repository-wide search finds these numbers only in the paper and the review file. The paper says
"every record behind the numbers in this paper [is] in the companion repository". So U1's
"16,844 instances, all pass" and the paper's test counts are **unverifiable here**. My own
recomputation reproduces the parts I could reach (see "Verified" below) but not the totals.

**3. The negative control is misattributed (unclear exposition / overclaim; paper line 853).**
"random degree-2 generators do exceed the bound, so the test could have failed" is attached to a
test on equal-depth systems with 2 ≤ k ≤ 5. In that population random generators do **not**
exceed k+2: I found exact LFD ≤ k+2 on 200 random chain-degree-2 generators at (p,k) = (13,2)
and 24 at (17,3) in the ordinary ring, and on 30/20/12/6 more at (13,2), (17,3), (17,4), (97,4)
in the quotient model (LFD = k+2 on most, never above). The prior review says the same thing
("equal depths alone would not discriminate"); its control exceeded the bound only at unequal
depths (1,3), (1,4), (2,4), which are outside the theorem as stated in the paper. As far as I
can test, the equal-depth experiment has no working negative control. Repair: say where the
control was run, or drop the sentence.

**4. Prior red-team findings: none resolved in U1; some resolved in the paper.**
The file hash equals the hash the red team recorded, so U1 was not edited after that review.

| prior finding | status in U1 now | status in paper |
|---|---|---|
| F1 "Replaces the dimension bound 2k+1", "every prime p" (U1 line 11) | unresolved | resolved (table keeps 2k+1 as Lean, k+2 as written, nonsingular, p ≥ 5) |
| F2 §2.6/§2.7 are not independent tests; "No instance exceeds k+2" could not have come out otherwise for k ≥ 3 | unresolved (text unchanged; scripts absent, finding 2) | paper cites the red-team test instead |
| F3 optimality wording (T1′) | unchanged; acceptable as stated | resolved ("not claimed to be sharp") |
| F4 unequal depths outside the kernel bridge; k = 1 outside the Lean presentation | disclosed in §1, still not proved formally; Corollary 1 still states min(k₁,k₂)+2 without the caveat | not claimed in paper |
| F5 Propositions 4, 5 stated for D′ ≥ D but used for D′ ≥ D−1 (lines 370, 386) | unresolved; I confirm the mismatch | n/a |
| F6 "(The true value is 5.)" line 579; "(4, 5, 6 for k ≤ 13)" line 275 | unresolved | resolved (line 857) |
| F7 E2 "first depth" wording (lines 17, 545) | unresolved | n/a |
| F9, F10, F11 (C1 loci, E1 label, formalization/hygiene) | unresolved | n/a |

None of F1–F11 affects the proof of the theorem.

**5. Corollary 1 is wider than what is kernel-supported (gap fixable; U1 lines 87–94).**
Two clauses rest on an unformalized bridge or lie outside the Lean presentation: unequal depths
(the bridge is stated for `chainInputs K m k`, one depth) and k = 1 (`PrimeM2General.inputs`
needs 2 ≤ k). The mathematics in A is fine for both; the label "proof" for min(k₁,k₂)+2 should
carry "bridge not formalized". The paper states neither clause, so `thm:selector` is unaffected.

**6. Theorem 1 hypothesis is representation-dependent at k₂ = 1 (unclear exposition; U1 lines
75–80).** With N₂ = 2, x₂² = 1 in A, so g₀, g₁, g₂ and hence q_α are not determined by g; the
hypothesis is about a chosen lift. The proof still works (q = c(X²−1) gives q̄ = 0, r = 2, all
corrections zero; I checked this case). Outside the paper's range k ≥ 2.

**7. Notation clash in Lemma 2.3 (typo; U1 lines 133–136).** ρ denotes both an element of R and
"the image in B" (ρ(u_j²), ρ(X²)). Also line 208 says the last step "is
`lfd_bound_of_local` with σ = k₁+1" while §2.5 actually concludes through (1.2); both routes are
valid but the text should pick one.

**8. Paper wording (typo; paper line 808).** "written entries are proved in the companion
manuscript, in Theorem `thm:selector`" – the theorem is in this manuscript and only sketched; the
proof is in `docs/proofs/U1-upper-bounds.md`. Line 849 "A separate review checked every step" is
a provenance claim I cannot verify beyond the existence of the review file.

**9. Paper statement vs U1 (no gap).** Paper: p ≥ 5, k ≥ 2, 2^k | p−1, nonsingular, equal depths,
coset x_i = c·u_{i,0} with one c. U1 Corollary 1 covers this strictly (k ≥ 1, independent scales
s₁, s₂). Lemma 2.2 uses exactly p ≠ 2, 3 and 4a³+27b² ≠ 0. The paper's generator
(`PrimeM2General.normalized`: top −2(t²+a)u₀w₀ − 2t(u₀w₁+u₁w₀) + u₁w₁, tail t²(u₁+w₁) −
(2at+4b)(u₀+w₀) + a² − 4bt) equals the chain normal form of the S₃ used in Lemma 2.2; I checked
the coefficients.

## Edge cases examined (all covered by the proof)

- x₁ = x₂, points of order 2, target with y = 0: no step depends on them. disc_X S₃(α,X,t) =
  16 f(α) f(t) (verified symbolically), so f(t) = 0 makes every fibre polynomial a square; a
  double root in H₂ gives gcd d of degree 1 (X^N − 1 is separable), handled by case r = 1.
- t ∈ H (leading coefficient (α−t)² vanishes at α = t): fibre polynomial of degree ≤ 1, or a
  nonzero constant when f(t) = 0; handled by the deg q = 1 and deg q = 0 cases. Nonsingularity is
  used exactly here and nowhere else.
- Unit ideal vs proper ideal: (C1) puts Σ_{empty fibres} e_α in V_D; if all fibres are empty this
  is 1. No case distinction is needed; the bound is for every h ∈ I, every budget max(D, deg h),
  via (1.2) for all D′ ≥ D.
- "A = functions on H₁×H₂, I = functions vanishing on Z": true (p ∤ N, μ_N ⊂ F_p) but the proof
  only uses A = ⊕ B·e_α and B/q̄B = F_p[X]/(gcd), which I checked directly.
- Small k and p: k = 2 at p = 5 (H = F₅^*) and k = 3 at p = 17 are covered and were computed.
- Singular curves: the proof stops (zero fibre at α = t = double root ∈ H). Not claimed.

## Steps verified, and how

| step | method | result |
|---|---|---|
| (1.1), Lemma 1.1, (1.2) | by hand, against the Lean definitions of `InClosure`, `HasLFDBound`, `InQuotClosure` | correct; the budget is the literal deg g + deg h ≤ D, applied to derived rows |
| Bridge (dimension form) | ordinary-ring closure vs quotient closure + chain part, 255 (instance, D) pairs, k = 2, 3 | all equal |
| Bridge (Lean) | statement read only | **unverified build** |
| Lemma 2.1 | by hand | correct |
| Lemma 2.2 | sympy: X²-coefficient (α−t)²; at α = t linear coefficient −4f(t), constant f′(t)² − 8t f(t); brute force over all (a,b,t,α), p = 5, 7, 11, 13 | identities hold; no zero fibre on a nonsingular curve |
| Lemma 2.3 | both recursions expanded by hand (j = 0, j = k₂, k₂ = 1 included); literal closure U(q) computed for every nonzero q of degree ≤ 2 up to scaling at p = 5, 13 (k₂ ≤ 2), 17 (k₂ ≤ 4), 41 (k₂ ≤ 3), plus samples at (97,5), (193,6): 7,825 polynomials | (a) and (b) hold in all |
| Lemma 2.4, (C1), (C2), induction, conclusion | by hand | correct; operand weights ≤ k₁+1 ≤ D−1 and ≤ max(r−1, k₁+1) ≤ max(D, r)−1 |
| Theorem, ordinary ring, k = 2 | own closure engine on all monomials of degree ≤ D, D = 0..5; ideal dimension from the Gröbner structure of the chains and an explicit basis of gA; every (a,b,t) at p = 5, 13, 17 (7,235 systems, singular included) and 1,400 mixed-coset systems | no LFD above 4; nonsingular subgroup histogram 3: 416, 4: 6,336, equal to U1 §2.7; LFD = 3 iff 3t²+a = 0 without exception |
| Theorem, ordinary ring, k = 3 | same, D = 0..7 (1,716 monomials), 120 systems at p = 17, 41: random, t ∈ H, t ∈ H with f(t) = 0, 3t²+a = 0, coset ratio 3, singular with vanishing fibre | LFD = 4 on all nonsingular (104), ≤ 4 on singular; bound 5 never reached |
| Theorem, quotient model, k = 2, 3, 4 | V_D = I_D for D = k+2..2k via the bridge; 53 + 59 + 53 + 29 nonsingular systems at (13,2), (17,3), (17,4), (97,4) | no violation; max LFD 4 |
| U1 §2.7 witness p = 13, y² = x³+x+2, t = 3 | ordinary dims 5, 25, 70, 126 against ideal 15, 35, 70, 126 | LFD = 4 confirmed |

Not verified: the Lean build and axioms; U1 §2.6 and §2.7 as run by the authors (scripts
absent); anything at k ≥ 5 (my k = 5 run did not finish in the time available); the paper's
totals 149,122, 132,303 and 2,570; Sections 3–7 of U1 beyond confirming prior finding F5.
A bug in my own Lemma 2.3 checker (duplicate dictionary key) produced false failures on the
first run; it was fixed and the rerun is the one reported.
