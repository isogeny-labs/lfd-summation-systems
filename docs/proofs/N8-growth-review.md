# Review of `N8-growth.md` (written growth theorem): main-session adversarial read

Date 2026-11-12. Reviewer: the main session (frontier tier), reading the draft at commit `9da01a3`
(the agent's last snapshot before this review). Independent second read: `N8-growth-review-redteam.md`
(red-team agent, same draft, written without sight of this file). Nothing in the draft is Lean-checked;
nothing here changes that. Labels as in the draft: **kernel**, **proof**, **exact**, **empirical**, **conjecture**.

## Verdict

The written proof of Theorem 5.1 (FEI_k(w) for k ≥ w + Λ(w) + 1, hence 1 ∉ V_D(k) for k ≥ k_D = D + Λ(D − 1))
**stands** on my reading, with two corrections to the text that do not change any statement, one citation that
turns the draft's "proved here" step into a quoted one, and one change of framing in the headline. Corollary 5.2
(LFD unbounded, slope ≥ 1/14 for unit instances with t ≠ 0) follows from it and the kernel bridge, with the
framing caveat below. The result is modest (first new value LFD ≥ 8 at k = 90) and, as the draft says, says
nothing about the measured slope 1/4; its content is that the last fall degree of this family cannot stay bounded.

Status after this review: **proof (reviewed, written)**. It may be cited as such once the four items under
"Required before citation" are in the draft. It must not be cited as kernel-checked.

## What I checked, step by step

| step | verdict | note |
|---|---|---|
| (1.1) lift identity | holds | `x_i^(a_i+2^j) − x_i^(a_i+2^j−N) = x_i^(a_i+2^j−N)(x_i^N − 1)`; the correction term is a polynomial because `a_i ≥ N − 2^j` on the top run |
| Prop. 2.1 (closed ⇔ FEI) | holds | uses only: 𝔭 prime, `x_i ∉ 𝔭`, `x_i^N − 1 ∉ 𝔭` (ξ_i non-constant) |
| Prop. 2.2 | holds | least-closed-subspace argument; `1 ∉ 𝔭` |
| Thm. 3.1 (a) Wronskian criterion | holds | valid for any derivation on any field, constants = ker δ |
| Thm. 3.1, kernel of δ = d/dz | **gap in the text, not in the mathematics** | the draft says ker δ is "a field between F^p and F, different from F". Hypothesis 2 is independence over F₀ = F^p, so the proof needs ker δ = F^p exactly. True: F ⊇ ker δ ⊇ F^p with [F : F^p] = p prime and δ ≠ 0, so ker δ = F^p. One sentence to add. |
| Thm. 3.1 (b), (c), (3.2), (3.3) | hold | (3.3): `δ_t(t^p) = 0` so a column of p-divisible order loses nothing; the uniformizer is separating because its order 1 is prime to p |
| Thm. 3.1 (d) | holds | Σ ord_P(W) = 0, Σ ord_P(dz) = 2γ − 2, Σ A_P = 0, −Σ μ_P = H(v); the excluded index i_P is the one of minimal order, so the column sum is A_P − μ_P |
| Thm. 3.1, novelty | **relabel** | this is the known characteristic-p form of the Mason–Stothers/Brownawell–Masser bound: the Wronskian criterion over F^p is Garcia–Voloch, "Wronskians and linear independence in fields of prime characteristic", Manuscripta Math. 59 (1987) 457–470; the S-unit sum bound in characteristic p with the "not all p-th powers" hypothesis is in Voloch, "Diagonal equations over function fields", Bol. Soc. Brasil. Mat. 16 (1985) 29–39 (and in Wang 1999, Hsia–Wang for several variables). The draft's proof is fine and may stay as a self-contained account; T2 should be labelled a known theorem with these references, not a new result. |
| Lemma 3.2 divisors | holds | poles of ξ₁ at O only, of ξ₂ at R only; `t ≠ 0` makes each a unit at the other's pole (`ξ₂(O) = x(R) = t`, `ξ₁(R) = t`); coincidences inside 𝒮 (R = ±2P₀, R 2-torsion, b = 0) only shrink 𝒮 |
| Thm. 3.3 | holds | ord_O(v_e/v_e′) ≡ −2(a − a′) mod p because ord_O(φ_e/φ_e′) ∈ pZ; the m = 2 case gives spread 0, which with `|2(a − a′)| < p` forces equal exponents, so no two distinct pairs in [0, N)² are parallel |
| Lemma 3.4 | holds | matroid of the vectors ξ^e over F₀; a dependent set contains a circuit; the circuit lies in neither independent B_i. K ⊂ k̄ = k̄^p ⊂ F^p = F₀ |
| Lemma 4.1 monotonicity | **typo, conclusion survives** | the two-step ratio `C(2ℓ+2, ≤w)/C(2ℓ, ≤w)` is bounded by `((2ℓ+1)/(2ℓ+1−w))²`, not `((2ℓ+2)/(2ℓ+2−w))²` (the one-step factor `(n+1)/(n+1−w)` decreases in n). With ℓ ≥ 4w the bound is still `(8/7)⁴ < 2` and the ratio of ρ still exceeds 1.17. |
| Lemma 4.1 numbers | holds | Λ(1) = 25 by hand (2^24 < 7776·51², 2^25 > 7776·53²); Λ(w) ≤ 13w + 23 by the stated evaluation |
| Lemma 4.2 window | holds | `|2(a − a₀)| < 2N < p` leaves only the multiples 0, ±p of p; the far case forces `a > N − σ`, `a₀ < σ` |
| Lemma 4.3 | holds | `m ≤ 36·C(2ℓ, ≤w)` from 6 × 6 block pairs and the low ℓ bits; `2^(ℓ−1) < 2σ + 1 ≤ 6m²` gives ℓ ≤ Λ(w) |
| Lemma 4.4 | holds | both digit arguments checked; needs λ − 1 ≤ j, provided in Thm. 5.1 by j ≥ k − w ≥ Λ(w) + 1 |
| Thm. 5.1 | holds | `w < k − λ` from k ≥ w + λ + 1; Lemma 3.4 needs "no circuit inside M meets both parts", which is what Lemmas 4.3/4.4 exclude |
| Cor. 5.2 | holds, **reframe** | p ≡ 1 mod 2^k and p ≠ 2^k + 1 give p ≥ 2^(k+1) + 1 > 2^(k+1); `⌊(k − 10)/14⌋ + 1` from k_D ≤ 14D + 10 needs D ≥ 2, i.e. k ≥ 38, as stated. Framing: for a fixed prime p the depth is bounded by v₂(p − 1), so "LFD unbounded in k" is a statement about the family (p, k) with 2^(k+1) | p − 1 (unit instances exist at every such depth by the kernel lemma `UnitInstance.unit_of_primitiveRoot`); the headline should say so. |
| Cor. 5.3 | holds | same locality argument at depth k + 1 with weight D; p > 2^(k+2) is the right hypothesis for Lemma 4.3 there |
| Prop. 7.1 | holds | direct |
| Sanity against data | consistent | the theorem is silent below k = 27 (D = 2) and k = 39 (D = 3); it never contradicts 1 ∈ V₇(14), 1 ∉ V₆(14) or LFD ≤ k + 2 (Λ(D − 1) ≥ 4(D − 1) gives D ≤ (k + 4)/5 < k + 2) |

Cheap attacks I tried and that fail, as they should: a direct wrap-around construction of 1 at small D and large k
(multiplying g by `x₁^(2^(k−1))` twice returns g; wrap-around alone is cheap, and the proof's content is precisely
that wrap-around *and* vanishing on the curve cannot be combined at low weight); and a Frobenius-type counterexample
to locality (`S₃(x₁^p, x₂^p)`), which the draft already handles by working over F^p and which is outside the window
`p > 2^(k+1)` anyway.

## Required before citation

1. State ker(d/dz) = F^p in the proof of Theorem 3.1 (prime-degree argument above).
2. Cite Garcia–Voloch 1987 and Voloch 1985 for Theorem 3.1 and relabel T2 from "proved here" to "known; proof included".
3. Correct the ratio in the proof of Lemma 4.1 to `(2ℓ+1)/(2ℓ+1−w)`.
4. Reframe T1/C1: unboundedness over the family of (p, k) with 2^(k+1) | p − 1, not for a fixed prime; name the
   kernel lemma that supplies the unit instances at every depth.

## Not settled by this review

* The draft's exact tables (Section 6) were not recomputed by me. The red-team read reruns the finite lemma checks
  (`ic_bounds_growth_lemmas.py lambda|window|straddle|semaev`); the closure tables at D = 3, 4 are the agent's runs
  (**exact** by its records, unverified by a second engine).
* The kernel bridge `ChainClosureBridge.chain_mem_closure_iff` (`formal/ICFallDegrees/ChainClosureBridge.lean:404`,
  audited at `Audit.lean:486`) assumes only `(2 : K) ≠ 0` and `2 ≤ D`: it is field-general, so the draft's remark that
  it is "stated for 2^k ∣ p − 1" is too cautious and Corollary 5.2 needs 2^k ∣ p − 1 only for its own unit-instance
  hypothesis, not for the bridge. (Read by me; the bridge's own correctness is the kernel's.)
* Knowledge-base check: the `research-kb` MCP server and its CLI fallback are not available in this session (no
  sibling checkout); the novelty assessment above rests on the web search of the two references. A `search_knowledge`
  pass with the negative form ("last fall degree lower bound unbounded prime chain") is still owed.
* The singular case, t = 0, more than two summands, and every slope below 1 for the upper bound remain open, as the
  draft says.

## Sources checked

Web: Garcia–Voloch, Manuscripta Math. 59 (1987) 457–470 (eudml.org/doc/155261; portal.mardi4nfdi.de Q1098891);
Voloch, Bol. Soc. Brasil. Mat. 16 (1985) 29–39 and Brownawell–Masser, Math. Proc. Cambridge Philos. Soc. 100 (1986)
(both as cited in the Mason–Stothers survey literature); Cherry, arXiv:0807.4332 (Hasse derivatives and generalized
Wronskians in arbitrary characteristic). Repository: `docs/proofs/N8-growth.md`, `docs/proofs/U1-upper-bounds.md` §1.

## Reconciliation with the independent red-team read (`N8-growth-review-redteam.md`, same draft, written blind)

Concordant on every step and on the four required changes above (its F2, F3, F4, F5, F6 are items 3, 4, 4, 2, 1
of "Required before citation"). It reran the finite checks (`lambda`, `window`, `straddle`, `semaev`: all `ok`,
counts equal to the draft's) and the D = 3 ideal table at k = 2..9 (every entry reproduced). It adds, and I adopt:

5. **Missing artifact (its F1).** The D = 4 main-instance table of Section 6.2 (p = 12289, a = 4021, b = 9917,
   t = 4408, k = 4..9) and the "33 computed cases" cite `ideal-D4-*.json`, but `experiments/ic-bounds-growth/`
   holds only the unit-instance file (k ≤ 8). Either commit the main-instance record or relabel those rows and the
   count as unrecorded. I am rerunning `ic_bounds_growth_ideal.py 12289 4021 9917 4408 4 4 9` in this session; the
   result is appended below when it finishes.
6. **Labelling (its F16).** The "circuits" row of Section 6.1 tests circuits over the constants; Theorem 3.3 in
   characteristic p is about F^p-circuits. The row is a consistency check, not a test of the statement.
7. **Slope constant (its F17).** The asymptotic slope of the proved bound is about 1/13.05, not 1/12.9; 1/14 stands.

Two cheap tests in the gap between data and theorem, neither of which can refute Theorem 5.1 but both of which bear
on the conjectured pattern FEI_k(w) ⇔ k ≥ 4w − 4 (E1) and on Theorem 3.3: FEI at (w, k) = (2, 10..12) and (3, 10..11)
on the p = 40961 instance (running now, results below), and an m = 5 F^p-circuit search at p = 151 or 157, where the
threshold 72 is below p/2 and E(F_{p²}) is enumerable (not run; a dependent pairwise-far 5-set would refute
Theorem 3.3 outright).

### Reruns by the reviewer (this session, 2026-11-12, `ulimit -v 3000000`, one thread)

* **Item 5 settled.** `python3 scripts/ic_bounds_growth_ideal.py 12289 4021 9917 4408 4 4 9` (1,028 s) reproduces
  every entry of the draft's D = 4 main-instance table: dim S_4(k) = 105, 270, 564, 1025, 1657, 2403 and
  dim V_4(k) = 163, 386, 605, 829, 1061, 1301 for k = 4..9; 1 ∈ V_4 at k = 4, 5 only; V_4 ⊂ S_4 from k = 7;
  S_4 closed ⇔ FEI_k(3) at k = 8, 9 only (k = 7 fails on 12 rows at the top shifts). The record is now committed as
  `experiments/ic-bounds-growth/ideal-D4-p12289-s12.json` (regenerated by the reviewer, not by the author agent).
* **FEI in the gap between data and theorem** (`ic_bounds_growth_fei.py 40961 21223 9887 4748 w k`, 1–8 s each):
  FEI_k(2) holds at k = 10, 11, 12 (dim S_2 = 2) and FEI_k(3) holds at k = 10, 11 (dim S_3 = 138, 142 = 4k + 98).
  Consistent with the conjectured pattern FEI_k(w) ⇔ k ≥ 4w − 4 (E1); records committed as
  `experiments/ic-bounds-growth/fei-w{2,3}-k{10,11,12}-p40961.json`. Not run: the m = 5 F^p-circuit search.

Both are **exact** on the instances named and `[pcc:unverified]` as certificates (one engine, one machine).
