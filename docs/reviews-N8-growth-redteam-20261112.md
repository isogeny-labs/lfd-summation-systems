# Red-team review of `docs/proofs/N8-growth.md` (2026-11-12)

Reviewer: breakthrough-tier red team, a separate session that neither produced nor approved N8.
Evidence labels as in the repository: **kernel**, **proof**, **exact**, **empirical**, **conjecture**.
The task was to find errors. The labels recommended here are recommendations only; the coordinator decides them.

## 0. Reviewed artifact (frozen by the reviewer)

The dispatch did not name a commit, so I froze the artifact myself:

* `docs/proofs/N8-growth.md`: git blob `6c98619ec725e856006bbb0991578414c09686b6`, SHA-256
  `51c303185ee1369b1c7f3ecedcf86c65b2b8a4464f1af82ca3ea21fa22cf79bb`. It was last changed in commit `9da01a31`,
  a WIP snapshot marked "not imported, not audited". During the review the branch head moved from `3c11e63d` to
  `9258f78b` (more snapshot commits), but this file did not change. **This review applies to that blob only.**
* `scripts/ic_bounds_growth_*.py` and `experiments/ic-bounds-growth/*.json` at the same head.
  `ideal-D4-p12289-s12.json` was untracked when the review began and is tracked now.
* Background files read: `formal/ICFallDegrees/Closure.lean`, `PrimeM2General.lean`, `ChainClosureBridge.lean`,
  `SelectorExtension.lean`, `docs/proofs/U1-upper-bounds.md` §1, `experiments/ic-bounds-k14/RESULTS.md`,
  `experiments/ic-bounds-freeext/table.md` and `records/D4-p12289-t0-X0.json`.

No existing file was edited and no Lean was built. My own test scripts are in the session scratchpad
(paths in §5). They are not part of the repository.

## 1. Verdict per claim

| claim | label claimed | recommended label | fatal error found | main caveat |
|---|---|---|---|---|
| T1 (Theorem 5.1) | proof | **proof** | none | The only non-elementary input is Theorem 3.1 in characteristic p. I checked it line by line and it is correct, but it is not Lean-checked (N8 §8 names the gap). |
| C1 (Corollary 5.2) | proof, given the kernel bridge | **proof** (given T1; the bridge itself is kernel) | none | The hypotheses are stated inaccurately: the bridge needs only `(2 : K) ≠ 0`, and the hypothesis the proof actually uses is p > 2^(k+1). The table row omits it (F4). |
| T4 (Corollary 5.3) | proof | **proof** | none | none |
| Prop 7.1 | proof | **proof** | none | elementary |
| Prop 7.2 | proof | **proof** | none | Part (c) is conditional. The "Reduction B is equivalent" sentence needs one more hypothesis (F11). |

**Is the artifact set sufficient?** Yes for review at the **proof** level: every step is written out, the bridge
statements are in kernel files, and the scripts and outputs are present. No for **kernel** level: no step of N8 has
Lean. T1 is therefore not a kernel closure. Its explicit gap is Theorem 3.1/3.3 (N8 §8).

**Consistency.** No N8 statement contradicts:
* the kernel bound `SelectorExtension.lfd_le_depth_add_two_of_two_ne` (LFD ≤ k + 2 for the same
  `inputs k depth a b t`, over any field with 2 ≠ 0);
* `exact_lfd_four_at_depth_two`;
* the recorded exact data, including 1 ∈ V_D(4D−11) and 1 ∉ V_D(4D−10) for D = 3..6 and LFD ≥ 7 at k = 14;
* N8's own refuted list.

No intermediate lemma is claimed in a range where the data refute it. FEI is proved only for k ≥ w + Λ(w) + 1, and
the data show FEI_k(w) failing only below 4w − 4. Theorem 3.3 and Lemma 3.4 hold at every depth but say nothing
about 1 ∈ V.

## 2. Findings

Severities: fatal to T1 / fatal to a secondary claim / repairable gap / presentation. **No finding is fatal.**

| # | location | severity | finding | repair; cheapest overturning test |
|---|---|---|---|---|
| F1 | Lemma 4.1 proof, lines 230–231 | repairable gap (minor) | The bound C(2ℓ+2,≤w)²/C(2ℓ,≤w)² ≤ ((2ℓ+2)/(2ℓ+2−w))⁴ does not follow from the binomial ratio quoted before it. Applying that ratio twice gives the factor (2ℓ+1)/(2ℓ+1−w)·(2ℓ+2)/(2ℓ+2−w), and the square of this exceeds the bound as written. The stated inequality held in every case I tested (w < 30, ℓ ∈ [4w, 4w+200)), but the proof does not establish it. | Use [(2ℓ+1)(2ℓ+2)/((2ℓ+1−w)(2ℓ+2−w))]². Each factor is < 8/7 for ℓ ≥ 4w, so the product is < (8/7)⁴ < 2, and ρ(ℓ+1)/ρ(ℓ) ≥ 2(7/8)⁴ still holds. Λ, k_D and the 13w+23 bound are unchanged (recomputed, §5). |
| F2 | §3.1 lines 120–123; T2 row (line 14) | presentation (attribution) | Brownawell–Masser Theorem B, as quoted in arXiv:2201.10964 (Thm 5, "follows directly from Theorem B"), assumes **no vanishing proper subsum**. It bounds max_i H(u_i) by C(k,2)·(\|S\| + max(0, 2g−2)) for 1 + u₁ + … + u_k = 0. N8's Theorem 3.1 instead assumes linear independence, bounds the projective height, and uses 2γ−2 in place of max(0, 2γ−2). For γ = 1, dividing the relation by v_{e′} (u_i = v_i/v_{e′}, i ≠ e′; a circuit has no vanishing proper subsum) turns the quoted BM statement into exactly the pairwise bound H(v_e/v_{e′}) ≤ 3(m−1)(m−2) that Theorem 3.3 uses. **So the constant 6(m−1)(m−2) agrees with the literature in characteristic 0.** The characteristic-p version (independence over K^p, ordinary derivatives) is attributed in the literature to Voloch (1985). It is not new, though the "orders ≡ 0 mod p outside 𝒮" hypothesis is a harmless generalization. I could not retrieve Voloch's exact statement. | Cite BM Theorem B in its actual form, and cite Voloch for characteristic p. Nothing downstream (Λ, k_D) changes. |
| F3 | §6.3 item 1, lines 461–465 | presentation | The fallback "for fixed (D, k), for all but finitely many p" is true for the **generic** system over F_p(a, b, t). It is false as stated for every specialized instance (a, b, t) ∈ F_p³: specialization can change rank decisions on a Zariski-closed set of parameters. | The point is moot because the characteristic-p proof stands. Reword it or delete it. |
| F4 | C1 row (line 13); Cor 5.2 (line 298); §6.3 last paragraph (lines 475–476) | presentation | **The bridge does not need 2^k ∣ p − 1.** `ChainClosureBridge.prime_m2_lfd_gt` / `not_one_of_chain_witness` need only `(2 : K) ≠ 0`, and the witness `W` = image of S_D(k) satisfies their hypotheses. The hypothesis the proof uses is char 0 or p > 2^(k+1). As written, the C1 row ("unit ideal, t ≠ 0 and p ≠ 2^k + 1") is correct only under the repository convention 2^k ∣ p − 1, which it does not state. Side remark: in C1's range k ≥ 27, p = 2^k + 1 prime forces k = 2^n with n ≥ 5, and F_5, …, F_32 are composite. **So the excluded case is empty for all 27 ≤ k < 2^33.** "Unbounded" also silently needs unit-ideal instances at every depth. They exist: each pair (α, β) ∈ μ_N² kills at most 2 targets t (by hand from the S_3 formula; S_3(α, α, t) = f′(α)² when f(α) = 0), so every prime p ≡ 1 mod 2^k with p > 2^(2k+1) + 2 has a unit-ideal target. The generic system over K(a, b, t) is another source. | State C1 for every field of characteristic 0 or p > 2^(k+1). Add the one-line existence argument. |
| F5 | §6 intro, lines 352–354 | presentation (false as written) | "dim V_4(k) = 163, 386, 605, 829, 1061, 1301 for k = 4..9 with low profile (0, 0, 2, 4k+87)" is false for k = 4, 5. There V_4 = M_≤4, with profiles (1, 9, 37, 93, 163) and (1, 11, 56, 176, 386) in `ideal-D4-p12289-s12.json`. The formula holds for k = 6..9. | Restrict the claim to k ≥ 6. |
| F6 | E4 row (line 20); §6.2 bullets (lines 429–436) | presentation | Only k = 10 is a low entry of a V_5 profile: p_4(V_5(10)) = 3163 = dim S_4(10). For k = 11 and 12 the profiles have p_4 = 3441 and 3727, while dim S_4 = 3569 and 3855. The row's wording reads as if all three were profile entries. The compared quantities also come from different instances: S_4 on p = 40961, S_3 on p = 12289, V_6(14) on p = 65537 and 163841. The **empirical** label is right. | Reword. Name the instance of each number. |
| F7 | §6.1 rows "Theorem 3.3 on computed circuits" and "Theorem 3.3 in characteristic p"; §6.3 item 1 (b), (c) | repairable gap (test coverage; does not affect the proof) | **Every script tests (in)dependence over the constants** (ranks of evaluation matrices), never over F₀ = F^p. The "control" S_3(x₁^p, x₂^p) is itself a constant-coefficient relation. So the one place where the characteristic-p argument differs from characteristic 0 was never exercised: coefficients in F^p, including non-constant ones, and condition (3.3). I built a control that shows the gap. Take the twisted S_3 relation with exponents e + p·u_e: my Wronskian test finds it dependent over F^p, while its rank over the constants is 9 of 9. **The N8 tests cannot see such dependences.** | Add an F^p (Hasse–Wronskian) test. Mine (§5.3) found no violation of Theorem 3.3 for m = 3 (exhaustive), 5, 6, 7, 8, including curve coefficients in F_{p²} ∖ F_p. |
| F8 | Refuted list row 2 (line 567); D = 4, k = 7 table cell (line 406) | presentation | "12 rows fail" depends on the basis. My independent basis gave 23/22 failing rows at D = 4, k = 6, against 99/28 in the record. At k = 7 both bases give 12/12. The verdicts agree. | Report the rank of the obstruction (the dimension of the image of the top-run map modulo 𝔭). |
| F9 | §0 lines 24–27 | presentation | The cited record `D4-p12289-t0-X0` is y² = x³ + x with t = 0, so b = 0 and R = (0,0) ∈ E[2]. There S_3(x₁, x₂, 0) = (x₁x₂ − 1)² (checked with `s3_coeffs`) and ξ₁ξ₂ = 1 (by hand: x(P + T) = a/x(P) for T = (0,0) on y² = x³ + ax). The two-term elements x₁^(2^j) x₂^(2^j) − 1 ∈ 𝔭 have weight 2 and violate the m = 2 case of Theorem 3.3, and FEI_k(2) fails at every k. This shows the method really needs to exclude this case, but the record is evidence about (t = 0, b = 0), not t = 0 alone. N8's own run at t = 0 with b ≠ 0 (`ideal-D3-p12289-t0.json`) is the relevant evidence. | Say so. Cite the t = 0, b ≠ 0 run instead. |
| F10 | C1 row "about 1/12.9" | presentation | The asymptotic slope is 1/(1+c), where 4·H₂(1/(2c)) = 1 gives c ≈ 11.99. That is about 1/12.99. The value 1/12.9 is the slope near D ≈ 60 (k_61/61 = 12.90). | Say "→ 1/13 asymptotically". |
| F11 | §7.2 Reduction B, lines 555–557 | presentation | Prop 7.2 shows that "T_D splits at k" implies the free extension from k. The converse also needs V_D(k) ⊂ T_D: given that, the free-extension law gives V_D(k′) ⊂ T_D and then the splitting. E3 supplies it only at D = 4, k = 7. The sentence "T_D splits at 4D − 9 **is equivalent to** the measured free-extension law" omits this. | Add "together with V_D(4D−9) ⊂ T_D". |
| F12 | T3 row (line 15) | presentation | Prop 2.1 proves FEI_k(D−1) ⇔ "S_D(k) closed" for every instance (I checked the proof). The "exact on the computed cases" qualifier undersells it. | Label the equivalence **proof**. |
| F13 | script docstrings | presentation | The cross-references have drifted. `lemmas.py` calls Lemma 4.4 "Lemma 4.3" and Theorem 3.3 "Theorem 3.2". `charp.py` says "Theorem 3.2". `lowpart.py` calls the Prop 7.1 criterion "Proposition 7.2". | Renumber. |
| F14 | process | procedural | No commit was named in the dispatch. I froze blob `6c98619e…` myself (§0). The file sits in a WIP snapshot that is "not audited". | Re-review if the blob changes. |

## 3. Critical steps re-derived by hand (verified, not just read)

* **Theorem 3.1 (a), the Wronskian criterion.** The induction holds in any characteristic. ker(d/dz) = F^p, because
  ker ⊇ F^p, ker ≠ F, and [F : F^p] = p is prime. An F^p-independent family has at most p members, so the rows
  δ^0..δ^(r−1) never reach δ^p = 0.
* **Theorem 3.1 (b), change of derivation.** δ^i = c^i δ_t^i + (lower terms), so
  W_δ = c^(r(r−1)/2) W_{δ_t} with ord_P c = −ord_P(dz).
* **Theorem 3.1 (c), local bounds.** (3.2) is the termwise estimate. (3.3) follows from
  δ_t(t^(pe) y) = t^(pe) δ_t(y) in characteristic p. A uniformizer at a smooth point is separating, and d/dt is the
  formal derivative on k̄((t)).
* **Theorem 3.1 (d), conclusion.** W = ±W(v_j : j ≠ i) for every i, by multilinearity. Choosing i_P with minimal order
  and summing gives 0 ≥ −n₀(2γ−2) + H − n₀\|𝒮\| with n₀ = (m−1)(m−2)/2. This holds with Hasse or ordinary
  derivatives alike, since only orders 0..m−2 < p occur. **The constant is right.**
* **Lemma 3.2, orders.**
  - ξ₁ = x(P) has a double pole at O and zeros at ±P₀ (a double zero if b = 0).
  - ξ₂ = x(R − P) has a double pole at P = R and zeros at R ∓ P₀.
  - ξ₂(O) = t and ξ₁(R) = t, so **t ≠ 0 is exactly what makes ord_O ξ₂ = ord_R ξ₁ = 0**. It also keeps R ∓ P₀ ≠ O.
  - \|𝒮\| ≤ 6, γ = 1, 2γ − 2 = 0, so H ≤ 3(m−1)(m−2).
  - (3.1) gives 2H ≥ \|ord_O\| + \|ord_R\|, and these orders are ≡ −2Δa and −2Δb mod p because φ_e ∈ F^p.
  - **Theorem 3.3 follows.**
* **Where p > 2^(k+1) is used.** Not in Theorems 3.1 and 3.3. It is used only in Lemmas 4.2–4.4 and Theorem 5.1:
  with p ≥ 2N + 1 the mod-p closeness of 2(a − a′) can fold only onto the two ends of [0, N). For N < p < 2N the fold
  lands inside the interval and Lemma 4.4 fails. Whether T1 itself fails for such p is open.
* **Lemma 3.4.** Correct: it is a circuit-exchange argument over F₀, and K ⊂ k̄ ⊂ F^p.
* **Lemmas 4.2, 4.3, 4.4 and Theorem 5.1.** Correct as written. I checked the three-interval argument, the count
  m ≤ 36·C(2ℓ, ≤w), the chain 2σ + 1 ≤ 6m² ≤ 7776·C², the two digit cases, and the conditions λ ≤ k and λ − 1 ≤ j.
* **Propositions 2.1 and 2.2.**
  - The lift identity (1.1) is correct, with no collisions because monomial multiplication is a bijection on
    exponents.
  - Primality of 𝔭 together with x_i, x_i^N − 1 ∉ 𝔭 gives the equivalence.
  - g̃ = S_3 for k ≥ 2.
  - The model matches the kernel definitions. `InQuotClosure` uses `w.totalDegree < D` on reduced (multilinear)
    polynomials, and that total degree is the binary weight. `PrimeM2General.normalized` is S_3 with x_i² ↦ u_{i,1}
    (`raw_normalization_identity`), weight 2. The chain tail u_{i,k−1}² − 1 gives x_i^N = 1.
* **Corollary 5.2.**
  - p ≡ 1 mod 2^k and p ≠ 2^k + 1 give p ≥ 2^(k+1) + 1.
  - The step 1 ∉ V_D ⇒ ¬HasLFDBound e for every e ≤ D is `prime_m2_lfd_gt`.
  - The bound LFD ≥ 1 + max{D : k_D ≤ k} ≥ ⌊(k−10)/14⌋ + 1 holds for k ≥ 38, and trivially below that.
  - "LFD ≥ 8 from k = 90" checks out (k_7 = 90).
* **Corollary 5.3.** Lemma 4.3 at depth k + 1 needs p > 2^(k+2), which is assumed. Both digit cases need
  k ≥ D + Λ(D). The dimension recursion follows weight by weight.
* **Proposition 7.1.** Correct (closedness of S^♯, and FEI ⇒ Crit).
* **Proposition 7.2.** Correct:
  - (b) ρ_k does not raise weights and sends x_i^(2^j) to u_{i,j} or to 1.
  - (c) The slice argument works because weight is additive across the cut. The low part is closed in A_k and is
    sandwiched by (b).

## 4. Literature

* **Brownawell–Masser 1986** (Math. Proc. Camb. Phil. Soc. 100, 427–434). The statement as quoted in
  [arXiv:2201.10964](https://arxiv.org/abs/2201.10964), Theorem 5: F/ℂ of genus g; S-units with 1 + u₁ + … + u_k = 0
  and no vanishing proper subsum; then max_i H(u_i) ≤ C(k,2)·(\|S\| + max(0, 2g−2)). Comparison: F2.
* **Voloch 1985** (Bol. Soc. Brasil. Mat. 16), and the Wronskian method in characteristic p (J. T.-Y. Wang). I found
  only secondary descriptions; the exact statement and constants were **not verified**.
* Neither is needed for correctness, because N8 writes the proof out in full and I checked it.

## 5. Computations

All runs used `ulimit -v 1500000` and one thread. No heavy lock was needed and no Lean was built.

### 5.1 Reruns of the N8 scripts (results identical to the records)

* `ic_bounds_growth_lemmas.py lambda | window | straddle | semaev`: ok. 60,920 window cases (at most 2 blocks);
  74,856 straddle cases; 2000 Semaev cases.
* `ic_bounds_growth_charp.py`: 266 + 268 + 269 sets, none dependent; control rank 8 of 9.
* `ic_bounds_growth_ideal.py 12289 4021 9917 4408 3 2 6` and `… 4 6 7`: every dimension and verdict as in §6.2.
* `ic_bounds_growth_circuits.py … 3 3 5`, `ic_bounds_growth_lowpart.py … 4 6 7`, `ic_bounds_growth_nowrap.py … 3 3 6`:
  as recorded.

**What these scripts test against what the text says:**
* `lambda`, `window`, `straddle` and `semaev` test exactly the finite statements claimed.
* `sparse` and `charp` test only a weak consequence of Theorem 3.3: independence over the constants.
* `circuits` tests the characteristic-0 inequality on circuits over the constants, which need not be F^p-circuits.
  Its docstring admits this; the §6.1 row title does not.

### 5.2 Independent recomputation (own code, own elimination, two disjoint point sets of size > 4(2^k − 1))

* Script: `/tmp/claude-0/-home-user/a6e79f2a-a5ec-5904-b377-343a14909fff/scratchpad/n8rt/my/indep_S.py`.
* dim S_3(k) = 3, 16, 41, 74, 111 for k = 2..6, and dim S_4(k) = 564, 1025, 1657 for k = 6, 7, 8. Both match the
  records.
* Closedness was tested from the definition (multiply in A_k, then test the lift on the second point set), separately
  from FEI. Both are false at D = 3, k = 2, 3 and true at k = 4..6. At D = 4 both are false at k = 6, 7 (12/12 failing
  rows at k = 7) and true at k = 8.
* Script: `/tmp/claude-0/-home-user/a6e79f2a-a5ec-5904-b377-343a14909fff/scratchpad/n8rt/my/lam.py`.
* Λ(D−1) = 25, 36, 48, 59, 71, 83, 95, 106, 118, 142, 189, 237 and k_D = 27, 39, 52, 64, 77, 90, 103, 115, 128, 154,
  205, 257 match the table.
* For w ≤ 60: 4w ≤ Λ(w) ≤ 13w + 23, and max over D ≤ 61 of k_D − (14D + 10) is −11.
* The proof's differences are 7.77, 7.73, 8.07, 8.55, 9.11, 9.71, 10.34, and 2·log₂(2e·16) = 12.885.

### 5.3 Falsification test of Theorem 3.3 over F^p (new; regime not tried by the author)

**Criterion.** For f₁..f_m ∈ F with m ≤ p, the family is F^p-dependent if and only if the Hasse Wronskian with respect
to x₁ (orders 0..m−1) vanishes identically. This is part (a) of Theorem 3.1, using D^i = i!·D^(i). The Wronskian is
evaluated at random points of C(F_{p²}) through the local expansion x₁ = α + s, x₂ = y(s), checked to satisfy
S_3 = 0 mod s^m. A nonzero value at one point **proves** F^p-independence.

**Exhaustiveness.** F^p-dependence of monomials depends only on the exponents mod p, because
ξ^(e+pu) = (ξ^u)^p·ξ^e. So scanning all p² residues covers every monomial of every degree.

Scripts in `/tmp/claude-0/-home-user/a6e79f2a-a5ec-5904-b377-343a14909fff/scratchpad/n8rt/my/`: `wronski.py`
(library), `wdrive.py` (run log `wdrive.log`), `w3fix.py`, `w8.py`, `wcheck.py`.

**Controls (p = 131, y² = x³ + 5x + 9, t = 2):**
* The 9 exponents of S_3 give W ≡ 0 at 12 points, and every 8-subset gives W ≠ 0 somewhere.
* Scanning all 17,161 residues with 8 of the S_3 terms recovers exactly the removed term (3 of 3).
* Twisted S_3 with exponents e + p·u_e, u_e ∈ [0,3]² (5 trials): W ≡ 0, while the rank over the constants is 9.
  **So the dependence uses non-constant F^p coefficients, which the N8 tests cannot detect.**
* A pair e, e + p·u gives W ≡ 0.

**Non-vacuity:**
* For 10 random sets B, the survivors of the scan are exactly B.
* Scan values equal direct 5×5 determinants in 50 of 50 cases.
* For the 3-circuit code, the cross-product test agrees with the direct 3×3 Wronskian.

**Results.** No violation; in fact no F^p-dependence at all outside the trivial ones:
* **Exhaustive 3-term F^p-circuits** (WLOG containing exponent 0, by translation invariance), with wildcard points
  handled: none for p = 61, 131 and 251 (t ∈ F_p), and none for p = 131 with a, b, t ∈ F_{p²} ∖ F_p. Rigorous: every
  pair was eliminated by a point where W ≠ 0.
* **Pairwise-far sets B**, distance ‖2Δa‖_p + ‖2Δb‖_p > 6(m−1)(m−2), so Theorem 3.3 predicts that B ∪ {f} is
  independent for every f ∉ B. Each B was scanned against all p² residues:

  | p | m | threshold | sets B | note |
  |---|---|---|---|---|
  | 131 | 5 | 72 | 150 | |
  | 131 | 5 | 72 | 100 | a, b, t ∈ F_{p²} ∖ F_p |
  | 311 | 6 | 120 | 40 | |
  | 401 | 7 | 180 | 25 | |

  **No F^p-dependent extension.**
* **p = 521, m = 8** (threshold 252 < p − 1): 20 random independent 7-sets, each scanned against all 271,441
  residues. No dependent extension.

**Limits of this test:**
* For m ≥ 5 the sets B are random, so this is a search for counterexamples, not a proof.
* The positive controls show that the method does detect both constant and non-constant F^p-dependences.

## 6. What I could not verify

* **Not re-run by me** (I read the records):
  - the D = 5 certificate FEI_12(4) on p = 40961 (heavy lock, 665 s) and FEI_k(4) for k = 8..11;
  - `ideal.py` for D = 3, k = 7..9 and D = 4, k = 9;
  - the unit and t = 0 instances;
  - `nowrap.py` for D = 4, k = 7, L = 10..12.
* I did check that the D = 5, k = 12 record uses disjoint kernel and check point sets, with 16,445 check points
  (more than 4(2^12 − 1)), and that its float64 arithmetic stays below 2^53.
* Voloch's characteristic-p statement and constants (§4).
* Lean side: nobody has formalized the identification of S_D(k) with a submodule `W` satisfying the hypotheses of
  `prime_m2_lfd_gt`. I checked only that the bridge's hypotheses are the ones the N8 argument supplies.

## 7. Labels after this review, and what would raise them

* **T1: proof** (no change). To raise it, kernel-formalize Props 2.1/2.2, Lemma 3.4, Lemmas 4.1–4.4 and
  Theorem 5.1, with Theorem 3.3 as an explicit hypothesis on 𝔭 (the N8 §8 plan). That would make it
  "kernel, conditional on Theorem 3.3 (proof)". Theorem 3.1 itself needs function-field Wronskian theory, which the
  pinned Mathlib lacks.
* **C1: proof** (given T1; the bridge is kernel). Fix the hypothesis statement (F4). To raise it, instantiate
  `ChainClosureBridge.prime_m2_lfd_gt` with W = S_D(k). This needs the same conditional formalization as T1.
* **T4: proof.** To raise it, use the same infrastructure as T1.
* **Prop 7.1: proof.** A small Lean file over an abstract prime ideal would raise it to kernel.
* **Prop 7.2: proof.** It lives in K[x₁, x₂] with digit-sum weights and needs no function fields, so kernel is
  feasible.
