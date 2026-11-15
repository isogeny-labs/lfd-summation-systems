# Red-team review of `docs/proofs/N8-growth.md` (Theorem 5.1 / Corollary 5.2)

Reviewed artifact: `/home/user/ml-cryptanalysis/docs/proofs/N8-growth.md` (603 lines, status line
2026-11-12), with `docs/proofs/U1-upper-bounds.md` §1, `experiments/ic-bounds-k14/RESULTS.md`,
`scripts/ic_bounds_growth_*.py` and `experiments/ic-bounds-growth/*.json` as context. Nothing was
edited; the only file written is this one. Verdict vocabulary: **holds** (argument correct as
written), **gap** (true but not proved by the text as written, or stated imprecisely), **wrong**
(a stated step is false), **unsupported** (claim without a frozen artifact), **unpriced** (cost or
hypothesis transferred without being stated). Line numbers refer to the draft.

## 1. Findings, ordered by severity

| # | item | verdict | draft line | reason (one sentence) |
|---|---|---|---|---|
| F1 | Section 6.2, D = 4 table on the main instance (k = 4..9), and the "33 computed cases" of Prop. 2.1 | **unsupported** (artifact missing) | 374, 401–408 | The table and the count cite `ideal-D4-*.json` for "D = 4: k = 4, …, 9 on the main instance" (p = 12289, a = 4021, b = 9917, t = 4408), but `experiments/ic-bounds-growth/` contains only `ideal-D4-p12289-unit-s10.json` (unit instance, k ≤ 8); no file in the directory records dim S_4(9) = 2403, "V_4(9) ⊂ S_4(9)" or "S_4(9) closed" on the main instance (the number 2403 occurs only in `fei-w4-k9-p40961.json`, a third instance), so six of the 33 "exact" cases and the k = 9 row have no frozen evidence; the unit file does reproduce rows k = 4..8 of the table (105/270/564/1025/1657; k = 7: `{'0,6': 12, '1,6': 12}` = "12 rows fail under u_{i,6}"). |
| F2 | Lemma 4.1, two-step ratio bound | **gap** (derivation wrong, inequality true, conclusion survives) | 230–231 | The text derives "C(2ℓ+2, ≤w)²/C(2ℓ, ≤w)² ≤ ((2ℓ+2)/(2ℓ+2−w))⁴" from binom(n+1,i)/binom(n,i) ≤ (n+1)/(n+1−w), but two applications of that one-step bound give the factor ((2ℓ+2)/(2ℓ+2−w))²·((2ℓ+1)/(2ℓ+1−w))², and (2ℓ+1)/(2ℓ+1−w) > (2ℓ+2)/(2ℓ+2−w), so the displayed exponent-4 bound with numerator 2ℓ+2 does not follow; the displayed inequality happens to be true numerically (no violation for 1 ≤ w ≤ 40, 4w ≤ ℓ < 14w+40, exact rationals) and the corrected bound ((2ℓ+1)/(2ℓ+1−w))⁴ < (8/7)⁴ = 1.706 < 2 for ℓ ≥ 4w (w/(2ℓ+1−w) ≤ w/(7w+1) < 1/7) gives the same monotonicity, with the exact minimum ratio ρ(ℓ+1)/ρ(ℓ) = 1.178 > 2·(7/8)⁴ = 1.172 as claimed. |
| F3 | C1 / Corollary 5.2: "the last fall degree is unbounded" | **gap** (statement needs a quantifier) | 13, 298–304, 311–313 | For a fixed prime the hypotheses 2^k ∣ p−1 and p > 2^{k+1} bound the depth by k ≤ v₂(p−1) and k ≤ log₂ p − 1, so over one prime only finitely many depths exist; what is proved is "for every D there is k_D such that every unit instance at every depth k ≥ k_D over every admissible prime (or over any field of characteristic 0) has LFD > D", a uniform lower bound over a family of primes, and C1 should say so rather than "unbounded in the depth". |
| F4 | Corollary 5.2 hypothesis "2^k ∣ p − 1" and the remark that the bridge is "stated in the repository for 2^k ∣ p − 1" | **gap** (misattributed hypothesis; corollary weaker than its own proof) | 298, 475–476 | The kernel bridge `chain_mem_closure_iff` (`formal/ICFallDegrees/ChainClosureBridge.lean:404`) is stated for an arbitrary `[Field K]` with the sole hypotheses `(2 : K) ≠ 0` and `2 ≤ D`; the draft uses 2^k ∣ p−1 only to turn "p ≠ 2^k+1" into p ≥ 2^{k+1}+1 (line 307, correct), so the corollary holds for every prime p > 2^{k+1} with unit ideal, and the "needs the kernel bridge … stated for 2^k ∣ p−1" sentence is inaccurate (I checked the Lean signature, not its semantics). |
| F5 | Theorem 3.1: novelty claim "the characteristic-p form with coefficients in F^p is proved here" | **holds** as mathematics; **gap** as attribution | 14, 101–127, 454 | The proof is correct (checked step by step below), but the statement with hypothesis 3 strengthened to "v_i are 𝒮-units" and hypothesis 2 "linearly independent over K^p" is the known characteristic-p n-term abc theorem: three terms Voloch, *Diagonal equations over function fields*, Bol. Soc. Brasil. Mat. 16 (1985) 29–39; the Wronskian criterion "W ≠ 0 ⇔ independent over K^p" Garcia–Voloch, *Wronskians and linear independence in fields of prime characteristic*, Manuscripta Math. 59 (1987) 457–469 (Thm 1, as cited by later papers); n terms J. T.-Y. Wang, *A note on Wronskians and the ABC theorem in function fields of prime characteristic*, Manuscripta Math. 98 (1999) 255–264, and Hsia–Wang, *The ABC theorem for higher-dimensional function fields*, Trans. AMS 356 (2004) 2871–2887; the only extra point in the draft is that hypothesis 3 asks for p ∣ ord_P(v_i) instead of ord_P(v_i) = 0 outside 𝒮, which is a one-line generalization of the same proof (its (3.3)), so the theorem should be cited as "Wang/Hsia–Wang, with hypothesis 3 relaxed" (I could not fetch the papers' texts through the proxy; theorem numbers unverified, bibliographic data verified). |
| F6 | Theorem 3.1 proof: ker δ "a field between F^p and F, different from F" | **holds** (make explicit) | 125–127 | ker δ is a subfield containing F^p and not equal to F, and [F : F^p] = p is prime, so ker δ = F^p exactly (no intermediate field); the draft states the three facts but never draws the conclusion, which is what makes hypothesis 2 (independence over F₀ = F^p) exactly the hypothesis of criterion (a); for p = 0 the constants of k̄(X) are k̄ = F₀. Note also hypothesis 2 forces m − 1 ≤ [F : F^p] = p, automatically satisfied for a circuit. |
| F7 | Theorem 3.1 (a)–(d): Wronskian criterion, change of derivation, (3.2), (3.3), the sum over places | **holds** | 129–166 | (a) is the standard induction and uses nothing about the characteristic; (b) δ = (dt/dz)·δ_t with ord_P(dt/dz) = −ord_P(dz) gives W_δ = c^{r(r−1)/2} W_{δ_t}; (c) δ_t is the formal derivative in k̄((t)) because t_P ∉ F^p is separating and derivations extend uniquely through separable extensions, (3.2) is the determinant expansion, and (3.3) is δ_t(t^{pe}) = 0 for every integer e (negative allowed) with ord δ_t^l(y) ≥ 0 for a unit y; (d) excludes the function of minimal order at each P (any choice gives ±W by Σ v_i = 0), applies (3.3) at P ∉ 𝒮 to the remaining m−1 functions, all of p-divisible order by hypothesis 3, and sums Σ ord_P(W) = 0, Σ ord_P(dz) = 2γ−2 (dz ≠ 0 as z is separating), Σ_P A_P = 0, −Σ μ_P = H(v): H ≤ n₀(2γ − 2 + \|𝒮\|). |
| F8 | Lemma 3.2 / Theorem 3.3: divisors, the set 𝒮, orders modulo p | **holds** | 172–197 | div(x) = (P₀)+(−P₀)−2(O) (with multiplicity 2(P₀) when b = 0) and ξ₂ = x∘τ, τ(P) = R−P an automorphism, so div ξ₂ = τ*div x = (R−P₀)+(R+P₀)−2(R); the only facts used are "ξ₂ is a unit at O" ⇔ O ∉ {R, R±P₀} ⇔ R ∉ {O, ±P₀} ⇔ t ≠ 0, and symmetrically "ξ₁ is a unit at R" ⇔ t ≠ 0; R being 2-torsion or R = ±2P₀ only merges points of 𝒮 (\|𝒮\| ≤ 6 still), and S₃ reducible (S₃ = c·h² when 2R = O) is irrelevant because everything is phrased through 𝔭; with φ_e = ψ_e^p, ord(φ_e) = p·ord(ψ_e), so ord_O(v_e/v_{e′}) ≡ −2(a−a′) and ord_R ≡ −2(b−b′) mod p, and (3.1) gives ‖2Δa‖_p + ‖2Δb‖_p ≤ 2H ≤ 6(m−1)(m−2). |
| F9 | Lemma 3.4: matroid argument, parallel elements, K ⊂ F₀ | **holds** | 199–206 | K ⊂ k̄ ⊂ F₀ because k̄ is perfect (every constant is a p-th power) so c_e ∈ F₀ and y lies in both F₀-spans; bases B_i ⊂ M_i exist, B₁ ∪ B₂ is dependent when y ≠ 0, and a minimal dependent subset is a circuit of the matroid on the ground set M (parallel elements, i.e. e ≠ e′ with ξ^e/ξ^{e′} ∈ F₀, are 2-circuits and are handled by exactly the same sentence); the hypothesis is about circuits J ⊂ M, which is what the proof uses. |
| F10 | Lemmas 4.2, 4.3, 4.4 | **holds** | 242–274 | 4.2: \|2(a−a₀)\| ≤ 2(N−1) < p so the distance condition is \|a−a₀\| ≤ σ or ≥ (p−2σ)/2 ≥ N+½−σ, the far set lies in [0, σ) ∪ (N−σ, N), each of the three intervals has ≤ 2σ+1 ≤ 2^ℓ elements and meets ≤ 2 aligned blocks (6 is a valid, loose bound; brute force finds 2); 4.3: m ≤ 6²·C(2ℓ, ≤w) because low ℓ bits determine a pair inside a block pair, 2^{ℓ−1} < 2σ+1 ≤ 6m² ≤ 7776·C² is the defining inequality of Λ so ℓ ≤ Λ(w) and σ < 2^{Λ(w)−1}, ℓ = 0 gives σ = 0, and ℓ > k only makes the count an overcount; 4.4: a′ = (N−2^{j+1}) + (2^j − e) with 1 ≤ e < 2^{λ−1} ≤ 2^j puts k−j−1 ones at bits j+1..k−1 and j−λ+1 ones at bits λ−1..j−1, and in the far case a > N − 2^{λ−1} has ones at bits λ−1..k−1 (k−λ+1 of them). |
| F11 | Theorem 5.1: j ≥ k − w, λ − 1 ≤ j, the use of p ≥ 2N+1, the circuit hypothesis | **holds** | 282–296 | A top-run exponent a ≥ N−2^j has bits j..k−1 set, so k−j ≤ wt(a) ≤ w and j ≥ k−w ≥ Λ(w)+1 ≥ λ−1 (two units of slack); the binding constraint is the third one, w < k−λ ⇔ k ≥ w+Λ(w)+1, and λ ≤ k follows from it; p > 2^{k+1} is p ≥ 2N+1, exactly what Lemmas 4.2–4.4 need; Lemma 3.4 needs "no circuit J ⊂ M meets both parts", and the proof excludes every F₀-circuit of pairs in [0,N)² of weight ≤ w containing a straddling pair, a superset of the circuits inside M, so the hypothesis is the right one and is established. |
| F12 | Corollary 5.2 arithmetic: p ≥ 2^{k+1}+1, LFD > D, the floor bound | **holds** | 302, 307–309 | p ≡ 1 (2^k) and p ≠ 2^k+1 give p ≥ 2·2^k+1; 1 ∉ V_D with 1 ∈ I and V_D ⊂ V_{D+1} gives LFD > D under U1 (1.2); k_D ≤ 14D+10 (Lemma 4.1) gives 1+max{D : k_D ≤ k} ≥ ⌊(k−10)/14⌋+1 as soon as ⌊(k−10)/14⌋ ≥ 2, i.e. k ≥ 38 (checked exactly for 38 ≤ k < 560 with the computed Λ; sample: k = 39 → 4 vs 3, k = 90 → 8 vs 6, k = 400 → 32 vs 28); for k < 38 the displayed bound is ≤ 2 and trivial, so the unqualified form in C1 is harmless. |
| F13 | Corollary 5.3 and Proposition 7.1 | **holds** | 315–330, 484–494 | 5.3 is Theorem 5.1's argument at depth k+1 for elements of weight ≤ D (hence λ = Λ(D), p ≥ 2·2^{k+1}+1 = hypothesis p > 2^{k+2}), with the near case a′ ∈ (N−2^{λ−1}, N) of weight ≥ k−λ+1 > D and the far case a > 2N−2^{λ−1} of weight ≥ k−λ+2, then the same for the second exponent and division by x^{Nq} (𝔭 prime, x_i ∉ 𝔭); 7.1 needs g ∈ S_{D−1}(k), i.e. D ≥ 3 as stated, and the implication FEI ⇒ Crit uses S_{D−1}(k) = S_D(k) ∩ M_{≤D−1}, which is immediate from the definition. |
| F14 | Hidden assumptions: model identification, weight = chain degree, D ≥ 2 edge cases | **holds** (with F4's caveat) | 37–50, 70–71 | The definition of V_D(k) at lines 46–47 is verbatim U1 §1's "quotient closure"; g is the image of NF(S₃) = image of S₃ in A_k (reduction modulo the chain relations is the quotient map), of weight 2 for k ≥ 2; weight = degree of the multilinear normal form, and carries are exactly the wrap of a top run, which (1.1) tracks correctly; the Lean bridge has `[Field K]`, `2 ≠ 0`, `2 ≤ D` (signature checked, semantics not); D = 2 is vacuous (V_₂ = span g, so k_₂ = 27 is a true but empty row) and Prop. 7.1 correctly requires D ≥ 3; the unit-ideal claim for the generic system over K(a,b,t) (line 304) is asserted without the one-line argument (no (α, β) ∈ μ_N² solves S₃(α, β, t) = 0: the t²-coefficient is (α−β)², and for α = β the coefficient of a·t is −4α ≠ 0). |
| F15 | Headline trivially true or false? Cheap attacks | **holds** (not trivial either way) | 12–13, 311–313 | (i) Wrap-around construction at large k, small D: for D = 3, V₃(k) = span{g, u_{i,j} g} (dim 2k+1 on the data; the only possible collapse is u_{1,0}g and u_{1,1}g sharing their weight-3 support (3,2),(3,1), which cancels only when a = −3t²), so no wrap ever produces 1 and the theorem's k₃ = 39 is far from the truth but not false; for D = 4, 5, 6 the certified exclusions 1 ∉ V_D at k = 6, 9/10, 14 plus the kernel depth-restriction theorem (U1 P3: 1 ∈ V_D(k+1) ⇒ 1 ∈ V_D(k)) already give 1 ∉ V_D(k) for all larger k on those instances, so Theorem 5.1 adds nothing testable below D = 7 (k ≥ 90) and cannot be contradicted by any present data; (ii) consistency: 1 ∈ V₇(14) is compatible because depth restriction only propagates exclusions upward and k₇ = 90 > 14; 1 ∉ V₆(14) vs k₆ = 77 is consistent; the lower bound ⌊(k−10)/14⌋+1 is ≤ k+2 = U1's upper bound for all k, and k_D ≥ 5D−4 > 4D−10 so the theorem never contradicts the measured base depth. Simpler proofs: none found; the statement is not implied by counting (unit ideal) or by depth restriction (which needs a seed at each D). |
| F16 | Section 6.1 row "Theorem 3.3 on computed circuits" | **gap** (labelling) | 376 | Reduced-echelon kernel vectors are circuits over the constants (pivot columns are independent, so each null vector's support is minimal), but Theorem 3.3 in characteristic p is about F^p-circuits; a constant-circuit of large spread would not refute the theorem unless it is also an F^p-circuit, so this row is a consistency check, not a test of the statement, and "spread ≤ m/2" is an observation about constant-circuits only. |
| F17 | Slope "about 1/12.9 with the exact Λ" | **gap** (minor numerics) | 13 | Asymptotically Λ(w)/w → c with H(1/(2c)) = 1/4 (binary entropy), c ≈ 12.05, so k_D/D → ≈ 13.05 and the exact slope tends to ≈ 1/13.05; the table gives (Λ(120)−Λ(100))/20 = 12; not load-bearing, the proved 1/14 stands. |
| F18 | Section 6 "exact" method for S_D(k) | **holds** (one unstated case) | 356–359 | "A polynomial of bidegree < 2^k vanishing on more than 4(2^k−1) affine F_p-points of C lies in 𝔭" is Bézout on P¹×P¹ against the (2,2)-curve V(𝔭) and is fine; when 2R = O (t a root of x³+ax+b) S₃ = c·h² with V(𝔭) = V(h) of bidegree (1,1), the point set is the same and the bound only improves, so the method is still exact; the docstring's "irreducible (2,2)-curve" should say "the curve V(𝔭) ⊂ V(S₃), of bidegree ≤ (2,2)". |

Finite checks named by the draft (item 10 of the dispatch), run from `/home/user/ml-cryptanalysis`
under `ulimit -v 1500000`, outputs as seen:

```
python3 scripts/ic_bounds_growth_lemmas.py lambda     (66.6 s)
  {"test": "lambda", "ok": true, "rows": [{"w": 1, "D": 2, "Lambda": 25, "k_D": 27}, {"w": 2, "D": 3, "Lambda": 36, "k_D": 39},
   {"w": 3, "D": 4, "Lambda": 48, "k_D": 52}, {"w": 4, "D": 5, "Lambda": 59, "k_D": 64}, {"w": 5, "D": 6, "Lambda": 71, "k_D": 77},
   {"w": 6, "D": 7, "Lambda": 83, "k_D": 90}, ..., {"w": 19, "D": 20, "Lambda": 237, "k_D": 257}, ..., {"w": 120, "D": 121, "Lambda": 1445, "k_D": 1566}]}
python3 scripts/ic_bounds_growth_lemmas.py window     (1.2 s)
  {"test": "window", "ok": true, "cases": 60920, "max_blocks": 2}
python3 scripts/ic_bounds_growth_lemmas.py straddle   (0.2 s)
  {"test": "straddle", "ok": true, "cases": 74856}
python3 scripts/ic_bounds_growth_lemmas.py semaev     (0.1 s)
  {"test": "semaev", "ok": true, "cases": 2000}
```

All four agree with the draft (lines 227–228, 240, 370–373) and with the frozen `lemmas.json`,
`semaev.json`. Independent direct computation (own script, exact integers): Λ(1) = 25 (the
inequality 2^24 < 7776·51² holds, 2^25 < 7776·53² fails), Λ(2) = 36 (2^35 < 7776·2629² holds,
2^36 < 7776·2776² fails), Λ(1..8) = 25, 36, 48, 59, 71, 83, 95, 106; the set of ℓ satisfying the
inequality is an initial segment of ℓ ≥ 1 for w ≤ 9 (not only on [4w, ∞)); ρ(13w+24) ≥ 1 for
1 ≤ w ≤ 200; the seven evaluated differences in the proof of Lemma 4.1 are 7.77, 7.72, 8.06,
8.54, 9.10, 9.71, 10.34 (draft: 7.7, 7.7, 8.0, 8.5, 9.1, 9.7, 10.3) and 2·log₂(5.437·16) = 12.886
< 13. Also rerun (allowed range D = 3, k ≤ 9, 47 s):
`python3 scripts/ic_bounds_growth_ideal.py 12289 4021 9917 4408 3 2 9` reproduces every entry of the
D = 3 table (dim V₃ = 5, 7, 9, …, 19; dim S₃ = 3, 16, 41, 74, 111, 126, 130, 134; V₃ ⊂ S₃ from k = 3;
S₃ closed ⇔ FEI_k(2) from k = 4, equal verdicts at all eight depths; 1 ∉ V₃ throughout).

## 2. Cheapest falsification not carried out

Theorem 5.1 itself cannot be falsified by computation for D ≤ 6 (its conclusion there is already
a kernel-checked consequence of the certified exclusions plus depth restriction, F15) and is out
of reach for D ≥ 7 (k ≥ 90). What the proof actually establishes is the stronger, non-monotone
statement FEI_k(w), which is testable in the gap between the data (k ≤ 9 for w = 2, 3; k ≤ 12 for
w = 4) and the theorem (k ≥ 39, 52, 64):

```
python3 scripts/ic_bounds_growth_fei.py 40961 21223 9887 4748 2 10   # 221 columns, 4092 points needed
python3 scripts/ic_bounds_growth_fei.py 40961 21223 9887 4748 2 11   # 265 columns, 8188 points
python3 scripts/ic_bounds_growth_fei.py 40961 21223 9887 4748 2 12   # 313 columns, 16380 of 40815 points
python3 scripts/ic_bounds_growth_fei.py 40961 21223 9887 4748 3 10   # 1351 columns
python3 scripts/ic_bounds_growth_fei.py 40961 21223 9887 4748 3 11   # 1794 columns
```

Each is a kernel computation of seconds to a few minutes (the recorded w = 4, k = 11 case took
244 s with 9109 columns); they were not run because the dispatch excludes k ≥ 10. A failure of
FEI_k(2) or FEI_k(3) at any k ≥ 10 would refute the conjectured pattern "FEI_k(w) ⇔ k ≥ 4w − 4"
(E1) and show that the splitting property is not monotone in k; it would not touch Theorem 5.1,
whose range starts at 39. The second cheap test, named but not done by the draft (line 460):
F^p-circuits with m = 5 terms, threshold 6·4·3 = 72 < p/2, need p ≥ 149; `ic_bounds_growth_charp.py`
enumerates E(F_{p²}) in pure Python and p = 151 or 157 (≈ 23k–25k x-values) is well within a
few minutes, so "not testable" at line 460–461 is a resource estimate that does not hold; a
dependent set of 5 pairwise-far monomials would refute Theorem 3.3 outright. The missing artifact
of F1 is settled by `python3 scripts/ic_bounds_growth_ideal.py 12289 4021 9917 4408 4 4 9
experiments/ic-bounds-growth/ideal-D4-p12289-s12.json` (D = 4, excluded here by the same rule).

## 3. Overall verdict

The written proof of Theorem 5.1 stands: every step from Proposition 2.1 through Lemma 4.4 was
checked and holds, including the one step the draft presents as new (Theorem 3.1 in
characteristic p with coefficients in F^p), whose Wronskian argument is correct and is in
substance the known characteristic-p abc theorem of Voloch / Garcia–Voloch / Wang / Hsia–Wang
with the 𝒮-unit hypothesis relaxed to "orders divisible by p outside 𝒮" (F5). The only
defect inside a proof is the mis-derived two-step ratio bound in Lemma 4.1 (F2), which is a
wrong justification of a true inequality and is repaired by replacing the numerator 2ℓ+2 by
2ℓ+1; the conclusions Λ(w) ≥ 4w, Λ(w) ≤ 13w+23 and the table are confirmed exactly. Before the
theorem is cited, four things must change: (1) fix Lemma 4.1 line 231 and state ker δ = F^p
explicitly (F2, F6); (2) cite Theorem 3.1 as a known theorem with the relaxation named, and drop
"is proved here" from T2 (F5); (3) restate C1/Corollary 5.2 as a uniform bound over all
admissible primes p > 2^{k+1} (and characteristic 0), remove the 2^k ∣ p−1 hypothesis or mark it
as inherited from the repository's setting rather than from the bridge, and say that for a fixed
prime the depth is bounded (F3, F4); (4) either commit `ideal-D4-p12289-s12.json` or relabel the
main-instance D = 4 rows and the "33 cases" count of Section 6 as unrecorded (F1). None of these
changes the mathematics; the result is correctly described by the draft itself as having no
content below D = 7 on any computed instance (line 311–313), and its value is the uniform
unboundedness and the slope 1/14, with the true threshold, on all data, near 4D − 8 rather than
D + Λ(D−1). Status and labels are for the Coordinator and Validator to decide; nothing here is a
status change.

Sources checked outside the repository (bibliographic data only; full texts not reachable
through the proxy): Garcia–Voloch, Manuscripta Math. 59 (1987) 457–469
(https://eudml.org/doc/155261); J. T.-Y. Wang, Manuscripta Math. 98 (1999) 255–264
(https://ir.sinica.edu.tw/handle/201000000A/35518); Hsia–Wang, Trans. AMS 356 (2004) 2871–2887
(https://scholar.lib.ntnu.edu.tw/en/publications/the-abc-theorem-for-higher-dimensional-function-fields-2/);
Brownawell–Masser, Math. Proc. Cambridge Philos. Soc. 100 (1986)
(https://www.cambridge.org/core/journals/mathematical-proceedings-of-the-cambridge-philosophical-society/article/vanishing-sums-in-function-fields/34E7F8E89B830D9CBE606324B7E254D9).
