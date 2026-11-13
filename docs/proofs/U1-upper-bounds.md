# U1: upper bounds for the two-summand prime chain LFD below the dimension bound

> **Status note (2026-11-12).** The proof of the `k+2` bound (Sections 1–2.5) was refereed again in `docs/reviews-20261112/referee-U1.md` with no gap found. Findings F1, F2, F5, F6, F7, F9, F10 and F11 of `docs/reviews-U1-upper-bound-redteam-20261110.md` are not resolved in this note; none affects the proof. The directories `experiments/ic-bounds-upper/` and `experiments/ic-bounds-upper-redteam/` named below are not in this repository, and their per-instance records were not kept; the counts are recorded only in the 2026-11-10 review. The paper's proof sketch now follows Section 2.5.

Status (2026-11-10). Evidence labels: **kernel** (Lean-checked), **proof** (written out in full
here), **exact** (certified computation, command and output named), **empirical** (pattern in
data), **conjecture**. Scripts and outputs: `experiments/ic-bounds-upper/`.

## 0. Results

| # | statement | label |
|---|---|---|
| T1 | **LFD ≤ min(k₁,k₂) + 2** for every two-summand prime chain system (depths k₁, k₂ of the two blocks; in the repository k₁ = k₂ = k, so **LFD ≤ k + 2**), every prime p with 2^k ∣ p−1, every nonsingular curve and every target t ∈ F_p; also for coset bases. Replaces the dimension bound 2k+1. | **proof** (Section 2); the derivation of the proof was executed with an enforced degree budget on 16,844 instances, k ≤ 6 (**exact**, Section 2.6) |
| T1' | The constant 2 cannot be lowered in a bound of the form k + c valid for all k: LFD = k+2 is attained at k = 1 and at k = 2. | **exact** (Section 2.7) |
| P3 | Depth restriction. The substitution u_{i,k} ↦ 1 maps degree-D derivations at depth k+1 to degree-D derivations at depth k. Hence 1 ∈ V_D(k+1) ⇒ 1 ∈ V_D(k); for unit ideals LFD(k) ≤ LFD(k+1); an exclusion 1 ∉ V_D(k₀) holds at every deeper depth with no further certificate. | **kernel** (`formal/ICFallDegrees/DepthRestriction.lean`, general m and further inputs; `DepthRestrictionM2.lean`, the normalized two-summand presentation, all k′ ≥ k); model statement **proof** (Section 3) |
| P4 | Filtered generation criterion: LFD ≤ D ⇔ V_D = I_D and I_{D'+1} = I_{D'} + Σ_v v·I_{D'} for all D' ≥ D; by T1 only D' ≤ k+1 has to be checked. | **proof** (Section 4) |
| P5 | Counting obstruction: LFD ≤ D ⇒ dim I_{D'+1} ≤ (2k+1)·dim I_{D'} for all D' ≥ D. | **proof**; its use as a lower bound of order 0.22·k for primes p ≈ 2^k is **conditional** (Section 6.2) |
| E1 | π(V_D(k+1)) = V_D(k) on all 63 tested pairs (D = 3, 4, 5, k ≤ 7, four instances); the filtered version fails at the transition depth; above the base both follow from the free extension. | **empirical** / **proof** above the base (Section 5.2) |
| E2 | For D = 4 and D = 5 the conjectured base depth 4D−10 is the first depth at which the free-extension polynomial f_D(k) fits below \|M_≤D(2k)\|: f_4(5) = 389 > 386, f_5(9) = 12649 > 12616. | **exact** arithmetic on the certified profiles (Section 6.1) |
| C1 | LFD(k) = ⌊(k+14)/4⌋ for generic systems: equivalent, for unit ideals, to 1 ∈ V_D(4D−11) and 1 ∉ V_D(4D−10) for every D ≥ 4. No inductive mechanism was found; the candidates tested are listed in Section 5. | **conjecture** |

What failed or is not claimed: no bound with slope below 1 is proved; no derivation of the
slope 1/4 is given (Section 6 explains what the dimension counts do and do not show); the m ≥ 3
analogue of T1 is only sketched (Remark 2.9).

---

## 1. Setting

Fix a prime p ≥ 5, depths k₁, k₂ ≥ 1 with 2^{k_i} ∣ p−1, N_i = 2^{k_i}, and H_i = μ_{N_i} ⊂ F_p^*.
The chain presentation has variables u_{i,j} (i = 1, 2; 0 ≤ j < k_i), chain equations
u_{i,j}² − u_{i,j+1} (j < k_i − 1), u_{i,k_i−1}² − 1, and one further input ĝ of degree δ.
In the repository ĝ is the normalized generator (the chain normal form of S_3(x₁, x₂, t),
x_i = u_{i,0}), δ = 2, and k₁ = k₂ = k.

**Quotient model.** Modulo the chains the ring is the group algebra

    A = F_p[G],  G = Z/N₁ × Z/N₂,   monomial x₁^{e₁} x₂^{e₂},

with **degree** wt(e₁) + wt(e₂) (binary Hamming weights); u_{i,j} = x_i^{2^j} acts by adding 2^j to
e_i modulo N_i. M_≤D is the span of the monomials of degree ≤ D. Multiplication by a variable
raises the degree by at most one (a carry never increases the Hamming weight by more than the bit
added), so

    wt(a·b) ≤ wt(a) + wt(b)   for all a, b ∈ A                                  (1.1)

where wt of an element is the largest degree in its support. Let g ∈ A be the image of ĝ and
I = gA. As p ∤ N₁N₂ and H_i ⊂ F_p, A is the algebra of functions on H₁ × H₂ and I is the ideal of
functions vanishing on Z = {(α, β) ∈ H₁ × H₂ : g(α, β) = 0}.

**Quotient closure.** For D ≥ δ let V_D ⊂ M_≤D be the smallest subspace containing g and closed
under w ↦ u_{i,j}·w for w ∈ V_D with wt(w) ≤ D − 1. Clearly V_D ⊂ V_{D+1} and V_D ⊂ I ∩ M_≤D =: I_D.

**Bridge (kernel, `ChainClosureBridge.chain_mem_closure_iff`).** For D ≥ 2 and F = chains ∪ {ĝ}, a
polynomial P lies in the ordinary degree-D closure of F iff deg P ≤ D and its chain normal form
lies in V_D (`InQuotClosure`; the identification of reduced polynomials with A is
`ChainTranslation.normalForm_X_mul_monomial`). A polynomial P lies in the ideal ⟨F⟩ iff its normal
form lies in I. (The kernel statement is for m blocks of equal depth k and arbitrary further
inputs; for unequal depths k₁ ≠ k₂ the same proof applies but is not formalized.) Hence, with LFD
as in `Closure.HasLFDBound`,

    LFD(F) ≤ D   ⇔   I_{D'} ⊂ V_{D'} for every D' ≥ D      (D ≥ max(2, δ)).           (1.2)

(⇐: P ∈ ⟨F⟩ of degree d has normal form in I of weight ≤ d, hence in V_{max(D,d)}. ⇒: an element
of I_{D'} is its own normal form, of degree ≤ D'.) Everything below is stated in A.

**Lemma 1.1 (budgeted products).** If w ∈ V_D, a ∈ A and wt(a) + wt(w) ≤ D, then a·w ∈ V_D.

*Proof.* By linearity a is a monomial, a product of r = wt(a) distinct variables v₁ ⋯ v_r. By
(1.1) wt(v_s ⋯ v₁ w) ≤ wt(w) + s ≤ D − 1 for s < r, so the closure rule applies r times. ∎
(Kernel counterpart: `ChainClosureBridge.quot_mul_mem`.)

---

## 2. Theorem 1: the selector bound

**Theorem 1.** Let g ∈ A have x₂-degree at most two,

    g = g₀(x₁) + g₁(x₁)·x₂ + g₂(x₁)·x₂²,     g_e ∈ F_p[Z/N₁],

such that for every α ∈ H₁ the fibre polynomial q_α(X) = g₀(α) + g₁(α)X + g₂(α)X² ∈ F_p[X] is not
the zero polynomial. Let ĝ be an input of degree δ ≥ wt(g) with normal form g. Then

    LFD(chains ∪ {ĝ}) ≤ k₁ + max(2, δ).

The roles of the two blocks may be exchanged when g has x₁-degree at most two and the
corresponding fibre condition holds.

**Corollary 1 (two-summand prime chains).** For every prime p ≥ 5, every k ≥ 1 with 2^k ∣ p−1,
every nonsingular curve y² = x³ + ax + b over F_p and every t ∈ F_p, the normalized presentation
{chains, NF(S_3(x₁, x₂, t))} satisfies **LFD ≤ k + 2**; with depths k₁, k₂ it satisfies
LFD ≤ min(k₁, k₂) + 2. The same holds for coset bases x_i = s_i·u_{i,0} (s_i ∈ F_p^*), and the
raw input S_3 (degree 4, unreduced) gives LFD ≤ min(k₁, k₂) + 4.

The corollary follows from Theorem 1, Lemma 2.2 (fibre condition) and the symmetry of S_3 in
x₁, x₂. There is no genericity, unit-ideal or target-order assumption.

### 2.1 Selectors

For α ∈ H₁ put

    e_α = (1/N₁) Σ_{e=0}^{N₁−1} α^{−e} x₁^e  =  Π_{j<k₁} (1 + α^{−2^j} u_{1,j}) / 2  ∈ A.

**Lemma 2.1.** (i) wt(e_α) ≤ k₁. (ii) x₁^e·e_α = α^e·e_α, hence f·e_α = f(α, x₂)·e_α for every
f ∈ A, in particular u_{1,j}·e_α = α^{2^j} e_α. (iii) e_α e_β = δ_{αβ} e_α and Σ_α e_α = 1.
(iv) A = ⊕_α B·e_α, where B = F_p[Z/N₂] is the block-2 algebra, and b ↦ b·e_α is injective on B.

*Proof.* e_α is the indicator function of {x₁ = α} on H₁ × H₂ (orthogonality of characters of
Z/N₁, valid because p ∤ N₁ and H₁ ⊂ F_p); (ii)–(iv) are immediate. (i): e_α involves block 1
only, which has k₁ variables. ∎

### 2.2 The fibre condition for S_3

**Lemma 2.2.** Let f(x) = x³ + ax + b with 4a³ + 27b² ≠ 0 in F_p, p ≥ 5, and t ∈ F_p. For every
α ∈ F_p the polynomial S_3(α, X, t) ∈ F_p[X] is nonzero, of degree ≤ 2.

*Proof.* S_3(x₁, x₂, t) = (x₁−x₂)²t² − 2((x₁+x₂)(x₁x₂+a) + 2b)t + (x₁x₂−a)² − 4b(x₁+x₂). The
coefficient of X² in S_3(α, X, t) is t² − 2αt + α² = (α − t)², nonzero unless α = t. For α = t the
coefficient of X is −2t³ − 2t(t² + a) − 2at − 4b = −4f(t) and the constant term is
t⁴ − 2at² − 8bt + a² = f′(t)² − 8t·f(t). If both vanish then f(t) = f′(t) = 0, so f has a double
root, contradicting 4a³ + 27b² ≠ 0. ∎

For a coset base the fibre polynomial is S_3(s₁α, s₂X, t), nonzero by the same lemma. The normal
form g of S_3 has x₂-exponents 0, 1, 2 and weight ≤ 2 (x_i² = u_{i,1}, or 1 when k_i = 1).

### 2.3 The univariate lemma (the block-2 fibre has LFD ≤ 2)

Let q ∈ F_p[X] be nonzero of degree ≤ 2, q̄ its image in B = F_p[X]/(X^{N₂} − 1) (so wt(q̄) ≤ 1),
and let U(q) ⊂ B be the smallest subspace containing q̄ and closed under w ↦ u_j·w for
w ∈ U(q) with wt(w) ≤ 1 (u_j = X^{2^j}, 0 ≤ j < k₂). Let d = gcd(q, X^{N₂} − 1), r = deg d ∈ {0,1,2},
R = {1, X, …, X^{r−1}} (so R = ∅, {1} or {1, X}; all of weight ≤ 1), and let π_q : B → span R be
reduction modulo d (B/q̄B = F_p[X]/(d)). The roots of d are the roots of q in H₂, they are simple,
and π_q is interpolation on them.

**Lemma 2.3.** (a) If r = 0 then 1 ∈ U(q). (b) For every j < k₂ and every ρ ∈ R:
u_j·ρ − π_q(u_j·ρ) ∈ U(q).

*Proof.* Write ρ(·) for the image in B, and note ρ(u_j²) = u_{j+1} for j < k₂−1, ρ(u_{k₂−1}²) = 1.

*deg q = 0.* 1 ∈ U(q); r = 0 and (b) is empty.

*deg q = 1,* q = c(X − β). Put w₀ = X − β ∈ U(q) and β_j = β^{2^j}. If w_j = u_j − β_j ∈ U(q)
(weight ≤ 1), then u_j·w_j + β_j·w_j = ρ(u_j²) − β_j² ∈ U(q). So w_j ∈ U(q) for all j < k₂ and
1 − β^{N₂} ∈ U(q). If β^{N₂} ≠ 1 then 1 ∈ U(q) and r = 0. Otherwise β ∈ H₂, d = X − β, R = {1},
π_q(u_j) = β_j, and (b) is w_j ∈ U(q).

*deg q = 2,* q = c(X² − sX + n). For h ∈ F_p[X] write [h] = a + bX for its remainder modulo q.
Let (a_j, b_j) = [X^{2^j}] for 0 ≤ j ≤ k₂, and

    Q₁ = ρ(X²) − sX + n = q̄/c ∈ U(q),
    w_j  = u_j − a_j − b_j X           (0 ≤ j < k₂),      w_{k₂} = 1 − a_{k₂} − b_{k₂} X,
    w′_j = X·u_j − [X·X^{2^j}]         (0 ≤ j < k₂),      w′_{k₂} = X − [X·X^{N₂}],

all of weight ≤ 2, the w_j of weight ≤ 1. We show by induction that all of them lie in U(q).
w₀ = 0. Given w_j ∈ U(q):

    X·w_j + b_j·Q₁ = ρ(X u_j) − a_j X − b_j(sX − n) = w′_j

(for j = k₂ read u_{k₂} as 1), because X·(a_j + b_j X) ≡ a_j X + b_j(sX − n) modulo q; the
product X·w_j is allowed as wt(w_j) ≤ 1. Next, for j < k₂,

    u_j·w_j + a_j·w_j + b_j·w′_j = ρ(u_j²) − a_j(a_j + b_j X) − b_j[X(a_j + b_j X)] = ρ(u_j²) − [(a_j + b_j X)²],

which is w_{j+1}. So w_{k₂} and w′_{k₂} lie in U(q); they are the remainders modulo q of
1 − X^{N₂} and X(1 − X^{N₂}), so they span the image of the ideal (X^{N₂} − 1) in F_p[X]/(q), which
is (d)/(q), of dimension 2 − r.

* r = 0: span{w_{k₂}, w′_{k₂}} = span{1, X} ∋ 1. This is (a).
* r = 1: U(q) contains d = X − β (β ∈ H₂). Then w_j + b_j(X − β) = u_j − (a_j + b_j β) = u_j − β^{2^j}
  (evaluate X^{2^j} ≡ a_j + b_j X at the root β of q). This is (b) with R = {1}.
* r = 2: q divides X^{N₂} − 1, R = {1, X}, π_q is reduction modulo q, and (b) is w_j, w′_j ∈ U(q). ∎

### 2.4 Lifting the fibre derivations

**Lemma 2.4.** Put D = k₁ + max(2, δ). For every α ∈ H₁: U(q_α)·e_α ⊂ V_D.

*Proof.* g·e_α = q̄_α·e_α by Lemma 2.1 (ii), and g·e_α ∈ V_D by Lemma 1.1 since
wt(e_α) + wt(g) ≤ k₁ + δ ≤ D (g enters at budget D because D ≥ δ). If w ∈ U(q_α), wt(w) ≤ 1 and
w·e_α ∈ V_D, then wt(w·e_α) ≤ k₁ + 1 ≤ D − 1, so u_{2,j}·(w e_α) = (u_j w)·e_α ∈ V_D. The set
{w ∈ B : w e_α ∈ V_D} is a subspace containing q̄_α and closed under the rule defining U(q_α). ∎

### 2.5 Proof of Theorem 1

Let D = k₁ + max(2, δ) ≥ k₁ + 2. For α ∈ H₁ let d_α, R_α, π_α be the data of Section 2.3 for q_α
(nonzero by hypothesis). By Lemma 2.1 (iv), I = gA = ⊕_α q̄_α B·e_α, hence

    A = I ⊕ S,   S = ⊕_α span(R_α)·e_α,   wt(S) ≤ k₁ + 1,

and the projector π : A → S along I is π(Σ_α b_α e_α) = Σ_α π_α(b_α) e_α. It satisfies
π(a·π(c)) = π(a·c) (I is an ideal) and I = ker π.

*Corrections.* (C1) 1 − π(1) = Σ_{α : r_α = 0} e_α ∈ V_D by Lemma 2.3 (a) and Lemma 2.4.
(C2) For every variable v and every basis element ρ·e_α of S (ρ ∈ R_α):
v·ρe_α − π(v·ρe_α) ∈ V_D. Indeed for v = u_{1,j} it is zero (Lemma 2.1 (ii)); for v = u_{2,j} it
is (u_jρ − π_α(u_jρ))·e_α ∈ U(q_α)e_α ⊂ V_D by Lemmas 2.3 (b) and 2.4.

*Induction.* Claim: for every monomial m of degree r, m − π(m) ∈ V_{max(D, r)}. For r = 0 this is
(C1). For r ≥ 1 write m = v·m′ with v a variable at a set bit of m, so m′ is a monomial of degree
r − 1. By induction m′ − π(m′) ∈ V_{max(D, r−1)} and its weight is at most
max(r − 1, k₁ + 1) ≤ max(D, r) − 1, so v·(m′ − π(m′)) ∈ V_{max(D, r)}. By (C2) and linearity
v·π(m′) − π(v·π(m′)) ∈ V_D, and π(v·π(m′)) = π(m). Adding, m − π(m) ∈ V_{max(D, r)}.

*Conclusion.* For a ∈ I_{D'} with D' ≥ D, a = a − π(a) = Σ_m c_m (m − π(m)) over monomials of
degree ≤ D', so a ∈ V_{D'}. By (1.2), LFD ≤ D. ∎

This is the fibre-selector argument of Theorem 1 of `paper/binary-lfd.tex` with two changes: the
Boolean selector χ_α of degree k is replaced by the character idempotent e_α (weight k₁ because
the chain block has k₁ variables), and Gaussian elimination of an affine fibre is replaced by
Lemma 2.3 (the univariate chain system of a polynomial of degree ≤ 2 has last fall degree ≤ 2).
The last step is `LocalReductionBound.lfd_bound_of_local` with σ = k₁ + 1.

### 2.6 Exact check of the proof

`experiments/ic-bounds-upper/selector_bound_check.py` executes the proof literally. Rows are
objects that can only be created by three operations (τ_m g with wt(m) + 2 ≤ budget; linear
combination; variable product of a row of weight ≤ budget − 1), each asserting its budget. The
script builds e_α g, the rows of Lemma 2.3 in each case, the corrections (C1), (C2), and then the
row of the induction for **every** monomial m of A (all degrees up to k₁ + k₂, budget
max(k₁+2, wt m)), and asserts that it equals m − π(m), where π is computed independently from the
roots of the fibre polynomials and is verified to interpolate on Z.

Result (**exact**; `results-selector.jsonl`, `results-selector-summary.json`): 16,844 instances,
all pass. They comprise every nonsingular (a, b) and every t for (p, k) = (5,1), (5,2), (7,1),
(13,2), (17,2), (17,3), (17,4); random instances at (41,3), (97,4), (97,5), (193,5), (193,6),
(257,6), (12289,5), (40961,4), one third of them with t ∈ H (a fibre polynomial of lower degree);
and unequal depths (k₁,k₂) = (1,4), (2,4), (2,5), (3,5), (3,6), (4,6). Every case of Lemma 2.3
occurs (constant, linear with and without a root in H, quadratic with 0, 1, 2 roots in H).

Command: `cd experiments/ic-bounds-upper && python3 selector_bound_check.py` (about 6 minutes,
under 300 MB).

### 2.7 Tightness and comparison with exact data

`lfd_structure.py small` computes the exact LFD (Proposition 4: closures verified closed, ideal
parts from evaluation ranks, generation checked for D' ≤ k+1) on 14,538 instances: every
nonsingular (a, b) and every t for (p, k) = (5,1), (7,1), (13,1), (5,2), (13,2), (17,2), (17,3), and
random instances at (41,3), (97,3), (17,4), (97,4), (113,4), (97,5), (193,5). **Exact**
(`results-structure-small-summary.json`):

| k | k+2 | LFD histogram | bound attained |
|---|---|---|---|
| 1 | 3 | ≤2: 252, 3: 2170 | yes (2170 instances) |
| 2 | 4 | 3: 416, 4: 6336 | yes (6336 instances) |
| 3 | 5 | 4: 4924 | no (max 4) |
| 4 | 6 | 4: 360 | no |
| 5 | 7 | 4: 80 | no |

No instance exceeds k + 2. At k = 1 and k = 2 the bound is attained, so no bound k + c with c < 2
holds for all k (for k = 2 this is also the kernel result LFD ≥ 4 for generic unit ideals). A
concrete witness of LFD = k+2 at k = 2: p = 13, y² = x³ + x + 2, t = 3 (unit ideal):
dim V_2 = 1, dim V_3 = 5 < 15 = |M_≤3|, dim V_4 = 16 = |M_≤4|, so 1 ∉ V_3 and LFD = 4.

Comparison with the certified data (`experiments/ic-curve-sweep`, `ic-chain-depth8`,
`ic-chain-depth9`; two-summand prime rows, subgroup and random-coset bases, both covered by
Corollary 1): the maximal exact LFD is 4, 4, 4, 4, 5, 5, 5, 5 at k = 2, …, 9, against
k + 2 = 4, …, 11, and LFD = 6 at k = 10 against 12. The bound is never violated and is far from
tight for k ≥ 4.

Cross-check of the engine and of Proposition 4 (`lfd_structure.py sweep`, **exact**): on the 311
subgroup two-summand cases of `ic-curve-sweep` with k ≤ 5 and on six cases at k = 6, the LFD
computed here (independent engine, generation criterion) and |Z| agree with the recorded
certified values in 317 of 317 cases (`results-structure-sweep-summary.json`).

### 2.8 Where the argument stops

The budget k₁ + 2 is used in exactly two places: the product e_α·g (weight k₁ times degree 2) and
the variable products of rows w·e_α with wt(w) = 1. Both are forced by the choice of the full
selector e_α, whose weight is k₁ for every α (its top coefficient is α^{−(N₁−1)}/N₁ ≠ 0). A proof
with slope below one cannot localize at single values of x₁:

* any element of A that vanishes on all but one value of x₁ is a multiple b(x₂)·e_α and has
  weight ≥ k₁ unless it is zero;
* partial selectors (indicators of cosets of μ_{2^j} in H₁, weight k₁ − j) leave a system of
  depth j in block 1, for which the same argument costs j + 2: total k₁ + 2 again.

So the selector method gives slope exactly 1, as in the binary case. The exact values
(4, 5, 6 for k ≤ 13) show that the true derivations do not localize.

### 2.9 Remark: more summands

For m summands the same construction with the selector Π_{i<m} e_{α_i} (weight (m−1)k) and a
univariate lemma for fibre polynomials of degree ≤ 2^{m−1} in x_m would give
LFD ≤ (m−1)k + c_m with c_m depending only on m, against the kernel bound mk + 1. This is **not
proved here**: it needs (i) the analogue of Lemma 2.2 (no fibre polynomial S_{m+1}(α₁,…,α_{m−1},X,t)
vanishes identically; where one does, the fibre is all of H and the section needs weight mk), and
(ii) the univariate lemma for degree 2^{m−1}, whose budget is about 2(m−1) because the remainder
basis 1, X, …, X^{2^{m−1}−1} has weight up to m−1. For m = 3 the exact data (LFD = 7 at k = 3, 4, 5,
`experiments/ic-bounds-m3/`, and LFD = 7 at k = 2 on 105 systems) forces c_3 ≥ 3.

---

## 3. Depth restriction

**Proposition 3.** Let π : A_{k+1} → A_k be the algebra map x_i ↦ x_i (restriction of functions
from H_{k+1}² to H_k², i.e. u_{i,k} ↦ 1; one block or both may be restricted). Then for every D ≥ δ

    π(V_D(k+1)) ⊂ V_D(k).

In the ordinary ring: the substitution u_{i,k} ↦ 1 sends each chain equation of depth k+1 to a
chain equation of depth k or to 0, sends the generator to the generator, and does not raise
degrees; so it maps the degree-D closure at depth k+1 into the degree-D closure at depth k.

*Proof.* π is a ring homomorphism, π(g_{k+1}) = g_k, and wt(π(a)) ≤ wt(a) (reducing an exponent
modulo 2^k deletes its bit k). The space {w ∈ M_≤D(k+1) : π(w) ∈ V_D(k)} contains g and is closed
under the closure rule: if wt(w) ≤ D − 1 then wt(π(w)) ≤ D − 1 and π(u_{i,j} w) is u_{i,j}π(w)
(j < k) or π(w) (j = k), both in V_D(k). ∎

**Kernel.** `formal/ICFallDegrees/DepthRestriction.lean` proves the ordinary-ring statement for
any number of summands m, any field and arbitrary further inputs R, R′ whose restrictions are
derivable: `restrict_chainPolynomial`, `restrict_closure`, `one_of_deeper`, `not_one_deeper`,
`lfd_bound_of_deeper_unit`. No root of unity, characteristic or bridge hypothesis is used.
`formal/ICFallDegrees/DepthRestrictionM2.lean` instantiates it for the normalized two-summand
presentation `PrimeM2General.inputs` (the restriction of the normalized generator of depth k+1 is
the one of depth k, `restrict_normalized`), for every field, a, b, t, k ≥ 2 and D ≥ 2:

    one_of_deeper   : InClosure (range (inputs (k+1) _ a b t)) D 1 → InClosure (range (inputs k _ a b t)) D 1
    not_one_of_le   : ¬ InClosure (range (inputs k _ a b t)) D 1 → ∀ k' ≥ k, ¬ InClosure (range (inputs k' _ a b t)) D 1
    lfd_gt_of_le    : the same hypothesis and 1 ∈ ⟨inputs k'⟩ give ∀ e ≤ D, ¬ HasLFDBound (range (inputs k' _ a b t)) e

So any kernel witness of "1 ∉ V_D" at one depth (for instance through
`ChainClosureBridge.prime_m2_not_one`) is a kernel lower bound LFD ≥ D+1 for the unit-ideal
instances with the same (a, b, t) at every larger depth.

**Corollary 3.1 (unit ideals are monotone in the depth; proof).** Fix (p, a, b, t) and D ≥ 2.
If 1 ∈ V_D(k+1) then 1 ∈ V_D(k). Hence, for instances that are unit at depth k+1,
LFD(k) ≤ LFD(k+1); and if 1 ∉ V_D(k₀) then 1 ∉ V_D(k) for every k ≥ k₀ with 2^k ∣ p−1 (in the
ordinary ring: for every k ≥ k₀, no condition on p).

Consequences.

1. A single exclusion certificate at depth k₀ is a lower bound LFD ≥ D+1 for the unit-ideal
   instances at **every** deeper depth. For D = 4 this gives the conclusion of N5
   (`ChainModelMap.n5_lfd_ge_five`) from "1 ∉ V_4(7)" (or already from "1 ∉ V_4(6)") alone: the
   closedness certificate of the depth-8 free extension is needed for the dimension formula
   4k² + 172k − 571, not for the lower bound. For D = 5, "1 ∉ V_5(10)" (certified on the
   p = 12289 instances) gives 1 ∉ V_5(11) and 1 ∉ V_5(12) for the same (a, b, t) with no further
   computation, hence LFD ≥ 6 there whenever the ideal is still the unit ideal.
2. For unit instances the set {k : LFD(k) ≤ D} is an initial segment. The conjecture C1 for
   unit ideals is therefore the pair of statements U_D, L_D of Section 5.1, one depth each.

**Corollary 3.2 (proper ideals; proof).** Suppose LFD(k+1) ≤ D, D ≥ 2.
(a) If Z_k = ∅ then LFD(k) ≤ D (even when Z_{k+1} ≠ ∅).
(b) In general I_d(k) ⊂ V_{max(D, d+2)}(k) for every d.

*Proof.* Let ε = (1 + u_{1,k})(1 + u_{2,k})/4 ∈ A_{k+1}, the indicator of H_k² (weight 2), and
h ∈ I_d(k) with a lift h̃ on the same monomials. Then ε·h̃ vanishes on Z_{k+1}, has weight ≤ d + 2
and π(ε h̃) = h; by hypothesis ε h̃ ∈ V_{max(D, d+2)}(k+1), and Proposition 3 gives (b). For (a) take
h = 1: ε ∈ I(k+1) has weight 2 ≤ D, so 1 = π(ε) ∈ V_D(k). ∎

A statement "LFD(k) ≤ LFD(k+1)" for proper ideals at both depths is not proved; it would need a
lift of weight d instead of d + 2. In the sweep no instance has LFD(k) > LFD(k+1)
(Section 5.2, test T4).

---

## 4. The filtered generation criterion

**Proposition 4.** For D ≥ max(2, δ):

    LFD ≤ D   ⇔   V_D = I_D   and   I_{D'+1} = I_{D'} + Σ_v v·I_{D'}  for every D' ≥ D,

the sum over the chain variables v. By Theorem 1 the second condition holds for every
D' ≥ min(k₁,k₂) + 2, so only finitely many D' have to be checked.

*Proof.* Let W = I_{D'} + Σ_v v I_{D'} ⊂ I_{D'+1}. W contains g, and W ∩ M_≤D' ⊂ I_{D'}, so W is
closed under the closure rule of degree D'+1; hence V_{D'+1} ⊂ W. (⇒) If LFD ≤ D then by (1.2)
V_{D'+1} = I_{D'+1}, so W = I_{D'+1}. (⇐) Induction on D': if V_{D'} = I_{D'} then
V_{D'+1} ⊃ V_{D'} + Σ_v v·V_{D'} = W = I_{D'+1} (the rule applies to all of V_{D'} ⊂ M_≤D'). ∎

So LFD = min{D ≥ γ(Z) : V_D = I_D}, where the **generation degree** γ(Z) is the least D₀ such that
the generation condition holds for all D' ≥ D₀. γ(Z) depends only on the point set Z ⊂ H₁ × H₂
(for a unit ideal γ = 0), and a failure of generation at D' forces LFD ≥ D' + 2.

**Proposition 5 (counting).** If LFD ≤ D then dim I_{D'+1} ≤ (2k+1)·dim I_{D'} for all D' ≥ D
(k₁ = k₂ = k; in general 1 + k₁ + k₂). With r_j = rank of evaluation M_≤j → F_p^Z:
|M_≤D'+1| − r_{D'+1} ≤ (2k+1)(|M_≤D'| − r_{D'}).

*Proof.* Proposition 4 and dim W ≤ (1 + 2k) dim I_{D'}. ∎

**Data (exact).** Generation was computed on all 14,538 instances of Section 2.7 (D' = 2, …, k+1)
and on Z-heavy instances (`slope_counts.py`):

* k ≤ 5, all instances: generation holds for every D' ≥ 2, so γ(Z) ≤ 2 and LFD is the first D with
  V_D = I_D.
* k = 8, p = 257 (H = F_p^*), two curves with |Z| = 262 and 238: dim I_2 = 2, dim I_3 = 435 and 459,
  but I_2 + Σ_v v·I_2 has dimension 34 = 17·2. Generation fails at D' = 2 and holds at D' = 3
  (I_4 = I_3 + Σ_v v I_3, dimensions 2255 and 2279); D' ≥ 4 was not computed (memory). So
  γ(Z) ≥ 3 there, and "LFD ≥ 4" follows from Proposition 5 alone (435 > 17·2), without any closure.
* k = 7, p = 257 (|Z| = 66) and k = 6, p = 193 (|Z| = 28): generation holds at every computed
  D' ≥ 2 (D' ≤ 4 and D' ≤ 5 respectively).

Where γ(Z) was computed in full (k ≤ 6: the 14,538 instances above and the 317 sweep cases of
Section 2.7) it is at most 2, below the LFD: there the binding condition for proper ideals is
V_D = I_D. At k = 8 the generation degree is at least 3, in line with Section 6.2.

---

## 5. The tight direction

### 5.1 What exactly has to be proved

Write k_D = 4D − 11. The conjecture C1, LFD(k) = ⌊(k+14)/4⌋ = min{D : k ≤ k_D}, for a generic
system at every depth k ≥ 2, splits as follows.

**Unit ideals.** By Corollary 3.1 it is equivalent to, for every D ≥ 4,

    (U_D)  1 ∈ V_D(4D − 11),        (L_D)  1 ∉ V_D(4D − 10),

for the generic instance (the instance over F_p(a, b, t), or a Zariski-dense set of
specializations with 2^{4D−10} ∣ p−1). Each is one statement at one depth; the other depths
follow from Proposition 3. Known: U_4, L_4 (kernel/exact, k = 5, 6), U_5, L_5 (exact, k = 9, 10).
First open cases: **U_6: 1 ∈ V_6(13)** and **L_6: 1 ∉ V_6(14)**. The certified "LFD ≥ 6 at
k = 11..13" is, for unit ideals, a consequence of L_5 at depth 10 (Corollary 3.1) and says nothing
about U_6.

**Proper ideals.** By Proposition 4, LFD ≤ D at depth k ≤ 4D − 11 needs

    (a_D) V_D(k) = I_D(k),      (b_D) γ(Z) ≤ D, i.e. I_{D'+1} = I_{D'} + Σ_v v I_{D'} for D ≤ D' ≤ k+1.

(b_D) is a statement about the point set Z alone (regularity of Z for the Hamming filtration); for
k ≤ 6 the computed γ(Z) is at most 2 (Section 4, data), so (a_D) is the binding condition there; by
Proposition 5, (b_D) is expected to fail for D below about 0.22·k when |Z| ≈ 2^k (conditional,
Section 6.2).
There is no monotonicity in k for (a_D) (Corollary 3.2 loses two degrees).

**Lower direction for all depths.** L_D at the single depth 4D − 10 suffices for unit ideals.
The free-extension description (W_k closed) is needed only for the exact dimension of V_D(k)
and for proper ideals.

### 5.2 Candidate inductive statements, tested at k ≤ 7

All tests are **exact** on the instances named; the scripts are in
`experiments/ic-bounds-upper/`. Instances: (p; a, b; t) = (12289; 3, 4; 5) (unit for k ≤ 7),
(40961; 21223, 9887; 4748) (the primary instance of `ic-bounds-fall`; unit for k ≤ 5, |Z| = 2 from
k = 6), (257; 1, 3; 190) (|Z| = 0, 12, 66, 262 at k = 5, 6, 7, 8), (193; 5, 7; 11) (|Z| = 4, 28 at
k = 5, 6).

**T1. Depth restriction (Proposition 3), `depth_restriction.py`.** π(V_D(k+1)) ⊂ V_D(k) in all 63
pairs (D = 3, 4, 5; k = 2, …, 7 for D ≤ 4, k ≤ 6 for D = 5, k ≤ 5 at p = 193; four instances); no pair has
1 ∈ V_D(k+1) and 1 ∉ V_D(k). This only confirms the proposition and the engine.

**T2. Surjectivity (E1).** In all 63 pairs **π(V_D(k+1)) = V_D(k)**: the image has the full
dimension of V_D(k). This includes the transitions where V_D(k) contains 1 and V_D(k+1) does not
(D = 4, k = 5 → 6: the 605-dimensional V_4(6) maps onto all of M_≤4 at depth 5, dimension 386) and
proper ideals on both sides (p = 257, D = 5, k = 6 → 7: 3407 → 1574 = |M_≤5| − 12).
Status of E1:

* above the base it is a consequence of the free extension (**proof**): if
  V_D(k+1) = ⊕_h τ_h(S ∩ M_≤D−wt h) with h supported on levels ≥ k′ for some k′ ≤ k, then deleting the
  level-k variables from h maps V_D(k+1) into itself, does not change the image under π and does
  not raise degrees; so even the filtered statement π(V_D(k+1) ∩ M_≤j) = V_D(k) ∩ M_≤j holds;
* below the transition, for unit ideals, both sides are everything;
* **at the transition it contains the upper statement**: for a unit ideal at depth k₀ − 1,
  π(V_D(k₀)) = V_D(k₀ − 1) = M_≤D says in particular that 1 has a lift in V_D(k₀). So E1 is not a
  reduction of U_D to something easier; the part of it that matters is the sufficient condition

      (S_D)  there is v ∈ V_D(4D − 10) with π(v) = 1, i.e. v ≡ 1 on H_{4D−11}²    ⇒   U_D   (**proof**, Prop. 3),

  and S_4 holds (**exact**, three instances, `lift_of_one.py`).

**T2′. The filtered version fails at the transition.** dim π(V_D(k+1) ∩ M_≤j) for j = 0, …, D
against dim(V_D(k) ∩ M_≤j): at D = 4, k = 5 → 6 it is (0, 0, 2, 107, 386) against (1, 11, 56, 176, 386).
The least degree of an element of V_4(6) that restricts to 1 at depth 5 is 4 = D on the three
instances tested (`results-lift-of-one.json`). So the wrap at depth 5 produces low-degree
elements (falls) whose lifts at depth 6 exist only in full degree: the candidate "low-degree parts
lift" is false, and an induction on the depth cannot proceed degree by degree.

**T3. Low-level part.** dim(V_D(k+1) ∩ M_{levels<k}) equals dim V_D(k) at D = 4, k = 7 (829; this is
`ic-bounds-fall` 5.2) but is smaller at k = 6 (563 < 605) and k = 5 (270 < 386): the candidate
"V_D(k) embeds in V_D(k+1) on the same monomials" is false at and one step above the transition.

**T4. Monotonicity in the depth on certified values, `sweep_monotone.py`.** 339 subgroup
instances of the sweeps, 35 of them computed at several depths, 254 pairs k < k′ (45 with a unit
ideal at k′, where LFD(k) ≤ LFD(k′) is Corollary 3.1, and 209 with a proper ideal at k′): no pair
with LFD(k) > LFD(k′). So "LFD is nondecreasing in k for fixed (p, a, b, t)" is **empirical** for
proper ideals.

**T5. Generation degree (Proposition 4).** Section 4, data: γ(Z) ≤ 2 on every instance with
k ≤ 6 where it was computed in full; the candidate "γ(Z) is the obstruction" is false in that
range (V_D = I_D is). At k = 8, p = 257, generation fails at D' = 2.

**T6. Candidate "LFD(k+1) ≤ LFD(k) + 1".** True on all certified values (it is implied by C1); it
would give slope ≤ 1, which Theorem 1 already gives. The candidate needed for C1 is (Step_D) below;
its only available instances are D = 4 (k = 5 → 9, true) and D = 5 (k = 9 → 13, open: it is U_6).

### 5.3 What an induction on the depth would need

The spanning lemma that would prove the upper half of C1 by induction is

    (Step_D)   1 ∈ V_D(k)  ⇒  1 ∈ V_{D+1}(k + 4)        (generic unit instance, k = 4D − 11).

Together with U_4 it gives U_D for all D. Three facts constrain any proof of (Step_D).

1. *Lifting costs two degrees per level with selectors.* A_{k+j} = ⊕_η A_k^η over the 4^j
   characters η of the top j levels, and each component is a twisted depth-k system
   (x_i^{2^k} = η_i). Multiplying a depth-k derivation by the selector of η (weight 2j) gives
   LFD(k+j) ≤ max_η LFD^η(k) + 2j for unit ideals (the argument of Lemma 2.4 with selectors in
   both blocks; not written out). (Step_D)
   asks for +1 where this gives +8; Theorem 1 gives +4.
2. *The wrap is not an independent source of rows (E1, empirical), but this does not reduce
   U_D.* By T2, E1 at the transition depth contains U_D, and above it E1 follows from
   the free extension. What remains is a concrete form of the target: V_D(4D − 10), which omits 1
   and is the base of the free extension for D = 4 and 5, must contain an element congruent to 1
   modulo (u_{1,4D−11} − 1, u_{2,4D−11} − 1) (condition S_D). An inductive proof would have to
   produce this element of V_{D+1}(4D − 6) from the one in V_D(4D − 10); the profiles (Section 6.1)
   show no simple relation between the two bases (Q_4 and Q_5 are not truncations of one series).
3. *The filtered version of E1 is false* (test T2′): elements of V_D(k) of low degree have lifts
   in V_D(k+1) only of higher degree. So an induction cannot be run degree by degree on the
   filtration; it has to control dim(V_D(k+1) ∩ M_≤j) for all j at once, which is the profile
   recursion p_j(k+1) = p_j(k) + 2p_{j−1}(k) + p_{j−2}(k) above the base and is unknown below it.

No statement of the form "V_D(k) = I_D(k) and V_D(k+1) ⊇ (shifts) ⇒ V_D(k+1) = I_D(k+1)" was
found that is both true on the data and inductive: the only general relation between consecutive
depths that was proved is Proposition 3, which goes downward.

---

## 6. Dimension counts and the slope

### 6.1 The free-extension polynomial crosses |M_≤D| at the conjectured base depth

If V_D(k+1) is the free extension of V_D(k), the profile p_j(k) = dim(V_D(k) ∩ M_≤j) satisfies
p_j(k+1) = p_j(k) + 2p_{j−1}(k) + p_{j−2}(k), i.e. with P_k(z) = Σ_j p_j(k) z^j

    P_k(z) ≡ (1 + z)^{2k}·Q_D(z)   mod z^{D+1}

for a fixed polynomial Q_D (a "virtual depth-0 profile"), and
f_D(k) := dim V_D(k) = Σ_i Q_{D,i}·C(2k, D − i). From the certified profiles (**exact**
arithmetic, `slope_counts.py`, `results-slope-counts.json`):

    Q_3 = z² + z³                              f_3(k) = 2k + 1
    Q_4 = 2z² + 87z³ − 571z⁴                   f_4(k) = 4k² + 172k − 571
    Q_5 = 2z² + 98z³ + 823z⁴ − 18791z⁵         f_5(k) = (8k³ + 576k² + 4648k − 56373)/3

(the four certified D = 5 profiles at k = 10, 11, 12, 13 give the same Q_5; three of them are
needed to determine it, the fourth is a check). Compared with |M_≤D(2k)| = Σ_{j≤D} C(2k, j):

| D | k | f_D(k) | \|M_≤D(2k)\| | |
|---|---|---|---|---|
| 4 | 4 | 181 | 163 | free extension impossible |
| 4 | 5 | 389 | 386 | free extension impossible |
| 4 | 6 | 605 | 794 | first possible depth = 4D − 10 |
| 5 | 8 | 7257 | 6885 | impossible |
| 5 | 9 | 12649 | 12616 | impossible |
| 5 | 10 | 18569 | 21700 | first possible depth = 4D − 10 |

So for D = 4 and D = 5 the closure is the free extension **from the first depth at which the
count allows it** (and it is everything, for unit ideals, at the depth before). In words: V_D(k)
is as large as possible until the free-extension polynomial fits. This is a consistency
observation on two degrees, not a derivation: Q_D is read off the data, and below k = 4, 8 the
polynomial is not monotone against |M_≤D| (f_4(3) < 0).

What it would take to turn this into a prediction of the slope: a formula for Q_D. With
"maximal growth" (**conjecture**: V_D(k) = M_≤D for unit ideals exactly while f_D(k) > |M_≤D(2k)|)
the base depth is 1 + max{k : f_D(k) > |M_≤D(2k)|}, and the conjectured value 4D − 10 is the
statement

    f_D(4D − 11) > |M_≤D(8D − 22)|   and   f_D(4D − 10) ≤ |M_≤D(8D − 20)|.

For large D, |M_≤D(8D)| = 2^{8D·H(1/8)(1+o(1))} ≈ 20.4^D, so slope exactly 1/4 requires
f_D(4D) to grow like 20.4^D; a slope 1/c requires growth 2^{2cD·H(1/(2c))}. The two available
values, f_4(6) = 605 and f_5(10) = 18569 (ratio 30.7, against |M| ratio 27.3), do not determine
a growth rate, and the leading term 2·C(2k, D−2) of f_D is a small part of |M_≤D| at the base
(17% at D = 4, 11% at D = 5, 2/49 asymptotically at k = 4D): the size of f_D at the base is carried
by the lower coefficients Q_{D,3}, …, which are not explained. **No heuristic derivation of the slope 1/4 from dimension counts was
found**; the count only says that the base depth is where it has to be given Q_D.

A direct test would be the D = 6 profile at one depth k ≥ 14 on one instance (it determines Q_6
together with the recursion checked at a second depth); the prediction of C1 is
f_6(13) > |M_≤6(26)| = 313,912 and f_6(14) ≤ |M_≤6(28)| = 499,178.

### 6.2 A counting obstruction of order 0.22·k for primes close to 2^k

Proposition 5 gives a lower bound that does not involve closures. Let D*(k) be the least D with
|M_≤D(2k)| ≥ |Z|. For D' < D* the evaluation map on M_≤D' cannot be onto, and if it is injective
up to the few low-degree ideal elements (dim I_{D'} small) while dim I_{D'+1} ≥ |M_≤D'+1| − |Z|
is large, generation fails at D' and LFD ≥ D' + 2.

*Exact instance.* k = 8, p = 257, |Z| = 262: |M_≤2| = 137 < |Z| ≤ |M_≤3| = 697, so D* = 3;
dim I_2 = 2 (the elements g and w of `ic-bounds-fall`), dim I_3 = 697 − 262 = 435 > 17·2, hence
LFD ≥ 4 = D* + 1. (The true value is 5.)

*Thresholds.* With |Z| = 2^k (the order of magnitude of |Z| ≈ 2^{2k}/p when p is close to 2^k;
always |Z| ≤ 2^{k+1}), D*(k) = min{D : |M_≤D(2k)| ≥ 2^k}:

| k | 8 | 16 | 24 | 32 | 48 | 64 | 96 | 128 |
|---|---|---|---|---|---|---|---|---|
| D*(k) | 3 | 5 | 7 | 8 | 12 | 16 | 23 | 30 |
| ⌊(k+14)/4⌋ | 5 | 7 | 9 | 11 | 15 | 19 | 27 | 35 |
| k + 2 (Theorem 1) | 10 | 18 | 26 | 34 | 50 | 66 | 98 | 130 |

Asymptotically |M_≤D(2k)| = 2^{2k·H(D/2k)(1+o(1))}, so D*(k)/k → 2·H^{−1}(1/2) = 0.2200.

*Conditional statement.* If dim I_{D*−1} < (|M_≤D*| − |Z|)/(2k+1), then LFD ≥ D* + 1. For
k = 16, p = 65537 (|Z| ≈ 2^16, D* = 5) the condition is dim I_4 < 5373; the free-extension values
would give dim V_4(16) = 3205 and dim(V_5(16) ∩ M_≤4) = 4951, both below 5373, but I_4 can be
larger than these and has not been computed (the evaluation matrix is 65,537 × 41,449). So
"LFD ≥ 6 at k = 16, p = 65537" is a **conjecture** supported by a count, and in general the
counting bound D* + 1 ≈ 0.22·k is **conditional** on the ideal having few elements below degree
D*. It stays below ⌊(k+14)/4⌋ by 1 to 4 for k ≤ 128, so it does not contradict C1; it shows that
for primes close to 2^k a slope of at least 0.22 is forced by the size of Z if the low-degree
part of the ideal is as small as at k = 8. For unit ideals (p ≫ 4^k) Proposition 5 is empty.

### 6.3 What a proof of a bound with slope below 1 needs

1. It cannot use selectors of single values of one coordinate (Section 2.8), nor the dimension
   of the quotient. It has to produce the section corrections v·b − π(v·b) (Section 2.5), or for
   unit ideals the element 1, through elements of weight < k: equivalently, by Proposition 4, it
   has to prove (a) V_D = I_D at one degree D = D(k) and (b) γ(Z) ≤ D.
2. For (b) Proposition 5 shows that generation fails at D' = D*(k) − 1, hence D ≥ D*(k) + 1, when
   the ideal has few elements below D*; for |Z| ≈ 2^k this is D ≳ 0.22·k. So a uniform bound (all primes with
   2^k ∣ p−1, all curves and targets) with slope below 0.22 is excluded if the smallness of
   dim I_{D'} below the threshold holds (it is exact at k = 8, p = 257; **conditional** beyond).
   For primes p ≫ 4^k (unit ideals) there is no counting obstruction at all, and the conjectured
   slope 1/4 there rests only on the closure computations.
3. For (a) in the unit case the statement is U_D of Section 5.1. By Proposition 3 it is enough to
   prove it at the largest depth 4D − 11, and it follows from S_D (Section 5.2), a statement about
   the single space V_D(4D − 10).
4. Linear growth of the true LFD with slope ≥ 1/4 for unit ideals would follow from L_D for all
   D, i.e. from one exclusion certificate per degree at depth 4D − 10 (Corollary 3.1); the
   free-extension recursion is not needed for that.

---

## 7. Lean formalization of Theorem 1 (outline, not started)

Target, in `ICFallDegrees.ChainQuotient` style (two summands, equal depths, K = ZMod p):

    theorem prime_m2_lfd_le_depth_add_two (p : ℕ) [Fact p.Prime] (large : 5 ≤ p)
        (k : ℕ) (depth : 2 ≤ k) (root : (2^k) ∣ p - 1) (a b t : ZMod p)
        (nonsingular : 4*a^3 + 27*b^2 ≠ 0) :
        HasLFDBound (Set.range (PrimeM2General.inputs k depth a b t)) (k + 2)

Route: `LocalReductionBound.lfd_bound_of_local` with the quotient map q : P (2k) K →ₐ A′
(A′ = K-valued functions on Z, or P/⟨F⟩), the section s onto span{ρ·e_α}, σ = k + 1, D = k + 2.
Its hypothesis `local_reductions` (p − s(q p) ∈ InClosure F (k+2) for deg p ≤ k+2) is the induction
of Section 2.5 restricted to degree ≤ k+2 and is itself obtained from the corrections by
`LocalReductionBound.global_reduction_of_local`-style induction on monomials.

Lemmas needed, with rough sizes (lines of Lean, by comparison with existing files):

| lemma | content | size |
|---|---|---|
| `selector` | e_α as the polynomial Π_j (1 + C(α^{-2^j}) X(var 0 j))/2; totalDegree ≤ k; X(var 0 j)·e_α − C(α^{2^j})·e_α ∈ chain ideal with degree ≤ k+1, hence in `InClosure` at budget k+1 (by `chain_mem_closure_iff`, normal form 0) | 150 |
| `selector_sum` | Σ_α e_α − 1 ∈ chain ideal, degree ≤ k (partition of unity; needs the 2^k-th roots of unity of ZMod p, `ZMod.pow_card_sub_one`, `IsPrimitiveRoot`) | 150 |
| `fibre_nonzero` | Lemma 2.2 (`ring`/`linear_combination` identities and a discriminant argument) | 60 |
| `univariate_closure` | Lemma 2.3 stated for rows multiplied by a fixed e of degree ≤ k with X(var 0 j)·e ≡ c_j e: the recursion w_j, w′_j as `InClosure` terms; the three cases r = 0, 1, 2 via gcd in `Polynomial (ZMod p)` or, avoiding gcd, via the rank of the 2×2 matrix (w_k, w′_k) | 350 |
| `section` | s, q, `section_id`, `bounded` (degree ≤ k+1), `vanish` (q kills chains and g) | 200 |
| `corrections` | (C1), (C2) | 150 |
| assembly | `lfd_bound_of_local` | 50 |

Estimated total: 1,000–1,300 lines, no heavy kernel computation (no `decide` on data). The main
risk is the bookkeeping between the polynomial ring and functions on H × H (the existing
`ChainQuotient.expand` and the grid vanishing theorem used for `reduced_chain_ideal_trivial` cover
"a reduced polynomial vanishing on the grid is zero", which gives `section_id` and the selector
identities). A cheaper intermediate target is the model statement in the exponent model of
`DepthReduction` / `ChainModelMap` (Theorem 1 in A), followed by the existing bridge.

Kernel-checked now: Proposition 3 and Corollary 3.1 only (`DepthRestriction.lean`,
`DepthRestrictionM2.lean`, Section 3).

---

## 8. Files

* `docs/proofs/U1-upper-bounds.md` (this file).
* `formal/ICFallDegrees/DepthRestriction.lean`, `formal/ICFallDegrees/DepthRestrictionM2.lean`
  (Proposition 3 and Corollary 3.1, kernel; `lake build ICFallDegrees.DepthRestrictionM2`; axioms
  `propext`, `Classical.choice`, `Quot.sound`).
* `experiments/ic-bounds-upper/`: `README.md` (table of scripts, outputs and run times),
  `chain_algebra.py` (independent exact engine), `selector_bound_check.py`, `lfd_structure.py`,
  `depth_restriction.py`, `lift_of_one.py`, `sweep_monotone.py`, `slope_counts.py`, and their
  `results-*` outputs.

All scripts run under `ulimit -v 1500000` with one thread.
