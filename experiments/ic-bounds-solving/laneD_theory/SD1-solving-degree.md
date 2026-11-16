# SD1: the static (XL) solving degree of the prime chain systems, beside d_ff, d_lf and sd^CG

Status (2026-11-15; written as lane D of the study `crypto:research/solving_degree_20261115/`, data lanes in `experiments/ic-bounds-solving/`). Evidence labels:
**kernel** (Lean-checked), **proof** (written out in full here), **exact** (certified computation,
command and output named), **empirical** (pattern in data), **conjecture**.

The proofs below have not been through independent review. Every result is
**[pcc:unverified]**: no pcc claim id has been issued. An assertion that is not a proof written
out here is also tagged [pcc:unverified] at the point where it is made. Red-team record:
`experiments/ic-bounds-solving/laneD_theory/REDTEAM.md`. Scripts: `experiments/ic-bounds-solving/laneD_theory/tools/check_small.py`, `tools/exceptional_count.py` (outputs `tools/*.out`).

**Naming.** PROTOCOL.md calls the static quantity "the Caminata–Gorla solving degree". That is
corrected here (REDTEAM A0). The Caminata–Gorla solving degree sd^CG is defined with mutant
rows: new lower-degree rows are multiplied again. For degree-compatible orders, CG Theorem 3.1
gives sd^CG = max(d_lf, max.GB.deg). The no-mutant quantity "least D with 1 ∈ Mac_D" is the
**static** or XL solving degree, written **sd^stat** throughout this note.

## 0. Results

| # | statement | label |
|---|---|---|
| L | **Lemma.** For chains plus one generator ĝ of degree δ, with g = NF(ĝ) ∈ A = K[(Z/2^k)^m]: (a) ker NF ∩ R_{≤D} = span of chain multiples of degree ≤ D, contained in Mac_D; (b) NF(Mac_D) = g·M_{≤D−δ}. With several generators, NF(Mac_D) = Σ_i g_i·M_{≤D−δ_i}. Any field K, any characteristic. | **proof** (§2) [pcc:unverified] |
| T | **Theorem.** If ⟨F⟩ = R (equivalently g ∈ A^×, equivalently ĝ has no zero on the chain variety), then **sd^stat = δ + maxwt(supp g^{-1})**, with δ = deg ĝ, the ordinary degree. | **proof** (§3) [pcc:unverified]; **exact** on 68/68 unit instances, k ≤ 3 (REDTEAM §B) |
| C1 | d_lf = sd^CG ≤ sd^stat ≤ δ + mk for unit ideals; sd^stat = δ + mk ⇔ the coefficient c_top of x^{(N−1,…,N−1)} in g^{-1} is nonzero. | **proof** (§4.1); sd^CG = d_lf cites CG Thm 3.1 |
| C2 | **c_top = N^{-m} Σ_{α∈H^m} (α₁⋯α_m)/ĝ(α)** (ĝ evaluated at u_{i,j} = α_i^{2^j}); for m = 2 with the normalized S_3: c_top(t) = N^{-2} Σ_{α₁,α₂∈H} α₁α₂ / S_3(α₁, α₂, t). Needs p ∤ N and H ⊂ K. | **proof** (§4.2); **exact** agreement with the linear-solve inverse on 68/68 unit rows |
| C3 | Fix (a, b), m = 2, normalized S_3. Either c_top(t) ≡ 0, or **#{t ∈ F_p : instance unit, c_top(t) = 0} ≤ N² − 1 = 4^k − 1**. For general m and an S_m-symmetric generator of t-degree d_t, the bound is d_t·C(N+m−1, m). Proportion O(1/p) for fixed k; vacuous once N² ≥ p. | **proof** (§4.3); the bound is attained at k = 1 (**exact**, `tools/exceptional_count.out`) |
| C4 | Generic gap: under U1 (d_lf ≤ k + 2), every unit instance with c_top ≠ 0 has **sd^stat − d_lf ≥ k** (m = 2, δ = 2). With the measured law, the gap is about 2k + 2 − ⌊(k+14)/4⌋ ≈ 7k/4. | **proof** given U1; measured-law form **conjecture** |
| P | Proper ideals: Mac_D ⊇ ⟨F⟩ ∩ R_{≤D} ⇔ g·M_{≤D−δ} ⊇ gA ∩ M_{≤D}; this holds for all D ≥ δ + mk. The static GB degree is max_{G∈𝒢} max(deg G, δ + w(G)), with w(G) = min{wt h : g h = NF G}. sd^CG = max(d_lf, max.GB.deg) ≤ max(d_lf, 2, mk). | **proof** (§5); sd^CG formula cites CG Thm 3.1 |
| B | Binary (Boolean ring with field equations): Lemma (a) and the several-generator form of (b) survive with squarefree weight. The Theorem does not transfer: the only unit of a Boolean ring is 1. sd^CG ≤ D* (the programme's certificate degree), hence sd^CG − LFD ≤ D* − LFD ∈ {0, 1, 3} on the recorded rows. | **proof** (§6); recorded counts **exact** in the source but not recomputed [pcc:unverified] |
| N | Measured: m = 2, sd^stat = 2k + 2 on every unit instance at k = 5..13 (29,616 generic rows, 159 exceptions, all at k = 3, 4 with c_top = 0 and sd^stat = 2k + 1, at a rate of about 1/p); brute-force Macaulay agreement with 2 + maxwt on 18,828/18,828 instances (k = 3..6); m = 3, sd^stat = 6 + 3k on generic units k = 2..8, brute force 40/40 at k = 2, 3; binary proper systems d_full^stat = k + 1 or k + 2 on 185/185 exact systems, k = 2..8. | **exact** (lanes A, B, C; `experiments/ic-bounds-solving/`) [pcc:unverified] |

What is not claimed: no lower bound sd^stat ≥ δ + f(k) without genericity (REDTEAM A8 gives a
counterexample with sd^stat = 2 at every k); no proof that c_top ≢ 0 for all nonsingular curves;
no proof that the exceptional proportion is O(1/p) in the repository regime p ≈ c·2^k; no novelty
claim (research-kb was unavailable, REDTEAM §C).

---

## 1. Setting

Let K be a field and k, m ≥ 1, N = 2^k. The ordinary ring is R = K[u_{i,j}], 1 ≤ i ≤ m,
0 ≤ j < k, with the total degree. The chains are

    c_{i,j} = u_{i,j}² − u_{i,j+1}   (j < k−1),        c_{i,k−1} = u_{i,k−1}² − 1.

There is one further generator ĝ ∈ R of total degree δ := deg ĝ, so F = {c_{i,j}} ∪ {ĝ}. In the
repository ĝ is the normalized generator: the chain normal form of S_3(x₁, x₂, t) for m = 2 with
δ = 2, or of S_{m+1} in general.

**Quotient model.** Let G = (Z/N)^m and A = K[G], with basis monomials x^e (e ∈ G) and degree
**wt(e)** = Σ_i popcount(e_i). For a ∈ A, wt(a) is the largest weight in its support (−∞ for 0).
M_{≤d} ⊂ A is the span of the x^e with wt(e) ≤ d. The ring map

    NF : R → A,   u_{i,j} ↦ x_i^{2^j}

is surjective with kernel the chain ideal C = ⟨c_{i,j}⟩ (§2). A reduced monomial is one that is
squarefree in every u_{i,j}. The reduced monomial Π u_{i,j}^{b_{ij}} maps to x^e with
e_i = Σ_j b_{ij} 2^j, and its degree is wt(e). So NF restricts to a degree-preserving linear
bijection ι^{-1} from reduced polynomials to A; write ι : A → R for the inverse (the reduced lift).

**Weight inequality.** For a ∈ A and an ordinary monomial μ, wt(NF(μ)·a) ≤ deg μ + wt(a), and in
particular

    wt(a·b) ≤ wt(a) + wt(b)                                                       (1.1)

(U1 §1, (1.1)). Proof: multiplication by u_{i,j} adds 2^j to e_i mod N. This clears a run of
r ≥ 0 ones and sets at most one bit, so the weight changes by at most 1 − r ≤ 1.

**Invariants (fixed for the study).** Write V_{F,D} for the iterated degree-D closure (IC_FALL_PROOFS
§1; CG Definition 1.5; HKYY Definition 2.1).

* **d_ff** (affine first strict fall): the least D with V_{F,D} ∩ R_{≤D−1} ≠ V_{F,D−1}.
* **d_lf** = d_F = min{d : f ∈ V_{F,max(d, deg f)} for every f ∈ ⟨F⟩}. For a unit ideal this is
  the least D with 1 ∈ V_{F,D}.
* **sd^CG** (Caminata–Gorla Definition 1.1): the least d such that rowsp(M_d) contains a Gröbner
  basis of ⟨F⟩, where M_d is the Macaulay matrix in degree d closed under the mutant step.
  For degree-compatible σ, CG equation (1) gives rowsp(M_d) = V_{F,d}, and **CG Theorem 3.1**
  gives sd^CG_σ = max(d_lf, max.GB.deg_σ). For a unit ideal, the reduced GB is {1}, so
  **sd^CG = d_lf**. (Cited. The quotes were checked against the arXiv PDF by the coordinator,
  not by lane D.)
* **Mac_D** = span{μ·f : f ∈ F, μ a monomial, deg(μ f) ≤ D}: the static Macaulay row space, with
  no mutant step. **sd^stat** = the least D with 1 ∈ Mac_D (unit ideals). This is the degree that
  plain XL, or a single Macaulay matrix without re-multiplication of new rows, must reach.

Mac_D ⊆ V_{F,D}, because V_{F,D} contains F ∩ R_{≤D} and is closed under multiplication by
variables within degree D. Hence **d_lf ≤ sd^stat** for unit ideals.

## 2. Lemma (the static Macaulay space in the quotient)

**Lemma 2.1.** For every D ≥ 0:

(a) ker NF ∩ R_{≤D} = span{μ·c_{i,j} : deg(μ c_{i,j}) ≤ D} ⊆ Mac_D.

(b) {NF(μ) : μ an ordinary monomial, deg μ ≤ d} = {x^e : wt(e) ≤ d}. Hence span NF(R_{≤d}) = M_{≤d},
and

    NF(Mac_D) = g · M_{≤D−δ}     (M_{≤d} = 0 for d < 0).

For several generators ĝ_1, …, ĝ_s of degrees δ_i, NF(Mac_D) = Σ_i g_i·M_{≤D−δ_i}.

*Proof.* (a) Fix a degree-compatible order with u_{i,j}² > u_{i,j+1}, for example DRL with [pcc:unverified]
u_{i,j} > u_{i,j+1}. The leading monomial of c_{i,j} is u_{i,j}². These leading monomials are
pairwise coprime, so the chains form a Gröbner basis of C (Buchberger's first criterion). Let
P ∈ R_{≤D}. One division step picks a term c·u_{i,j}²·μ′ of the current remainder and subtracts
c·μ′·c_{i,j}. The subtracted multiple has degree deg(u_{i,j}² μ′), which is at most the degree
of the current remainder. That degree never increases, so it stays ≤ deg P ≤ D. Division ends
at a reduced polynomial r with P − r ∈ span{μ c_{i,j} : deg ≤ D}. NF is a ring map that kills
each c_{i,j}, so NF(P) = NF(r), and NF is injective on reduced polynomials. If NF(P) = 0 then
r = 0, which proves "⊆". "⊇" holds because NF(c_{i,j}) = x_i^{2^{j+1}} − x_i^{2^{j+1}} = 0 for
j < k−1, and x_i^{N} − 1 = 0 for j = k−1. The chain multiples are rows of Mac_D by definition.

(b) Each chain rewrite u_{i,j}² ↦ u_{i,j+1} or ↦ 1 maps a monomial to a monomial with coefficient
1. So NF(μ) = x^e for one e, and wt(e) ≤ deg μ by (1.1), applied deg μ times starting from 1.
Conversely, x^e = NF(ι(x^e)) with deg ι(x^e) = wt(e). For the second statement, R is a domain,
so deg(μ ĝ) = deg μ + δ, and the rows μ ĝ of Mac_D are exactly those with deg μ ≤ D − δ. NF is
multiplicative and kills the chain rows. Therefore

    NF(Mac_D) = span{NF(μ)·g : deg μ ≤ D − δ} = g·M_{≤D−δ}.

The several-generator case is the same computation, term by term. ∎

*Remarks.* (i) The proof uses no property of K: any field and any characteristic will do, and
the same argument works over any commutative ring because the leading coefficients are 1.
(ii) δ is the **ordinary** degree of ĝ. If ĝ is not reduced, δ > wt(g) is possible (raw S_3:
δ = 4, wt g = 2), and the formula needs δ, not wt g (REDTEAM A3, checked exactly).
(iii) Kernel counterpart for the iterated closure: `ChainClosureBridge.chain_mem_closure_iff`
(ml-cryptanalysis, `formal/ICFallDegrees/ChainClosureBridge.lean:404`). There is no Lean
counterpart yet for the static statement.

## 3. Theorem (static solving degree of a unit instance)

**Theorem 3.1.** The following are equivalent: (i) ⟨F⟩ = R; (ii) g ∈ A^×; (iii) over K̄, ĝ has
no zero on the chain variety {u_{i,j} = α_i^{2^j} : α_i^N = 1}. When they hold,

    sd^stat = δ + maxwt(supp g^{-1}).

*Proof.* ⟨F⟩ ⊇ C = ker NF, so ⟨F⟩ = NF^{-1}(gA). Hence ⟨F⟩ = R ⇔ gA = A ⇔ g is a unit, which
gives (i) ⇔ (ii). (i) ⇔ (iii) is the weak Nullstellensatz: V(F) is the zero set of ĝ on V(C),
and V(C) is the stated set over K̄.

Let w = maxwt(supp g^{-1}). We show 1 ∈ Mac_D ⇔ D ≥ δ + w.

⇒: 1 ∈ Mac_D gives 1 = NF(1) ∈ NF(Mac_D) = g·M_{≤D−δ} by Lemma 2.1(b). So g h = 1 for some
h ∈ M_{≤D−δ}. Inverses are unique, so h = g^{-1}, and w ≤ D − δ.

⇐: Let D ≥ δ + w and P = ι(g^{-1})·ĝ. Expanding ι(g^{-1}) into monomials of degree ≤ w ≤ D − δ
writes P as a combination of rows μ ĝ of Mac_D. Also deg(P − 1) ≤ D and
NF(P − 1) = g g^{-1} − 1 = 0. By Lemma 2.1(a), P − 1 is a combination of chain rows of degree
≤ D. (If D < 2 there are no such rows, but P − 1 is then reduced with NF 0, so P − 1 = 0.)
Hence 1 = P − (P − 1) ∈ Mac_D.

Mac_D increases with D, so the least D is δ + w. ∎

*Exact check.* `tools/check_small.py` (own sparse elimination over the ordinary ring; no
repository engine). It covers 174 instances with (m, k) ∈ {(2,1), (2,2), (2,3), (3,1)} and
p ∈ {3, 5, 7, 13, 17, 41, 73}, with generators: normalized S_3, raw S_3, random reduced
polynomials of degree 1–3, and u_{1,k−1}. On all **68 unit instances** sd^stat equals
δ + maxwt supp g^{-1}. On all **106 non-unit instances**, 1 ∉ Mac_D up to δ + mk + 1. There
were **0 mismatches** (REDTEAM §B). [pcc:unverified]

## 4. Corollaries

### 4.1 Bounds and the top coefficient

**Corollary 4.1.** For a unit instance, d_lf = sd^CG ≤ sd^stat ≤ δ + mk. Moreover
sd^stat = δ + mk iff c_top := [x^{e*}] g^{-1} ≠ 0, where e* = (N−1, …, N−1).

*Proof.* d_lf ≤ sd^stat by Mac_D ⊆ V_{F,D}. sd^CG = d_lf is CG Theorem 3.1 with GB = {1}
(cited). wt(e) ≤ mk with equality iff e = e*, because every bit is set only at e*. Theorem 3.1
finishes the proof. ∎

For the two-summand normalized system (δ = 2) this reads sd^stat ≤ 2k + 2, with equality iff
c_top ≠ 0. Three summands give δ_3 + 3k, with δ_3 = deg NF(S_4) = 6 generically (lane B: δ_3 = 6 on every generic unit instance tested at k = 2..8, 5 on the locus 5t³ + 3at + 6b = 0) [pcc:unverified].

### 4.2 The top coefficient as a character sum

Assume char K = p ∤ N and H := μ_N ⊂ K (for K = F_p: N | p − 1). Then A ≅ K^{H^m} via
a ↦ (a(α))_α, where a(α) = Σ_e a_e α^e and α^e = Π α_i^{e_i}. This is the discrete Fourier
transform on G, an isomorphism because the characters e ↦ α^e (α ∈ H^m) are the |G| = N^m
distinct characters of G and N is invertible in K. The inversion formula is

    a_e = N^{-m} Σ_{α∈H^m} a(α) α^{-e}.

Apply it to a = g^{-1}. Then a(α) = 1/g(α), and α^{-e*} = Π α_i^{-(N−1)} = Π α_i. So:

**Corollary 4.2.**

    c_top = N^{-m} Σ_{α∈H^m} (α₁ α₂ ⋯ α_m) / g(α).

Here g(α) = ĝ(u_{i,j} = α_i^{2^j}), because NF is evaluation-compatible. For m = 2 and the
normalized generator,

    c_top(t) = N^{-2} Σ_{α₁,α₂∈H} α₁ α₂ / S_3(α₁, α₂, t),
    S_3(x₁,x₂,t) = (x₁−x₂)² t² − 2((x₁+x₂)(x₁x₂+a) + 2b) t + (x₁x₂ − a)² − 4b(x₁+x₂).

Equivalently, sd^stat < δ + mk iff Σ_α (Π α_i)/g(α) = 0.

*Proof.* As above. The value of g at α is the value of ĝ at the corresponding point of V(C)
because ĝ − ι(g) ∈ C. ∎

*Exact check.* On all 68 unit rows, c_top from the character sum equals the coefficient of the
linearly solved inverse, and c_top ≠ 0 ⇔ sd^stat = δ + mk (REDTEAM §B). [pcc:unverified]

### 4.3 How often the top coefficient vanishes

Fix m = 2, p odd with N | p − 1, the curve coefficients (a, b), and let t range over F_p. For
α ∈ H² put q_α(t) := S_3(α₁, α₂, t) ∈ F_p[t]. Then deg_t q_α ≤ 2, and deg_t q_α ≤ 1 on the
diagonal α₁ = α₂, because the t² coefficient (α₁ − α₂)² vanishes there. S_3 is symmetric, so
q_α depends only on the orbit O = {α, ᾱ} of the swap α ↦ ᾱ = (α₂, α₁), and so does the weight
α₁α₂. Write U ⊂ F_p for the set of unit targets (q_α(t) ≠ 0 for all α), and
E = {t ∈ U : c_top(t) = 0} for the exceptional set.

**Corollary 4.3.** Either c_top vanishes on all of U, or |E| ≤ N² − 1 = 4^k − 1. One unit target
t₀ with c_top(t₀) ≠ 0 certifies the second case.

*Proof.* Group the sum by orbits: N² c_top(t) = Σ_O |O| w_O / q_O(t), with |O| ∈ {1, 2},
nonzero mod p since p is odd. Put Q = Π_O q_O and P = Σ_O |O| w_O Π_{O′≠O} q_{O′}. Then
N² c_top = P/Q on U, and Q ≠ 0 on U, so E = U ∩ {P = 0}.

Degrees: there are N diagonal orbits, each of t-degree ≤ 1, and N(N−1)/2 off-diagonal orbits, each
of t-degree ≤ 2. So deg Q ≤ N + N(N−1) = N², and deg P ≤ deg Q − min_O deg q_O. If some diagonal
q_O is constant, then deg Q ≤ N² − 1 and deg P ≤ deg Q ≤ N² − 1. Otherwise min_O deg q_O ≤ 1
(a diagonal orbit), so again deg P ≤ N² − 1.

If P ≡ 0, then c_top ≡ 0 on U. If not, P has at most N² − 1 roots. A unit t₀ with
c_top(t₀) ≠ 0 has P(t₀) ≠ 0, so P ≢ 0. ∎

**General m.** Let ĝ be the normal form of S_{m+1}(x₁, …, x_m, t), of t-degree d_t (d_t ≤ 2^{m−1}),
symmetric under S_m. The same argument over the S_m-orbits of H^m (multisets, C(N+m−1, m) of
them) gives |E| ≤ d_t·C(N+m−1, m) in the case P ≢ 0, and |E| ≤ d_t N^m − 1 without using the
symmetry. For m = 3, d_t ≤ 4.

*Proportion.* For fixed k (and m) and p → ∞ over primes p ≡ 1 mod N, the proportion is
|E|/|U| = O(4^k/p) = O(1/p). Here |U| ≥ p − 2N² because each q_α ≢ 0 has at most two roots.
In the repository regime, p = c·2^k with small c, the bound N² − 1 exceeds p and says
**nothing**. The preregistered P2 (proportion ≥ 1 − 3/p) is therefore **not** a consequence of
this note. It is empirical: lane A, and `tools/exceptional_count.out`.

| k | p | unit t (all curves) | exceptional t | fraction | 1/p | max per curve / N² − 1 |
|---|---|---|---|---|---|---|
| 1 | 101 | 590 | 10 | 0.0169 | 0.0099 | 3 / 3 (attained) |
| 2 | 97 | 502 | 3 | 0.0060 | 0.0103 | 1 / 15 |
| 2 | 113 | 411 | 2 | 0.0049 | 0.0088 | 2 / 15 |
| 3 | 97 | 199 | 1 | 0.0050 | 0.0103 | 1 / 63 |
| 3 | 73 | 3060 | 40 | 0.0131 | 0.0137 | lane A |
| 3 | 41..401 (9 primes) | 18,492 | 129 | about 1/p | | lane A, all brute-forced |
| 4 | 97..2017 | 282 | 2 | about 1/p | | lane A |
| 5..13 | 193 ≤ p < 2^31 | 11,000+ | 0 | 0 | ≤ 0.005 | lane A (1/p too small to see) |

The c_top ≢ 0 certificate was present on every curve tested. A proof that c_top ≢ 0 for every
nonsingular (a, b) is open (§8). [pcc:unverified]

### 4.4 Static versus iterated: the gap

**Corollary 4.4 (generic gap, m = 2).** For p ≥ 5 with N | p − 1, a nonsingular curve, and a unit
instance with c_top ≠ 0:

    d_lf = sd^CG ≤ k + 2 < 2k + 2 = sd^stat,     sd^stat − d_lf ≥ k   (k ≥ 1).

*Proof.* sd^stat = 2k + 2 is Corollary 4.1 with δ = 2. d_lf ≤ k + 2 is U1 Theorem 1 / [pcc:unverified]
Corollary 1 (ml-cryptanalysis `docs/proofs/U1-upper-bounds.md`, **proof**, refereed twice
without a gap; not re-proved here). ∎

The other recorded bounds for comparison (none re-proved here): [pcc:unverified]

| bound on d_lf (m = 2 unless stated) | label in its source | source |
|---|---|---|
| d_lf ≤ mk + 1 (K1, every m) | proof | ml-cryptanalysis `docs/IC_BOUNDS_PROGRAM.md` K1; `docs/IC_FALL_PROOFS.md` |
| d_lf ≤ k + 2 (U1) | proof (k + 2 attained at k = 1, 2: exact) | `docs/proofs/U1-upper-bounds.md` T1, T1' |
| d_lf ≥ ⌊(k − 10)/14⌋ + 1 for k ≥ 38 (char 0 or p > 2^{k+1}, t ≠ 0) | proof given T1 of N8 | `docs/proofs/N8-growth.md` §0, C1 |
| d_lf = ⌊(k + 14)/4⌋ generically | conjecture (measured) | U1 C1 |
| sd^stat = 2k + 2 generically; exact when c_top ≠ 0 | proof (this note) | Corollary 4.1 |

Under the measured law, sd^stat − d_lf = 2k + 2 − ⌊(k+14)/4⌋, about 7k/4 for large k [pcc:unverified]
(**conjecture**). On exceptional instances only d_lf ≤ sd^stat ≤ 2k + 1 is proved, and the gap
can be smaller (REDTEAM A8, A10).

| k | fraction sd^stat = 2k+2 | exact d_lf | sd^stat − d_lf |
|---|---|---|---|
| 3 | 1 − 40/3060 at p = 73; 1 at p ≥ 193 | 4 | 4 (3 on exceptions) |
| 4 | 1 − 2/282; 1 at p ≥ 193 | 4 | 6 (5 on exceptions) |
| 5 | 1 (48 + 11,000 scanned) | 4 | 8 |
| 6 | 1 | 5 | 9 |
| 7 | 1 | 5 | 11 |
| 8 | 1 | 5 | 13 |
| 9..13 | 1 | 5, 6, 6, 6, 6 (law, not measured) | 15, 16, 18, 20, 22 |
| m = 3, k = 2 | 1 at p ≥ 29; 156/163 at p = 13 | 7 (measured) | 5 |
| m = 3, k = 3 | 1 at p ≥ 137 | 7 (measured) | 8 |
| m = 3, k = 4, 5 | 1 at p ≥ 1153 | 7 (recorded) | 11, 14 |
| m = 3, k = 6..8 | 1 | not measured (≥ 8 recorded at k = 6, one proper instance) | ≤ 16 at k = 6 |

All rows [pcc:unverified] (lanes A and B, `experiments/ic-bounds-solving/`).

## 5. Proper ideals (Z ≠ ∅)

Let I = gA, a proper ideal: g vanishes on Z ⊂ H^m with Z ≠ ∅. Then ⟨F⟩ = NF^{-1}(I), and
I = {a ∈ A : a|_Z = 0} when p ∤ N and H ⊂ K. Lemma 2.1 still gives NF(Mac_D) = g·M_{≤D−δ}. But
1 ∉ ⟨F⟩, and the right target is now a Gröbner basis or all of I ∩ R_{≤D}.

**Proposition 5.1 (transfer of the full-containment degree).** For every D ≥ 0,

    Mac_D ⊇ ⟨F⟩ ∩ R_{≤D}   ⇔   g·M_{≤D−δ} ⊇ I ∩ M_{≤D}.

Both sides hold for every D ≥ δ + mk. Define d_full^stat as the least D₀ such that the
containment holds for every D ≥ D₀. Then d_full^stat ≤ δ + mk.

*Proof.* ⇒: Let a ∈ I ∩ M_{≤D}. Then ι(a) ∈ ⟨F⟩ has degree wt(a) ≤ D, so ι(a) ∈ Mac_D, and
a = NF(ι(a)) ∈ g M_{≤D−δ}. ⇐: Let P ∈ ⟨F⟩ ∩ R_{≤D}. Then NF(P) ∈ I, and wt NF(P) ≤ deg P ≤ D by
Lemma 2.1(b), so NF(P) = g h with h ∈ M_{≤D−δ}. Now ι(h) ĝ ∈ Mac_D, and P − ι(h) ĝ ∈ ker NF ∩ R_{≤D}
⊆ Mac_D by Lemma 2.1(a). For D − δ ≥ mk, M_{≤D−δ} = A, so g M_{≤D−δ} = I. ∎

The protocol's d_full also asks that the monomials of degree < D map onto A/I (functions on Z).
In A that condition is M_{≤D−1} + I = A, which holds for D ≥ mk + 1.

**Proposition 5.2 (static Gröbner degree).** Fix DRL on R and let 𝒢 be the reduced Gröbner basis
of ⟨F⟩. For G ∈ 𝒢 with NF(G) ≠ 0 put w(G) := min{wt(h) : g h = NF(G)}; the minimum runs over the
coset h₀ + Ann(g), and Ann(g) = functions supported on Z. Put δ + w(G) := −∞ if NF(G) = 0. Let
sd^stat_GB be the least D such that the reduced row echelon form of Mac_D (columns in DRL order)
contains 𝒢. Then

    sd^stat_GB = max_{G∈𝒢} max(deg G, δ + w(G)) ≤ max(2, δ + mk).

*Proof.* Let W ⊆ ⟨F⟩ be a subspace and G ∈ 𝒢 ∩ W with leading monomial L. The other monomials of
G are standard, so none of them is a leading monomial of an element of W ⊆ ⟨F⟩. Hence G is the
RREF row of W with pivot L. Conversely, the RREF of W contains G only if G ∈ W. So the condition
is 𝒢 ⊆ Mac_D. For G ∈ 𝒢 and deg G ≤ D, the argument of Proposition 5.1 gives: G ∈ Mac_D iff
NF(G) ∈ g M_{≤D−δ}, that is, iff NF(G) = 0 or w(G) ≤ D − δ. If deg G > D, then G ∉ Mac_D ⊆ R_{≤D}.

Bound: w(G) ≤ mk. A leading monomial of 𝒢 is minimal non-standard, so it is either u_{i,j}² or
squarefree of degree ≤ mk. In DRL the tail has degree ≤ the leading degree. ∎

Elements G ∈ 𝒢 with NF(G) = 0 are exactly the elements of 𝒢 that lie in the chain ideal; for
example c_{i,j} itself, when u_{i,j}² is a minimal leading monomial and u_{i,j+1} is standard.
For a unit ideal, 𝒢 = {1} and Proposition 5.2 gives back Theorem 3.1, with w(1) = maxwt g^{-1}.

**Corollary 5.3 (the iterated solving degree for proper ideals).** For DRL,
sd^CG = max(d_lf, max.GB.deg) (CG Theorem 3.1, cited), and

    max.GB.deg ≤ max(2, min(σ + 1, mk)) ≤ max(2, mk),

where σ is the largest degree of a DRL standard monomial of ⟨F⟩. In particular
sd^CG ≤ max(d_lf, 2, mk).

*Proof of the GB degree bound.* A leading monomial L of a reduced GB element is a minimal
non-standard monomial. Either L = u_{i,j}², or L is squarefree, hence of degree ≤ mk, and
L = u·L′ with L′ standard, hence deg L ≤ σ + 1. In a degree-compatible order the tail has degree
≤ deg L. ∎

If the monomials of degree < D* map onto A/I (the binary programme's certificate condition, §6),
then σ ≤ D* − 1: images of monomials of degree ≤ d span the same space as the standard monomials
of degree ≤ d, and standard monomials are independent mod I. Hence max.GB.deg ≤ max(2, D*).

## 6. The binary case

Let B = F₂[x₁, …, x_n]/(x_i² + x_i) (the Boolean ring, functions F₂^n → F₂), with descended
generators f₁, …, f_s of degrees δ_s, and let F = {x_i² + x_i} ∪ {f_l} in R = F₂[x].

**What survives (proof).**

* *Lemma (a).* The field equations have pairwise coprime leading monomials x_i², so they form a
  Gröbner basis. Division is degree-nonincreasing, with reducer multiples of degree ≤ deg P. So
  ker(R → B) ∩ R_{≤D} = span of field-equation multiples of degree ≤ D ⊆ Mac_D. The proof of §2
  goes through verbatim.
* *Lemma (b), several generators.* The normal form of a monomial is the squarefree monomial on
  its support, of degree ≤ deg μ, and every squarefree monomial of degree ≤ d is reached. So
  NF(Mac_D) = Σ_l NF(f_l)·M^B_{≤D−δ_l}, with M^B_{≤d} the squarefree monomials of degree ≤ d.
  Proposition 5.1 transfers in the same way, with A replaced by B and mk by n.
* *Upper bound.* M^B_{≤n} = B, so for a unit ideal sd^stat ≤ max_l δ_l + n.

**What does not survive.**

* The chains u² − u′ and the field equations x² − x are both degree-lowering. The difference is
  in the quotient. In A, every variable u_{i,j} is a **unit** (x_i^{2^j} ∈ G), A is a group
  algebra, and one generator can generate the unit ideal. In B, every element is idempotent
  (f² = f), so **the only unit is 1**. A single generator gives the unit ideal only if it is 1,
  and there is no inverse formula. (B is not a group algebra over F₂ either:
  F₂[Z/2] = F₂[y]/(y+1)² is local, not Boolean.)
* With several generators, the unit condition 1 ∈ Σ_l NF(f_l) M^B_{≤D−δ_l} is a degree-bounded
  Nullstellensatz certificate, a linear system with no closed form. The ideal is principal,
  ⟨f₁, …, f_s⟩ = (1 + Π_l(1 + f_l))B, but the Macaulay space is not the one of that single
  generator: its degree is up to Σ δ_l.

**Relation to the programme's certificate degree D*.** The binary programme
(ml-cryptanalysis `experiments/ic-bounds-binary/README.md`, "Certificates") certifies a proper
system at the first D* such that M_{<D*} maps onto the functions on Z and codim V_{D*} = |Z|. The
argument at the end of §5 gives max.GB.deg ≤ max(2, D*) in DRL. With CG Theorem 3.1,

    LFD ≤ sd^CG = max(LFD, max.GB.deg) ≤ max(LFD, 2, D*).

On the recorded rows D* ≥ LFD, so sd^CG ∈ [LFD, D*]. The recorded D* − LFD is 0 in 1476, 1 in 8
and 3 in 2 of 1486 proper systems (`experiments/ic-bounds-binary-law/README.md` line 93; not
recomputed here, [pcc:unverified]). Hence **sd^CG − LFD ∈ {0, 1, 3}** at most, row by row. The
static d_full^stat on the same rows (lane C, `experiments/ic-bounds-solving/laneC_binary/`): d_full^stat = k + 1 or k + 2 on all 185 proper systems with an exact value (143 at k + 1, 42 at k + 2; k = 2..8), lower bounds ≥ 8 at k = 9 and ≥ 7 at k = 10, never below LFD, and above the Bardet semi-regular degree by 4 to 5 at k = 8; LFD on the same systems stays in {3, 4, 5} [pcc:unverified].

## 7. What this means for the Gröbner route in index calculus

[pcc:unverified] for the cost remarks; the degree statements are the proofs above.

* An algorithm that builds **one** Macaulay matrix in degree D from the original generators,
  with no re-multiplication of new low-degree rows (plain XL, the protocol's sd^stat), must
  reach **D = 2k + 2** on a generic two-summand unit instance (δ + 3k for three summands). It
  cannot stop earlier. The certificate of 1 is ι(g^{-1})·ĝ, and generically g^{-1} has full
  support, top monomial included. In the ordinary ring that matrix has
  C(2k + 2k + 2, 2k) columns at D = 2k + 2, compared with C(2k + k + 2, 2k) at D = k + 2.
  [pcc:unverified]
* An algorithm that re-multiplies the rows it finds (MutantXL, F4/F5 in the degree-driven
  normal-strategy sense, the iterated closure V_{F,D}) stops at d_lf = sd^CG ≤ **k + 2**
  (U1, proof), and at about ⌊(k + 14)/4⌋ on the measured instances (conjecture). The mutant step
  is worth at least k degrees on every generic instance (Corollary 4.4).
* The first fall degree d_ff and degree-of-regularity heuristics are not a substitute for either
  invariant. HKYY warn against replacing LFD by FFD, and the exact counterexample with
  d_ff = 4 < d_lf = 5 is in IC_FALL_PROOFS §7. The static degree is pinned by the single
  algebraic quantity maxwt supp g^{-1}, which a semi-regularity heuristic in the style of
  Bardet–Faugère–Salvy does not see: the system is never semi-regular in the ordinary ring,
  because the chains have a degree fall at degree 3. [pcc:unverified]
* None of this is a speedup claim. The relation search still needs a solve per target, and the
  comparison with rho is AGENTS.md §5a's job.

## 8. Open problems

1. **c_top ≢ 0 for every nonsingular curve.** Prove that for every (a, b) with 4a³ + 27b² ≠ 0, and [pcc:unverified]
   every k, the rational function c_top(t) is not identically zero. Candidate: its Laurent
   expansion at t = ∞. The diagonal terms α₁ = α₂ = α contribute
   −(4N²)^{-1} Σ_{α∈H} α²/(α³ + aα + b) to the 1/t coefficient (when no α³ + aα + b vanishes);
   the off-diagonal terms are O(1/t²). So a sufficient condition is
   Σ_{α∈H} α²/(α³ + aα + b) ≠ 0. [pcc:unverified; derivation sketched, not checked.]
2. **A bound on the exceptional proportion in the regime p ≈ c·2^k**, where the root count is vacuous:
   for example an equidistribution or Weil-type bound for the character sum c_top(t) as t
   varies.
3. A partial-fraction or residue evaluation of the inner sum Σ_{β∈H} β/q(β). It reduces
   c_top to N terms over F_{p²}, in terms of the x-coordinates of P₁ ± T (T the target point).
   This would give a closed criterion for sd^stat < 2k + 2.
4. Lean formalization of Lemma 2.1 and Theorem 3.1. This should be short on top of
   `ChainTranslation` and `ChainClosureBridge`.
5. The static three-summand statement with the exact δ_3 and exceptional count (lane B data).
6. For proper ideals, the exact distribution of w(G_j) relative to mk, and whether
   d_full^stat = δ + mk generically (lane C for binary; untested for prime).

## References

* M. Caminata, E. Gorla, *Solving degree, last fall degree, and related invariants*,
  arXiv:2112.05579v2: Section 1 (mutant Macaulay construction), Definitions 1.1 and 1.5, eq. (1),
  Theorem 3.1. **Quoted by the coordinator from the PDF; not re-read by lane D.**
* M.-D. Huang, M. Kosters, Y. Yang, S. L. Yeo, *On the last fall degree of zero-dimensional Weil
  descent systems*, arXiv:1505.02532: Definitions 2.1 and 2.5, Proposition 2.6. **Not re-read
  here**; cited through IC_FALL_PROOFS §1.
* I. Semaev, *Summation polynomials and the discrete logarithm problem on elliptic curves*,
  IACR ePrint 2004/031: definition of S_3, S_{m+1}. **Not re-read here**; the S_3 formula used in
  the scripts is the standard one, checked only through the exact runs.
* C. Petit, J.-J. Quisquater, *On polynomial systems arising from a Weil descent*, ASIACRYPT 2012:
  first-fall-degree heuristics for binary descent. **Not re-read here.**
* M. Bardet, J.-C. Faugère, B. Salvy, *On the complexity of Gröbner basis computation of
  semi-regular overdetermined algebraic equations* (ICPSS 2004) and related work on the degree of
  regularity of semi-regular sequences. **Not re-read here.**
* ml-cryptanalysis: `docs/proofs/U1-upper-bounds.md` (§0–2 read), `docs/IC_FALL_PROOFS.md` (§1–2
  read), `docs/proofs/N8-growth.md` (§0 read), `docs/IC_BOUNDS_PROGRAM.md` (K1, line 30),
  `experiments/ic-bounds-binary/README.md`, `experiments/ic-bounds-binary-law/README.md`.

## Files

* `experiments/ic-bounds-solving/laneD_theory/SD1-solving-degree.md`: this note.
* `experiments/ic-bounds-solving/laneD_theory/REDTEAM.md`: refutation attempts A0–A12, the exact check tables, sources checked.
* `experiments/ic-bounds-solving/laneD_theory/tools/check_small.py`, `tools/check_small.out`: ordinary-ring static Macaulay brute
  force against δ + maxwt supp g^{-1}, with the character-sum check (174 cases, 0 mismatches).
* `experiments/ic-bounds-solving/laneD_theory/tools/exceptional_count.py`, `tools/exceptional_count.out`: exceptional targets
  against the N² − 1 bound.
