# Lane D red team: attempts to refute the Lemma and Theorem of PROTOCOL.md

Study: `research/solving_degree_20261115/` (protocol frozen 2026-11-15T02:05Z). Lane D, theory.
Date 2026-11-15. Everything here is [pcc:unverified] until an independent validator runs it.
No pcc claim id has been issued for any statement below.

**Verdict.** The Lemma and the Theorem **hold**. No gap was found in the mathematics. One issue
with how the result was named was found and fixed (A0). One hypothesis has to be stated
explicitly (A3: δ is the ordinary degree of ĝ, not the weight of g). Three of the Corollaries need
qualifiers (A9, A10, A11). P6 survives in a stronger form than preregistered: the Lemma and the
Theorem need no roots of unity and no condition on the characteristic (A6).

Small exact check: `tools/check_small.py` → `tools/check_small.out`.
**174 instances, 68 unit, 106 non-unit, 0 mismatches** (section B).

---

## A. Attacks

Each entry gives the step attacked, the attack, and the outcome.

### A0. Attribution: "sd = Caminata–Gorla solving degree" (the naming was refuted and corrected)

*Step.* PROTOCOL.md calls "least D with 1 ∈ Mac_D" the "static solving degree (Caminata–Gorla
for a unit ideal)".

*Attack (from the coordinator, 2026-11-15, and checked against the text the coordinator quoted).*
Caminata–Gorla, arXiv:2112.05579, Section 1, build M_d by Gaussian elimination and then add
rows u·f_ℓ for the new lower-degree rows f_ℓ, repeating until no new rows appear (a mutant
step). Definition 1.1: sd_σ(F) is the least d such that rowsp(M_d) contains a Gröbner basis.
Equation (1): V_{F,d} = rowsp(M_d) for degree-compatible σ. Theorem 3.1:
sd_σ(F) = max{d_F, max.GB.deg_σ(F)}. For a unit ideal the reduced Gröbner basis is {1}, so
**sd^CG = d_lf**. The coordinator checked this exactly at m = 2, k = 1, p = 101. On 8/8 random
unit instances the static least D with 1 ∈ Mac_D was 4 = 2 + maxwt supp g^{-1}, and the iterated
d_lf was 3. So the static quantity is not the Caminata–Gorla solving degree.

*Outcome.* **The attribution is refuted. The mathematics is unaffected.** SD1 renames the static
quantity **sd^stat** (the XL / no-mutant Macaulay degree) and states the Theorem for sd^stat. It
lists d_ff, d_lf = sd^CG (unit case) and sd^stat side by side, and cites CG Theorem 3.1 as
re-read by the coordinator. Lane D did not fetch the PDF itself; the quotes come from the
coordinator's message. My own k = 1 rows agree with the coordinator: for the normalized S_3
instances, sd^stat = 4 = 2 + 2k at k = 1 (section B). The exact d_lf = 3 at k = 1 is the U1
attainment T1' (k + 2 = 3).

### A1. Lemma (a): chain multiples cover the kernel within the degree budget

*Step.* ker NF ∩ R_{≤D} ⊆ span{m·c_{i,j} : deg(m c_{i,j}) ≤ D}.

*Attack.* Look for a P of degree ≤ D with NF(P) = 0 whose every chain representation needs a
multiple of degree > D. That would be a "degree fall inside the chain ideal".

*Outcome.* **No such P exists.** The leading monomials u_{i,j}² of the chains, in any
degree-compatible order with u_{i,j}² > u_{i,j+1}, are pairwise coprime. Buchberger's first
criterion then makes the chains a Gröbner basis. Division of P by the chains rewrites a term
c·u_{i,j}²·m′ using the multiple c·m′·c_{i,j}, whose degree equals deg(u_{i,j}² m′) ≤ deg P. The
remainder is the unique reduced (squarefree in each u_{i,j}) representative, so it is 0 when
NF(P) = 0. This holds over any commutative ring, because the leading coefficients are 1.

### A2. Lemma (b): the reachable set is exactly {wt e ≤ d}

*Step.* {NF(m) : m an ordinary monomial, deg m ≤ d} = {x^e : wt(e) ≤ d}.

*Attacks.* (i) A carry could raise the weight by more than one. (ii) The wrap-around
u_{i,k−1}·x^{e} with bit k−1 set (2^{k−1} + 2^{k−1} ≡ 0 mod N) could make NF(m) a non-monomial or
introduce a coefficient other than 1. (iii) Some e with wt(e) ≤ d might not be reachable.

*Outcome.* **All three attacks fail.** (i) Adding 2^j to a k-bit word clears a run of r ≥ 0 ones
and sets at most one zero bit, so the weight changes by 1 − r ≤ 1. A carry out of bit k−1 is
dropped mod N, and the weight then falls by r. (ii) Every chain rewrite maps a monomial to a
monomial with coefficient +1, because the chains are binomials u² − u′ and u² − 1. So NF(m) is a
monomial x^e. (iii) The reduced monomial Π u_{i,j}^{bit_j(e_i)} has degree wt(e) and NF x^e.

### A3. Theorem: δ versus wt(g) (a hypothesis that must be explicit)

*Step.* deg(m·ĝ) = deg m + δ is used to convert "deg(m ĝ) ≤ D" into "deg m ≤ D − δ".

*Attack.* If ĝ is not reduced, deg ĝ can exceed wt(g). For the raw S_3 (degree 4), NF has weight 2.
If the formula used wt(g), it would predict sd = 2 + maxwt and be wrong.

*Outcome.* **The Theorem holds with δ := deg ĝ, the ordinary total degree. It is false with
δ := wt(g) when ĝ is not reduced.** Exact check: on the raw S_3 rows, sd^stat = 4 + maxwt supp g^{-1}
(for example m = 2, k = 3, p = 41, a = 2, b = 3, t = 12 gives sd^stat = 10 = 4 + 6), whereas
2 + maxwt would predict 8. SD1 states δ = deg ĝ and notes that δ = wt g for the normalized
(reduced) generator. Products in R do not drop degree because R is a domain, so
deg(m ĝ) = deg m + deg ĝ holds exactly.

### A4. Theorem: the lift step "1 ∈ g·M_{≤D−δ} ⇒ 1 ∈ Mac_D"

*Step.* The proof lifts h = g^{-1} to the reduced polynomial ĥ, forms P = ĥ·ĝ ∈ Mac_D, and puts
P − 1 into the chain span.

*Attack.* deg(P − 1) might exceed D, or P − 1 might need chain multiples of degree > D.

*Outcome.* **The step holds.** deg ĥ = wt(h) ≤ D − δ, so deg P ≤ D, and deg(P − 1) ≤ max(D, 0).
NF(P − 1) = g h − 1 = 0, so A1 puts P − 1 into the chain multiples of degree ≤ D. For D < 2
the chain span is empty in degree ≤ D, but then P − 1 is reduced with NF 0, hence P − 1 = 0.

### A5. Theorem: the converse "1 ∈ Mac_D ⇒ g^{-1} ∈ M_{≤D−δ}"

*Step.* 1 ∈ Mac_D ⇒ NF(1) = 1 ∈ NF(Mac_D) = g·M_{≤D−δ}. So g h = 1 with h ∈ M_{≤D−δ}, and h = g^{-1}
by uniqueness of inverses.

*Attack.* A non-unit g could have several h with g h = 1. Alternatively, Mac_D could contain 1
through a chain-only combination.

*Outcome.* **The step holds.** If g h = 1 then g is a unit and h is unique. A chain-only
combination has NF 0 ≠ 1.

### A6. P6: field and characteristic hypotheses

*Attack.* Does any step need 2^k | p − 1, p odd, or H ⊂ K?

*Outcome.* **No step does, for the Lemma and the Theorem.** They hold over any field K, including
characteristic 2 and fields without the N-th roots of unity: the proofs use only the Gröbner
property of the chains, the weight inequality and uniqueness of inverses in the commutative ring
K[G]. The hypotheses enter only in three places:
(i) "g unit ⇔ g has no zero on H^m" needs H ⊂ K. Over K̄ it is the Nullstellensatz
statement V(F) = ∅, which holds in general.
(ii) The character-sum formula needs p ∤ N and H ⊂ K. In characteristic 2, K[G] is local and
g is a unit iff g(1, …, 1) ≠ 0.
(iii) The O(1/p) count needs p odd, so that orbit sizes 1 and 2 are nonzero.
**P6 survives**, more strongly than preregistered.

### A7. Generic sd^stat = δ + mk: might the top coefficient always vanish by symmetry?

*Attack.* Look for a symmetry of S_3 (x₁ ↔ x₂, or α ↦ −α) that forces
c_top = N^{-2} Σ α₁α₂ / S_3(α₁, α₂, t) to vanish identically.

*Outcome.* **No such symmetry exists.** The x₁ ↔ x₂ swap fixes the summand. The map α ↦ −α does [pcc:unverified]
not preserve S_3 (its odd-degree terms change sign) unless a special relation holds. Exact data:
c_top ≠ 0 on most unit rows (section B). In `tools/exceptional_count.out`, the per-curve count of
targets t with c_top = 0 was 0 to 3 out of 63 to 105 unit targets at p = 97, 101, 113. The
character sum was checked equal to the top coefficient of the linearly solved inverse on all
68 unit rows.

### A8. Low-weight inverses: is sd^stat ≥ some function of k always?

*Attack.* Choose g with g² = 1, for example ĝ = u_{1,k−1} (g = x₁^{N/2}).

*Outcome.* **There is no lower bound of the form sd^stat ≥ δ + f(k) without genericity.** This is
consistent with the Theorem: sd^stat = 2 = 1 + 1 at every k (rows `ghat=u_{1,k-1}`, all
matching). Consequence for the corollaries: the gap sd^stat − d_lf ≥ k is a statement about
**generic** instances (c_top ≠ 0), not about all unit instances. SD1 states it that way.

### A9. Corollary "exceptional set has proportion O(1/p)"

*Attack.* (i) The rational function c_top(t) might vanish identically, making every unit t
exceptional. (ii) The degree bound might be too large to mean anything.

*Outcome.* **The statement is true but must be qualified.** (i) The bound is a dichotomy: either
c_top ≡ 0 for that (a, b), or #exceptional ≤ N² − 1 (m = 2). A single non-exceptional unit target
certifies c_top ≢ 0. Every curve in `tools/exceptional_count.out` has one. A proof that c_top ≢ 0
for every nonsingular (a, b) is **open**: the expansion at t = ∞ gives a candidate certificate,
but it was not completed. (ii) The bound N² − 1 = 4^k − 1 is O(1/p) only for fixed k as p → ∞.
In the repository regime p ≈ c·2^k with c < 2^k it exceeds p and is **vacuous**. P2's
1 − 3/p is therefore **not** implied by a proof here. It remains empirical (lane A). The bound
is attained: at k = 1, p = 101, two of six curves have exactly 3 = N² − 1 exceptional targets.

### A10. Corollary "sd^stat − d_lf ≥ k"

*Attack.* The corollary combines U1 (d_lf ≤ k + 2, proof) with sd^stat = 2k + 2. The second holds
only when c_top ≠ 0.

*Outcome.* **It holds on the generic set only.** sd^stat − d_lf ≥ (2k + 2) − (k + 2) = k when
c_top ≠ 0. On exceptional instances sd^stat ≤ 2k + 1 and only sd^stat ≥ d_lf survives. One
exceptional row appears in the small check (m = 2, k = 2, p = 13: sd^stat < 6).

### A11. "D* − LFD ∈ {0, 1}" (coordinator's paraphrase of the binary README)

*Attack.* Read `ml-cryptanalysis/experiments/ic-bounds-binary-law/README.md` line 93.

*Outcome.* **The paraphrase is slightly wrong.** The README says D* − LFD is 0 in 1476 systems,
1 in 8, and **3 in 2** of 1486 proper systems. SD1 uses {0, 1, 3} and labels the counts as
recorded data that were not recomputed here.

### A12. Proper ideals: is "Mac_D ⊇ ⟨F⟩ ∩ R_{≤D}" monotone in D?

*Outcome.* **It is not shown to be monotone.** SD1 therefore defines d_full^stat as the least D
from which the containment holds for every D′ ≥ D. The transfer to A is exact at every D, and the
containment holds for every D ≥ δ + mk.

---

## B. Small exact check (ordinary ring, brute force, own code)

Command (single thread, `ulimit -v 6000000`, 49 s):

    python3 tools/check_small.py > tools/check_small.out

The script builds the ordinary-ring Macaulay rows m·f (f a chain or ĝ, deg ≤ D) over F_p in
F_p[u_{i,j}] by sparse elimination with pivots on the highest monomial. It reports the least D
with the constant monomial as a pivot, and compares it with deg ĝ + maxwt supp g^{-1}. The
inverse g^{-1} is computed twice, by a linear solve with the multiplication matrix of F_p[G] and
by Fourier inversion over H^m; the script asserts that the two agree. It also checks
c_top = N^{-m} Σ_α Π α_i / g(α) and "sd^stat = δ + mk ⇔ c_top ≠ 0". No repository engine is used.

| (m, k, p) | generators | unit / total | mismatches |
|---|---|---|---|
| (2,1,3), (2,1,7) | normalized S_3, raw S_3, random reduced δ ∈ {1,2}, u_{1,0} | see `.out` | 0 |
| (2,2,5), (2,2,13), (2,2,17) | normalized S_3, raw S_3, random reduced δ ∈ {1,2,3}, u_{1,1} | see `.out` | 0 |
| (2,3,17), (2,3,41), (2,3,73) | normalized S_3, raw S_3 (p = 41), random reduced, u_{1,2} | 14 unit at k = 3 | 0 |
| (3,1,5) | random reduced, u_{1,0} | see `.out` | 0 |
| **all** | | **68 unit / 174** | **0** |

The non-unit rows (106) were checked too: on each of them, 1 ∉ Mac_D up to D = δ + mk + 1 and
neither inverse exists.

Normalized S_3 unit rows split by whether sd^stat = 2k + 2:

| k | p | sd^stat = 2k+2 | sd^stat < 2k+2 |
|---|---|---|---|
| 1 | 3 | 2 | 1 |
| 1 | 7 | 8 | 2 |
| 2 | 13 | 11 | 1 |
| 2 | 17 | 2 | 0 |
| 3 | 17, 41, 73 | 1 + 5 + 3 | 0 |

At k = 3 the prediction is the full brute force: sd^stat = 8 = 2 + 6 on all 9 normalized unit
instances, against a 3003-column ordinary Macaulay matrix in 6 variables at D = 8.

Exceptional-target count (`tools/exceptional_count.py` → `.out`, under 1 s), with c_top evaluated
by the character sum for every t ∈ F_p:

| k | p | curves | unit t | exceptional t | fraction | 1/p | per-curve max / bound N²−1 |
|---|---|---|---|---|---|---|---|
| 1 | 101 | 6 | 590 | 10 | 0.0169 | 0.0099 | 3 / 3 |
| 2 | 97 | 6 | 502 | 3 | 0.0060 | 0.0103 | 1 / 15 |
| 2 | 113 | 4 | 411 | 2 | 0.0049 | 0.0088 | 2 / 15 |
| 3 | 97 | 3 | 199 | 1 | 0.0050 | 0.0103 | 1 / 63 |

## C. Sources checked

* research-kb: **unavailable** in this lane. There is no MCP tool, and
  `cryptanalysis/scripts/research-kb status` fails with "no executable found". No source_id is
  recorded. The novelty of the inverse formula is **not assessed**.
* Read: `ml-cryptanalysis:docs/proofs/U1-upper-bounds.md` §0–2; `ml-cryptanalysis:docs/IC_FALL_PROOFS.md`
  §1–2; `ml-cryptanalysis:docs/proofs/N8-growth.md` §0; `ml-cryptanalysis:docs/IC_BOUNDS_PROGRAM.md`
  line 30 (K1); `ml-cryptanalysis:experiments/ic-bounds-binary/README.md` (certificates);
  `ml-cryptanalysis:experiments/ic-bounds-binary-law/README.md` line 93;
  `ml-cryptanalysis:formal/ICFallDegrees/ChainClosureBridge.lean:404` (name only).
* Caminata–Gorla arXiv:2112.05579: quoted by the coordinator from the PDF (Section 1,
  Definitions 1.1 and 1.5, eq. (1), Theorem 3.1). **Not fetched by lane D.**
