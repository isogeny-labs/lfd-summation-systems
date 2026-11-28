<!-- Agent report saved by the parent session from the workflow journal. Status: mixed. Independent skeptic verification: NOT YET RUN (failed on usage limit). -->

# Breakthrough red team: lead 1 (principal polarizations in the isogeny class of A_131)

**Bottom line.** At ell = 131 I found no error in the mathematics: Theorems A and B hold up under my own derivation and my own code. The committed artifact set is still not enough for a closure sign-off, because:
- the written proof is missing (only a sketch is committed);
- there is no Lean proof and no named Lean gap;
- the routine that produced the small-ell controls has a false-negative bug.

These are review findings. I am not marking anything supported or closed.

Reviewer slug: `redteam-lead1`. This is a separate session from the lead-1 producer and from the coordinator.

## 1. What was reviewed

- **Lead-1 files.** `experiments/ecc2k130-weil-factor/lead1_polarization/*` at commit `6450ea9`. They are byte-identical at HEAD `1159028` and at the WIP snapshot `87fae43` (`git diff` is empty).
  - `REPORT.md`: sha256 `6a4982dd…`
  - `decide_ppav.gp`: `dababaab…`
  - `run_17_orders.out`: `89b51566…`
- **What lead 1 builds on.**
  - Lead 2, `lead2_index/PROOF.md` (commit `688c805`, sha256 `a17b8c30…`).
  - The PR #57 certificate, `explicit_ppav.json` (sha256 `6b91b86d…`).
- **Validation record.** No `VALIDATION.json` applies to a mathematical claim. There is no Lean proof, and the container has no Lean toolchain.

**Is the artifact set enough at this tier?**
- To review the mathematics: yes, but only because I re-derived the argument and recomputed the controls myself.
- To close the question: no. `REPORT.md` contains a five-line proof sketch, and `decide_ppav.gp` and `run_17_orders.gp` cite "REPORT.md sect. 2" and "sect. 4.3", which do not exist. The parent session says it checked "consistency, not every line".

## 2. Verdict per claim

| Claim | Verdict | Basis |
|---|---|---|
| The CM type depends only on the embedding of sqrt-7 | **upheld** (proved) | Lemma L1. Checked for every Galois element at ell = 3–29, and 2-adically at 131 |
| Positivity reduces to r being totally positive or totally negative | **upheld** (proved) | Lemma L2. Running both signs is redundant (lambda -> -lambda swaps them) |
| Principal means lambda * conj(T) = T^dual | **upheld** (proved) | Lemma L3 |
| Hilbert-reciprocity bookkeeping | **upheld** (proved); two omissions in the write-up | L4–L6. Places above 7 and inert places where End is maximal are missing from the list; both give +1. Checked with `nfhilbert` place by place on explicit ppavs |
| Non-invertible modules | **upheld** (proved) | Every intermediate order is Gorenstein, so Bass. Lattice-verified at ell = 11 and 17 |
| Theorem A at 131: End non-maximal at v0, no ppav contains zeta, A not principally polarizable | **upheld**, no GRH needed | Section 3 |
| Theorem B at 131: every conjugation-stable order non-maximal at v0 occurs | **upheld** (proved via class field theory plus a local twist) | Section 3. PR #57 module re-verified |
| General rule for every odd ell ≠ 7 | **unsupported as stated** | Needs v_p(U_ell) = 1 for inert-type p |
| Recorded counts "4 of 5" (ell = 11), "43/45" (ell = 17), "146/148" (ell = 19) | **refuted**; the conclusions survive | GP zero-vector bug; see F1 |
| Pre-registered prediction at ell = 29 | **confirmed** (computed, GRH) | `rt29.out` |

## 3. The derivation I checked

Notation: K = Q(zeta_ell, sqrt-7), c is complex conjugation, K+ is its fixed field, and Gal(K/Q) = (Z/ell)^× × {±1} with c = (-1,-1). Let eps = (-7/ell). By Deligne and Howe (1995, §4), a ppav in the class is a pair (T, lambda): T a lattice in K stable under pi and pibar, lambda with conj(lambda) = -lambda, Im phi(lambda) > 0 on Phi, and Tr(lambda x ybar) unimodular on T.

- **L1 (CM type).** Phi = {phi : v_2(phi(pi)) > 0}, because Frobenius is topologically nilpotent on the Lie algebra of the canonical lift. Since phi(pi) = phi(alpha)phi(zeta), zeta is a unit, and the 2-adic roots of x^2+x+2 have valuations 1 and 0, Phi = {phi : phi(sqrt-7) = s0·sqrt-7}. So Phi is induced from Q(sqrt-7).
- **L2 (positivity).** Write lambda = sqrt-7 · r with r in K+. For phi in Phi_+, Im phi(lambda) = sqrt7 · phi(r). So positivity means r is totally positive (Phi_+) or totally negative (Phi_-).
- **L3 (principal condition).** The form is unimodular iff lambda · Tbar = T^dual.
- **L4 (K/K+ and its places).**
  - No inertia group (I_ell = (Z/ell)^× × 1, I_7 = 1 × {±1}, I_2 = 1) contains c, so K/K+ is unramified at every finite place.
  - The place v_ell is tamely and totally ramified, so v_ell(D_K+) = ell - 2, which is odd. It is inert in K/K+ iff eps = -1.
  - Places above 2 split in K/K+.
  - At places above 7: e = 2 and v(sqrt-7) = 1.
  - Local symbols: (r,-7)_v = (-1)^v(r) at inert v, and +1 at split v.
- **L5 (local structure, assuming v_p(U_ell) = 1).** Every p dividing U_ell is either of split type (p ≡ 1 mod ell and (-7/p) = 1) or of inert type (p ≡ -1 mod ell and (-7/p) = -1).
  - Inert type: with beta/alpha ≡ omega^j, the gluing is (i, a) ~ (i - j, b). This gives exactly one self-glued place v0 (i = j/2), one free place, and (ell - 3)/2 glued place-pairs.
  - The local components are Z_p + pZ_{p^2} at v0 (index p) and {(x, y) : x ≡ psi(y) mod p} at each glued pair (index p^2). The total index is p^(ell-2), which matches lead 2.
  - All component orders are monogenic, hence Gorenstein. So every order between Z[pi,pibar] and O_K is Bass, and equals Z[pi,pibar] + c for some conductor c.
- **L6 (local parity).** At a non-maximal component, R^dual = delta · R with v(delta) = -1 at each of its places. So v(r) = -1 - 2v(t) is odd at every inert place of a non-maximal component. At v_ell, v(r) = -(ell - 2) - 2v(t) is odd. At places above 7, v(r) = -2 - 2v(t) is even. At inert places where End is maximal, v(r) is even.
- **Theorem A.** The product of all local symbols is eps · (-1)^n(T), and Hilbert reciprocity forces it to be +1. At ell = 131:
  - eps = -1;
  - U_131 = 263 · p2, squarefree, with both factors proven prime;
  - 263 is of split type and p2 is of inert type.
  So n(T) = 1: End is non-maximal at v0. Hence zeta is not in End, since otherwise End would contain Z[pi, pibar, zeta] = O_K. In particular A itself has no principal polarization.
- **Theorem B.**
  - Necessity: R^dual · O_K = (D_K c)^-1 (lattice-verified), so a ppav gives lambda · S · Sbar = D^-1 c^-1, that is, r · N(S) = b+^-1 with b+ · O_K = sqrt-7 · D_K · c.
  - Sufficiency: the local unit condition can always be met. At v0 it is automatic. At a glued pair, twist the local generator by an element whose norm fixes a ratio in F_p^×. At a split pair, choose the generator at the conjugate component freely.
  - So a ppav exists iff the class of b+ lies in N(Cl K). By class field theory this is the index-2 kernel of the Artin map Cl+(K+) -> Gal(K/K+), and the Artin symbol of b+ is eps · (-1)^n.
- **Natural polarization.** lambda_nat = (2 - zeta - zeta^-1)/(ell · sqrt-7) satisfies lambda_nat · O_K = D^-1 · w_ell, so it has degree ell^2. Its valuation at w_ell is -(ell - 3), which is even, so all its local symbols are +1. A principal polarization on an O_K-module would have to remove w_ell. That flips the v_ell parity, and an O_K-module has no other inert place where End is non-maximal to flip it back. This is consistent with the degree-131^2 polarization.

## 4. Independent computations

All files are in `/home/user/lfd-summation-systems/experiments/ecc2k130-weil-factor/round2/redteam-lead1/`. My code shares nothing with `decide_ppav.gp`:
- K is built as K+(sqrt-7) with `rnfinit`, from a resultant (no `polcompositum`).
- **Method A:** class-by-class search over Cl(K), with positivity decided by the units of K+.
- **Method B:** narrow class group Cl+(K+) plus the image of the norm map. It does no unit or generator search in K.
- Every "exists" answer comes with an exact integer Gram determinant on an explicit lattice.

| Run | Result |
|---|---|
| `bugcheck_decide.gp/.out` (lead 1's `decide()` unchanged vs. a one-line fix) | ell = 11, End = O_K: 4 becomes 5 of 5, and all 5 have det 1. ell = 17: the antisymmetric counts 43 / 42 / 44 become 45. The unchanged code gives different counts on different runs |
| `rt_run_3_17.out`, `rt_run_19.out` | ell = 3, 5, 13: no ppav. ell = 11: 5/5 for O_K and for both split-pair orders, with det-1 lattices. ell = 17: O_K 0, v0 45/45, one pair 0, v0 + pair 45/45, two pairs 0. ell = 19: 0/148 for every order. Methods A and B agree everywhere, and \|W\| = \|ker N\| (coset check) |
| `rt_twist.out` (ell = 17) | All 256 orders checked by method B: the rule agrees in 256/256 (128 with ppav, 128 without, split exactly by n(v0)). v0 + 2 pairs: det-1 lattice; the symbol is -1 at v_17, at v0 and at the 4 pair places. Two pairs: W = 0. v0 + 3 pairs: W = 45, but the lattice search hit the 28-minute cap |
| `rt_gor.out` | Gorenstein, dual-hull, multiplier-ring and index checks all pass at ell = 11 and 17, including Z[pi,pibar] itself |
| `rt29.out` (ell = 29, eps = +1) | O_K has a ppav containing zeta (56×56 lattice, det 1). R(v0) and R(v0 + pair) have none. R(pair) has one |
| `rt131.py/.out/.json`, `rt131_w131.py/.json` | lambda_nat has Gram [[0, -G], [G, 0]], det 131^2. The PR #57 module, rebuilt from its recorded u, m and v: det on O_K is D^2, det on I is 1 (flint); stable under pi and pibar, not under zeta; positive at all 130 embeddings (700 digits, min Im about 8e-58). v_p2(m) = 1 only at v0, and v/D is non-integral only at v0. v_w131(lambda) = -129 versus -128 for lambda_nat |
| `rt131_int.out` | Integer facts at 131 confirmed. No prime of a third type and no square factor among the fully factored primes of U_ell for ell ≤ 113 |

## 5. Findings

Each finding gives its severity and the cheapest test that would overturn it.

- **F1 (high: evidence integrity).** `decide_ppav.gp` uses `if (sol == 0, next)` and `if (x == 0, next)`. In GP, `[0,..,0]~ == 0` is true, so whenever `matsolvemod` returns the zero solution the case is discarded.
  - The recorded counts 4/5, 43/45 and 146/148 are impossible: the set W of working classes is a coset of ker(N : Cl K -> Cl+ K+), so its size is 0 or a subgroup order.
  - Where h(K) = 1 (ell = 3, 5, 13), a single miss would have flipped the negative control.
  - The conclusions hold after correction.
  - Test: replace both tests with `type(.) == "t_INT"` and rerun ell = 11 (23 s). Done.
- **F2 (high: proof obligation).** No committed written proof, no Lean proof, no named Lean gap. Section 3 above is my own derivation; it is unreviewed and needs an independent review before it can stand as the written proof.
  - Test: falsify any one of L4–L6. Examples: the place over 131 splits in K/K+, v_131(D_K+) is even, or K/K+ is ramified above 7. Each was checked at ell ≤ 29, and the integer facts at 131.
- **F3 (medium: scope).** The general rule silently needs v_p(U_ell) = 1 for every inert-type p.
  - Test: find an ell with p^2 dividing U_ell for an inert-type p. Then an order non-maximal at v0 with even conductor exponent breaks the rule as stated.
- **F4 (low).** The bookkeeping list omits the places above 7 and the maximal inert places. Both contribute +1.
- **F5 (low).** "Both signs were run" is not an independent check.
- **F6 (low).** The heading "non-maximal at one place" is misleading: End may also be non-maximal at pairs above p2 or 263. Separately, the report table (45/45 and 43/45) and the main doc (45/45) disagree because of F1.
- **F7 (information, outside lead 1).** The main doc says Jacobian candidates are reached "only through isogenies of p2-power degree". Theorem A gives only that p2 divides the degree. The PR #57 ppav sits at index 263 · p2 over O_K. I did not review this further.

## 6. Lessons

I did not write these to the ledger, because this task restricts writes to the slug directory.

- **Candidate:** in GP, a zero vector compares equal to 0. Test "no solution" returns from `matsolvemod` with `type(x) == "t_INT"`. Evidence: `bugcheck_decide.out`.
- **Candidate:** a count of solutions over a class group must be 0 or the order of a subgroup. A count like 4/5 or 43/45 signals a bug.
- I ran `lessons.py brief` for `red-team-breakthrough`. Its adopted lessons are about PR and journal workflow; none applied to this review.

## 7. Not addressed

- Whether any ppav in the class is a Jacobian.
- Automorphisms defined only over extensions of F_2.
- A Lean formalization.

## Caveats (agent)

- **GRH.** Class groups and units at ell ≤ 29 use PARI's `bnfinit` under GRH; none were certified.
- **ell = 29 negatives.** The no-ppav results for R(v0) and R(v0+pair) also rely on the class-field-theory fact [Cl+(K+) : N Cl(K)] = 2. The subgroup generated by norms of primes up to 233 reaches that index.
- **ell = 131 positivity.** Positivity of the PR #57 certificate was checked by 700-digit floating-point evaluation, not interval arithmetic. PR #57 itself carries a rigorous enclosure.
- **ell = 17 with three pairs.** The explicit lattice for v0 + 3 glued pairs was not finished within the 28-minute cap. Method A found W = 45 for it, and method B agrees.
- **Superseded outputs.** I deleted the earlier outputs `rt_run_small.out` and the partial `rt_run_17_19.out`. They came from a version whose random Hilbert control drew large coefficients, and factoring a big norm stalled at ell = 19. The final drivers reproduce every result (`rt_run_3_17.out`, `rt_run_19.out`).
- **ell = 3 genus-2 check.** I did not rerun lead 1's genus-2 enumeration. I checked the Howe–Nart–Ritzenthaler criterion by hand instead: a = -1, b = -1, q = 2.
- **Proof status.** My derivation in section 3 is itself an agent proof and is unreviewed. It needs independent review and a Lean formalization or a named Lean gap; there is no Lean toolchain here.
- **Lessons ledger.** Nothing was written to the ledger, because the task restricts writes to the slug directory.
- **Commit pin.** The reviewed lead-1 files are pinned at commit 6450ea9, and are identical at 1159028 and 87fae43.
- **Out of scope.** Whether any ppav is a Jacobian; automorphisms defined only over extensions; lead 3's wording about isogeny degrees.
