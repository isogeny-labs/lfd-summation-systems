# Red-team review: `paper/curve-lfd-bounds.tex` (written arguments)

Evidence reviewed: commit `553e01c7` (HEAD of `/home/user/ml-cryptanalysis`). The probe evidence comes
from commit `7299db81`. The SHA-256 values in `experiments/ic-chain-depth10/verification.json` match the
committed probe script, the csd sources, `_elim.c` and `probe.json`. The working tree has an
uncommitted edit to `scripts/run_ic_chain_depth10_probe.py`. That edit is not durable evidence and is
not reviewed. No repository file was edited.

Independent small checks I ran are in the scratchpad (k <= 7, each under 1 minute):

- `check_ordinary.py`: my own ordinary-ring closure in the 2k chain variables, with the chain equations
  as inputs. It is compared with `chain_degree_closure` and with `shift_closure` (non-reaching,
  verified), at k = 2, 3, D = 2..6, p in {17, 97, 193, 257}, including order-two and (0,0) targets.
  Result: 0 mismatches in ordinary rank, quotient rank or "1 in V_D".
- `check_shift_deficient.py`: `shift_closure` in the deficient (non-reaching) regime at k = 6, 7,
  D = 4, with a small product cap (17-35 rounds) and a large one. Its rank equals the dense engine's
  every time, and verification passes. It reproduces dim V_4 = 603, 607, 605, 730, 857 and 829.
- `check_nocarry.py`: the same generator closed in the carry-free algebra K[u]/(u_v^2) (see Finding C).

Verdict scale: **sound**, **gap with fix**, **wrong**. Finding labels: wrong, unsupported, unpriced/unreplayed.

---

## Item 1. Lemma `lem:degree` and Corollary `cor:allm`

**The weight claim.** NF is linear. Each monomial x_i^e reduces to the squarefree monomial of
e mod 2^k: the chain ideal contains x_i^{2^k} - 1 and x_i^{2^j} - u_{i,j}, and squarefree
representatives are unique by Theorem `thm:reduced`. Reducing mod 2^k drops the binary digits at
positions >= k, so wt(e mod 2^k) <= min(k, wt(e)). For 0 <= e <= 2^{m-1}, wt(e) <= m-1 (the maximum,
m-1, is reached at 2^{m-1}-1; 2^{m-1} has weight 1). Collecting terms can cancel monomials but cannot
create new ones. So chain degree <= m*min(k, m-1). For cosets (x_i = c*u_{i,0}) only scalars change.
The fact that S_{m+1} has degree 2^{m-1} in each variable is Semaev's standard degree.

**Lean applicability.** `ChainQuotient.chain_lfd_bound` needs `Field K`, `(2:K) != 0`,
`chainInputs K m k ⊆ F` and a uniform input-degree bound delta. The direct presentation
F = G_{m,k} ∪ {NF(S_{m+1})} lives in `P (m*k) K` and contains `chainInputs` exactly. The sign
convention is the same, which `chainInputs_two_subset` checks for m = 2. Its degrees are <= max(2, mk).
So it applies, with LFD <= max(mk+1, delta) = mk+1. Neither 2^k | p-1 nor nonsingularity is needed.

**Remark (understatement, not an error).** `cor:allm` does not need `lem:degree`. Any NF modulo
G_{m,k} is squarefree in mk variables, so by Lemma `lem:leading` its degree is <= mk. The corollary is
therefore the kernel theorem `ChainQuotient.normalized_chain_lfd_bound` plus the definition of NF.
Calling it "written for m >= 3" understates it. The finer bound m*min(k, m-1) is true but plays no
role in the LFD bound.

**Verdict: sound.** Optional change: say that `cor:allm` follows from `lem:leading` and
`normalized_chain_lfd_bound` alone. Add the one-line uniqueness argument (via `thm:reduced`) to the
proof of `lem:degree`.

---

## Item 2. Soundness of `csd/ic_chain_shift_closure.py` as a computation of the ordinary V_D

### (a) Representation modulo the chains

Let G be the chain equations, J = <G>, N = NF and U_D the engine's space. U_D is the least subspace of
span(M_{<=D}) that contains N(mu*g) for wt(mu) + deg g <= D and is closed under w -> N(u_v w) for
wt(w) < D. Then V_{F,D} = (J ∩ R_{<=D}) ⊕ U_D, for D >= 2:

- G has pairwise coprime leading squares, so it is a Gröbner basis for a degree-compatible order in
  every characteristic. Every h in J of degree <= D therefore has a standard representation
  sum a_v c_v with deg(a_v) + 2 <= D. So J ∩ R_{<=D} ⊆ V_D, and h - N(h) ∈ V_{deg h} with
  deg N(h) <= deg h.
- U_D ⊆ V_D. A squarefree w ∈ U_D has ordinary degree wt(w), so u_v w is admissible when wt(w) < D,
  and N(u_v w) = u_v w - (element of J ∩ R_{<=D}).
- Conversely, W = (J ∩ R_{<=D}) + U_D contains the inputs and is closed under admissible products.
  For f = j + w with w = N(f), deg w <= deg f, so deg a + deg f <= D implies deg a + wt(w) <= D, and
  monomial multipliers iterate through intermediate weights < D. The engine's rule ("weight of the
  NF < D") is never more permissive than the ordinary rule, because N(h) is itself in V_D with
  ordinary degree wt(N(h)).

`generator_rows` (N(mu*g), wt(mu) <= D-2) lies in U_D, because a shift raises weight by at most 1
(carries only lower it). Carries are implemented as exponent addition mod 2^k (`shift_permutation`),
which is exactly u_{i,j} * x_i^a = x_i^{a+2^j} with x_i^{2^k} = 1. A row of pivot weight D is never
shifted, so a degree-D product is never "made affordable" by a carry. This matches the argument in
`docs/IC_FALL_COUNTEREXAMPLES.md` §3 that justifies `chain_degree_closure`. My independent
ordinary-ring implementation agrees at k = 2, 3.

**Verdict (a): sound.**

### (b) "Every low pivot processed once" plus the verification pass

- **Graded columns.** Columns are sorted by (weight, mask) descending, and `fp_echelon` returns row
  echelon form (not reduced). A combination of REF rows has leading column equal to the smallest
  pivot used. So any v in S ∩ M_{<D} uses only rows with pivot weight < D, and those rows lie entirely
  in M_{<D}. Rows with pivot weight < D therefore span exactly S ∩ M_{<D}. This is correct.
- **"Processed once" is not a closure proof.** Pivots are keyed by column and the basis is REF. The
  row at a processed pivot can be replaced by a different vector, and the docstring argument ("later
  tail reductions ... are processed in turn") does not cover a difference vector whose pivot column is
  already marked processed. The module docstring's "so when no low pivot is left the row space is
  closed" is therefore heuristic.
- **The verification pass is a valid closure proof.** It re-reduces the 2k shifts of every final low
  row against the final basis, and passing means rank is unchanged. The final S contains the seeds and
  every derived row is in U_D, so S ⊆ U_D. S is closed under the defining rule of U_D, so U_D ⊆ S.
  Hence S = U_D, and V_D follows by (a). This covers products of newly derived lower-degree rows
  (they are in S ∩ M_{<D} = span of low rows) and the budget rule (deg a + deg h <= D, by (a)).
- **The "reached" path needs no closure.** Every row lies in U_D, and reaching rank = dim(I ∩ M_{<=D})
  is a sound upper certificate on its own.

My non-reaching runs at k = 6, 7 (17-35 rounds with a small cap) matched the dense engine in every
case, including the degenerate instances with dim V_4 = 603, 730 and 857.

**Verdict (b): sound, with a wording gap.** The manuscript's depth-10 paragraph correctly rests
exactness on the verification. The module docstring and the sweep-section phrase "each low pivot
processed once" should not suggest that processing alone implies closure.

---

## Item 3. Prime local-to-global criterion (`prime_local_profile`, `shift_profile`)

**The instantiation of `LocalReductionBound.lfd_bound_of_local`.**

- A = K^Z, q = evaluation at Z (a K-algebra map).
- s = a linear section with values in span(M_{<D}). It exists because M_{<D} evaluates onto K^Z, and
  gives sigma = D-1 and room sigma+1 <= D.
- vanish: every input (chain equations and g) vanishes on Z.
- local_reductions: for deg p <= D, p - s(q p) ∈ I(Z)_{<=D}. Here N(V_D) ⊆ I(Z) ∩ span(M_{<=D}), which
  has dimension |M_{<=D}| - |Z| under onto. So the rank equality gives N(V_D) = I(Z) ∩ span(M_{<=D}).
  Then p - s(q p) = N(.) + (J ∩ R_{<=D}) ∈ V_D.

So HasLFDBound F D holds.

**Radicality and rationality.** These are correct but not needed for the upper bound. With
2^k | p-1 and p odd, K[u]/J ≅ K[(Z/2^k)^m] ≅ K^{2^{mk}} (Maschke, split). So every ideal containing
J is the vanishing ideal of its zero set, and <F> = I(Z). This is needed for the unit-ideal claim and
for converting "element of I(Z) of degree <= D not in V_D" into "element of <F>".

**The evaluation formula.** g0^{a1 E1 + a2 E2} is right. `subgroup_2power` returns an element of
order exactly 2^k, and the solution set is enumerated in chain coordinates.

**The exact-value rule is correct.** LFD <= D-1 iff I_{<=D-1} ⊆ V_{D-1} (given LFD <= D). With
V_D = I ∩ M_{<=D} and onto, dim(I ∩ M_{<D}) = |M_{<D}| - |Z|. So a strict fall at D is decided by
dim V_{D-1} against that number, and when there is no fall at D, descending gives LFD = last strict
fall. Quotient and ordinary dimension comparisons agree by item 2(a).

**Robustness.** An incomplete Z makes the target unreachable, because U_D ⊆ I(Z_true). The code
asserts that every supplied point vanishes on g, so false positives are excluded.

**Verdict: sound** as a written instantiation. The manuscript already says it is not a kernel
instance. Minor fix: state that the upper bound uses only "F vanishes on Z" plus rank/onto, and that
radicality is used for the unit-ideal and lower-bound statements.

---

## Item 4. Depth-10 deduction (`run_ic_chain_depth10_probe.py`, `probe.json`)

**The logic is sound.** A verified-closed space S equals U_5 = N(V_5) (item 2b). If 1 ∈ U_5, closing
1 under shifts (intermediate weights < 5) generates all of M_{<=5}, giving rank 21,700. So rank
18,569 means 1 ∉ U_5. The probe also checks the constant pivot directly; in REF with the constant
column last, "1 ∈ span" is equivalent to "last column is a pivot". 1 ∉ V_5 while 1 ∈ <F> gives
LFD >= 6. The arithmetic checks: binom(20, <=5) = 21,700 and binom(20, 6) partial sum = 60,460.
uint16 checkpoints are lossless for p = 12289. `rank_history` has 17 strictly increasing entries, so
there was no closure-phase resume (a resume would duplicate an entry). The 11 chunks match
ceil(3163/300).

**The conclusion does not depend on the unit-ideal claim (a strengthening the paper should state).**
For each x1 in H, S_3(x1, ., t) has degree <= 2 in x2. Its x2^2 coefficient is (x1 - t)^2. At x1 = t it
is linear with x2 coefficient -4(t^3 + a t + b) = -4 y_t^2 != 0 (y_t = 1126). So |Z| <= 2*1024 = 2048,
and dim(I ∩ M_{<=5}) >= 21,700 - 2048 = 19,652 > 18,569. LFD >= 6 holds even if the solution
enumeration were wrong, provided <F> = I(Z) (item 3).

**Evidence gap (unreplayed).** Exactness of V_5 rests on the verification pass. Its only durable
record is the status string in `probe.json`, a producer-written `verification.json` ("passed"), and
the hash of an uncommitted 770 MB basis. The per-chunk log is not committed. An independent check
requires the ~2 h / ~8 GB rerun; no cheaper certificate exists. Until a second party replays it and
the basis hash matches `d4fe1334...`, the result is a producer observation.

**Verdict: sound (logic), unreplayed (evidence).**

---

## Item 5. Claims beyond the evidence (findings ordered by severity)

### A. WRONG: "genericity" of dim V_3 and dim V_4, and the formal items built on it

The manuscript says (Section "Why the prime value steps at depth six", and Section `sec:missing`,
"Generic closure dimensions" and item 1):

- "At every depth dim V_3 = 2k+1".
- "For k >= 6 ... dim V_4 and dim(V_4 ∩ M_{<=3}) = 4k+87 no longer depend on the prime, curve, target
  or solution count".
- "dim V_4 = 4k^2+172k-571 on every curve and target tested".

The paper's own sweep (`experiments/ic-curve-sweep/results.jsonl`, field `quotient_ranks`)
contradicts this:

- **86 of 333 distinct prime systems at k = 6, 7 have dim V_4 != 4k^2+172k-571.**
  - All 81 two-torsion targets (y_t = 0) have 607 (k=6) or 835 (k=7).
  - The j=1728, target (0,0) systems have 730, 794, 1343 and 1471.
  - p=193, (a,b) = (23,169), target (182,69), subgroup: dim V_3 = 34 and dim V_4 = 603. Its target
    abscissa lies in H.
  - p=257, (34,176), target (17,28), coset c=206: dim V_4 = 857.
  - I reproduced 603, 607, 730 and 857 with both engines (`check_shift_deficient.py`).
- dim V_3 != 2k+1 for 3 systems at k=2 and 1 at k=6.
- **A kernel theorem in the paper's own Lean development refutes the proposed formal items.**
  `PrimeM2Sharp.one_at_four` (via `PrimeFieldConclusions.sharp_lfd`) says: for y^2 = x^3 + a x,
  target (0,0), a != 0 and a^(2^k) != 1, 1 ∈ V_4 at every depth k >= 2. So dim V_4 = binom(2k, <=4),
  which is quartic. Consequences:
  - "A proof of a quadratic upper bound on dim V_4" (depth-step paragraph) cannot exist.
  - `sec:missing` item 3 cannot hold as stated: "an upper bound dim V_4 <= f(k) with f quadratic,
    proved for all instances, gives LFD >= 5 for every system with |M_{<=4}| - |Z| > f(k)". That family
    is a unit ideal with LFD exactly 4.
  - `sec:missing` item 1, "never pivots on a quantity that vanishes for an actual curve and target, so
    that dim V_D = d_D(k) for every instance", is false at every threshold.
- **The proposed route is based on a false premise.** "The pivots are monomials in the top
  coefficients of S_3, which are nonzero for every nonsingular curve" is wrong. The top coefficients
  of the normalized generator are 1, -2(t^2+a) and -2t (`PrimeM2General.normalized`,
  `top (...) (-2*(t^2+a)) (-2*t)`). They vanish for t = 0 or t^2 = -a, independently of whether the
  curve is nonsingular.
- **The account of the (0,0) family is wrong.** "The (0,0)-target family ..., whose degree-4 ideal
  part is smaller than generic". For its unit-ideal members (8 coset systems at k=6, 4 at k=7) the
  degree-4 ideal part is all of M_{<=4}, the largest possible. For its proper members (|Z| = 2^k) it
  is 730 > 605. What is non-generic is the closure, not a smaller ideal part.
- Cheapest settling test: one-line scan of `results.jsonl` (k >= 6, `quotient_ranks[4]` against the
  fit), and reading `PrimeM2Sharp.one_at_four`.
- Fix: restrict every genericity statement to the five depth-step instances and the 30 depth-8/9
  instances (all random-rational or planted targets). Name the exceptional loci that the data already
  show (y_t = 0; t = 0 on b = 0; t ∈ cH). Restate items 1-3 of `sec:missing` for instances outside an
  explicit exceptional set.

### B. WRONG: abstract and Observations misstate the sweep counts

- **Abstract:** "prime chains have LFD 4 for k <= 5 and 5 for k = 6,...,9 independently of the prime,
  curve, target and of whether the base is the subgroup or a random coset."
  - 3 systems at k=2 have LFD 3: p=97 cosets on (27,42) and (95,81), and p=193 coset on (151,81).
  - 21 systems at k=6, 7 have LFD 4: all j=1728 with target (0,0); 4+8 at k=6, 3+6 at k=7.
  - k = 8, 9 come from separate 20- and 8-system arms that sample no order-two targets, and k=9 uses a
    single prime (7681). "Independently of the prime" is supported at k=9 by one prime only.
  - The values for k = 8, 9 are also attributed to the 3005 sweep systems, which stop at k = 7.
- **Observations:** "step ... to 5 at k = 6 and 7, for every prime, curve, target class and both base
  classes". This is false even as a per-cell maximum: every system in the cell (p, j1728, k, order-two)
  has value 4. "The step is a property of the depth" should read "of the depth, for targets other than
  (0,0) on b = 0 curves".
- "Consistent with the three systems of value 4 at each of k=6,7": the table in the same PDF shows 12
  (k=6) and 9 (k=7).
- The kernel family `sharp_lfd` requires a^(2^k) != 1, so it covers only the unit-ideal members. The
  subgroup systems with a = 1 ∈ H (64 or 128 points, LFD 4) are measured, not kernel-covered. The
  table entry "=4 at every k >= 2 for y^2=x^3+ax, target (0,0) (Lean)" must carry that hypothesis.
- Cheapest test: the paper's own `tab:prime` (`curve-sweep-data.tex`).

### C. WRONG (minor): "shifts that carry ... are the only source of strict falls"

The same generator (p=257, (2,3), t=58, k=6) closed in the carry-free algebra K[u]/(u_v^2) has a
strict fall at 4 and contains 1 at degree 4 (`check_nocarry.py`). The falls come from cancelling the
weight-2 top form and leaving g's lower-weight tail, not from carries. Here g is a unit in that local
algebra because its constant term a^2 - 4bt != 0. In the chain algebra a fall element is a sum of tail
and carry contributions. Fix: drop "only source", or prove it.

### D. UNSUPPORTED: "d_5(10) <= 18,569", "predicted by a dimension count", "a step ... is then forced"

- **"d_5(10) <= 18,569".** Specialising a generic derivation gives dim V_D(c0) >= d_D only on the open
  set U where the generic denominators and pivots stay nonzero. Outside U the instance can fall below
  the generic value (603 < 605 shows that this happens). The single k=10 instance is on the j=1728
  curve, where non-generic behaviour has already been seen, and is not shown to lie in U.
  - Cheapest test: a second k=10 instance on a random curve with a random target with y != 0 and
    t ∉ H. The modified probe script supports --a/--b/--seed. Equal rank 18,569 supports genericity;
    a different rank refutes "d_5(10) <= 18,569" as a generic statement.
- **"The next step predicted by a dimension count"** (abstract, Conclusion: "as the dimension count
  predicts"). The pre-run text (commit `cced8188`) predicted no k_0. It said "a step to 6 at some depth
  k_0 is then forced" and called k_0 "the single most informative missing measurement". No
  falsification target was declared, so this is post hoc. Change "predicted" to "consistent with".
- **The forced-step argument.** "a universal bound LFD <= 5 holds exactly when d_5(k) keeps pace with
  binom(2k, <=5), which a polynomial of degree below 5 cannot do; a step to 6 ... is then forced" has
  three problems:
  - It presupposes that d_5 is a polynomial of degree < 5, for which there is no evidence. For k <= 9
    the data show d_5(k) >= binom(2k, <=5) - |Z|, which grows with degree 5.
  - "Exactly when" conflates the generic d_5 with the all-instance statement.
  - It omits the onto condition and |Z|.

  The depth-10 refutation of "LFD <= 5 for all k" stands on its own; this argument adds nothing to it.
- **"LFD(k) = min{D : d_D(k) >= |M_{<=D}| - |Z|}"** needs the onto condition at D, and should use
  "=" rather than ">=".

### E. UNSUPPORTED: "No residual effect of invariance is detectable" and the Conclusion's "do not change the exact value beyond what their solution structure accounts for"

- No test is declared or reported. The paper's own stratified summaries are:
  - binary ordinary: invariant lower in 13 strata, higher in 5, weighted mean -0.282 (two-sided sign
    test p ≈ 0.096);
  - Koblitz: lower in 12, higher in 6, mean -0.18 (p ≈ 0.24).

  Both point estimates are nonzero and point the same way. Low power is not absence of effect.
- The manuscript omits the strongest raw contrast: on Koblitz proper ideals, random versus invariant
  non-subfield is higher in 2 cells and lower in 20 of 90 (p ≈ 1e-4). This is explained, not refuted,
  by the solution count.
- "Accounts for" is a mediation claim. The base class changes the solution structure (70% unit ideals
  on subfields; twice the solutions on invariant Koblitz bases). So symmetric bases do change the
  value, through Z.
- The prime "symmetric vs random" arm cannot test symmetry at all:
  - S_3 is weighted-homogeneous (x weight 1, a weight 2, b weight 3, total 4), so
    S_3^{a,b}(c x1, c x2, c t') = c^4 S_3^{a/c^2, b/c^3}(x1, x2, t'). A random-coset system is exactly
    the subgroup system of a scaled or twisted curve with target t/c.
  - cH is also closed under negation (-1 ∈ H) and under multiplication by H, the two properties the
    design lists as "symmetric".
- Fix: report the test and its interval, write "not established" instead of "no residual effect",
  state the mediation explicitly, and describe the prime coset arm as a reparametrisation control.

### F. UNSUPPORTED (wording): "Beyond the crossing ... the missing elements require degree-5 derivations, and the exact value is 5"

The dimension count gives only LFD >= 5, and only for instances with the generic dim V_4. "Exactly 5"
needs the degree-5 certificate, which was measured for k <= 9 and fails at k = 10. Fix: say so.

### G. UNREPLAYED: depth-10 evidence and "verified-complete" in the abstract

See item 4. The abstract and Conclusion present a producer-side verification as settled. Before it
counts as validated evidence, a validator replay with a matching basis hash is required. The abstract
also calls the sweep values "certified" next to "kernel-checked" without saying that the prime
local-to-global step is a written instantiation.

### Not applicable

The paper makes no cost, yield or rho comparison and says so. The §5a pricing band, the
fastest-external-implementation rule and the m=83 / ICV1 rules are not triggered.

---

## Verdict summary

| Item | Verdict |
|:--|:--|
| 1. `lem:degree`, `cor:allm`, Lean applicability | sound. Optional: `cor:allm` follows from `lem:leading` and `normalized_chain_lfd_bound` without `lem:degree` |
| 2a. Squarefree representation, shift = exponent addition mod 2^k | sound |
| 2b. Processed-once loop plus verification implies exactly V_D | sound, because of the verification. Gap with fix: "processed once" alone is not a closure argument |
| 3. Local-to-global criterion and exact-value rule | sound (written instantiation). Radicality is needed only for the unit and lower-bound statements |
| 4. Depth-10: closed V_5 of rank 18,569 implies 1 ∉ V_5 and LFD >= 6 | sound and robust to the unit-ideal claim (|Z| <= 2048 < 3131). Unreplayed evidence |
| 5. Abstract, depth-step section, `sec:missing`, Conclusion | wrong in A-C, unsupported in D-F, unreplayed in G |

## Required manuscript changes

1. **Abstract and Observations.** Replace "LFD 4 for k <= 5 and 5 for k = 6,...,9 independently of the
   prime, curve, target ..." with the actual counts: three LFD-3 systems at k=2; 21 LFD-4 systems at
   k = 6, 7, all j=1728 with target (0,0); k = 8, 9 from separate arms without order-two targets, k=9
   from one prime. Fix "three systems of value 4 at each of k=6,7" (12 and 9). Fix "for every ...
   target class" and "the step is a property of the depth".
2. **Depth-step paragraph and `sec:missing`.** Withdraw "dim V_3 = 2k+1 at every depth" and the claim
   that dim V_4 is independent of curve, target and solution count "on every curve and target tested".
   Restrict both to the depth-step and depth-8/9 instances, and list the sweep exceptions (86 of 333 at
   k = 6, 7).
3. **(0,0) family.** Replace "whose degree-4 ideal part is smaller than generic" with "whose degree-4
   closure is non-generic (1 ∈ V_4 when a^(2^k) != 1, kernel-checked)". Add the hypothesis
   a^(2^k) != 1 to the Lean entry in `tab:twosided`.
4. **Open formal items.** Withdraw the open formal item "a quadratic upper bound on dim V_4" and
   `sec:missing` items 1 and 3 as stated, since they are refuted by `PrimeM2Sharp.one_at_four`.
   Restate them for instances outside an explicit exceptional locus. Correct "top coefficients ...
   nonzero for every nonsingular curve" (they are 1, -2(t^2+a), -2t).
5. **Mechanism paragraph.** Delete or prove "shifts that carry ... are the only source of strict falls".
6. **Depth-ten paragraph and conclusion.** Replace "d_5(10) <= 18,569" with "dim V_5 = 18,569 for this
   instance; equal to d_5(10) if the instance is generic". Replace "predicted by a dimension count" and
   "as the dimension count predicts" with "consistent with". Remove or mark as conjecture the "forced
   step" argument that presupposes deg d_5 < 5. Add the onto condition and "=" to the conjectural
   formula.
7. **Depth-ten robustness.** State that LFD >= 6 does not depend on the unit-ideal claim
   (|Z| <= 2048 < 21,700 - 18,569). Label the verification as producer-side until a validator replay
   matches the basis hash.
8. **Base classes.** Replace "No residual effect of invariance is detectable" and the Conclusion's
   "do not change the exact value beyond what their solution structure accounts for" with the
   stratified estimates, a declared test, and the explicit mediation statement. Report the Koblitz
   proper-only 2:20 contrast. Describe the prime coset arm as a reparametrisation of the subgroup
   system (S_3 weighted homogeneity; cH is closed under negation and H-multiplication).
9. **Wording.** "Beyond the crossing ... the exact value is 5" should say that the upper bound 5 is a
   separate measured certificate (k <= 9) and fails at k=10. The sweep-section phrase "each low pivot
   processed once" should say that exactness, when claimed, rests on the verification pass.
10. **Optional (item 1).** Note that `cor:allm` is `lem:leading` plus `normalized_chain_lfd_bound`
    (only the definitional "NF is squarefree" step is unformalized). Add the uniqueness argument for
    NF of a monomial to `lem:degree`.
