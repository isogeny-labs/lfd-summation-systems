# Red-team audit of `paper/curve-lfd-bounds.tex` (2026-11-12)

Frozen evidence reviewed: commit `1c5120db796f3b162dde3bb0c0353ce53b0d6460` (working tree clean for
`paper/` and `scripts/make_paper_figures.py`), `paper/curve-lfd-bounds.tex` (2143 lines),
`paper/curve-sweep-data.tex`, `paper/figures/*.pdf`, and the records named in Appendix C. No file
was edited and no experiment was rerun. Line numbers refer to `paper/curve-lfd-bounds.tex`.

Classes used: **wrong** (contradicted by the record, the cited source or a computation),
**unsupported** (no record or argument backs it as stated), **inconsistency** (the paper disagrees
with itself), **reference**, **prose**. "Unpriced" does not arise: the paper makes no cost claim.

## What was checked

* Every number in the abstract and contribution list, Tables 1-8 and Figures 1-6 (numbering of the built PDF), and the
  Limitations, against `experiments/ic-bounds-k14/RESULTS.md` and `runs/*.json` (all closure and lift
  records), `ic-curve-sweep/{summary,cases,results,base-class-analysis}.json*`,
  `ic-chain-depth8/9/10`, `ic-bounds-d5/verify/*.json`, `ic-bounds-m3/README.md`,
  `ic-bounds-m3-depth/closures.jsonl`, `ic-bounds-msummands/{README.md,pod-scores.json,records/}`,
  `ic-bounds-binary{,-law,-deep}/*`, `ic-bounds-freeext/{predictions,outcomes}-2026-11-10.json`,
  `ic-bounds-locus/`, `ic-bounds-cert/`, `ic-bounds-exact-k6/`, `docs/proofs/N6-genericity.md`,
  `docs/reviews-U1-upper-bound-redteam-20261110.md`, `docs/IC_BOUNDS_PROGRAM.md`.
* All 75 `\code{}` names in Table 7 (Appendix A) against `formal/ICFallDegrees/*.lean`,
  `Audit.lean` (582 `#print axioms` lines) and `experiments/ic-lean/verification.json`
  (582 declarations, axiom sets only `propext`/`Classical.choice`/`Quot.sound`, all 124 source
  hashes equal to the checkout). Lean statements of the selector, three-summand, first-fall,
  sharp-family, bridge, restriction, depth-14 and growth theorems were read for hypotheses.
* Identities with sympy 1.14: the selector identity (line 686), symmetry of S_3, the unit-lemma
  identity (873-874) and its x1 = x2 case, the step to line 882; Λ(1..7), k_D, k_D <= 14D+10;
  all |M_{<=D}| counts, the recursion and its unrolled cubic, 829+230(k-7)+2C(2k-14,2), the
  lift arithmetic, every "ideal part" subtraction, sign tests (p-values), and the numeric facts of
  the 65537, 163841, 114689, 1153, 40961, 12289 instances.
* `scripts/make_paper_figures.py` data functions were executed in memory (no figure written): all
  asserts pass; the rendered PDFs were inspected.
* References: Petit-Kosters-Messeng (PKC 2016 PDF, Tables 4-5, Section 3.3, System (4)), Kosters-Yeo
  arXiv:1503.08001v3 (full text), HKYY arXiv:1505.02532 (Def. 2.1, 2.5, Prop. 2.6, 2.8),
  Caminata-Gorla arXiv:2112.05579 (Thm 2.9, Ex. 4.3), arXiv entries for Huang 2103.07282,
  Kousidis-Wiemers 1906.05594, Steiner 2310.03637, Guo et al. 2602.12772, ePrint 2015/573, and the
  JSC record of Gorla-Mueller-Petit.

## Verified findings, by severity

### F1 (wrong; must fix). The PKM comparison is not "the same presentation"
Lines 297-300, 330-333, 1102-1104, 1475-1477, 1768-1772.
The paper's system is `{NF(S_3)} ∪ G_{2,k}` with a chain-reduced generator of degree 2 (lines
376-382). PKM's System (4) (PKM Section 3.1, Algorithm 1; Section 3.3 sets L_j(x) = x^{p_j}, L_{n'}(x) = 1 - x^{p_{n'}})
uses the unreduced `S_{m+1}(x_{11},...,x_{m1},X)`, which has degree 4 for m = 2 (PKM Section 4.2, after Table 5, compares
with "n equations of degree 2 and 1 equation of degree 4"). The paper itself stresses that the LFD
depends on the presentation (line 258) and the Lean development proves exactly this for raw versus
normalized `S_3` (`PrimeCounterexample.normalized_lfd_three`, `raw_lfd_four`,
`lfd_is_not_ideal_invariant`). The matching is still true for every value >= 4: S_3 - NF(S_3) is
derived from the chain equations within degree 4 (x1^2 x2^2 - u11 u21 = x1^2(x2^2-u21) + u21(x1^2-u11),
etc.), so V_D(raw) = V_D(normalized) for all D >= 4 and V_3(raw) ⊆ V_3(normalized); hence the two LFDs
agree whenever the normalized LFD is >= 4. That argument is not in the paper.
Minimal fix: replace "the same presentation" by "PKM's presentation (unreduced S_3)" and add the
two-line transport argument above, with the exception at value 3.
Cheapest test: none needed beyond the argument; optionally rerun one k = 6 PKM-form system.

### F2 (inconsistency / overclaim; must fix). Uncertified value 7 described as exact or certified
* Line 903-906: "All values in this section are exact on the instances named, certified as in
  Section 3.3" — the section contains the lift value 7 (no certificate, line 1196-1197).
* Lines 655-657 (Table 2 caption): "measured values are exact and certified on the seeded grids";
  the same table lists "7 (k = 14, three instances, by degree lift)", which is neither certified nor on
  a seeded grid.
* Line 1165-1166: "The exact value was obtained by a degree lift"; Table 6 column "exact LFD" row
  14 shows 7. By the paper's own label definition (lines 93-94: exact = "ends with a certificate"),
  the lift result is not exact.
* Line 232-233: "the LFD reaches 6 at k = 10 and 7 at k = 14" with no qualifier; line 1594-1595
  "the measured values (7 at k = 14)".
Minimal fix: restrict line 903 and the Table 2 caption to "except the upper bound at k = 14";
introduce one label for the lift value (e.g. "computed by the degree lift, written proof, no
certificate") and use it at 232, 1165, 1461 and 1595.

### F3 (overclaim; must fix). "LFD = deg g + 1 for four or five summands" is mostly forced by proved bounds
Abstract line 62-63; lines 1355-1363; Figure 4 caption.
In 4 of the 5 four/five-summand cells (m = 4, k = 2, 3; m = 5, k = 2, 3) deg g = mk (Lemma 8 gives
m·min(k, m-1) = mk; records: deg 8, 12, 10, 15), so LFD <= mk+1 (Corollary 9, kernel) and
LFD >= deg g + 1 (the carry criterion of Theorem 13, hypotheses checked per instance) force
LFD = mk+1 = deg g + 1. `experiments/ic-bounds-msummands/README.md` says so: "the runs confirm the
implementation and the hypotheses, they are not an independent measurement of an unknown number".
Only m = 4, k = 4 (deg g = 12 < 16) carries information on the upper side. The README also discloses
two smoke-test closures (m = 4, 5 at k = 2) run before the prediction freeze; the paper says only
"the predictions were recorded before the runs" (1361-1362).
Minimal fix: state that at deg g = mk the value is forced, that m = 4, k = 4 is the one informative
cell, and carry the smoke-test disclosure; drop or qualify the abstract sentence.

### F4 (wrong). Two false generalizations about three summands and the m-summand figure
* Line 1330-1332: "So for three summands, at k <= 5, the only fall is at deg g + 1". On the seeded
  m = 3 grid, the locus 5t^3+3at+6b = 0 has deg NF(S_4) = 5 and LFD = 7 = deg g + 2 at k = 3 and 4
  (10 and 11 cases; `experiments/ic-bounds-m3/README.md`, "Reading", and its k-table: "direct,
  deg NF(S_4) = 6: 7: 107" of 117 at k = 3, 106 of 117 at k = 4, all 117 with value 7).
  Fix: "for the twelve generic instances (5t^3+3at+6b != 0)".
* Lines 1377-1378 (Figure 4 caption): "m = 4, k = 4 is the cell where deg g = 12 < mk separates the
  two". m = 3 at k = 3, 4 also has deg g = 6 < mk (bounds 10 and 13 against value 7); the figure
  itself draws open markers there. Fix: "the cells with deg g < mk (m = 3, k = 3, 4; m = 4, k = 4)".
* Line 1319-1320: "Theorem 13 proves the lower bound 7 at every depth outside this locus" —
  Theorem 13 is stated for k >= 3; k = 2 rests on the per-instance check. Fix: "at every depth k >= 3".

### F5 (overclaim). The dimension identity at depth >= 7 is attributed to the kernel theorem
Lines 1025-1031. "V_4(k) is the free high-level extension of V_4(7), of dimension
4k^2+172k-571, by the depth reduction (Theorem 15, kernel-checked)". Theorem 15 gives that
W_k(S) is closed and contains g, hence only V_4(k) ⊆ W_k(S), i.e. dim V_4(k) <= 4k^2+172k-571.
Equality (W_k(S) ⊆ V_4(k)) is not in the theorem; it is an exact computation at k = 8..12 on three
instances (`experiments/ic-bounds-freeext/README.md`: canonical RREF equality) and dimension
agreement to k = 13. The base facts are kernel-checked only for p = 40961 and p = 11, not for the
"depth-step curves" this paragraph discusses. Minimal fix: "the kernel theorem gives the upper bound
dim V_4(k) <= 4k^2+172k-571 per instance whose base facts hold; equality is an exact computation at
k <= 12 (dimensions to 13)".

### F6 (mathematical; definition does not match the theorem). Depth reduction
Lines 780-792. W_k(S) is defined only for a subspace of "the depth-7 chain algebra (degree at most
4)" with high patterns "on levels 7, ..., k-1" and S ∩ M_{<=4-wt q}; Theorem 15 then quantifies
"for every subspace S and all k_0, D, m". As written the theorem uses an undefined object. The Lean
definition is general (`DepthReduction.W k0 D S k`, weight budget `lowWt l + lowWt q <= D`).
Minimal fix: define W^{k_0,D}_k(S) with levels k_0..k-1 and S ∩ M_{<=D-wt q} for m blocks, then
specialise to k_0 = 7, D = 4.

### F7 (wrong). Appendix C promises a README in every data directory
Lines 2033-2035: "Each directory contains a README.md". Seven listed directories have none:
`experiments/ic-curve-sweep/`, `ic-chain-depth8/`, `ic-chain-depth9/`, `ic-depth-step/`,
`ic-chain-depth10/`, `ic-bounds-binary-deep/` (has RESULTS.md), `ic-curve-validation/` (has
REPORT.md). Fix: add the READMEs or reword ("each directory contains a README.md, RESULTS.md or
REPORT.md, or a summary.json ...").

### F8 (inconsistency). "The bound exceeds LFD >= 4 only from k = 52"
Line 1594. For the displayed bound floor((k-10)/14)+1 (line 1539) and Figure 6(b), the value exceeds
4 only from k = 66 (it equals 4 on k = 52..65). k = 52 is right only for the k_D form
(k_4 = 4 + Λ(3) = 52 gives LFD >= 5). Fix: "the k_D form of the theorem exceeds ... from k = 52
(the linear form from k = 66)".

### F9 (inconsistency). "No prediction recorded in advance has failed"
Line 1774-1775 (Conclusion). The paper reports preregistered binary laws refuted (lines 1390-1414,
abstract line 65) and the m-summand point rule had 3 misses (`ic-bounds-msummands/README.md`, G5).
Fix: "none of the six predictions recorded for the two-summand prime law has failed".

### F10 (reference, wrong). PKM Table 5 misreported; "Magma's corrected degree"
* Lines 304-305: "Table 5 ... shows 5 up to kernels of order 2^9". PKM Table 5 has D^corr = 4.00 for
  N = 2^1·4 (p = 11), 5.00 for 2^2 ... 2^9, 6.00 for 2^10, 2^11. Fix: "4 at 2^1, 5 from 2^2 to 2^9".
* Lines 202 and 298-299: "Magma's corrected degree of regularity". The correction (removing steps
  with no pair reduced) is PKM's own post-processing of Magma's verbose output (PKM Section 4.2); the
  uncorrected D_av does not show the clean breakpoints (4.48 at n1 = 5, 5.36 at 7 and 8, 5.78 at 9,
  6.36/6.29 at 10/11). The interpretive sentence at 331-333 ("reflects the LFD ... and not a property
  of one solver") should mention that the agreement holds for the corrected measure only.
* Verified correct: PKM Table 4 values (4.00 for n1 <= 5, 5.00 for 6..9, 6.00 for 10, 11), 100 runs,
  smallest prime with 2^{n1} | p-1 and p >= 2^{2 n1}, nbsols 0.59-1.34, Section 3.3.

### F11 (reference / unsupported). The binary written bounds rest on an uncited internal draft
Lines 648, 654-655, 662-665, 1743-1745. "The companion manuscript" is cited three times and supports
the written bounds min(k, l)+2 and floor(k/2)+1. It is not in the bibliography or in Appendix C;
the repository copy `paper/binary-lfd.tex` has no author and is dated "Internal review draft, 30
October 2026". Fix: cite it (with a stable location) or label these entries "written, unpublished".

### F12 (reference). Kosters-Yeo claims overstated
Lines 227-228 and 315-318. KY (arXiv:1503.08001v3) prove NP-completeness "under an assumption"
(their Assumption 3.3; Theorem 3.5(i)); the unconditional part, Theorem 3.5(ii), is for the singular curve y^2 = x^3. The first fall
degree 2 is shown "in general" ("in practice, this almost always happens", Section 4). The growth of
the degree of regularity with n is in Section 4 (Remark 4.12 table), the dependence on the subspace
in Section 5. Fix: "under an assumption"; "is 2 in general".

### F13 (overclaim). The k = 14 control on p = 114689
Lines 1239-1242: "a lift that overshoots would be caught". The bound dim V_7 <= |M_{<=7}| - |Z| detects
only lifted elements outside the ideal; an unsound lift that adds ideal elements not derivable in
degree 7 lands below or at the bound and is not caught. `experiments/ic-bounds-k14/RESULTS.md`
words it correctly ("A lift that produced anything outside the ideal would have exceeded it"). Since
the unit-instance value 7 is exactly "V_7 = M_{<=7}", the control rules out the "always full" failure
mode but not every unsoundness. Fix: use the RESULTS.md wording.

### F14 (scope). The p = 11 "instance" is never a subgroup system
Lines 110-111, 799-805, 815-818. Over F_11 there is no subgroup of order 2^k for k >= 2 (4 ∤ 10), so
"LFD >= 5 at every depth k >= 8 for p = 11" concerns a polynomial system that never encodes a factor
base; likewise p = 40961 for k >= 14 (2^14 ∤ 40960). The paper makes this caveat for 12289 (lines
1107-1110) but not here, and the contribution list calls both "concrete instances" of the curve
systems. Fix: one sentence stating where each is a subgroup system.

### F15 (inconsistency). Evidence labels are not on every statement
Line 90: "Every statement carries one of five labels". Lemma 4, Theorems 5 and 6, Corollaries 7
and 9, and Theorem 12 carry no label in their headers (their status is only in Appendix A); the
claims "the same generator closed in the carry-free algebra K[u]/(u^2) still has a fall at degree 4"
(405-407) and "the method cannot give a slope below 1" (719-720) carry no label and no record. Fix:
label each theorem header (kernel / written) and label or remove the two remarks.

### F16 (inconsistency). Appendix A names
Line 1808: `QuotientLowerBound.closure_le_quotient_span` does not exist; the declaration is
`FiniteQuotient.closure_le_quotient_span` (namespace `ICFallDegrees.FiniteQuotient`, file
`QuotientLowerBound.lean`, `Audit.lean` line 70). Line 1807 lists `HomogeneousFirstFall` (a module)
as a declaration, followed by `degree_three_syzygy_zero`, which is
`GeneralQuadraticClassification.degree_three_syzygy_zero`, not in `HomogeneousFirstFall`. Several
other short names (`standard_le_one`, `chain_lfd_bound`, `prime_m2_lfd_bound`, `lfd_ge_four_of_unit`)
are ambiguous under the stated convention "names relative to the project namespace". All other names
exist and are audited. Fix: fully qualified names.

### F17 (figure). Figure 2 (sweep)
`paper/figures/curve-sweep.pdf`: the legend says "bound 2k+1 (right axis)" but no right axis is drawn;
the dashed line is plotted on the count axis (5 to 15 systems), where it has no meaning. The caption
(981) says "Capped systems are shown separately"; there are 0 capped systems and none are shown.
Fix: draw the bound on a labelled right axis (or as text) and drop the capped sentence.

### F18 (title framing). "Random and symmetric factor bases" for the prime family
Title (32), abstract (49), contribution 3 (116), auto-table header "subgroup (symmetric)"
(`curve-sweep-data.tex` line 10). Lines 951-955 state that the prime arm "does not test symmetry"
because a coset cH is also closed under negation and multiplication by H. Fix: describe the prime
arm as "subgroup and coset bases"; keep "random versus Frobenius-invariant" for binary/Koblitz.

### F19 (minor errors of fact and wording)
* Line 1391-1392: "80 new systems" — the record has 80 cases, 60 distinct systems
  (`ic-bounds-binary-law/README.md` line 7). Line 1396: "unit ideals reach 4 at k = 7" comes from the
  earlier grid ((n,k) = (12,7)), not from the preregistered test (P2, P3 have units with 4 at k = 9, 10).
* Lines 648 and 1384-1385: "<= 5 on the grids through k = 10", while the scaled subfield has 6 at
  k = 10 (Section 5.8). Qualify as "the earlier grids".
* Line 1057-1058: invariance "fails for μ ∉ H"; it fails only for μ with μ^2 ∉ H
  (`docs/proofs/N6-genericity.md` line 280); for μ ∉ H with μ^2 ∈ H the coefficients a stay in the
  same class.
* Lines 1182, 1221-1222, 1242 ("kept" basis) against 1665-1666 ("no k = 14 degree-6 basis ... was
  kept"): say "kept for the checker and the lift, not archived".
* Line 1675-1676: "the free-extension structure of V_6 starts at k = 14" — by the paper's own
  convention for V_5 (1286-1288: structure from base 11, recursion from 10) and prediction
  `d6-base-15`, the base would be 15 and the recursion from 14.
* Line 312: Caminata-Gorla Example 4.3 is the chain subsystem {x1^2-1, x2^2-1} alone, not "the
  case k = 1 of the chain system" (with the generator, k = 1 has LFD up to 3; U1 review table).
* Lines 184-185: Frobenius is a purely inseparable endomorphism of degree 2, not "a curve
  automorphism"; it acts as an automorphism of the group E(F_{2^n}).
* Line 1087: "Two further instances give the same degree-5 closure" — same dimension and profile;
  the spaces differ (different basis hashes in `ic-chain-depth10/*.json`).
* Line 838: "closure dimensions (0,0,2,111,605)" is the profile dim(V_4 ∩ M_{<=j}) (111 = 4k+87 is the weight-<=3 part, not dim V_3, which is generically 2k+1 = 13). The
  definition of "low profile" (349-350, j <= D) does not match its use (1152, 1342: j < D).
* Line 1720: "a fall of dimension 2k+86" for the fixed space V_4(7): give 100, or state it for V_4(k).
* Line 1768-1772: PKM data stop at n1 = 11; values at k = 12-14 extend beyond their table and do not
  "reproduce" it.
* Figure 3(b) (closure-dims): the D = 3 codimension 10 is not read from a record; the script takes it
  from the text of RESULTS.md and checks it against 2k+1 (`scripts/make_paper_figures.py`,
  `closure_dims`). The caption of 3(a) does not say that V_4 and V_5 at k <= 9 come from p = 12289 and
  the rest from p = 65537.
* Growth window (1580-1583): correct, but for 24 <= k <= 37 the Lean conclusion is d >= 2 (vacuous);
  one clause would explain why Lean says k >= 24 and the theorem k >= 38.

### F20 (reference, bibliographic)
* `silv` (Silverman) is in the bibliography and never cited.
* `gmp21`: add J. Symbolic Comput. 109 (2022) 386-398.
* `kw19`: published in J. Math. Cryptol. 13(3-4) (2019) 229-237 (doi 10.1515/jmc-2017-0022).
* `hky15`: pages 581-600 missing.
* `guo26`: authors are Junyu Guo, Hao Shen, Junqi Liu and Lihong Zhi (initials missing; format
  differs from all other entries).
* `ky` has a v3 (June 2015); cite the version used.

## Suspicions (not fully verified)

* **S1. Semaev attribution (lines 182-183, 282-284).** Semaev's ePrint 2004/031 introduces summation
  polynomials; the descent to F_2-linear coordinate spaces is usually credited to Gaudry (2009),
  Diem (2011), Faugere-Perret-Petit-Renault (2012) and Petit-Quisquater (2012). I did not read
  Semaev 2004 in full. Cheapest test: read its Section on factor bases.
* **S2. Brownawell-Masser and Voloch (1555-1562, 1589-1590, 1683-1684).** The constant
  6(m-1)(m-2) is consistent with H <= (n-1)(n-2)(|S|+2g-2)/2 for g = 1, |S| = 6 and the pole orders at
  O and R; I did not verify the cited statements, and the paper says neither review could. The
  analytic input of the growth theorem is therefore attributed to sources no one in the chain has
  checked. Cheapest test: obtain the two papers and compare the exact hypotheses (char p: p-th power
  coefficients).
* **S3. Steiner (313-314).** The arXiv abstract mentions iterated systems but not coprime leading
  terms; the attribution may still be right in the body. Not checked.
* **S4. CLO citation (549).** "[Ch. 2, §9]" for coprime leading terms matches the 4th edition's
  "Refinements of the Buchberger criterion" as far as I recall; not checked against the book.
* **S5. HKYY numbering.** Prop. 2.8 and Def. 2.1/2.5 match the arXiv version; the JSC version cited
  may renumber. Not checked.
* **S6. First binary grid "28 systems with value 5" (1385-1387).** `ic-bounds-binary/README.md`
  says 30 systems; `ic-bounds-binary/VERIFY.md` line 152 says "28 systems / 30 cases". The paper
  follows VERIFY; the README table (184 cases = 184 systems at k = 10) disagrees with it. One of the
  two records is wrong.

## Prose: the ten worst places

1. 794-796: "not up to 2 as in an earlier written statement; the two versions agree on the computed
   instances" — history of the project, not mathematics.
2. 1993-2010 ("Proof reviews"): "found the overclaims corrected in an earlier version", "pinned to
   its commit" — process talk about reviews and repository states.
3. 2012-2014 ("Prior art. A literature check located ...") — says nothing a reader can use.
4. 2016, 2024-2025: dated replay log ("passed on 4 November 2026", "replayed in full on 8 November
   2026").
5. 1880-1885, 1899-1900, 1928-1929, 1955-1959: run times, gigabytes, "checkpointed so that it
   survives a restart", "spills tiles to host memory, on a GPU with 141 GB" — operations, not methods.
6. 1943-1944: "its record also states which parts of the engine its author had read" — agent
   process; state the independence boundary instead.
7. 1164-1165 and 1753: "on the hardware we used", "out of reach on our hardware".
8. 1066-1067: "about 52 CPU-hours" inside the genericity discussion.
9. 1979-1980: "An independent validator wrote its own exact closure" — an actor, not a method
   ("an independent implementation").
10. 1699-1702: "the upper bounds of Table 2 hold for every instance without genericity assumptions,
    the lower-bound families are exact" — unclear ("families are exact"?) and wrong for the written
    entries of that table.

Also: notation clashes — A is the polynomial ring (190) and the group algebra (396); G is a subset
of F (354), the chain set G_{m,k} and the group (Z/2^k)^2 (396); K is the field and the lift subspace
(1173-1178); g is the generator and the subgroup generator (1858); a is the curve coefficient and a
multiplier in the same theorem (747-749); N is 2^k (367) and the projector (690); p is the prime and
a polynomial (772-773). Engine names drift: "independent engine" (1118), "second engine" (1100),
"reference engine" (1664, 1931), "independent reference engine" (1902), "third engine" / "GPU
engine". Informal phrases: "with no hypothesis left" (110, 813-814), "The value does not move with
the depth in the range we can reach" (1324).

## Checked and consistent (coverage only, not endorsement)

Abstract and Table 6 numbers (18,569 ... 66,073; 60,460 ... 313,912; 448,922; 499,178; 1,683,218;
codimensions 10/189/3131/50,256; part ranks; lift counts 586,501/5519/30,655/30,639 and falls per
round; 1,680,848; 46,188/47,886/43,431); sweep counts (3163/3005, histogram, the 3 systems at k = 2
and 21 at k = 6, 7, the 86 of 333 deviations, the four dim V_3 != 2k+1 systems, all with 3t^2+a = 0);
sign-test p-values; depth-8/9 arm counts; m = 3 counts (234, 105/116, 1092, k = 6 bound with
|Z| = 48, 19,958, 62,956); binary-deep counts (57, 24, 45/52, trace-zero 5/6/6, subfield 6/7);
selector test counts (149,122, 132,303, 2570); prediction file timing (committed 03:06 UTC, before
the k = 14 runs) and outcomes; Λ(1..5) = 25, 36, 48, 59, 71 and k_D = 27, 39, 52, 64, 77; the worked
example (LFD 3, FFD 3, LFD 2 for F ∪ {x}); the LFD definition matches HKYY Def. 2.5; unit lemma
and existence argument; Lean hypotheses of Corollary 7, Theorems 11, 12, 13, 17, Lemma 18, the depth-14 summary
(four recorded closures as hypotheses) and the growth theorem (CurveModel hypothesis).

## Verdict

Must fix before submission: F1 (PKM presentation), F2 (the uncertified 7 called exact/certified in
four places), F3 (the four/five-summand "finding" is mostly forced by proved bounds; disclose), F4
(two false generalizations and the Figure 4 caption), F5 (kernel theorem credited with an equality
it does not give), F6 (depth-reduction definition), F7 (missing READMEs promised in Appendix C), F9,
F10, F11 (uncited companion manuscript). Should fix: F8, F12-F18 and the bibliography (F20). The
Lean side is clean: every audited declaration exists, the audit record matches the sources, and the
kernel-checked statements carry the hypotheses the paper states, apart from the naming slips in F16.

## Resolution (applied to the working tree after `39d6a8bf`)

Each finding was re-checked against the record, the Lean source or the cited paper before it was
applied. Theorem numbers below are those of the audit; the rebuilt PDF renumbers from Remark 11 on.

* **F1 applied.** PKM System (4) (PKC 2016 PDF, Algorithm 1, Section 3.3: `L_j(x) = x^2`,
  `L_n'(x) = 1 - x^2`) uses the unreduced `S_{m+1}`. Section 2.5 now says so and adds a written
  transport argument, checked by hand: `S_3 - NF(S_3)` lies in the degree-4 closure of the chain
  equations, so `V_D(raw) = V_D(normalized)` for `D >= 4`, while for `D <= 3` the raw closure is the
  chain closure and misses `NF(S_3) != 0` (Theorem 5); hence `LFD(raw) = max(4, LFD(normalized))`,
  equal at every value 4-7 of Table 6, different only for the three value-3 systems at `k = 2`. The
  Lean counterexample `PrimeCounterexample.lfd_is_not_ideal_invariant` (F_5, depth 1, target (0,0) on
  `y^2 = x^3 + 2x`: normalized 3, raw 4 = max(4, 3)) is cited (name in Appendix A). "Same
  presentation" is replaced in the abstract, contributions, Section 5.4, Table 6 caption and the
  Conclusion.
* **F2 applied.** New label *lifted* for the `k = 14` upper bound (Evidence labels); qualified in the
  abstract, Section 2.3, the Section 5 opening sentence, Table 2 (cell and caption), Table 6 (header
  "LFD (exact unless marked)", row 14, caption), the three "exactly 7" sentences of Section 5.5,
  the conjecture paragraph, Section 6 Scope and the Conclusion.
* **F3 applied.** Abstract, Section 5.7, Figure 4 caption and Limitations now say that in four of
  the five m = 4, 5 cells `deg g = mk` and the value is forced (kernel bound plus the carry
  criterion with per-instance hypotheses), name the informative cells (m = 4, k = 4; and m = 3,
  k = 3, 4, where `deg g = 6 < mk` per `ic-bounds-m3/README.md` and `ic-bounds-msummands` records),
  and disclose the two smoke-test closures exactly as `predictions-2026-11-12.json` does.
* **F4 applied** (all three corrections; the locus values 7 at k = 3, 4 on 10 and 11 cases are stated).
* **F5 applied, with one correction to the audit.** Theorem 15 now gives only
  `dim V_4(k) <= 4k^2+172k-571`. The space equality is recorded in `ic-bounds-freeext/summary.json`
  only at k = 8 (three generic instances) and k = 9 (two), not k = 8..12 as the audit states; the
  paper says k = 8, 9, and dimensions at k = 6..13 (reference engine, `ic-bounds-d5/VERIFY.md`).
* **F6 applied.** `W^{k_0,D}_k(S)` defined for m blocks, levels k_0..k-1, matching `DepthReduction.W`.
* **F7 applied** by rewording Appendix C (no README files were added).
* **F8 applied** (k_D form from k = 52, linear form from k = 66).
* **F9 applied** in the Conclusion; the abstract sentence already concerns the law only.
* **F10 applied.** Verified against the PKM PDF: Table 4 `D_av` 4.00/4.00/4.00/4.01/4.48/5.00/5.36/
  5.36/5.78/6.36/6.29 and `D^corr` 4.00 (n1 <= 5), 5.00 (6-9), 6.00 (10, 11); Table 5 `D^corr` 4.00 at
  N = 2·4 (p = 11), 5.00 at 2^2..2^9 (times a cofactor, not always odd: 2^4·22), 6.00 at 2^10, 2^11.
  The correction ("removing any step in which no pair was reduced") is PKM's own post-processing
  of Magma's verbose output (Section 4.2). Sections 2.2 and 2.5 report both measures.
* **F11 applied, option "move out of the bounds table".** `paper/binary-lfd.tex` has an empty
  author field and is dated "Internal review draft, 30 October 2026". Labelling its statements as
  written bounds of this paper would overstate them, so the binary `min(k, l)+2` and
  `floor(k/2)+1` entries were removed from Table 2 and stated in Remark 11 as unpublished, unrefereed
  statements that this paper neither claims nor verifies; Appendix C gives the file location.
* **F12 applied** (verified against arXiv:1503.08001v3: Theorem 3.5(i) under Assumption 3.3,
  3.5(ii) unconditional for the singular curve; first fall degree 2 "in general"; growth of the
  degree of regularity in Remark 4.12; coordinate-space dependence in Section 5).
* **F13 applied** (RESULTS.md wording; the undetected failure mode is stated).
* **F14 applied** (Section 4.5 and the contribution list).
* **F15 applied.** Kernel-checked labels on Lemma 4, Theorems 5, 6, 12, Corollaries 7, 9 ("written
  proof" on Lemma 8); the carry-free remark is labelled empirical with its instance (p = 257,
  (a,b) = (2,3), t = 58, k = 6; first red-team review, Finding C); the slope remark was removed.
* **F16 applied.** Names qualified (`ChainQuotient.*`, `FiniteQuotient.closure_le_quotient_span`,
  `GeneralQuadraticClassification.degree_three_syzygy_zero`, `GrowthMain.growth_lfd_window`), with
  the convention stated; every name was checked against `Audit.lean`.
* **F17 applied.** `scripts/make_paper_figures.py` now redraws `curve-sweep.pdf` with the bound on a
  labelled right axis and asserts that no system is capped; the caption says so.
* **F18 partly applied.** Abstract and contribution list say "subgroup and coset bases" for the prime
  arm. Not changed: the title (frozen by the task) and the auto-generated header "subgroup
  (symmetric)" in `paper/curve-sweep-data.tex` (written by `scripts/build_ic_curve_paper.py`, which
  this pass could not edit).
* **F19 applied** (every item; Caminata-Gorla Theorem 2.9 and Example 4.3 checked in
  arXiv:2112.05579; the V_6 prediction is now "base depth 15, recursion from 14" per `d6-base-15`).
* **F20 applied.** `gmp21`, `kw19`, `hky15` (pages and DOI), `guo26` (full names), `ky` (v3) completed
  from Crossref and arXiv metadata. `silv` is kept: it is now cited for the Frobenius (Ch. II §2).
  `gaudry09` and `diem11` added from Crossref (see S1).
* **S1 resolved.** The ePrint 2004/031 abstract speaks only of index calculus with bounded solutions
  of summation-polynomial equations (the full PDF could not be fetched); the descent over extension
  fields is credited to Gaudry (JSC 44, 2009) and Diem (Compositio 147, 2011), the F_2-linear
  coordinate spaces to Faugère-Perret-Petit-Renault.
* **S6 resolved.** `ic-bounds-binary/results-m2.jsonl` has 30 LFD-5 cases with 28 distinct
  `system_hash` values, all at k = 10; the README's "30 systems" conflates cases and systems. The
  paper's 28 is correct and now also says "30 cases".
* **S2-S5 not checked** in this pass.
* **Prose.** All ten places rewritten or removed; run times, gigabytes, hardware and dates removed from
  the body and from Appendices A and B except the statement that two engines run on a graphics
  processor; engine names fixed to dense, shift, reference, third and bit-packed engine; notation
  clashes resolved locally (group algebra `\mathcal A = K[\Gamma]`, polynomial `h`, projector `\Pi`,
  lift subspace `\mathcal K`, multiplier `\lambda`, subgroup generator `\zeta`).
