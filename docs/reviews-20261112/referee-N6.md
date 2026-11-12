# Referee report: N6 "genericity" (`docs/proofs/N6-genericity.md`) and its use in `paper/curve-lfd-bounds.tex`

Referee: independent, adversarial. Repository read-only at commit `705279f`. Date 2026-11-12.
All computations below are my own (scripts in
`/tmp/claude-0/-home-user/c25f3f7b-b111-516d-a596-5c5b49ca2dae/scratchpad/n6/`; they import nothing
from the repository; the repository's engines were read for conventions only).

## 0. Verdict

| Claim | Where | Verdict |
|---|---|---|
| Genericity of the depth-7 base ("Cert passes off a proper closed set") | N6 Conjecture 6.2 | **Not proved, and N6 says so.** Correctly labelled a conjecture. No theorem of the paper depends on it. |
| Theorem 2.4 (specialisation of a parametric echelon basis) | N6 §2.3 | Proved as written. Its exceptional set `{Δ = 0}` is **not canonical** (finding 3). |
| Proposition 3.1, Corollary 3.2 (twist correspondence, H-invariance) | N6 §3.1 | Proved as written; identities re-derived and checked numerically. |
| Theorem 3.3 (twisted closedness of `W_k(S)`), Theorem 6.1 (all depths `k >= 8` on the `G_m`-orbit) | N6 §3.3, §6 | Proved, modulo N5 and the bridge (not my object) and one index slip (finding 7). The restriction to `k >= 8` is unnecessary: the same hypothesis gives `k = 7` (finding 4). |
| Proposition 2.2 (constructibility) | N6 §2.2 | Proved for the closure; the argument does not cover `Cert` as defined (existential in `S`) without finding 4. |
| Propositions 2.1, 2.3, 3.3 ("neither direction", "not open", "dimension not invariant") | N6 §2, §3.2 | The point values are correct (replicated). The inference "special > generic" and every use of "generic = 605/829" is **conditional on the unproved conjecture** (finding 2). Only "special < generic" is rigorous (on the line `b = t = 0`). |
| Section 4 sizes (`deg D >= 480`, 68,803 entries) | N6 §4.2 | Replicated at `k = 6` under the repository's column order. Labelled "exact" although it is an inference from samples, and it is a property of one echelon form, not of the exceptional set (finding 3). `k = 7` not replicated. |
| `thm:three` (`3t^2+a != 0`), `thm:carry` (`5t^3+3at+6b != 0`) | paper | Do not rely on N6. The stated conditions are exactly the ones the proofs use; both are sufficient, not necessary, and the paper states them as sufficient. Verified (finding 11). |
| Paper's summary of N6 (Section "What is missing", item 1) | paper l.1441–1458 | Faithful to N6 on the main point (genericity open). Four wording overclaims (findings 1, 2, 8, 9). |

No fatal gap in anything claimed as proved. The single most important point for a reader: **nothing in N6
or the paper proves a statement about generic `(a, b, t)`**; every lower bound above `LFD >= 4` is per
certified instance (and its scaling orbit).

## 1. Findings

### 1. Fixable gap (artifact): N6's own evidence is not in the repository
* Location: N6 header, §3.4, §4.2, §7; paper l.1445–1447, l.1456–1458, l.1622–1625.
* N6 cites `experiments/ic-bounds-genericity/` (README, `genericity_engine.py`, `twist_test.json`,
  `coset_test.json`, `extension_twist_test.json`, `line_degree_k6.json`, `line_degree_k7.json`) and
  `docs/IC_BOUNDS_PROGRAM.md`. None exists in the repository (`grep -r` finds the names only in N6
  itself). "item 2", "item 5", "N3", "N4a", "K4" are never defined in this repository.
* The paper says "every record behind the numbers in this paper are in the companion repository" and
  quotes N6 numbers as "exact": degree `>= 480` at `k = 6` and `7`, `10^10` coefficients, 52 CPU-hours.
  No record supports them here.
* Repair: commit the directory, or relabel these numbers. My replications (section 2) cover the `k = 6`
  degree data, depth-7 wrap-insensitivity, the twisted depth-8 check and the orbit statement at fresh
  points. Not replicated: `k = 7` line degrees, the 192 + 96 coset sweep, "four primes, identical support".

### 2. Fixable gap (logic): "generic value 605 / 829" is used as a fact; it is the conjecture
* Location: N6 Prop. 2.1 ("Generic `dim V_4 = 605`", label *exact*); §0 item 5 and §5 item 1
  ("C2 ⊂ E … by Corollary 2.5", "so `D` must vanish there by Theorem 2.4"); Cor. 2.5 parenthetical
  ("so `(A+3) | Δ`"); §4.2 ("the pivot set is the generic one at all 1200"); §2.4; paper l.1441–1443
  ("the closure dimension is semicontinuous in neither direction (… `607 > 605`, … `603 < 605`)").
* What is wrong: `dim V_4^K` over `K = F_p(A, B)` is not computed anywhere. N6 §2.4 itself explains that
  no finite set of points determines it. So "607 > generic" and Corollary 2.5 applied to C2
  (`857 > 829`) presuppose `d_gen = 605` resp. `829`, i.e. Conjecture 6.2. Over `F_193` every subset of
  the plane is Zariski closed, so the plane data alone cannot refute any semicontinuity statement.
* What is rigorous: on the line `b = t = 0` the generic value is `|M_{<=4}|` (kernel-checked
  `1 ∈ V_4` for `a^{2^k} != 1`, cofinite) and the special points `a ∈ H` have a proper ideal, hence
  `1 ∉ V_4`. That proves "special < generic" and "`{1 ∉ V_4}` is not open" (Prop. 2.3). There is no
  rigorous instance of "special > generic" in N6. "Not closed" in Prop. 2.3 is also conditional
  ("if the generic statement holds"), though the sentence ends unconditionally.
* Repair: state Prop. 2.1 and Cor. 2.5's applications as "relative to the value attained at 36,486 of
  37,249 points of the `F_193` plane", or as conditional on Conjecture 6.2. In the paper replace
  "is semicontinuous in neither direction" by "is not known to be semicontinuous, and the computed
  values 603 < 605 < 607 show that it is not constant; no argument transfers the value from points".
  The conclusion N6 draws (points cannot prove the generic facts) stands, for the simpler reason that
  no semicontinuity is available.

### 3. Unclear exposition / overclaim: `{Δ = 0}` and "the degree of the exceptional divisor" depend on the column order
* Location: N6 §0 item 3, §2.3 ("the exceptional set of the parametric certificate is `{Δ = 0}`",
  "equivalently the product … of the pivots met by any elimination"), §4.2 ("the exceptional divisor of
  the base does not grow with the depth", "settles question 2 negatively"), §8 items 1 and 4;
  paper l.1444–1447.
* What is wrong. `Δ` is the common denominator of the reduced echelon basis for one ordering of the
  monomials inside each weight class. A different weight-sorted order gives other pivots and another
  `Δ`. Measured (section 2, item R6): with the repository's order (descending mask) I reproduce N6
  exactly at `k = 6` (lcm degree 480, numerators up to 538, 68,801 nonzero entries). With ascending
  `(e_1, e_2)` inside each weight class, the same subspace has 79,169 nonzero entries, 43 of 192 sampled
  weight-4 entries have `deg num + deg den >= 1260` (not reconstructible from 1260 points), and the lcm
  of the 149 reconstructed denominators already has degree 667. So "480" is a property of one basis.
  Also "equivalently the product of the pivots" is false: the common denominator divides a product of
  pivots; the zero sets differ in general.
* The canonical statement. `V^K` is a `K`-point of a Grassmannian, i.e. a rational map from the
  `(A, B)`-plane, defined outside finitely many points (codimension 2 on a smooth surface). Where it is
  defined and the limit space `W_P` meets each `M_{<=e}`, `e <= 3`, in the generic dimension,
  `Gen(P) ⊆ W_P`, closedness and `1 ∉ W_P` hold by continuity (containment and closedness are closed
  conditions inside the open Schubert stratum). So the canonical exceptional set is
  `{indeterminacy points} ∪ {P : dim(W_P ∩ M_{<=e}) jumps for some e <= 3}`, contained in `{Δ = 0}` for
  every admissible order, and possibly of far smaller degree. N6's "explicit `E` is out of reach" is
  shown only for one presentation of an over-approximation.
* Repair: define `E` canonically; state 480 as "for the echelon form in the order of
  `csd/ic_bounds_fall.monomials`"; drop "equivalently the product of the pivots".
* Also (label): `deg D >= 480` is tagged "exact". It is an inference: a rational function of degree
  `(538, 480)` matches 1200 values. If the true entry were different it would have
  `deg num >= 720` or `deg den >= 662`. So "large" is safe, "denominator `>= 480`" is not strictly
  proved, and "the pivot set is the generic one" is circular.

### 4. Misleading remark, and a strengthening N6 missed: `Cert` is a depth-7 statement and gives `k = 7`
* Location: N6 Def. 1.4, §3.3 Remark "why `k >= 8` and not 7", §3.4 first bullet, §4.1 ("for a theorem
  that also covers `k = 7` one adds the `14 × 115` depth-7 identities"), Prop. 2.2.
* Claim (short proof). Let `v ∈ S ∩ M_{<=3}` and `u = 2^j e_i`, `j <= 6`. In `A_8`,
  `τ_u v = v_0 + τ_{128 e_i} v_1` with `v_0, v_1` supported on exponents `< 128`. Closedness of `W_8(S)`
  and uniqueness of the high decomposition give `v_0 ∈ S` and `v_1 ∈ S ∩ M_{<=3}` (this is the step N6
  uses in Theorem 3.3). In `A_7^{(c)}`, `τ^{(c)}_u v = v_0 + c·v_1`. Hence `S` is `τ^{(c)}`-closed at
  depth 7 for **every** `c`, so `V_4^{(c)}(7)(P) ⊆ S`, and by Prop. 3.1(c)
  `1 ∉ V_4(7)(μ·P)`, `dim V_4(7)(μ·P) <= dim S` for every `μ`. Conversely, if `S ⊇ Gen_7`,
  `1 ∉ S` and `S` is closed under `τ^{(0)}` and `τ^{(1)}` at depth 7, then `W_8(S)` is closed (shifts at
  level 7 never need closedness). Therefore
  **`Cert(P)` holds iff `1 ∉ S_min(P)`, the least subspace of `M_{<=4} ⊆ A_7` containing `g` and closed
  under both `τ^{(0)}` and `τ^{(1)}`**, and Theorems 3.3 and 6.1 hold for `k >= 7`.
* Consequences. (a) The Remark is wrong in spirit: the hypothesis already contains depth-7 closedness in
  every twisted algebra. (b) "Depth 7 is wrap-insensitive" (§3.4) is not an extra exact datum on three
  instances: `V_4^{(c)}(7) ⊆ V_4(7)` follows from `Cert` with `S = V_4(7)`. (c) The least `S` exists
  (intersections of admissible `S` are admissible), which is what Prop. 2.2 needs to call `Cert`
  constructible; its proof as written covers only `S = V_4(7)`. (d) The Lean-facing certificate needs no
  depth-8 object.
* Checked numerically (R4): at six points `W_8(S_min)` is closed for `c ∈ {1, 0, 2}`, and
  `S_min = V_4(7)` whenever `W_8(V_4(7))` is closed.

### 5. New facts that settle two items N6 leaves open, and narrow Conjecture 6.2
* Location: N6 §5 item 1 ("whether the free extension of the 857-dimensional space is closed is a finite
  check per point, not done"), §8 item 2 ("decides whether C1, the tangent line and the special points
  pass"), Conjecture 6.2.
* `Cert` **fails for every `S`** at the C2-class target `p = 257`, `(a, b, t) = (148, 227, 143)`, `k = 7`:
  `dim V_4(7) = 857`, the spaces `V_4^{(c)}(7)` differ for `c = 0, 1, 2, 5`, `W_8(V_4(7))` is not closed,
  and `S_min = M_{<=4}` contains `1`. Yet at depth 8 (an honest prime-field instance, `H = F_257^*`)
  `1 ∉ V_4(8)`, `dim V_4(8) = 1093 = 857 + 2·117 + 2`, `|Z| = 275`, so `LFD >= 5` there. At `(3, 0, 116)` the depth-8 closure
  also has dimension 1093 and excludes `1` (`|Z| = 286`; `Cert` not run there). So `Cert` is strictly
  sufficient, not only `{Δ != 0}`.
* `Cert` **passes** at points that Conjecture 6.2 places in `E`: tangent `(148, 227, 76)` with profile
  `(0,0,3,117,828)`; order two `(0, 1, 256)` with `(0,0,2,117,835)`; C1 `(0, 1, 116)` with
  `(0,0,2,115,826)`; `t = 0` `(148, 227, 0)` with `(0,0,3,123,849)` (all `p = 257`). Theorem 6.1 then
  gives `LFD >= 5` for `k >= 7` on their orbits with the corresponding bounds
  (e.g. `828 + 234(k-7) + 3·C(2k-14, 2)`). One point per class; not a statement about the classes.
* Repair: record these; state Conjecture 6.2 for `Cert` with unspecified profile, and say that C2 is where
  the method (not only the generic profile) fails.

### 6. Unclear exposition: which `(p, k)` are actual prime-field instances; small primes
* Location: N6 §0 item 1, Theorem 6.1 ("`p` odd, `k >= 8`"), §8 item 2; paper l.1455–1458.
* A prime-field subgroup chain of depth `k` needs `2^k | p - 1`. At `p = 193` (`192 = 2^6·3`) there is no
  depth-7 or depth-8 prime-field instance, and none at `p = 11`. "A complete exact statement for that
  prime at every depth `k >= 8`" is a statement about systems over `F_{193^e}` whose curve and target are
  defined over `F_193`. N6 §1.1 says this once; §8 and the paper's "third" conclusion do not.
* For small `p` a dense open set need not have many `F_p`-points. Measured at `k = 6`, `t = 1`, all
  `(A, B)`: `p = 3`: 1 of 9 points has the profile `(0,0,2,111,605)`; `p = 5`: 8 of 25. N6 §5 item 3 says
  the same in principle (`deg Δ > p` for 193, 257, 449, 577).
* Repair: add the hypothesis `2^k | p^e - 1` and the phrase "with `a, b, t ∈ F_p`" wherever "every depth"
  is used for a fixed prime.

### 7. Typo with content: `wt h <= 2` in the restatement of N5
* Location: N6 §1.2 ("`h` supported on levels `>= 7`, `wt h <= 2`") against Theorem 6.1's bound, which
  has the terms `C(2k-14, 3) s_1 + C(2k-14, 4) s_0`.
* For `k >= 9` the proof of Theorem 3.3 raises `wt h` from 2 to 3; with patterns of weight `<= 2` only,
  `W_k(S)` is closed only when `S ∩ M_{<=1} = 0`. Use all `h` with `wt h <= 4` (as in
  `docs/proofs/N5-reduction.md`), or assume `s_0 = s_1 = 0`.

### 8. Overclaim in the paper: "the exceptional set is known completely" at `k = 6`, `p = 193`
* Location: paper l.1430–1437; `experiments/ic-bounds-locus/README.md` §3.
* What the data contain (re-analysed independently, R7): 1276 targets on 13 curves with certified exact
  LFD, and all 37,249 points of the sheet `t = 1` with `D <= 4` only, run with `--no-certify`
  (`certified: false` on every row). The rule "LFD 4 iff `t = 0, ab = 0` or `(a/t^2, b/t^3) ∈ C2`" has 0
  violations on the 1276 targets. But: (i) the plane gives saturation of `V_4`, not exact LFD;
  (ii) the transfer from the sheet `t = 1` to other cosets is empirical, and N6 Prop. 3.3 proves the
  dimension is not invariant; off the sheet, 768 of 74,498 `(A, B, coset)` classes were computed;
  (iii) the plane `t = 0` is represented by 13 curves.
* "C2 lies on no curve of degree at most 18" carries no information: there are exactly 190 monomials of
  degree `<= 18`, so 190 points in general position lie on no such curve, and every 190 points lie on a
  20-dimensional family of degree-19 curves (rank 190 in both cases, checked). The README's conjecture
  "C2 is an absolutely irreducible curve of degree `>= 19`" has no support; N6 §5 item 3 agrees.
* Repair: "is known on the 1276 targets computed and, for the degree-4 closure, on the whole sheet
  `t = 1`".

### 9. Paper wording: "generic" is not defined
* Location: abstract ("generic targets"), l.800 table, l.373 growth conjecture ("generic two-summand
  prime chain systems"), l.1090–1093 (`d_D(k)` "attained on a Zariski-open set of coefficient vectors of
  `g`, if such a value exists"), l.1096, l.1191.
* The parameter space is `(a, b, t)` (weighted, 2-dimensional after scaling plus the coset), not the
  coefficient vectors of `g`. A generic value always exists over the algebraic closure
  (constructibility), so "if such a value exists" hedges the wrong thing; what is unknown is its value
  and whether the open set has `F_p`-points for the `p` in use. The exceptional set depends on `k`
  (C2: saturated at `k = 6`, 857 at `k = 7`) and on `p`. Every use of "generic" in the paper is in fact
  "random curve and random target of large order, on the instances computed". The paper mostly says
  this ("on the instances where it was determined", "instances we regard as generic"); the table entry
  "5 (`k = 6, …, 9`, generic targets)" and the growth conjecture do not.
* l.1444: "what is needed is the generic depth-7 basis" should be "what would suffice"; finding 4 gives a
  smaller sufficient object (`S_min^K` excluding `1`).

### 10. Minor: Conjecture 6.2 "Equivalently (Theorem 2.4)"
* Theorem 2.4 yields `S = W_ρ ⊇ V_4(7)(P)`, possibly strict (N6's own second remark). To get
  `S = V_4(7)(P)` with profile `(0,0,2,115,829)` one must also avoid the locus where the specialised
  closure is smaller. Either say `S = W_ρ` or enlarge `Δ`. The converse direction needs Prop. 2.2.

### 11. `thm:three` and `thm:carry`: conditions are the ones needed (no finding against the theorems)
* `α - β^2 = -2(3t^2 + a)` for `α = -2(t^2+a)`, `β = -2t`: verified symbolically.
* Coefficient of `x_1^3 x_2^3 x_3^3` in `S_4 = Res_X(S_3(x_1,x_2,X), S_3(x_3,t,X))` is
  `-8(5t^3 + 3at + 6b)`; `S_4` has 125 terms, degree 4 in each variable, and `(3,3,3)` is the only
  exponent of binary weight 6: verified symbolically. The step `V_d ⊆ (⟨G⟩ ∩ M_{<=d}) + K g` is correct
  (a leading monomial of an ideal element is non-standard, so top parts cannot cancel against reduced `g`).
* Brute force of `thm:three` in the ordinary polynomial ring (own code, chains included, no quotient):
  all `(a, b, t) ∈ F_p^3` for `k = 2`, `p ∈ {3,5,7,11,13,17}` and `k = 3`, `p ∈ {3,5,7}` (9,431
  systems): 0 systems with `1 ∈ V_3` and `3t^2 + a != 0`. The condition is not necessary (e.g.
  `p = 17, k = 2`: 177 of 289 tangent systems have `1 ∈ V_3`); the paper claims sufficiency only.
* Neither theorem cites or needs N6. The abstract and l.526–528 state them correctly ("unit ideal and
  `3t^2+a != 0`"; "`k >= 3` … `5t^3+3at+6b != 0`").
* Typo: paper l.718–719, "separates ordinary / ordinary first fall degree" (word repeated across the line break).

### 12. Minor
* N6 §0/§6: "the universal bound `|Z| <= 2^{k+1}`" needs the target not to be a singular point:
  `S_3(t, x_2, t)` vanishes identically iff `t^3+at+b = 0` and `t^4 - 2at^2 - 8bt + a^2 = 0`. At the cusp
  point `(A, B) = (-3, 2)` the plane data have `|Z| = 127` (the two full lines), so the stated bound still
  holds there, but not by the two-roots argument. N4a is not in this repository; unverified.
* N6 Prop. 2.1 writes "`low_3 = 111`, `low_2 = 1` at degree 3", mixing the degree-4 and degree-3 closures.
* N6 §4.2: 68,803 nonzero entries "at two points over `p = 40961` and `p = 33554393`"; I get 68,801 at two
  points over `p = 40961` and 68,799 at one over `p = 12289`. Within the variation N6 itself describes.
* "Orbit of `p - 1` scalings" (paper l.1454): `p - 1` points only when `t != 0`.

## 2. Steps verified, and how

R1. `g` of N6 §1.1 equals Semaev's `S_3` (sympy); vanishes on `x(P), x(Q), x(P+Q)` on random curves.
R2. Weighted homogeneity, Prop. 3.1(a) (sympy). `S_3 = q^2` on `t^3+at+b = 0` (sympy).
R3. Prop. 3.1(b),(c) by hand; and numerically at `p = 193`, `k = 6`: `φ_μ(V_4^{(c)}(P)) = V_4(μ·P)` as
    subspaces at four points (generic, `j = 1728`, C2, tangent), `c = μ^{-64}`.
R4. Theorem 3.3 / 6.1, by hand (all cases of the carry argument), and numerically with my own engine
    (group algebra, twisted shifts, exact arithmetic mod `p`):
    * `p = 257`, `(88, 32, 5)` and `p = 12289`, `(4321, 987, 55)`: `dim V_4(7) = 829`, profile
      `(0,0,2,115,829)`, same subspace for `c ∈ {1,0,2,5}`, `W_8` of dimension 1061 closed for
      `c ∈ {1,0,2}`, `1 ∉`.
    * `p = 12289`: `V_4(8)(P) = W_8` (1061, profile `(0,0,2,119,1061)`); for `μ = 3, 7, 11`
      (`μ^256 != 1`): `dim V_4(8)(μ·P) = 1061`, `1 ∉`.
R5. Point values: `p = 193`, `k = 6`: 605, 607, 603, 602, 600, 599, 681; `(13,0,0)`: 730, `1 ∉`;
    `(52,0,0)`: 794, `1 ∈`; `(20,0,1)`: 784, `|Z| = 10`; `(5,0,96)`: 762, `|Z| = 32`. 29 plane rows
    (all classes) re-computed: 29 of 29 agree in `dim`, `low_4`, `|Z|`, `1 ∈`. Depth 6 with `c = 0, 2`:
    same dimension and pivots, neither inclusion; the joint `(τ^{(0)}, τ^{(1)})` closure at depth 6 is
    all of `M_{<=4}`.
R6. Line degrees at `k = 6`: 1260 closures along my own line over `p = 40961`, 288 entries, rational
    reconstruction by Newton interpolation and extended Euclid, every reconstruction re-evaluated.
    Repository order: all 288 reconstructed, slack `>= 241`; denominators: weight 2 `<= 1`, weight 3
    `<= 112`, weight 4 `<= 480` (30 of 192 at `>= 400`); numerators `<= 538`; lcm degree 480; lcm
    vanishes at the order-two and tangent points of the line. This matches N6 §4.2. Other order: see
    finding 3.
R7. Locus data re-analysed without the repository's summaries: 37,249 plane rows; classes are exactly
    `1+A+B = 0` (192 + cusp), `A = -3` (191 + 2), the cubic C1 (186 + 4 special) and 190 saturated
    points off all three; 1276 targets, 12 with LFD 4, rule has 0 violations; target class equals plane
    class at the normalized point for all 1268 targets with `t != 0`.
R8. Counts: 106 generator rows, 119 low rows of `W_8`, 1904 identities, 1471 and 2517 columns,
    `4k^2 + 172k - 571`, the table `2005, 3024, 4148, 5013, 4759` against `1061, …, 2069`, 320, 192, 96.
R9. Theorem 2.4 and Lemmas 1.1, 1.3, column-order paragraph, §4.3 spanning claim: by hand, correct.

## 3. Not verified

* N5 and the bridge themselves (kernel-checked per the paper; outside this review).
* N6 §3.4 coset sweep (192 + 96), "identical support across four primes", `line_degree_k7` (480 at
  `k = 7`), the `k = 7` entry counts (93,011; 8,893), timings and cost estimates of §8.
* The closed form `w = Q_w g` with denominator `a + 3t^2` (not given in N6); I only see that the
  weight-2 denominators have degree `<= 1` and vanish at `A = -3` on my line.
* "K4" (kernel-checked `1 ∈ V_4` on `b = t = 0`, `a^{2^k} != 1`): taken from the paper; my engine
  agrees at one point.
* The `p = 11` and `p = 40961` Lean instances.

## 4. Paper versus N6: mismatches

1. l.1441–1443 "semicontinuous in neither direction": conditional on the unproved generic value
   (finding 2).
2. l.1444–1447 "what is needed … common denominator has degree at least 480 … (exact, `k = 6` and 7)":
   sufficient, not needed; order-dependent; no record in the repository; `k = 7` unreplicated
   (findings 1, 3, 9).
3. l.1430 "known completely": finding 8.
4. l.1455–1458 per-prime description at `p = 193`: not about prime-field instances of depth `>= 8`
   (finding 6).
5. l.1447–1455 (scaling orbit, twisted closedness, unit and proper ideals): faithful to N6, and proved.
6. l.1049 "What remains open is the genericity of the depth-7 base in `(a,b,t)`": accurate.
7. No theorem, proposition or lemma of the paper marked "written" or "kernel-checked" cites N6 or
   depends on genericity. `thm:three`, `thm:carry`, `thm:reduction`, `thm:selector` are independent of it.

## 5. Scripts (all in `/tmp/claude-0/-home-user/c25f3f7b-b111-516d-a596-5c5b49ca2dae/scratchpad/n6/`)

`sym.py` (identities), `eng.py` (closure engine), `cert.py` (depth-7/8 certificate, `S_min`, orbit),
`twist.py` (Prop. 3.1(c)), `d6.py` (depth 6, C2 rank), `spot.py` (plane spot check), `locus.py`
(locus data), `three.py` (`thm:three` brute force), `tiny.py` (`p = 3, 5`), `line.py` + `recon.py`
(line degrees; `ORDER=repo` for the repository's column order; samples in `line_part*.pkl` and
`ord_mine/`), `nnz.py`. Outputs: `cert_*.out`, `three.out`.
