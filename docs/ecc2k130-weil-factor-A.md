# The 130-dimensional Weil-restriction factor A of ECC2K-130: exact invariants

Date: 2026-11-23. Follows the transfer survey (`docs/ecc2k130-transfer-survey.md`, PR #47), which left
one open route: is `A`, the 130-dimensional factor of the Weil restriction of ECC2K-130, isogenous to
the Jacobian of a curve of genus about 130 to 300, by a computable map? This note computes `A`'s exact
invariants directly (not by citation), as the first rigorous step. Reproduce with
`python3 experiments/ecc2k130-weil-factor/compute_A.py`; the polynomial is in `P_A.json`.

## Setup

`E / F_2 : y^2 + xy = x^3 + 1` (the Koblitz curve `a = 0, b = 1`), `#E(F_2) = 4`, Frobenius
`pi^2 + pi + 2 = 0` (trace `t = -1`), so `E` has CM by `Q(sqrt(-7))`. The Weil restriction
`W = Res_{F_2^131 / F_2} E` has dimension 131 and splits `W ~ E x A` with `dim A = 130`. The
2-Frobenius eigenvalues of `W` are `{alpha zeta^j, beta zeta^j : j = 0..130}`, `zeta` a primitive
131st root of unity; `j = 0` is `E`, the other 130 pairs are `A`.

## Exact result

The characteristic polynomial of 2-Frobenius on `W` is
`P_W(x) = (x^131 - alpha^131)(x^131 - beta^131) = x^262 - s_131 x^131 + 2^131`,
with `s_131 = alpha^131 + beta^131 = -22283658519494248867`. Dividing by `E`'s Weil polynomial
`x^2 + x + 2` (exact, remainder zero) gives the degree-260 integer polynomial `P_A`.

Everything below is computed exactly and cross-checked:

| Invariant | Value / verdict |
|---|---|
| `#E(F_2^131) = #W(F_2)` | `2722258935367507707729280517973639940516`, cofactor exactly 4 |
| `#A(F_2) = P_A(1)` | `680564733841876926932320129493409985129` — **prime, 130 bits** |
| `deg P_A` | 260 (= 2·dim) |
| constant term | `2^130` (functional equation) |
| middle coefficient | **odd ⇒ `A` is ordinary** (p-rank 130 = dim) |
| `P_A` squarefree | yes (`gcd(P_A, P_A') = 1`) |
| CM field | `K = Q(zeta_131, sqrt(-7))`, degree 260 (`sqrt(-7) ∉ Q(zeta_131)` since `7 ∤ 131`) |
| `P_A` mod 3, 5, 13 | two irreducible factors of degree 130 each |

**`A` is simple** (two independent proofs):
- *Direct* (not Honda-Tate by citation): the only `sigma ∈ Gal(K/Q)` that could fix `alpha·zeta` would
  need `alpha/beta = (-3 - sqrt(-7))/4` to be a 131st root of unity; it has trace `-3/2`, so it is not an
  algebraic integer, hence not a root of unity. So `alpha·zeta` has trivial stabiliser,
  `[Q(alpha·zeta):Q] = 260 = deg P_A`, and `P_A` is irreducible.
- *By CAS*: `python-flint` factors `P_A` over `Z` into a single irreducible factor of degree 260.

(Consistent with the mod-p pattern: `Gal(K/Q) ≅ C_130 x C_2` is non-cyclic, so `P_A` never factors into
a single piece mod `p`, but always into equal-degree pieces — degree 130 at `p = 3, 5, 13`.)

`P_A` is a genuine Weil polynomial for `q = 2`: all 260 complex roots satisfy `|root| = sqrt(2)`
(checked numerically), and the coefficients obey the functional equation `a_i = 2^(g-i) a_{2g-i}`
(`g = 130`), verified exactly.

### What this settles

- The easiest hope — that `A` splits into a smaller transfer target — is **false**: `A` is simple.
- `#A(F_2) = r`, the 130-bit ECDLP prime, so `A` genuinely carries the discrete logarithm: a group
  isomorphism `A(F_2) ≅ E(F_2^131)[r]` follows from the Weil-restriction structure. Any Jacobian
  search must match `P_A` exactly; the polynomial is now pinned (`P_A.json`).
- `A` is **ordinary**, so there is no supersingular/p-rank obstruction to being a Jacobian.

## Where the window question stands

`A` carries a rational order-131 automorphism `rho` (the `Gal(F_2^131 / F_2)` action on the Weil
restriction, defined over `F_2`). This is **consistent** with `A` being a Jacobian factor: a genus-130
cyclic cover of `P^1` of degree 131 exists by Riemann–Hurwitz (`2·130 - 2 = 131·(-2) + 130·r` gives
`r = 4` branch points), and its Jacobian carries exactly such a `rho`. So `rho` is not by itself an
obstruction — `A` even has one.

The open, discriminating question is therefore a **CM-type match**: does the CM type of `A` (which 130
of the 260 embeddings `K -> C` are the Frobenius half) coincide with the CM type of a degree-131 cyclic
cover's Jacobian factor (fixed by the branch exponents `a_1 + ... + a_4 ≡ 0 mod 131` via
Shimura–Taniyama)? `A`'s CM field `Q(zeta_131, sqrt(-7))` strictly contains `Q(zeta_131)`, so the extra
`sqrt(-7)` (from `E`'s CM) must also be produced by the cover — which is exactly the constraint that, in
`crypto/research/.../RESEARCH_ECC2K130_CURVE_CONSTRUCTION.md`, pushes any carrier with an order-131
symmetry to genus ≥ 394, outside the usable window. Reconciling that bound against this CM-type
computation is the next concrete step; it is a finite check, not a run.

## Jacobian necessary-condition tests (passed)

If `A` were isogenous to the Jacobian of a genus-130 curve `C / F_2`, then `C`'s zeta numerator is the
reciprocal of `P_A`, and both the point counts `#C(F_{2^k}) = 2^k + 1 - p_k` (`p_k` = power sums of
`P_A`'s roots) and the closed-place counts `a_d = (1/d) sum_{e|d} mu(d/e) #C(F_{2^e})` must be
non-negative integers. A single negative value would prove no curve has this zeta function, hence `A`
is not a Jacobian. Computed in `experiments/ecc2k130-weil-factor/jacobian_tests.py` (`p_k` by Newton's
identities, cross-checked against the exact `p_k = -s_k^E` for `131 ∤ k`, `130 s_k^E` for `131 | k`):

- `#C(F_2) = 2`, `#C(F_4) = 2`, `#C(F_8) = 14`, `#C(F_16) = 18`; **all `#C(F_{2^k}) ≥ 0`** (min 2,
  checked `k = 1..1599`).
- closed-place counts `a_1 = 2, a_2 = 0, a_3 = 4, a_4 = 4, ...`; **all `a_d ≥ 0` integers** (min 0,
  checked `d = 1..400`).

So **`A` passes every cheap decisive test**: its virtual zeta function is that of a genuine genus-130
curve candidate. There is no point-count or place-count obstruction to `A` being a Jacobian.

## What is at stake: index-calculus cost if A is a Jacobian

To make the window question's payoff precise, `experiments/ecc2k130-weil-factor/ic_cost_model.py`
computes the standard smooth-divisor index-calculus cost for the *hypothetical* genus-130 curve `C`
(whose place counts `a_d` and effective-divisor counts `B_m` are fixed by `P_A`). Factor base = places
of degree `<= t`; a relation is a `t`-smooth reduced divisor of degree `~ g`; cost = `N(t)/P_smooth`
trials plus sparse linear algebra `~ N(t) g`.

- Effective-divisor counts `B_m` are all `>= 0` (a further consistency check).
- Optimal `t = 24`, `|FB| ~ 1.46e6`: relation collection `~ 2^33` trials (`~ 2^47` operations with a
  conservative `g^2 ~ 2^14` per-trial smoothness-test factor), sparse linear algebra `~ 2^27.5`.
- So, **conditional on `A ~ Jac(C)`, ECDLP on ECC2K-130 would cost about `2^33`-`2^47`, i.e. 14 to 28
  bits below the matched Pollard rho reference (`~ 2^61`).** This refines the `2^37.17` figure in the
  crypto notes with exact counts.

Two hard caveats: (1) this is **conditional** on the open window question; (2) even if `A ~ Jac(C)`,
running index calculus needs the *explicit* curve `C` and an *evaluatable* isogeny `A -> Jac(C)`, which
existence alone does not provide. This is a cost model, not a validated run (no replay/validation-gate
evidence), and claims nothing about `A` being a Jacobian.

## CORRECTION / RETRACTION (2026-11-24, red-team catch)

An independent red-team review found a genuine error in the `sqrt(-7)` argument used below and in
`not_a_jacobian.py`. **It is RETRACTED.** The error: I claimed the Frobenius eigenvalues of a `mu_131`
cyclic cover of `P^1` lie in `Q(zeta_131)`. That holds only for 1-dimensional character eigenspaces
(3 branch points). With 4 branch points (the genus-130 case, forced by Riemann-Hurwitz), each nontrivial
character eigenspace is **2-dimensional** (Chevalley-Weil), Frobenius acts as a 2x2 block, and only the
*characteristic polynomial coefficients* (trace, determinant) are forced into `Q(zeta_131)` -- the
eigenvalues may lie in a quadratic extension. And `A` satisfies exactly this: its eigenvalue pairs
`{alpha^130 zeta^-m, beta^130 zeta^-m}` have `sum = s_130 zeta^-m` and `prod = 2^130 zeta^-2m` with
`s_130 = alpha^130+beta^130 = 73000192456347156057` an integer, both in `Q(zeta_131)`; the discriminant
is `-7 (c zeta^-m)^2`, a non-square in `Q(zeta_131)`, so the eigenvalues carry `sqrt(-7)` -- perfectly
compatible with a 2-dim eigenspace of a 4-branch-point `P^1` cover. **So the `sqrt(-7)` "obstruction" is
illusory**, the P^1-cover exclusion and the genus `>= 458`/`917 | N` bound below are NOT established,
and the not-a-Jacobian theorem's step 5 fails. `A` being the Jacobian of a 4-branch-point cyclic cover
of `P^1` is NOT ruled out by the field argument -- it is attack-relevant and open.

A possible *eigenvalue-independent* rescue (class field theory): an `F_2`-rational order-131 automorphism
makes the deck group the constant group scheme `Z/131`, and a constant-`Z/131` cover of `P^1` with branch
total degree 4 needs `131 | |Cl^0_m(P^1)|`, but `|Cl^0_m|` is coprime to 131 for every degree-4 branch
configuration (it needs a degree-130 place, `ord_131(2) = 130`). If sound, this would re-close the
isomorphism case -- but it is a CANDIDATE pending independent review, not a claim, and it does not touch
the isogeny case. The ray-class cover-of-E result (`rayclass_cover_E.py`) uses the same mechanism and is
likewise downgraded to "candidate, pending re-review". The decisive test the red-team recommends: compute
the 2x2 Frobenius block of an explicit 4-branch-point cover via Jacobi sums and check whether any matches
`A`'s `{alpha^130 zeta^-m, beta^130 zeta^-m}`.

---

## (RETRACTED) Sub-result: cyclic covers OF P^1 cannot realize A in-window

**Scope correction (2026-11-24).** An earlier version of this note claimed this result closed "every
construction with abelian symmetry". That was an overstatement: the argument below rules out cyclic /
abelian covers of **`P^1`** only. Covers of `E` itself (the CM base) escape it entirely and are treated
in the next section, where the route is *reopened*.

The standard way to build a Jacobian with prescribed CM is a cyclic (superelliptic) cover of `P^1`;
this is the family the crypto notes' genus >= 394 bound targets. For covers of `P^1` it is **rigorously
excluded** (`experiments/ecc2k130-weil-factor/cyclic_cover_exclusion.py`):

- *Mechanism.* For an abelian cover `C -> P^1` with Galois group of exponent `N` (char 2, `N` odd),
  every Frobenius eigenvalue on `Jac(C)` is a Jacobi sum, so lies in `Q(zeta_N)` (Weil). If `A` were an
  isogeny factor, its eigenvalues would lie in `Q(zeta_N)`, forcing the CM field
  `K = Q(zeta_131, sqrt(-7))` (which `A`'s eigenvalues generate, `A` being simple) into `Q(zeta_N)`,
  i.e. `conductor(K) = lcm(131, 7) = 917 | N`.
- *The obstruction is `sqrt(-7)`.* `A`'s base-change eigenvalues to `F_{2^130}` are
  `alpha^130 * zeta_131^{-j}` with `alpha = (-1 + sqrt(-7))/2`. Computed exactly,
  `alpha^130 = (73000192456347156057 - 4061839345336951189 sqrt(-7))/2`, which has a **nonzero
  `sqrt(-7)` component**, so `alpha^130 notin Q(zeta_131)` (because `Q(sqrt(-7)) cap Q(zeta_131) = Q`,
  as `7 nmid 131`). A pure `Z/131` cover has all eigenvalues in `Q(zeta_131)` and therefore cannot carry
  `A`. The `sqrt(-7)` that `A` inherits from `E`'s CM is exactly what a `Z/131` cover lacks.
- *Genus consequence.* Supplying `sqrt(-7)` needs `7 | N` (since `sqrt(-7) in Q(zeta_7)`, `7 ≡ 3 mod 4`),
  so `917 | N`. A tame `Z/917` cover of `P^1` has genus `>= 458` (`r = 3` branch points; more branch
  points only raise it) -- far above the usable window (`~130`-`300`). So **no cyclic cover realizes `A`
  at usable genus.** This gives a clean mechanism for, and sharpens, the crypto notes' `>= 394` bound.

For contrast, a `Z/131` cover with 4 branch points has genus exactly 130 (in-window) -- but its
eigenvalues live in `Q(zeta_131)`, so it carries a *different* (no-`sqrt(-7)`) abelian variety, not `A`.

This closes the `P^1`-based cyclic/abelian covers. It does **not** close covers of `E`, nor the fully
general case.

## Cover of E / Prym route: CANDIDATE closure (class field theory), PENDING RE-REVIEW

The natural carrier for `A` is a cyclic cover of `E` itself, not of `P^1`. A degree-131 cover
`C -> E` totally ramified at a degree-2 conductor has genus 131 (`2g-2 = 131(2*1-2) + 130*2 = 260`),
in-window, with `Jac(C) ~ E x Prym`, `dim Prym = 130`. For `Prym ~ A` we need `Jac(C) ~ E x A = W`,
where `W = Res_{F_2^131/F_2}(E)` has `P_W = (x^2+x+2) P_A = x^262 - s_131 x^131 + 2^131`. The `P^1` field
obstruction does **not** apply here (`E` already has CM by `sqrt(-7)`), and `W` passes the place-count
necessary conditions (`cover_of_E_reduction.py`): all `a_d >= 0`, and strikingly `#C(F_2^k) = 2^k + 1`
for all `k = 1..130` (the virtual curve is `P^1`-like up to degree 130). So the cheap tests do not close
it.

But **class field theory does** (`experiments/ecc2k130-weil-factor/rayclass_cover_E.py`). A
geometrically connected cyclic degree-131 cover of `E` with conductor `m` corresponds to an order-131
character of the ray class group `Cl_m(E)` that is *nontrivial on the degree-0 part* `Cl^0_m`
(characters through the degree map give the disconnected constant-field extension `Res(E)` itself, not a
connected curve). Since `E` has a degree-1 place, `Cl_m = Z (+) Cl^0_m`, so such a character exists
**iff `131 | |Cl^0_m|`**, and `|Cl^0_m| = #E(F_2) * |ray part|`, `#E(F_2) = 4`, ray part (tame, since
`131` is coprime to char 2) `= prod_{v|m}(2^deg(v) - 1)` up to a 2-power.

`A` as the Prym forces `g_C = 131`, hence conductor of total degree 2. Both cases give
`|Cl^0_m| in {4, 12}`, **coprime to 131** -- so there is no such character and **no connected cyclic
degree-131 cover of `E` has genus 131: `A` is not its Prym.** 131-torsion in the ray class group first
appears when the conductor contains a degree-130 place (`ord_131(2) = 130`, so `131 | 2^d - 1` iff
`130 | d`), which forces genus `>= 8451` -- astronomically outside the window. (A connected degree-131
cover is Galois iff cyclic since 131 is prime, and `A`'s `zeta_131`-eigenvalue structure requires the
cyclic one.)

So the cover-of-E route is a **candidate** closure (pending review -- same mechanism as the retracted-and-
candidate P^1 ray-class rescue), by the same `ord_131(2) = 130` prime-degree obstruction that
governs every other level of this problem -- now realized as the ray class group lacking 131-torsion
below a degree-130 conductor. Correcting the record: an earlier version called this route "open"; class
field theory closes it. What remains genuinely open is only the fully general case -- `A` as the
Jacobian of a curve with *no* cover structure or useful symmetry -- which has no known decision
procedure at dimension 130 and, lacking any construction, carries no attack.

**Literature grounding.** Scholten realizes `Res_{L/K}(E)` as a genus-2 Jacobian for quadratic `L/K`;
genus-3 cover attacks do the analogue for cubic extensions. But Freeman (ePrint 2009/103) notes the Weil
restriction is *generally not itself a Jacobian* for higher degree -- one must find a curve whose
Jacobian is isogenous to it -- and the converse direction (a curve with a nonconstant map to `E` has
Jacobian isogenous to `Res(E)`) is exactly our reduction. The standard GHS cover of `P^1` for ECC2K-130
is forced to genus `~2^130` or the zero map (the crypto survey), so GHS does not produce the genus-131
object; whether a *non-GHS* degree-131 cover of `E` does is the open window question, right at the genus
floor (`>= 130`).

**Why no cheap obstruction (structure).** The Prym eigenvalues come from *ramified* order-131 twists of
`E`: a nontrivial *unramified* rank-1 local system `L` on the genus-1 curve `E` has `H^1(E, L) = 0`
(`chi = rank * (2 - 2g) = 0`), so it contributes nothing; a twist ramified at 2 points has
`chi_c = -2`, contributing exactly 2 eigenvalues per character, `260` total `= dim Prym * 2`. So `A`'s
eigenvalue pairs `{alpha zeta^t, beta zeta^t}` would have to be the `H^1(E, rho^t)` eigenvalues of a
*ramified* character -- a concrete existence condition with no a-priori contradiction. (This corrects an
earlier mistaken intuition that `A`'s eigenvalues carried an "unramified signature".)

Still genuinely open, needing cohomology / polarization theory or a specialized search (the field
`F_2^131` is too big to brute-force): whether such a ramified character / genus-131 cover exists, and
whether `W` admits an *indecomposable* principal polarization (a Jacobian's theta divisor is irreducible;
`W` is only *isogenous*, not isomorphic, to `E x A`, so this is not settled by the isogeny
decomposition).

## Honest status

What is established by direct computation: `A` is simple, ordinary, carries the ECDLP prime, has CM
field `Q(zeta_131, sqrt(-7))`, and is **not ruled out** by any computable zeta-function obstruction
(point counts, place counts, effective-divisor counts all pass, for both `A` and `W = E x A`). The easy
escapes (splitting, a p-rank obstruction, a negative virtual count) are closed. Cyclic/abelian covers of
`P^1` are rigorously excluded (the `sqrt(-7)` obstruction).

**The "no automorphism" case turns out to be constrained.** The natural model of `A` (with its full CM
ring `End(A) = O_K`) carries a rational order-131 automorphism `rho = zeta_131` (the Weil-restriction
deck transformation), and `rho` preserves every polarization (Rosati = complex conjugation on `K`,
`zeta * conj(zeta) = 1`). The earlier `sqrt(-7)` proof that **`A` is not isomorphic to a Jacobian is
RETRACTED** (see CORRECTION above; the eigenvalue-field step was wrong). A candidate ray-class rescue of
the isomorphism case is noted there, pending review. The sketch below is kept only to show the structure; if
`(A, theta) ~= (Jac C, theta_C)`, Torelli forces an order-131 automorphism of the genus-130 curve `C`,
Riemann-Hurwitz forces `C -> P^1` cyclic degree-131 (the unique solution is `g_D = 0`, 4 branch points),
and then `Jac(C)`'s eigenvalues lie in `Q(zeta_131)` -- contradicting `A`'s `sqrt(-7)`. The "curve with
no useful automorphism group" cannot occur for the natural model: `A`'s own automorphism forces one onto
any isomorphic curve, collapsing the general case back onto the (excluded) `P^1`-cover case.

**What is genuinely left (isogeny residue).** The ECDLP transfer needs only `A` *isogenous* to a
Jacobian. A different model `A'` in the isogeny class, with `End(A')` a smaller order NOT containing
`zeta_131`, would have no forced automorphism and is not excluded by the argument above. So the open
residue is precisely: **does the minimal endomorphism ring `Z[pi, pi-bar]` contain `zeta_131`?** If yes,
every isogenous model carries the automorphism and `A` is not isogenous to a Jacobian (full closure); if
no, a `no-automorphism` model `A'` is not excluded. This is a concrete endomorphism-ring question, not
resolved here. (History: an earlier summary wrongly called the whole window "closed"; it was corrected
to "cover-of-E open"; class field theory then closed the cover-of-E route; and the automorphism/Torelli
argument now closes the isomorphism case, leaving only the isogeny residue above. This result should be
independently reviewed.)

## Cheap obstructions for the cover-of-E route: exhausted (none closes it)

Recording which obstructions were tried and FAIL, so they are not re-attempted:

- **Place / point / divisor counts for W:** all pass (above). Not an obstruction.
- **The sqrt(-7) field argument:** does not apply to covers of E (E has the CM). Not an obstruction.
- **A per-character c_1 trace argument** ("the degree-1 trace of L(E, rho^t) must equal zeta^t, a sum of
  roots of unity that cannot hold"): **INVALID.** It assumes an F_2-rational order-131 Hecke character
  rho^t, but `ord_131(2) = 130` means the 130 nontrivial order-131 characters form a single Frobenius
  orbit, so there is NO nontrivial F_2-rational order-131 character; the per-piece trace is only defined
  over F_{2^130} (a sum over ~2^130 points, with no cheap form). This attempt does not close the route.
- **E[131] as a Galois module:** `x^2 + x + 2` is irreducible mod 131 (disc `-7` is a non-residue), so
  E[131] is irreducible over F_2 and **E has no rational 131-isogeny**. This is suggestive *toward* the
  route staying open: the E<->A gluing inside W cannot split rationally, so W is plausibly
  *indecomposable* as a ppav (a necessary condition to be a Jacobian) -- but indecomposability is not
  thereby proven, and it is only necessary, not sufficient.

Net: the *cheap counting* obstructions do not close the cover-of-E route -- but **class field theory
does** (see "Cover of E / Prym route: CLOSED for cyclic covers" above). The ray class group `Cl^0_m` of
`E` for a genus-131 conductor has order `4` or `12`, coprime to 131, so no connected cyclic degree-131
cover of genus 131 exists; 131-torsion needs a degree-130 place in the conductor, forcing genus
`>= 8451`. The `E[131]`-irreducibility note above was *suggestive toward* openness; the ray class group
computation is the actual decision, and it closes the route. (The per-character `c_1` attempt remains
recorded as invalid.)

## Cover-of-P^1 route: small-ell test (corrected 2026-11-26)

By Torelli, if `A` (natural model) were isomorphic to a Jacobian it would be the Jacobian of a 4-branch-point
cyclic degree-`ell` cover of `P^1`. Over `F_q = F_{2^{ell-1}}`, `A_ell`'s eigenspace traces are `s*zeta^-j`
with `s` = E's trace over `F_q`. For a cover `y^ell = f(x)`, `#C(F_q) = q+1+sum_m (S_m + eps)` with
`S_m = sum_{f(x)!=0} chi^m(f(x))` and `eps = 1` for the unramified point(s) at infinity, so the trace is
`-(S_m+eps)`; `S_m` is the Galois conjugate of `S_1`, so a match needs `-(S_1+eps) = s*zeta^k`, i.e. the
`dlog`-histogram `h' = h + eps*e_0` constant except one entry offset by exactly `-s`
(`experiments/ecc2k130-weil-factor/cover_trace_test.py`).

| ell | q | A needs | covers (4 rational branch pts) | spread min/max | coherent-trace covers: trace integers | match |
|---|---|---|---|---|---|---|
| 5 | 16 | 1 | 169 | 2/5 | -3, 2 | 0 |
| 11 | 1024 | 57 | 92,911 | 4/64 | -64, -53, -20, 13, 24, 35, 46 | 0 |
| 13 | 4096 | -47 | 544,369 | 4/128 | -128, -102, -76, -37, -24, -11, 80 | 0 |

No cover with four `F_q`-rational branch points matches `A_ell` at `ell = 5, 11, 13`. Some covers do have
the coherent "integer times root of unity" trace structure, but with other integers (e.g. `-64`, `-128`
equal traces of supersingular curves over `F_2`), never `s`.

**All branch configurations (`cover_trace_nonrational.py`).** Conjugate branch points must carry equal
exponents, so a 4-point branch locus over `F_q` is: 4 rational; 2 rational + a conjugate pair; two
conjugate pairs; 1 rational + a cubic orbit (a single degree-4 orbit would need `4a = 0 mod ell`, impossible).
Each family was enumerated exhaustively up to PGL_2(F_q) (sharply 3-transitive on rational points,
transitive on conjugate pairs) and exponent scaling, with the point count `#C = q+1+sum_m(S_m+eps)`
validated against direct counting on random covers of every family, the degenerate `f = q1^ell` case
excluded, and every coherent hit checked against the Weil bound `2 sqrt(q)`:

| configuration | ell=5 | ell=11 | ell=13 | matches |
|---|---|---|---|---|
| 4 rational | 169 | 92,911 | 544,369 | 0 |
| 2 rational + pair | 24 | 4,608 | 22,528 | 0 |
| two conjugate pairs | 119 | 523,775 | 8,386,559 | 0 |
| rational + cubic orbit | 85 | 349,525 | 5,592,405 | 0 |

`ell = 19` (`F_2^18`), sampled: 20,000 random covers with four rational branch points, max spread 972
(Weil bound 1024), no coherent-trace cover, no match (`cover_trace_sample19.py`; a random sample,
not exhaustive).

**Correction.** The first version of this section (PR #70 draft) said covers' histograms have spread
`<= 4` and are "incompatible by a wide margin". That was wrong: 4 was the *minimum* spread; the maximum is
64 / 128, which contains the needed 57 / 47. The detector also looked only for a spike up and omitted the
point at infinity. The tables above are the corrected result. The conclusion (no match) stands, as an
exhaustive check over every branch configuration at `ell = 5, 11, 13`, not as a margin argument.
`ell = 5` is also confirmed independently by brute-force point counting over `F_16^k`. The trace is a
necessary condition, so no determinant check is needed: nothing passes the trace test.

Scope: so at `ell = 5, 11, 13` no 4-branch-point cyclic cover of `P^1` over `F_{2^{ell-1}}` has a Jacobian
with `A_ell`'s eigenspace traces -- by Torelli, `A_ell` (natural model, with its order-`ell` automorphism)
is not isomorphic to a Jacobian at these `ell`. Not covered: `ell = 131` itself (an extrapolation from
small `ell`), and the isogeny case (a model `A'` in the isogeny class without the order-`ell`
automorphism).

## Isogeny residue: zeta_131 is NOT in Z[pi, pi-bar] (2026-11-27)

The cover-trace evidence above applies to models of A that carry the order-131 automorphism zeta.
By Deligne's equivalence for ordinary abelian varieties, every order R with Z[pi, pi-bar] <= R <= O_K
is End(A') for some A' isogenous to A. So the question is whether zeta lies in Z[pi, pi-bar].

Set pi = alpha*zeta and pi-bar = beta*zeta^-1, with pi*pi-bar = 2. Then Z[pi, pi-bar] is the Z-span of
pi^0..pi^259 and pi-bar^1..pi-bar^259, written in the basis of O_K = Z[zeta_131] (x) Z[(1+sqrt-7)/2].

| ell | U_ell (Lucas, x^2+x+2) | [O_K : Z[pi,pi-bar]] | zeta in Z[pi,pi-bar]? |
|---|---|---|---|
| 5 | -1 | 1 | yes |
| 11 | 23 | 23^9 | **no** |
| 13 | -1 | 1 | yes |
| 19 | -457 | 457^17 | **no** |
| 131 | 263 * 146505763881528721 | (263 * 146505763881528721)^129 | **no** |

Two independent methods give these answers:
- `zeta_in_order.gp`: exact Hermite normal form over Z, with the index and an exact membership solve.
  The denominator of zeta's coordinates is |U_ell|.
- `zeta_in_order_modp.py`: pure Python, no PARI. Mod p = 263, the 519 generators span only 131 of 260
  dimensions, and adding zeta raises the rank to 132. So zeta is not in Z[pi, pi-bar] + 263 O_K, and hence
  not in Z[pi, pi-bar]. The same test reports the rank drop at ell = 11 (p = 23) and ell = 19 (p = 457).
  The controls (ell = 5 and 13, and p = 3 at ell = 131) give full rank, with no false positive.

**Structure (now PROVED, see `lead2_index/PROOF.md`; originally computed at four ell).** At these four ell, Z[pi, pi-bar] fails to be
maximal exactly at the primes dividing the Lucas number U_ell = (alpha^ell - beta^ell)/(alpha - beta),
and the index is |U_ell|^(ell-2). The mechanism: modulo a prime p dividing U_ell, alpha/beta is an
ell-th root of unity, so pi-bar collapses onto pi's span and the rank halves (131 of 260). So zeta lies
in the minimal order iff U_ell = +-1. Among primes ell <= 37, that happens only for ell = 2, 3, 5, 13.

**Consequence: the isogeny residue is a real gap.** Some varieties A' isogenous to A have
End(A') = Z[pi, pi-bar]. Such an A' contains no nontrivial power of zeta, so it has no order-131
automorphism. The Torelli reduction (a Jacobian in the isogeny class must be a cyclic degree-131 cover of
P^1) does not apply to these models, and neither does the cover-trace evidence.

The gap is localized:
- The models differ from the natural one only at p = 263 and p = 146505763881528721. They are reached
  from A by isogenies whose degree is a power of these primes.
- That degree is coprime to N = #A(F_2), so the ECDLP would transfer to any of them for free.
- Whether one of them is principally polarized and a Jacobian stays open. No construction is known.

This is a reversal of the hoped-for closure, and the record should say so: the cover-of-P^1 evidence
covers the natural model only, not the isogeny class.

**Index formula: proof (lead 2, 2026-11-27).** `experiments/ecc2k130-weil-factor/lead2_index/PROOF.md` proves, for every
odd prime ell, that [O_F (x) Z[zeta_ell] : Z[pi, pi-bar]] = |U_ell|^(ell-2), that the quotient is (Z/|U_ell|)^(ell-2),
that |U_ell| O_K is contained in Z[pi, pi-bar], and that zeta lies in the order iff U_ell = +-1. For ell != 7 the ring
O_F (x) Z[zeta_ell] is O_K.

How the proof works:
- Lift to O_F[T]/(T^ell - 1). There Z[alpha T, beta T^-1] is graded, and each graded piece has index |U_ell| by the
  pairing D(alpha^i, beta^j) = U_{i+j}.
- Divide by the norm line O_F * N, which meets the ring in f O_F * N (index f^2). This gives f^ell / f^2.

Every prime dividing U_ell is congruent to +-1 mod ell: 263 = 2*131 + 1, and 146505763881528721 = -1 mod 131.
The formula was checked by exact HNF and SNF for every prime ell <= 113, by an independent pure-Python HNF for ell <= 61,
and here independently by SNF at ell = 11, 17, 19, 23, 29. No Lean proof yet; the container has no Lean toolchain.

## Principal polarizations: every ppav in the class lacks zeta_131 (lead 1, 2026-11-27)

`experiments/ecc2k130-weil-factor/lead1_polarization/REPORT.md` derives a parity rule from Deligne modules, Howe's
polarization criterion and Hilbert reciprocity over K+. Let eps = (-7/ell), and let n(T) count the
self-glued inert places v0 (above primes dividing U_ell that are -1 mod ell and inert in Q(sqrt-7)) where
End(T) is non-maximal. A Deligne module T carries a principal polarization only if eps * (-1)^n(T) = +1,
and every conjugation-stable order with that parity occurs.

**At ell = 131.** eps = -1, 263 is of split type, and p2 = 146505763881528721 is of inert type. So:
- **A_131 itself (End = O_K) is NOT principally polarizable.** This is consistent with the natural
  restricted-product polarization, whose degree is 131^2 (docs/ecc2k130-explicit-ppav.md).
- **Every ppav in the isogeny class has End non-maximal at v0 above p2,** so zeta_131 is not in End.
- **Such ppavs exist:** the explicit module of PR #57, re-audited independently, has det 1, positive type,
  and is non-maximal exactly at v0.

**Consequence for the record.** No Jacobian isogenous to A_131 can have an F_2-rational automorphism of
order 131. So the Torelli -> cyclic-cover-of-P^1 reduction, and the cover-trace evidence above (PR #70),
apply to NO Jacobian in this class. That evidence remains correct as a statement about cyclic covers,
but it does not bear on the window question. The window question is exactly: is one of the ppavs
(End non-maximal at v0, no order-131 automorphism) a Jacobian? This is the "no useful automorphism"
case, with no known decision procedure at dimension 130.

**Controls.** The rule matches every computed case:
- ell = 3, 5, 13, 19: no ppav at all. At ell = 3 this is also confirmed unconditionally by enumerating
  all 768 genus-2 curves over F_2.
- ell = 11: O_K is principally polarizable.
- ell = 17, the same type as 131: a ppav exists iff End is non-maximal at v0 (45/45 vs 0/45 classes).

The small-ell enumerations assume GRH. The ell = 131 statement does not.

## Explicit models without zeta_131, and where the Jacobian candidates live (lead 3, 2026-11-27)

`experiments/ecc2k130-weil-factor/lead3_deligne/REPORT.md` constructs A' explicitly as a Deligne module, checked by two
independent methods (PARI and pure Python).

- **Minimal degree.** The smallest isogeny A -> A' with zeta not in End(A') has degree 263, and this is proved
  minimal. The smallest degree among those of p2-power order is p2.
- **End and Aut.** End(A') is a Gorenstein order of index p in O_K, and Aut(A') = {+-1}.
- **Kernel.** The isogeny kernel is cyclic, and Frobenius acts on it by lambda with lambda^2 = 2.

**Combined with lead 1.** 263 is of split type, so no index-263 model is principally polarizable. Every
possible Jacobian in the class (a ppav, End non-maximal at v0 above p2) is reached from A only through
isogenies of p2-power degree. On the minimal such isogeny, Frobenius acts on the kernel with order
(p2-1)/12, about 1.2e16. So the kernel points are defined only over F_{2^(1.2e16)}.

**Practical note: SUPERSEDED (2026-11-28, see 'Round 2' below). The kernel-field fact holds, but it is not an evaluation barrier.** Even if such a model were a
Jacobian, carrying the ECDLP across that isogeny would need an isogeny evaluation in dimension 130 with
kernel points over a field of degree about 1.2e16. No known method does this. The window question is
therefore open in principle (is one of these ppavs a Jacobian?), but any affirmative answer would also
need a way to cross this isogeny before it gave an attack.

## Round 2: all open leads in parallel (2026-11-27/28)

These are agent reports in `experiments/ecc2k130-weil-factor/round2/*/REPORT.md`. The independent skeptic pass failed on a usage
limit and is being re-run. Until it finishes, everything below is an agent result reviewed only by the parent session.

**Red-team of lead 1 (`redteam-lead1`).**
- **Upheld at ell=131.** Theorems A and B were re-derived and recomputed, unconditionally. No ppav in the class has
  zeta_131 in End, and A itself has no principal polarization. The independent cross-check: lambda_nat has valuation
  -128 (even) at the place above 131, while the PR #57 ppav has odd valuations exactly at that place and at v0.
- **Bug found.** A false-negative bug in `decide_ppav.gp` changed the small-ell counts but no verdict; corrected in
  `lead1_polarization/REPORT.md`.
- **Hypothesis.** The general rule needs v_p(U_ell) = 1.

**ell = 29 and 37 (`prediction-29-37`).** The parity rule is confirmed by unconditional det-1 certificates: ppavs with
End = O_K exist, and at ell=37 so does one on the order non-maximal at both self-glued places. Orders non-maximal at
exactly one v0 are excluded. The finer point: the specific module O_K (A_ell itself) is not principally polarizable at
29 or 37, for a narrow-class reason (the ramified prime above ell has no totally definite generator).

**Jacobian obstructions and expected count (`jacobian-obstructions`).** No known obstruction that varies within the
isogeny class rules out a Jacobian. Newly settled for the whole class:
- **Serre's twist obstruction is vacuous.** P_A(-x) would need #C(F_4) < #C(F_2).
- **No hyperelliptic curve.** a_1 = 2 and a_2 = 0 are incompatible.
- **Every ppav is geometrically indecomposable.**
- **A rational symmetric theta divisor always exists.** P_A = x^130 Phi_131 mod 2.

A curve with Jacobian in the class would be non-hyperelliptic, with Aut_{F_2}(C) = 1. The random-model heuristic gives
log2 E[#Jacobians] about -3697 (arithmetic band [-3700, -3567]). It rests on Lang-Weil leading terms that are not
justified at q=2, g=130, so read it as "no reason to expect one", not as a bound. The computed inputs include
h^-(K_131), about 2^209.46, from the analytic formula, cross-checked against the cyclotomic tables and bnfinit at
small ell.

**Isogeny evaluation (`isogeny-evaluation`). This corrects the previous section.**
- **(a) The kernel-field fact is a theorem, conditional on lead 1.** Every F_2-isogeny between A and a principally
  polarizable variety has kernel points defined only over F_{2^m} with (p2-1)/12 | m.
- **(b) It is NOT an evaluation barrier.** If some ppav were Jac(C) over F_2, there would be a nonconstant map
  f: C -> E over F_{2^131}. Then Phi = Tr o f^* is an F_2-isogeny A -> Jac(C) that can be evaluated in
  poly(deg f, g) group operations, with no kernel points.
- **The least deg f** is the minimum of an explicit Hermitian lattice. Blichfeldt gives <= 58 for the PR #57 module.
- **So the real barrier is the existence of C,** not kernel fields. The earlier "practical barrier" note is withdrawn.

**Cover of E (`cover-of-E-review`). This upgrades the old candidate status.** The old class-field-theory argument
(`rayclass_cover_E.py`) had a gap. It ignored twisted deck groups, and mu_131 Kummer covers of genus 131 do exist.
New results:
- **Twisted forms excluded by a trace lemma.** Frobenius permutes the eigenspaces with no fixed point, so the Prym
  trace is 0, but tr(pi|A) = 1.
- **Every geometrically Galois cover is closed** (any group, any twisted form, genus 131 with Jac(C) ~ E x A'). A
  gluing-degree constraint (Kani) forces deg(C -> E) to divide 131 * 263 * p2. Together with Riemann-Hurwitz, this
  leaves only 131, which the trace lemma and class field theory close.
- **Still open:**
  - degree-131 maps with AGL-type monodromy (an explicit consistent configuration exists, pending 131 | h(E') for a
    genus-11 E');
  - A_131 / S_131 monodromy and the other divisors n.
- **Lead 1 does not exclude this route.** The Prym is not principally polarized; it has type (1^128, n, n).

## Round 2: skeptic pass and the remaining leads (2026-11-28)

An independent skeptic re-ran the decisive computation of each lead with its own code. The core of all five leads is
upheld. Scope and label corrections:
- **jacobian-obstructions.** "Aut = {+-1} for every ppav" is conditional on lead 1. Note that lead 1 has since been
  upheld by the red-team, but only as an agent derivation. The ppav count of 2^4431 is heuristic and conditional on
  Theorem B.
- **isogeny-evaluation.**
  - "Any useful transfer must be an isogeny" holds only for algebraic maps into abelian varieties; other transfers are
    excluded only heuristically.
  - The ell = 11 and 17 sampling, and the Hermitian d_min formula, could not be reproduced (a script load failure), so
    they are unsupported.
  - The kernel-field lemma (a) stands, modulo lead 1.
- **redteam-lead1.** The skeptic's from-scratch End = O_K decider and its place structure above p2 (one self-glued v0,
  64 glued pairs) agree. Not re-checked: the ell = 29 lattice (timed out) and the Theorem B existence direction
  (rests on the PR #57 certificate).
- **prediction-29-37.** The narrow-class obstruction for T = O_K was re-proved by an independent route. "Every order
  non-maximal at exactly one v0 is excluded" is checked only for orders maximal away from v0.
- **cover-of-E-review.** The trace lemma, the gluing-degree constraint and the Galois exclusion are upheld; the content
  of P_A(alpha) = 131 * U_131 was recomputed two ways. The AGL arithmetic-equals-geometric monodromy step is
  unsupported (that case is listed as open anyway).

**Loose ends (`round2/loose-ends`, skeptic-upheld).**
- p2 and 263 are proved prime by independent pure-Python Lehmer (N-1) and Morrison (N+1) certificates, with controls.
- Pic of the index-p orders follows a closed formula, h_K (p -+ 1)/(ell Q), which matches ell = 11, 17, 19 (GRH).
- At ell = 131, the number of isomorphism classes with End = the order at v0 is h_K (p2+1)/262, about
  6.3e77 * h^+(K+). Here h^-(K_131), about 1.13e63, is unconditional, and h^+ is unknown.
- **Lean: blocked.** Running the elan installer was denied by the session's permission classifier, and the agent did
  not work around it. Lemmas 1 and 2 of the index proof were checked symbolically only.

**Extension fields and A as a factor (`round2/extension-and-factor-routes`).**
- **For 131 not dividing k: proved.** Q(pi^k) = K, so End and Aut do not grow over F_{2^k}, and no order-131
  automorphism appears for any Jacobian in the class.
- **At 131 | k: heuristic dead end.** A ~ E^130 there, but transferring into a genus-130 Jacobian over F_{2^131}
  (order about 2^17030) costs about 2^680 by large-genus index calculus.
- **A as a factor of a larger Jacobian.** GHS covers of P^1 and all geometrically Galois covers of E are closed.
  Non-Galois maps to E, and curves with no map to E, remain open.
