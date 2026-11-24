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

## Sub-result: cyclic covers OF P^1 cannot realize A in-window

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

## Cover of E / Prym route (open, not closed)

The natural carrier for `A` is a cyclic cover of `E` itself, not of `P^1`. A degree-131 cover
`C -> E` totally ramified at 2 points has genus 131 (`2g-2 = 131(2*1-2) + 2*130 = 260`), in-window,
with `Jac(C) ~ E x Prym`, `dim Prym = 130`. For `Prym ~ A` we need `Jac(C) ~ E x A = W`, where
`W = Res_{F_2^131/F_2}(E)` has `P_W = (x^2+x+2) P_A = x^262 - s_131 x^131 + 2^131`. So:

> **A is the Prym of a degree-131 cover of E  ==>  W = Res(E) is a genus-131 Jacobian (carrying a
> degree-131 map to E).**

The `P^1` field obstruction does **not** apply here: `E` already has CM by `sqrt(-7)`, so a cover of `E`
can supply it. Testing the place-count necessary conditions for `W`
(`experiments/ecc2k130-weil-factor/cover_of_E_reduction.py`):

- All place counts `a_d >= 0` (`d = 1..399`).
- Strikingly, **`#C(F_2^k) = 2^k + 1` for all `k = 1..130`** -- the virtual curve looks exactly like
  `P^1` up to degree 130, revealing its genus only at degree 131.

So **`W` passes the cheap necessary conditions and the cover-of-E route is open** -- it is not ruled out.
Whether `W` is actually a genus-131 Jacobian (hence whether the cover exists) is the real open question,
and it is attack-relevant: a genus-131 curve `C` with `Jac(C) ~ W` would transfer the ECC2K-130 DLP to
genus-131 index calculus (cost on the order of the model in the previous section). This is the correct
state of the window question; the earlier "every attackable path is closed" summary was wrong.

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

**What remains open.** Covers of `E` itself (the Prym route, reducing to "is `W = Res(E)` a genus-131
Jacobian") are NOT closed and are attack-relevant; and the fully general case (a Jacobian with no useful
automorphism group) has no known decision procedure at dimension 130. These are the live questions. An
earlier summary that called the window "closed / solved in the negative" was wrong and is retracted
here: the Prym-of-E route in particular is open.

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

Net: the cheap obstructions are exhausted and none closes the cover-of-E route. A genuine verdict needs
explicit cohomology over F_{2^130} (the ramified order-131 twisted L-functions of E) or a
polarization-lattice computation of W's principal polarization -- real tools beyond cheap counts. The
route is open.
