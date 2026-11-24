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

## Honest status

This is a **necessary, not sufficient**, result. Passing the point- and place-count tests does not
construct a curve; for dimension `g ≥ 4` there is no known decision procedure for "is this simple
ordinary abelian variety a Jacobian", so the window question stays open. What is now established, by
direct computation: `A` is simple, ordinary, carries the ECDLP prime, has CM field `Q(zeta_131,
sqrt(-7))`, and is **not ruled out** by any computable zeta-function obstruction. The easy escapes
(splitting, a p-rank obstruction, a negative virtual point/place count) are all closed off. The
remaining route to a verdict is a CM-type / Jacobi-sum comparison against explicit cyclic-cover
families (`PARI/GP` is now available in-container for this), reconciled with the crypto repo's
genus ≥ 394 bound.
