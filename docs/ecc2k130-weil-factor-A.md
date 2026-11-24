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

The actual trace-zero representative `A` is **not F2-isomorphic to a genus-130 Jacobian**.
Its rational integral order-131 action preserves every F2-polarization (Rosati is complex
conjugation on its CM endomorphism field). Torelli would force a rational order-131 action
on the curve. Riemann–Hurwitz then requires four branch points, while Frobenius orbits
of branch points must have length divisible by `ord_131(2)=130`. This is impossible.
The earlier claim that Riemann–Hurwitz alone made this cover consistent over F2
omitted this descent constraint.

The isogeny question is different. The isogeny class **contains a principally polarized
variety**, by Howe's Theorem 1.1 with an explicit inert prime above
`146505763881528721` dividing `F-V`. The order-131 symmetry is not forced integrally
by `Z[F,V]`; a reduction modulo 263 proves this directly. Thus the symmetry argument
does not exclude a Jacobian elsewhere in the class. Proofs and reproducible arithmetic
certificates are in [the integral-symmetry follow-up](ecc2k130-integral-symmetry.md).

## Honest status

The exact representative `A` is excluded as a genus-130 Jacobian. Its isogeny class is
principally polarizable; whether that class contains a genus-130 Jacobian, or occurs
as a factor in a usable larger-genus Jacobian, remains unresolved. Checking only
cyclic-cover CM types cannot decide members whose integral cyclotomic action is absent.
