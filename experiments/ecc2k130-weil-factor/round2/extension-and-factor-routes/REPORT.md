# Item 4: routes outside "Jacobian isogenous to A over F_2"

The workflow agent for this lead was stopped by a safety classifier before producing anything. The parent session did
sub-route (a) directly and summarized (b) from the repository's existing results.

## (a) Automorphisms over extension fields
**Lemma (proved).** For 131 not dividing k:
- Q(pi^k) = K, because the 260 conjugates alpha^k zeta^{jk}, beta^k zeta^{jk} are distinct (see `ext_fields.py`; the
  numeric control shows distinctness fails exactly at k = 131, 262).
- Hence, by Deligne's equivalence with base change (T, pi^k) (standard), End_{F_{2^k}}(A') = End(A') for every A' in
  the class.
- No new endomorphisms or automorphisms appear over F_{2^k}. In particular every ppav keeps Aut = {+-1}, so no
  Jacobian in the class acquires an order-131 automorphism over any such extension.

**131 | k.** Here A x F_{2^131} ~ E^130 and End^0 = M_130(Q(sqrt-7)), which leaves plenty of room for automorphisms.
Using this means transferring A(F_2) (order N) into the Jacobian of a genus-130 curve over F_{2^131}. That group has
order about 2^17030.
- **Heuristic cost.** Large-genus index calculus there (Enge-Gaudry type, L_{q^g}[1/2, sqrt 2]) costs about 2^680
  operations, far above rho (about 2^65). Moreover g/log2(q) is about 1, the boundary of that algorithm's regime.
- **Verdict.** No gain; this route is only a heuristic dead end.

## (b) A as a proper factor of a larger Jacobian (g ~ 131-300)
| family | status | source |
|---|---|---|
| GHS Weil-descent covers of P^1 | closed: genus forced to about 2^130 (or zero map) | docs/ecc2k130-transfer-survey.md |
| cyclic covers of P^1 carrying A (natural model) | cover-trace evidence (PR #70); moot for Jacobians in A's class, still relevant for factors | docs/ecc2k130-weil-factor-A.md |
| geometrically Galois covers C -> E with Jac(C) ~ E x A' (g = 131) | closed (round 2, skeptic-upheld; the AGL monodromy step is only sketched) | round2/cover-of-E-review |
| geometrically cyclic degree-131 covers of E, A in the Prym, g <= 300 | closed (round 2) | round2/cover-of-E-review |
| non-Galois maps C -> E (AGL with d > 1, A_131/S_131 monodromy, other n dividing 131*263*p2) | OPEN | round2/cover-of-E-review |
| A a factor of Jac(C) with C not mapping to E (no structure) | OPEN; no decision procedure | — |

Nothing here was computed beyond what the cited reports computed.
