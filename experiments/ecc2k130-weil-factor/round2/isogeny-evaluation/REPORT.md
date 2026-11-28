<!-- Agent report saved by the parent session from the workflow journal. Status: mixed. Independent skeptic verification: NOT YET RUN (failed on usage limit). -->

# Round 2, isogeny-evaluation: is the transfer barrier a theorem?

**Verdict.** There are two separate statements here.

1. **Kernel-field barrier: a theorem**, conditional on Lead 1 Theorem A at ell=131. Every F_2-isogeny, in either direction, between A and a principally polarizable variety has a kernel part at the self-glued ideal m0 above p2. Every nonzero point of that part needs F_{2^m} with **(p2-1)/12 | m**, and (p2-1)/12 = 12208813656794060.
2. **Evaluation barrier: not a theorem, and refuted as one.** Kernel points over a huge field do not prevent cheap evaluation. If any ppav in the class is a Jacobian Jac(C), then an F_2-isogeny A -> Jac(C) exists that is evaluable in poly(deg f, g) through a map f: C -> E over F_{2^131}, without any kernel point. Its kernel still meets m0.

So the real obstruction is whether a curve C exists and can be found together with f. Kernel fields are not the obstruction.

## (a) The kernel theorem

**Setup.** R = Z[pi, pi-bar]. Deligne modules are covariant, and for an isogeny x: T_C -> T_B the prime-to-2 kernel is the finite R-module G = T_B / x T_C, with pi acting as Frobenius. A finite R-module splits as the sum of G_m over the maximal ideals m of R.

**Residue data at p2 (computed).** Above p2, R has 66 maximal ideals: one self-glued m0, 64 glued pairs and 1 unglued place. R/m0 = F_p2, and pi = pi-bar = lambda there, with lambda^2 = 2. The value is lambda = 14975904886999252, and ord(lambda) = ord(-lambda) = (p2-1)/12. Only the places of K over v0 contract to m0.

**Lemma (proved).** Let C have End(C) maximal at v0 (e.g. A with End = O_K, or A^dual), and let B be principally polarizable. Then for any isogeny C->B or B->C, G_{m0} != 0.

*Proof.*
- If G_{m0} = 0, the two lattices agree after completing at m0.
- (T_C)_{m0} is a module over (O_K)_{m0} = O_{v0}, so End(B) would be maximal at v0.
- That contradicts Lead 1 Theorem A: eps = (-7/131) = -1 forces n(T) odd, so End(B) is non-maximal at v0.

**Consequence (proved).**
- Let P be a nonzero point of G_{m0}. Then M = R·P satisfies M/m0M != 0 (Nakayama), and pi acts on M/m0M as the scalar lambda. So pi^m P = P forces lambda^m = 1.
- A point of G with a nonzero m0-component carries that component in Z[pi]·P, so it needs the same field.

**What the lemma covers.** Compositions (lengths add, so some step has an m0-part), higher p2-power kernels (unipotent extensions only multiply m by powers of p2), duals (the Cartier dual carries pi-bar, which is congruent to lambda mod m0) and isogenies *to* A.

**A kernel cannot be defined over a smaller field.** The subgroup scheme is F_2-rational, but its m0-points are not rational over any F_{2^m} with ord(lambda) not dividing m.

**Correction to Lead 3.** The degree must be divisible by p2, but need not be a p2-power. The PR #57 module is reached by an isogeny of degree 263·p2.

### Checks

- **Method A** (`residues.gp`, PARI finite fields) and **Method B** (`residues_check.py`, hand-built F_{p^2}, Miller-Rabin) agree on lambda and its order for all four (ell, p). They also match Lead 3's table:

| ell | p | lambda | ord(lambda) |
|---|---|---|---|
| 11 | 23 | 18 | 11 |
| 17 | 271 | 96 | 270 |
| 19 | 457 | 377 | 152 |
| 131 | 263 | 110 | 262 |
| 131 | p2 | 14975904886999252 | (p2-1)/12 |

- **ell=131 instance** (`ppav131_kernel.py`, own arithmetic in Z[zeta][alpha], reading `explicit_ppav.json`):
  - The kernel I/O_K = Z/(263·p2).
  - pi acts by c with c = 110 mod 263 and c = lambda mod p2. This equals the JSON's `F_and_V_on_cyclic_quotient`: three-way agreement.
  - Kernel points need 1599354589040021860 | m. The 263-part alone needs only F_{2^262} (negative control).
- **ell=17 positive test** (`isog_small.gp`, `run17a.out`, `run17b.out`; GRH class group). Sources are two O_K-classes; targets are 4 ppav classes with End = R + v0.
  - Into ppav targets: 571 of 571 isogenies have m0-socle > 0, and 270 | m_G for all.
  - **Control**, O_K targets: 0 of 249 have m0-socle > 0.
  - **Control**, R-stable lattice non-maximal only at a glued pair (not a ppav): 0 of 140.
  - 140 of 140 glued-pair-control isogenies have 270 | m_G with no socle, so plain order divisibility is a weak detector. The socle is the right one.
- **ell=11 negative control** (eps=+1, `run11.out`). ppavs with End = O_K exist, and most isogenies into them (e.g. 46 of 49) have no lambda-socle. The barrier is absent, as Lead 1 predicts.
- **ell=19.** No ppav exists, so the theorem is vacuous. The lambda order of 152 matches Lead 3.

## (b) Routes that avoid kernel points

**1. Homomorphism given only on the N-subgroup: no (proved, trivial).**
- A homomorphism Z/N -> Z/N known only by the image of a generator is DLP-equivalent.
- Interpolating an isogeny on A[N] needs A[N]. The Weil-dual line to A(F_2) has eigenvalue 2, so this needs F_{2^k} with ord_N(2) | k, and ord_N(2) = (N-1)/24 (computed, `ntorsion.out`).
- A useful transfer must therefore be an algebraic F_2-map. By rigidity it is a homomorphism, and since the target is simple it is 0 or an isogeny, so (a) applies to it.

**2. Endomorphisms first, then a quotient: no escape.**
- x in O_K acting A -> A, followed by A -> A/H, is the isogeny with kernel x^{-1}(H). (a) applies.
- Endomorphisms themselves are cheap even when their kernel meets m0: [p2] is the trivial example. So "kernel in a huge field" never implied "hard to evaluate".

**3. Through E over F_{2^131}: the decisive route (proved).**
- Suppose B = Jac(C), with C/F_2 of genus 130.
- J_{F_{2^131}} ~ E^130, so Hom(J, E) != 0. u composed with Abel-Jacobi gives a nonconstant f: C -> E over F_{2^131}; C(F_{2^131}) is non-empty by Weil.
- Define Phi(P) = sum_{i=0}^{130} Frob^i(f^*(P)) for P in E(F_{2^131}) = Res E(F_2). Phi is Galois-equivariant, so it is defined over F_2.
- Hom_{F_2}(E, J) = 0 and Phi != 0, because Phi^dual Phi is circulant with diagonal deg f > 0. So Phi restricted to A is an F_2-isogeny A -> J.
- Its kernel meets m0 by (a), yet evaluation costs one pullback of a degree-d divisor plus 131 Frobenius-additions in J: poly(d, g, 131).
- Caveat: Phi may kill A(F_2) only if N divides its degree. This is avoidable by changing u, which is argued here but not proved.

**Cover degree (proved formula, computed values).** The least d is the minimum of the lattice L = sqrt(-7)·(Z[alpha]T)^dual with form Q(c) = |Tr_{K/Q(sqrt-7)}(c c-bar / (sqrt(-7) lambda))|.
- **Formula:** det G = 7^(ell-1) · [Z[alpha]T : T]^2.
- **ell=11:** det G = 7^10, Q integral, d_min = 3, 4, 4.
- **ell=17:** I = 271, det G = 7^16·271^2, d_min = 5 for 4 ppav classes.
- **ell=131, PR #57 module:** I = p2, so det(Q)^(1/260) = 1.793. The Blichfeldt bound gives **d_min <= 58**; a Minkowski-Hlawka typical value is about 27 (heuristic).

So if that ppav were a Jacobian, a degree-at-most-58 map C -> E would carry the ECDLP across cheaply. A byproduct: d_min = 1 would split B over F_{2^131} and prove "not a Jacobian". That never occurred at ell=11 or 17.

**4. Generic high-dimensional algorithms (literature, cost heuristic, not verified here).**
- Theta-based methods (Lubicz-Robert, Cosset-Robert) need n^g >= 2^130 coordinates per point, plus kernel points.
- HD/Kani embeddings (Robert 2022+; the SIDH attacks) work in dimension 2g-8g, at 2^260 or more per step, and need torsion images of phi. The 2-torsion is étale rank only g, because A is ordinary.
- Analytic/canonical-lift methods need 2^130 theta constants.
- All of these are far above rho (about 2^65).
- GHS descent is degenerate for E defined over F_2 (Koblitz curves are GHS-immune). This is a literature claim.

## Status summary

| Status | Claims |
|---|---|
| Proved | Lemma and field bound (conditional on Lead 1 Thm A at 131); the cover-transfer proposition; det formula; N-subgroup triviality |
| Computed | lambda and orders (3 ways); the 131 PR #57 kernel; ell=11/17 lattice tests with controls; d_min at 11 and 17; ord_N(2) |
| Heuristic | typical d_min of about 27 at 131; generic-algorithm costs |
| Open | whether any ppav in the class is a Jacobian; constructing C and f; Lead 1 review |

**Precise reduction.** An efficient transfer to a Jacobian exists if and only if (up to the N-divisibility caveat) there is a C/F_2 of genus 130 with Jac C ~ A together with an explicit F_{2^131}-map C -> E. Lead 1's Theorem B ppavs exist with d_min <= 58. Whether one is a Jacobian is open.

## Caveats (agent)

- (a) depends on Lead 1 Theorem A at ell=131, an agent proof that has not been independently reviewed. The ell=17 analogue was checked by enumeration under GRH plus the sampling here.
- The ell=11 and ell=17 lattice tests sample isogenies (LLL basis plus random small combinations); they do not enumerate them.
- `run17a.out` came from an earlier script version. Its ppav and O_K sections are valid; its glued-pair control crashed because the old construction was not R-stable, and was redone in `run17b.out`.
- The ell=131 d_min value is only a Blichfeldt upper bound; no LLL was run in dimension 260.
- The cover-transfer proposition assumes a Jacobian exists. "Changing u avoids N | deg Phi" is argued, not proved.
- Algorithm-cost statements (theta, HD/Kani, GHS degeneracy) come from the literature and were not verified here.
- Primality of p2 and N rests on PARI isprime and Miller-Rabin.
- No git commits were made. No files were written outside round2/isogeny-evaluation/.
