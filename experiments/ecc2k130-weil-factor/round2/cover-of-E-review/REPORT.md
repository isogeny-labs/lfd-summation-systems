<!-- Agent report saved by the parent session from the workflow journal. Status: mixed. Independent skeptic verification: NOT YET RUN (failed on usage limit). -->

# Re-review of the cover-of-E results (slug: cover-of-E-review)

**Verdict: closed except for one case (X).** At genus 131, every geometrically Galois cover C -> E with Jac(C) ~ E x A' is excluded (proved), whatever the group and whatever the twisted form of the deck group. Geometrically cyclic degree-131 covers whose Prym contains A are excluded at every genus up to 300 (proved).

**X (open):** maps C -> E that are not geometrically Galois. Their degree must be a divisor of 131 * 263 * p2 greater than 1. At degree 131 the monodromy must be one of:

- Z/131 x| U_d with d > 1, where arithmetic monodromy equals geometric monodromy (forced);
- A_131 or S_131.

The old statement ("A's zeta-eigenvalue structure requires the cyclic one", with degree 131 forced) was not justified. Its CFT computation is correct, but it missed whole classes of covers (below).

Files: everything is under `experiments/ecc2k130-weil-factor/round2/cover-of-E-review/`. The scripts are `content_gp.gp`, `content_py.py`, `kummer_control.py` and `constant_count_control.py`, each with a matching `.out` file.

## 1. What the old argument got right and what it missed

**Correct (re-derived).** For a cover with constant deck group Z/131 over F_2, the following all hold:

- The cover is geometrically connected exactly when its character is nontrivial on Cl^0_m.
- Because E has degree-1 places, Cl_m = Z (+) Cl^0_m.
- From the exact sequence 1 -> F_2^* (trivial) -> prod (O_v/v)^* -> Cl^0_m -> Pic^0 -> 0, we get |Cl^0_m| = 4 * prod(2^{deg v} - 1).
- Ramification is tame (131 is odd), so the conductor exponent is 1.
- Riemann-Hurwitz gives 2g - 2 = 130 * deg(m).

At genus 131 the conductor has degree 2, so |Cl^0_m| is 4 or 12, and neither is divisible by 131. 131-torsion first needs a place of degree 130, which means genus 8451.

**Missed, part 1: twisted forms.** The deck group can be a twisted form of Z/131 (Frobenius acting by u != 1). The mu_131 Kummer covers y^131 = f exist at genus 131: since 131 does not divide #E(F_2) = 4, multiplication by 131 is a bijection on Pic^0(E)(F_2), so there is an f with div f = a(P - Q) + 131D. The constant-group CFT says nothing about these covers.

**Missed, part 2: maps that are not Galois.** "Galois iff cyclic for prime degree" is true, but a degree-131 map need not be Galois. Nothing forces the degree to be 131 either; see section 2.

## 2. New rigorous constraints (proved)

**(a) Point count.** If Jac(C) ~ E x A', then:

- #C(F_2) = 3 - Tr(E) - Tr(A) = 3 + 1 - 1 = 3, where Tr(pi|A) = 1 was computed.
- #C(F_{2^k}) = 2^k + 1 for k <= 130.

**(b) Gluing degree.** Any map C -> E factors through an optimal map C -> E' followed by an isogeny E' -> E. Here E' is E itself: End = O_F, class number 1, so E is alone in its isogeny class.

- For an optimal map of degree n, Kani's result gives E cap P = f^*(E[n]), isomorphic to E[n]. So E[n] embeds in P, and P_A(pi) kills it.
- Hence n divides P_A(alpha) in Z[alpha].
- P_A(alpha) = 131 * alpha^130 * U_131. Here alpha^130 contributes no integer factor, and 131 does not divide U_131 (U_131 = -1 mod 131).
- So the content is 131 * |U_131| = 131 * 263 * p2. This was computed by PARI, and independently by Python via P_A(alpha) = P_W'(alpha)/sqrt(-7). The two agree at ell = 3, 5, 11, 13, 17, 19, 23, 29, 31, 37 and 131, and the content equals ell * |U_ell| every time.
- Conclusion: n is a divisor of 131 * 263 * p2 greater than 1, and n is odd. These are the same primes as lead 2's index and lead 3's 263 and p2 isogenies.

**(c) Twisted cyclic covers.** Let g be the geometric deck generator and u the twist, with F g = g^u F. Then pi maps the eigenspace V_k onto V_{uk}, and for u != 1 no V_k (k != 0) is fixed. So Tr(pi|Prym) = 0, which is not Tr(pi|A) = 1.

A second, independent proof: over each unramified rational point, Frobenius acts on the fibre torsor as x -> ux + b, which has exactly one fixed point. A ramified rational point has exactly one point above it. So #C(F_2) = #E(F_2) = 4, which is not 3. The full eigenvalue statement is that P_Prym lies in Z[x^m], with m = ord(u). P_A is not in any Z[x^m] (checked; its x^259 coefficient is -1).

**(d) Constant cyclic covers.** Besides the CFT argument: a rational point has 0 or 131 preimages if unramified and 1 if ramified, so #C(F_2) = 131a + r with r <= 2. That cannot equal 3. This gives two independent proofs.

**(e) All geometrically Galois covers at genus 131.** Take the optimal factor, which is Galois of degree n. It is tame because n is odd, so each branch point contributes at least 2n/3. Then 260 >= 2n/3 gives n <= 390. Combined with (b), n is 131 or 263.

- n = 263: each branch point contributes 262, which does not divide 260. Excluded.
- n = 131: cyclic, so (c) and (d) apply.

The whole Galois family is closed.

**(f) Window up to genus 300, cyclic degree 131 with A inside the Prym.**

- Eigenvalues closed under mu_m, and A's 260 eigenvalues stay pairwise distinct under those shifts. So dim Prym >= 130m.
- With g = 65r + 1 <= 300, r <= 4, so m <= 2.
- m = 1: CFT over F_2 would need 131 to divide some 2^d - 1 with d <= 4. It does not.
- m = 2: CFT over F_4 would need 131 to divide #E(F_4) = 8 or some 4^d - 1 with d <= 4 (values 3, 15, 63, 255). It does not.
- More generally, 131 does not divide #E(F_{2^m}) for any m | 130 (ord(alpha mod 131) = 5720 = 130 * 44), so an unramified 131-part never helps.

**(g) AGL type (degree 131, d > 1).** The Hecke field Q(zeta)^{U_d} acts on the Prym, and Frobenius acts on it by sigma_u. If u is not in U_d, the trace is 0, a contradiction. So G_arith = G_geom, and D -> E' is a constant Z/131 cover over F_2. CFT then needs 131 to divide h(E') or d = 130. Riemann-Hurwitz allows, for example, d = 10:

- C -> E has two branch points with contributions 104 (tame, e = 5) and 156 (wild Z/10, jump 3);
- E' -> E has 2g' - 2 = 8 + 12 = 20 = 2d, so g' = 11;
- D -> E' is etale, and the remaining condition is 131 | h(E').

Not excluded. **Open.**

## 3. Reconciling with lead 1

Lead 1 is about principal polarizations on varieties in A's isogeny class (dimension 130). In this route the Prym of an optimal degree-n map carries the restricted polarization, of type (1^128, n, n) with kernel isomorphic to E[n]. That is never principal, and the ppav here is Jac(C), in W's isogeny class (dimension 131).

The natural A, with End = O_K containing zeta, carries exactly such a polarization, with kernel E[131], from W = Res(E). Its natural principal polarization is geometrically a product, so W itself is not a Jacobian. So lead 1 neither excludes the route nor is used anywhere above. The exclusions come from Frobenius traces, point counts, CFT and Riemann-Hurwitz.

## 4. Controls

- **Twisted forms (positive control for the lemma; negative control for "Prym = A").** Covers C_ell: z^ell = x over E (genus ell) over F_2 with ell = 3, 5, 7, 11:
  - the Prym polynomial lies in Z[x^m] with m = ord_ell(2);
  - its trace is 0 and #C(F_2) = 4 = #E(F_2);
  - the Prym is never A_ell;
  - an independent count at degree g+1 matches the zeta function every time.
  - Over F_{2^13} at ell = 131: #C = #E = 8012.
- **Constant group (m = 1).** Bases F_4, F_8, F_16, F_64 and F_4096:
  - #C = 2 mod ell, as the count rule predicts;
  - the Prym traces are nonzero (-16, -12, -1344, -72), so the trace-0 signature is specific to twists. (One case, ell = 3 over F_4, happens to have trace 0.)
- **A bug caught by the controls:** the first version of my g+1 sanity check had a Newton sign error. It produced mismatches, I fixed it, and all checks pass now. Separately, an exponent bug for F_2 caused a hang and was fixed.

## 5. Status by claim

| Status | Claims |
|---|---|
| Proved | (a) to (f); the reduction in (g) to G_arith = G_geom; the lead-1 reconciliation |
| Computed | content = 131 * 263 * p2; Tr(A) = 1; ord(alpha mod 131) = 5720; all controls |
| Open | non-Galois maps of degree n dividing 131 * 263 * p2, n > 1, which are the AGL types Z/131 x U_d (d > 1) and A_131/S_131 at n = 131, and the maps of other degrees; the general window question (Jac(C) ~ E^a x A x B with g <= 300) beyond cyclic degree-131 covers |
| Heuristic only | the count of Jacobians against Weil classes makes existence implausible, but W's class is special |

**Suggested doc change (not made, as instructed):**

1. Replace "closed for cyclic covers" with the statement proved in (e) and (f).
2. Remove "A's zeta-eigenvalue structure requires the cyclic one" and "A as the Prym forces degree 131".
3. Add the gluing-degree constraint n | 131 * 263 * p2.
4. Record the AGL / S_131 residue as open.

The stakes are as the task says: a curve of this kind would give a genus-131 DLP target. The open residue is a search over non-Galois covers, for example looking for E' of genus d+1 with 131 | h(E') at d = 10 and checking the eigenvalue match.

## Caveats (agent)

No other agent has reviewed these proofs. Kani's lemma (E cap P is isomorphic to E[n] for optimal maps) is cited from the literature, not re-proved. The AGL step (g) relies on the Hecke-algebra/Frobenius twisting argument, which is only sketched. Maps that are not geometrically Galois remain fully open, and that includes A_131/S_131 monodromy, where no structure is available to exploit. For the window beyond genus 131, only geometrically cyclic degree-131 covers are closed; Galois covers of other degrees, and non-Galois covers, with A as a proper factor of the Prym are not addressed. The controls run only at small ell (3 to 13) plus one ell = 131 count over F_{2^13}. Nothing was committed, and no doc was edited.
