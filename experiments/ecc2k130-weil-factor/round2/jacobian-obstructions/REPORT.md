<!-- Agent report saved by the parent session from the workflow journal. Status: mixed. Independent skeptic verification: NOT YET RUN (failed on usage limit). -->

# Window question: is any ppav in the class a Jacobian? (slug: jacobian-obstructions)

Everything below lives in `experiments/ecc2k130-weil-factor/round2/jacobian-obstructions/`. Nothing was committed and no file outside that folder was changed.

## Bottom line
- **No known non-isogeny-invariant necessary condition decides the question.** Each condition I could find either holds for every ppav in the class, or is ruled out for the whole class, or cannot be tested without curve equations.
- **Four conditions are now settled for the whole class:**
  - Serre's twist obstruction cannot block anything: the twist class contains no curve.
  - The curve cannot be hyperelliptic.
  - Every ppav in the class is geometrically indecomposable.
  - The obstruction to a rational symmetric theta divisor vanishes.
- **What is left:** a non-hyperelliptic genus-130 curve C over F_2 with Aut_{F_2}(C)=1 and Jac(C) isomorphic to a ppav in the class.
- **Heuristic count:** the random-model estimate is log2 E[#Jac] ≈ **-3697** (band [-3700, -3567] on the arithmetic inputs). It carries a large modelling caveat, set out in section 2.
- **Status:** open, with the precise reduction given in section 1.

## 1. Survey of non-isogeny-invariant necessary conditions

| Condition | What it says | Testable here? | Result |
|---|---|---|---|
| N1. A principal polarization exists | — | yes (lead 1) | Holds only for orders non-maximal at v0 over p2. Not re-proved here. |
| N2. Geometric indecomposability of (A',λ) (theta divisor irreducible) | — | **yes, for the whole class** | **Holds for every ppav (proved).** If (A',λ) split over F_2-bar, Galois would permute the factors transitively, because A is simple over F_2. Uniqueness of the decomposition into indecomposables (cited) then gives A' ≅ Res_{F_{2^k}/F_2} B with k\|130, k>1, so P_A(x) = P_B(x^k). But the x^259 coefficient of P_A is -1 and gcd(259,130)=1. The detector flags a toy Q(x^2), so it can fire. |
| N3. Serre's obstruction | A ppav that is a Jacobian over F_2-bar (non-hyperelliptic curve) is, over F_2, either Jac(C) or the -1 twist of Jac(C). The twist has Frobenius -π and Weil polynomial P_A(-x). (Serre, appendix to Lauter 2001; Beauville–Ritzenthaler. Cited, not re-derived.) | **yes, for the whole class** | **The twist class contains no curve (proved).** It would need N_1=4 and N_2=2, but C(F_2) ⊂ C(F_4); a_2=-1 as well. So any ppav in the class that is a Jacobian over F_2-bar is a Jacobian over F_2, and Serre's sign cannot block anything. (For ECDLP purposes a twisted Jacobian would have been usable over F_4 anyway, but none exists.) |
| N4. Hyperelliptic sub-case | Serre's sign is vacuous for hyperelliptic curves. Over F_2 the g^1_2 descends to C → P^1. | **yes** | **Excluded (proved, elementary).** Each of the 3 points of P^1(F_2) contributes 2 rational points, 1 rational point (ramified), or 1 degree-2 place. a_2(C)=0 forces a_1 ≥ 3, but a_1=2. The genus-2 brute-force control (1920 weighted models over F_2, 320 of them with a_2=0) is never rejected. Independent consistency check: the étale A[2] has Frobenius characteristic polynomial Φ_131, which would force h(x) irreducible of degree 131, so the étale 2-torsion alone does not exclude this case. The fiber count does. |
| N5. Aut(C) and morphisms | Aut(Jac,θ) = Aut(C)×{±1}, or Aut(C) when C is hyperelliptic. | yes (lead 3) | Aut(A',λ)={±1}, so Aut_{F_2}(C)=1 (C is non-hyperelliptic by N4). A is simple over F_2, so C has no non-constant F_2-map to a curve of genus 1 to 129 (proved). Consistent with a Jacobian, not an obstruction. |
| N6. Rational symmetric theta divisor | A Jacobian in characteristic 2 has one over the base field: W_{g-1} minus the canonical theta characteristic. The torsor class lies in H^1(F_2, A'[2]). | **yes, on the explicit module** | **The obstruction vanishes.** P_A mod 2 = x^130·Φ_131 (computed), with Φ_131 irreducible mod 2. On the explicit lattice mod 2 (T/2T = O_K/2O_K because the index 263·p2 is odd), the étale dimension is 130 and F-1 is invertible on the étale part. So H^0 = H^1 = 0. A toy negative control with rational 2-torsion is detected as H^1 ≠ 0. |
| N7. Geometric structure | — | proved (Deligne over F_{2^131}, where Frobenius α^131 lies in Q(√-7); torsion-free modules over a quadratic order are sums of ideals, Bass / Borevich–Faddeev) | A'_{F_{2^131}} ≅ E_1×…×E_130 as an abelian variety, with E_i ordinary and isogenous to E. So C would be a genus-130 curve whose Jacobian is geometrically isomorphic to a product of elliptic curves, with an indecomposable principal polarization. That is an Ekedahl–Serre-type completely decomposable Jacobian; whether ordinary examples exist in genus 130 is not known to me. |
| N8. Geometric automorphism bound | For an ordinary curve, Nakajima (1987) gives \|Aut(C_{F_2-bar})\| ≤ 84g(g-1), so \|Aut(A'_{F_2-bar},λ)\| ≤ 2,817,360. A geometric automorphism of order 131 would also force a cyclic cover of P^1 with 4 branch points (Riemann–Hurwitz); that is the earlier cover-of-P^1 route, analysed only for F_2-rational symmetry. | **in principle (finite computation)** | **Not run.** Aut(A'_{F_2-bar},λ) is the centralizer of δ = π^131 − π̄^131 in O(T,Q), with Q(x,y) = b(x,δy). End over F_2-bar is the commutant of Q(√-7), so it is already reached over F_{2^131}. Q is positive definite of rank 260 with det 7^130·U^260; Plesken–Souvignier is beyond the budget. The mass-formula expectation is {±1}. |
| N9. Schottky / Andreotti–Mayer (dim Sing Θ ≥ g-4), theta-null relations, Serre's sign via χ18-type modular forms (explicit only for g=3) | — | no (needs equations or theta constants at g=130) | — |
| N10. Howe–Lauter resultant/gluing; Howe–Nart–Ritzenthaler | — | not applicable | The real Weil polynomial is irreducible over Q, and HNR is for g=2 only. |
| N11. Zarhin-type results | — | not applicable | Zarhin's trick (A^4×A^∨4) concerns the factor case. His endomorphism results are for characteristic 0 or supersingular varieties. |
| N12. Canonical lift / Coleman | — | — | No condition: the canonical lift of a Jacobian need not be a Jacobian (Oort–Sekiguchi, Dwork–Ogus). |

**Precise residue.** A non-hyperelliptic genus-130 curve C/F_2 with Aut_{F_2}(C)=1 and Jac(C) ≅ (A',λ), for some ppav A' non-maximal at v0. Geometrically, Jac(C) is isomorphic to a product of 130 elliptic curves isogenous to E. Over F_2 there is no descent obstruction. The only remaining computable necessary test I know of is N8, the geometric automorphism bound.

## 2. Quantitative heuristic
E[#Jac] ≈ N_ppav(class) · #M_g(F_2)/#A_g(F_2), with:
- #M_g(F_2) ≈ 2^{3g-3} = 2^387;
- #A_g(F_2) ≈ 2^{g(g+1)/2} = 2^8515;
- both are Lang–Weil leading terms for the stack counts.

**Estimating N_ppav.** Use N_ppav ≈ Σ_R \|Pic R\|/\|Pic R^+\| · \|S^+(R)\| over the conjugation-stable orders R that admit a ppav. This assumes every ppav with End=R occurs (lead 1, Theorem B) and that all orders are Bass. The sum is dominated by R0 = Z[π,π̄]: every conjugation-stable un-gluing above 263 or p2 divides the minus index by at least p. The R0 term is:

\|Pic R0\|/\|Pic R0^+\| ≈ h^-(K) · [O_K:R0]/[O_K+:R0^+] / ([O_K^×:R0^×]/[O_K+^×:R0+^×])

**Inputs, each cross-checked:**
- **h^-(K_131) is exact and unconditional:** 1128247071505617895582943540743356840380563482194685535482598575 ≈ 2^209.456.
  - Formula: h^- = Q·w·∏_{χ odd}(-B_{1,χ}/2), with w=262. Q=2 is forced because Q=1 gives a half-integer.
  - Controls: Q(ζ_p) for p = 23, 29, 31, 37, 41, 43, 47, 59 reproduces the table values 3, 8, 9, 37, 121, 211, 695, 41241.
  - The same code gives h^-(K_ell) = bnfinit h_K/h_{K+} for ell = 3, 5, 11, 13, 17, 19, namely 1, 1, 5, 7, 45, 148 (GRH; h^+=1 in each case).
  - h^-(Q(ζ_131)) = 28496379729272136525 divides h^-(K).
- **Indices, by discriminants (independent of lead 2's HNF):**
  - [O_K:R0]^2 = U^258 exactly;
  - [O_K+:Z[π+π̄]] = (263·p2)^64, so the minus index is U^65 (log2 4229.07);
  - N(π−π̄)^2 = 263^2·p2^2;
  - P_A recomputed as Res_z(Φ_131, x^2+xz+2z^2) equals P_A.json.
- **Unit ratio:** 262/[R0^×:±R0+^×] ∈ {131, 262}, i.e. -7.0 to -8.0 bits.
- **Check of the local minus factor (ell=17, GRH):** the ray class group gives \|Pic R\| = 360 for R glued at a degree-2 prime above 271. The formula h^-·(p+1)/(2ℓ) also gives 360.
- **Polarization units mod norms, \|S^+(O_K)\|:** 2^0, 2^0, 2^2 at ell = 11, 17, 19 (bnfinit). Allowed band [0, g-1] bits.

**Result.**
- log2 N_ppav(class) ≈ **4431.0**, band [4428, 4562]. Band components: Euler factors -1.5 to 0, norm cokernel ≤1, sum over orders ≤ +1, S^+ up to +129.
- Cross-check: the DiPippo–Howe count of q-Weil polynomials, v_g·q^{g(g+1)/4}, gives an average of 2^4388.6 ppav per class. Ours is 42 bits above that, so this class is typical in log scale.
- Equivalently, log2 N_ppav ≈ g(g+1)/4 + 173.5.

**log2 E[#Jac] = 4431 + 387 − 8515 = -3697**, band **[-3700, -3567]**.

| ell | g | log2 N_ppav | DiPippo–Howe average | log2 E |
|---|---|---|---|---|
| 11 | 10 | 21.0 | 22.1 | -7.0 [-10, 3.5] |
| 17 | 16 | 65.6 | 63.1 | -25.4 [-28, -9] |
| 19 | 18 | **0** (parity rule: no ppav) | 81.2 | **0 Jacobians** (negative control: the 82-bit class-number estimate must be gated by existence) |
| 131 | 130 | 4431 | 4388.6 | **-3697 [-3700, -3567]** |

The two methods agree at ell=11 and 17 to within 1 to 3 bits. Lead 1 found 4 ppavs with End=O_K at ell=11, against an O_K term of h^-=5.

### Honest error assessment
1. **The arithmetic inputs are solid** to within the band above.
2. **The M_g/A_g factor is the weak point.** Lang–Weil with q=2 and g=130 is not justified. The Betti numbers of M_130 are astronomically large: \|χ(M_g)\| ≈ \|B_{2g}\|/(2g(2g-2)) is about 2^1000. So the point count of M_130(F_2) could in principle differ from 2^387 by a huge factor. I have no rigorous upper bound on #M_130(F_2) better than ones far above 2^8515.
3. **The independence assumption is unproved.** The assumption is that Jacobians are spread across ppav isogeny classes like random points. This class is special: geometrically it is CM, isomorphic to a product of elliptic curves. Its ppavs, however, have End non-maximal and Aut={±1}, so they look generic within the class.
4. **How robust the conclusion is:** "E << 1" fails only if Jacobians are over-represented in this class by more than about 2^3566 relative to random.

**Classification:** a heuristic, not evidence of non-existence. The scaling is log2 E ≈ 3g − g(g+1)/4 + O(g).

## Files
- `zeta_tests.py`, `zeta_tests.json`: the N2 coefficient test with its detector control; point and place counts for P_A and P_A(-x) (the Serre twist); the hyperelliptic fiber test with the genus-2 control. Power sums are checked two ways, Newton identities against the Lucas closed form, for k ≤ 600.
- `two_torsion.py`, `two_torsion.json`: the N6 test on the explicit lattice mod 2, with a negative control.
- `hminus.gp`, `hminus.out`, `hminus_defs.gp`, `hminus131.gp`, `hminus131.out`, `hplus_fix.gp`, `hplus_fix.out`: analytic h^- with its controls.
- `index_plus.gp`, `index_plus.out`, `index131_exact.gp`, `index131_exact.out`, `PA131.txt`: plus and minus indices.
- `picR17.gp`, `picR17.out`: Pic(R) at ell=17.
- `splus2.gp`, `splus2.out`: \|S^+(O_K)\| at ell=11, 17, 19.
- `heuristic.py`, `heuristic.json`: the final numbers.

## Caveats (agent)

- **Cited, not re-derived:** Serre's descent theorem; uniqueness of the decomposition of a ppav into indecomposables in characteristic p; Nakajima's bound for ordinary curves; Bass/Borevich–Faddeev for modules over quadratic orders; the DiPippo–Howe constant v_g (formula as recalled).
- **Depends on lead 1:** Theorem B (every conjugation-stable order of the right parity occurs) is used for the N_ppav sum. It has not been independently reviewed.
- **GRH-conditional controls:** the bnfinit comparisons (h_K, h_{K+}, ray class group, S^+) assume GRH. The analytic h^- formula itself is unconditional.
- **Pic(R) at ell=17:** the value assumes F_p^* lies in the image of the units, which is consistent with the computed orders.
- **Untested necessary condition:** the geometric automorphism bound N8 was not computed.
- **Weakest part of the heuristic:** the Lang–Weil/M_g/A_g factor, which is not justified at q=2, g=130, together with the assumption that Jacobians are independent of isogeny class. The -3697 figure is a random-model estimate, not evidence that no Jacobian exists.
- **Scope:** the hyperelliptic and Serre-twist exclusions hold for curves over F_2, which by descent covers curves over F_2-bar too. Curves that contain A only as a Jacobian factor were not considered.
