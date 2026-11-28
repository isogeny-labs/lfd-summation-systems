<!-- Agent report saved by the parent session from the workflow journal. Status: mixed. Independent skeptic verification: NOT YET RUN (failed on usage limit). -->

# Testing lead 1's prediction at ell = 29 and ell = 37

**Summary.** The lead-1 parity rule holds at both primes, and the constructive positives are unconditional. One part of the prediction as stated was wrong: A_ell itself (the module T = O_K) is **not** principally polarizable at either ell. Some other module with End = O_K is. The cause is a finer narrow-class obstruction that the parity rule does not see.

All files are in `experiments/ecc2k130-weil-factor/round2/prediction-29-37/`.

## 0. Integer data (`u_factor.gp/.out`)
Two independent computations of U_ell (a recurrence and a computation in Q(alpha)) agree.

| ell | eps = (-7/ell) | U_ell | primes and types |
|---|---|---|---|
| 11 | +1 | 23 | split |
| 17 | -1 | 271 | inert |
| 29 | +1 | -8641 | 8641: inert (8641 = -1 mod 29, (-7/8641) = -1) |
| 37 | +1 | 194399 = 73 * 2663 | both inert |

So ell = 37 has two self-glued places, v0 over 73 and v0' over 2663. The rule allows n in {0, 2} there and forbids n = 1.

## 1. Reduction that needs no class group of K
- **K+ is cheap to work with.** K+ is cyclic of degree ell-1 and conductor 7*ell. Its discriminant is 7^((ell-1)/2) * ell^(ell-2); this matches the computation and nfcertify returned nothing to check. K/K+ is unramified at every finite prime.
- **Shape of lambda.** On T, lambda * T-bar = T^dual. Writing lambda = sqrt-7 * r with r in K+, the CM type for the canonical lift forces r to be totally definite.
- **The different.** D = (sqrt-7 * (z - 1/z)^(ell-2)), and qO_K = (z - 1/z), where q is the totally ramified prime of K+ above ell.
- **T = O_K.** The ideal of r is 7^-1 (c^2-4)^(-(ell-3)/2) q^-1, where c = z + 1/z. Since c^2-4 < 0 at every real place, **A_ell itself is a ppav iff [q] = 0 in Cl+(K+)** (proved).
- **General T.** T = S with End = O_K, or an order non-maximal at a set of places v0, needs q * prod(v0) * N(S) to be narrowly principal. By class field theory, N(Cl_K) is the kernel of the Artin map of K/K+, so some S works iff the Artin symbol eps*(-1)^n is trivial. This is lead 1's rule, Theorems A and B, now with a one-line proof.
- **Local computation at v0.** Let R_v = Z_p + p O_v. Then R_v^dual = (sqrt-7/p) R_v, so v(r) = -1 at v0. That gives (r, -7)_v0 = (p, -7)_p = -1. At a maximal place the symbol is +1.

## 2. Results

| test | prediction | result | status |
|---|---|---|---|
| ell=11, T=O_K (control) | ppav | [q]=0 in Cl+=Z/2; explicit det 1, all positive | computed, unconditional certificate |
| ell=29, End=O_K, some T | ppav | T = P above split p=7193 with [qv]=0: det 1, integral, alternating, conj(lambda) = -lambda, stable under pi and pibar, r totally positive (nfeltsign and numeric Im phi(lambda) on all 28 phi) | **confirmed constructively, unconditional** |
| ell=37, End=O_K, some T | ppav | T = P above p=149: same checks, 36 phi | **confirmed constructively, unconditional** |
| ell=37, n=2 (both v0) | ppav | T = R itself (index 194399): det 1, positive, zeta not in End | **confirmed constructively, unconditional** |
| ell=29, n=1 | no ppav | Artin/Hilbert obstruction; nfhilbert gives -1 only at v0, and the generator has 11 negative embeddings (odd); [q v0] is not in the norm group (GRH) | **proved (unconditional), checked locally** |
| ell=37, n=1 (v0 over 73, or over 2663) | no ppav | same; -1 only at v0; 13 and 23 negative embeddings | **proved, checked locally** |
| ell=29, 37, **T = O_K (A_ell itself)** | the task said ppav | **NOT a ppav**: [q] = [1,1,0,1,1] in Cl+ = (Z/2)^5 at 29, and [1,1] in (Z/2)^2 at 37 (GRH). Unconditional: bnrclassfield gives quadratic L/K+ with rnfdisc of norm 1 in which q does not split (4 of 5 at ell=29, 2 of 2 at ell=37). | **proved unconditionally**; it corrects the prediction's wording, not the rule |

## 3. Notes
- **Negative controls.** In every negative case, the lambda built from a non-definite generator still gives an integral, alternating Gram with det 1; only positivity fails. This shows the detector separates the two conditions.
- **Product-formula check.** The (r, -7) symbols multiply to 1 in every case (`run_local.out`). For the n=2 construction, the symbol is -1 at both v0 and v0' and r is totally positive.
- **The rule needs a caveat.** Lead 1's rule speaks about End(T) only. It does not decide whether the particular module A_ell (T = O_K) is a ppav. That depends on the narrow class of q, which is trivial at ell=11 and nontrivial at ell=29 and 37. The parity rule is *consistent* at 29 and 37 because some S does work; it is not contradicted.
- **Bug in lead 1's `decide_ppav.gp`.** `matsolvemod` returns 0 when there is no solution, but an all-zero solution vector also tests `== 0` in PARI, so valid classes were skipped. The fix tests `type(sol) == "t_INT"`, in `decide_ppav_fixed.gp`.
  - Reruns: `rerun_fixed_run_small/17_19/17_orders.out`, plus `ctrl11_decide*.out`.
  - The ell=11 count becomes 5/5 (not 4/5).
  - The ell=17 order counts become 45/45 (not 43/45).
  - The ell=19 antisym counts become 148/148.
  - Every existence verdict and the 'agreement with prediction: 1' line are unchanged.
  - The lead-1 table's '4 of 5 classes' and '43/45' should be corrected.

## 4. Proof status
- **Proved:** the reduction in section 1; T=O_K not a ppav at 29 and 37 (via the unramified quadratic L, verified by rnfdisc and rnfidealprimedec); no ppav for n=1 at 29 and 37 (Hilbert reciprocity).
- **Computed, unconditional:** the constructive certificates. These are exact rational Gram determinants; positivity uses nfeltsign, cross-checked numerically at 300 digits.
- **GRH:** only the bnfinit, bnrinit and narrow-class computations used to *find* generators, the Cl+ structures, and the 'unit correction impossible' outputs.
- **Not done:** an independent human review of the Artin/CFT argument; whether any of these ppavs is a Jacobian. Ideal classes S at ell=29 and 37 were not enumerated beyond the single constructed one.

## Caveats (agent)

The CFT and Hilbert-reciprocity arguments were written by this agent and have not been independently reviewed. Positivity is certified by PARI's nfeltsign and a numeric embedding check at 300 digits, not by interval arithmetic. The narrow class groups, the unit-sign impossibility claims, and the search that found the generators are GRH-conditional. The final certificates (ideal equalities, exact det = 1, rnfdisc) do not depend on GRH. The statement that A_ell has Deligne module O_K is taken from the PR #57 docs argument. Nothing was committed, and no files outside round2/prediction-29-37 were edited. The decide_ppav.gp bug is reported here, not fixed in place.
