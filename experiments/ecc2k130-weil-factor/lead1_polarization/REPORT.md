# Lead 1: principal polarizations in the isogeny class of A_131

The lead-1 agent wrote this analysis; the harness blocked it from saving the report file itself.
The parent session reviewed it and added the integer-fact checks and the consistency table below.

## Result (ell = 131; proved without GRH, from Deligne + Howe + Hilbert reciprocity + class field theory)
1. **No principally polarized member contains zeta_131.** No ppav isogenous to A_131 has zeta_131 in End over F_2.
   In particular A_131 itself (End = O_K) is NOT principally polarizable.
2. **Every ppav is non-maximal at one place.** Its End is non-maximal at the unique self-glued place v0 of K+
   above p2 = 146505763881528721. Locally End = Z_p2 + p2*O_v0, so p2 divides [O_K : End].
3. **Ppavs exist, and the orders that occur are known.** Every conjugation-stable order R with
   Z[pi, pi-bar] <= R <= O_K that is non-maximal at v0 occurs as End of a ppav. Orders not stable under
   conjugation never occur. All orders here are Bass, so every Deligne module is invertible over its End.
   The impossibility proof does not assume invertibility.

**Consequence.** Every Jacobian isogenous to A_131 has no F_2-rational automorphism of order 131. So the
Torelli -> cyclic-cover-of-P^1 reduction (and the cover-trace evidence of PR #70) applies to no Jacobian
in this class.

## Criterion (re-derived)
- **Deligne modules:** lattices T in K stable under pi and pi-bar; End(T) = {x : xT <= T}.
- **Polarizations (Howe 1995, section 4):** b(x, y) = Tr(lambda * x * conj(y)), with conj(lambda) = -lambda
  and Im phi(lambda) > 0 for phi in Phi.
- **Principal:** lambda * conj(T) equals the trace dual of T (checked as det Gram = 1).
- **CM type:** for the canonical lift, Phi depends only on the embedding of sqrt(-7). Writing
  lambda = sqrt(-7) * r with r in K+, positivity means r is totally positive or totally negative
  (both signs were run).

## General rule (odd prime ell != 7)
Let eps = (-7/ell). Each prime p dividing U_ell has one of two types:
- **inert type:** p = -1 mod ell and (-7/p) = -1;
- **split type:** p = 1 mod ell and (-7/p) = 1.

Each inert-type p has exactly one self-glued place v0. Let n(T) be the number of such places where End(T)
is non-maximal.
- **Theorem A (necessary):** eps * (-1)^n(T) = +1.
- **Theorem B (sufficient):** every conjugation-stable order with that parity occurs as End of a ppav.

**Proof of A** (Hilbert reciprocity for (r, -7)_v over K+). The contributions are:
- the real places together: +1;
- places split in K/K+: +1;
- the inert place above ell, where the different has odd exponent: eps;
- each self-glued v0: -1 exactly when End is non-maximal there;
- each glued pair of inert places: +1.

At ell = 131, eps = -1, so End must be non-maximal at v0. If zeta were in End, End would contain
Z[pi, pi-bar, zeta] = O_K, a contradiction.

**Proof of B.** Class field theory for K/K+ (unramified at every finite prime), followed by a local twist
at the primes above U_ell.

## Integer facts at ell = 131 (rechecked by the parent session in gp)
- 263 and p2 are prime.
- p2 = 130 = -1 mod 131 and (-7/p2) = -1: p2 is of inert type.
- 263 = 1 mod 131 and (-7/263) = 1: 263 is of split type.
- (-7/131) = -1.

## Consistency of the rule with every control
| ell | eps | U_ell primes (type) | prediction | computed (GRH enumeration) |
|---|---|---|---|---|
| 3 | -1 | none | no ppav | none; unconditional: none of the 768 genus-2 curves over F_2 has the Weil polynomial |
| 5, 13 | -1 | none | no ppav | none |
| 11 | +1 | 23 (split) | O_K ppav | 5 of 5 classes (corrected; was recorded 4 of 5, see below) |
| 17 | -1 | 271 (inert) | ppav iff non-maximal at v0 | non-maximal at v0: 45/45 and 45/45 (corrected from 43/45); v0 maximal: 0 |
| 19 | -1 | 457 (split) | no ppav | 0 of 148 (O_K), 0 of 148 (glued) |
| 131 | -1 | 263 (split), p2 (inert) | ppav iff non-maximal at v0 | explicit module (PR #57) re-audited, see below |

## Audit of the existing ell = 131 module (PR #57, `../construct_ppav.py`)
Rerun in memory and checked with independent arithmetic:
- the form is Tr(lambda * x * conj(y)) on all 67,600 entries;
- conj(lambda) = -lambda, and positivity holds at all 130 embeddings;
- det = 1 exactly;
- the lattice is stable under pi and pi-bar but not under zeta;
- it is non-maximal exactly at v0 above p2, and maximal at 263.

Reproduce with `dump_prior_ppav131.py` and `audit_prior_ppav131.gp`.

## Scope
- **Proved:** the ell = 131 statements above (written proof by the agent; the parent session checked its
  consistency, not every line).
- **Computed under GRH:** the small-ell enumerations.
- **Not addressed:**
  - whether any ppav in the class is a Jacobian;
  - automorphisms defined only over extensions;
  - curves whose Jacobian contains A_131 only as a factor.

## Untested prediction
At ell = 29 and 37, eps = +1 and there are inert-type primes. The rule says O_K is principally
polarizable there, but orders non-maximal at an odd number of v0 places are not. Those fields are too
large to compute here.

## Files
- `local_data.*`, `local_structure.*`
- `decide_ppav.gp`, `run_small.*`, `run_17_19.*`, `run_17_orders.*`
- `genus2_F2_control.*`
- `dump_prior_ppav131.py`, `audit_prior_ppav131.*`

## CORRECTIONS (2026-11-28, from round-2 red-team and the ell=29/37 lead)
- **Bug in `decide_ppav.gp`.** Lines 68 and 76 (`if (sol == 0, next)`, `if (x == 0, next)`) discard a valid all-zero
  solution, because in GP `[0,..,0]~ == 0` is true. This undercounts.
- **Corrected counts.** Two independent agents found the bug. With it fixed:
  - ell=11: O_K 5/5 (was 4/5).
  - ell=17: 45/45 (was 43/45).
  - ell=19: the antisymmetric counts are 148/148.
  No existence verdict changes.
- **Hypothesis of the general rule.** As stated, the rule assumes v_p(U_ell) = 1 for the inert-type primes. If p^2 divides
  U_ell, n(T) must count self-glued places with odd local conductor exponent. The hypothesis holds at ell=131
  (U_131 is squarefree).
- **Upheld at ell=131.** The independent red-team re-derivation (`../round2/redteam-lead1/REPORT.md`) upholds Theorems A
  and B. Its process note: this file is a sketch, not a full written proof. There is no Lean proof.
- **What "A_ell itself" means at ell=29, 37** (`../round2/prediction-29-37/REPORT.md`):
  - Principally polarized varieties with End = O_K exist there, as the rule predicts.
  - The specific module T = O_K (the natural A_ell) has none, for a narrow-class reason: the ramified prime above ell in
    K+ has no totally definite generator.
  - So "End = O_K admits a ppav" and "A_ell itself admits one" are different statements.
