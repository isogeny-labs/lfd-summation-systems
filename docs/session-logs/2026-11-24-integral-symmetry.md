# 2026-11-24: ECC2K-130 integral symmetry (Codex, PR #57)

## Objective and scope

Continue the Jacobian question beyond PR #53's bounded zeta tests.
Acceptance for the original objective remains a curve with a usable
transfer or a proof excluding Jacobians in the relevant isogeny class.
Neither has been obtained. This draft records verified intermediate results.

## Results

- The concrete trace-zero kernel is not F2-isomorphic to a genus-130
  Jacobian: Rosati and Torelli would give a rational order-131 curve
  automorphism; Riemann--Hurwitz requires four branch points, incompatible
  with Frobenius branch orbits of length divisible by 130.
- The isogeny class contains a principally polarized variety by Howe's
  Theorem 1.1. The certificate gives an inert prime above
  146505763881528721 dividing F-V. This is existence of a polarization
  somewhere in the class, not a curve construction.
- Two reductions modulo 263 agree on F,V but not zeta_131. Thus the
  cyclotomic symmetry is not an automatic integral consequence of F,V.
- Point, place and effective-divisor consistency holds in every degree,
  by finite exact checks and proved tail inequalities/identities.
- PR #53's N*g estimate for linear algebra counts one sparse multiply,
  not a complete generic solve; the report identifies the correction.

## Verification

`integral_symmetry.py` uses standard-library exact integers and a
deterministic 64-bit primality test. The saved JSON reproduces exactly.
`verify_integral_symmetry.py` independently checks the field witness,
factorization and Frobenius-polynomial gcd with SymPy, and verifies the
small-degree counts with Newton/Mobius identities and a different
effective-divisor recurrence. Both pass.

`python -m pytest -q tests`: 53 passed, 1 skipped (16.88 seconds).
Python compilation and `git diff --check` pass. The existing CI workflow
now replays both arithmetic scripts and compares the saved JSON.
No new Lean theorem was added; the geometric proofs are written proofs
using the cited Rosati, Torelli, and Howe results.

## Remaining obligation and successor

The automorphism argument closes the concrete representative, not every
variety with its Weil polynomial. Principal polarizability is positively
settled and cannot be the missing exclusion. The remaining task is to
construct or exclude principally polarized Deligne modules with smaller
endomorphism orders lying on the Torelli locus. No algorithm handling
that dimension has been implemented in this repository. The exact next
construction target is a principally polarized module exploiting the
inert-prime witness, then an actual curve or a Torelli obstruction.
This is a mathematical/construction gap, not a request for user permission
or an assertion that all other routes have been exhausted.
