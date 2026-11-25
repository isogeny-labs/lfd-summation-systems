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

## Continuation: solve the isogeny-class problem

The user explicitly requested a solution of the remaining isogeny-class
question. The full Jacobian verdict has not been obtained. The following
successor work was executed, rather than left as a proposed next step.

1. Proved an isogeny-degree obstruction: for every principally polarized
   target B/F2 and every isogeny f:A->B, the valuation of deg(f) at
   ell=146505763881528721 is odd. The proof uses local polarization
   kernel lengths and Artin reciprocity, with exact conductor inputs.
   Maximal-order representatives cannot be principally polarized.
2. Proved that every principally polarized member is geometrically
   indecomposable, so decomposability cannot close the remainder.
3. Excluded every hyperelliptic curve in the class by #C(F4)=2, since
   each rational point of its degree-two genus-zero quotient has a
   point above it over F4, requiring at least three points.
4. Constructed an explicit principally polarized Deligne module. A
   circular-unit signature computation has rank 129; an explicit
   product of 69 units corrects all 130 signs of the trial form.
   Rational interval arithmetic certifies the signs. Polynomial
   expansion gives the actual 260 unit coefficients. A cyclic quotient
   lattice of degree 263*ell carries a positive integral alternating
   form with determinant one by an exact norm/index identity. Frobenius
   compatibility and lattice stability are checked. Its endomorphism
   order has index ell in O_K.

Artifacts: docs/ecc2k130-isogeny-class-obstruction.md,
docs/ecc2k130-explicit-ppav.md, local_orders.py / local_orders.json,
ppav_unit.py, construct_ppav.py / explicit_ppav.json in the existing
ecc2k130-weil-factor experiment directory. The full principal matrix
is reconstructible; its SHA-256 (decimal rows) is
4ab7dd31504657e5ad4f180a96dbd9e5838d4cfaeb9c7fe7fccde26b988b94a4.

Verification: both new producers pass; independent full FLINT modular
factorizations match the Frobenius-gcd factor-pattern checks. FLINT's
integer polynomial discriminant matches the stated index identity.
The positive-unit signs are exact; an off-kernel generator negative
control is rejected. Existing independent SymPy certificate passes.
Repository tests after installing FLINT: 53 passed, 1 skipped, 16.09s.
No new Lean formalization was made.

Two optional exploratory computations were stopped without results:
a dense alternating determinant after about 84 seconds and a generic
HNF basis computation after about 63 seconds. The final certificate
uses exact unit/norm/index identities and an explicit cyclic basis;
it does not claim either abandoned computation succeeded.

The remaining failure to solve the original objective is specific:
no characteristic-two Torelli-locus certificate or curve realization
has been implemented for this module, and a verdict for this one
candidate would not classify all polarized representatives. ThetAV's
documented genus-2 conversions and the complex numerical Fay method
were investigated and do not provide the needed finite-field verdict.
No curve, pointwise cryptanalytic transfer, or ECDLP improvement is claimed.

### Additional verification and elliptic-map probe

The actual principal matrix has determinant one modulo five independently
computed prime moduli, supplementing the exact norm/index proof. A CSV
export is available through --matrix-out; its SHA-256 matches the record.
The isogeny degree is checked coprime to the target subgroup order.

A prime-degree elliptic-map obstruction was then tried on the constructed
polarization. Acceptance and a 120-second exact-Gram LLL cap were written
before running. That worker timed out without returning a result. After
recording this, an adaptive approximate-Gram mode with exact final
vector checks was given a separate 60-second cap and also timed out.
Both records are preserved, with no claim that the lattice has no short
vectors. No further unchanged LLL retry was made. The protocol explains
why a returned degree in {1,2,3,5,7,11} would have excluded this one
polarization, and why even that would not settle the full isogeny class.
