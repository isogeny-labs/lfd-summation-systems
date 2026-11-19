# E6 growth beyond the unit-ideal argument

18 November 2026. Continuation after merged PR #31, tracked in draft PR #33, local-only work.

## New result and exact scope

`GrowthJetPolynomial.e6_prime_lfd_linear` proves:

For every prime p, depth k>=38 and 2^(k+1) dividing p-1, no
last-fall-degree bound e<=floor((k-10)/14) holds for the actual normalized
E6 chain system of y^2=x^3+6x-7, target (1,0).

The signature has NO unit-ideal, CurveModel, FEI or spread premise. The
proof composes the native prime-field bounded-box spread theorem already
in the repository with the new polynomial first-derivative invariant and
an actual Frobenius ideal witness. This includes exceptional proper-ideal
primes satisfying the stated divisibility. It is a fixed E6 result, not
an extension to all rational parameters or to non-two-torsion targets.

`e6_prime_lfd_gt` gives the degree form: D>=2, k>=14D+10,
p>=2*2^k+1, and 2^k dividing p-1 imply that no e<=D is an LFD bound.
The linear corollary derives both arithmetic guards from admissibility.

This file records a local Lean check. A final CI run on the exact pushed
commit must pass before reporting independent Linux verification. The
new source is automatically included by scripts/lean_check_all.py even
though it is not imported by the older default root.

## Why the unit-ideal restriction is removed

Let h=xy-x-y-8. The unwrapped summation polynomial is exactly h^2.
Consider rows whose lift and whose polynomial Euler derivative both
vanish in the rational-function model y=(x+8)/(x-1).

The polynomial Euler map is X_i times the formal partial derivative. Its
coefficient on exponent d is d_i times the original coefficient, so it
preserves support and the boxed binary-weight bound. It commutes with
the existing topPart projection. The new proof uses the existing cleared
lift identity directly; differentiating it and applying FEI to both the
original and differentiated lift shows that the jet space is closed under
every admissible chain shift. Thus actual InClosure rows lie in that space.
This is the concrete normal-form bridge missing from the earlier array
prototype; a separate wrapping-vanishing premise is no longer supplied.

In the finite chain quotient, 2^k dividing p-1 makes each active coordinate
Frobenius-fixed. The factor is therefore fixed too: h^p=h. The normalized
input and chain relations imply h^2=0 in the input quotient, hence h=0
there. This proves membership of the degree-2 witness h in the ACTUAL
input ideal. It does not claim a low-degree derivation of that membership.

On the rational-function model, the derivative of h is y*(x-1)=x+8,
which is nonzero. Consequently h cannot occur in a closure satisfying the
jet invariant. If any e<=D were an LFD bound, its degree-2 ideal membership
would derive h by degree max(e,2)<=D, a contradiction. This argument works
whether the finite input ideal is proper or the unit ideal.

## Proof layers and obligations

- `coeff_euler`, `boxed_euler`, `euler_topPart`, `euler_mul` establish the
  concrete polynomial identities.
- `derivative_shift_mem` and `jetSpace_closed` derive the shifted jet
  invariant from FEI for the actual polynomial normal form.
- `closure_normalForm_mem_jetSpace` applies the existing quotient bridge
  to actual InClosure.
- `lfd_gt_of_derivative_witness` is a general conditional witness theorem.
- `e6Witness_mem_ideal` establishes witness membership, rather than
  assuming it in the final E6 theorem.
- `e6_prime_lfd_gt` discharges FEI using native bounded-box spread.
- `e6_prime_lfd_linear` is the final arithmetic corollary.

The general witness theorem retains FEI and witness hypotheses. The final
E6 theorem discharges them. Do not describe the general lemma as
unconditional, or the E6 corollary as an arbitrary-curve theorem.

The old small-depth counterexamples do not contradict the new result:
they are below its depth threshold and demonstrate why FEI cannot be
removed from the intermediate jet-closure statement.

## Verification

The saved local Lean log has no errors, placeholders or extra axioms.
All printed declarations use only propext, Classical.choice and Quot.sound.
The exact source digest and parsed axiom report are in local-verification.json.

check_algebra.py independently checks the symbolic S3=h^2 identity,
the nonzero rational derivative and product identity, plus five finite
Frobenius/proper-grid examples and two deliberate algebra changes.
These finite regressions are not an asymptotic closure replay. The
asymptotic statement is justified by the Lean proof, subject to final CI.

Reproduction from the repository root:

```
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python3 -m pytest -q tests
python3 experiments/lfd-jet-polynomial-20261118/check_algebra.py
cd formal
lake env lean ICFallDegrees/GrowthJetPolynomial.lean
```

The script rewrites its output record; use a separate checkout to preserve
original records. No full depth-14 basis, split degree-8 or depth-18 closure
was run. Those computational obligations remain open. The submitted
paper/PDF package has not been regenerated to include this new theorem.

## Integration with current main

After this proof was written, current main was integrated through #32.
The integrated local suite has 50 passing tests and four failures in the
root-finding tests: the C kernel includes x86-only intrinsics and does not
compile on this ARM Mac. This is not a Lean failure. The full Linux CI
suite is still required; do not report a locally green integrated suite.
The 439 submission consistency checks pass on the integrated tree.
The new proof does not regenerate the newer paper package from #32.
