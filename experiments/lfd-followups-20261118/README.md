# Concrete first-derivative operators, 18 November 2026

Continuation of draft PR #31; all work remains local. No paid resources
were launched. This work does not resolve the full depth-14 certificate,
split degree-8 cell, depth-18 breakpoint or asymptotic proper-ideal bound.

## New Lean results

`formal/ICFallDegrees/GrowthJetClosure.lean` defines actual linear operators
on coefficient arrays A -> Fin N -> K:

- Euler derivative: multiply the distinguished exponent-b coefficient by b;
- cyclic shift: pull the coefficient from (b+N-s) mod N;
- wrapping projection: keep source exponents at least N-s.

It proves the exact commutator at every N>0 and s<N, over any commutative
ring, then specializes it to two exponent blocks N=2^k, s=2^j, j<k.
This replaces the finite commutator tests with a kernel-checked identity.
No inverse of N or characteristic restriction is needed for this identity.

It also proves that Euler differentiation and wrapping preserve any
support-weight budget, and that shifts in the other coordinate commute
with the Euler operator. Taking the weight to be the sum of binary
popcounts gives the intended chain degree; the theorem is stronger in
that it accepts any weight function without increasing support.

For a submodule S, define its first-jet subspace by f in S and delta(f) in S.
The conditional `jetSpace_shift` theorem proves shift stability under
explicit ordinary-shift and wrapping-shift stability assumptions on
low-weight rows. A derivative outside S excludes a witness from that
first-jet subspace. None of these statements supplies the analytic
wrapping hypotheses or the embedding into the existing polynomial
normal-form/FEI development.

## Falsification of the unqualified statement

`check_jet_invariant.py` checks completed small closures by evaluating
unwrapped rows and their Euler derivatives at rational curve samples.
It finds counterexamples at k=2,D=3 and at D=4 for k=2..6. A nonzero
sample certifies failure; a list of zero samples proves no identity.
The old numerical exact values are not revised by these tests.

The k=2,p=17,D=3 counterexample is accompanied by a derivation DAG.
`check_certified_falsifier.py` checks all output derivations with the
separate pure-integer checker, then identifies a derived row with
phi=8 and phi(delta)=14 at x=2,y=10. The input is the normalized summation
polynomial for E6 and target (1,0); the chain relations are included via
the previously justified degree-compatible normal-form model.

`certified-falsifier.json` fixes the expected output digest outside the
certificate. The CI test reconstructs the curve input, checks that digest,
checks the derivation DAG and recomputes both nonzero evaluations.
This refutes dropping the wrapping hypotheses. It does not refute a
large-depth theorem with those hypotheses established.

## Remaining bridge

The next proof must identify these coefficient operators with expansion
and chain normal form in GrowthReduction, then derive the ordinary and
wrapping shift hypotheses from FEI at the required characteristic/depth
threshold. It must also bind the finite-chain Frobenius cofactor witness
to the ordinary input ideal. The new conditional subspace theorem alone
is not an LFD lower bound. These gaps remain explicit in PR #31.

## Reproduction

From the repository root:

```
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python3 -m pytest -q tests
python3 experiments/lfd-followups-20261118/check_certified_falsifier.py
python3 experiments/lfd-followups-20261118/check_jet_invariant.py
cd formal
lake env lean ICFallDegrees/GrowthJetClosure.lean
```

The two scripts rewrite their local output records; use a separate checkout
when preserving the saved records. The direct Lean check prints every new
theorem's axioms. The all-module CI checker automatically discovers this
module even though the older default root does not import it.
