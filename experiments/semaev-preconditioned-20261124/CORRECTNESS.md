# Trace constraints and exact affine preprocessing

**Status:** written proof and exhaustive finite checks, not kernel-checked.
These are standard trace and affine-elimination ideas, not a new asymptotic
attack or a new root-finding algorithm.

## Trace equations on the complete PDP domain

Let k=GF(2^n), E:y²+xy=x³+ax²+b nonsingular, and define
phi(O)=0 and phi(P)=Tr(x(P)+a) for affine P in E(k).

**Lemma.** phi is a group homomorphism E(k) -> GF(2).

**Proof.** For distinct abscissas, the addition formula is
x(P+Q)=lambda²+lambda+x(P)+x(Q)+a. Since Tr(lambda²+lambda)=0,
it gives phi(P+Q)=phi(P)+phi(Q). If P=Q and 2P is affine,
x(2P)=lambda²+lambda+a, so phi(2P)=0=2phi(P). If Q=-P,
then P+Q=O and the two abscissas agree, so both sides are zero.
The point at infinity is immediate. These cases exhaust the group law. ∎

For P1+P2+P3=R, the lemma gives

    Tr(x1+x2+x3+r)=0.

For an affine summation-polynomial triple, it gives

    Tr(x+y+z+a)=0.

Thus the finite-u split can receive
Tr(x1+x2+u+a)=Tr(u+x3+r+a)=0. The first two coordinates lift
to rational points by membership, and the S3 property forces the point
over u to be rational too; the sign choice does not change phi.

These statements use rational-point membership. They are not asserted for
the raw resultant zero set over arbitrary coordinates. The complete PDP
correspondence, including the infinity branch, is proved in the
[prior experiment](../semaev-matched-20261124/CORRESPONDENCE.md).
The split's infinity branch is recovered separately and is never subjected
to the finite-u trace equations.

## Boolean row span and substitution

Work in B=GF(2)[z1,...,zN]/(zi²+zi), representing polynomials as sets of
squarefree monomials. For every Boolean assignment, row operations over
GF(2) preserve the common zero set of a list of equations. Echelonization
with high-degree columns first exposes affine members of their row span.

**Proposition.** Replacing a system containing zj+g=0, with g affine and
independent of zj, by the remaining equations after zj=g gives a bijection
of Boolean solution sets. The inverse map inserts zj=g.

**Proof.** Every original solution has zj=g and satisfies the substituted
equations. Conversely, every residual solution extends uniquely by that
formula and satisfies every original equation. Multiplication of monomials
uses union of supports, exactly the identities zi²=zi, so the sparse
substitution is equality in B. ∎

The implementation repeats this operation, composing an affine recovery
map for every original variable, and compacts the remaining indices. It
eliminates one variable per successful iteration, so terminates in at most
N eliminations. A constant 1 proves inconsistency. If no equations remain,
the residual Boolean cube is free and is enumerated in recovery. No oracle
or assignment search is used to obtain the preprocessing equations.

Affine substitution does not increase the Boolean degree of an equation:
a monomial containing zj becomes its other factors times an affine
polynomial. This does not prove a decrease in solving degree or runtime.
The latter is tested with charged preprocessing and recovery.

The residual equations go to the same ordinary GF(2) Buchberger backend
with Boolean field equations. The bijection therefore also justifies
recovery from its residual Gröbner basis. Projection and the separately
recovered infinity branch give the full ordered PDP x-coordinate tuples.

## Independent controls

- 420 seeded generic Boolean systems, 1 to 7 variables: enumerate every
  original and residual assignment and check exact bijective recovery.
- 140 satisfiable systems with planted affine maps hidden by invertible row
  mixing, 2 to 8 variables: force nontrivial elimination and check full
  bijective recovery, rather than relying only on random dense systems.
- Explicit hidden-affine and inconsistent systems.
- Curve addition trace identity, including doubling and inverse cases.
- Complete small-field PDP oracles versus trace augmentation and reduction.
- Fast shared input builder versus the audited builder: identical equation
  lists, and independent full coordinate scans versus lift-generated points.
- All 45 frozen tuple oracles versus existing LRF plus the explicit infinity
  correction; LRF is an existing comparator, not the candidate contribution.
