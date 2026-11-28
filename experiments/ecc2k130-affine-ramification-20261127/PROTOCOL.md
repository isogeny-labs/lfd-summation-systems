# Frozen screen: tame ramification in the degree-131 dihedral route

## Question and scope

Let `E/F_2` be `y^2+xy=x^3+1`. Suppose `D/E` is a geometrically
connected Galois cover **over `F_2`** with full arithmetic and geometric
group the nontrivial dihedral group `G=C_131 semidirect C_2`. Let
`C=D/C_2` be the degree-131 intermediate curve, with `g(C)=131`, and
`E'=D/C_131`. Can the normal `C_131` layer ramify?

The full arithmetic group hypothesis is essential. A merely geometric
dihedral group can have a larger arithmetic normalizer, changing the
Frobenius action at branch places.

## Predeclared gate

1. At any place where inertia contains `C_131`, determine the possible
   inertia subgroups using the wild-inertia normality constraint.
2. Use the local tame Frobenius relation to constrain the degree of that
   place from the order of `2` modulo `131` and the conjugation action in
   `G`.
3. Compute the local Artin conductor of the 131-sheet permutation
   representation and compare its degree-weighted contribution with
   `deg Diff(C/E)=2g(C)-2=260`.
4. If ramification is excluded, combine only with the prior, separately
   stated unramified quotient result. Recheck whether its deck-group and
   constant-field hypotheses follow from the present hypotheses before
   making a combined claim.

Record a conditional exclusion if every implication holds. If a step
fails, retain the ramified case and state the countercondition. Do not
extend the result to other affine complements, larger arithmetic
monodromy, or arbitrary Jacobians.

## Verification

Check `ord_131(2)` and `2^65 mod 131` independently by integer arithmetic.
The written local-field and conductor proof is the actual gate; a finite
calculation only checks its small numerical inputs.
