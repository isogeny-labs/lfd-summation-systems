# Compact closure candidate: not yet implemented or certified

The current exact-pivot engines store a separate row for every retained
monomial pivot, including many monomial multiples of earlier low-degree
rows. Both split runs reached roughly 300 million terms at only 107,024
retained rows. Increasing that budget alone is unlikely to be sufficient
for the larger depth inputs on this machine.

The next candidate is a degree-bounded, divisor-based polynomial basis.
This is a proposed algorithm, not evidence of any open cell's exact value.

Represent normal monomials by chain masks and auxiliary exponents. Reduce
a leading monomial by a prior polynomial when its chain mask contains the
prior leading mask and its auxiliary exponents dominate the prior ones.
The multiplier is then disjoint from the prior leading chain mask; its
product has degree at most D. Reconstruct the full product in the cyclic
quotient, rather than simply shifting term indices without carries.

Candidate closure obligations:

1. Every retained polynomial has an independently replayable derivation in
   V_D. Reducers must record their monomial multiplier, coefficient and
   source row; the verifier checks deg(source)+deg(multiplier) <= D.
2. Reduce critical pairs whose leading monomial lcm has degree <= D.
3. Include carry pairs: multiply a retained polynomial by every chain
   variable present in its leading mask, when its degree is below D.
4. Prove that these obligations imply closure of the *whole truncated
   span*, including low-degree cancellations. Do not infer this from the
   producer reaching an empty queue or from ordinary ideal membership.
5. For a nonunit result, count the leading-monomial upward closure through
   D and prove it equals the span dimension. For a unit result, replay the
   constant polynomial directly. An exact lower bound still needs the
   compact closedness criterion to be proved and independently checked.

Before scaling, compare against the direct finite-space closure on saved
chain and split gates and exhaustive small/random polynomial controls.
Include cancellations that lower degree, the final cyclic wrap, auxiliary
exponents, pairs crossing the degree boundary, and a deliberately missing
critical/carry pair. A counterexample blocks use of this method for claims.

For depth14, bind any verified replacement space to the original input and
lift. A possible rank sandwich additionally needs the *actual archived*
upper space's checked containment of V_6, its exact input identity and
basis hashes. Matching two rank numbers alone does not establish equality.
The second depth14 prime and the frozen depth18 D7/D8 instance remain full
requirements; a compact prototype or a capped search does not replace them.
