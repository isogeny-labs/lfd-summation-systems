# Exact counterexamples to identifying first and last fall degree

## What is disproved

There is no unspecified numerical "first/last fall degree of elliptic-curve
index calculus": the field, polynomial encoding, factor base, target, and
definition of first fall must be fixed. Here we retain this repository's
**direct prime-field PKM coordinate-subgroup chains**, with one normalized
S3 equation and no additional filters or auxiliary variables.

For this encoding, the following universal statements are **false**:

1. First fall degree always equals last fall degree.
2. Last fall degree is always at most first fall degree.
3. Normalized two-summand systems always have last fall degree at most 4.

An exact counterexample has homogeneous-syzygy FFD = 4, affine first strict
fall = 4, and LFD = 5. Thus the distinction between those two FFD conventions
does not avoid the counterexample. It uses a nonsingular curve and a prime-order
target outside the factor base. A second example has FFD = 3 and LFD = 4.

These are **computer-assisted finite disproofs**, with exact finite-field
certificates and independently implemented ordinary-ring checks. They do not
disprove the possibility of a larger constant LFD bound, establish an
unbounded family, or settle a probabilistic/asymptotic heuristic. No novelty
claim is made. Binary curves and other encodings are not covered.

## 1. Definitions and the system

We use the definitions and symbolic theorems in
[IC_FALL_PROOFS.md](IC_FALL_PROOFS.md). In particular, V_(F,D) is the smallest
vector space containing input polynomials of degree at most D and closed
under multiplications whose ordinary polynomial degree is at most D.
**Newly fallen rows can be multiplied.** The last fall degree is the least d
such that every h in (F) belongs to V_(F,max(d,deg h)).

The affine first strict fall is the least D for which

\[
 \dim(V_{F,D}\cap R_{\le D-1})>\dim V_{F,D-1}.
\]

Homogeneous-syzygy FFD means the first non-Koszul homogeneous syzygy of
the top forms. The finite-field truncated-ring convention agrees with the
ordinary-ring calculation in the degrees used here because p>4.

Set

\[
 p=2689,\qquad E:y^2=x^3+1739x+1574,\qquad R_0=(797,1184),\qquad k=6.
\]

The curve is nonsingular. Its group has order 2614=2*1307. The target is
nonzero and [1307]R_0=O, so it has exact prime order 1307.

Take the coordinate subgroup V={x in F_p^*: x^64=1}. It has size 64 because
64 divides p-1. The target abscissa 797 does not belong to V. Introduce twelve
variables x_0,...,x_5,y_0,...,y_5 and the input chains

\[
 x_j^2-x_{j+1},\quad y_j^2-y_{j+1}\quad(0\le j<5),\qquad
 x_5^2-1,\quad y_5^2-1.
\]

Let u=x_0, v=x_1, w=y_0, z=y_1. The thirteenth input equation is exactly
the chain normal form of S3(x_0,y_0,797):

\[
 g=vz+1095(uz+vw)+690uw
              +605(v+z)+2164(u+w)+1447.                       \tag{1}
\]

All coefficients and equations are over F_2689. All thirteen input
polynomials are quadratic. Independent evaluation on the complete 64-by-64
coordinate grid finds no zeros of S3, in agreement with the unit-ideal
certificate below. A failed two-summand relation attempt is a valid instance;
its solving cost is part of relation collection.

## 2. Exact first-fall degree

Section 4.4 of [the proof note](IC_FALL_PROOFS.md) proves that for normalized
m=2 systems with k>=2,

\[
 \operatorname{FFD}_{\rm hom}=3\text{ if }3t^2+a=0,
 \qquad \operatorname{FFD}_{\rm hom}=4\text{ otherwise}.
\]

Here (3*797^2+1739) mod 2689 = 865, which is nonzero. Hence

\[
                         \operatorname{FFD}_{\rm hom}=4.        \tag{2}
\]

This is a proved symbolic formula, not a prediction from the earlier ML
model. Its proof uses the active four-variable square-zero quotient: the
linear-annihilator map for q=690uw+1095(uz+vw)+vz is invertible, while uv is
a nonzero quadratic annihilator. Thus a nontrivial syzygy first occurs in
total degree 2+2=4.

## 3. Why the quotient certificate computes ordinary-ring closures

Let G be the twelve input chain equations, J=(G), and N their normal-form
map. Their pairwise-coprime leading squares make G a degree-compatible
Groebner basis. In particular, for every h,

\[
 \deg N(h)\le\deg h,\qquad h-N(h)\in V_{G,\deg h}.              \tag{3}
\]

The standard monomials are squarefree. Write A_(<=D) for their span through
degree D, and let U_D be the smallest subspace of A_(<=D) that contains g
when deg(g)<=D and satisfies

\[
 u\in U_D,\quad\deg u<D
       \quad\Longrightarrow\quad N(z_i u)\in U_D               \tag{4}
\]

for every chain variable z_i. The condition is checked **before** reduction;
a degree-D row is never multiplied and then made affordable by a carry.

**Lemma.** In the ordinary polynomial ring,

\[
 V_{F,D}=(J\cap R_{\le D})\mathbin\oplus U_D.                  \tag{5}
\]

**Proof.** Equation (3) puts all of J intersect R_(<=D) in V_(F,D). Every
normal-form generator and every step (4) can also be derived in V_(F,D),
because reduction uses only input chain multiples of degree at most D.
This proves the right-hand side of (5) is contained in V_(F,D).

Conversely let W be that right-hand side. It contains every available input.
For h in W of degree less than D, write h=j+u using u=N(h). Both j and u
have degree less than D by degree compatibility. Then z_i*j lies in J
intersect R_(<=D); (4) puts N(z_i*u) in U_D; and the difference
z_i*u-N(z_i*u) lies in J intersect R_(<=D) by (3). Thus z_i*h belongs to W.
Closure under all allowed polynomial multipliers follows by iteration and
linearity. Minimality of V_(F,D) proves the reverse inclusion. The sum is
direct because a normal form in J is zero. QED.

A degree-compatible echelon basis computes the lower-degree subspace U_D
by its lower-degree pivot rows. A combination of higher-degree pivot rows
cannot cancel all their distinct leading pivots. Consequently the algorithm
can apply (4) to those rows, then echelonize and repeat. A strictly increasing
rank is bounded by the number of standard columns, so no heuristic stopping
rule or fixed iteration cap is needed.

There are c_D=sum_(j=0)^min(12,D) binomial(12,j) standard columns. Formula
(5) gives the ordinary-ring dimension

\[
 \dim V_{F,D}=\binom{12+D}{D}-c_D+\dim U_D.                    \tag{6}
\]

This proves the method computes the full V space, rather than only a static
Macaulay space or an affine-linear stopping proxy.

## 4. Lower and upper certificates: exact LFD 5

Exact elimination yields:

| D | Standard columns | dim U_D | dim V_(F,D) | dim(V_(F,D) intersect R_(<=D-1)) | 1 in V_(F,D)? |
|---:|---:|---:|---:|---:|:---:|
| 0 | 1 | 0 | 0 | 0 | no |
| 1 | 13 | 0 | 0 | 0 | no |
| 2 | 79 | 1 | 13 | 0 | no |
| 3 | 299 | 13 | 169 | 13 | no |
| 4 | 794 | 605 | 1631 | 267 | no |
| 5 | 1586 | 1586 | 6188 | 1820 | yes |

The **lower-bound certificate** contains a 605-by-794 matrix B spanning a
subspace of normal forms of degree at most 4, and a vector ell of length 794.
The checker verifies, independently of how B was generated:

1. Its columns are all degree-at-most-4 squarefree monomials, in the specified
   degree-compatible order; its rows are independent.
2. g belongs to rowspan(B).
3. For every lower-degree echelon basis row b and every variable z_i,
   N(z_i*b) belongs to rowspan(B).
4. B*ell=0 and ell applied to the coefficient vector of 1 equals 1.

By (3)--(5), conditions 2 and 3 imply N(V_(F,4)) is contained in rowspan(B).
Condition 4 therefore proves 1 is **not** in V_(F,4). This argument requires
only a containing closed subspace, not trust that B is the minimal one.
A separate ordinary-ring implementation, without chain-quotient elimination,
also obtains rank 1631 and exclusion of 1 at degree 4.

For the **upper bound**, the checker replays exact echelonization starting
from g and applying only (4), recording these ranks at degree 5:

\[
 1,13,79,299,472,732,940,1168,1349,1503,1573,1585,1586.
\]

Every row constructed in this process belongs to V_(F,5), by (3). The final
space has all 1586 standard pivots, including the constant pivot, proving
1 belongs to V_(F,5). This is a finite arithmetic derivation; solution counts
and ML predictions are not used to insert any ideal consequence.

For a unit ideal, the least degree at which 1 is derived is exactly LFD:
once 1 is available at degree D, every polynomial h can be multiplied from
it within degree max(D,deg h). If 1 is unavailable below D, the defining
LFD condition fails there for h=1. Thus

\[
                 \boxed{\operatorname{LFD}(F)=5>4
                         =\operatorname{FFD}_{\rm hom}(F).}     \tag{7}
\]

The table also proves affine first strict fall = 4: there are no strict falls
through degree 3, while at degree 4 the lower-degree dimension is 267>169.
Another strict fall occurs at 5 because 1820>1631. Hence both FFD conventions
equal 4 in this counterexample, and neither equals LFD.

## 5. A second, smaller counterexample

Over F_257 take E:y^2=x^3+105x+196, target (42,56), k=3, and the same
normalized two-summand construction. Here 3*42^2+105=0 mod 257, giving
homogeneous FFD=3. Exact ordinary-ring calculations give ranks

\[
 \dim V_{F,2}=7,\qquad\dim V_{F,3}=58,\qquad\dim V_{F,4}=210.
\]

The last number is the full dimension binomial(6+4,4). The constant 1 is
absent at degree 3 and present at 4. Thus exact LFD=4, and checking the
lower-degree intersections gives affine first strict fall=3. This example
has fewer variables and can be checked entirely in the ordinary ring.

The normalized F_257 system now also has a complete Lean certificate path:
`Prime257System.affine_first_three_lfd_four` proves affine first strict fall
3 and exact LFD 4. A checked 58-row containing space excludes 1 at degree 3,
a checked 177-node ordinary-ring DAG derives 1 at degree 4, and a separate
quadratic witness establishes the first affine strict fall. Lean also checks
the actual curve equation, nonsingularity, S3 normalization identity, and
input-list semantics. See [the certificate guide](IC_CHECKED_PRIME_GAP.md).
The ordinary homogeneous FFD=3 statement now has an independent Lean proof,
modulo the full polynomial Koszul submodule, using the same indexed input
list. It is not inferred from the affine strict-fall theorem. See
[the homogeneous proof guide](IC_HOMOGENEOUS_FFD_LEAN.md). Equivalence with the
truncated finite-field quotient convention remains a written argument.

## 6. Reproduction and remaining questions

From the repository root:

```sh
python -m pip install numpy sympy pytest
python scripts/verify_ic_fall_counterexamples.py --ordinary-check
python scripts/verify_ic_fall_counterexamples.py --check-only
python -m pytest -q tests/test_ic_fall_certificates.py tests/test_ic_fall_proofs.py
```

`experiments/ic-fall-counterexamples/counterexamples.json` records every
coefficient, exact rank, curve/target check, version, and source hash. The
two compressed `*-lower-bound.npz` files contain the exclusion matrices and
dual functionals. `--check-only` checks those saved certificates and replays
the upper derivations, rather than regenerating the certificate files.
The code uses exact integers modulo the stated prime; floating-point ranks
and numerical tolerances play no part.

The earlier general theorem LFD<=m*k+1 still holds. The exact normalized
m=2 FFD theorem also still holds; it is precisely the gap to LFD that these
instances establish. Raw versus normalized presentation dependence remains
a separate issue.

An exact generic normalized FFD formula for m>=3, a sharp normalized LFD
formula, a constant LFD bound independent of k, or an unbounded-LFD family
are **not** established by these two counterexamples. Nor does a finite
counterexample disprove an average-case claim that permits exceptional
instances. No cryptographic-size index-calculus speedup follows.
