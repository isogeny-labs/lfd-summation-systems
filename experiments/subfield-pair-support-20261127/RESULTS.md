# Exact subfield pair-support screen

## Frozen question

Does a subfield-coordinate factor base improve two-summand relation support
and complete single-target cost, compared with an equal-size coordinate
subspace on the **same** non-descending binary curve and target subgroup?
This screen uses exact point addition to verify the relations behind the
binary third summation polynomial. It benchmarks both an exhaustive pair
table and a direct subfield-quadratic S3 oracle. The latter is a specialized
two-summand construction; it does not measure a general S3/LRF solver.

## Fixture and validity

`replay.py` fixes `K=F_(2^16)` with modulus
`z^16+z^5+z^3+z^2+1` and the curve
`E: y^2+xy=x^3+z^11*x^2+z`. Its group order is `65986=2*32993`,
with prime `32993`. Since `j(E)=1/z` is outside `F_(2^8)`, this curve
cannot descend to that subfield even after a `K`-isomorphism.

The candidate x-space is the embedded `F_(2^8)`. The control x-space is
`<1,z,...,z^7>`. Both lie in the absolute trace-zero hyperplane and
their nonzero-x points occupy the same non-prime 2-coset. The candidate
has 128 nonzero x-values with two points over each; the control has 129,
and the last sorted x-class is dropped to match **128 unknown factor-base
logs and 256 signed points**. The product-space dimensions are 8 and 15,
respectively. The order-two point at x=0 is excluded from both bases.

Each table entry is an actual equality `P_i+P_j=Q` on `E(K)`.
For every nonzero `Q`, the replay checks the corresponding binary S3
identity

    (x_i*x_j+x_i*x_Q+x_j*x_Q)^2+x_i*x_j*x_Q+z = 0.

It projects the equality into the prime-order subgroup by doubling and
solves the resulting linear system for all base logs and the target log.
The recovered scalar is checked independently by scalar multiplication.

## Frobenius-eliminated S3 oracle

Here `q=256` and `K/F_q` is quadratic. For a nonzero query point `Q`, put
`t=x(Q)`, and for candidate abscissas `x1,x2 in F_q`, put
`s=x1+x2` and `p=x1*x2`. The binary S3 equation becomes

    p^2 + t^2*s^2 + t*p + b6 = 0.

The exact base enumeration proves every point with `x in F_q` lies outside
the odd-order subgroup, apart from the excluded order-two point. Hence
every nonzero subgroup query has `t not in F_q`. Write
`A=t+t^q`, `B=b6+b6^q`, `N=t*t^q`, and `C=B/A`. Taking the
`q`-Frobenius difference of S3 yields `p=A*s^2+C`. Substituting
`u=s^2` back into S3 gives the **quadratic over F_q**

    A^2*u^2 + N*u + D = 0,    D=C^2+t*C+b6 in F_q.

For each of its at most two roots, `s` is uniquely determined because
squaring is bijective in `F_q`; `p=A*u+C`. The two abscissas are the roots
of the second quadratic `X^2+s*X+p`. Each unordered abscissa pair gives
at most one signed representation of a fixed nonzero `Q`, because the
base excludes order-two summands. Thus this oracle returns at most two
signed relations per query. It precomputes the Artin-Schreier map
`v -> v^2+v` on `F_q` to solve the two quadratics, and verifies every
returned point relation. A separate exhaustive pair-table audit confirms
exact set equality on all 1,024 queries for each target.

This lowers the **univariate degree of this two-summand relation query**
to two. It is not a bound on the solving degree of other arities or a
generic binary S3 presentation.
Subfield-coordinate factor bases are established index-calculus objects
([Diem, Sections 1–2](https://msp.org/ant/2013/7-6/ant-v7-n6-p01-s.pdf));
this record makes no novelty claim for the construction or elimination.

## Results

Exact enumeration over all 32,993 subgroup elements gives:

| Equal-size base | Distinct two-summand sums, including zero | Nonzero query coverage | Maximum multiplicity away from zero |
|---|---:|---:|---:|
| Subfield | 24,605 | 24,604 / 32,992 = 74.58% | 2 |
| Monomial control | 20,959 | 20,958 / 32,992 = 63.52% | 6 |

Both bases have exactly 32,896 unordered signed pairs. The subfield
base therefore has **17.4% greater distinct-query coverage** at fixed
base size and pair count. This is a measured component benefit, not a
lower total solve cost.

For a single target, collecting every base logarithm is unnecessary.
Write the relation matrix as `[M_B | c]`, with `c` the target column.
The target log is unique exactly when
`rank([M_B | c]) > rank(M_B)` for the consistent system. The replay
finds the first such prefix and solves that underdetermined system;
it also retains the later first-full-rank prefix as a separate record.

The following are cold **single-target** seconds from the frozen
SageMath 10.10.rc0 replay on this host. Pair-table runs charge base
construction, all 32,896 pair additions, queries through target
determination, and elimination. Direct-oracle runs charge base
construction, its Artin-Schreier lookup, queries, and elimination; the
pair table used only to audit completeness is excluded. The scalar is
planted for verification; the algorithms receive only its public target
point. Query counts are prefixes of the same deterministic 1,024-query
stream for each candidate and control.

| Target scalar | Direct subfield queries / seconds | Subfield pair-table queries / seconds | Control pair-table queries / seconds | Cold BSGS seconds | Cold rho seconds |
|---:|---:|---:|---:|---:|---:|
| 3425 | 43 / 0.396 | 43 / 1.501 | 80 / 1.154 | 0.00136 | 0.00103 |
| 17123 | 82 / 0.487 | 82 / 1.551 | 56 / 1.424 | 0.00186 | 0.01637 |
| 29113 | 78 / 0.457 | 78 / 1.307 | 62 / 1.150 | 0.00232 | 0.00170 |

The subfield base raises the chance that one query yields a relation,
but it does not consistently reduce the number of queries that
determine the target. The direct quadratic oracle saves the pair-table
work, but its cold method still takes 0.40–0.49 seconds. Cold BSGS uses
182 baby steps and finishes in roughly two milliseconds. Even this
specialized oracle remains behind the matched generic methods on these
three targets. The result is a bounded negative cost screen for this
fixture, not a universal complexity bound. The specialization requires
a proper degree-eight subfield of the degree-16 field, so it does not
transfer directly to prime extension degree 131.

`cold_target_seconds` is optimistic: it uses the known first prefix
that determines the target and excludes the cost of detecting that rank
condition during live collection. The old full-base-log cost remains in
`cold_total_seconds`; using it as the single-target cost would overcharge
this method.
Rho timing varies with the deterministic walk seed. All three targets,
query-hit counts, row counts and rank gates remain in `results.json`;
`replay.py` regenerates the individual queries and relation rows.

## Lazy single-target walk

`lazy_walk.py` tests the cheaper cold path on the same curve, generator and
three targets. It prepares eight deterministic point increments
`da*G+db*T` and generates each query with one group addition. The S3
quadratic roots are tested against the curve, and an abscissa is lifted
only when a returned relation uses it. Every admitted signed equality
and its S3 identity are checked. The target stops at the first prefix
where its column is independent of the base columns; the replay also
verifies the recovered scalar by multiplication.

| Target scalar | Queries to target | Relation rows | Distinct base abscissas used | Optimistic cold seconds | Cold BSGS seconds | Cold rho seconds |
|---:|---:|---:|---:|---:|---:|---:|
| 3425 | 85 | 91 | 97 | 0.210 | 0.00136 | 0.00110 |
| 17123 | 58 | 60 | 78 | 0.136 | 0.00197 | 0.01638 |
| 29113 | 75 | 74 | 89 | 0.172 | 0.00227 | 0.00172 |

These are measured on one host in SageMath 10.10.rc0. The candidate
charges its increment setup, lazy lifts, walk queries and solving. It
still trails cold BSGS by roughly 69–154 times across these targets.
The timing remains optimistic because it identifies the first target
determination offline and omits the cost of checking that condition
during collection. The 1,024-query stream, failed queries and point
lifts are replayable from `lazy_walk.py`; the counts and target checks
are frozen in `lazy_results.json`. This screen does not establish a
scaling law or a win on a larger field.

## Replay

```bash
sage -python experiments/subfield-pair-support-20261127/replay.py \
  > /tmp/subfield-pair-support-replay.json
python3 experiments/subfield-pair-support-20261127/verify.py \
  /tmp/subfield-pair-support-replay.json
sage -python experiments/subfield-pair-support-20261127/lazy_walk.py \
  > /tmp/subfield-lazy-walk-replay.json
python3 experiments/subfield-pair-support-20261127/verify.py \
  /tmp/subfield-lazy-walk-replay.json lazy_results.json
```

The verifier compares every discrete field, support, query, rank, and
recovered-scalar record with `results.json`; wall-time fields are omitted
from equality because they vary by host and load. The committed result
SHA-256 is `4f55a14e595b7765c79cd7caaf3bb34edc03db1027b6264d8daa3aa725575ff2`.
The committed lazy-walk result SHA-256 is
`a7c7b9a4cf11bff54bfc408770b56f63d2dd5d6ca23207f2783084f34a3bd956`.
