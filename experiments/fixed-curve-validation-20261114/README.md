# Fixed-curve validation and two follow-ups

Status: exact computations and written arguments, **not new Lean theorems**.
Source baseline: companion commit `49585b5`. That baseline already proves the
complex growth window for `|a| > 5` and the Fermat-prime obstruction in Lean.
The earlier depth-2 determinant is also already in
`scripts/check_n8_prime_specialization.py`. Those results are not claimed as new here.

## Reproduce

From the repository root, with SymPy installed:

```sh
python3 scripts/check_fixed_curve_obstructions.py
```

`--write` regenerates `results.json`. Without it the script recomputes and
compares the entire record, including its own SHA-256. It imports no project
closure engine, Semaev generator, or subgroup builder. No timings are reported.

The completed checks are:

* 4,759 cells: every prime `p < 100000` with `p = 1 mod 2^(k+1)`, at
  `k = 2..10`. There are 4,636 unit cells and 123 proper cells.
* Direct raw-S3 evaluation on 16,368 ordered grid pairs agrees with the
  one-variable graph test. Changing the graph numerator from `x+8` to
  `x+7` disagrees on 12 of the 15 controls.
* At `k = 2,3,4`, an integer resultant predicts proper/unit status at every
  scanned prime. Below 20,000 a separate polynomial gcd gives the exact
  ordered solution count. The k=2 resultant equals an independently built
  16-by-16 integral multiplication determinant.
* Computed integral inverses satisfy `q*v=d` at k=2,3,4. The record includes
  every coefficient. Cyclic convolution separately checks `q^2*v^2=d^2`.
  These are unit certificates away from the scalar's prime divisors, not
  last-fall-degree certificates. The replay verifies `d | R_N` and
  `R_N | d^(N^2)`, establishing the same prime support without factoring R_N.
* 260 exact rational-parameter gcd checks: 65 distinct rational parameters,
  `k = 2..5`. A root occurs exactly at the parameters `-3,-1,1` in these tests.

The resultant factorizations are complete at k=2,3. At k=4 only trial
division through 20,000 was performed; the remaining cofactor is explicitly
labelled unfactored. This avoids an unbounded factorization task.

| k | Digits in resultant | Digits in inverse denominator d |
|---|---|---|
| 2 | 15 | 6 |
| 3 | 58 | 24 |
| 4 | 232 | 91 |

The first two denominators are 283815 and 159329046096905711288535.
The third is stored in `results.json`; it is computed exactly, while its
complete factorization is not needed. The unique rational inverse is solved
by exact linear algebra and its common integer denominator is reduced by
the gcd of all numerators and the denominator.

All these finite depths are below the linear-growth theorem's k>=24 range.
They test its algebraic ingredients and obstructions, not the large-depth
closure conclusion. These bounded prime counts imply no prime-density law.

## 1. Explicit exceptional-prime criterion and threshold

**Written argument.** Let N=2^k, k>=2, and

    q(x,z) = (x-1)z-x-8,
    f(X) = X^N-1,
    g(X) = (X+8)^N-(X-1)^N,
    R_N = Res_X(f,g) in Z.

On the complex unit-circle grid, q=0 would imply
`(x-1)(z-1)=9`, whereas its left side has absolute value at most 4.
Thus R_N is nonzero: for any complex N-th root x, x=1 gives g(1)=9^N,
and otherwise g(x)=0 would give the forbidden N-th root
z=(x+8)/(x-1).

Since f is monic, the product formula for the resultant gives

    0 < |R_N| <= (9^N+2^N)^N <= 2^N * 9^(N^2).

For any admissible prime p>2N with N dividing p-1, f splits into distinct
linear factors over F_p, and x=1 is not a root of g (p != 3). Consequently

    p divides R_N  <=>  the chain system has a geometric grid zero
                   <=>  its ideal is proper.

The equality between grid and geometric zeros uses that every chain solution
is determined by its initial N-th roots and their successive squares. Raw
and normalized S3 differ by chain-ideal terms; S3=q^2 has the same zeros as q.
The chain quotient is a product of copies of F_p, so a grid function with
no zero is invertible, giving the unit ideal without a radical-ideal
assumption on the full S3 system.

In particular every admissible prime

    p > B_k := 2^(2^k) * 9^(4^k)

gives a unit ideal. This is an explicit but very loose threshold, replacing
the existential unit threshold in `GrowthFixedE6NativeUnit`. Together with
the existing large-prime nonmembership and selector bound it gives a written
route to the same LFD window with an explicit sufficient prime threshold.
It is exponential in 4^k and provides no efficient prime-search or solve-cost
bound. The all-depth resultant argument and threshold are not formalized here.

The complete admissible bad-prime lists for the first two depths are:

| k | All admissible bad primes |
|---|---|
| 2 | 17 |
| 3 | 17, 97, 257, 2417, 4481 |

For example k=3,p=97 has the ordered zero (75,47), separate from the
Fermat-factor construction. This is consistent with the theorem's stated
prime dependence.

### Diagonal versus off-diagonal obstructions

There is a useful written factorization, checked exactly at k=2,3,4:

    R_N = (4^N-1)(2^N-1) * S_N^2,   S_N in Z.

Indeed R_N is the determinant of multiplication by q in the grid algebra,
and over C it is the product of q(x,z) over all ordered pairs of N-th roots.
The diagonal part is the product of
`q(x,x)=(x-4)(x+2)`, giving `(4^N-1)(2^N-1)` since N is even. Each unordered
off-diagonal pair contributes the square of q(x,z), by symmetry of q.
The product S_N over unordered pairs is Galois invariant and an algebraic
integer, hence an integer. This is an all-depth written argument, not a
formalized determinant decomposition.

Over each scanned prime, the diagonal roots can only be 4 and -2; they
are distinct because p>8. The replay counts their membership in the subgroup
separately and verifies that the number of off-diagonal ordered solutions
is even. Among the 123 proper cells, 107 have no diagonal solution and 117
have an off-diagonal solution. The record stores both counts. Fermat-factor obstructions are
diagonal; the S_N factor is where further off-diagonal obstructions can
occur. The two integer factors need not have disjoint prime support.

## 2. Rational parameters: an exact complex unit criterion

**Written argument.** For rational a, target t=1, and
`E_a: y^2=x^3+a*x-a-1`, the raw S3 polynomial equals

    ((x-1)z-x-2-a)^2.

For every k>=2, its complex two-power grid has a zero **if and only if**
`a` is one of `-3,-1,1`.

Proof. For a=-3 the pair x=z=1 is a zero. Otherwise x=1 is impossible.
Let x,z be roots of unity satisfying the displayed linear factor. As a is
real and |x|=|z|=1, taking squared absolute values in
`(x-1)z=x+2+a` gives

    (a+3)(a+1+x+x^-1)=0.

Therefore u=x+x^-1=-a-1 is rational. It is an algebraic integer, since
x and x^-1 are roots of unity; a rational algebraic integer is an integer.
Also u is real with |u|<=2. If u=-1 then x^2+x+1=0, so x has order 3.
If u=1 then x^2-x+1=0, so x has order 6. Neither order divides 2^k.
The remaining possibilities u=-2,0,2 give a=1,-1,-3. Conversely the pairs
(-1,-1), (i,-i), (1,1) witness those parameters, respectively, on every
grid with k>=2. This proves the criterion.

The singular parameters of the elliptic model are -3 and -3/4, because
`4a^3+27(a+1)^2=(a+3)^2(4a+3)`. Thus, among nonsingular rational members,
only a=-1 and a=1 obstruct the complex unit ideal.

This can extend the existing complex growth family beyond |a|>5. The
current `GrowthTwoTorsionFamily.model` and `actual_spread_bound` already
need only a+3 != 0. Once this new unit criterion is formalized, the existing
growth and selector arguments can be composed for every rational a outside
{-3,-3/4,-1,1}. The gcd checks supply finite regressions, not that formal proof.
At a=-3 the elimination gcd degree is N, not the full two-dimensional
solution count; only the presence of a root is used in the parameter tests.

## Follow-up order and falsifiers

1. Formalize the rational trace-of-root-of-unity argument, then transport
   the no-zero result through the existing chain-grid unit theorem. Check
   a=-1 and a=1 as negative controls and a=0,2,-2 as members newly covered
   beyond the norm-five condition. Do not retain a norm-five assumption in
   the final generalized statement.
2. Formalize the resultant elimination and size bound. Its useful output
   is a concrete B_k; reusing the existential Bezout scalar alone does not
   settle this follow-up. Compare the specialization at k=2,3 with the
   complete prime lists above.
3. Seek a smaller uniform threshold using the denominator of q^-1 in the
   finite integral quotient instead of its full determinant. The record
   establishes a reduction at three depths with identical prime support.
   Every column of the inverse multiplication matrix is a translate of
   q^-1, so its least common denominator d divides det(M), while the integer
   matrix d*M^-1 gives det(M) dividing d^(N^2). This explains the prime-support
   equality. A useful all-depth size bound or recurrence for d remains open;
   three examples do not establish one.
4. Keep the non-torsion, genus-one target separate. None of these arguments
   discharges its function-field hypothesis or establishes the conjectured
   slope 1/4. The present rational-model extension remains a two-torsion result.
