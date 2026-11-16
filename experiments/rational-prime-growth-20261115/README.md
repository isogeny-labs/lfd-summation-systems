# Rational-parameter prime-field growth with an explicit cutoff

The parameter family is `y^2=x^3+a*x-(1+a)`, with target abscissa 1.
For every rational `a` outside `{-3,-3/4,-1,1}`, write `a=m/q`, where
`m` is an integer and `q>0` is an integer. Lowest terms are not required.
Put

    N = 2^k
    H = 3q + |m+2q|
    C(k,m,q) = 2^N * H^(N*N) + 2N + q + |m+3q| + |4m+3q|.

For every k>=24 and every prime p>C(k,m,q) with `2^k | p-1`, the
actual normalized chain system over F_p, at parameter `m/q` reduced
modulo p, has a least last fall degree d with

    floor((k-10)/14) < d <= k+2.

The theorem proves the denominator is invertible; it does not assume a
field homomorphism Q -> F_p. The canonical numerator and positive denominator
give a specific cutoff for every allowed rational parameter.

## Lean statements

`formal/ICFallDegrees/GrowthRationalPrimeNative.lean` generalizes the native
finite-box proof to the section `(X,(X+c)/(X-1))`, for arbitrary c. Common
denominator clearing, degree bounds, rational Wronskian nonvanishing,
cofactor ratios and pole-gap estimates are proved in the native field.
The spread result needs `c+1 != 0`. With c=a+2, this gives nonmembership
of 1 in the actual closure whenever `a+3 != 0` and `p>2^(k+1)`.
It assumes no transfer of circuits from characteristic zero.

`formal/ICFallDegrees/GrowthRationalPrimeCutoff.lean` proves:

- `matrix_det_bound`: the denominator-cleared integer eliminant matrix has
  determinant with absolute value at most `2^N * H^(N*N)`.
- `value_ne_zero_complex` and `matrix_det_ne_zero`: its eigenvalues cannot
  vanish, by the proved rational root-grid classification. The actual
  Vandermonde matrix of complex Nth roots diagonalizes it and is invertible.
- `matrix_det_mod_prime`: the nonzero determinant survives reduction for
  every prime above the displayed determinant bound.
- `unit_prime`: nonzero reduced determinant excludes S3 grid zeros and
  gives the unit ideal for the actual normalized inputs.
- `denominator_and_parameter_guards`: the cutoff proves q, m+3q and
  4m+3q remain nonzero modulo p, so division is valid and the reduced curve
  is nonsingular. The integer factors are nonzero by the rational exclusions.
- `prime_lfd_window`: the unit, native nonmembership and selector results
  give the exact least-LFD window for any representation a=m/q with q>0.
- `rational_prime_growth`: the public version uses a.num and a.den.

The proofs use neither a supplied unit certificate nor a CurveModel nor an
unspecified exceptional-prime integer. The cutoff is intentionally
conservative. The sharper fixed-E6 bound `9^(4^k)` is retained; this generic
family cutoff is not asserted to improve that specialized result.

## Finite regression

`scripts/check_n8_rational_prime_growth.py` recomputes `regression.json`:

- 12 parameters at depths 2 and 3, including all excluded parameters and
  the singular parameter -3/4 as polynomial controls;
- exact integer matrix determinants and independent univariate resultants;
- 4,720 raw S3 checks against the denominator-cleared linear factor;
- agreement of determinant divisibility with root-grid solutions;
- two rejected denominator cases and wrong-constant negative controls;
- seven above-cutoff unit instances at depth 2, at primes below 2^64
  (where SymPy primality testing is deterministic).

The depth-2 examples test the unit and arithmetic components. They do not
numerically evaluate a depth-24 closure. The all-depth growth statement is
symbolic Lean evidence. The finite replay uses the project's raw S3
coefficient routine, while the determinant/resultant computations use SymPy.
The singular controls are polynomial systems, not nonsingular curves.

The result concerns target order two and rational function fields. It does
not formalize the general genus-one growth theorem or resolve the numerical
LFD lower bound at the proper-ideal exceptional primes below the cutoff.
It concerns admissible (p,k) pairs, with p allowed to grow with k.

## Reproduce

    python3 scripts/check_n8_rational_prime_growth.py
    python3 scripts/verify_lean_fall_degrees.py
    python3 scripts/lean_check_all.py --out /tmp/rational-prime-all
    python3 scripts/build_ic_curve_paper.py
    python3 scripts/check_publication_artifacts.py

The existing CI script glob runs the finite replay. Linux CI builds and audits
all Lean modules using Lean 4.19.0 and pinned Mathlib
`c44e0c8ee63ca166450922a373c7409c5d26b00b`. The default audit and source hashes
are retained in `experiments/ic-lean/`.
