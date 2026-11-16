# Explicit exceptional-prime cutoff, checked in Lean

The source is `formal/ICFallDegrees/GrowthExplicitPrimeBound.lean`.
The default target imports it, and `Audit.lean` prints its nine headline
statements. The cutoff is the closed natural-number expression

    B_k = 2^(2^k) * 9^(4^k).

The curve coefficients are fixed at `(a,b,t)=(6,-7,1)`. The public results
refer to the actual normalized `PrimeM2General.inputs`.

## Statements and quantifiers

* `unit_prime`: for every k>=2, every prime p>B_k, and `2^k | p-1`, the
  chain ideal over F_p is the unit ideal. There is no supplied unit witness.
* `proper_prime_le_bound`: any admissible prime at which that ideal is
  proper satisfies p<=B_k. No converse is claimed.
* `prime_lfd_window`: for every k>=24 and every admissible prime p>B_k,
  the actual system has a last fall degree d with

      floor((k-10)/14) < d <= k+2.

* `exists_prime_window_above`: for every k>=24 and every external B there
  is a prime p>max(B,B_k), congruent to 1 mod `2^(k+1)`, with that window.
  The sufficient cutoff is explicit; no least prime or numerical
  growing-prime sequence is evaluated.

Admissibility in the first three statements is only `2^k | p-1`, so the
universal conclusion also covers the smaller progression used by the
factor base. The last statement keeps the stronger congruence of the
earlier existence theorem. The window assumes neither CurveModel nor a
unit certificate nor an unspecified good-prime cutoff.

## Algebraic proof

Let N=2^k. The integer matrix S has entries `S_ij=1` when j is the next
index after i in the cyclic order, and zero otherwise. Set

    A_k = (S+8I)^N - (S-I)^N,     D_k = det(A_k).

1. `eliminant_mod_two` evaluates the already proved polynomial identity
   `(X+8)^N-(X-1)^N=1` over F_2 at S. Therefore A_k reduces to I, D_k is
   odd and, in particular, nonzero. No complex nonvanishing or existential
   integral Bezout identity is required for this step.
2. Each column of S has absolute sum 1. The module proves the additive,
   multiplicative and power bounds for column absolute sums, so columns
   of A_k have sum at most `9^N+2^N`.
3. The determinant expansion is bounded by the sum over permutations of
   products of absolute entries. Enlarging to all index functions bounds
   this by the product of column absolute sums. Consequently

       0 < |D_k| <= (9^N+2^N)^N <= 2^N * 9^(N*N) = B_k.

4. A prime p>B_k cannot divide D_k, so A_k has nonzero determinant modulo p.
5. An N-th root x gives the nonzero vector `v=(1,x,...,x^(N-1))` and
   `S*v=x*v`. Thus `A_k*v=((x+8)^N-(x-1)^N)*v`. An S3 grid zero satisfies
   `(x-1)z=x+8` and z^N=1; it would give A_k*v=0, contradicting the
   nonzero determinant.
6. The finite-grid theorem and raw-to-normalized chain identity turn this
   exclusion into membership of 1 in the actual chain ideal. The module
   proves B_k dominates both characteristic thresholds needed by the
   existing native prime nonmembership and selector results, and composes
   them to obtain the least-LFD window.

The matrix is an explicit cyclic eliminant matrix. Its identification with
the usual univariate resultant is checked numerically in the regression;
the Lean proof does not require a general resultant API or assume that
identification.

The cutoff is deliberately conservative. It supplies a sufficient condition,
not a practical prime-search bound or a computation of large-depth closures.
The rational-complex family is unchanged; this prime-field statement is for
the fixed E6 coefficients. The non-torsion, genus-one case remains separate.

## Reproduce and inspect evidence

```sh
python3 scripts/verify_lean_fall_degrees.py
python3 scripts/lean_check_all.py --out /tmp/explicit-prime-bound-all-modules
python3 scripts/check_n8_explicit_prime_bound.py
python3 scripts/check_publication_artifacts.py
```

The pinned Lean version is 4.19.0 and the Mathlib revision is
`c44e0c8ee63ca166450922a373c7409c5d26b00b`. The local default build reuses
compiled dependencies. `build.log` is the standalone module elaboration;
`default-build.log` records the full default-target verifier. The full
axiom output and source hashes are in `experiments/ic-lean/`.

`regression.json` is recomputed by the script, with `--write` required to
replace it. It includes k=0..4 integer matrices and their determinants,
agreement with independently computed univariate resultants, mod-two
identity checks, column-sum bounds, 139 finite-field eigenvector checks,
and five deliberately wrong-constant negative controls.

It also checks a real unit instance above B_2 at
`p=29648323021629577`, congruent to 1 mod 8. SymPy primality testing is
deterministic in this p<2^64 range. That example checks the unit criterion
at k=2, not the growth conclusion whose range begins at k=24.
The existing Python CI loop runs the regression automatically. Linux Lean
CI builds all modules and audits every file that prints axioms.
