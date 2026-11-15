# N8 fixed-E6: one admissible bad characteristic at every depth

**Evidence:** elementary exact proof below; **kernel** for the statement about
the polynomial system itself, in `formal/ICFallDegrees/GrowthFermatDiagonalObstruction.lean`
(Lean 4.19.0, Mathlib v4.19.0; axioms `propext`, `Classical.choice`, `Quot.sound`):

* `e6_diagonal_grid_root`: if `2^(2^k) = −1` in a commutative ring then `4^(2^k) = 1`
  and `h(4,4) = 0` (the statement of ml-cryptanalysis `0712c0d5`);
* `e6_not_unit_of_fermat`: over a field with `2^(2^k) = −1`, `1` is not in the ideal of
  `PrimeM2General.inputs k _ 6 (−7) 1` (evaluation at the point whose chain variable of
  level `j` is `4^(2^j)`, for both summands);
* `e6_not_unit_of_fermat_factor`: for a prime `p ∣ 2^(2^k)+1` and `k ≥ 2`:
  `p ≡ 1 (mod 2^(k+1))`, `p > 2^(k+1)`, and the system over `F_p` is a proper ideal;
* `exists_bad_admissible_prime`: such a prime exists at every depth `k ≥ 2`.

The last three were added in the companion repository on 2026-11-14; upstream has the
first. Exact integer regression: `scripts/check_n8_fermat_diagonal_obstruction.py`
(`k = 2, …, 5`, including the composite `F_5`).

For each integer k ≥ 2 let N=2^k and F_k=2^N+1. Choose **any prime divisor**
p of F_k. In F_p we have 2^N=−1, hence 2^(2N)=1 but 2^N≠1. Because the
multiplicative order of 2 divides 2N and does not divide N, it is exactly
2N=2^(k+1). Therefore 2^(k+1)|(p−1), so p meets the N8 admissibility
condition. The multiplicative order of 4 is N, giving 4^N=1.

For the fixed E6 curve y²=x³+6x−7 the normalized eliminant is

    h(x,z)=(x−1)z−x−8.

There is the *integer identity* h(4,4)=3·4−4−8=0. Thus the torus system
x^N=z^N=1, h(x,z)=0 has an F_p-point (4,4), and its chain extensions
obtained by successive squaring give a common root of the corresponding
E6 chain equations. Consequently the prime-field E6 chain ideal is
**not** the unit ideal at this (p,k).

The E6 discriminant's nonzero factor is 2187=3^7. For k≥2, N is even,
so F_k=2^N+1 ≡ 2 (mod 3), and p≠3. Hence E6 remains nonsingular in
these counterexamples. The conclusion applies **at every depth k≥2**,
not just the known p=17,k=2 and p=257,k=3 examples.

Exact sample factors:
- k=2: F_2=17;
- k=3: F_3=257;
- k=4: F_4=65537;
- k=5: F_5=641·6700417.

This is **not an obstruction to the almost-all-primes theorem**: for each
fixed k the divisors of F_k constitute a finite exceptional set. To
build a growing sequence of good primes exceeding the separate rank and
Bézout exceptional bound B(k), choose factors of F_m for sufficiently
large m relative to k and B(k), rather than automatically choosing
factors of F_k. The existing #249 prime-supply lemma gives the
appropriate unbounded choice. No claim about all elliptic curves,
uniform prime-field LFD growth, or ECDLP speedups follows.
