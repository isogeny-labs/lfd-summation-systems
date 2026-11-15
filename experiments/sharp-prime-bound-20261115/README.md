# Sharper prime-field cutoff and follow-up evidence

`GrowthSharpPrimeBound.lean` improves the sufficient exceptional-prime cutoff
for the actual normalized E6 chain inputs from

    2^(2^k) * 9^(4^k)

to

    9^(4^k).

The improvement is strict at every depth, by the factor `2^(2^k)`.
It does not change the growth window or the admissibility condition.
For every k>=24 and every prime p above the new cutoff with `2^k | p-1`,
the actual least last fall degree d satisfies

    floor((k-10)/14) < d <= k+2.

## What removes the factor

Let N=2^k and let S be the cyclic shift matrix. For
`U=S+8I`, `V=S-I`, entrywise `|V| <= U`; multiplication and induction give
entrywise `|V^N| <= U^N`. Thus `A=U^N-V^N` has nonnegative entries.
Every column of U sums to 9 and every column of V sums to zero, so every
column of A sums to exactly `9^N`. The determinant bound by the product of
absolute column sums gives `|det A| <= 9^(N*N)` without the previous
triangle-inequality loss. The earlier mod-two identity still proves
the integer determinant is nonzero.

The new module also exposes `unit_of_not_dvd` and `prime_lfd_window_of_det`.
A prime that does not divide the explicit integer determinant has the unit
ideal, even below the size cutoff. The growth window additionally requires
the displayed characteristic bounds and k>=24. These are sufficient tests;
the Lean module does not claim a converse to the determinant criterion.

## Rational-parameter probe

The regression clears denominators for a=m/q and computes

    Res(X^N-1, (q*X+m+2*q)^N - (q*(X-1))^N).

It checks 12 rational parameters at depths 2 and 3, with finite-field grid
agreement whenever the denominator is invertible. The determinant vanishes
in these cases exactly at a=-3,-1,1, agreeing with the complex classification.
The singular a=-3/4 has a nonzero eliminant, illustrating that unit exclusion
and curve nonsingularity are distinct requirements.

These finite checks are not a formal all-rational prime-field growth theorem.
That still needs denominator-sensitive nonvanishing/size estimates and a
generalization of the native characteristic-p nonmembership proof, whose
current composition is specialized to E6. No such extension is claimed here.

## Proper-ideal exceptional primes

The exact closure replay gives:

| k | p | solutions | exact LFD |
|---|---|-----------|-----------|
| 2 | 17 | 1 | 4 |
| 3 | 257 | 1 | 4 |
| 4 | 65537 | 1 | 4 |
| 5 | 641 | 3 | 4 |
| 5 | 6700417 | 1 | 4 |
| 6 | 274177 | 1 | 5 |

Every listed prime divides `2^(2^k)+1`. The profiles check closed containing
spaces and the evaluation-rank/codimension criterion for an exact local
upper bound. This uses the project's closure engine, not an independent
implementation or Lean replay. These depths lie below the k>=24 growth
range; they neither establish nor refute that asymptotic lower bound.
The depth-6 case rules out extrapolating a constant exact value 4 from the
smaller cases. The general exceptional-prime LFD problem remains open.

## Reproduce

    python3 scripts/check_n8_sharp_prime_bound.py
    python3 scripts/verify_lean_fall_degrees.py
    python3 scripts/lean_check_all.py --out /tmp/sharp-prime-all
    python3 scripts/check_publication_artifacts.py

The deterministic JSON record is recomputed on every replay. `--write`
explicitly replaces it. The existing CI script glob includes this replay.
