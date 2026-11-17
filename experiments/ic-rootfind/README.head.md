# Point decomposition by root finding instead of a Gröbner basis (binary curves)

**Question.** The relation search of index calculus asks for every `R = P_1 + ... + P_m` with `x(P_i)` in a
subspace `V` of `F_(2^n)` (dimension `k`). The usual route is the Weil descent of `S_(m+1)` and F4/F5. Is there a
faster route that finds the roots directly?

**Answer.** Yes, when the factor base is `V = <1, w, ..., w^(k-1)>` (or any `V` whose product space `V2 = span(V·V)`
is small). For such `V` the two-summand problem is one F_2-linear system in `k + 1` unknowns and one quadratic, so
`m` summands cost `|FB|^(m-2)` small linear systems. Enumeration costs `|FB|^(m-1)` curve additions. For a random
`V` with `k(k+1)/2 >= n` the method gives nothing; the table below shows the cliff.

Code: `csd/ic_rootfind.py` (mathematics, tables, Python reference), `csd/_rootfind.c` (kernels, `n <= 63`),
`tests/test_ic_rootfind.py`, driver `scripts/run_ic_rootfind.py`.

## The algorithm (LRF, "linearized root finding")

`P_1 + P_2 = T` iff `P_1, P_2, -T` lie on a line `y = λ(x + xT) + yT + xT`. Substituting into
`y^2 + xy = x^3 + a x^2 + b`:

    e1 = x_1 + x_2 = λ^2 + λ + a + xT          e2 = x_1 x_2 = xT λ^2 + ((yT + xT)^2 + b)/xT

Both are F_2-affine in `λ`. `x_1, x_2 ∈ V` forces `e1 ∈ V` and `e2 ∈ V2`. Writing `e1 = Σ t_i v_i` and
`λ = Q(e1 + a + xT) + ε`, with `Q` a fixed linear right inverse of `z ↦ z^2 + z`, the condition `e2 ∈ V2` plus one
trace condition is a linear system in `(t, ε)`: `k + 1` unknowns, `n - dim V2 + 1` equations. Each solution gives
`x^2 + e1 x + e2`; it is kept when its roots lie in `V`. The points come with it (`y_i = λ x_i + μ`), so only
rational decompositions with the right signs are produced.

For `m` summands, enumerate `m - 2` factor-base points and solve the two-summand problem for
`T = R - P_3 - ... - P_m`.

Per target the kernel does `k + 1` field multiplications, a `(k + 1)`-column elimination on 64-bit words, and about
six more multiplications per candidate; inversions are shared across blocks of 256 targets and 256 candidates
(Montgomery's trick). The enumeration baseline shares its inversions the same way.

The linear structure of `S_3` in the symmetric functions is not new (it is behind the symmetrized systems of
Faugère–Perret–Petit–Renault and Huang–Petit–Shinohara–Takagi, and Semaev's splitting into `S_3`'s). What is here
is the solver built on it and its measurement against the baselines.

## Checks

- `tests/test_ic_rootfind.py`: on 480 instances (`n` from 7 to 17, `m` = 2, 3, 4, standard and random `V`, random
  and Koblitz curves; 87 of them have decompositions) LRF, the C enumeration, the C meet in the middle and a
  plain Python enumeration on `csd.ic_dreg`'s group law return the same sets. The C field arithmetic and factor base
  are compared with `csd.ic_dreg`'s. The targets with `x(T) = 0` (the one case the line cannot take) are tested
  separately.
- Every row of every record below: all methods run on the instance return the same set, and every decomposition is
  re-added to `R` with `csd.ic_dreg`'s group law before the row is written.

## Results

Single thread, wall clock, minimum of three runs for the C kernels, median over three seeded instances per size
(random `a ∈ {0,1}`, random `b`, random `R`). Machine: 4-core Intel Xeon (PCLMULQDQ), gcc 13.3 `-O3 -march=native`.
Set-up of the tables is excluded (it is below 0.2 s everywhere).

