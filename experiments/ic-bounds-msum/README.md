# Three to eight summands: degree of NF(S_{m+1}) and the last fall degree it forces

Script: `scripts/ic_bounds_msum.py` (numpy only). Records: `results.jsonl` (one line per instance),
`summary.json`, `table.md`. Evidence labels: **kernel** = Lean-checked statement; **exact** = computation in
exact arithmetic modulo p, command named here.

## What is computed

For the direct presentation G_{m,k} + {g}, g = NF(S_{m+1}(x_1, ..., x_m, t)), over F_p with 2^k | p - 1:

- **kernel** upper bound: LFD <= mk + 1 (`ChainQuotient.reduced_inputs_lfd_bound`);
- **kernel** lower bound: LFD >= deg g + 1 when h = NF(x_1 g) has degree at most deg g and is not a
  constant multiple of g (`CarryLower.chain_carry_lower`, with a = x_1).

So LFD = mk + 1 exactly when deg g = mk and the two hypotheses on h hold. The script computes g, h, their
degrees and the proportionality test (**exact**). It does not expand S_{m+1}: modulo the chains a polynomial is
a function on H^m, so g is the inverse discrete Fourier transform of the values of S_{m+1} on H^m, and the values
come from Semaev's recursion S_n = Res_Y(S_{n-1}(.., Y), S_3(x_{n-1}, X, Y)) carried out with numbers
(interpolation in X, resultant with a quadratic as a norm).

    python3 scripts/ic_bounds_msum.py selftest                 # must print SELFTEST PASSED
    python3 scripts/ic_bounds_msum.py run --max-log2-points 24 # writes this directory (about 5 minutes)

## Self-test (**exact**)

- S_3 and S_4 on the grid equal `csd.ic_dreg.s3_prime_value` and `csd.ic_bounds_m3.s4_value`; NF(S_4) equals
  `csd.ic_bounds_m3.direct_generator` coefficient by coefficient (three instances).
- The top coefficient of NF(S_4) at k = 2 equals -8(5t^3 + 3at + 6b) (three instances).
- For m = 2, ..., 8: S_{m+1}(x(P_1), ..., x(P_m), x(P_1 + ... + P_m)) = 0 for random rational points P_i
  (group law), and the value is invariant under a permutation of the x_i.
- Grid values equal a separate scalar evaluation at 40 random grid points for (m, k) = (4, 2), (5, 2), (4, 3).
  The scalar evaluation shares the resultant routine with the grid code; the planted-sum test does not.
- All 344 certified three-summand records of `experiments/ic-bounds-m3/` with k <= 3 (direct encoding) are
  consistent: same generator degree, same number of solutions, exact LFD inside [lower, upper], and equal to
  mk + 1 on the 205 records where this script claims an exact value.

## Results (**exact**; 1,132 instances; primes 97, 193, 257, 7681, 12289, 65537; curves y^2 = x^3 + x,
y^2 = x^3 + 1 and two random curves per prime; two random rational targets per curve)

See `table.md`. For m >= 3 and 2 <= k <= m - 1:

| m | S_{m+1} | k = 2 | k = 3 | k = 4 |
|---|---|---|---|---|
| 3 | S_4 | LFD = 7 on 48 of 48 | | |
| 4 | S_5 | 9 on 48 of 48 | 13 on 47 of 48 | |
| 5 | S_6 | 11 on 48 of 48 | 16 on 48 of 48 | 21 on 8 of 8 |
| 6 | S_7 | 13 on 48 of 48 | 19 on 8 of 8 | 25 on 2 of 2 |
| 7 | S_8 | 15 on 47 of 48 | 22 on 8 of 8 | not computed (2^28 points) |
| 8 | S_9 | 17 on 47 of 48 | 25 on 2 of 2 | not computed |

- 409 of 412 instances have deg g = mk and satisfy the hypotheses on h: LFD = mk + 1 exactly there.
- The three exceptions are on y^2 = x^3 + x ((m, k, p, t) = (4, 3, 193, 133), (7, 2, 97, 46), (8, 2, 193, 154)):
  deg g = mk - 1 and deg h = mk, the ideal is the unit ideal (no solution), and this script decides nothing for them.
- For k >= m: deg g = m(m - 1) on every instance, so LFD >= m(m - 1) + 1: 7 for m = 3 at k = 3, 4 (96 instances,
  as known), 13 for m = 4 at k = 4, 5 (56 instances; upper bounds 17 and 21). No exact value there.
- m = 2 is in the records for completeness: the hypotheses on h fail (deg h > deg g), as they must, since the
  two-summand values are above deg g + 1.
- k = 1 (H = {1, -1}) is in the records and not used in the paper.

## Limits

- Instances are few at the largest sizes (8 at 2^20 and 2^21 points, 2 at 2^24).
- Depths with more than 2^24 grid points were not run: (m, k) = (6, 5), (7, 4), (8, 4) and beyond.
- The exact values are for the direct presentation only. Nothing here concerns split encodings.
- The lower bound for m >= 4 uses the Lean statement with hypotheses established by this computation, not by
  the kernel; for m = 3 and 5t^3 + 3at + 6b != 0 the kernel proves them (`M3CarryLower.m3_lfd_ge_seven`).
