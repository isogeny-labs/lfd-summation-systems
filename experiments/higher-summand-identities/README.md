# Checks for growth with every number of summands (paper, Section `sec:growthm`)

Added 2026-11-23 (PR #49). Pure Python, no dependencies; the resultant evaluator
is written from scratch and does not import the PR #45 engine.

| script | what it checks | result (`*-20261123.txt`) |
|---|---|---|
| `check_identities.py 6 8` | symmetry of the recursive `S_n` (n = 4..7); the order-two identities (A) `S_n(U,1,1) = c^(2^(n-3)) S_(n-2)(U)^4` and (B) `S_4(x,y,1,t) = (t-1)^4 S_3(x,y,u)^2` on `E_a`; both parity formulas for m = 2..8 (Remark `rem:ordertwoslice`); a mutation control (wrong scalar exponent) and the specialize-first trap (must lose `(x-y)^4`); a small unit-instance sanity check on `y^2 = x^3 - 1` | all pass: SYM 96, A 432, B 72, EVEN 288, ODD 216, mutant detected 360/360, trap detected 24/24 |
| `check_slice_membership.py` | the new step of Corollary `cor:growthm`: `S_(m+1)(x(P), x(R'-P), q_3..q_m, t) = 0` with `R' = R - Q_3 - ... - Q_m`, m = 3..7, four primes, random curves and rational points; cancelling pairs, one unpaired summand, arbitrary `Q_i` | all pass: PAIR 96, ODD 144, GEN 240; control `R' + W` (random `W`) nonzero 160/160 |
| `thresholds.py` | `k_D`, and the first depth at which the slope-1/14 bound exceeds the constant floor `m(m-1)+1` | 90, 167, 270, 399, 554, 736 for m = 3..8 |

These are exact evaluations at random points (Schwartz-Zippel style evidence for
polynomial identities), not proofs and not closure computations. A first version
of the membership control used `R' + 2Q_3`; that is the target of another sign
pattern, where `S_(m+1)` also vanishes, so it was replaced by a random offset.
