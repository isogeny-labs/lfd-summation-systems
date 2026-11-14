# Depth step: dim V_3 and dim V_4 against the ideal

`growth.json` has one row per (curve, target, depth k = 3, ..., 10) for five curves over p = 257, 7681 and
12289: dim V_3, dim V_4, dim(V_4 ∩ M_≤3), |M_≤4|, the number of scalar solutions and the dimension of the
degree-4 part of the ideal. Observation recorded with the table: dim V_3 = 2k + 1, and for k ≥ 6 dim V_4 does
not depend on the prime, curve or target; dim V_4 grows quadratically while |M_≤4| grows quartically, so V_4 is
the full degree-4 part of the ideal for k ≤ 5 and is deficient from k = 6 on.

Recompute (about 30 seconds; overwrites `growth.json` and `verification.json`, and only the timing field
changes):

    python scripts/run_ic_depth_step.py
