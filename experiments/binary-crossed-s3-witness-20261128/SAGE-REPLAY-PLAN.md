# Independent Sage replay

Filed 2026-11-28 after `result.json`, before the Sage replay. This is a
post-run validation of the frozen binary screen, not a new selection rule.

Use Sage's finite field and elliptic-curve implementations, with the same
eight explicit `(n, modulus, a, b)` cells. Reconstruct polynomial-basis field
elements from their bit vectors; do not import `run.py` or use its group law.
For each cell independently enumerate Sage rational points, compare group
order and attained x count, compare direct pair-sum labels with roots of the
binary S3 formula plus O, and repeat every quadruple-target direct-versus-tree
comparison. Match the exact counts and both zero-mismatch results in the
committed `result.json`. Save all mismatches before the final assertion. The
replay output binds the source and original result by SHA-256.

This checks the instrument and finite cells. It does not prove a general
algebraic theorem or measure a solver.
