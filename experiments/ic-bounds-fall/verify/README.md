# Verifier outputs for item 5 (fall-space structure)

These are the outputs of `scripts/ic_bounds_fall_verify.py` and
`scripts/ic_bounds_fall_verify_npz.py`. The engine is independent and does not use
`csd/ic_bounds_fall.py`. The commands, the per-claim verdicts and the findings are in
`../VERIFY.md`.

| file | contents |
|---|---|
| `closures.json` | Primary instance, D=4, k=6..12: dimensions, profiles, rank histories, low birth counts, closed check, whether 1 is in V_4. |
| `free-p40961.json` | Free extensions from base 7 (k=8..10) and from base 6 (k=7..10). Certificates, comparison with V_4(k), the stable low part, top-level ranks. |
| `k11/free-k11.txt` | Log excerpt of the k=11 free-extension check. |
| `k12/free-p40961.json` | The k=12 free-extension check. |
| `levels.json` | Möbius and graded level counts of F_4 at k=6,7,8, and the fall-space structure at k=8,9,10. |
| `w.json` | Symbolic 7x9 kernel and the check of w on the instance. |
| `cross.json` | p12289r3 and p12289j1728 at k=6..9: pivot sets of V_4(7) and the free-extension checks. |
| `mult.json` | Unit instance p=1073153: multiplier filtration and kappa_4, kappa_3 at k=6,7. |
| `d5-k{6,7,8}.json` | Degree-5 closures and the local-to-global certificate. |
| `npz.json` | The producer's exported bases compared with the recomputed bases, plus a check of the producer's file hashes. |
