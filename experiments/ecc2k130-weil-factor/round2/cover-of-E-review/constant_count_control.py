#!/usr/bin/env python3
"""Positive control for the constant-group count rule (m = 1): over F_Q with ell | Q-1 the cover
z^ell = x of E is a CONSTANT Z/ell cover; the rule predicts #C(F_Q) = ell*a + 2 (2 rational branch
points), i.e. #C(F_Q) = 2 mod ell, and a nonzero Prym trace is allowed. Twisted bases (ell nmid Q-1)
predict #C(F_Q) = #E(F_Q). Cheap: only F_Q itself is counted."""
from kummer_control import count
rows = []
for ell, b in [(3, 2), (5, 4), (7, 3), (3, 4), (13, 12), (7, 6), (5, 1), (5, 2), (7, 2), (11, 5), (131, 13)]:
    Q = 2 ** b
    const = (Q - 1) % ell == 0
    nC, nE = count(ell, b)
    ok = ((nC - 2) % ell == 0) if const else (nC == nE)
    rows.append((ell, Q, const, nC, nE, ok))
    print(f"ell={ell:3d} Q=2^{b:<2d} constant={const!s:5}  #C(F_Q)={nC:5d}  #E(F_Q)={nE:5d}  "
          f"{'#C=2 mod ell' if const else '#C=#E'}: {ok}   trace_on_Prym(F_Q) = {Q+1-nC - (Q+1-nE)}")
assert all(r[-1] for r in rows)
print("ALL CONSTANT/TWISTED COUNT CONTROLS PASSED")
