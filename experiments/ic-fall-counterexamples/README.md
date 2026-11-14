# Exact counterexamples: first fall degree does not bound last fall degree

The [full mathematical argument](../../docs/IC_FALL_COUNTEREXAMPLES.md) fixes
the prime-field PKM chain encoding and explains the exact certificate checks.

| p | Curve coefficients (a,b) | Target | k | Homogeneous FFD | Affine first strict fall | Exact LFD |
|---:|:---|:---|---:|---:|---:|---:|
| 2689 | (1739,1574) | (797,1184) | 6 | 4 | 4 | 5 |
| 257 | (105,196) | (42,56) | 3 | 3 | 3 | 4 |

The first target has prime order 1307, on a nonsingular curve of order 2614,
and lies outside the coordinate factor base. The polynomial input is the
chains plus one normalized S3 equation. Its degree-4 closure excludes 1;
its degree-5 closure derives 1. For a unit ideal the first degree deriving
1 is exactly LFD. Both FFD definitions equal 4 in this example.

This disproves universal FFD=LFD, LFD<=FFD, and the normalized two-summand
bound LFD<=4. It does not prove an unbounded-LFD family or decide whether a
larger constant bound, or an average-case heuristic, holds.

## Reproduce

From the repository root:

```sh
python -m pip install -r experiments/ic-fall-counterexamples/requirements.lock.txt
python scripts/verify_ic_fall_counterexamples.py --ordinary-check
python scripts/verify_ic_fall_counterexamples.py --check-only
python -m pytest -q tests/test_ic_fall_certificates.py tests/test_ic_fall_proofs.py tests/test_ic_dreg.py tests/test_ic_degree_ml.py
```

`counterexamples.json` records explicit polynomial coefficients, exact closure
dimensions and rank histories, group checks, package versions, source hashes,
and certificate hashes. The two `*-lower-bound.npz` files contain normal-form
row-space bases and dual functionals that annihilate them but send 1 to 1.

The lower-bound checker verifies generator inclusion and stability under
every permitted variable multiplication. It does not trust how the saved
matrix was generated. The upper bound is checked by replaying exact linear
combinations and permitted multiplications, including newly fallen rows.
Every multiplication is charged at its degree before chain reduction; the
proved quotient-closure equivalence makes this an ordinary-ring V-space
computation. There is no fixed iteration cap or solution-oracle stopping rule.

The optional ordinary-ring check uses a separately implemented algorithm on
all ordinary monomials, including 1820 columns at degree 4 in the larger
example. It obtains the same rank 1631 and excludes 1. Both implementations
use exact finite-field arithmetic, not floating-point ranks.

The earlier ML and group-search results are not used as proof. These are
finite computer-assisted counterexamples. The F_257 affine strict-fall/LFD
claim also has the separate Lean verification described below; this does not
formalize the larger case or establish a cryptographic-size speedup.

A separate [Lean project](../../formal/README.md) proves the F_5 presentation
counterexample and now the actual normalized F_257 system's affine first
strict fall 3 and exact LFD 4. Its 58-row containing space and 177-node upper
DAG are checked by the kernel, with the actual curve/S3/input identities.
See [the checked-prime guide](../../docs/IC_CHECKED_PRIME_GAP.md).
An independent Lean proof now verifies ordinary Koszul FFD=3 on the F_257
system. The truncated-quotient equivalence and larger F_2689 gap remain
outside Lean.
