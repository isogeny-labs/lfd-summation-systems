# Verified Lean build and axiom audit

> **Note (2026-11-12).** Links to `docs/IC_*.md` point to notes of the original development that are not included in this repository; only `docs/IC_FALL_COUNTEREXAMPLES.md` is. The statements they describe are in the Lean sources and in `paper/`.

The [verification report](verification.json) records a successful Lean 4.19.0
build with pinned mathlib revision
`c44e0c8ee63ca166450922a373c7409c5d26b00b`. The
[build log](build.log) and [axiom log](axioms.log) are retained. The 513 audited (339 when this README was first written)
declarations use only `propext`, `Classical.choice`, and `Quot.sound`, with no
project axioms or proof placeholders. The verifier also checks source tokens,
toolchain version, source hashes, and the mathlib revision.

The prime-field result here is an exact
**presentation** counterexample: for E/F_5:y^2=x^3+2x, target (0,0), and
coordinate subgroup +/-1, the actual raw S3 presentation has LFD 4 and its
normalized presentation has LFD 3. The generated ideals are equal. This is
distinct from the larger, Python-checked FFD=4/LFD=5 counterexample.

The general closure/ideal soundness theorem, conditional normal-form
criterion, and graded-lex GB-to-normal-form bridge are formalized. This proves
the general finite-quotient theorem, explicit grid and Boolean bounds, and a
matching proper-ideal Boolean family in every dimension r>=2. The
[finite-quotient note](../../docs/IC_FINITE_QUOTIENT_LEAN.md) states the exact
hypotheses. The executable sparse proof-tree checker and actual F_5 certificate
replay are now verified; [sparse-demo.log](sparse-demo.log) records the runtime
diagnostic. The shared-DAG and containing-space checkers now verify the actual
F_257 normalized curve system: affine first strict fall 3 and exact LFD 4.
An independent ordinary Koszul FFD proof now gives FFD 3, including a
parameterized m=2,k=3 derivative-zero theorem for arbitrary coefficients.
See [the complete certificate explanation](../../docs/IC_CHECKED_PRIME_GAP.md).
The 58-row lower space and 177-node upper DAG are kernel-replayed. Concrete
arbitrary-parameter PKM quotient instantiation and the larger prime/binary
certificate imports remain obligations. No
generic sharp ECC bound or publication novelty is established by this build.

The new [parameterized binary/Koblitz family](../../docs/IC_BINARY_SHARP_FAMILY.md)
is also Lean-checked: normalized LFD 3, raw LFD 4, and equality of proper
ideals after actual linear coefficient descent. It has two summands, a
dimension-one coordinate base, and an order-two target. It is distinct from
the three older dimension-two finite binary certificates.

```sh
python scripts/verify_lean_fall_degrees.py
```

See [the project instructions](../../formal/README.md) for a cold build and
the exact coverage of each Lean module.

Two actual dimension-two descents are now included: ordinary F_8 and
non-invariant Koblitz F_16, each with exact normalized LFD 3. Lean proves
the finite-field constructions, power bases, coefficient maps, target facts,
factor-base dimension 2, and non-invariance under Frobenius. The
[satisfiable invariant case remains open](../../docs/IC_BINARY_DIMENSION_TWO_LEAN.md).


The [general prime-field follow-up](../../docs/IC_GENERAL_PRIME_LEAN.md) adds
all-k m=2 ordinary/truncated FFD classification and exact normalized LFD 4
for a growing-depth order-two-target curve family. Its separate
[verification report](../ic-general-prime/verification.json) checks current
audit coverage and three finite closure regressions.
