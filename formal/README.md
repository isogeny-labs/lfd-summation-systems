# Lean foundation for degree-bounded relation systems

This project formalizes full iterated degree closure, an explicit bounded
normal-form criterion, the graded-lex finite-quotient theorem, explicit grid
bounds, and exact last fall degrees for a quadratic family.
It includes an actual prime-field summation-polynomial specialization with
LFD 3 after normalization and LFD 4 before normalization, despite equality of
the generated ideals. A subsequent parameterized binary/Koblitz family now
has Lean-checked normalized LFD 3 and raw LFD 4 under coefficient descent.
The general Boolean r+1 bound and a matching proper-ideal family in every
dimension r>=2 are also formalized. An actual normalized F_257 curve system now has kernel-checked affine first
strict fall 3, ordinary Koszul FFD 3, and exact LFD 4, using independent
homogeneous arguments and exported upper/lower certificates.
For every k>=2, the normalized m=2 prime-field presentation now has exact
ordinary and truncated FFD 3/4 according to `3*t^2+a=0`. A growing-depth
order-two-target prime curve family has exact normalized LFD 4. Broad sharp
bounds for arbitrary curves/targets and larger m remain open.

## Pinned build

The toolchain is Lean 4.19.0. The committed Lake manifest pins mathlib to
`c44e0c8ee63ca166450922a373c7409c5d26b00b` and pins its transitive dependencies.
Install the toolchain through elan, or use a preinstalled matching Lean/Lake.
From the repository root:

```sh
elan toolchain install leanprover/lean4:v4.19.0
cd formal
MATHLIB_NO_CACHE_ON_UPDATE=1 lake update
lake --no-cache build ICFallDegrees
lake env lean ICFallDegrees/Audit.lean
```

The cache-independent build can take substantially longer on first use.
For the build, source scan, revision checks, and theorem axiom audit together,
run from the repository root:

```sh
python scripts/verify_lean_fall_degrees.py
# Or: python scripts/verify_lean_fall_degrees.py --lake /absolute/path/to/lake
```

The verifier writes [the report](../experiments/ic-lean/verification.json),
[build log](../experiments/ic-lean/build.log), and
[axiom log](../experiments/ic-lean/axioms.log). The checked run covers 513
public declarations, whose proof
dependencies are subsets of mathlib's standard foundational axioms
`propext`, `Classical.choice`, and `Quot.sound`. No project axioms,
`sorry`/`admit`, unsafe implementation overrides, or `native_decide` occur.
Finite computations use kernel-checked reduction, `decide`, or proved `norm_num` steps.

## What is formalized

| Module | Mathematical content |
|---|---|
| [Closure](ICFallDegrees/Closure.lean) | Inductive V_(F,D), ideal and degree soundness, least closed-submodule property, transport with explicit preprocessing budget, and exact LFD from first derivation of 1 |
| [CyclicCarryWindow](ICFallDegrees/CyclicCarryWindow.lean) | Bounded periodic carries constructed from modular equality; actual odd uniqueness and alternating rigidity without a supplied witness |
| [CyclicProductBridge](ICFallDegrees/CyclicProductBridge.lean) | All-length telescoping products and the even alternating bridge at max(3,d+1) |
| [AffineTransport](ICFallDegrees/AffineTransport.lean) | Affine degree bounds and full-closure transport with checked input budgets |
| [BinaryTriple](ICFallDegrees/BinaryTriple.lean) | Exact LFD4 of the three-bit S3 obstruction, including full degree-three exclusion |
| [SplitChainLower](ICFallDegrees/SplitChainLower.lean) | Every normalized base-F2 split chain with m>=4 and either target bit has LFD>=4; the all-m upper4 specialization remains written |
| [CyclicNAF](ICFallDegrees/CyclicNAF.lean) | Periodic bounded-carry rigidity, alternating exceptions, and odd-cycle uniqueness under an explicit carry witness |
| [CyclicInverseRewrite](ICFallDegrees/CyclicInverseRewrite.lean) | Degree-three carry/adjacency derivations and quadratic-row multiplication with actual degree budgets |
| [BinaryFiber](ICFallDegrees/BinaryFiber.lean) | Actual characteristic-two S3 fiber identity, nonzero polynomial, and at-most-two-root theorem |
| [LocalReductionBound](ICFallDegrees/LocalReductionBound.lean) | A degree-bounded linear right section and local reductions through degree D imply global reductions and full-ideal LFD<=D; covers proper finite ideals |
| [BinaryS3Bounds](ICFallDegrees/BinaryS3Bounds.lean) | Actual raw/normalized S3 Boolean correction, basis extraction degree bounds, and arbitrary-target/graph split-system LFD<=max(r+1,4/3) |
| [NormalFormBound](ICFallDegrees/NormalFormBound.lean) | A linear section of an algebra quotient, bounded representatives, and explicit degree-bounded reductions imply an LFD bound |
| [LowerBounds](ICFallDegrees/LowerBounds.lean) | Top-degree containing-space and evaluation exclusion principles |
| [FirstFall](ICFallDegrees/FirstFall.lean) | Affine strict falls cannot exceed an LFD bound; this is distinct from a homogeneous-syzygy FFD definition |
| [QuadraticFamily](ICFallDegrees/QuadraticFamily.lean) | Over every field, c^2 != 1 implies exact LFD 3 and affine first strict fall 3 for {x^2-1,y^2-1,xy-c} |
| [RawPresentation](ICFallDegrees/RawPresentation.lean) | An explicit raw polynomial lift has exact LFD 4, with the same ideal as the normalized family |
| [PrimeCounterexample](ICFallDegrees/PrimeCounterexample.lean) | Over F_5, E:y^2=x^3+2x, target (0,0), and coordinate subgroup +/-1, the raw lift is the actual S3 polynomial; its LFD is 4 versus normalized LFD 3 |
| [Certificate](ICFallDegrees/Certificate.lean) | A degree-aware semantic certificate checker implies V-space and ideal membership |
| [PointIdeal](ICFallDegrees/PointIdeal.lean) | Budgeted coordinate differences imply degree-bounded reduction for every polynomial, a full-ideal bound, and equality with an evaluation kernel |
| [BooleanPair](ICFallDegrees/BooleanPair.lean), [BooleanRaw](ICFallDegrees/BooleanRaw.lean) | For any field and coefficient family with a nonzero member, the Boolean systems have exact normalized LFD 3 and raw LFD 4, with the same proper point ideal |
| [BinaryFamily](ICFallDegrees/BinaryFamily.lean) | Actual characteristic-two S3 identity, linear coefficient descent, exact degrees for a dimension-one coordinate base in every extension, Koblitz specialization, and a separate Frobenius criterion |
| [FiniteQuotient](ICFallDegrees/FiniteQuotient.lean) | Constructs unique normal forms and their degree-bounded derivations from graded-lex division and the GB property; proves LFD<=max(sigma+1,input degree) |
| [GridQuotient](ICFallDegrees/GridQuotient.lean) | Explicit nonempty finite grids and Boolean field-equation corollaries, with no assumed GB or normal-form data for these systems |
| [QuotientLowerBound](ICFallDegrees/QuotientLowerBound.lean), [BooleanSharp](ICFallDegrees/BooleanSharp.lean) | Closed quotient-span exclusion invariants and the exact proper-ideal Boolean family LFD=r+1 for all r>=2 |
| [SparsePolynomial](ICFallDegrees/SparsePolynomial.lean), [SparseCertificate](ICFallDegrees/SparseCertificate.lean) | Executable sparse arithmetic and a checker with proved interpretation, degree-budget, and full V-space soundness |
| [SparsePrimeCertificate](ICFallDegrees/SparsePrimeCertificate.lean) | Actual F_5 input semantics and kernel-replayed normalized/raw unit certificates, combined with independent lower proofs to recover exact LFD 3/4 |
| [VariableClosure](ICFallDegrees/VariableClosure.lean), [FilteredSpan](ICFallDegrees/FilteredSpan.lean) | Variable-wise closure implies closure under every allowed polynomial multiplier; filtered RREF spaces give general exclusion proofs |
| [SparseCoefficient](ICFallDegrees/SparseCoefficient.lean), [SparseLowerCertificate](ICFallDegrees/SparseLowerCertificate.lean), [LowerCertificateSupport](ICFallDegrees/LowerCertificateSupport.lean) | Executable containing-space checker with input inclusion, pivot identities, and all low-row variable products checked; coefficient and low-degree exclusions |
| [SparseDAG](ICFallDegrees/SparseDAG.lean), [SparseTrustChecks](ICFallDegrees/SparseTrustChecks.lean) | Shared derivation DAG and claimed-row verification, proved sound to the full V-space; mutation checks reject bad pivots, rows, references, and budgets |
| [GeneratedPrimeLower](ICFallDegrees/GeneratedPrimeLower.lean), [GeneratedPrimeUpper](ICFallDegrees/GeneratedPrimeUpper.lean) | The actual F_257 input list, 58-row degree-3 containing space excluding 1, and 177-node degree-4 derivation of 1, all kernel-replayed |
| [Prime257System](ICFallDegrees/Prime257System.lean), [GeneratedPrimeFirstFall](ICFallDegrees/GeneratedPrimeFirstFall.lean), [Prime257Counterexample](ICFallDegrees/Prime257Counterexample.lean) | Actual curve/S3 normalization identity, exact LFD 4, and affine first strict fall 3 on the same fixed normalized presentation |
| [ChainQuotient](ICFallDegrees/ChainQuotient.lean) | The m*k subgroup-chain equations of m summands at depth k are a reduced-ideal-trivial subsystem over every field with 2≠0 (two-power expansion of chain variables, binary-digit injectivity, roots of unity in an algebraic closure and the grid vanishing theorem); hence every input list containing them has LFD ≤ max(m*k+1, delta), and the normalized two-summand prime presentation of `PrimeM2General` has LFD ≤ 2k+1 at every depth k ≥ 2, prime p ≥ 5, curve and target |
| [CurveBounds](ICFallDegrees/CurveBounds.lean) | One kernel-checked entry point per family: prime chains (general m), prime two-summand at every depth, binary descended left-chain S3 systems (LFD ≤ r+1), and the Koblitz specialization b = 1 |

The [homogeneous FFD guide](../docs/IC_HOMOGENEOUS_FFD_LEAN.md) covers four
additional modules: `HomogeneousFirstFall` defines ordinary syzygies modulo
the full polynomial Koszul submodule; `QuadraticHomogeneousFamily` proves
matching degree-3 bounds with arbitrary lower tails; `PrimeM2Homogeneous`
connects the m=2,k=3 derivative-zero S3 branch for arbitrary coefficients;
`Prime257Homogeneous` identifies the actual F_257 inputs.

The [general prime-field guide](../docs/IC_GENERAL_PRIME_LEAN.md) covers
`PowerIdeals`, `LowDegreeCoefficients`, `GeneralQuadraticTop`,
`GeneralQuadraticClassification`, and `PrimeM2General` for arbitrary-k FFD;
`TruncatedFirstFall` and `PrimeM2Truncated` for the quotient convention;
`SquareChainRefutation`, `QuadraticNoLow`, and `PrimeM2Sharp` for matching
LFD-4 bounds at every depth; and `PrimeFieldConclusions` for prime-field
entry points, nonsingularity, and three concrete admissible instances.

The [dimension-two binary guide](../docs/IC_BINARY_DIMENSION_TWO_LEAN.md)
covers actual F_8/F_16 constructions, proved irreducibility and power bases,
coefficient extraction, the S3 Boolean correction, target equations,
dimension-2 bases and Frobenius witnesses, and exact LFD-3 unit certificates.
These are the `BinaryExtensionFields`, `BinaryDescent`,
`BinaryDimensionTwoSystems`, `BinaryDimensionTwoCurves`,
`BinaryDimensionTwoBases`, and two `GeneratedBinary*` modules.

Closure is defined over any commutative coefficient ring. The multiplier rule
charges `totalDegree a + totalDegree f <= D` before multiplication, including
when f is a newly derived lower-degree polynomial. Over a field this agrees
with the usual product-degree budget for nonzero factors. Arbitrary rings may
have zero divisors, so no equivalence with a product-degree convention is
claimed for them.

## Explicit remaining obligations

- The graded-lex GB-to-normal-form bridge is formalized and now instantiated
  for the direct prime chains: `ChainQuotient.reduced_chain_ideal_trivial`
  proves the reduced-ideal hypothesis for every m and k, and
  `ChainQuotient.chain_lfd_bound` gives LFD ≤ max(m*k+1, delta) for every
  input list containing the chains; `prime_field_m2_lfd_bound` is the
  all-depth two-summand corollary 2k+1. The bound for m ≥ 3 is stated in
  terms of the input degree delta; a Lean statement of the normalized
  S_{m+1} degree for m ≥ 3 is not provided. See
  [the finite-quotient proof note](../docs/IC_FINITE_QUOTIENT_LEAN.md).
- The raw lift in `RawPresentation` is synthetic for general K,c. The separate
  F_5 identity in `PrimeCounterexample` is what connects one specialization
  to the curve's S3. The group law is not formalized here.
- The F_2689 FFD=4/LFD=5 certificate and the satisfiable Frobenius-invariant
  F_8 dimension-two profile remain exact Python results. The ordinary F_8
  and non-invariant Koblitz F_16 dimension-two cases now have curve-connected
  Lean proofs of exact LFD 3. The new dimension-one family is a separate parameterized Lean
  theorem; see [its precise scope](../docs/IC_BINARY_SHARP_FAMILY.md).
- The legacy `checkCertificate` is noncomputable. The separate sparse checker
  is executable and sound, with actual prime certificates checked by kernel
  reduction. Shared DAG traces, deterministic Python exporters, and a general
  executable containing-space checker now replay the F_257 gap. Larger
  F_2689 and the satisfiable invariant binary import remain open. See
  [the sparse-certificate guide](../docs/IC_SPARSE_CERTIFICATES.md).
- Prime m=2 FFD is classified for every k>=2 in both ordinary and truncated
  conventions. The growing-depth prime LFD-4 family is formalized but uses
  an order-two target and inconsistent systems. Binary filtered FFD,
  arbitrary-target/m>=3 sharp LFD, and generic mk+1 sharpness remain open.
  The growing Boolean sharpness family is also formalized; no arbitrary-r
  realization of that Boolean family as a curve system is asserted.
- The new binary family has an order-two target already in the factor base
  and permits repeated coordinates. Rational lifting needs Tr(a)=0. The trace
  and group-law interpretation are written/finite-checked, not Lean group-law
  theorems.

See the [theorem ledger](../docs/IC_THEOREM_LEDGER.md) for all hypotheses and
the [novelty audit](../docs/IC_NOVELTY_AUDIT.md) for relevant prior work,
including existing Lean Gröbner-basis and polynomial-certificate projects.
This artifact is a foundation for further research, not a publication-ready
proof of broad sharp ECDLP degree bounds.

The [broad campaign](../docs/IC_BROAD_CAMPAIGN.md) and [inverse-family note](../docs/IC_BINARY_INVERSE_FAMILY.md) separate exact Python data and written candidate theorems from these kernel-checked general bounds. The inverse-family Moore bridge, growing lower bound, and concrete local matrix instances are not yet formalized.

The [proof-chain follow-up](../docs/IC_PROOF_CHAIN_REVIEW.md) adds 44 audited
results: arbitrary-dimensional inverse coordinate identities, actual Moore
matrix invertibility, exact finite coefficient-field extension, and exact
affine/input transport. The universal polynomial normalization, lower
obstruction and concrete filtered-input equivalence specialization remain
open. The aggregate reproduction command builds the pinned project alongside
the complete saved arithmetic and manuscript.
