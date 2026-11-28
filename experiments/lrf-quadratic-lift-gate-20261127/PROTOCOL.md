# Frozen pilot: quadratic root-coordinate lift against native LRF

## Question and evidence boundary

The existing [quadratic-lift screen](../quadratic-root-lift-screen-20261125/SCREEN.md)
proves that the residual two-point LRF condition has equations of input
degree at most two in kernel and root coordinates. This pilot asks whether
solving that presentation is actually cheaper than the native C kernel
on the **same residual queries**. Input degree is not solver degree, and
residual-query cost is not complete ECDLP cost.

The active index's statement that high projected ANF degree proves the
good set cannot be found without scanning is inconsistent with the
linked screen; correct that wording when reporting this pilot.

## Frozen source and fixture

- `csd/ic_rootfind.py` SHA-256
  `727a14d638964fa19c487294a8d281c2cf87d3ef47fecefa6cefe7469fd58cc5`.
- `csd/ic_dreg.py` SHA-256
  `0e5a650d3717ac029e3cbbc10c8e7f1525c11f46461d651187ce40d89c25afe1`.
- `experiments/ic-rootfind-yield/goodset_structure.py` SHA-256
  `dc906ad70c54e63b6f324307ffd1cbfb46ad4c4658566a38686329f6eb2cf762`.
- Field and curve construction follow `goodset_structure.py` with seed 1:
  binary fields `(n,k,delta)=(23,10,6)` and `(23,11,9)`, modulus
  `w^23+w^5+1`, `b=Random(1).randrange(1,2^23)`, and
  `E:y²+xy=x³+b`. The base is `<1,w,...,w^(k-1)>`.
- Generate residual targets in the same deterministic sequence as that
  script. Among valid nondegenerate, solvable systems with the declared
  kernel dimension, retain the first two with a relation and first two
  without one per cell, using complete native-kernel enumeration solely
  for this predeclared selection and independent audit. Stop after 5,000
  attempted targets per cell; record a missing category as censored.

These conditioned small panels are a diagnostic, not a population sample
or scaling result.

## Predeclared checks and costs

1. Build the n Boolean coordinate equations
   `x(z)^2+e1(g)x(z)+e2(g)=0`, plus Boolean field equations for all
   `g` and `z`. Check input degree at most two and compare the **entire**
   projected `(g,z)` solution set by direct enumeration with the native
   root-membership set. Lift every returned relation to actual points
   and verify group addition and factor-base membership.
2. Solve each lifted ideal with SageMath 10.10.rc0 / Singular over GF(2)
   in degree-reverse-lexicographic order. Record maximal degree seen if
   the backend exposes it; otherwise report it as unmeasured, never as
   two. Cap each solver call at 30 seconds, retain timeout and memory
   failures, and record wall and CPU. Include equation construction,
   Gröbner computation, solution extraction, point recovery, and final
   checks. Record solver setup and shared LRF setup separately.
3. Time the native C `Solver.lrf(2,T)` on the same targets and host,
   including all failed queries and relation recovery. Record its wall
   and CPU time. Do not pool the positive and empty categories without
   reporting each.

The pilot passes only if the lift has identical verified relations and
lower complete residual-query time on the frozen cells. Even a pass does
not establish increased relation yield or a same-target ECDLP win;
that would require a cold relation-collection, linear-algebra, and
generic-solver comparison on matched targets and factor bases.
