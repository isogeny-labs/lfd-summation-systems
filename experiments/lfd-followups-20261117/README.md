# LFD validation follow-ups, 17 November 2026

Draft PR #31 tracks all five requested follow-ups. This directory contains
new experiments and an internal code review; it does not replace the
submission package merged in #24. No paid worker has been launched; the owner selected local-only.

## Results so far

| Track | Completed work | Remaining obligation |
|---|---|---|
| Degree-6 input certificate | Pure-integer DAG checker; six small closures checked against the direct engine; saved replayable certificate and malformed-input controls | Emit and check derivations of the actual depth-14 input basis; prototype is not that certificate |
| Independent lift review | Confirmed high-prime row-reduction defect, corrected it, added CI tests; four-prime arithmetic and 16 fresh lift runs pass | Separate external review / formal algorithm proof; input containment still required |
| Five summands, depth 3 | Tiled engine matches 16 control closures; bounded D8 retry stopped after 183.433 s at rank 1677 | Local ladder checked through D6, capped at D7; D8 needs a method that fits local resources |
| Proper-ideal exceptional primes | Exact LFD 5 at (k,p)=(7,257),(7,769),(8,7681),(9,12289); exact LFD 6 at (10,12289) | General proper-ideal lower bound remains open |
| Depth-18 breakpoint | Unit instance fixed and closed V4 computed: rank 3821, 1 excluded, hence LFD at least 5 | Checked exclusion at D7 and upper witness at D8; prediction 8 is untested |

The depth-7 p=257 profile agrees between the ordinary chain closure and
the compressed reference engine. At depth 10, p=12289, both the reference
and optimized engines completed checked closures and give exact LFD 6. The other new exact profiles currently
have the compressed reference engine only. Exactness uses the checked
closed space plus evaluation-rank/codimension criterion; these new values
are finite computational evidence, not an all-depth result or Lean replay.

The initial depth-10 reference attempt stopped at 240 seconds with no
completed D6 closure. Its record is retained separately. The extended
attempt completed in 977.99 seconds, within its 1800-second cap; the
optimized engine completed and checked the same exact value in 409.831
seconds. These local timings are not a controlled performance comparison.
An independent involution enumeration agrees with all five complete zero
sets, including the 86 points at depth 10.

## Certificate scope

The saved certificate is `certificate-k2-D4.json`. Its checker accepts an
externally supplied generator list and parameters, then checks input,
linear-combination and admissible-shift nodes without importing producer
code. It also verifies every supplied output row against the DAG. It
certifies containment of those rows, not completeness or missing large
basis rows. The checker requires an expected output digest from the caller;
an empty or substituted output set cannot choose its own target identity.
For a full depth-14 check, that target must come from the actual basis,
not from the certificate being checked. The direct comparison additionally checks row independence
and equality of rank with the reference closure.

These scripts use the normalized generator as the declared input. They do
not silently admit a higher-degree raw generator at a lower budget.
The inadmissible-carry control specifically rejects a multiplication whose
normal form drops degree but whose source degree already equals D.

## Reproduction

From the repository root:

```
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python3 -m pytest -q tests
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python3 experiments/lfd-followups-20261117/validate.py
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python3 experiments/lfd-followups-20261117/split_retry.py
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python3 experiments/lfd-followups-20261117/exceptional_reference.py
python3 scripts/check_publication_artifacts.py
```

`validate.py` and the other experiment scripts rewrite their own output
records and timing fields. Replay in a separate checkout to retain the
saved originals. The archival `calibration/split_retry.py.txt` is the exact
runner source of the initial three-minute retry, before the full-ladder
option was added. Restore it to its original path in a separate checkout
if replaying that calibration. Its hash matches the calibration record.

See [PROTOCOL.md](PROTOCOL.md) for identities, predictions, stopping rules,
chronology and the declined US$100 split-worker proposal and local-only amendment; see
[REVIEW.md](REVIEW.md) for the arithmetic counterexample and review limits.

## Resource boundary

The tiled full-column basis at D8 has a worst-case size of about 191 GiB,
versus about 382 GiB for one float64 full-column basis, before the old
engine's full-size update/copy temporaries. The tiled version also removes
the queued dense row copies. These are storage estimates, not measured
peak usage or a completion guarantee.

The depth-18 direct D8 space has 40,999,516 columns. A compressed split
basis has a worst-case storage estimate near 764 TiB before buffers. A
structured method is needed before allocating a larger run. No paid budget is approved. Current free disk space also rules out storing
these large bases locally.

## A candidate proof route

The explicit witness h=xy-x-y-8 belongs to the finite chain ideal, and
the new probe finds it outside V4 at depths 7..9 and outside V5 at depth
10. The first-derivative commutator passed 61,888 monomial checks at
depths 2..6. Two elementary supporting identities compile in Lean with
only propext. See PROPER-IDEAL-ROUTE.md for the proposed invariant and
the unproved bridge to the growth splitting theorem. The witnesses were
selected after observing the exact profiles; they are structural probes,
not preregistered predictions or an asymptotic theorem.

Native macOS RSS fields in the historical optimized engine use a Linux
divisor and are mislabeled. local-memory-units.json preserves the raw
fields and records the correct MiB interpretation for these local runs.

## Final local split outcome

The 1800-second local ladder cap produced checked closed exclusions at
D=3,4,5,6 with ranks 1,39,711,5243. At D7 the derived span reached rank
3644 and hit the deadline before closure; total measured time was
1806.759 seconds because the current batch finished before the stop.
This fresh local ladder establishes only a lower bound of 7. It does not
replace the older log-supported lower bound 8, and it does not establish
D8 membership. Partial profile fields describe a derived span, not V_D.
