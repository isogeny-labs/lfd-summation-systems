# Matched torsion crossover result

The frozen protocol is in `PROTOCOL.md`; `run.py` is the timed producer and
`results.json` contains all 24 cells. `validate.py` provides a separate,
uncapped correctness replay and writes `validation.json`. Run from the repo
root with `python3 experiments/torsion-crossover-20261125/run.py` and then
`python3 experiments/torsion-crossover-20261125/validate.py`.

For this curve, target panel, factor base, and implementation, there is no
three-repeat crossover at N=8 through 1024. Independent relation rank is zero
at N=8 after sign reduction and positive at the other seven rungs. The torsion
arm completed and beat both direct arms in one of three repeats at each of
N=64 and N=256. The other cells were either slower when complete or
stopped under the predeclared dominance rule. A censored cell is only a lower
bound on the torsion arm's full cost; it is not a full timing.

The direct controls use identical points and targets. The shared fixture and
target setup is reported and charged equally to every arm. `validation.json`
records that uncapped torsion enumeration recovers exactly the same point-pair
sets as direct group scanning for all 288 targets at each N. The rank is for
the random subgroup-target relations after cofactor projection; planted
targets are correctness controls. The result is limited to two summands and
this curve, range, and Python implementation.

`invalid-pre-rational-base-results.json` and
`invalid-pre-rank-correction-results.json` and
`invalid-pre-sign-reduction-results.json` preserve discarded development
receipts. The first included abscissae outside the rational factor base in the
torsion arm. The second used two formal challenge columns even though Q=113G.
The third treated opposite points as separate rank columns.
Neither is used for the conclusion.
