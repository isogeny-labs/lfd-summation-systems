# Producer-local two-place kernel evidence

The tested source commit is recorded in receipt.json. The subsequent commit only
adds this evidence directory. Run from formal/ with the pinned Lean toolchain:

```sh
lake build ICFallDegrees.GrowthNonTorsionFunctionField ICFallDegrees.GrowthNonTorsionTargetSecant
lake env lean N8NonTwoTorsionPolesAudit.lean > poles.log
python3 ../scripts/check_non_two_poles_axioms.py poles.log
```

statements-and-axioms.log includes all 42 exact statements and transitive axiom
answers. audit-checks.txt records their acceptance, the earlier 16 checks, and
negative parser checks. build.log is a successful targeted build, not a full
default-library rebuild. The regression logs are from independent Lean commands.

failed-elaborations.log and failed-index.json preserve draft compiler failures and
one intentionally interrupted dependency build. They are debugging history, have
no retained draft source snapshots, and are excluded from every successful audit.
They must not be interpreted as kernel evidence. No numerical study ran.
