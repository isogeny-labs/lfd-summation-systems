# Frozen source custody after the rootfinder extension

The original preconditioning measurements pinned SHA-256 hashes for `csd/_rootfind.c` and `csd/ic_rootfind.py` in `environment.json`. PR #64 later extended those modules for the conditional-image experiment. The historical source bytes are preserved verbatim under `frozen-source/csd/`, recovered from the pre-extension main commit `ea196668897341d389c9bc99e45bfa89aa830100`. Their hashes match the unchanged frozen environment record:

| Path | Frozen SHA-256 |
| --- | --- |
| `csd/_rootfind.c` | `7b09ef3b21f2489699008d6086dce15e3632f6f4824ccee32d02c7365fc80fee` |
| `csd/ic_rootfind.py` | `ef0335bf39e1c54ce68292a2c471e6e585bd16ace0d09a1ab27c80490da2acc7` |

`scripts/check_semaev_preconditioned.py` now checks those immutable snapshots when the current source differs. It still regenerates the input cases and independently enumerates every point oracle on the current code. This is a source/record integrity check and a current-code semantic replay; it does **not** claim to reproduce historical timings or execute the historical native binary. The frozen environment, measurements and summary were not rewritten.
