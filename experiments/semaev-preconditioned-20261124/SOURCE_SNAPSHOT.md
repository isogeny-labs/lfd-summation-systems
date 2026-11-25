# Frozen source snapshot for the preconditioning record

The saved `environment.json` hashes the source used for this experiment.
Subsequent root-finder changes altered two files in the current source tree,
so checking those paths only at HEAD made the certificate self-test fail even
though the saved record was unchanged. The original bytes are archived under
`source_snapshot/` and are checked against the unchanged SHA-256 values in
`environment.json`.

| Path | SHA-256 | Git source |
|:--|:--|:--|
| `csd/_rootfind.c` | `7b09ef3b21f2489699008d6086dce15e3632f6f4824ccee32d02c7365fc80fee` | `5d2da6f10652e5d2099ba478f9616bec26f6a058` |
| `csd/ic_rootfind.py` | `ef0335bf39e1c54ce68292a2c471e6e585bd16ace0d09a1ab27c80490da2acc7` | `5d2da6f10652e5d2099ba478f9616bec26f6a058` |

`scripts/check_semaev_preconditioned.py` verifies these archived bytes when
the current source differs. It still regenerates the saved inputs using the
current implementation and checks the independent point oracle, recorded
rows, and summaries. This source snapshot does not replay the original
timings or transfer them to the newer implementation.
