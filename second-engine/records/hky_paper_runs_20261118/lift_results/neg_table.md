| control | exit | outcome | dim V_D read | dim V_{D+1} | 1 in V_{D+1} | random closedness: shifted combos leaving | generator check: not in |
|---|---|---|---|---|---|---|---|
| unit clean | 0 | done | 829 | 3473 | true | 0 | 0 |
| unit clean, include-low | 0 | done | 829 | 3473 | true | 0 | 0 |
| unit flip-byte:2000 (raw) | 2 | refused: mask 0x20a9 has weight 5 > D | - | - | - | - | - |
| unit drop-row:0 | 2 | refused: the stored space is not closed: shifted random low elements leave it | - | - | - | 3 | 0 |
| unit drop-row:0 --no-refuse | 0 | done | 828 | 3473 | true | 3 | 0 |
| unit perturb:0 | 2 | refused: the stored space is not closed: shifted random low elements leave it | - | - | - | 3 | 0 |
| unit perturb:0 --no-refuse | 0 | done | 829 | 3473 | true | 3 | 0 |
| unit drop-row:100 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 12 | 3 |
| unit drop-row:100 --no-refuse | 0 | done | 828 | 3473 | true | 12 | 3 |
| unit perturb:100 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 12 | 3 |
| unit perturb:100 --no-refuse | 0 | done | 829 | 3473 | true | 12 | 3 |
| unit drop-row:400 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 36 | 3 |
| unit drop-row:400 --no-refuse | 0 | done | 828 | 3473 | true | 36 | 3 |
| unit perturb:400 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 36 | 3 |
| unit perturb:400 --no-refuse | 0 | done | 829 | 3473 | true | 36 | 3 |
| unit drop-row:700 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 42 | 3 |
| unit drop-row:700 --no-refuse | 0 | done | 828 | 3473 | true | 42 | 3 |
| unit perturb:700 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 42 | 3 |
| unit perturb:700 --no-refuse | 0 | done | 829 | 3473 | true | 42 | 3 |
| unit drop-row:828 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 42 | 4 |
| unit drop-row:828 --no-refuse | 0 | done | 828 | 3473 | true | 42 | 4 |
| unit perturb:828 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 42 | 4 |
| unit perturb:828 --no-refuse | 0 | done | 829 | 3473 | true | 42 | 4 |
| unit hook: no fall shifts | 0 | done | 829 | 2844 | - | 0 | 0 |
| unit hook: drop shift 0 | 0 | done | 829 | 3473 | true | 0 | 0 |
| unit hook: drop shift 13 | 0 | done | 829 | 3473 | true | 0 | 0 |
| proper clean | 0 | done | 829 | 3471 | - | 0 | 0 |
| proper clean, include-low | 0 | done | 829 | 3471 | - | 0 | 0 |
| proper flip-byte:2000 (raw) | 2 | refused: mask 0x20a9 has weight 5 > D | - | - | - | - | - |
| proper drop-row:0 | 2 | refused: the stored space is not closed: shifted random low elements leave it | - | - | - | 3 | 0 |
| proper drop-row:0 --no-refuse | 0 | done | 828 | 3471 | - | 3 | 0 |
| proper perturb:0 | 2 | refused: the stored space is not closed: shifted random low elements leave it | - | - | - | 3 | 0 |
| proper perturb:0 --no-refuse | 0 | done | 829 | 3473 | true | 3 | 0 |
| proper drop-row:100 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 12 | 3 |
| proper drop-row:100 --no-refuse | 0 | done | 828 | 3471 | - | 12 | 3 |
| proper perturb:100 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 12 | 3 |
| proper perturb:100 --no-refuse | 0 | done | 829 | 3473 | true | 12 | 3 |
| proper drop-row:400 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 36 | 3 |
| proper drop-row:400 --no-refuse | 0 | done | 828 | 3471 | - | 36 | 3 |
| proper perturb:400 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 36 | 3 |
| proper perturb:400 --no-refuse | 0 | done | 829 | 3473 | true | 36 | 3 |
| proper drop-row:700 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 42 | 3 |
| proper drop-row:700 --no-refuse | 0 | done | 828 | 3471 | - | 42 | 3 |
| proper perturb:700 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 42 | 3 |
| proper perturb:700 --no-refuse | 0 | done | 829 | 3473 | true | 42 | 3 |
| proper drop-row:828 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 42 | 4 |
| proper drop-row:828 --no-refuse | 0 | done | 828 | 3471 | - | 42 | 4 |
| proper perturb:828 | 2 | refused: the generator g (or some m g) is not in the stored space | - | - | - | 42 | 4 |
| proper perturb:828 --no-refuse | 0 | done | 829 | 3473 | true | 42 | 4 |
| proper hook: no fall shifts | 0 | done | 829 | 2844 | - | 0 | 0 |
| proper hook: drop shift 0 | 0 | done | 829 | 3471 | - | 0 | 0 |
| proper hook: drop shift 13 | 0 | done | 829 | 3471 | - | 0 | 0 |
