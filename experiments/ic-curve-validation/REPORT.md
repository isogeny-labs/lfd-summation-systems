# Validator report: independent re-verification of exact last fall degrees (LFD)

Scope: `/home/user/ml-cryptanalysis`, evidence `experiments/ic-curve-sweep/{systems,results,cases}.jsonl` and
`experiments/ic-depth-step/growth.json`. Nothing in the repository was edited, committed or pushed.
All code is in `/tmp/claude-0/-home-user/a6e79f2a-a5ec-5904-b377-343a14909fff/scratchpad/` (list and hashes at the end).

## Disposition: `validated` (scoped to the claims below)

* 1440 of 1440 systems recomputed (every prime k=2 system, every prime k=3 system, every Boolean system with 2k <= 8,
  both `binary` and `koblitz` cases): **0 mismatches** in exact LFD, zero-set size, unit/proper status, or the list of
  strict falls. Required minimum was 12 + 12; all of them were done, with 17+18 prime and 16+16+14 Boolean in the stratified
  tables of section 3 and the full populations in section 3.4.
* Item 3: ordinary-ring dimension of V_3 is 49 = 42 + 7 for k = 3, consistent with "dim V_3 = 2k+1 modulo the chains".
  The two dimension counts are related by an exact identity that I verified (section 5). All 30 `growth.json` rows
  (5 curves, k = 3..8) reproduce dim V_3 = 2k+1 and the solution counts; all 15 rows with k <= 5 reproduce dim V_4.
* Not claimed: anything about k >= 7 (not run) and about k = 4..6 beyond the 14 bonus systems of section 6; anything about performance, rho, the IC leaderboard
  or the validation gate (none of these is in the dispatch, see section 7).

## 1. Binding

| item | value |
|---|---|
| commit at start of validation | `553e01c7fcda4e4d45107f5694687ef9a56dc32e` |
| commit at end of validation | `7e4f2191e62db43959989d7bbf24d495f863d3ed` (another session committed paper/Lean/script work meanwhile) |
| evidence blobs, identical at both commits AND in the work tree (`git rev-parse <c>:<f>` == `git hash-object <f>`) | `systems.jsonl eb683b31...`, `results.jsonl 70e366cc...`, `cases.jsonl ded53d7f...`, `growth.json c87d492c...` |
| `git diff 553e01c7 HEAD -- experiments/ic-curve-sweep experiments/ic-depth-step formal/ICFallDegrees/Closure.lean` | empty |
| sha256 of data files | equal to the values recorded in `experiments/ic-curve-sweep/verification.json` and `experiments/ic-depth-step/verification.json` (systems `99248b43...`, results `b7a88e88...`, cases.jsonl `8143967b...`, growth `d6e52ea6...`) |
| tree state | At start the tree was dirty (`scripts/run_ic_chain_depth10_probe.py`, `paper/curve-lfd-bounds.{tex,pdf}`, `experiments/ic-curve-paper/verification.json`: none is evidence I used). At the end `git status` was clean. The dispatch named no commit; I bound to the byte-identical blobs above, not to a moving tree. |
| environment | Python 3.13.18, numpy 2.5.3, Linux 6.18.44, 4 shared cores (two other long jobs running), `nice -n 10..15`, 1-2 BLAS threads. Peak RSS: 79 MB for the required parts, 97 MB for the k=5 bonus, 867 MB for the k=6 bonus. Longest required run: 0.7 s. |
| seeds | system selection `random.Random(20261105)`; control `Random(5)`; naive sample `Random(99)`. Exhaustive audit needs no seed. |

## 2. What I implemented (independent of the repository's closure engines)

Code: `lfd_engine.py`. It imports nothing from `csd/`. The only repository code used anywhere is, in `item3.py` growth mode,
`csd.ic_dreg.PrimeCurve` and `csd.ic_dreg.pkm_generator` to build the growth.json inputs (as allowed). The main audit
reads the JSON files only. I did NOT use `ic_fall_certificates`, `ic_chain_shift_closure`, `ic_boolean_closure`, or the
`ic_curve_sweep` profiles. (I did read `ic_curve_sweep.py` and `ic_dreg.py` to learn the file format: bit j of block i of a
generator mask = variable u_{i,j}; Boolean equations = lists of squarefree masks XORed.)

Definitions, from `formal/ICFallDegrees/Closure.lean`, in the ORDINARY ring K[z_1..z_r] (full monomial space of degree <= D,
r = 2k, no reduction by the chains or by z^2 = z):

* Inputs, prime: `u_{i,j}^2 - u_{i,j+1}` (j < k-1), `u_{i,k-1}^2 - 1` for i = 0,1, plus the stored normalised generator G rebuilt
  as a sum of coefficient times product of u over the set bits. Binary/Koblitz: `z^2 - z` for each of the 2k variables plus every
  stored equation.
* `V_D`: start with the span of inputs of degree <= D, repeat "multiply every element of V of degree <= D-1 by every variable"
  to a fixed point (dense exact elimination mod p or GF(2), columns ordered by degree descending, so the pivot is the leading
  degree and the elements of degree <= D-1 are exactly the rows with pivot degree <= D-1).
* `I_D = I cap P_{<=D}`. Zero set Z found by brute force from G and the chain (all x_i with x_i^{2^k} = 1, u_{i,j} = x_i^{2^j}) or
  from the 2^{2k} Boolean points. The ideal is radical and every point is rational: K[u]/(chain) is isomorphic to
  K[x_1,x_2]/(x_1^N - 1, x_2^N - 1) with N | p-1, which is K^{N^2}, reduced, and any quotient of a product of copies of K is again
  one; likewise K[z]/(z^2 - z). So `dim I_D = |M_{<=D}| - rank(evaluation of M_{<=D} on Z)`.
* `LFD = 1 + max{ D : V_D != I_D }` (equivalent to the Lean `HasLFDBound`: every h in I lies in `V_{max(d, deg h)}`, because
  `V_D` is inside `I_D`). Cross-check: the largest "fall" D, where `dim(V_D cap P_{<=D-1}) > dim V_{D-1}`, which must give the same number.
* Certificate that LFD <= D0 (proof in the `lfd_engine.py` docstring): `dim V_{D0} = |M_{<=D0}| - |Z|` and the monomials of degree
  < D0 evaluate onto K^Z. Then V_D = I_D for every D >= D0. For the unit ideal the certificate is `1 in V_{D0}`.
  The upper bound needs no radicality (V <= I <= ker eval and the dimensions agree); radicality is only used for I_D at the degrees below D0.

### Self-tests and controls (all pass; `selftest.py`, `control.py`, `naive_check.py`, `naive_lfd.py`, `soundness_check.py`)

1. Row-echelon rank equals a pure-Python exact rank on random low-rank matrices for p in {2, 3, 97, 257, 12289}.
2. Hand-computed toy ideals: `{x^2-x}` LFD 0; `{x^2-1, x^2-x}` over F_3 LFD 2 (x-1 appears at degree 2); `{x, x-1}` LFD 1; `{x^2, x^2-1}` LFD 2.
3. Variable-wise closure equals the literal Lean rule (pure Python integers, own echelon, every monomial multiplier `a` with
   `deg a + deg f <= D` applied to every basis row): equal (dim V_D, dim V_D cap P_{<=D-1}, 1 in V_D) for D = 0..5 on 12 k=2 / 4-variable
   systems, and for D = 0..4 on two proper k=3 prime systems (one only to D = 3), one unit k=3 prime system and two 8-variable Boolean systems (log `log_naive.txt`).
4. A second, fully separate pipeline: `exact_lfd` driven by the pure-Python literal-rule closure on a stratified sample of 60
   systems (10 prime k=2, 20 prime k=3, 10+10+10 Boolean): 60/60 equal to the repository in LFD, point count and falls (`log_naive_lfd.txt`).
5. Soundness at the certificate: for all 534 proper systems of the audit, at the certified degree D0 every basis row of V_{D0}
   evaluates to 0 at every point of Z, `dim V_{D0} = |M| - |Z|`, and evaluation of degree < D0 monomials has rank |Z|:
   534/534 pass (`log_soundness.txt`). `dim V_D <= dim I_D` was asserted at every D of every run and never fired.
6. Sensitivity control (`control.py`): the engine does not return a constant. Dense random generators over p = 97, k = 3:
   degree 1 gives LFD 4, degrees 2 and 3 give LFD 5, degree 4 gives LFD 6 (4 samples each, all certified). So "LFD 4 for every sweep system with k = 3"
   is a property of the sweep generators (degree <= 2, only u_{i,0}, u_{i,1}), not an artefact of my code.
7. Input integrity (supplement, own S_3 formulas): for every audited prime system the stored normalised generator equals
   `S_3(c*u_1, c*u_2, x_R)` at every point of mu_N x mu_N (c = `base_scale`, curve and target from `cases.jsonl`) and its zero set equals the stored
   `coordinate_solutions` decoded through `generator_value`; for every Boolean system the stored equations vanish on exactly the
   Boolean assignments where `S_3` (GF(2^n) with the case's modulus, basis, b, x_R) vanishes, and the zero set equals the stored list. 0 failures in 1440.
   Caveat: these use the textbook S_3 and the case metadata; they show the stored systems mean what the repository says, not that the metadata is the intended problem.

## 3. Items 1 and 2: stratified tables (system_hash, mine / repo)

Every row: certificate found at the listed degree D0 (V_{D0} = I_{D0} with surjective evaluation, or 1 in V_{D0}); all zero-set sizes,
statuses, LFD and fall lists equal. "falls" is my list of D with `dim(V_D cap P_{<=D-1}) > dim V_{D-1}`; the repository calls these `strict_falls`.

#### Item 1a: prime, k = 2 (4 variables) (17 systems)

| system_hash | family | status | zero-set size (mine / repo) | exact LFD (mine / repo) | falls (mine / repo) | certificate at D | result |
|---|---|---|---|---|---|---|---|
| 44e2542d2ce965ec1c3c1fb8c0cb7ac86454f78f7cd0aca03fa3809cefa6b39d | prime k=2 p=193 | proper-local | 2 / 2 | 3 / 3 | [3] / [3] | 3 | match |
| 60534f932de0d996ecf9626d15f76e7b6efd8220bc736c3c389693375dc9a7b7 | prime k=2 p=97 | unit | 0 / 0 | 3 / 3 | [3] / [3] | 3 | match |
| e135e8cdb667da9ec542618f14248fb98dd36bd502bb1e464d422a397d19710b | prime k=2 p=97 | unit | 0 / 0 | 3 / 3 | [3] / [3] | 3 | match |
| d6b3997ae5bd392ed06a9f60112b69574d8b72bf7c1f9298b9887bcf67fdd58c | prime k=2 p=97 | proper-local | 2 / 2 | 4 / 4 | [4] / [4] | 4 | match |
| 600022908abcd7ed333c118a1fa5b2cafc63602e6bdb5029384bae78cce3a38e | prime k=2 p=97 | proper-local | 2 / 2 | 4 / 4 | [4] / [4] | 4 | match |
| 31e9b4bcf5d1dcff6186a382162b33086c6395d5a6107a8ad71cd0a52d5189c8 | prime k=2 p=193 | proper-local | 2 / 2 | 4 / 4 | [4] / [4] | 4 | match |
| de2209ac604084e48e7c55b673a7c5d1618f7b65c691cabb9867a708271af636 | prime k=2 p=193 | proper-local | 4 / 4 | 4 / 4 | [4] / [4] | 4 | match |
| c6d36621550034f80da562d08bca7861424dcd9f2b9efbdad876f63666354fe3 | prime k=2 p=257 | proper-local | 2 / 2 | 4 / 4 | [4] / [4] | 4 | match |
| fb2ab4d82ca2c293b2e9c1f91f3ee805ae41b114f793e875ad9e88c050aba6f8 | prime k=2 p=257 | proper-local | 4 / 4 | 4 / 4 | [4] / [4] | 4 | match |
| 29ca1bec6ebe47379ac0eaffe064c4b10c1e4898b82bbb8414fb5af61747554b | prime k=2 p=641 | proper-local | 2 / 2 | 4 / 4 | [4] / [4] | 4 | match |
| 099f0d4df65b54db0b83303e8f8a38e114704649f1e7e2ed982495dc53018360 | prime k=2 p=97 | unit | 0 / 0 | 4 / 4 | [4] / [4] | 4 | match |
| 5551ea27a3b1cc79260d8a6fa85313c8cee4f66f3aa8049a3d3d6a6027149b11 | prime k=2 p=97 | unit | 0 / 0 | 4 / 4 | [4] / [4] | 4 | match |
| b386a6469ecfbc48f0a9e291b7bc70e0a5d5d490d63e0a99eda027959ee2d727 | prime k=2 p=193 | unit | 0 / 0 | 4 / 4 | [4] / [4] | 4 | match |
| 1bf0149ad1ed03d44e859106c8afc92f049de83cf4bff45c2906cc1cc12d55c8 | prime k=2 p=193 | unit | 0 / 0 | 4 / 4 | [4] / [4] | 4 | match |
| 5db841657f16c3c90dcd2b0c5bf708e96b1731bcef731474427eaf66ef62f084 | prime k=2 p=257 | unit | 0 / 0 | 4 / 4 | [4] / [4] | 4 | match |
| 155b7e9bd8b6eb9575be1fa7945d15b779e33e812ff97eb4540b19c66e188f79 | prime k=2 p=257 | unit | 0 / 0 | 4 / 4 | [4] / [4] | 4 | match |
| 22ccadea0fbbaf049f244ef6b0ff15bd76f6d4c05527938d1a5968ae5f0cdcaf | prime k=2 p=641 | unit | 0 / 0 | 4 / 4 | [4] / [4] | 4 | match |

#### Item 1b: prime, k = 3 (6 variables) (18 systems)

| system_hash | family | status | zero-set size (mine / repo) | exact LFD (mine / repo) | falls (mine / repo) | certificate at D | result |
|---|---|---|---|---|---|---|---|
| 2dceabe91fba9ccf94619c5384e6c83298005182de98bd28e84c56e9590b476c | prime k=3 p=97 | proper-local | 2 / 2 | 4 / 4 | [4] / [4] | 4 | match |
| 8261f9d5e21e0492f7f70a5ad1f092e003f0cd10d8a8593c36795e80e0cc767f | prime k=3 p=97 | proper-local | 2 / 2 | 4 / 4 | [4] / [4] | 4 | match |
| 3ca5b4813c8f931818dec360efb0fe93a80d341602e25f3e3ef071a9eb9620e9 | prime k=3 p=193 | proper-local | 2 / 2 | 4 / 4 | [4] / [4] | 4 | match |
| ddd557fafbbc03fee48944eda541fb8de1faa2a5215bc821b983b6b322b0745a | prime k=3 p=193 | proper-local | 2 / 2 | 4 / 4 | [4] / [4] | 4 | match |
| abbf0a53b8745de8cb318d62740559e920ad1063148b971aef434c8fc9dde273 | prime k=3 p=257 | proper-local | 2 / 2 | 4 / 4 | [4] / [4] | 4 | match |
| 3738ef21f81c17605cbc17f87ba3d2ea2993679ebaeb8f1842555c816a8f7bb9 | prime k=3 p=257 | proper-local | 4 / 4 | 4 / 4 | [4] / [4] | 4 | match |
| 212717d28f7e1e6a76783847306fac1b5f4ecfaa3fae8bc3a058a251b3844c33 | prime k=3 p=641 | proper-local | 2 / 2 | 4 / 4 | [4] / [4] | 4 | match |
| 9edda4ab69dad428110ef7a7da11ae5f96294a21f9e388677b885703d917cc4b | prime k=3 p=641 | proper-local | 4 / 4 | 4 / 4 | [4] / [4] | 4 | match |
| 99247aa917040fd20c3888a777329e4ec3118b084962e2190ba80b657f1e4ba9 | prime k=3 p=7681 | proper-local | 2 / 2 | 4 / 4 | [4] / [4] | 4 | match |
| 44eb174d4df52648945a3bfa0e444213345e11550c3f356e7d636e1aac9c33ca | prime k=3 p=97 | unit | 0 / 0 | 4 / 4 | [4] / [4] | 4 | match |
| 368b0271bfcd6eb46b2cc776267e9805c12261cda1aae7675c00979167834f2b | prime k=3 p=97 | unit | 0 / 0 | 4 / 4 | [4] / [4] | 4 | match |
| d569e3e1392fc8ff17c63e311a19c3935d79426926fce00abc17ef7aa5e7fc2f | prime k=3 p=193 | unit | 0 / 0 | 4 / 4 | [4] / [4] | 4 | match |
| fc489857694fea4c1521bb3fabc68f2557518e872fde10afb60351d765f563c2 | prime k=3 p=193 | unit | 0 / 0 | 4 / 4 | [4] / [4] | 4 | match |
| b34d90542a9caadee3eb7864d288f54ef159289d39727c55c09325073bf463a9 | prime k=3 p=257 | unit | 0 / 0 | 4 / 4 | [4] / [4] | 4 | match |
| 8ba44093bf73bd28bf9daf86b77df0f32779693a2dd9be77b2fce84cde7ca447 | prime k=3 p=257 | unit | 0 / 0 | 4 / 4 | [4] / [4] | 4 | match |
| 6ccf36ee7cbd6fcc6908189eec73c6f3aa1d18236ef18e78ee4025aae3dc4a6a | prime k=3 p=641 | unit | 0 / 0 | 4 / 4 | [4] / [4] | 4 | match |
| 42ced8d20e585312248c1b1f695bdb70f963b6d9418382675d02434808a178a7 | prime k=3 p=641 | unit | 0 / 0 | 4 / 4 | [4] / [4] | 4 | match |
| c80c20420e1eb79881244b18999646bb426bd8725c445312c1e702fff891edc4 | prime k=3 p=7681 | unit | 0 / 0 | 4 / 4 | [4] / [4] | 4 | match |

#### Item 2a: Boolean, 4 variables (16 systems)

| system_hash | family | status | zero-set size (mine / repo) | exact LFD (mine / repo) | falls (mine / repo) | certificate at D | result |
|---|---|---|---|---|---|---|---|
| 45d9bc85fb75ce7e7a0b97f523a4016903f376d465d126f38aa2a1729c698531 | binary vars=4 | proper-local | 2 / 2 | 2 / 2 | [2] / [2] | 2 | match |
| 0a137d8469cd2882d75b0432d2ce54d68769945aa04e00b39aceaa8a6fe37913 | binary vars=4 | proper-local | 4 / 4 | 3 / 3 | [2, 3] / [2, 3] | 3 | match |
| 9dd942e95b1cf4d77418e29b5b7829b7423c1f075a75bb4c1f6b0d67dba37d6b | binary vars=4 | unit | 0 / 0 | 0 / 0 | [] / [] | 0 | match |
| 0a365999efec3e987444800a57ff75a098fb9236e9b8bff41898dd7dc60caa93 | binary vars=4 | unit | 0 / 0 | 1 / 1 | [1] / [1] | 1 | match |
| c66e40c7095c67d381f48d8a7b3fb9454f6e1bfb7e7cf5694fc4cc6932b2d87e | binary vars=4 | unit | 0 / 0 | 2 / 2 | [2] / [2] | 2 | match |
| 1cb8c6577b377e723cfe8afbd5f5d9a64a61c4e28c35d095c0bfd36c06dc51cf | koblitz vars=4 | proper-local | 6 / 6 | 0 / 0 | [] / [] | 3 | match |
| f9ae7b464bd8f9594cd3a4db2041aa231962cb52e98c74bea4e01f3bf020d909 | koblitz vars=4 | proper-local | 6 / 6 | 2 / 2 | [2] / [2] | 3 | match |
| b09951d509586bd76a205a3151bcf36bd5d735f92fa99b477b2eb605a32cf824 | koblitz vars=4 | proper-local | 3 / 3 | 3 / 3 | [3] / [3] | 3 | match |
| 3a766d21628c309c82eece6020308a292da4b7b4b49bcc9cd76b336c5754cc5c | koblitz vars=4 | unit | 0 / 0 | 0 / 0 | [] / [] | 0 | match |
| 91d0925384e905fee149950607345db2fc2e4c373237fa4de3d5545a1e178904 | koblitz vars=4 | unit | 0 / 0 | 1 / 1 | [1] / [1] | 1 | match |
| 50422d4a15b895c5b5b4343c37c4a648d76e47a4131cc0a597176386f6d14891 | koblitz vars=4 | unit | 0 / 0 | 2 / 2 | [2] / [2] | 2 | match |
| cd3b59cd8d108b38dd62a6f0c37b79c58f5b1ecc22b2cc2003f12674a7440889 | koblitz vars=4 | unit | 0 / 0 | 3 / 3 | [3] / [3] | 3 | match |
| 24ca263c117f23673733655c49c0488883e3d6591954b7def92031030fa5d49e | binary vars=4 | proper-local | 2 / 2 | 2 / 2 | [2] / [2] | 2 | match |
| c15e1d8ee56c955dc5717571bf9b6fef2db97f276f995d21a2427b91ce8c4e80 | binary vars=4 | proper-local | 4 / 4 | 3 / 3 | [2, 3] / [2, 3] | 3 | match |
| 6bd3476fa00170f91d308d68796ee004898ad8264e9cc041c69bad0bde57101f | binary vars=4 | unit | 0 / 0 | 0 / 0 | [] / [] | 0 | match |
| bf4c5909e3e7239f76c4a6cad8b2bc61a1686d93b05808732b1295d677412838 | binary vars=4 | unit | 0 / 0 | 1 / 1 | [1] / [1] | 1 | match |

#### Item 2b: Boolean, 6 variables (16 systems)

| system_hash | family | status | zero-set size (mine / repo) | exact LFD (mine / repo) | falls (mine / repo) | certificate at D | result |
|---|---|---|---|---|---|---|---|
| bf41b9cca7b3091a1c2caca5c2381452528c75b4b76f8c988791d36ef65765ee | binary vars=6 | proper-local | 2 / 2 | 2 / 2 | [2] / [2] | 2 | match |
| 79ffc649e4ea26e7b1446a9b87652ae76b618df73d4e17f14672f9be9cfec480 | binary vars=6 | proper-local | 5 / 5 | 3 / 3 | [3] / [3] | 3 | match |
| a07755fbf7982ed9a8cb7b08cb003c57bdbd9c533205b0a3811bbaf97b5ef973 | binary vars=6 | unit | 0 / 0 | 0 / 0 | [] / [] | 0 | match |
| 2534fd65de386f01d03aa7b6be2bde0381220b8234a47f707df6dcf59a964ee8 | binary vars=6 | unit | 0 / 0 | 2 / 2 | [2] / [2] | 2 | match |
| b87ba5d1ba4fe444cfe1cf3a095a93d8ba56d56f29475c26bea90ccce902389f | binary vars=6 | unit | 0 / 0 | 3 / 3 | [3] / [3] | 3 | match |
| 73bbec6e6e6a0684a52615f60c3be05d3c0a47e3b6deb9992a1c103c638fa2fa | koblitz vars=6 | proper-local | 2 / 2 | 2 / 2 | [2] / [2] | 2 | match |
| 493ec536ba485cbf2291b060ae5d90df002a74fc7c89b04bd6f83e9feb1d275d | koblitz vars=6 | proper-local | 2 / 2 | 3 / 3 | [2, 3] / [2, 3] | 3 | match |
| f75c4294318cc4dea30437b774cfc6e54401628e5a1bafa1add3df6cc442a419 | koblitz vars=6 | unit | 0 / 0 | 0 / 0 | [] / [] | 0 | match |
| 819aeb93e164e292975a3f1d24b304cfe17f8ff8bfd6c8d480189538d7e06a0d | koblitz vars=6 | unit | 0 / 0 | 1 / 1 | [1] / [1] | 1 | match |
| 379551bd07c0ffee57a16a4ce841b999006c4e16c70e845902e6cce372a7130e | koblitz vars=6 | unit | 0 / 0 | 2 / 2 | [2] / [2] | 2 | match |
| 2ff9acff41da6a5672a531f7c7cc61b4b2b0f283ad6910686eab2e3ef4710d25 | koblitz vars=6 | unit | 0 / 0 | 3 / 3 | [3] / [3] | 3 | match |
| c029140aad33bca62a4de4638003c2222dce9d7cfdaba1382516240d3f964a9f | binary vars=6 | proper-local | 2 / 2 | 2 / 2 | [2] / [2] | 2 | match |
| 2a87bb28de86a7982a71df6832c98f3d3e94ee4ca317132c4be5eb8ef02ff9c1 | binary vars=6 | proper-local | 2 / 2 | 3 / 3 | [3] / [3] | 3 | match |
| 1547e4e0ef07e4d2354a4345d0eb25dce315eb8d2fceecf44df24eb92550304b | binary vars=6 | unit | 0 / 0 | 0 / 0 | [] / [] | 0 | match |
| 97cdbb677d47ea8b08fe1811d2775cf5a0cf92c1368e3bd5daeb9bd9b86a2b8b | binary vars=6 | unit | 0 / 0 | 2 / 2 | [2] / [2] | 2 | match |
| 999601f185c54ac4d594a2b28ac1373e7859d8fda89975eeca6a372ce49db0f4 | binary vars=6 | unit | 0 / 0 | 3 / 3 | [3] / [3] | 3 | match |

#### Item 2c: Boolean, 8 variables (14 systems)

| system_hash | family | status | zero-set size (mine / repo) | exact LFD (mine / repo) | falls (mine / repo) | certificate at D | result |
|---|---|---|---|---|---|---|---|
| ed64ecd0cb5cb4a755f6b53615e443f8001b40b904e51416f55f334bf7b10bb1 | binary vars=8 | proper-local | 2 / 2 | 2 / 2 | [2] / [2] | 2 | match |
| f44266bac9f3eaae93ad8bdd51b986d587cedb542f89c080b6b95e97ab90edf8 | binary vars=8 | proper-local | 5 / 5 | 3 / 3 | [3] / [3] | 3 | match |
| 4f19eb4483d6ee4797a6971203ab4aa376a41367e94f3aa6918a819d7adafec8 | binary vars=8 | proper-local | 14 / 14 | 4 / 4 | [3, 4] / [3, 4] | 4 | match |
| 61c9f4cb7e817998a7cb9c6b12cbfd27560a9e7fb6c844df655a71280e1f84a2 | binary vars=8 | unit | 0 / 0 | 0 / 0 | [] / [] | 0 | match |
| 9f946b65a207a5c8c8d5f1043c1a0abeaa0997f0a4026a8d757622bac5176c66 | binary vars=8 | unit | 0 / 0 | 2 / 2 | [2] / [2] | 2 | match |
| 3d6a488f67e0171b0a158dce08122712b613f572fb016bc248850843f439869a | binary vars=8 | unit | 0 / 0 | 3 / 3 | [2, 3] / [2, 3] | 3 | match |
| 650974aeb829e2dd52f14f2bc31d8a71df0647061ae16f7e628091529ce26e3a | koblitz vars=8 | proper-local | 2 / 2 | 2 / 2 | [2] / [2] | 2 | match |
| e3b587c6afa4c6d28d037127a383970c7197bca4c7a49bab6cc44b32b289e54d | koblitz vars=8 | proper-local | 15 / 15 | 3 / 3 | [3] / [3] | 3 | match |
| 42591ae614563251d5dfc46bfda66bb21668aecbb316694aaf678900caa4e490 | koblitz vars=8 | proper-local | 14 / 14 | 4 / 4 | [2, 3, 4] / [2, 3, 4] | 4 | match |
| 2c41edd1c487b60d95a105f7befc1dae7551ecfc95093b19f744fdf37d61ce83 | koblitz vars=8 | unit | 0 / 0 | 2 / 2 | [2] / [2] | 2 | match |
| 342411005bf13bcb0de11e4a91052ab59279b217aad99d76f9c016e842d5dd57 | koblitz vars=8 | unit | 0 / 0 | 3 / 3 | [2, 3] / [2, 3] | 3 | match |
| afb67182ba8e29a6acdeb415232c8762fee68e558a81ff4fdde5ff391a9d3dc8 | binary vars=8 | proper-local | 2 / 2 | 2 / 2 | [2] / [2] | 2 | match |
| 1faf7a191a7fd89136541c2f39f766b6c64d266afa66ba77450fdfaad73634c1 | binary vars=8 | proper-local | 4 / 4 | 3 / 3 | [2, 3] / [2, 3] | 3 | match |
| f6c46fea2d5237dbbeca95ce45975abbb3c0cec6603c048c03e043cf159c8096 | binary vars=8 | proper-local | 13 / 13 | 4 / 4 | [3, 4] / [3, 4] | 4 | match |

### 3.4 Full populations (not a sample)

Three complete runs of all 1440 systems: `all_*.jsonl` (first engine version), `final_*.jsonl` (after a memory-only patch of the elimination merge) and `rerun_*.jsonl` (final engine file, sha256 `75b34ce1...`). LFD, falls and zero-set sizes are identical across the three; per-degree dimensions are identical between `final` and `rerun`. The table is built from `final_*.jsonl` via `all_compare.csv` (which also contains the bonus runs, group column `bonus`).

| population | systems | mismatches (LFD / zero set / status / falls) | LFD histogram (mine = repo) |
|---|---|---|---|
| prime k=2 (4 variables) | 220 | 0 | LFD 3: 3, LFD 4: 217 |
| prime k=3 (6 variables) | 235 | 0 | LFD 4: 235 |
| `binary` kind, 4 / 6 / 8 variables | 158 / 179 / 198 | 0 | see per-size below |
| `koblitz` kind, 4 / 6 / 8 variables | 128 / 152 / 170 | 0 | |
| Boolean, 4 variables (286) | | 0 | 0: 41, 1: 10, 2: 224, 3: 11 |
| Boolean, 6 variables (331) | | 0 | 0: 19, 1: 1, 2: 220, 3: 91 |
| Boolean, 8 variables (368) | | 0 | 0: 3, 2: 117, 3: 239, 4: 9 |
| total | 1440 | **0** | |

Also compared for all 1440: the repository's `quotient_ranks` list equals (my ordinary-ring dim V_D) minus (number of non-squarefree monomials of degree <= D),
for every scanned D and with the same list length (5458 entries compared, 0 differences). That is the same identity as in section 5, checked on every system.

**Mismatches: none. There is no system hash to report.**

## 4. Reading of the numbers

* Prime k = 2 and 3: the zero set has 0, 2, or 4 points; LFD is 4 for all k = 3 systems and for all but three k = 2 systems (LFD 3).
  The fall is a single fall at D = LFD: `dim V_3 = 49` against `dim I_3 = 80..84` for k = 3.
* Boolean: falls lists such as [2, 3, 4] (koblitz, 8 variables, 14 points) reproduce exactly, which is a non-trivial agreement.

## 5. Item 3: relating the two dimension counts for k = 3

Identity used (all pieces computed by my own code, `item3.py`):

    dim V_D (ordinary)  =  dim( C cap P_{<=D} )  +  dim NF(V_D),

where C is the chain ideal and NF is reduction by u_{i,j}^2 -> u_{i,j+1}, u_{i,k-1}^2 -> 1, a linear map that sends each monomial to a squarefree monomial
without raising its degree. The chain equations are a Groebner basis (pure powers of distinct variables), so ker NF = C and
`dim(C cap P_{<=D}) = |M_{<=D}| - #squarefree(<=D)` ("nonstandard count"). I verified, for each system and D = 1..4, that (a) the rank of NF on `P_{<=D}` equals
`#squarefree`, (b) the closure of the chain equations alone, `V^chain_D`, has exactly this dimension, so `C cap P_{<=D}` lies in `V_D`
(the chain equations are inputs), and (c) `dim V_D = dim ker NF + rank(NF applied to V_D basis)`. "Modulo the chains" means `q_D = dim NF(V_D) = dim V_D - dim(C cap P_{<=D})`.
Repository quantity: `quotient_ranks[D]` (sweep) and `dim_V3`, `dim_V4` (growth.json). Two systems in the ordinary ring for k = 3 (r = 6, |M_{<=3}| = 84, squarefree 42, nonstandard 42):


System `2dceabe91fba9ccf94619c5384e6c83298005182de98bd28e84c56e9590b476c` (repo `quotient_ranks` [0, 0, 1, 7, 55]):

| D | monomials M<=D | squarefree | nonstandard | dim ker NF = dim(C cap P<=D) | dim V^chain_D | dim V_D ordinary | dim NF(V_D) (modulo chains) | repo quotient rank | dim I_D ordinary |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 7 | 7 | 0 | 0 | 0 | 0 | 0 | 0 | 5 |
| 2 | 28 | 22 | 6 | 6 | 6 | 7 | 1 | 1 | 26 |
| 3 | 84 | 42 | 42 | 42 | 42 | 49 | 7 | 7 | 82 |
| 4 | 210 | 57 | 153 | 153 | 153 | 208 | 55 | 55 | 208 |

System `3738ef21f81c17605cbc17f87ba3d2ea2993679ebaeb8f1842555c816a8f7bb9` (repo `quotient_ranks` [0, 0, 1, 7, 53]):

| D | monomials M<=D | squarefree | nonstandard | dim ker NF = dim(C cap P<=D) | dim V^chain_D | dim V_D ordinary | dim NF(V_D) (modulo chains) | repo quotient rank | dim I_D ordinary |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 7 | 7 | 0 | 0 | 0 | 0 | 0 | 0 | 3 |
| 2 | 28 | 22 | 6 | 6 | 6 | 7 | 1 | 1 | 24 |
| 3 | 84 | 42 | 42 | 42 | 42 | 49 | 7 | 7 | 80 |
| 4 | 210 | 57 | 153 | 153 | 153 | 206 | 53 | 53 | 206 |

System `44eb174d4df52648945a3bfa0e444213345e11550c3f356e7d636e1aac9c33ca` (repo `quotient_ranks` [0, 0, 1, 7, 57]):

| D | monomials M<=D | squarefree | nonstandard | dim ker NF = dim(C cap P<=D) | dim V^chain_D | dim V_D ordinary | dim NF(V_D) (modulo chains) | repo quotient rank | dim I_D ordinary |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 7 | 7 | 0 | 0 | 0 | 0 | 0 | 0 | 7 |
| 2 | 28 | 22 | 6 | 6 | 6 | 7 | 1 | 1 | 28 |
| 3 | 84 | 42 | 42 | 42 | 42 | 49 | 7 | 7 | 84 |
| 4 | 210 | 57 | 153 | 153 | 153 | 210 | 57 | 57 | 210 |

System `b34d90542a9caadee3eb7864d288f54ef159289d39727c55c09325073bf463a9` (repo `quotient_ranks` [0, 0, 1, 7, 57]):

| D | monomials M<=D | squarefree | nonstandard | dim ker NF = dim(C cap P<=D) | dim V^chain_D | dim V_D ordinary | dim NF(V_D) (modulo chains) | repo quotient rank | dim I_D ordinary |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 7 | 7 | 0 | 0 | 0 | 0 | 0 | 0 | 7 |
| 2 | 28 | 22 | 6 | 6 | 6 | 7 | 1 | 1 | 28 |
| 3 | 84 | 42 | 42 | 42 | 42 | 49 | 7 | 7 | 84 |
| 4 | 210 | 57 | 153 | 153 | 153 | 210 | 57 | 57 | 210 |

Result: dim V_3 (ordinary ring) = 42 + 7 = 49 for both proper and both unit systems; modulo the chains it is 7 = 2k + 1, equal to the repository's
`quotient_ranks[3]` and to `growth.json` `dim_V3`. Reason for the count: V_3 modulo the chains is spanned by G (degree 2) and the six products u_{i,j} G, and
no cancellation occurs. In the ordinary ring the same space is bigger by the 42 chain consequences of degree <= 3.
At D = 4 the proper systems give `dim V_4 = dim I_4` (208 for 2 points, 206 for 4 points), i.e. modulo the chains 55 and 53 = 57 - |Z|; the unit systems give 210 (= everything, 1 in V_4).

`growth.json` reproduction (all 30 rows, curves (257,2,3), (7681,1,0), (7681,5,7), (12289,1,0), (12289,3,11), k = 3..8, generator `pkm_generator`):
`dim_V3` = 2k+1 in 30/30 rows (my NF-image dimension in the ordinary ring), `solutions` equal in 30/30, and `dim_V4` equal in 15/15 rows with k <= 5
(k = 3: 57, 57, 57, 57, 57 for the five curves; k = 4: 161, 163, 163, 163, 163; k = 5: 384, 386, 386, 386, 386). Rows k = 6..8 at D = 4 and `deficiency` were not recomputed.
The relation `dim V = dim ker NF + dim NF(V)` and `V^chain = C cap P_{<=D}` held in all 45 profiled (system, D) rows.

## 6. Bonus checks beyond the dispatch (k >= 4)

These are outside the dispatch (k in {2,3}); I ran them because the repository's non-trivial claims (LFD 5, falls [4, 5]) only occur for k >= 6.

* Prime k = 4 (8 variables, D <= 5 scanned): 6 systems stratified over primes and proper/unit (hashes in `bonus_k4.jsonl`): LFD 4, falls [4], zero-set size equal, 6/6 equal to the repository, under 1 s each.
* Prime k = 5 (10 variables): 4 systems (`bonus_k5.jsonl`): LFD 4, falls [4], 4/4 equal, about 4 s each.
* Prime k = 6 (12 variables, ordinary ring with 6188 monomials at D = 5, 893 MB peak). One proper and one unit system from each of the two LFD classes the repository reports at k = 6 (LFD 5 with falls [4, 5], and the minority LFD 4). 4/4 equal:

| system_hash | prime | status | zero set | LFD mine / repo | falls mine / repo | certificate D0 | repo `quotient_ranks` vs mine | time / process peak RSS (cumulative high-water mark) |
|---|---|---|---|---|---|---|---|---|
| 00ba644009d0398fb83ec46158f8c4f4a4bf75cab96579494d4a4d60f6d6fc67 | p=641 | proper-local | 6 | 5 / 5 | [4, 5] / [4, 5] | 5 | equal | 415 s / 867 MB |
| 09345f6decff8608f4575814c0a2595094aa19181046c9045f7f7a2c8521fc95 | p=7681 | unit | 0 | 5 / 5 | [4, 5] / [4, 5] | 5 | equal | 406 s / 893 MB |
| cebac5121b4a52884f40bd8c3309b5b7a8401a5354843687d65403475b4b61fa | p=257 | proper-local | 64 | 4 / 4 | [4] / [4] | 4 | equal | 19 s / 893 MB |
| 060d360f192a92ee12bd59d9ad26e7b15f39a13c534ca6d4d253a13bfff60948 | p=257 | unit | 0 | 4 / 4 | [4] / [4] | 4 | equal | 19 s / 893 MB |

Per-degree data for the first k = 6 system (proper, 6 points, p = 641): D = 3: dim V = 169, dim I = 449; D = 4: dim V = 1631 < dim I = 1814 (fall at 4: dim(V_4 cap P_{<=3}) = 267 > dim V_3 = 169); D = 5: dim V = 6182 = dim I = |M_{<=5}| - 6 = 6188 - 6 with evaluation of degree <= 4 monomials of rank 6 (fall at 5: dim(V_5 cap P_{<=4}) = 1814 > dim V_4 = 1631). Hence LFD = 5 exactly, as in the repository (`quotient_ranks` [0, 0, 1, 13, 605, 1580] equal to my ordinary dimensions minus the nonstandard counts: 1631 - 1026 = 605, 6182 - 4602 = 1580).


## 7. Requirement table

| # | requirement | status | evidence |
|---|---|---|---|
| 1 | >= 12 prime systems, k in {2,3}, recompute exact LFD in the ordinary ring, compare by system_hash | verified complete: 455 of 455 (17 + 18 in the stratified tables) | `final_prime2.jsonl`, `final_prime3.jsonl`, `all_compare.csv`, `out_prime{2,3}.jsonl` |
| 2 | >= 12 binary/Koblitz systems with 2k <= 8 | verified complete: 985 of 985 (46 in the stratified tables) | `final_bin{4,6,8}.jsonl`, `out_bin{4,6,8}.jsonl` |
| 3 | dim V_3 consistency for k = 3 and the explanation of how the counts relate | verified complete: 4 sweep systems plus 30 growth rows | `item3_sweep.json`, `item3_growth.json`, `log_item3_growth.txt` |
| 4 | no closure engine of the repository used | verified: only `PrimeCurve`, `pkm_generator` (growth mode of item 3) | `lfd_engine.py`, `run_checks.py` |
| 5 | memory < 2 GB, each computation < 15 min | verified: max RSS 79 MB (required), 867 MB (largest bonus); longest required run 0.7 s, longest bonus run 415 s | run logs |
| 6 | report every mismatch with hash | none to report | `all_compare.csv`, column `match_all` all True |
| 7 | certificate discipline | each LFD <= D0 claim carries the explicit certificate of section 2, re-verified directly for the 534 proper systems and via `1 in V` for the 906 unit ones | `soundness_check.py`, `log_soundness.txt` |
| 8 | controls: null / cheaper baseline / reference arm | not applicable to an exact-algebra recomputation; the analogues are the self-tests, the literal-rule second implementation and the sensitivity control (section 2) | |
| 9 | reporting contract, matched rho (section 5a band), validation gate (`VALIDATION.json`), m = 83 / ECC2K-130 disclosure | not applicable: no cost, speedup, rho comparison or ECC2K-130 claim is made or tested here; this report does not supply a validator report for the IC leaderboard or `benchmarks/README.md` | |

## 8. Limits and things I did not do

* The 1440 systems are exactly the k in {2,3} prime and 2k <= 8 Boolean systems of the sweep. For prime k = 4..7 and Boolean 2k > 8 the repository's
  claims (including the 312 LFD-5 results at k = 6, 7 (177 + 135) and the Boolean 10..14-variable results) are covered only by the 14-system bonus sample of section 6 (k = 4, 5, 6) and not at all for k = 7 or Boolean 2k > 8.
* The recomputation relies on the stored generators/equations (cross-checked against S_3 with my own formulas and the case metadata, section 2 item 7).
* Exact LFD below the certificate degree uses `I_D = ker(evaluation)`, which uses radicality (argued in section 2, not machine-checked here). The upper bound `LFD <= D0` does not.
* The unit-ideal certificate `1 in V_D` is tested by my closure code; its correctness rests on the cross-implementation checks (items 3 and 4 of section 2), not on an explicit cofactor expression for 1.
* The repository's Lean file `Closure.lean` was used only as the definition; no Lean was run.

## 9. Code and logs (all in the scratchpad directory)

`/tmp/claude-0/-home-user/a6e79f2a-a5ec-5904-b377-343a14909fff/scratchpad/`

| file | role | sha256 |
|---|---|---|
| `lfd_engine.py` | engine: monomial spaces, exact echelon, V_D closure (variable-wise and literal), eval/rank, `exact_lfd`, system builders | `75b34ce1e31d7b677409a3bc153691de66c97cf4125dfc646d8b7b7b127510ea` |
| `selftest.py` | self-tests (rank, toys, literal-rule equality) | `0589f69d8d92977d233bb88e293bfcfb423de3aae5f0279568ef9557d08d56e3` |
| `run_checks.py` | driver: selectors `allprime2`, `allprime3`, `allbin4/6/8`, `prime2/3`, `bin4/6/8`, `sample:k:n`, `hash:<h>`; integrity checks | `cf420057b9990c96740671fb58bd1b30ce7858cee2b42fea756f5494c5722d7c` (changed once after the first runs, adding the `all*` and `sample:` selectors and the `LFD_DMAX` override; the final audit and rerun used this file) |
| `item3.py` | NF, chain-only closure, dimension identity, growth.json reproduction | `6d24033065032c687732734da34431a450b982178019f1a0e4e6a503526f8349` |
| `naive_check.py`, `naive_lfd.py` | literal-rule pure-Python cross-checks | `89431abd...`, `a0439bfe...` |
| `soundness_check.py` | vanishing / dimension / surjectivity re-check at the certificate | `2bcba8fea3aae8f2d15f2a02ea3c15c1566814f444924887c13e7208420d4e10` |
| `control.py` | sensitivity control | `926229f3...` |
| `summarize.py` | builds `all_compare.csv` | `923bca6a...` |

Reproduce (from the scratchpad directory, no repository writes):

    export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 PYTHONDONTWRITEBYTECODE=1
    python3 -I selftest.py
    python3 -I run_checks.py out.jsonl allprime3 allprime2 allbin4 allbin6 allbin8        # 1440 systems, about 2 minutes
    python3 -I soundness_check.py                                                         # needs final_*.jsonl
    python3 -I item3.py sweep 2dceabe91fba 3738ef21f81c 44eb174d4df5 b34d90542a9c
    CSD_CACHE=<scratch> python3 -I item3.py growth 8 5
    LFD_DMAX=5 python3 -I run_checks.py bonus_k6.jsonl $(cat k6_picks.txt)               # about 7 minutes per LFD-5 system, 870 MB
