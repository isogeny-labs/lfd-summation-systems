# Free-extension structure of the two-summand closures (D = 3, 4, 5)

Status 2026-11-10. Labels: **kernel** (Lean-checked), **proof** (written), **exact** (certified
computation, command and record named), **empirical**, **conjecture**.

## What is checked

For a subspace S of the depth-k0 chain algebra (S ⊆ M_≤D), the free high extension is
`W_{k0+1}(S) = ⊕_q τ_q(S ∩ M_{≤D−wt q})`, q ∈ {0,1}² the pattern of the new level k0
(`docs/proofs/N5-reduction.md` §5, `DepthReduction.W`). The kernel theorem
`DepthReduction.W_closed` / `ChainModelMap.not_one_of_base` / `lfd_gt_of_base` (any k0, any D)
needs three base facts:

1. `g ∈ S`;
2. `W_{k0+1}(S)` closed at depth k0+1 (every u_{i,j}, j ≤ k0, applied to every element of weight < D lands in W);
3. `1 ∉ S`;

and then gives `1 ∉ V_D(k)` and `V_D(k) ⊆ W_k(S)` for every k ≥ k0+1.

For S = V_D(k0) (the closure of `csd/ic_bounds_d5_verify.py`) each record checks, at every
computed k0:

* fact 2 by the **deterministic full sweep** of that engine run on W (every block-1 shift of
  every basis row of weight < D, both swap parts, must reduce to zero; with the swap split
  closedness under u_{0,j} gives closedness under u_{1,j} = σu_{0,j}σ since W is σ-invariant;
  D = 3 and D = 4 are repeated without the split, same verdicts), plus a random sweep;
* fact 1 as `g` reducing to zero against V_D(k0), and additionally all 1794 (D=5, k0+1=11)
  generator rows `NF(m g)`, wt m ≤ D−2, against W;
* fact 3 as "no pivot at the constant column";
* `W_{k0+1}(V_D(k0)) = V_D(k0+1)`: canonical RREFs identical (`same_space`), or, when W was
  freed before computing V_D(k0+1) (`--drop-w`), equal SHA-256 of the canonical RREFs; at the
  top depth against the recorded V_D(k0+1) of `experiments/ic-bounds-d5/verify/` (same engine
  source hash `c3d123ed…`; record path and file hash are stored);
* the profile `p_j = dim(V ∩ M_≤j)`, the recursion `p_j(k+1) = p_j(k) + 2p_{j−1}(k) + p_{j−2}(k)`
  (this is exactly the profile of W), and `dim(I_k ∩ M_≤D) = |M_≤D| − rank(eval of M_≤D on Z)`
  for depths with 2^k | p−1 ("deficient" = V_D(k) ≠ I_k ∩ M_≤D);
* the shifted-generator witness search: every low pair e with e_i + 2 < 2^k (τ_e g has the same
  support at every depth ≥ k) and all weights of τ_e g ≤ D, in order of decreasing wt e,
  tested for τ_e g ∉ V_D(k). The search is exhaustive over that set.

W is assembled by `free_extension_direct` (component RREFs placed side by side; the mixed
components q = (1,0),(0,1) are the RREF of S ∩ M_≤D−1 in full coordinates). It was checked equal
(identical canonical RREF) to the generic elimination `free_extension` on 14 cases (D = 3, 4, 5,
k0 = 2..8, generic, order-two and tangent targets) before the D = 5 runs.

## Commands

```sh
# light (D = 3, 4): ulimit -v 1500000, one thread
R="python3 scripts/ic_bounds_freeext_run.py chain"
$R --p 193 --a 1 --b 0 --seed 10 --D 3 --kmin 2 --kmax 6 --ext-top --witness --name D3-p193-unit-s10
$R --p 257 --a 3 --b 7 --seed 21 --want --D 3 --kmin 2 --kmax 8 --ext-top --witness --name D3-p257-s21
$R --p 12289 --a 4021 --b 9917 --seed 12 --want --draw-k 10 --D 3 --kmin 2 --kmax 8 --ext-top --witness --name D3-p12289-s12
$R --p 193 --a 1 --b 0 --seed 10 --D 3 --kmin 2 --kmax 6 --ext-top --nosplit --name D3-p193-unit-s10-nosplit
$R --p 257 --a 1 --b 0 --seed 10 --D 4 --kmin 4 --kmax 8 --ext-top --witness --name D4-p257-unit-s10
$R --p 7681 --a 3 --b 7 --seed 21 --want --D 4 --kmin 4 --kmax 9 --witness --name D4-p7681-s21
$R --p 12289 --a 4021 --b 9917 --seed 12 --want --draw-k 10 --D 4 --kmin 4 --kmax 9 --witness --name D4-p12289-s12
$R --p 12289 --a 4021 --b 9917 --seed 12 --want --draw-k 10 --D 4 --kmin 5 --kmax 8 --nosplit --name D4-p12289-s12-nosplit
$R --p 12289 --a 1 --b 0 --X 1479 --D 4 --kmin 4 --kmax 9 --witness --name D4-p12289-order2-X1479
$R --p 12289 --a 1 --b 0 --X 64   --D 4 --kmin 4 --kmax 9 --witness --name D4-p12289-tangent-X64
$R --p 12289 --a 1 --b 0 --X 0    --D 4 --kmin 4 --kmax 9 --witness --name D4-p12289-t0-X0
# heavy (D = 5), sequential under the heavy lock:
#   flock /tmp/lfd-locks/heavy.lock sh -c 'ulimit -v 7000000; OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 ...'
R5="python3 scripts/ic_bounds_freeext_run.py chain --D 5 --witness --drop-w --ext-top"
$R5 --p 12289 --a 1 --b 0 --seed 10 --kmin 9 --kmax 11 --ref-name 'k{k}-p12289-s10-D5.json' --name D5-p12289-unit-s10
$R5 --p 12289 --a 4021 --b 9917 --seed 12 --want --draw-k 10 --kmin 9 --kmax 11 --ref-name 'k{k}-p12289-s12-D5.json' --name D5-p12289-s12
$R5 --p 13313 --a 1 --b 0 --seed 11 --want --draw-k 10 --kmin 9 --kmax 11 --name D5-p13313-s11
$R5 --p 15361 --a 4021 --b 9917 --seed 12 --want --draw-k 10 --kmin 9 --kmax 11 --name D5-p15361-s12
python3 scripts/ic_bounds_freeext_run.py chain --D 5 --witness --drop-w --p 40961 --a 1 --b 0 --seed 10 --kmin 9 --kmax 12 --ref-name 'k{k}-p40961-s10-D5.json' --name D5-p40961-unit-s10
python3 scripts/ic_bounds_freeext_table.py      # table.md, summary.json
```

Exceptional abscissas: `python3 scripts/ic_bounds_freeext_run.py exceptional --p 12289 --a 1 --b 0 --kind order2|tangent`.
Records: `records/<name>.json` (command, source SHA-256 of the four sources, ranks, profiles,
basis SHA-256 per swap part, sweep counts, generator-row counts, witness search, timings).

## Results

The outcome of each frozen prediction, with the records that test it, is in
[`outcomes-2026-11-10.json`](outcomes-2026-11-10.json); the free-extension table is [`table.md`](table.md)
and [`summary.json`](summary.json). The results are reported in `paper/curve-lfd-bounds.tex`, Section 7
("The recursion" and "The growth law, as a conjecture").
