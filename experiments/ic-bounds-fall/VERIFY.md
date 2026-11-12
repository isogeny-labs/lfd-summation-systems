# Independent verification of item 5, "Fall-space structure"

Verifier report for the producer result in this directory (`README.md`, `p40961-s7/`,
`csd/ic_bounds_fall.py`, `scripts/ic_bounds_fall_{export,analyze}.py`). Date: 2026-11-08.

## How the check was done

I wrote a separate engine, `scripts/ic_bounds_fall_verify.py`. It works from the
definitions only and does not import the producer's module or scripts, or any closure engine
in the repository. It reuses one thing: the dense mod-p row-echelon class `Ech` from the
independent validator `experiments/ic-curve-validation/code/lfd_engine.py`. No closure logic
comes from that file.

- **Model.** Modulo the chains, a monomial is an exponent pair in `(Z/2^k)^2` and its degree
  is the Hamming weight. The variable `u_{i,j}` acts as translation by `2^j` in block i,
  mod 2^k. The generator `g = S_3(x1, x2, t)` is rebuilt from `(a, b, t)` with the textbook
  formula `S_3 = (x1-x2)^2 t^2 - 2((x1+x2)(x1x2+a)+2b) t + (x1x2-a)^2 - 4b(x1+x2)`. It matches
  the producer's `instance.json` generator coefficient by coefficient
  (`verify/closures.json: generator_equals_producer = true`). The target is on the curve.
- **Closure.** The closure is breadth-first. Round 0 is the generator rows `tau_mu g`,
  wt mu <= D-2. Round r+1 adds every translate of `V^(r) ∩ M_<=D-1`. The RREF sorts columns
  by weight, highest first, so the rows whose pivot has weight <= d span `V ∩ M_<=d`.
  Closedness is checked separately at the end: every translate of every row of weight < D
  is re-reduced against the basis, and all of them reduce to zero.
- **Free extension and certificate.** Built from my own `V_4(k0)`. The certificate checks
  three things against the spanning set alone, without using `V_4(k)`: it contains all
  generator rows, it is closed under every translate of its weight-<=3 part, and it does not
  contain 1.
- **Everything else.** Level pieces come from coordinate-subspace intersections plus Möbius
  inversion over level sets. As a stronger check, I also computed the graded dimensions
  `dim(F∩M_L) - dim Σ_{L'⊊L} F∩M_{L'}`. Multipliers come from my own 2-D DFT over F_p, with
  a round-trip check. The kappa values are computed as `dim(B_w ∩ M_<=d·g^{-1})`. The
  symbolic kernel uses sympy, as an independent 7x9 system.

Resources: every run used `ulimit -v 3000000` and one thread. The longest run was the
degree-5 closure at k=8: 490 s, about 0.9 GB. New files total about 60 kB. I did not touch
`~/.cache/ic-chain-depth10`, and I ran no Lean builds.

## Commands

    ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    python3 scripts/ic_bounds_fall_verify.py closures --k 6 7 8 9 10 11 12     # verify/closures.json  (~5 min)
    python3 scripts/ic_bounds_fall_verify.py free --k0 7 6 --k 7 8 9 10 --top 6 7   # verify/free-p40961.json (~16 min)
    python3 scripts/ic_bounds_fall_verify.py free --k0 7 --k 11 --out experiments/ic-bounds-fall/verify/k11   # k=11. The result recorded in verify/k11/free-k11.txt is a log excerpt from a combined '--k 11 12' run. That run passed k=11, then hit the 3 GB cap at k=12, so k=12 was re-run on its own.
    python3 scripts/ic_bounds_fall_verify.py free --k0 7 --k 12 --out experiments/ic-bounds-fall/verify/k12   # verify/k12/free-p40961.json
    python3 scripts/ic_bounds_fall_verify.py levels --k 6 7 8 --kstruct 8 9 10   # verify/levels.json
    python3 scripts/ic_bounds_fall_verify.py w --k 6 7                          # verify/w.json
    python3 scripts/ic_bounds_fall_verify.py cross --k 6 7 8 9 --kfree 8 9      # verify/cross.json
    python3 scripts/ic_bounds_fall_verify.py mult --k 6 7                       # verify/mult.json (unit instance p=1073153)
    python3 scripts/ic_bounds_fall_verify.py d5 --k 6 7 8                       # verify/d5-k{6,7,8}.json
    python3 scripts/ic_bounds_fall_verify_npz.py --k 6 7 8 9 10 --k5 6          # verify/npz.json (producer bases vs own)

Instances, all with the target x-coordinate as t:

- p40961: a=21223, b=9887, target (4748, 524).
- p12289r3: a=3899, b=9710, target (9895, 11610).
- p12289j1728: a=1, b=0, target (9362, 1126). This is identical to
  `experiments/ic-chain-depth10/probe.json`.
- p1073153: a=281783, b=132345, target (247294, 569236).

The subgroup generator is `z = c^((p-1)/2^k)` for the least c >= 2 with `z^(2^(k-1)) ≠ 1`.
The subspaces do not depend on this choice. Only the solution coordinates do.

## Results, claim by claim

| # | claim (short) | producer label | verdict |
|---|---|---|---|
| 1 | dims of V_4, V_4∩M_<=3, V_4∩M_<=2, 1 ∉ V_4, closed, k=6..12 | exact | **confirmed** |
| 2 | V_4(k) = W_k, certificate, stable low part, k=7..12 | exact | **confirmed** |
| 3 | reduction theorem (iii) and its consequences | proof | **confirmed** (proof checked; hypothesis re-verified) |
| 4 | depth-6 base fails (42 each way, 563, ranks 37 and 4) | exact | **confirmed** |
| 5 | level-set counts and F_4(k) = F_4(7) ⊕ translates | exact | **confirmed** |
| 6 | w = Q_w g, symbolic kernel of dimension 2, instance check | proof | **confirmed** |
| 7 | cross-instance dims, pivot set, free extension | exact | **confirmed** |
| 8 | multiplier filtration and kappa_4, kappa_3 on the unit instance | exact | **confirmed** |
| 9 | degree-5 dims, local-to-global certificate, profiles | exact | **confirmed** (round count uses a different convention, see below) |
| 10 | 9 rounds at every depth; low birth histograms | exact | **confirmed**; the producer's reproduce one-liner fails |
| 11 | degree-5 free-extension prediction | conjecture | unverifiable (label is appropriate) |
| 12 | no closed form for V_4(7); pivot set instance-independent | empirical | **confirmed** for the pivot-set part; "not found" cannot be checked |

### Claim 1 (exact): confirmed

`verify/closures.json`, own engine:

| k | dim V_4 | profile (d=0..4) | 4k^2+172k-571 | 1 ∈ V_4 | closed (independent re-check) | solutions |Z| |
|---|---|---|---|---|---|---|
| 6 | 605 | 0,0,2,111,605 | 605 | no | yes | 2 |
| 7 | 829 | 0,0,2,115,829 | 829 | no | yes | 2 |
| 8 | 1061 | 0,0,2,119,1061 | 1061 | no | yes | 2 |
| 9 | 1301 | 0,0,2,123,1301 | 1301 | no | yes | 6 |
| 10 | 1549 | 0,0,2,127,1549 | 1549 | no | yes | 26 |
| 11 | 1805 | 0,0,2,131,1805 | 1805 | no | yes | – |
| 12 | 2069 | 0,0,2,135,2069 | 2069 | no | yes | – |

The rank histories are identical to the producer's, entry by entry, for k=6..10.

I also compared the producer's exported bases `k{6..10}-D4.npz` with my bases. They have the
same row spaces at every k (`verify/npz.json`), and each file's sha256 matches the hash
recorded in its `k*-D4.json`. I re-ran the producer's own `ic_bounds_fall_export.py` for
k=6,7,8 into a scratch directory. It reproduces the recorded `basis_sha256` and file sha256
exactly.

Evidence: exact. The basis is a certified closed subspace, and 1 ∉ V_4 is a reduction of the
constant against that closed basis.

### Claim 2 (exact): confirmed

From `verify/free-p40961.json` (k=8,9,10), `verify/k11/free-k11.txt` and
`verify/k12/free-p40961.json`, at every k=8..12:

- W_k, built from my own V_4(7), has exactly 829 + 230(k-7) + 2C(2k-14,2) rows.
- `rank(W_k) = rank(V_4(k)) = rank(W_k + V_4(k))`, so W_k = V_4(k).
- The certificate, computed on W_k alone, finds that W_k contains all generator rows, is
  closed under all 2k translates of its weight-<=3 part, and does not contain 1.
- `V_4(k) ∩ M_{levels<=6}` has dimension 829 and equals V_4(7) embedded, so it lies in the
  base.

k=7 is trivial: W_7 = V_4(7).

### Claim 3 (proof): confirmed

I read the proof in README §4.2 line by line. It holds:

- **Direct sum.** The decomposition by high pattern is a direct sum: low exponents are below
  128, and h ≡ 0 mod 128 blockwise.
- **(i)** uses no wrap-around: `mu_low + s <= 96 + 2 < 128`. It also uses wt(a+b) <= wt a + wt b.
- **(iii), case j >= 7.** Adding one bit to h raises the weight by at most 1, even with
  carries or a wrap at the top.
- **(iii), case j <= 6.** There is exactly one possible carry, into bit 7. The carried
  monomials lose at least 7-j >= 1 bits. Closedness of W_8 in A_8, together with uniqueness
  of the h-component in A_8, gives `v_0 ∈ S ∩ M_<=4` and `v_1 ∈ S ∩ M_<=3`. The weight
  budgets `wt v_0 <= 4 - wt h` and `wt v_1 + wt(h + 128 e_i) <= 3` both follow.
- **Weight-<=3 elements.** An element of W_k of weight <= 3 decomposes into components of
  weight <= 3, because the supports are disjoint and weights add. So the reduction to single
  `tau_h v` is valid.
- **Base not needed closed.** The proof never needs S to be closed in A_7, which is correct:
  S is read inside A_k without wrap.

The hypothesis of (iii), that W_8 is closed in A_8, is re-verified by my certificate at k=8.

Two scope notes:

- "Instance-free" applies to the reduction only. The conclusion 1 ∉ V_4(k) and the bound on
  dim V_4(k) hold for an instance whose own depth-8 certificate passes. The numbers 829,
  230 and 2 are that instance's profile of S.
- This is a written proof, not a kernel proof. Nothing here is Lean-checked.

### Claim 4 (exact): confirmed

From `verify/free-p40961.json`, for k=7,8,9 (and k=10 as well):

- Free_k^(6)(V_4(6)) has the predicted dimension.
- `dim(W + V_4) - dim V_4 = dim(W + V_4) - dim W = 42`.
- The certificate finds the space is not closed, although it does contain the generator rows
  and does not contain 1.
- `dim V_4(k) ∩ M_{levels<=5} = 563`, and that intersection lies in V_4(6).

The fall space of V_4(6) has rank 37 on the bit-5 columns. The fall space of V_4(7) has
rank 4 on the bit-6 columns.

### Claim 5 (exact): confirmed

Level-set analysis (`verify/levels.json`), k=6,7,8:

- The Möbius counts equal the producer's table key for key (`equals_producer_counts = true`).
  None are negative, and they sum to 4k+87.
- The graded dimensions agree with the Möbius counts and also sum to dim F. This is the
  direct-sum certificate.
- The span of the level pieces is the whole of F.

Caveat: at k=6 the 33-dimensional piece {0..5} is the full level set, so that piece holds
by definition.

Fall-space structure: for k=8,9,10, `F_4(7) ⊕ span(tau_{2^j e_i}(V_4(7) ∩ M_<=2))` for
j=7..k-1 equals F_4(k) as row spaces, with ranks 119, 123 and 127.

### Claim 6 (proof): confirmed

`verify/w.json`, own sympy computation:

- The weight-≥3 exponent pairs in {0..4}^2 are exactly the 7 listed.
- The 7x9 system over Q(a,b,t) has a kernel of dimension 2.
- The producer's vector Q_w lies in it. Together with delta_0, it spans the kernel
  (rank 2).
- On the instance (a+3t^2 = 25124 mod p), w has weight 2 and 17 monomials, and it lies in
  V_4. At k=6 and k=7, span(g, w) = V_4 ∩ M_<=2.

The proof is complete as a statement about Q(a,b,t), that is, generic (a,b,t) with
a+3t^2 ≠ 0. That the kernel stays 2-dimensional under specialization is only checked on the
instance.

### Claim 7 (exact): confirmed

`verify/cross.json`:

- p12289r3 and p12289j1728 both have dim V_4 = 605, 829, 1061, 1301 at k=6..9, with the same
  profiles as the primary instance.
- Neither contains 1, and all closures are closed.
- The pivot set of V_4(7) in my column order is identical on all three instances
  (symmetric difference 0).
- At k=8 and k=9, W_k = V_4(k) and the certificate passes on both instances.

The unit-ideal instance p1073153 has dims 605 and 829 at k=6,7 (`verify/mult.json`). g has no
zero on H^2 there, so 0 solutions. j1728 also has 0 solutions at k=6..8. For p12289r3 the
target order is 6151 (own point count: #E = 12302 = 2·6151).

### Claim 8 (exact): confirmed

`verify/mult.json`, own 2-D DFT, with a round-trip check `(v g^{-1}) g = v`:

- **k=6.** `dim(Q∩B_w)` = 1,13,79,123,212,280,377,465,546,588,605. kappa_4 =
  1,13,79,123,212,280,377,465,583,717,781,793,794. The fall filtration is
  1,13,18,34,43,59,72,88,103,111. kappa_3 = …,103,222,….
- **k=7.** The filtration is 1,15,106,160,282,368,497,611,724,796,829, and kappa_4(w<=6)
  equals it.

The maximum multiplier weight is 10 at both depths. All values match the producer's numbers.

### Claim 9 (exact): confirmed, with a counting-convention note

`verify/d5-k{6,7,8}.json`, own degree-5 closures:

- **Dimensions.** dim V_5 = 1584, 3471 and 6883. Each equals `|M_<=5| - |Z|`, with |Z| = 2.
- **Certificate.** Evaluation of M_<=4 onto Z has rank 2, so the local-to-global certificate
  holds and V_5 = I ∩ M_<=5.
- **Profiles.** The low profiles equal `|M_<=d| - rank eval` and match the claim.
- **Cross-check.** My V_5 ∩ M_<=4 at k=6 has the same row space as `k6-D5-low.npz`.

The rank histories match the producer's entry by entry. The stated "9, 11, 13 rounds",
however, counts a final round that adds nothing: each D=5 history ends with a repeated value.
By the same convention as the D=4 "9 rounds", where the last round adds 7 rows, the
productive rounds are 8, 10 and 12. This is cosmetic, but the round counts at D=4 and D=5 are
not comparable as stated.

The profiles are forced by |Z| = 2 once V_5 = I ∩ M_<=5, so they carry no closure structure.
The producer's README says this itself.

### Claim 10 (exact): data confirmed, reproduce command broken

My breadth-first closure reproduces the producer's rank histories for k=6..10. It needs 9
productive rounds at every k=6..12. The low birth histograms are
`[6+2k, 4+2k, 9, 16, 12, 13, 16, 9, 2, 0]`, identical across k except rounds 0 and 1, which
grow by 2 per depth.

The claim's reproduce one-liner (`cat k*-D4.json | python3 -c ...split('}\n{')...`) fails with
`JSONDecodeError: Extra data`. The JSON files have no trailing newline, so `cat` joins them
as `}{` and the split never happens. The command is broken; the data are not.

### Claim 11 (conjecture): unverifiable; the label is appropriate

The claim is explicitly a prediction, and it needs the heavy runs listed. One wording
problem: "1 ∈ V_5(k) for k <= 9" refers to the unit-ideal (probe-family) instances, not to
the primary instance, where 1 ∉ V_5 trivially.

### Claim 12 (empirical): pivot-set part confirmed; "not found" cannot be checked

The pivot-set coincidence on three instances is confirmed under a different column order.
It is correctly labelled empirical: three instances are not a genericity theorem.

## Other findings (not in the claims list)

- **Wrong primality statement in README §2.** The table says the target order is
  "40817 (= curve order, prime)". My own point count gives #E = 40817 = 7^4 · 17, which is
  not prime. The target order is indeed 40817, so the target generates the group. The
  "prime" is a factual error, harmless for the fall-space results.
- **"10 = 2k-2" holds only at k=6.** README §5.4 writes the needed multiplier weight as
  "10 = 2k-2". At k=7 the weight is also 10, not 2k-2 = 12, so the identity holds only at
  k=6.
- **README §9 commands run.** I ran `ic_bounds_fall_export.py --k 6 7 8` and
  `ic_bounds_fall_analyze.py --kmax 8 --moebius-kmax 7 --unit-prime 0 --extra-instances ''`
  with `--out` pointing at a scratch directory, so that no producer file was overwritten.
  Both completed. They reproduced the stored hashes, the free-extension results and the w
  checks. I did not re-run the full analysis at `--kmax 10` with the unit and extra
  instances; my own engine covered the same keys.
- **No label is upgraded.** No stall is treated as a bound: every dimension comes from a
  closure verified closed, and every "1 ∉" is checked against a closed basis or a passing
  exclusion certificate. The fits (`4k^2+172k-571`, `4k+87`) are explained by the proved
  reduction plus the computed base profile. For k <= 12 they are exact computations; for
  k > 12 they follow from the reduction theorem on this instance. They are not presented as
  a law beyond that.

## Heavy runs not done

None of these is needed for the verdicts above.

- The degree-5 conjecture (claim 11) needs V_5(11). The producer's estimate is about
  32 GB and 24 h, which exceeds the limits here.
- A full genericity sweep of the depth-7 base is the producer's third heavy run.
