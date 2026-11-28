"""Frozen high-delta quadratic-lift pilot against native two-point LRF."""

import hashlib
import json
import random
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from csd import ic_dreg as I
from csd import ic_rootfind as rf

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SEED = 1
PLAN = ((23, 10, 6), (23, 11, 9))
MAX_ATTEMPTS = 5000
EXPECTED_HASHES = {
    "csd/ic_rootfind.py": "727a14d638964fa19c487294a8d281c2cf87d3ef47fecefa6cefe7469fd58cc5",
    "csd/ic_dreg.py": "0e5a650d3717ac029e3cbbc10c8e7f1525c11f46461d651187ce40d89c25afe1",
    "experiments/ic-rootfind-yield/goodset_structure.py": "dc906ad70c54e63b6f324307ffd1cbfb46ad4c4658566a38686329f6eb2cf762",
}


def source_check():
    for relative, expected in EXPECTED_HASHES.items():
        actual = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        if actual != expected:
            raise RuntimeError(f"frozen source changed: {relative}: {actual}")


def affine_data(S, T):
    F, E = S.F, S.E
    xT, yT = T
    if xT == 0:
        return None
    ym = yT ^ xT
    c1 = E.a ^ xT
    c2 = F.mul(F.sq(ym) ^ E.b, F.inv(xT))
    Qf = lambda x: rf._apply(S.Q, x)
    pif = lambda x: rf._apply(S.pi, x)
    top = S.w - 1
    cols = [pif(F.mul(xT, Qf(v) ^ v)) | (F.trace(v) << top)
            for v in S.basis] + [pif(xT)]
    rhs = pif(F.mul(xT, Qf(c1) ^ c1) ^ c2) | (F.trace(c1) << top)
    sol = I._f2_solve(cols, rhs, S.w)
    if sol is None:
        return None
    part, kernel = sol

    def evaluation(mask):
        s = part
        for j, vector in enumerate(kernel):
            if mask >> j & 1:
                s ^= vector
        e1 = 0
        for i, value in enumerate(S.basis):
            if s >> i & 1:
                e1 ^= value
        eps = (s >> S.k) & 1
        lam = Qf(e1 ^ c1) ^ eps
        e2 = F.mul(xT, F.sq(lam)) ^ c2
        assert pif(e2) == 0
        return e1, e2, lam

    base = evaluation(0)
    increments = [tuple(a ^ b for a, b in zip(evaluation(1 << j), base))
                  for j in range(len(kernel))]
    if len(kernel) >= 2:
        check = evaluation(3)
        assert check == tuple(base[i] ^ increments[0][i] ^ increments[1][i]
                              for i in range(3))
    return dict(base=base, increments=increments, evaluation=evaluation,
                xT=xT, ym=ym, kernel_dim=len(kernel))


def exact_pairs(S, T, data):
    F, E = S.F, S.E
    pairs, rows = [], set()
    for gmask in range(1 << data["kernel_dim"]):
        e1, e2, lam = data["evaluation"](gmask)
        if e1 == 0:
            roots = [F.sqrt(e2)]
        else:
            quotient = F.mul(e2, F.inv(F.sq(e1)))
            if F.trace(quotient):
                continue
            first = F.mul(e1, rf._apply(S.Q, quotient))
            roots = [first, first ^ e1]
        for x in roots:
            if rf._apply(S.chk, x):
                continue
            assert (F.sq(x) ^ F.mul(e1, x) ^ e2) == 0
            zmask = x  # The frozen monomial basis is the low-k bit basis.
            assert zmask < 1 << S.k
            pairs.append([gmask, zmask])
            x2 = x ^ e1
            mu = F.mul(lam, data["xT"]) ^ data["ym"]
            P1 = (x, F.mul(lam, x) ^ mu)
            P2 = (x2, F.mul(lam, x2) ^ mu)
            assert rf._apply(S.chk, x2) == 0
            assert E.on_curve(P1) and E.on_curve(P2)
            assert E.add(P1, P2) == T
            rows.add(tuple(sorted((P1, P2))))
    return sorted(pairs), sorted(rows)


def worker_payload(S, data):
    F = S.F
    e1, e2, _ = data["base"]
    return {
        "n": S.n, "k": S.k, "delta": data["kernel_dim"],
        "e1_base": e1, "e2_base": e2,
        "e2_delta": [triple[1] for triple in data["increments"]],
        "basis_square": [F.sq(x) for x in S.basis],
        "e1_base_mul": [F.mul(e1, x) for x in S.basis],
        "e1_delta_mul": [[F.mul(triple[0], x) for x in S.basis]
                         for triple in data["increments"]],
    }


def run_cell(n, k, delta):
    F = I.GF2n(n, (1 << n) | (1 << 5) | 1)
    rng = random.Random(SEED)
    b = rng.randrange(1, 1 << n)
    E = I.BinaryCurve(F, 0, b)
    basis = [1 << i for i in range(k)]
    setup_start = time.perf_counter()
    setup_cpu_start = time.process_time()
    S = rf.Solver(E, basis)
    setup = dict(wall=time.perf_counter() - setup_start,
                 cpu=time.process_time() - setup_cpu_start)
    chosen = {"positive": [], "empty": []}
    attempts = 0
    while attempts < MAX_ATTEMPTS and (len(chosen["positive"]) < 2 or
                                       len(chosen["empty"]) < 2):
        attempts += 1
        R = E.random_point(rng)
        x3 = rng.choice(list(range(1, 1 << k)))
        P3 = E.lift_x(x3)
        if P3 is None:
            continue
        T = E.add(R, E.neg(P3))
        data = affine_data(S, T)
        if data is None or data["kernel_dim"] != delta:
            continue
        pairs, rows = exact_pairs(S, T, data)
        category = "positive" if rows else "empty"
        if len(chosen[category]) >= 2:
            continue
        native_start = time.perf_counter()
        native_cpu_start = time.process_time()
        native_rows, stats, _ = S.lrf(2, T)
        native_wall = time.perf_counter() - native_start
        native_cpu = time.process_time() - native_cpu_start
        assert native_rows == rows == S.d2_reference(T)

        payload = worker_payload(S, data)
        worker_start = time.perf_counter()
        try:
            call = subprocess.run(
                ["sage", "-python", str(HERE / "worker.py")],
                input=json.dumps(payload), text=True, capture_output=True,
                timeout=30, check=False, cwd=ROOT,
            )
            worker_total = time.perf_counter() - worker_start
            if call.returncode:
                solver = dict(status="error", returncode=call.returncode,
                              stderr=call.stderr[-2000:], worker_total_wall=worker_total)
            else:
                solver = json.loads(call.stdout)
                assert solver["pairs"] == pairs
                assert solver["input_degree_max"] <= 2
                solver["status"] = "exact_match"
                solver["worker_total_wall"] = worker_total
        except subprocess.TimeoutExpired:
            solver = dict(status="timeout_30s", worker_total_wall=time.perf_counter()-worker_start,
                          peak_rss_bytes=None)
        chosen[category].append({
            "attempt": attempts, "T": list(T), "pairs": pairs, "rows": rows,
            "native": dict(wall=native_wall, cpu=native_cpu, stats=stats),
            "solver": solver,
        })
        print(f"n={n} k={k} {category} attempt={attempts} "
              f"native={native_wall:.6f}s solver={solver['status']}", flush=True)
    return dict(n=n, k=k, delta=delta, b=b, attempts=attempts,
                setup=setup, selected=chosen,
                censored=[kind for kind in chosen if len(chosen[kind]) < 2])


def main():
    source_check()
    results = dict(protocol="PROTOCOL.md", seed=SEED, source_hashes=EXPECTED_HASHES,
                   sage_version="10.10.rc0", cells=[])
    for n, k, delta in PLAN:
        results["cells"].append(run_cell(n, k, delta))
        (HERE / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
