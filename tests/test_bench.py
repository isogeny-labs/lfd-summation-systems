"""Tests of the benchmark harness: instances, isolation (timeout, crash), the cross-checks that must catch a wrong
solver, the statistics, and the cost model against the figures recorded in experiments/ic-rootfind/README.md."""
from __future__ import annotations

import importlib.util
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from csd.bench import grids, model, runner, solvers, stats      # noqa: E402
from csd.bench.spec import Spec, build                           # noqa: E402


# ------------------------------------------------------------------------------------------------ instances
def test_spec_parse_key_and_validation():
    s = Spec.parse("n=25,k=9,m=3,seed=2,plant=1")
    assert (s.n, s.k, s.m, s.seed, s.plant) == (25, 9, 3, 2, True)
    assert s.key() == "n25-k9-m3-standard-random-planted-s2"
    assert Spec.parse("n=25,k=9,m=3").cell() == Spec(n=25, k=9, m=3).cell()
    with pytest.raises(ValueError):
        Spec(n=10, k=3, m=3, basis="mixed")           # mixed needs j
    with pytest.raises(ValueError):
        Spec(n=10, k=3, m=1)


def test_build_is_deterministic_and_matches_the_rootfind_driver():
    a, b = build(Spec(n=19, k=7, m=3, seed=4)), build(Spec(n=19, k=7, m=3, seed=4))
    assert (a.E.a, a.E.b, a.basis, a.R) == (b.E.a, b.E.b, b.basis, b.R)
    spec = importlib.util.spec_from_file_location("run_ic_rootfind", ROOT / "scripts" / "run_ic_rootfind.py")
    drv = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(drv)
    for kind, j in (("standard", None), ("random", None), ("mixed", 2)):
        E, basis, R = drv.instance(16, 6, 3, kind, j)
        c = build(Spec(n=16, k=6, m=3, seed=3, basis=kind, j=j))
        assert (E.a, E.b, list(basis), tuple(R)) == (c.E.a, c.E.b, c.basis, c.R), kind


def test_planted_target_has_a_decomposition_and_is_a_different_stream():
    for seed in range(4):
        inst = build(Spec(n=13, k=5, m=3, seed=seed, plant=True))
        rows = runner.run_instance(inst.spec, ["py_ref"], reps=1)
        assert rows[0]["status"] == "ok" and rows[0]["n_decs"] >= 1


def test_koblitz_curve_and_invariant_basis():
    inst = build(Spec(n=11, k=1, m=2, curve="koblitz", basis="invariant"))
    assert inst.E.b == 1 and len(inst.basis) == 1
    with pytest.raises(ValueError):
        build(Spec(n=11, k=5, m=2, curve="koblitz", basis="invariant"))   # n = 11 prime: dims 0, 1, 10, 11 only


# -------------------------------------------------------------------------------- runs, isolation, cross-checks
def test_solvers_agree_and_rows_are_complete():
    spec = Spec(n=13, k=5, m=3, seed=1, plant=True)
    rows = runner.run_instance(spec, ["py_ref", "enum", "lrf", "mitm4"], reps=2, aa=True, order_seed=7)
    by = {(r["solver"], r["arm"]): r for r in rows}
    assert by[("mitm4", "A")]["status"] == "skipped"
    for s in ("py_ref", "enum", "lrf"):
        for arm in "AB":
            r = by[(s, arm)]
            assert r["status"] == "ok", r
            assert r["n_decs"] >= 1 and r["t_min"] > 0 and len(r["times"]) == 2 and r["rss_kib"] > 0
            assert r["decs_sha"] == by[("py_ref", "A")]["decs_sha"]
    assert by[("lrf", "A")]["checked_against"] == "py_ref"
    # the order is a function of (order_seed, instance), and positions are all distinct
    again = runner.run_instance(spec, ["py_ref", "enum", "lrf", "mitm4"], reps=1, aa=True, order_seed=7)
    assert [(r["solver"], r["arm"]) for r in again] == [(r["solver"], r["arm"]) for r in rows]
    assert sorted(r["pos"] for r in rows) == list(range(len(rows)))


def test_mitm4_runs_for_m4():
    rows = runner.run_instance(Spec(n=13, k=4, m=4, seed=0, plant=True), ["py_ref", "enum", "lrf", "mitm4"], reps=1)
    assert {r["solver"]: r["status"] for r in rows} == {"py_ref": "ok", "enum": "ok", "lrf": "ok", "mitm4": "ok"}


class _Faulty(solvers.Solver):
    name = "faulty"
    oracle_rank = 50
    mode = "empty"

    def run(self, ctx, inst, threads=1):
        if self.mode == "empty":
            return [], {}
        if self.mode == "wrong_sum":                 # a point on the curve, but not a decomposition of R
            import random
            P = inst.E.random_point(random.Random(5))
            Q = (P[0], P[0] ^ P[1])                     # -P: m copies of P or -P add up to R only by accident
            return [tuple([P] * (inst.spec.m - 1) + [Q])], {}
        if self.mode == "sleep":
            time.sleep(30)
        if self.mode == "die":
            os._exit(3)
        if self.mode == "raise":
            raise RuntimeError("boom")


@pytest.fixture
def faulty(monkeypatch):
    f = _Faulty()
    monkeypatch.setitem(solvers.REGISTRY, "faulty", f)
    return f


def _status(rows, name):
    return next(r for r in rows if r["solver"] == name)


def test_cross_check_catches_a_solver_that_misses_decompositions(faulty):
    faulty.mode = "empty"
    rows = runner.run_instance(Spec(n=13, k=5, m=3, seed=1, plant=True), ["py_ref", "faulty"], reps=1)
    assert _status(rows, "py_ref")["status"] == "ok"
    r = _status(rows, "faulty")
    assert r["status"] == "mismatch" and "py_ref" in r["reason"]


def test_verification_catches_a_decomposition_that_does_not_sum_to_R(faulty):
    faulty.mode = "wrong_sum"
    rows = runner.run_instance(Spec(n=13, k=5, m=3, seed=1, plant=True), ["faulty"], reps=1)
    assert _status(rows, "faulty")["status"] == "invalid"


def test_timeout_is_enforced_by_the_parent(faulty):
    faulty.mode = "sleep"
    t = time.time()
    rows = runner.run_instance(Spec(n=13, k=5, m=3, seed=1), ["faulty"], reps=1, timeout=1.0)
    assert _status(rows, "faulty")["status"] == "timeout" and time.time() - t < 10


def test_crash_and_exception_cost_one_row(faulty):
    faulty.mode = "die"
    assert _status(runner.run_instance(Spec(n=13, k=5, m=3), ["faulty", "lrf"], reps=1), "faulty")["status"] == "crash"
    faulty.mode = "raise"
    rows = runner.run_instance(Spec(n=13, k=5, m=3), ["faulty", "lrf"], reps=1)
    assert _status(rows, "faulty")["status"] == "error" and "boom" in _status(rows, "faulty")["reason"]
    assert _status(rows, "lrf")["status"] == "ok"             # the sweep goes on


def test_inapplicable_solvers_are_skipped_with_a_reason():
    rows = runner.run_instance(Spec(n=52, k=18, m=3), ["py_ref", "enum", "mutantxl"], reps=1, timeout=5)
    for r in rows:
        assert r["status"] == "skipped" and r["reason"], r


def test_mutantxl_closes_on_a_small_system():
    rows = runner.run_instance(Spec(n=10, k=4, m=3, seed=0), ["mutantxl", "lrf"], reps=1, timeout=120)
    r = _status(rows, "mutantxl")
    assert r["status"] == "ok" and r["counters"]["closed"] and r["counters"]["d_mxl"] >= 3


def test_run_grid_writes_meta_and_resumes(tmp_path):
    specs = grids.smoke(seeds=1)[:2]
    out = tmp_path / "r.jsonl"
    assert runner.run_grid(specs, ["py_ref", "lrf"], out, reps=1, log=lambda *_: None) == 0
    meta = json.loads(out.with_suffix(".meta.json").read_text())
    assert meta["commit"] and "cpu" in meta and set(meta["sources"]) >= {"_rootfind.c", "runner.py"}
    n = len(out.read_text().splitlines())
    runner.run_grid(specs, ["py_ref", "lrf"], out, reps=1, resume=True, log=lambda *_: None)
    assert len(out.read_text().splitlines()) == n


# --------------------------------------------------------------------------------------------------- stats
def test_fit_log2_recovers_a_doubling_per_step_and_ci_covers_it():
    rng = np.random.default_rng(1)
    ks = np.arange(5, 16)
    ys = 3e-6 * 2.0 ** ks * np.exp(rng.normal(0, 0.05, ks.size))
    f = stats.fit_log2(ks, ys)
    assert abs(f["slope"] - 1) < 0.05 and f["lo"] <= 1 <= f["hi"] and f["r2"] > 0.99


def test_boot_ci_paired_ratio_and_noise_floor():
    xs = np.random.default_rng(0).normal(10, 1, 200)
    lo, hi = stats.boot_ci(xs)
    assert lo < 10 < hi
    pr = stats.paired_ratio({1: 4.0, 2: 9.0, 3: 16.0}, {1: 2.0, 2: 3.0, 3: 4.0})
    assert pr["n"] == 3 and abs(pr["gmean"] - (2 * 3 * 4) ** (1 / 3)) < 1e-9
    nf = stats.noise_floor([(1.0, 1.02), (2.0, 1.9), (3.0, 3.3)])
    assert nf["p05"] < nf["median"] < nf["p95"]


# --------------------------------------------------------------------------------------------------- model
def test_model_reproduces_the_readme_figures():
    p = model.ecc2k130(m=3)
    r = model.evaluate(p, 44)
    assert abs(r["mu_log2"] - math.log2(1 / 3)) < 0.1                         # "about a third of a relation"
    scan_years = 2.0 ** 44 * 1850e-9 / (365.25 * 86400)
    assert abs(scan_years - 1.0) < 0.05                                       # "about 1.0 core-years"
    assert abs(r["log2_relations_s"] - math.log2(3e13 * 365.25 * 86400)) < 0.3  # "about 3e13 core-years"
    assert abs(model.rho_iterations_log2(p) - 60.8) < 0.1                    # "about 2^61 group operations"
    assert model.optimize_k(p)["k"] in (44, 45)


def test_model_floor_and_required_speedup_are_consistent():
    p = model.ecc2k130(m=3)
    sp = model.required_speedup(p)
    assert math.isinf(sp["factor_log2"]) and sp["floor"]["log2_total_s"] > model.rho_stage_log2s(p)
    # a cheaper relation finder never raises the total; the total is monotone in t_lin
    from dataclasses import replace
    t = [model.optimize_k(replace(p, t_lin_ns=x, t_cand_ns=x / 7))["log2_total_s"] for x in (1850, 185, 18.5, 1.85)]
    assert t == sorted(t, reverse=True)
    sp4 = model.required_speedup(model.ecc2k130(m=4))
    assert not math.isinf(sp4["factor_log2"]) and sp4["factor_log2"] > 0
    q = model.ecc2k130(m=4, t_lin_ns=model.ecc2k130(m=4).t_lin_ns * 2.0 ** -sp4["factor_log2"],
                       t_cand_ns=model.ecc2k130(m=4).t_cand_ns * 2.0 ** -sp4["factor_log2"])
    assert abs(model.optimize_k(q)["log2_total_s"] - model.rho_stage_log2s(q)) < 0.05


def test_oracle_budget_is_solver_agnostic_and_consistent():
    b3 = model.oracle_budget(model.ecc2k130(m=3))
    assert not b3["feasible"]                                  # below the cost of forming the target
    b4 = model.oracle_budget(model.ecc2k130(m=4))
    assert b4["feasible"] and 1e-3 < b4["scan_s"] < 1.0       # tens of milliseconds per target
    # spending exactly the budget at the best k gives total == rho
    from dataclasses import replace
    p = model.ecc2k130(m=4)
    k = b4["k"]
    need = 2.0 ** b4["relations_log2"]
    rel = need * b4["scan_s"] / 2.0 ** b4["mu_log2"]
    assert abs(math.log2(rel + 2.0 ** model.la_stage_log2s(p, k)) - model.rho_stage_log2s(p)) < 1e-6


def test_model_from_records_reads_the_n131_record():
    kw = model.from_records(ROOT / "experiments" / "ic-rootfind" / "n131.jsonl")
    assert 1700 < kw["t_lin_ns"] < 2000 and 100 < kw["t_add_ns"] < 200


# --------------------------------------------------------------------------------------------------- driver
def test_driver_run_and_summarize(tmp_path):
    out = tmp_path / "s.jsonl"
    cmd = [sys.executable, str(ROOT / "scripts" / "bench.py")]
    r = subprocess.run(cmd + ["run", "--grid", "smoke", "--seeds", "1", "--solvers", "py_ref,enum,lrf", "--reps", "1",
                              "--aa", "--out", str(out)], capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stderr
    s = subprocess.run(cmd + ["summarize", str(out), "--baseline", "lrf"], capture_output=True, text=True, timeout=120)
    assert s.returncode == 0 and "A/A noise floor" in s.stdout and "Speed-up of lrf" in s.stdout, s.stdout + s.stderr
    m = subprocess.run(cmd + ["model", "--preset", "ecc2k130", "--window", "1"], capture_output=True, text=True, timeout=60)
    assert m.returncode == 0 and "floor with a FREE relation finder" in m.stdout


# ------------------------------------------------------------------------------------------------- fixtures
from csd.bench import fixtures                                   # noqa: E402


def test_fixture_manifest_and_instances():
    man = json.loads((fixtures.DIR / "MANIFEST.json").read_text())
    assert set(man["files"]) == set(fixtures.available()) == set(fixtures.SET_NAMES)
    for name, meta in man["files"].items():
        assert len(fixtures.load(name)) == meta["rows"]
    assert fixtures.verify(["smoke", "m3-planted"], log=lambda *_: None) == 0
    for r in fixtures.load("smoke"):
        inst = fixtures.instance(r["id"])
        assert inst.E.on_curve(inst.R)


def test_fixture_tamper_is_detected(tmp_path, monkeypatch):
    import shutil
    shutil.copytree(fixtures.DIR, tmp_path / "f")
    monkeypatch.setattr(fixtures, "DIR", tmp_path / "f")
    p = tmp_path / "f" / "smoke.jsonl"
    p.write_text(p.read_text().replace('"seed": 0', '"seed": 9', 1))
    msgs = []
    assert fixtures.verify(["smoke"], log=msgs.append) == 1
    assert "MANIFEST" in msgs[0]


@pytest.mark.parametrize("name", ["smoke", "m3-planted", "m4", "basis-cliff", "koblitz-invariant"])
def test_fixture_lrf_matches_frozen_set(name):
    solver = solvers.REGISTRY["lrf"]
    rows = fixtures.load(name)
    assert rows
    for r in rows:
        spec = Spec(**{**r["spec"], "fixture": r["id"]})
        if not solver.applicable(spec)[0]:
            continue
        inst = fixtures.instance(spec)
        decs, _ = solver.run(solver.prepare(inst), inst, 1)
        assert fixtures.check_result(r["id"], solvers.canon(decs)) is None, r["id"]


def test_fixture_check_catches_a_wrong_set():
    r = next(r for r in fixtures.load("m3-planted") if r["expected"]["count"] > 0)
    assert fixtures.check_result(r["id"], []) is not None


def test_fixture_run_through_runner(tmp_path):
    specs = fixtures.specs("smoke")[:2]
    out = tmp_path / "fx.jsonl"
    runner.run_grid(specs, ["lrf", "enum"], out, reps=1, timeout=60, log=lambda *_: None)
    rows = [json.loads(l) for l in out.read_text().splitlines()]
    assert rows and all(r["status"] == "ok" for r in rows)
    assert all(r.get("fixture_check") == "match" for r in rows)
