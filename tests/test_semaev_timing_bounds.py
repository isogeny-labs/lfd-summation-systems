"""Censored bounds must never turn partial wall time into exact solver time."""
import importlib.util
import math
from pathlib import Path

path=Path(__file__).resolve().parents[1]/'scripts/analyze_semaev_timing_bounds.py'
spec=importlib.util.spec_from_file_location('semaev_bounds',path)
B=importlib.util.module_from_spec(spec)
spec.loader.exec_module(B)


def test_dominated_censored_baseline_is_resolved():
    original=B.median_bounds([dict(status='ok',total_seconds=t) for t in (6,7,8)])
    split=B.median_bounds([dict(status='gb_timeout',partial_wall_seconds=30)]*3)
    best=(min(original[0],split[0]),min(original[1],split[1]))
    assert best==(7,7)
    assert B.gain_bounds(best,(1,1))==(7,7)


def test_both_gb_timeouts_only_give_a_lower_bound():
    base=B.median_bounds([dict(status='gb_timeout')]*3)
    assert base==(10,math.inf)
    assert B.gain_bounds(base,(2,2))==(5,math.inf)


def test_process_timeout_or_missing_run_is_not_a_solver_lower_bound():
    assert B.median_bounds([dict(status='process_timeout',process_wall_seconds=40)]*3)==(0,math.inf)
    assert B.median_bounds([])==(0,math.inf)
    assert B.gain_bounds((0,math.inf),(1,1))==(0,math.inf)
