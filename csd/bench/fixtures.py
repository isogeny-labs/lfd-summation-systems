"""Frozen benchmark fixtures: instances with their expected answers, independent of the random stream.

A fixture row is one instance written out in full (field modulus, curve, factor-base basis, target R, hex strings so
that other languages can read it) together with what every correct solver must return: the full set of
decompositions when two independent solvers agreed on it (`complete`), or the planted summands that must be in it
(`must_contain`) when only the fast kernel can reach the instance.  Rows may carry a `reference` (a number recorded
in experiments/ic-rootfind/ for the same instance) so that a run can print its ratio to the record.

    experiments/ic-bench/fixtures/<set>.jsonl   one row per instance
    experiments/ic-bench/fixtures/MANIFEST.json sha256 of every set file and the row counts
    experiments/ic-bench/fixtures/SCHEMA.md     the row format

`make` writes a set from its generator (deterministic, solvers cross-checked while writing); `verify` recomputes
from the stored data only (no random stream) and compares.  Use through scripts/bench.py (`fixtures ...`,
`run --fixture SET`).
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

from csd import ic_dreg as I

from . import solvers as S
from .spec import Instance, Spec, build, planted_points

ROOT = Path(__file__).resolve().parents[2]
DIR = ROOT / "experiments" / "ic-bench" / "fixtures"
RECORDS = ROOT / "experiments" / "ic-rootfind"
FORMAT = 1
MAX_STORED_DECS = 64          # larger sets are stored as count + hash only


# ------------------------------------------------------------------------------------------------ set generators
def _tag(spec: Spec) -> str:
    return (f"n{spec.n}-k{spec.k:02d}-m{spec.m}-{spec.basis}{'' if spec.j is None else spec.j}"
            f"{'-' + spec.curve if spec.curve != 'random' else ''}{'-planted' if spec.plant else ''}-s{spec.seed}")


def _sets() -> dict:
    """name -> (description, [Spec], reference file or None)."""
    sp = Spec
    out = {}
    out["smoke"] = ("Seconds: tiny instances with planted decompositions; every solver applies.",
                    [sp(n=3 * k - 2, k=k, m=3, seed=s, plant=True) for k in (4, 5, 6) for s in (0, 1)]
                    + [sp(n=13, k=4, m=4, seed=s, plant=True) for s in (0, 1)], None)
    out["m3-natural"] = ("m = 3, V = <1, w, .., w^(k-1)>, n = 3k - 2, k = 5..12, random target: the instances of "
                         "experiments/ic-rootfind/scaling-m3.jsonl (most have no decomposition).",
                         [sp(n=3 * k - 2, k=k, m=3, seed=s) for k in range(5, 13) for s in range(3)], "scaling-m3.jsonl")
    out["m3-planted"] = ("Same sizes with a planted decomposition: at least one hit per instance.",
                         [sp(n=3 * k - 2, k=k, m=3, seed=s, plant=True) for k in range(5, 13) for s in range(3)], None)
    out["m3-large-planted"] = ("k = 13..16, planted: LRF takes milliseconds to a second; enumeration (2^26 to 2^32 "
                               "additions, about 3 minutes for the set) was run once to freeze the full expected sets.",
                               [sp(n=3 * k - 2, k=k, m=3, seed=s, plant=True) for k in range(13, 17) for s in range(2)], None)
    out["m4"] = ("m = 4, n = 4k - 3, k = 4..8, random targets (instances of experiments/ic-rootfind/m4.jsonl) and planted.",
                 [sp(n=4 * k - 3, k=k, m=4, seed=s, plant=p) for k in range(4, 9) for s in range(3) for p in (False, True)],
                 "m4.jsonl")
    out["basis-cliff"] = ("n = 31, k = 10, m = 3: <1..w^(j-1)> plus k - j random vectors, j = 10, 9, 8, 0; planted. "
                          "The factor base decides how much LRF gains (ic-rootfind README).",
                          [sp(n=31, k=10, m=3, seed=s, basis="standard" if j == 10 else "mixed", j=None if j == 10 else j,
                              plant=True) for j in (10, 9, 8, 0) for s in range(3)], None)
    out["regime-k12"] = ("k = 12, m = 3, n = 27..45: from free unknowns in the linear system to none "
                         "(experiments/ic-rootfind/regime-k12.jsonl).",
                         [sp(n=n, k=12, m=3, seed=s) for n in range(27, 46, 2) for s in range(3)], "regime-k12.jsonl")
    out["koblitz-invariant"] = ("Koblitz curves y^2 + xy = x^3 + a x^2 + 1 with a Frobenius-invariant factor base, planted: "
                                "the trace-zero hyperplane (k = n - 1, n = 11, 13, m = 2; there V2 is the whole field and "
                                "LRF degenerates to many candidates) and sums of the irreducible factors of x^n - 1 for "
                                "n = 15, 21 (m = 3).",
                                [sp(n=n, k=n - 1, m=2, seed=s, basis="invariant", curve="koblitz", plant=True)
                                 for n in (11, 13) for s in range(3)]
                                + [sp(n=n, k=k, m=3, seed=s, basis="invariant", curve="koblitz", plant=True)
                                   for n, k in ((15, 5), (15, 6), (21, 7), (21, 9)) for s in range(3)], None)
    return out


SET_NAMES = ("smoke", "m3-natural", "m3-planted", "m3-large-planted", "m4", "basis-cliff", "regime-k12",
             "koblitz-invariant", "rust-records")


# ---------------------------------------------------------------------------------------------------- row format
def _hx(v: int) -> str:
    return format(int(v), "x")


def _hp(P) -> list:
    return [_hx(P[0]), _hx(P[1])]


def _unp(P) -> tuple:
    return (int(P[0], 16), int(P[1], 16))


def _row_from_instance(set_name, spec: Spec, inst: Instance, expected, reference, provenance) -> dict:
    F = inst.E.F
    return dict(
        id=f"{set_name}/{_tag(spec)}", format=FORMAT, spec={**spec.as_dict(), "fixture": None},
        n=spec.n, k=spec.k, m=spec.m, modulus=_hx(F.modulus), a=inst.E.a, b=_hx(inst.E.b),
        basis=[_hx(v) for v in inst.basis], R=_hp(inst.R), expected=expected, reference=reference,
        provenance=provenance)


def _expected_block(complete, decs, verified_by, planted=None) -> dict:
    canon = S.canon(decs) if decs is not None else None
    blk = dict(complete=complete, verified_by=verified_by)
    if canon is not None:
        blk.update(count=len(canon), sha=hashlib.sha1(repr(canon).encode()).hexdigest()[:16])
        if len(canon) <= MAX_STORED_DECS:
            blk["decs"] = [[_hp(P) for P in d] for d in canon]
    if planted is not None:
        blk["must_contain"] = [_hp(P) for P in sorted(planted)]
    return blk


def _solve_all(inst: Instance):
    """Run every applicable decomposition solver in-process; return {name: canon decs}."""
    out = {}
    for name, solver in S.REGISTRY.items():
        if solver.kind != "decomposition":
            continue
        ok, _ = solver.available()
        if ok:
            ok, _ = solver.applicable(inst.spec)
        if not ok:
            continue
        ctx = solver.prepare(inst)
        decs, _ = solver.run(ctx, inst, 1)
        out[name] = S.canon(decs)
    return out


def _reference(ref_file, spec: Spec):
    if not ref_file:
        return None
    path = RECORDS / ref_file
    for line in path.read_text().splitlines():
        r = json.loads(line)
        if (r["n"], r["k"], r["seed"], r.get("m", 3)) == (spec.n, spec.k, spec.seed, spec.m) and r.get("basis", "standard") == "standard" and not spec.plant:
            keep = {kk: r[kk] for kk in ("t_lrf", "t_brute", "t_mitm", "decompositions", "nfb", "dimV2") if kk in r}
            keep.update(file=f"experiments/ic-rootfind/{ref_file}", note="single thread, min of runs, 4-core Xeon, Python call included")
            return keep
    return None


def generate(name: str, log=print) -> list[dict]:
    if name == "rust-records":
        return _rust_records()
    desc, specs, ref_file = _sets()[name]
    rows = []
    for spec in specs:
        inst = build(spec)
        results = _solve_all(inst)
        vals = {tuple(v) for v in results.values()}
        if len(vals) != 1:
            raise RuntimeError(f"solvers disagree on {_tag(spec)}: {sorted(results)}")
        decs = next(iter(results.values()))
        independent = [s for s in results if s in ("py_ref", "enum", "mitm4")]
        complete = bool(independent)                # two kernels of different design agreed, or the reference did
        planted = planted_points(spec) if spec.plant else None
        if planted is not None and tuple(sorted(planted)) not in set(decs):
            raise RuntimeError(f"planted decomposition not returned for {_tag(spec)}")
        exp = _expected_block(complete, decs, sorted(results), planted)
        rows.append(_row_from_instance(name, spec, inst, exp, _reference(ref_file, spec),
                                       dict(generator="csd.bench.fixtures.generate", spec_key=spec.key(),
                                            python_random_stream=True, description=desc)))
        log(f"  {rows[-1]['id']}: {exp.get('count')} decompositions, verified by {','.join(sorted(results))}")
    return rows


def _rust_records() -> list[dict]:
    """The instances of the Rust kernel's records (n up to 163): data and the recorded cost per target; no expected set
    (the Rust crate plants one decomposition and scans a range; the planted summands are not recorded)."""
    rows, seen = [], set()
    for fname in ("n131.jsonl", "rust-scaling.jsonl"):
        for line in (RECORDS / fname).read_text().splitlines():
            r = json.loads(line)
            if r.get("threads", 1) != 1 or r.get("form") != "rows":
                continue
            key = (r["n"], r["k"], r["seed"])
            if key in seen:
                continue
            seen.add(key)
            F = I.binary_field(r["n"])
            low = sum(1 << t for t in r["modulus_low_terms"])
            if F.modulus != (1 << r["n"]) | low:
                raise RuntimeError(f"field modulus of n = {r['n']} differs from the record's")
            spec = Spec(n=r["n"], k=r["k"], m=3, seed=r["seed"])
            rows.append(dict(
                id=f"rust-records/n{r['n']}-k{r['k']:02d}-s{r['seed']}", format=FORMAT, spec=spec.as_dict(), n=r["n"], k=r["k"],
                m=3, modulus=_hx(F.modulus), a=r["a"], b=r["b"], basis=[_hx(v) for v in I.standard_subspace(F, r["k"])],
                R=[r["xR"], r["yR"]],
                expected=dict(complete=False, verified_by=[], note="one decomposition is planted in the scanned range "
                              "(planted_found is true in the record); the summands are not recorded"),
                reference=dict(file=f"experiments/ic-rootfind/{fname}",
                               **{kk: r[kk] for kk in ("ns_per_target", "ns_per_abscissa", "log2_range", "range_lo", "unknowns",
                                                       "equations", "dimV2", "lin_solvable", "targets") if kk in r},
                               note="Rust kernel, single thread, rows form"),
                provenance=dict(generator="csd.bench.fixtures._rust_records", source=fname)))
    return rows


# --------------------------------------------------------------------------------------------------- files
def path(name: str) -> Path:
    return DIR / f"{name}.jsonl"


def write(name: str, rows: list[dict]):
    DIR.mkdir(parents=True, exist_ok=True)
    path(name).write_text("".join(json.dumps(r, sort_keys=True) + "\n" for r in rows))
    _update_manifest()


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _update_manifest():
    files = {p.stem: dict(sha256=_sha(p), rows=len(p.read_text().splitlines())) for p in sorted(DIR.glob("*.jsonl"))}
    (DIR / "MANIFEST.json").write_text(json.dumps(dict(format=FORMAT, files=files), indent=1, sort_keys=True) + "\n")


def load(name: str) -> list[dict]:
    return [json.loads(l) for l in path(name).read_text().splitlines() if l.strip()]


def available() -> list[str]:
    return [p.stem for p in sorted(DIR.glob("*.jsonl"))]


_ROW_CACHE: dict[str, dict] = {}


def row(fid: str) -> dict:
    set_name = fid.split("/", 1)[0]
    if set_name not in _ROW_CACHE or fid not in _ROW_CACHE.get(set_name, {}):
        _ROW_CACHE.setdefault(set_name, {}).update({r["id"]: r for r in load(set_name)})
    return _ROW_CACHE[set_name][fid]


def specs(name: str) -> list[Spec]:
    """Specs that build their instance from the fixture row (not from the random stream)."""
    return [Spec(**{**r["spec"], "fixture": r["id"]}) for r in load(name)]


def instance(spec_or_id) -> Instance:
    fid = spec_or_id if isinstance(spec_or_id, str) else spec_or_id.fixture
    r = row(fid)
    F = I.binary_field(r["n"])
    if F.modulus != int(r["modulus"], 16):
        raise RuntimeError(f"{fid}: the field modulus of csd.ic_dreg changed ({F.modulus:x} vs {r['modulus']})")
    E = I.BinaryCurve(F, r["a"], int(r["b"], 16))
    spec = spec_or_id if not isinstance(spec_or_id, str) else Spec(**{**r["spec"], "fixture": fid})
    return Instance(spec, E, [int(v, 16) for v in r["basis"]], _unp(r["R"]))


def expected(fid: str):
    """(canonical expected set or None, must_contain tuple or None, complete flag)."""
    e = row(fid)["expected"]
    decs = S.canon([[_unp(P) for P in d] for d in e["decs"]]) if "decs" in e else None
    must = tuple(sorted(_unp(P) for P in e.get("must_contain", []))) or None
    return decs, must, bool(e.get("complete")), e


def check_result(fid: str, decs) -> str | None:
    """None if `decs` (canonical) is what the fixture requires, else a message."""
    exp, must, complete, e = expected(fid)
    if "sha" in e and complete:
        sha = hashlib.sha1(repr(decs).encode()).hexdigest()[:16]
        if sha != e["sha"] or len(decs) != e["count"]:
            return f"differs from the frozen set ({len(decs)} decompositions, sha {sha}; frozen {e['count']}, {e['sha']})"
    if must:
        flat = {tuple(sorted(d)) for d in decs}
        if must not in flat and not any(set(must) == set(d) for d in decs):
            return "does not contain the planted decomposition"
    return None


# ------------------------------------------------------------------------------------------------- verify
def verify(names=None, deep=False, log=print) -> int:
    """Re-derive what can be re-derived from the stored data alone.  Returns the number of problems."""
    problems = 0
    man = json.loads((DIR / "MANIFEST.json").read_text())
    for name in names or available():
        p = path(name)
        if name not in man["files"] or _sha(p) != man["files"][name]["sha256"]:
            log(f"{name}: file does not match MANIFEST.json")
            problems += 1
            continue
        for r in load(name):
            fid = r["id"]
            inst = instance(fid)
            if not inst.E.on_curve(inst.R):
                log(f"{fid}: R is not on the curve")
                problems += 1
                continue
            e = r["expected"]
            if not e.get("verified_by"):
                continue
            results = {}
            for sname in e["verified_by"] if deep else [s for s in ("lrf",) if s in e["verified_by"]]:
                solver = S.REGISTRY[sname]
                if not (solver.available()[0] and solver.applicable(inst.spec)[0]):
                    continue
                decs, _ = solver.run(solver.prepare(inst), inst, 1)
                results[sname] = S.canon(decs)
            for sname, decs in results.items():
                bad = check_result(fid, decs)
                if bad:
                    log(f"{fid}: {sname} {bad}")
                    problems += 1
    log("fixtures verified" if not problems else f"{problems} problem(s)")
    return problems


def rng_drift(names=None) -> list[str]:
    """Instances whose stored data differ from what the current random stream builds (informational: the fixture
    does not depend on it, but a difference means Python's `random` or the field choice changed)."""
    out = []
    for name in names or available():
        for r in load(name):
            if not r["provenance"].get("python_random_stream"):
                continue
            spec = Spec(**{**r["spec"], "fixture": None})
            inst = build(spec)
            if (_hx(inst.E.b), inst.E.a, [_hx(v) for v in inst.basis], _hp(inst.R)) != (r["b"], r["a"], r["basis"], r["R"]):
                out.append(r["id"])
    return out
