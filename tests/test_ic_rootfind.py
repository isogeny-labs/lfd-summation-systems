"""csd.ic_rootfind: the C kernels against csd.ic_dreg's field, group law and plain enumeration."""
import random

from csd import ic_dreg as I
from csd.ic_rootfind import Solver, _lib


def run(ns, trials, seed):
    rng = random.Random(seed)
    instances = nonempty = 0
    for n in ns:
        F = I.binary_field(n)
        for trial in range(trials):
            a = rng.getrandbits(1)
            b = rng.randrange(1, 1 << n) if trial % 2 else 1
            E = I.BinaryCurve(F, a, b)
            for k in range(2, min(6, n - 1)):
                basis = I.standard_subspace(F, k) if trial % 3 == 0 else I.random_subspace(F, k, rng)
                S = Solver(E, basis)
                for _ in range(20):
                    x, y = rng.randrange(1, 1 << n), rng.randrange(1 << n)
                    assert _lib().rf_mul(S._ctx, x, y) == F.mul(x, y)
                    assert _lib().rf_inv(S._ctx, x) == F.inv(x)
                fb = sorted(P for x in I.span(basis) for P0 in [E.lift_x(x)] if P0 for P in {P0, E.neg(P0)})
                assert sorted(zip(map(int, S.fx), map(int, S.fy))) == fb
                for m in (2, 3) + ((4,) if k <= 3 else ()):
                    R = E.random_point(rng)
                    ref = S.decompose_reference(m, R)
                    for P in (pt for dec in ref for pt in dec):
                        assert E.on_curve(P)
                    assert S.lrf(m, R)[0] == ref
                    assert S.lrf_image(m, R)[0] == ref
                    assert S.brute(m, R)[0] == ref
                    if m == 2 and R[0]:
                        assert S.d2_reference(R) == ref
                    if m == 4:
                        assert S.mitm4(R)[0] == ref
                    instances += 1
                    nonempty += bool(ref)
    return instances, nonempty


def test_small_instances():
    instances, nonempty = run([7, 8, 9, 11, 12, 13, 16, 17], 6, 1)
    assert instances == 480 and nonempty == 87


def test_conditional_image_table_and_diagonal():
    F = I.binary_field(7)
    E = I.BinaryCurve(F, 0, 1)
    basis = I.standard_subspace(F, 3)
    S = Solver(E, basis)
    fb = sorted(zip(map(int, S.fx), map(int, S.fy)))
    P = next(p for p in fb if E.mul(2, p) is not None)
    R = E.mul(2, P)
    assert (P, P) in S.lrf_image(2, R)[0]
    assert S.lrf_image(2, R)[0] == S.lrf(2, R)[0]
    for mask in range(1 << S.k):
        e1 = 0
        for j, v in enumerate(basis):
            if mask >> j & 1:
                e1 ^= v
        image = {F.sq(x) ^ F.mul(e1, x) for x in I.span(basis)}
        pivots = S._image_table[mask*S.k:(mask+1)*S.k]
        for e2 in range(1 << F.N):
            value = e2
            for pivot in pivots:
                value = min(value, value ^ int(pivot))
            assert (value == 0) == (e2 in image)


def test_zero_abscissa_target():
    """x(T) = 0 (the 2-torsion point) is the one target the line parametrization cannot take."""
    F = I.binary_field(11)
    E = I.BinaryCurve(F, 1, 5)
    S = Solver(E, I.standard_subspace(F, 4))
    T = E.lift_x(0)
    assert S.lrf(2, T)[0] == S.decompose_reference(2, T)
    # R = P + P' + T with P, P' in the factor base: the guess P leaves the target P' + T, the guess P' leaves
    # P + T, and T itself is in the factor base (0 is in V)
    P, P2 = (int(S.fx[3]), int(S.fy[3])), (int(S.fx[6]), int(S.fy[6]))
    R = E.add(E.add(T, P), P2)
    ref = S.decompose_reference(3, R)
    assert tuple(sorted([T, P, P2])) in ref
    assert S.lrf(3, R)[0] == ref
    R = E.add(P, P2)                      # the guess T leaves a generic target; the guess P leaves P2 + T - ...
    assert S.lrf(3, E.add(R, T))[0] == S.decompose_reference(3, E.add(R, T))
    assert S.lrf(2, E.add(P, T))[0] == S.decompose_reference(2, E.add(P, T)) != []


def test_threads_and_larger_instance():
    """n = 28, k = 10, m = 3: 2^20 enumerated pairs against 2^10 linear systems, one and four threads."""
    rng = random.Random(5)
    F = I.binary_field(28)
    E = I.BinaryCurve(F, 0, 1)
    S = Solver(E, I.standard_subspace(F, 10))
    hits = 0
    for _ in range(12):
        R = E.random_point(rng)
        ref = S.brute(3, R)[0]
        for dec in ref:
            acc = None
            for P in dec:
                acc = E.add(acc, P)
            assert acc == R
        assert S.lrf(3, R)[0] == ref == S.lrf(3, R, threads=4)[0] == S.brute(3, R, threads=4)[0]
        hits += len(ref)
    assert hits > 0


def test_rust_kernel_against_c_kernel():
    """rust/ic_rootfind (multi-word fields, row form of the linear system) on the instances of the C kernel."""
    import json
    import shutil
    import subprocess
    from pathlib import Path

    import pytest

    if shutil.which("cargo") is None:
        pytest.skip("no cargo")
    crate = Path(__file__).resolve().parents[1] / "rust" / "ic_rootfind"
    subprocess.run(["cargo", "build", "--release", "--quiet"], cwd=crate, check=True)
    exe = crate / "target" / "release" / "ic_rootfind"
    assert "SELFTEST PASSED" in subprocess.run([exe, "selftest"], capture_output=True, text=True, check=True).stdout
    rng = random.Random(11)
    hits = 0
    for n, k in [(13, 5), (16, 6), (19, 7), (25, 9), (28, 10), (31, 11), (40, 14), (61, 12)]:
        F = I.binary_field(n)
        for _ in range(3):
            E = I.BinaryCurve(F, rng.getrandbits(1), rng.randrange(1, 1 << n))
            S = Solver(E, I.standard_subspace(F, k))
            for m in (2, 3):
                R = E.random_point(rng)
                if m == 2:                   # a target that has a decomposition
                    R = E.add((int(S.fx[1]), int(S.fy[1])), (int(S.fx[-1]), int(S.fy[-1])))
                out = subprocess.run([exe, "solve", "--n", str(n), "--k", str(k), "--m", str(m), "--a", str(E.a),
                                      "--b", f"{E.b:x}", "--xr", f"{R[0]:x}", "--yr", f"{R[1]:x}"],
                                     capture_output=True, text=True, check=True).stdout
                got = sorted(tuple(sorted((int(x, 16), int(y, 16)) for x, y in dec)) for dec in json.loads(out))
                want = S.lrf(m, R)[0]
                assert got == want, (n, k, m)
                hits += len(want)
    assert hits > 20
