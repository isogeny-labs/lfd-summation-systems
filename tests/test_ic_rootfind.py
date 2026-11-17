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
