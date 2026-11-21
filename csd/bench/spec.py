"""Instance specifications.  One Spec is one seeded instance: a binary curve, a factor base V and a target R.

The random stream of the ``random`` curve kind is the one of scripts/run_ic_rootfind.py (`instance`), so a Spec with
the same (n, k, seed, basis, j) rebuilds the instance of a record in experiments/ic-rootfind/.
"""
from __future__ import annotations

import random
from dataclasses import asdict, dataclass

from csd import ic_dreg as I

BASES = ("standard", "random", "mixed", "invariant")
CURVES = ("random", "koblitz")


@dataclass(frozen=True)
class Spec:
    n: int
    k: int
    m: int
    seed: int = 0
    basis: str = "standard"      # standard <1,w,..>; random; mixed <1..w^(j-1)> + (k-j) random; invariant (Frobenius-stable)
    j: int | None = None          # only for basis == "mixed"
    curve: str = "random"         # random: a in {0,1}, b random; koblitz: a in {0,1}, b = 1
    plant: bool = False           # R = sum of m random factor-base points: at least one decomposition exists, so the
                                  # cross-check of the solvers is not an agreement on the empty set
    fixture: str | None = None    # "set/id": take the curve, V and R from the frozen row (csd.bench.fixtures), not from
                                  # the random stream; the other fields are then labels

    def __post_init__(self):
        if self.basis not in BASES:
            raise ValueError(f"basis {self.basis!r} not in {BASES}")
        if self.curve not in CURVES:
            raise ValueError(f"curve {self.curve!r} not in {CURVES}")
        if self.basis == "mixed" and self.j is None:
            raise ValueError("basis 'mixed' needs j")
        if not 0 < self.k <= self.n:
            raise ValueError("need 0 < k <= n")
        if self.m < 2:
            raise ValueError("need m >= 2")

    def cell(self) -> tuple:
        """Everything but the seed: the unit that is aggregated over seeds."""
        return (self.n, self.k, self.m, self.basis, self.j, self.curve, self.plant)

    def key(self) -> str:
        return f"n{self.n}-k{self.k}-m{self.m}-{self.basis}{'' if self.j is None else self.j}-{self.curve}{'-planted' if self.plant else ''}-s{self.seed}{'' if self.fixture is None else '@' + self.fixture}"

    def as_dict(self) -> dict:
        return asdict(self)

    @staticmethod
    def parse(text: str) -> "Spec":
        """'n=25,k=9,m=3,basis=standard,seed=0' (j and curve optional)."""
        kv = dict(p.split("=", 1) for p in text.split(",") if p)
        if "fixture" in kv and kv["fixture"] == "":
            del kv["fixture"]
        if "plant" in kv:
            kv["plant"] = kv["plant"].lower() in ("1", "true", "yes")
        for f in ("n", "k", "m", "seed", "j"):
            if kv.get(f) not in (None, ""):
                kv[f] = int(kv[f])
            else:
                kv.pop(f, None)
        return Spec(**kv)


@dataclass
class Instance:
    spec: Spec
    E: object
    basis: list
    R: tuple


def build(spec: Spec) -> Instance:
    if spec.fixture:
        from . import fixtures
        return fixtures.instance(spec)
    n, k = spec.n, spec.k
    F = I.binary_field(n)
    if spec.curve == "random" and spec.basis != "invariant":
        rng = random.Random(f"ic-rootfind/{n}/{k}/{spec.seed}/{spec.basis}/{spec.j}")
        a, b = rng.getrandbits(1), rng.randrange(1, 1 << n)
        E = I.BinaryCurve(F, a, b)
        if spec.basis == "standard":
            basis = I.standard_subspace(F, k)
        elif spec.basis == "random":
            basis = I.random_subspace(F, k, rng)
        else:
            while True:
                basis = I.standard_subspace(F, spec.j) + [rng.getrandbits(n) for _ in range(k - spec.j)]
                if I._f2_rank(basis) == k:
                    break
    else:
        rng = random.Random(f"ic-bench/{spec.curve}/{n}/{k}/{spec.seed}/{spec.basis}/{spec.j}")
        a, b = rng.getrandbits(1), 1 if spec.curve == "koblitz" else rng.randrange(1, 1 << n)
        E = I.BinaryCurve(F, a, b)
        if spec.basis == "invariant":
            subs = I.invariant_subspaces(F, k)
            if not subs:
                raise ValueError(f"no Frobenius-invariant subspace of dimension {k} for n = {n}")
            basis = subs[rng.randrange(len(subs))]
        elif spec.basis == "standard":
            basis = I.standard_subspace(F, k)
        elif spec.basis == "random":
            basis = I.random_subspace(F, k, rng)
        else:
            while True:
                basis = I.standard_subspace(F, spec.j) + [rng.getrandbits(n) for _ in range(k - spec.j)]
                if I._f2_rank(basis) == k:
                    break
    R = E.random_point(rng)
    if spec.plant:
        R = _planted(E, list(basis), spec)
    return Instance(spec, E, list(basis), tuple(R))


def _planted_with_points(E, basis, spec):
    """(R, [P_1..P_m]) with R = P_1 + ... + P_m, x(P_i) uniform in V and a random sign; resampled if the sum is O."""
    rng = random.Random(f"ic-bench/plant/{spec.key()}")
    while True:
        acc, pts = None, []
        for _ in range(spec.m):
            P = None
            while P is None:
                x = 0
                for v in basis:
                    if rng.getrandbits(1):
                        x ^= v
                P = E.lift_x(x)
            P = P if rng.getrandbits(1) else E.neg(P)
            pts.append(tuple(P))
            acc = E.add(acc, P)
        if acc is not None:
            return tuple(acc), pts


def _planted(E, basis, spec):
    return _planted_with_points(E, basis, spec)[0]


def planted_points(spec: Spec):
    """The planted summands of a spec built from the random stream (used to write fixtures)."""
    inst = build(Spec(**{**spec.as_dict(), "plant": False, "fixture": None}))
    return _planted_with_points(inst.E, inst.basis, spec)[1]
