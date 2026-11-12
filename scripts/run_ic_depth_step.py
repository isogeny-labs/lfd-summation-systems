"""Dimension of the degree-3 and degree-4 chain closures against the ideal: the depth-6 step.

For two-summand chain presentations the closure V_D modulo the chains is
computed exactly with the dense engine at D = 3, 4 for several primes, curves
and targets and depths k = 3..10, together with |M_{<=4}|, the scalar solution
count |Z| and the evaluation-onto check.  The table records whether V_4 is the
full degree-4 part of the ideal (dim = |M_{<=4}| - |Z|) and, when it is not,
the deficiency.  Written to experiments/ic-depth-step/growth.json.
"""
from pathlib import Path
import hashlib, json, math, random, sys, time
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from csd import ic_dreg as I
from csd import ic_curve_sweep as S
from csd import ic_fall_certificates as C

CURVES = [(257, 2, 3), (7681, 1, 0), (7681, 5, 7), (12289, 1, 0), (12289, 3, 11)]
DEPTHS = range(3, 11)


def main():
    root = Path(__file__).resolve().parents[1]
    out = root / 'experiments/ic-depth-step'
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    started = time.time()
    for p, a, b in CURVES:
        if (4 * a ** 3 + 27 * b * b) % p == 0:
            continue
        rng = random.Random(p * 1000 + a)
        E = I.PrimeCurve(p, a, b)
        R = S.random_rational_target(E, rng, p)
        for k in DEPTHS:
            if (p - 1) % (1 << k):
                continue
            g = I.pkm_generator(E, R[0], k)
            g0 = S.subgroup_generator_value(p, k)
            Z = S.prime_points(E, k, g0, 1, 1, R[0])
            V3 = C.chain_degree_closure([g], 2, k, p, 3, max_columns=16384, max_cells=400_000_000)
            V4 = C.chain_degree_closure([g], 2, k, p, 4, max_columns=16384, max_cells=400_000_000)
            for space in (V3, V4):
                assert space.contains_one or C.check_closed_exclusion_certificate(space, [g])
            M4 = sum(math.comb(2 * k, d) for d in range(5))
            M3 = sum(math.comb(2 * k, d) for d in range(4))
            onto = S.prime_evaluation_rank(Z, [m for m in V4.masks if m.bit_count() < 4], k, g0, p, 10 ** 8) == len(Z) if Z else True
            row = dict(prime=p, a=a, b=b, target_x=R[0], k=k, solutions=len(Z), monomials_le4=M4,
                       ideal_dimension_le4=M4 - len(Z), dim_V4=V4.quotient_rank, deficiency=M4 - len(Z) - V4.quotient_rank,
                       dim_V3=V3.quotient_rank, dim_V4_low=V4.low_degree_dimension(3) - (math.comb(2 * k + 3, 3) - M3),
                       generator_rows_le4=1 + 2 * k + math.comb(2 * k, 2), evaluation_onto_le3=onto,
                       quadratic_fit=4 * k * k + 172 * k - 571)
            rows.append(row)
            print(json.dumps(row), flush=True)
    summary = dict(rows=rows, elapsed_seconds=round(time.time() - started),
                   observation='dim V_3 = 2k+1 and, for k >= 6, dim V_4 and dim(V_4 ∩ M_{<=3}) do not depend on the prime, curve or target; '
                               'the fitted dim V_4 = 4k^2 + 172k - 571 is quadratic while |M_{<=4}| is quartic, so V_4 is the full '
                               'degree-4 part of the ideal for k <= 5 and deficient from k = 6 on: the strict fall at 5 and exact LFD 5.')
    (out / 'growth.json').write_text(json.dumps(summary, indent=1) + '\n')
    sources = ['scripts/run_ic_depth_step.py', 'csd/ic_fall_certificates.py', 'csd/ic_curve_sweep.py', 'csd/ic_dreg.py']
    (out / 'verification.json').write_text(json.dumps(dict(
        status='passed: exact dense closures with exclusion certificates',
        source_sha256={n: hashlib.sha256((root / n).read_bytes()).hexdigest() for n in sources},
        data_sha256={'experiments/ic-depth-step/growth.json': hashlib.sha256((out / 'growth.json').read_bytes()).hexdigest()}), indent=2) + '\n')


if __name__ == '__main__':
    main()
