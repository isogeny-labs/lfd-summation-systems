"""Independent arithmetic and exhaustive controls for the matched Semaev experiment."""
from itertools import permutations, product

from csd import ic_dreg as I
from csd import semaev_matched as M


def independent_mul(x, y, n, modulus):
    value = 0
    for j in range(n):
        if y >> j & 1:
            value ^= x << j
    for j in range(2*n-2, n-1, -1):
        if value >> j & 1:
            value ^= modulus << (j-n)
    return value


def test_field_products_and_sylvester_determinant():
    """S4 closed form against the 4x4 Sylvester determinant, not another formula."""
    F = I.binary_field(3)
    for x, y in product(range(8), repeat=2):
        assert F.mul(x, y) == independent_mul(x, y, 3, F.modulus)
    for x, y, z, r in product(range(4), repeat=4):
        C, B, A = M.quadratic_coefficients(F, 3, x, y)
        c, b, a = M.quadratic_coefficients(F, 3, z, r)
        matrix = [[A, B, C, 0], [0, A, B, C], [a, b, c, 0], [0, a, b, c]]
        det = 0
        for perm in permutations(range(4)):
            term = 1
            for i, j in enumerate(perm):
                term = independent_mul(term, matrix[i][j], 3, F.modulus)
            det ^= term  # signs coincide in characteristic two
        assert det == I.s4_binary_value(F, 3, x, y, z, r)


def test_infinity_counterexample_and_complete_pdp_correction():
    """The affine split really drops tuples; the correction is checked against points."""
    E = I.BinaryCurve(I.binary_field(2), 0, 1)
    V = list(range(4))
    oracle = M.decomposition_map(E, V)
    saw_missing = False
    for R in M.rational_points(E):
        sets = M.coordinate_sets(E, V, R[0])
        assert sets['pdp_direct'] == sets['pdp_split'] == oracle.get(R, set())
        missing = sets['pdp_direct'] - sets['pdp_finite']
        saw_missing |= bool(missing)
        assert missing <= sets['infinity']
    assert saw_missing


def test_descent_on_every_boolean_assignment():
    """Both descended systems versus direct extension-field evaluation, including membership."""
    E = I.BinaryCurve(I.binary_field(3), 1, 3)
    basis = [1, 2]
    points = M.rational_points(E)
    r = points[-1][0]
    liftable = {P[0] for P in points}
    for formulation in ('original', 'split'):
        eqs, nv = M.build_system(E, basis, r, formulation)
        for mask in range(1 << nv):
            x, y, z = M.decode(mask, basis)
            if formulation == 'original':
                zero = I.s4_binary_value(E.F, E.b, x, y, z, r) == 0
            else:
                u = mask >> 6
                zero = (I.s3_binary_value(E.F, E.b, x, y, u) == 0 and
                        I.s3_binary_value(E.F, E.b, u, z, r) == 0)
            want = zero and all(t in liftable for t in (x, y, z))
            assert all(M.eval_boolean(p, mask) == 0 for p in eqs) == want


def test_saved_matched_record():
    """Measured-source hashes, exact audit replay, case selection and saved summaries."""
    from pathlib import Path
    import subprocess
    import sys
    root = Path(__file__).resolve().parents[1]
    subprocess.run([sys.executable, str(root/'scripts/check_semaev_matched.py')],
                   cwd=root, check=True, capture_output=True, text=True)
