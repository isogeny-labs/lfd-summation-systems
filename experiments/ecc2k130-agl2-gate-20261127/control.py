#!/usr/bin/env sage -python
"""Small arithmetic controls for the affine degree-two quotient screen."""

from sage.all import EllipticCurve, GF


def count_cover(m):
    field = GF(2**m, "u")
    curve = EllipticCurve(field, [1, 0, 0, 0, 1])
    # The Artin-Schreier cover w^2+w=y has one point above infinity.
    return 1 + sum(
        sum(w*w+w == point[1] for w in field)
        for point in curve.points() if not point.is_zero()
    )


def main():
    base = EllipticCurve(GF(2), [1, 0, 0, 0, 1])
    assert base.cardinality() == 4
    assert sorted(str(point) for point in base.points()) == [
        "(0 : 1 : 0)", "(0 : 1 : 1)", "(1 : 0 : 1)", "(1 : 1 : 1)",
    ]
    assert (sum((-x) % 131 == x for x in range(131)),
            sum(x < (-x) % 131 for x in range(131))) == (1, 65)
    assert 70**2 * 2 < 99**2

    points = [count_cover(m) for m in (1, 2, 3)]
    assert points == [3, 11, 3]
    traces = [2**m + 1 - points[m-1] for m in (1, 2, 3)]
    a1 = -traces[0]
    a2 = (traces[0]**2-traces[1]) // 2
    a3 = (-traces[0]**3+3*traces[0]*traces[1]-2*traces[2]) // 6
    assert (a1, a2, a3) == (0, 3, -2)
    class_number = 9 + 5*a1 + 3*a2 + a3
    assert class_number == 16 and class_number % 4 == 0
    print("base order 4; affine involution (1 fixed, 65 pairs);")
    print("genus-three cover point counts", points,
          "; Jacobian order", class_number)


if __name__ == "__main__":
    main()
