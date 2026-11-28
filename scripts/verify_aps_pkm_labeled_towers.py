"""Reproduce the two small factor-base calculations discussed in chat.

Run with: python3 scripts/verify_aps_pkm_labeled_towers.py
Only the Python standard library is required.
"""

from itertools import product


def add(p, q, modulus, a):
    if p is None:
        return q
    if q is None:
        return p
    x1, y1 = p
    x2, y2 = q
    if x1 == x2 and (y1 + y2) % modulus == 0:
        return None
    if p == q:
        slope = (3 * x1 * x1 + a) * pow(2 * y1, -1, modulus)
    else:
        slope = (y2 - y1) * pow(x2 - x1, -1, modulus)
    slope %= modulus
    x3 = (slope * slope - x1 - x2) % modulus
    return x3, (slope * (x1 - x3) - y1) % modulus


def multiply(k, p, modulus, a):
    if k < 0:
        return multiply(-k, (p[0], -p[1] % modulus), modulus, a)
    result = None
    while k:
        if k & 1:
            result = add(result, p, modulus, a)
        p = add(p, p, modulus, a)
        k >>= 1
    return result


def curve_points(modulus, a, b):
    return [
        (x, y)
        for x, y in product(range(modulus), repeat=2)
        if (y * y - x * x * x - a * x - b) % modulus == 0
    ]


def polynomial_from_roots(roots, modulus):
    # Coefficients run from constant term upward.
    coefficients = [1]
    for root in roots:
        next_coefficients = [0] * (len(coefficients) + 1)
        for degree, coefficient in enumerate(coefficients):
            next_coefficients[degree] -= root * coefficient
            next_coefficients[degree + 1] += coefficient
        coefficients = [value % modulus for value in next_coefficients]
    return coefficients


def evaluate(coefficients, x, modulus):
    result = 0
    for coefficient in reversed(coefficients):
        result = (result * x + coefficient) % modulus
    return result


def square(coefficients, modulus):
    result = [0] * (2 * len(coefficients) - 1)
    for i, left in enumerate(coefficients):
        for j, right in enumerate(coefficients):
            result[i + j] = (result[i + j] + left * right) % modulus
    return result


def composite_order_example():
    modulus, a, b = 101, 1, 2
    points = curve_points(modulus, a, b)
    assert len(points) + 1 == 100
    p = (6, 27)
    assert multiply(100, p, modulus, a) is None
    assert multiply(50, p, modulus, a) is not None
    assert multiply(20, p, modulus, a) is not None
    q = multiply(17, p, modulus, a)
    assert q == (9, 29)

    d, quotient_order = 4, 25
    kernel_generator = multiply(quotient_order, p, modulus, a)
    assert kernel_generator == (1, 99)
    kernel = [multiply(j, kernel_generator, modulus, a) for j in range(d)]
    labels = [(30, 38), (13, 92), (50, 61), (19, 11)]
    representatives = [
        add(multiply(u, p, modulus, a), multiply(v, q, modulus, a), modulus, a)
        for u, v in labels
    ]
    assert representatives == [(18, 14), (41, 53), (85, 43), (20, 12)]

    # Velu's x-coordinate formula for the order-four kernel. Here Phi_x=A/B.
    def numerator(x):
        return (x**4 - x**3 + 11 * x * x + 9 * x + 12) % modulus

    def denominator(x):
        return (x**3 - x * x - x + 1) % modulus

    # The kernel chain <2K> < <K> yields two degree-two x maps.
    # The intermediate curve has equation y^2=x^3+82x+30.
    def first_isogeny_x(x):
        return ((x * x + x + 4) * pow(x + 1, -1, modulus)) % modulus

    def second_isogeny_x(x):
        return ((x * x - 3 * x + 8) * pow(x - 3, -1, modulus)) % modulus

    def velu_x(point):
        if point in kernel:
            return None
        x = point[0]
        return (
            x
            + sum(add(point, s, modulus, a)[0] - s[0] for s in kernel[1:])
        ) % modulus

    for point in points:
        if point not in kernel:
            x = point[0]
            assert velu_x(point) == (
                numerator(x) * pow(denominator(x), -1, modulus)
            ) % modulus
            assert velu_x(point) == second_isogeny_x(first_isogeny_x(x))

    targets = [velu_x(r) for r in representatives]
    assert targets == [100, 73, 0, 7]
    factor_base = [
        add(r, s, modulus, a) for r in representatives for s in kernel
    ]
    factor_x = {point[0] for point in factor_base}
    polynomial_roots = {
        x
        for x in range(modulus)
        if denominator(x) != 0
        and any(
            (numerator(x) - c * denominator(x)) % modulus == 0
            for c in targets
        )
    }
    # The condition above is L(x) = product_i (A(x)-c_i*B(x)) = 0.
    assert polynomial_roots == factor_x and len(factor_x) == 16

    relation = add(
        multiply(-1, representatives[2], modulus, a),
        add(
            multiply(2, representatives[3], modulus, a),
            multiply(3, kernel_generator, modulus, a),
            modulus,
            a,
        ),
        modulus,
        a,
    )
    assert relation is None
    shifted_r4 = add(representatives[3], multiply(3, kernel_generator, modulus, a), modulus, a)
    assert shifted_r4 == (26, 38)
    assert add(representatives[3], shifted_r4, modulus, a) == representatives[2]
    coefficient_p = (-labels[2][0] + 2 * labels[3][0] + 3 * quotient_order) % 100
    coefficient_q = (-labels[2][1] + 2 * labels[3][1]) % 100
    assert (coefficient_p, coefficient_q) == (63, 61)
    recovered = -coefficient_p * pow(coefficient_q, -1, 100) % 100
    assert recovered == 17
    print("Composite-order example:")
    print("  #E(F_101) = 100; P =", p, "; Q =", q, "; K =", kernel_generator)
    print("  representatives =", representatives)
    print("  tower: (x^2+x+4)/(x+1), then (z^2-3z+8)/(z-3)")
    print("  x(phi(R_i)) =", targets)
    print("  factor-base x roots =", sorted(factor_x))
    print("  factor-base relation: R_4 + (R_4+3K) = R_3")
    print("  point check: (20,12) + (26,38) = (85,43)")
    print("  labels: 63P + 61Q = O; recovered w =", recovered)


def prime_order_orbit_example():
    modulus, a, b = 101, 1, 32
    assert len(curve_points(modulus, a, b)) + 1 == 101
    p = (4, 10)
    q = multiply(3, p, modulus, a)
    roots = [
        multiply(2**j, r, modulus, a)[0] for r in (p, q) for j in range(4)
    ]
    assert roots == [4, 23, 49, 27, 43, 34, 24, 28]
    f = polynomial_from_roots(roots, modulus)
    assert list(reversed(f)) == [1, 71, 66, 38, 25, 97, 29, 34, 3]

    # If f=g(h) with deg h=2, normalize h=x^2+c*x. The x^7
    # coefficient forces c=71/4=43; then f(-c-x)=f(x) is necessary.
    c = 71 * pow(4, -1, modulus) % modulus
    assert c == 43
    assert evaluate(f, 0, modulus) == 3
    assert evaluate(f, -c, modulus) == 93

    # If deg h=4, normalize h=x^4+b3*x^3+b2*x^2+b1*x.
    # The top three coefficients determine b3,b2,b1. Then f-h^2
    # must be a scalar multiple of h plus a constant.
    inv2 = pow(2, -1, modulus)
    b3 = 71 * inv2 % modulus
    b2 = (66 - b3 * b3) * inv2 % modulus
    b1 = (38 - 2 * b3 * b2) * inv2 % modulus
    assert (b3, b2, b1) == (86, 72, 89)
    h = [0, b1, b2, b3, 1]
    h_squared = square(h, modulus)
    remainder = [(f[i] - h_squared[i]) % modulus for i in range(9)]
    assert remainder[4] == 36
    assert remainder[3] == 7
    assert remainder[4] * b3 % modulus == 66 != remainder[3]

    print("Prime-order orbit example:")
    print("  #E(F_101) = 101; P =", p, "; Q =", q)
    print("  orbit x roots =", roots)
    print("  root polynomial coefficients, high to low =", list(reversed(f)))
    print("  degree-2 right component: required symmetry fails (3 != 93)")
    print("  degree-4 right component: x^3 residual is 7, must be 66")


if __name__ == "__main__":
    composite_order_example()
    prime_order_orbit_example()
