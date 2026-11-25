#!/usr/bin/env python3
"""Independent symbolic checks of the saved integral-symmetry certificate."""
import json
from pathlib import Path

import sympy as s

here = Path(__file__).parent
record = json.loads((here / "integral_symmetry.json").read_text())
coeffs = json.loads((here / "P_A.json").read_text())["P_A_coeffs_low_to_high"]
x, a = s.symbols("x a")
ell = record["howe_theorem_1_1_witness"]["ell"]
assert s.factorint(record["alpha_power_131"]["alpha_coefficient"]) == {263: 1, ell: 1}
assert s.isprime(ell) and s.legendre_symbol(-7, ell) == -1
assert s.n_order(2, 131) == 130
f = s.Poly.from_list(list(reversed(coeffs)), x, modulus=ell)
pi = record["howe_theorem_1_1_witness"]["pi_equals_V_coordinates"][0]
assert s.gcd(f, s.Poly(x*x - 2, x, modulus=ell)) == s.Poly(x - pi, x, modulus=ell)
z0, z1 = record["howe_theorem_1_1_witness"]["zeta_coordinates"]
z = s.Poly(z0 + z1*a, a, modulus=ell)
quadratic = s.Poly(a*a + a + 2, a, modulus=ell)
assert s.rem(z**131, quadratic) == s.Poly(1, a, modulus=ell)
assert s.rem(s.Poly(a, a, modulus=ell)*z, quadratic) == s.Poly(pi, a, modulus=ell)

# Newton sums from the saved degree-260 polynomial, independently of the
# structural elliptic trace formula used by the producer.
descending = list(reversed(coeffs))
power_sums = [260]
for k in range(1, 32):
    power_sums.append(-k*descending[k] - sum(descending[i]*power_sums[k-i]
                                           for i in range(1, k)))
points = {k: 2**k + 1 - power_sums[k] for k in range(1, 32)}
for d in range(1, 32):
    places = sum(s.mobius(d//e)*points[e] for e in s.divisors(d))/d
    assert places.is_integer and places >= 0
assert min(points.values()) == 2

# Rational-series recurrence, independently of the producer's coefficient sum.
effective = []
for m in range(261):
    val = descending[m]
    if m >= 1:
        val += 3*effective[m-1]
    if m >= 2:
        val -= 2*effective[m-2]
    effective.append(val)
assert min(effective) >= 0
assert effective[259] == sum(coeffs)*(2**130 - 1)
assert effective[260] == sum(coeffs)*(2**131 - 1)
print("Independent SymPy checks passed: factorization, field witness, polynomial gcd,")
print("Newton/Mobius point and place counts, and effective-divisor recurrence.")
