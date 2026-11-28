#!/usr/bin/env python3
"""Lemma (extension fields): for every k with 131 not dividing k, Q(pi^k) = K, so End_{F_{2^k}}(A') = End(A') for every A'
in the class, and there are no new automorphisms (in particular no order-131 automorphism) over F_{2^k}.

Proof. The conjugates of pi^k are alpha^k zeta^{jk}, beta^k zeta^{jk} (j = 1..130). They are pairwise distinct iff
(i) zeta^k != 1, i.e. 131 does not divide k, and (ii) (alpha/beta)^k is not a 131st root of unity. (ii) always holds:
(alpha/beta)^k lies in Q(sqrt-7), whose roots of unity are +-1, and alpha/beta is not a root of unity (|alpha| = |beta|,
alpha != +-beta). So pi^k has 260 distinct conjugates, i.e. degree 260 = [K:Q], and Q(pi^k) = K.

For ordinary varieties, base change to F_{2^k} corresponds to (T, pi^k) in Deligne's category (standard, via the
canonical lift). End_{F_{2^k}} = {x in End^0 commuting with pi^k} intersected with the multiplier ring of T. The
commutant of pi^k is Q(pi^k) = K, so this is the multiplier ring of T, the same as over F_2.

At 131 | k, pi^k = alpha^k is in Q(sqrt-7): A x F_{2^k} ~ E^130 and End^0 = M_130(Q(sqrt-7)).

This script checks distinctness numerically for k = 1..262 (a control, not the proof), and shows that it fails exactly
at 131 | k."""
import cmath, math
a = complex(-0.5, math.sqrt(7)/2); b = a.conjugate()
bad = []
for k in range(1, 263):
    conj = [a**k*cmath.exp(2j*math.pi*j*k/131) for j in range(1,131)] + [b**k*cmath.exp(2j*math.pi*j*k/131) for j in range(1,131)]
    scale = abs(a)**k
    conj.sort(key=lambda z: (round(z.real/scale, 9), round(z.imag/scale, 9)))
    mind = min(abs(conj[i]-conj[i+1]) for i in range(len(conj)-1))
    # also full pairwise on a coarse grid
    distinct = len({(round(z.real/scale, 7), round(z.imag/scale, 7)) for z in conj}) == 260
    if not distinct: bad.append(k)
print("k in 1..262 where the conjugates of pi^k are NOT distinct:", bad)
assert bad == [131, 262], bad
print("OK: distinct exactly when 131 does not divide k (numeric control of the lemma).")
