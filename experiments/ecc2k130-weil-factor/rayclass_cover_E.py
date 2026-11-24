#!/usr/bin/env python3
# CANDIDATE, PENDING RE-REVIEW 2026-11-24: this uses the same class-field-theory mechanism as the
# (now candidate) P^1 ray-class rescue. The computation (|Cl^0_m| coprime to 131) is elementary and
# correct; whether it rigorously excludes every relevant cover (constant-Z/131 vs mu_131 deck group,
# isomorphism vs isogeny) needs independent review. Not a settled claim. See docs CORRECTION section.
"""Does a genus-131 cyclic cover C -> E realize A as its Prym? Settled by class field theory.

A as the Prym of a degree-131 cyclic cover C -> E (E = ECC2K-130 base, genus 1) forces dim Prym = 130,
hence g_C = 131. By function-field class field theory, a geometrically connected cyclic degree-131
cover of E with conductor m corresponds to an order-131 character of the ray class group Cl_m(E) that
is NONTRIVIAL on the degree-0 part Cl^0_m (characters that factor through the degree map Cl_m -> Z give
the disconnected constant-field extension, i.e. Res(E) itself, not a connected curve). Since E has a
degree-1 place, Cl_m = Z (+) Cl^0_m, so such a character exists IFF 131 | |Cl^0_m|.

   |Cl^0_m| = #Pic^0(E)(F_2) * |R_m|,   R_m = (O/m)^* / F_2^*  (the ray part).

#Pic^0(E)(F_2) = #E(F_2) = 4. For a tame modulus (131 coprime to char 2, so all ramification is tame),
R_m = prod over places v|m of (2^deg(v) - 1) times a 2-power (no odd torsion). So 131 | |Cl^0_m| iff some
place v|m has 131 | 2^deg(v) - 1 iff ord_131(2) | deg(v) iff 130 | deg(v).

Riemann-Hurwitz: a degree-131 cover of E of genus 131 has conductor of total degree 2 (2g-2 = 260 =
130*sum deg(v)). Both options -- two degree-1 places (|Cl^0_m| = 4) and one degree-2 place
(|Cl^0_m| = 12) -- have |Cl^0_m| coprime to 131. So NO connected cyclic degree-131 cover of E has
genus 131: A is NOT the Prym of such a cover. 131-torsion in the ray class group first appears when the
conductor contains a degree-130 place, forcing genus >= 8451 -- astronomically outside the window.

(A connected degree-131 cover is Galois iff cyclic since 131 is prime; A's zeta_131-eigenvalue structure
requires the cyclic one. The fully general question -- A a Jacobian of a curve with no cover structure --
is separate and remains open.)
"""
def ord_mod(a, n):
    e, x = 1, a % n
    while x != 1:
        x = (x*a) % n; e += 1
    return e

def num_E_F2(a, b):
    pts = 1
    for x in (0, 1):
        for y in (0, 1):
            if (y*y ^ x*y) == (x*x*x ^ (a*x*x) ^ b): pts += 1
    return pts

def ray_part_order(place_degs):           # tame, multiplicity 1; odd part is prod (2^d - 1)
    o = 1
    for d in place_degs: o *= (2**d - 1)
    return o

def main():
    ell = 131
    o2 = ord_mod(2, ell)
    hE = num_E_F2(0, 1)
    print(f"ord_131(2) = {o2}  => 131 | 2^d - 1 iff {o2} | d; smallest such d = {o2}")
    print(f"#Pic^0(E)(F_2) = #E(F_2) = {hE}; 131 | {hE}? {hE % ell == 0}")
    print("genus-131 cover: conductor total degree 2. |Cl^0_m| = #E(F2) * ray part:")
    cases = {"two degree-1 places": [1, 1], "one degree-2 place": [2]}
    any_131 = False
    for name, degs in cases.items():
        Cl0 = hE * ray_part_order(degs)
        div = (Cl0 % ell == 0)
        any_131 = any_131 or div
        print(f"  {name:22s}: |Cl^0_m| = {Cl0:3d}; 131 | |Cl^0_m|? {div}")
    g_force = 130*o2//2 + 1
    print(f"\n131-torsion needs a degree-{o2} place in the conductor; a cyclic 131-cover totally ramified")
    print(f"there has genus 2g-2 = 130*{o2} = {130*o2}, i.e. g = {g_force} >> window (~130-300).")
    assert not any_131 and hE % ell != 0, "unexpected 131-torsion"
    print("\nVERDICT: no connected cyclic degree-131 cover of E has genus 131; A is NOT its Prym."
          " The cover-of-E route is closed (genus forced to >= 8451 to carry the needed 131-structure).")

if __name__ == "__main__":
    main()
