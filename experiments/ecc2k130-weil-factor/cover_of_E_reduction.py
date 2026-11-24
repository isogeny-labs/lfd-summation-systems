#!/usr/bin/env python3
"""The cover-of-E / Prym route for A, and its reduction.

A degree-131 cyclic cover C -> E (E the ECC2K-130 base curve, genus 1) totally ramified at 2 points has
genus 131 (Riemann-Hurwitz: 2g-2 = 131*(2*1-2) + 2*(131-1) = 260). Then Jac(C) ~ E x Prym with
dim Prym = 130. For Prym ~ A we need Jac(C) ~ E x A = W, where W = Res_{F_2^131/F_2}(E) has
characteristic polynomial P_W(x) = (x^2+x+2) P_A(x) = x^262 - s_131 x^131 + 2^131.

So: (A is the Prym of a degree-131 cover of E) ==> (W is isogenous to a genus-131 Jacobian carrying a
degree-131 map to E). This script tests the Howe-Lauter place-count necessary conditions for the
isogeny class of W (what such a curve C would have to look like). Unlike the P^1-cyclic-cover case, the
sqrt(-7) field obstruction does NOT apply here (the base E already has CM by sqrt(-7)), so this route is
not closed by that argument.

Result: W passes all place-count tests, and strikingly #C(F_2^k) = 2^k + 1 for k = 1..130 (the virtual
curve looks exactly like P^1 up to degree 130, revealing genus only at degree 131). So the cover-of-E
route is OPEN: it is not ruled out by the cheap necessary conditions.
"""
def main():
    t = -1
    s = [2, t]
    for k in range(2, 2000): s.append(t*s[-1] - 2*s[-2])      # s_k^E = trace of Frob^k on E
    # power sums of P_W roots: s_k^W = s_k^E * (131 if 131|k else 0)
    sW = lambda k: 131*s[k] if k % 131 == 0 else 0
    C = lambda k: 2**k + 1 - sW(k)                            # virtual #C(F_2^k)
    def mu(m):
        if m == 1: return 1
        r, mm, p = 1, m, 2
        while p*p <= mm:
            if mm % p == 0:
                mm //= p
                if mm % p == 0: return 0
                r = -r
            p += 1
        if mm > 1: r = -r
        return r
    neg = []
    for d in range(1, 400):
        a = sum(mu(d//e)*C(e) for e in range(1, d+1) if d % e == 0)
        assert a % d == 0
        a //= d
        if a < 0: neg.append((d, a))
    p1like = all(C(k) == 2**k + 1 for k in range(1, 131))
    print("W char poly degree 262 (genus 131 if a Jacobian).")
    print("place counts a_1..a_6:", [sum(mu(d//e)*C(e) for e in range(1, d+1) if d % e == 0)//d
                                      for d in range(1, 7)])
    print("#C(F_2^k) = 2^k + 1 for all k=1..130 (P^1-like up to degree 130):", p1like)
    print("negative place counts a_d, d=1..399:", neg if neg else "NONE")
    print("VERDICT: W passes the place-count necessary conditions; the cover-of-E route is NOT closed by"
          " them. Whether W is actually a genus-131 Jacobian (hence whether the cover exists) is open.")
    assert not neg and p1like

if __name__ == "__main__":
    main()
