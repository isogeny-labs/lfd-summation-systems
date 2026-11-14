"""Deterministic finite-field sanity checks for N8's Semaev secant and Wronskian inputs.

Not a mathematical proof and not a Lean build; explicit counts in JSON output.
"""
from __future__ import annotations

from math import prod
from random import Random
import json


def s3(a: int, b: int, x: int, z: int, t: int, p: int) -> int:
    return ((x-z)**2*t**2 - 2*((x+z)*(x*z+a)+2*b)*t
            + (x*z-a)**2 - 4*b*(x+z)) % p


def check_secant(p: int, trials: int = 5000) -> dict[str, int]:
    rng = Random(110012 + p)
    sqrt = {(v*v) % p: v for v in range(p)}
    accepted = 0
    for _ in range(trials):
        a, b = rng.randrange(p), rng.randrange(p)
        if (4*a*a*a + 27*b*b) % p == 0:
            continue
        x, t = rng.sample(range(p), 2)
        y = sqrt.get((x*x*x+a*x+b) % p)
        Y = sqrt.get((t*t*t+a*t+b) % p)
        if y is None or Y is None:
            continue
        A = ((x+t)*(x*t+a)+2*b) % p
        z = ((A+2*y*Y)*pow((x-t) % p, -2, p)) % p
        assert s3(a,b,x,z,t,p) == 0, (p,a,b,x,t,y,Y,z)
        accepted += 1
    return {"trials":trials,"valid_point_pairs":accepted}


def check_charp_monomial_wronskian(p: int, trials: int = 3000) -> dict[str, int]:
    """For X^a over F_p(X^p), independent iff exponent residues are distinct.

    Ordinary Wronskian determinant coefficient is a Vandermonde in exponents.
    This checks its characteristic-p behavior, not the N8 function-field estimate.
    """
    rng=Random(220024+p)
    dependent=0
    for _ in range(trials):
        r=rng.randrange(1,min(4,p)+1)
        exponents=[rng.randrange(10*p) for _ in range(r)]
        independent=len({v%p for v in exponents})==r
        vandermonde=prod(exponents[j]-exponents[i] for i in range(r)
                        for j in range(i+1,r))%p
        assert (vandermonde != 0) == independent, (p,exponents)
        dependent += not independent
    return {"trials":trials,"dependent":dependent}


def main() -> None:
    results={"secant":{},"wronskian_monomials":{}}
    for p in (5,7,11,61,101,131,251):
        results["secant"][str(p)] = check_secant(p)
    for p in (3,5,7,11,13):
        results["wronskian_monomials"][str(p)] = check_charp_monomial_wronskian(p)
    print(json.dumps(results,indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
