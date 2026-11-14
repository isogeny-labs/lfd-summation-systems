#!/usr/bin/env python3
"""Independent exact Q(X) inverse, quotient and product differentiation tests."""
from fractions import Fraction as F
from random import Random

def norm(a):
    t=list(a)
    while len(t)>1 and t[-1]==0: t.pop()
    return tuple(t)
def add(a,b):
    return norm(tuple((a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0)
                     for i in range(max(len(a),len(b)))))
def neg(a): return tuple(-t for t in a)
def sub(a,b): return add(a,neg(b))
def mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j]+=x*y
    return norm(out)
def der(a): return norm(tuple(i*a[i] for i in range(1,len(a))) or (F(0),))
def power(a,n):
    out=(F(1),)
    for _ in range(n):out=mul(out,a)
    return out
class RF:
    def __init__(self,n,d):
        self.n,self.d=norm(n),norm(d)
        assert self.d!=(F(0),)
    def __mul__(self,o): return RF(mul(self.n,o.n),mul(self.d,o.d))
    def __truediv__(self,o):
        assert o.n!=(F(0),)
        return RF(mul(self.n,o.d),mul(self.d,o.n))
    def __neg__(self): return RF(neg(self.n),self.d)
    def __add__(self,o): return RF(add(mul(self.n,o.d),mul(o.n,self.d)),mul(self.d,o.d))
    def __sub__(self,o): return self+(-o)
    def __eq__(self,o): return mul(self.n,o.d)==mul(o.n,self.d)
    def derivative(self):
        return RF(sub(mul(der(self.n),self.d),mul(self.n,der(self.d))),power(self.d,2))
    def inverse(self):
        assert self.n!=(F(0),)
        return RF(self.d,self.n)

r=Random(20261113)
def random_poly():
    return norm(tuple(F(r.randint(-5,5),r.choice([1,1,2,3,5]))
                      for _ in range(r.randint(1,5))))
def random_nonzero_poly():
    while True:
        q=random_poly()
        if q!=(F(0),):return q

for _ in range(4096):
    f=RF(random_nonzero_poly(),random_nonzero_poly())
    g=RF(random_nonzero_poly(),random_nonzero_poly())
    assert f.inverse().derivative()==-(f.derivative()/(f*f))
    assert (f/g).derivative()==(f.derivative()*g-f*g.derivative())/(g*g)
    assert (f*g).derivative()==f*g.derivative()+g*f.derivative()
print("OK: 4096 exact Q(X) pairs; 12288 inverse, quotient, product identities")
