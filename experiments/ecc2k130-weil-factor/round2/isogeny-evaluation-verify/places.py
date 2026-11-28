#!/usr/bin/env python3
"""Necessary condition for some C/F_2 with Jac(C) ~ A: degree-d place counts a_d >= 0.
#C(F_{2^k}) = 2^k + 1 - t_k * c_k, t_k = alpha^k+beta^k, c_k = -1 (131 !| k) or 130.
Control: E itself (genus 1) gives #E(F_2^k) counts; a fake 'negative' trace vector is rejected."""
import math
def mob(n):
    r,m,p=1,n,2
    while p*p<=m:
        if m%p==0:
            m//=p
            if m%p==0: return 0
            r=-r
        p+=1
    return -r if m>1 else r
t=[2,-1]
for k in range(2,600): t.append(-t[-1]-2*t[-2])
def Ncount(k): return 2**k+1-t[k]*(130 if k%131==0 else -1)
neg=[];mins=None
for d in range(1,400):
    a=sum(mob(d//e)*Ncount(e) for e in range(1,d+1) if d%e==0)
    assert a%d==0; a//=d
    if a<0: neg.append(d)
    if d<=12: print("d=%d places=%d"%(d,a))
print("degrees d<400 with negative place count:",neg)
# Ihara bound check for #C(F_2), #C(F_4): N_1 <= q+1+ (sqrt((8q+1)g^2+4(q^2-q)g)-g)/2
g=130;q=2
print("#C(F_2)=%d  Ihara bound=%.1f ; #C(F_4)=%d bound(q=4)=%.1f"%(Ncount(1),q+1+(math.sqrt((8*q+1)*g*g+4*(q*q-q)*g)-g)/2,Ncount(2),5+(math.sqrt(33*g*g+48*g)-g)/2))
