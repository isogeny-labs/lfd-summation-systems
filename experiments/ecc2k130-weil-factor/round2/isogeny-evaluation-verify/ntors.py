# independent ord_N(2): factor N-1 partially via PARI-free trial + check (N-1)/24 exponent candidates
import subprocess
t=[2,-1]
for k in range(2,132): t.append(-t[-1]-2*t[-2])
nE=2**131+1-t[131]; N=nE//4; assert nE%4==0
print("N=",N, "nE%4",nE%4)
o=(N-1)//24; print("2^((N-1)/24) mod N == 1:", pow(2,o,N)==1)
for q in (2,3):
    print("2^((N-1)/%d) == 1:"%(24*q), (N-1)%(24*q)==0 and pow(2,(N-1)//(24*q),N)==1)
