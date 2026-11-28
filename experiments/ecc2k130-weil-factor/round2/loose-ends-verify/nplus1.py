from math import gcd
N=146505763881528721; P,Q=1,-29; D=P*P-4*Q
def jac(a,n):
    a%=n;r=1
    while a:
        while a%2==0:
            a//=2
            if n%8 in(3,5): r=-r
        a,n=n,a
        if a%4==3 and n%4==3: r=-r
        a%=n
    return r if n==1 else 0
def U(k):  # matrix power
    def mul(A,B): return [[(A[0][0]*B[0][0]+A[0][1]*B[1][0])%N,(A[0][0]*B[0][1]+A[0][1]*B[1][1])%N],[(A[1][0]*B[0][0]+A[1][1]*B[1][0])%N,(A[1][0]*B[0][1]+A[1][1]*B[1][1])%N]]
    R=[[1,0],[0,1]];M=[[P,-Q],[1,0]]
    while k:
        if k&1:R=mul(R,M)
        M=mul(M,M);k>>=1
    return R[1][0]
qs=[2,11,131,1777,28607065273]
prod=1
for q in qs: prod*=q
print("factorization ok",prod==N+1, "jacobi",jac(D,N))
print("U_{N+1}=0", U(N+1)==0, [gcd(U((N+1)//q),N) for q in qs])
def isp(n): return n>1 and all(n%d for d in range(2,int(n**.5)+1))
print("q primes", [isp(q) for q in qs])
