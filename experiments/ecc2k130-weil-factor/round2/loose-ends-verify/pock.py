# independent Pocklington/Lehmer full N-1 check with own factoring
def tf(n):
    f={};d=2
    while d*d<=n:
        while n%d==0: f[d]=f.get(d,0)+1; n//=d
        d+=1 if d==2 else 2
        if d>10**6: break
    if n>1: f[n]=f.get(n,0)+1
    return f
def rho(n):
    import math
    if n%2==0: return 2
    for c in range(1,50):
        x=y=2;d=1
        f=lambda v:(v*v+c)%n
        while d==1:
            x=f(x);y=f(f(y));d=math.gcd(abs(x-y),n)
        if d!=n: return d
def factor(n):
    if n==1: return {}
    if isprime(n): return {n:1}
    if n<10**12: return tf(n)
    d=rho(n); a=factor(d); b=factor(n//d)
    for k,v in b.items(): a[k]=a.get(k,0)+v
    return a
def isprime(n):  # Lehmer, recursive
    if n<2: return False
    if n<10**6: return all(n%d for d in range(2,int(n**.5)+1))
    if pow(2,n-1,n)!=1: return False
    F=factor(n-1)
    assert all(isprime(q) for q in F)
    for q in F:
        if not any(pow(a,n-1,n)==1 and pow(a,(n-1)//q,n)!=1 for a in range(2,200)): return False
    return True
p2=146505763881528721
print(isprime(p2), factor(p2-1), isprime(263))
print(isprime(2**61-3), isprime(2**61-1), isprime(p2*263))
