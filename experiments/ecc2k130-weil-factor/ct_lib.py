import numpy as np, sys
PRIM={4:0b10011,10:0b10000001001,12:0b1000001010011,18:0b1000000000010000001,28:0b10000000000000000000000001001}
def build(n):
    mod=PRIM[n]; q=1<<n
    def mul(a,b):
        r=0
        while b:
            if b&1: r^=a
            b>>=1; a<<=1
            if a>>n&1: a^=mod
        return r
    # generator test via factorization of q-1
    m=q-1; fs=[]; d=2
    while d*d<=m:
        if m%d==0:
            fs.append(d)
            while m%d==0: m//=d
        d+=1
    if m>1: fs.append(m)
    def pw(a,e):
        r=1
        while e:
            if e&1: r=mul(r,a)
            a=mul(a,a); e>>=1
        return r
    for g in range(2,q):
        if all(pw(g,(q-1)//p)!=1 for p in fs): break
    dlog=np.full(q,-1,dtype=np.int64); x=1
    for i in range(q-1): dlog[x]=i; x=mul(x,g)
    return q,g,dlog,mul,pw
def s_of(ell):
    s=[2,-1]
    for k in range(2,2*ell): s.append(-s[-1]-2*s[-2])
    return s[ell-1]
def classify(hp, s):
    # hp: corrected histogram h' (length ell). match iff constant except one entry offset by -s.
    srt=np.sort(hp)
    if s>0: match = srt[1]==srt[-1] and srt[0]==srt[1]-s
    elif s<0: match = srt[0]==srt[-2] and srt[-1]==srt[0]-s
    else: match = srt[0]==srt[-1]
    flat1 = (srt[1]==srt[-1]) or (srt[0]==srt[-2])   # "constant except one", any offset
    return bool(match), bool(flat1), int(srt[-1]-srt[0])
