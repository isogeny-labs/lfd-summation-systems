import sys, numpy as np
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ct_lib import build, s_of
from collections import Counter
def tables(n):
    q,g,dlog,mul,pw=build(n)
    exp=np.zeros(2*q,dtype=np.int64); x=1
    for i in range(q-1): exp[i]=x; x=mul(x,g)
    exp[q-1:2*q-2]=exp[0:q-1]
    return q,g,dlog,exp
def vmul_scalar(c,X,dlog,exp,q):
    # c*X elementwise (c scalar field elt, X array)
    out=np.zeros_like(X)
    if c==0: return out
    nz=X!=0
    out[nz]=exp[(dlog[c]+dlog[X[nz]])%(q-1)]
    return out
def tally(hist_rows, s, eps, coh):
    # hist_rows: (N, ell) histograms; add eps at col 0
    h=hist_rows.copy(); h[:,0]+=eps
    srt=np.sort(h,axis=1)
    lowone=(srt[:,1]==srt[:,-1]); highone=(srt[:,0]==srt[:,-2])
    for i in np.nonzero(lowone|highone)[0]:
        delta = int(srt[i,0]-srt[i,1]) if lowone[i] else int(srt[i,-1]-srt[i,0])
        coh[-delta]+=1
        assert abs(delta) <= 2*int(round((h[i].sum()+0.0)**0.5))+2, ("Weil bound violated", delta)
def hist2d(D,mask,ell):
    # D,mask shape (N,q): per-row bincount of D[mask]
    N=D.shape[0]; out=np.zeros((N,ell),dtype=np.int64)
    off=(np.arange(N)[:,None]*ell + D)
    out_flat=np.bincount(off[mask],minlength=N*ell)
    return out_flat.reshape(N,ell)
def run(ell,n,which,sample=None,seed=1):
    q,g,dlog,exp=tables(n); s=s_of(ell)
    X=np.arange(q,dtype=np.int64)
    X2=np.array([0]+[exp[(2*dlog[x])%(q-1)] for x in range(1,q)],dtype=np.int64)
    X3=np.array([0]+[exp[(3*dlog[x])%(q-1)] for x in range(1,q)],dtype=np.int64)
    rng=np.random.default_rng(seed)
    coh=Counter(); count=0
    V=np.arange(q,dtype=np.int64)[:,None]
    if which=='ii':   # branch 0, inf, pair roots of x^2+x+n ; f=x * (x^2+x+n)^a3, eps=0
        base=(X2^X)[None,:]^V              # x^2+x+v, rows v
        Dq=dlog[base]; Dx=dlog[X]
        irr=(Dq>=0).all(axis=1)            # no root in F_q
        rows=np.nonzero(irr)[0]
        for a3 in range(1,ell):
            if (-(1+2*a3))%ell==0: continue
            D=(Dx[None,:]+a3*Dq[rows])%ell
            mask=(Dx[None,:]>=0)&np.ones((len(rows),1),bool)
            H=hist2d(D,mask,ell); tally(H,s,0,coh); count+=len(rows)
    elif which=='iv': # branch inf + cubic orbit; f=(x^3+c1 x+c0), eps=0
        us=range(q) if sample is None else rng.choice(q,size=sample,replace=False)
        for c1 in us:
            base=(X3^vmul_scalar(int(c1),X,dlog,exp,q))[None,:]^V   # rows c0
            Dc=dlog[base]; irr=(Dc>=0).all(axis=1); rows=np.nonzero(irr)[0]
            if len(rows)==0: continue
            D=Dc[rows]%ell; mask=np.ones_like(D,bool)
            H=hist2d(D,mask,ell); tally(H,s,0,coh); count+=len(rows)
    elif which=='iii': # two pairs; f=(x^2+x+n1)*(x^2+u x+v)^(ell-1), eps=1
        # representative n1 with x^2+x+n1 irreducible
        Dq1=None
        for n1 in range(q):
            dq=dlog[(X2^X)^n1]
            if (dq>=0).all(): Dq1=dq; N1=n1; break
        us=range(1,q) if sample is None else rng.choice(np.arange(1,q),size=sample,replace=False)
        for u in us:
            base=(X2^vmul_scalar(int(u),X,dlog,exp,q))[None,:]^V
            Dq2=dlog[base]; irr=(Dq2>=0).all(axis=1)
            if int(u)==1: irr[N1]=False            # q2 == q1 gives f = q1^ell: not a curve
            rows=np.nonzero(irr)[0]
            D=(Dq1[None,:]+(ell-1)*Dq2[rows])%ell; mask=np.ones_like(D,bool)
            H=hist2d(D,mask,ell); tally(H,s,1,coh); count+=len(rows)
    return s,count,coh
if __name__=="__main__":
    ell,n,which=int(sys.argv[1]),int(sys.argv[2]),sys.argv[3]
    sample=int(sys.argv[4]) if len(sys.argv)>4 else None
    s,count,coh=run(ell,n,which,sample)
    print(f"ell={ell} config={which}{' sample u='+str(sample) if sample else ''}: covers={count}  A needs {s}  "
          f"coherent trace ints={dict(sorted(coh.items()))}  MATCH={coh.get(s,0)}")
