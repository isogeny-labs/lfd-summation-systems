"""Independent integer replay of sparse-producer row derivations.

Imports no producer algebra. A ledger proves containment of each supplied
row in V_D. Exactness additionally needs a closedness check; a stopped
partial basis is never treated as a closure. Output SHA is supplied by caller.
"""
import hashlib,math
from pathlib import Path


def read_input(path):
    data=list(map(int,Path(path).read_text().split()));p,m,k,a,D,seconds,cap,terms=data[:8]
    j=8;count=data[j];j+=1;gens=[]
    for _ in range(count):
        n=data[j];j+=1;g={}
        for _ in range(n):
            x,c=data[j:j+2];j+=2
            if x in g or not 0<c<p:raise ValueError('generator encoding')
            g[x]=c
        gens.append(g)
    if j!=len(data):raise ValueError('trailing input')
    if p<=3 or any(p%d==0 for d in range(2,math.isqrt(p)+1)) or (p-1)%(1<<k):raise ValueError('field/chains')
    return (p,m,k,a,D),gens


def replay(input_path,ledger_path,basis_path,expected_basis_sha,closed=False):
    if hashlib.sha256(Path(basis_path).read_bytes()).hexdigest()!=expected_basis_sha:raise ValueError('basis target digest')
    (p,m,k,a,D),gens=read_input(input_path);r=m*k;N=1<<k
    def degree(x):return (x%2**r).bit_count()+sum((x>>(r+6*j))%64 for j in range(a))
    def shift(row,v):
        result={}
        for code,c in row.items():
            chain=code%2**r;aux=code>>r
            if v<r:
                b,j=divmod(v,k);exps=[(chain>>(i*k))%N for i in range(m)]
                exps[b]=(exps[b]+2**j)%N
                dest=sum(e<<(i*k) for i,e in enumerate(exps))+(aux<<r)
            else:dest=code+(1<<(r+6*(v-r)))
            result[dest]=(result.get(dest,0)+c)%p
        return {x:c for x,c in result.items() if c}
    basis=[];pivots={}
    with open(ledger_path) as logs,open(basis_path) as rows:
        for ident,line in enumerate(logs):
            parts=line.split();id_,op,source,var,inv,count=parts[:6];id_,source,var,inv,count=map(int,[id_,source,var,inv,count])
            if id_!=ident or not 0<inv<p or len(parts)!=6+2*count:raise ValueError('ledger framing')
            if op=='input':
                if not 0<=source<len(gens):raise ValueError('input reference')
                row=dict(gens[source])
            elif op=='shift':
                if not 0<=source<ident or not 0<=var<r+a:raise ValueError('shift reference')
                if max(map(degree,basis[source]),default=0)>=D:raise ValueError('inadmissible source degree')
                row=shift(basis[source],var)
            else:raise ValueError('operation')
            if max(map(degree,row),default=0)>D:raise ValueError('input degree')
            nums=list(map(int,parts[6:]))
            for j in range(0,len(nums),2):
                old,c=nums[j:j+2]
                if not 0<=old<ident or not 0<c<p:raise ValueError('reducer reference')
                for x,v in basis[old].items():row[x]=(row.get(x,0)-c*v)%p
            row={x:c*inv%p for x,c in row.items() if c%p}
            target=list(map(int,rows.readline().split()))
            if len(target)<2 or target[0]!=ident or len(target)!=2+2*target[1]:raise ValueError('basis framing')
            pairs=list(zip(target[2::2],target[3::2]));out=dict(pairs)
            if len(out)!=len(pairs) or out!=row or not row:raise ValueError('derived row mismatch')
            lead=max(row,key=lambda x:(degree(x),x))
            if row[lead]!=1 or lead in pivots:raise ValueError('echelon pivot')
            pivots[lead]=ident;basis.append(row)
        if rows.read().strip():raise ValueError('extra basis rows')
    def reduce(row):
        row=dict(row)
        while row:
            x=max(row,key=lambda t:(degree(t),t))
            if x not in pivots:return row
            c=row[x]
            for y,v in basis[pivots[x]].items():row[y]=(row.get(y,0)-c*v)%p
            row={y:v for y,v in row.items() if v}
        return row
    exclusions=None
    if closed:
        if any(reduce(g) for g in gens if max(map(degree,g),default=0)<=D):raise ValueError('missing input')
        for row in basis:
            if max(map(degree,row))>=D:continue
            for v in range(r+a):
                if reduce(shift(row,v)):raise ValueError('space not closed')
        exclusions=not any(0==max(row,key=lambda x:(degree(x),x)) for row in basis)
    return dict(derived_rank=len(basis),contains_one=0 in pivots,closed_check=closed,one_excluded=exclusions,
                scope='derived rows contained in V_D'+('; exact closure checked' if closed else '; no completeness claim'))
