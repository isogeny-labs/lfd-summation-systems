"""Independent membership replay for compact derivations, not closedness."""
import hashlib,json
from pathlib import Path
from replay import read_input

def verify(inp,proof):
    inp=Path(inp);proof=json.loads(Path(proof).read_text())
    if hashlib.sha256(inp.read_bytes()).hexdigest()!=proof['input_sha256']:raise ValueError('input digest')
    (p,m,k,a,D),gens=read_input(inp);r=m*k;N=2**k;basis=[]
    def unpack(code):
        if code<0 or code>>(r+6*a):raise ValueError('monomial range')
        xs=[(code//2**(j*k))%N for j in range(m)]
        es=[(code//2**(r+6*j))%64 for j in range(a)]
        return xs,es
    def degree(code):
        xs,es=unpack(code);return sum(x.bit_count() for x in xs)+sum(es)
    def leading(poly):return max(poly,key=lambda x:(degree(x),x))
    def multiply(poly,mon):
        ys,fs=unpack(mon);out={}
        if max(map(degree,poly),default=0)+degree(mon)>D:raise ValueError('degree budget')
        for code,c in poly.items():
            xs,es=unpack(code);dest=0
            for j,(x,y) in enumerate(zip(xs,ys)):dest+=((x+y)%N)*2**(j*k)
            for j,(x,y) in enumerate(zip(es,fs)):
                if x+y>=64:raise ValueError('auxiliary exponent')
                dest+=(x+y)*2**(r+6*j)
            out[dest]=(out.get(dest,0)+c)%p
        return {x:c for x,c in out.items() if c}
    def lcm_and_multipliers(x,y):
        mx=x%2**r;my=y%2**r;mask=mx|my
        _,ex=unpack(x);_,ey=unpack(y)
        lcm=mask+sum(max(u,v)*2**(r+6*j) for j,(u,v) in enumerate(zip(ex,ey)))
        hx=(mask^mx)+sum((max(u,v)-u)*2**(r+6*j) for j,(u,v) in enumerate(zip(ex,ey)))
        hy=(mask^my)+sum((max(u,v)-v)*2**(r+6*j) for j,(u,v) in enumerate(zip(ex,ey)))
        if degree(lcm)>D:raise ValueError('pair degree budget')
        return hx,hy
    def subtract(row,poly,c):
        for x,v in poly.items():row[x]=(row.get(x,0)-c*v)%p
        return {x:v for x,v in row.items() if v}
    for ident,node in enumerate(proof['trace']):
        src=node['source'];kind=src[0]
        if kind=='input':
            i=src[1]
            if not 0<=i<len(gens):raise ValueError('input reference')
            row=dict(gens[i])
            if max(map(degree,row),default=0)>D:raise ValueError('input degree')
        elif kind=='carry':
            i,v=src[1:]
            if not 0<=i<ident or not 0<=v<r:raise ValueError('carry reference')
            row=multiply(basis[i],2**v)
        elif kind=='pair':
            i,j=src[1:]
            if not 0<=i<ident or not 0<=j<ident:raise ValueError('pair reference')
            hi,hj=lcm_and_multipliers(leading(basis[i]),leading(basis[j]))
            row=subtract(multiply(basis[i],hi),multiply(basis[j],hj),1)
        else:raise ValueError('source kind')
        for j,h,c in node['reductions']:
            if not 0<=j<ident or not 0<c<p:raise ValueError('reducer reference')
            row=subtract(row,multiply(basis[j],h),c)
        inv=node['inverse']
        if not 0<inv<p:raise ValueError('normalization')
        row={x:c*inv%p for x,c in row.items() if c}
        pairs=node['row'];target=dict(pairs)
        if len(pairs)!=len(target) or target!=row or not row:raise ValueError('row identity')
        if row[leading(row)]!=1:raise ValueError('monic row')
        basis.append(row)
    return {'status':'passed','verified_rows':len(basis),'contains_one':any(row=={0:1} for row in basis),'closed_check':False,'scope':'row membership only; no compact closedness claim'}
