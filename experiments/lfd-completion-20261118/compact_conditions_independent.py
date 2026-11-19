#!/usr/bin/env python3
"""Independent finite-condition checks. The bounded criterion remains draft."""
import argparse,functools,hashlib,json,time
from pathlib import Path
from replay import read_input
HERE=Path(__file__).resolve().parent

def check(run_dir,seconds=30,snapshot=None):
    d=Path(run_dir);inp=d/'input.txt';raw=Path(snapshot) if snapshot else d/'candidate.json';data=json.loads(raw.read_text());input_sha=hashlib.sha256(inp.read_bytes()).hexdigest()
    if input_sha!=data['input_sha256']:raise ValueError('input digest')
    (p,m,k,a,D),gens=read_input(inp);r=m*k;n=r+a;end=time.monotonic()+seconds
    @functools.lru_cache(maxsize=262144)
    def powers(code):
        if code<0 or code>>(r+6*a):raise ValueError('monomial range')
        return tuple((code>>j)&1 for j in range(r))+tuple((code>>(r+6*j))%64 for j in range(a))
    def degree(code):return sum(powers(code))
    def pack(es):return sum(es[j]<<j for j in range(r))+sum(es[r+j]<<(r+6*j) for j in range(a))
    def multiplier(x,y):
        ex,ey=powers(x),powers(y)
        if any(v<w for v,w in zip(ex,ey)):return None
        return pack(tuple(v-w for v,w in zip(ex,ey)))
    def lcm(x,y):return pack(tuple(max(v,w) for v,w in zip(powers(x),powers(y))))
    def times(row,code):
        if max(map(degree,row),default=0)+degree(code)>D:raise ValueError('product degree')
        ec=powers(code);out={}
        for x,c in row.items():
            ex=powers(x);dest=0
            for b in range(m):
                u=sum(ex[b*k+j]<<j for j in range(k));v=sum(ec[b*k+j]<<j for j in range(k));dest+=((u+v)%(1<<k))<<(b*k)
            for j in range(a):
                e=ex[r+j]+ec[r+j]
                if e>=64:raise ValueError('auxiliary overflow')
                dest+=e<<(r+6*j)
            out[dest]=(out.get(dest,0)+c)%p
        return {x:c for x,c in out.items() if c}
    def leading(row):return max(row,key=lambda x:(degree(x),x))
    G=[];leads=[]
    for pairs in data['polynomials']:
        row=dict(pairs)
        if not row or len(row)!=len(pairs) or any(not 0<c<p or degree(x)>D for x,c in row.items()):raise ValueError('candidate polynomial')
        L=leading(row);inv=pow(row[L],-1,p);G.append({x:c*inv%p for x,c in row.items()});leads.append(L)
    def normal_form(poly):
        row={x:c%p for x,c in poly.items() if c%p};rem={}
        while row:
            if time.monotonic()>end:raise TimeoutError
            x=leading(row);c=row[x]
            for j,L in enumerate(leads):
                h=multiplier(x,L)
                if h is None:continue
                product=times(G[j],h)
                if leading(product)!=x or product[x]!=1:raise ValueError('leading product')
                for y,v in product.items():
                    value=(row.get(y,0)-c*v)%p
                    if value:row[y]=value
                    else:row.pop(y,None)
                break
            else:rem[x]=row.pop(x)
        return rem
    counts={'inputs':0,'carries':0,'pairs':0};missing=None
    def test(poly,kind,*ids):
        nonlocal missing
        rem=normal_form(poly);counts[kind]+=1
        if rem:missing={'kind':kind,'source':ids,'remainder':sorted(rem.items())};return False
        return True
    try:
        for i,g in enumerate(gens):
            if max(map(degree,g),default=0)<=D and not test(g,'inputs',i):break
        if missing is None:
            for i,g in enumerate(G):
                if degree(leads[i])<D:
                    for v in range(r):
                        if powers(leads[i])[v] and not test(times(g,1<<v),'carries',i,v):break
                if missing:break
        if missing is None:
            for i,g in enumerate(G):
                for j in range(i):
                    if time.monotonic()>end:raise TimeoutError
                    L=lcm(leads[i],leads[j])
                    if degree(L)>D:continue
                    left=times(g,multiplier(L,leads[i]));right=times(G[j],multiplier(L,leads[j]))
                    for x,c in right.items():left[x]=(left.get(x,0)-c)%p
                    if not test({x:c for x,c in left.items() if c},'pairs',i,j):break
                if missing:break
        status='conditions-failed' if missing else 'finite-conditions-passed'
    except TimeoutError:status='capped'
    result={'status':status,'scope':'independent finite conditions only; bounded-criterion proof remains draft','closed_check':False,'checked':counts,'missing':missing,'contains_constant_candidate':0 in leads,'input_sha256':input_sha,'candidate_sha256':hashlib.sha256(raw.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'results'/(d.name+'-conditions-independent.json')).write_text(json.dumps(result,indent=2)+'\n');return result
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('run_dir',type=Path);ap.add_argument('--snapshot',type=Path);ap.add_argument('--seconds',type=int,default=30);a=ap.parse_args();out=check(a.run_dir,a.seconds,a.snapshot)
    print(json.dumps({'status':out['status'],'checked':out['checked'],'missing_source':None if out['missing'] is None else {'kind':out['missing']['kind'],'source':out['missing']['source']},'scope':out['scope']},indent=2))
