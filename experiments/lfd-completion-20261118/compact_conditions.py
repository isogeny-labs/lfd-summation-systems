#!/usr/bin/env python3
"""Falsify finite compact-pair conditions; this is not an independent proof."""
import argparse,hashlib,json
from pathlib import Path
from compact_candidate import Ring,Compact
from replay import read_input
HERE=Path(__file__).resolve().parent

def check(run_dir,seconds=30):
    d=Path(run_dir);inp=d/'input.txt';raw=d/'candidate.json';candidate=json.loads(raw.read_text())
    if hashlib.sha256(inp.read_bytes()).hexdigest()!=candidate['input_sha256']:raise ValueError('input digest')
    params,gens=read_input(inp);R=Ring(*params);polys=[]
    for pairs in candidate['polynomials']:
        row=dict(pairs)
        if len(row)!=len(pairs) or any(not 0<c<R.p or R.degree(x)>R.D for x,c in row.items()):raise ValueError('candidate polynomial')
        lead=max(row,key=R.key);inv=pow(row[lead],-1,R.p);polys.append({x:c*inv%R.p for x,c in row.items()})
    C=Compact(R,seconds=seconds,max_rows=len(polys)+1);C.G=polys;C.leads=[max(g,key=R.key) for g in polys];checked={'inputs':0,'carries':0,'pairs':0};missing=None
    def test(poly,kind,*ids):
        nonlocal missing
        rem,_=C.reduce(poly);checked[kind]+=1
        if rem:
            missing={'kind':kind,'source':ids,'remainder':sorted(rem.items())};return False
        return True
    try:
        for i,g in enumerate(gens):
            if max(map(R.degree,g),default=0)<=R.D and not test(g,'inputs',i):break
        if missing is None:
            for i,g in enumerate(polys):
                L=C.leads[i]
                if R.degree(L)<R.D:
                    for bit in range(R.r):
                        if L>>bit&1 and not test(R.scale_shift(g,1<<bit),'carries',i,bit):break
                if missing:break
        if missing is None:
            for i,g in enumerate(polys):
                for j in range(i):
                    lcm=R.lcm(C.leads[i],C.leads[j])
                    if R.degree(lcm)>R.D:continue
                    row=R.scale_shift(g,R.quotient(lcm,C.leads[i]));other=R.scale_shift(polys[j],R.quotient(lcm,C.leads[j]))
                    for x,c in other.items():row[x]=(row.get(x,0)-c)%R.p
                    if not test({x:c for x,c in row.items() if c},'pairs',i,j):break
                if missing:break
        status='conditions-failed' if missing else 'finite-conditions-passed'
    except TimeoutError:status='capped'
    out={'status':status,'scope':'prototype finite-condition check only; not an independent bounded-closedness certificate','checked':checked,'missing':missing,'candidate_sha256':hashlib.sha256(raw.read_bytes()).hexdigest(),'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),HERE/'compact_candidate.py']}}
    (HERE/'results'/(d.name+'-conditions.json')).write_text(json.dumps(out,indent=2)+'\n');return out
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('run_dir',type=Path);a=ap.parse_args();print(json.dumps(check(a.run_dir),indent=2))
