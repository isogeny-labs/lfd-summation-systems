#!/usr/bin/env python3
"""Independent compressed closure and local-to-global checks, bounded per case."""
import hashlib,json,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from csd import ic_bounds_d5_verify as R
HERE=Path(__file__).resolve().parent

def run(k,p,cap,seconds=240):
    start=time.time();deadline=start+seconds
    g=R.chain_generator(p,6,p-7,1,k);pts=R.solutions(p,6,p-7,1,k)
    assert pts and (p-1)%(2**(k+1))==0
    rows=[];previous_rank=0;falls=[];value=None
    for D in range(2,cap+1):
        C=R.Closure(g,k,p,D,split=True)
        complete=C.run(batch=48,deadline=deadline)
        sweep=C.sweep_full(deadline=deadline) if complete else None
        checked=bool(complete and sweep['complete'] and sweep['passed'])
        profile,_=C.profile();rank=C.rank()
        entry=dict(D=D,rank=rank,columns=C.L.c,closed_check=sweep,complete=checked,contains_one=C.contains_one(),profile=profile)
        rows.append(entry)
        if not checked: break
        assert not C.contains_one()
        if profile[D-1]>previous_rank: falls.append(D)
        previous_rank=rank
        ev=R.eval_rank(C.L,pts,p,k,D-1)
        entry.update(evaluation_rank_below_D=ev,codimension=C.L.c-rank)
        if ev==len(pts)==C.L.c-rank:
            value=max(falls,default=0);break
        if time.time()>deadline: break
    out=dict(k=k,p=p,curve=[6,-7],target=[1,0],points=len(pts),exact_lfd=value,
             status='exact-proper-local' if value is not None else 'capped',strict_falls=falls,
             lfd_lower_bound=max(falls,default=0),profile=rows,seconds=round(time.time()-start,3),
             source_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [Path(__file__),ROOT/'csd/ic_bounds_d5_verify.py']})
    (HERE/f'exceptional-reference-k{k}-p{p}.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['profile','source_sha256']}),flush=True)

if __name__=='__main__':
    for k,p,cap in [(7,257,5),(7,769,5),(8,7681,5),(9,12289,5),(10,12289,6)]:run(k,p,cap)
