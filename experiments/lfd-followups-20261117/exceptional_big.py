#!/usr/bin/env python3
"""Bounded optimized-engine cross-check of the deeper proper-ideal case."""
import hashlib,json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from csd import ic_bounds_big as B
from csd import ic_bounds_d5_verify as R
HERE=Path(__file__).resolve().parent

def main():
    k,p,D=10,12289,6
    Z=R.solutions(p,6,p-7,1,k);g=B.chain_generator(p,6,p-7,1,k)
    assert g==R.chain_generator(p,6,p-7,1,k) and len(Z)==86
    inst=dict(k=k,prime=p,a=6,b=p-7,seed=None,target=[1,0],points=len(Z),Z=Z,
              generator={f'{i},{j}':c for (i,j),c in g.items()})
    src=HERE/'exceptional-k10-instance.json';src.write_text(json.dumps(inst,indent=2)+'\n')
    out=HERE/'exceptional-big-k10-p12289.json';log=HERE/'exceptional-big-k10.log'
    cmd=[sys.executable,str(ROOT/'scripts/ic_bounds_big_run.py'),'closure','--instance',str(src),
         '--D','6','--device','cpu','--work-gb','1','--tile-mb','32','--batch','128','--threads','4',
         '--no-eval','--out',str(out)]
    start=time.time()
    with log.open('w') as stream:
        try: result=subprocess.run(cmd,stdout=stream,stderr=subprocess.STDOUT,timeout=900);status='completed';code=result.returncode
        except subprocess.TimeoutExpired:status='capped';code=None
    receipt=dict(status=status,exit_code=code,seconds=round(time.time()-start,3),cap_seconds=900,
                 source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    if status=='completed':
        assert code==0
        j=json.loads(out.read_text());assert j['closed'] and j['closure_complete'] and not j['contains_one']
        ev=R.eval_rank(B.Layout(k,D,False),Z,p,k,D-1)
        # Lower-degree space comes from a completed independent reference run.
        prior=json.loads((HERE/'exceptional-reference-k10-p12289-240s.json').read_text())
        prev=next(s for s in prior['profile'] if s['D']==5)
        assert prev['complete'] and prev['closed_check']['passed']
        strict=j['profile'][D-1]>prev['rank']
        exact=bool(j['columns']-j['rank']==len(Z)==ev and strict)
        receipt.update(points=len(Z),evaluation_rank_below_D=ev,codimension=j['columns']-j['rank'],
                       strict_fall_at_D6=strict,exact_lfd=6 if exact else None,
                       scope='closed-space/local-to-global check; full reference D6 may remain capped')
        assert exact
    (HERE/'exceptional-big-k10-reading.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
