#!/usr/bin/env python3
"""Bind compressed production receipts before independently replaying rows."""
import argparse,json,subprocess
from pathlib import Path
from run_heap import HERE,digest

def verify(run_dir,*,closed=False,reuse=False,engine='baseline'):
    d=Path(run_dir);receipt=json.loads((d/'result.json').read_text())
    if digest(d/'input.txt')!=receipt['input_sha256']:raise ValueError('input digest mismatch')
    for name in ['derivations.txt.gz','basis.txt.gz']:
        record=receipt['artifacts'][name]
        if digest(d/name)!=record['sha256']:raise ValueError('compressed digest mismatch')
        if digest(d/name,True)!=record['decoded_sha256']:raise ValueError('decoded digest mismatch')
    if engine not in ['baseline','indexed']:raise ValueError('unknown replay engine')
    suffix='' if engine=='baseline' else '-indexed'
    source=HERE/('replay.rs' if engine=='baseline' else 'replay_indexed.rs');source_sha=digest(source);binary=d/('independent-replay'+suffix)
    cached=d/('replay'+suffix+'.json')
    if reuse and cached.exists():
        prior=json.loads(cached.read_text())
        if (prior.get('status')=='passed' and prior.get('verifier_source_sha256')==source_sha
            and prior.get('input_sha256')==receipt['input_sha256'] and prior.get('target')==receipt['artifacts']
            and prior.get('derived_rank')==receipt['derived_rank'] and prior.get('contains_one')==receipt['contains_one']
            and (not closed or prior.get('closed_check'))):
            return dict(prior,cache_reused=True)
    subprocess.run(['rustc','-O',str(source),'-o',str(binary)],check=True)
    cmd=[str(binary),str(d/'input.txt'),str(d/'derivations.txt.gz'),str(d/'basis.txt.gz')]
    if closed:cmd.append('closed')
    with (d/('replay'+suffix+'.stderr')).open('w') as err:
        result=json.loads(subprocess.run(cmd,stdout=subprocess.PIPE,stderr=err,text=True,check=True).stdout)
    assert result['derived_rank']==receipt['derived_rank']
    assert result['contains_one']==receipt['contains_one']
    if closed:assert receipt['status']=='closed'
    result.update(verifier_engine=engine,verifier_source_sha256=source_sha,input_sha256=receipt['input_sha256'],target=receipt['artifacts'])
    cached.write_text(json.dumps(result,indent=2)+'\n')
    (HERE/'results'/(d.name+'-replay'+suffix+'.json')).write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run_dir',type=Path);p.add_argument('--closed',action='store_true');p.add_argument('--engine',choices=['baseline','indexed'],default='baseline');a=p.parse_args()
    print(json.dumps(verify(a.run_dir,closed=a.closed,engine=a.engine),indent=2))
