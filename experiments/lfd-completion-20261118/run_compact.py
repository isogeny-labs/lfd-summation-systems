#!/usr/bin/env python3
"""Bounded compact search; replays every retained derivation independently."""
import argparse,hashlib,json,os,shutil,signal,subprocess,time
from pathlib import Path
from replay import read_input
from compact_candidate import Ring,Compact
from compact_replay import verify
HERE=Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser();ap.add_argument('input',type=Path);ap.add_argument('tag');ap.add_argument('--seconds',type=int,default=30);ap.add_argument('--max-rows',type=int,default=500);ap.add_argument('--engine',choices=['python','rust'],default='python');ap.add_argument('--max-terms',type=int,default=1000000);ap.add_argument('--max-queue',type=int,default=2000000);ap.add_argument('--replay-engine',choices=['python','rust'],default='python');ap.add_argument('--prune-coprime',action='store_true');a=ap.parse_args()
    if a.prune_coprime and a.engine!='rust':raise ValueError('pair pruning requires the Rust producer')
    if min(a.seconds,a.max_rows,a.max_terms,a.max_queue)<1:raise ValueError('positive resource budgets required')
    if shutil.disk_usage(HERE).free<3*2**30:raise RuntimeError('need 3 GiB free')
    dest=HERE/'runs'/a.tag
    if dest.exists():raise FileExistsError(dest)
    params,gens=read_input(a.input);original_sha=hashlib.sha256(a.input.read_bytes()).hexdigest();dest.mkdir();inp=dest/'input.txt';inp.write_bytes(a.input.read_bytes())
    producer=HERE/('compact_producer.rs' if a.engine=='rust' else 'compact_candidate.py')
    if a.replay_engine=='rust':
        from compact_replay_rust import verify as replay_rows
        checkers=[HERE/'compact_replay_rust.py',HERE/'compact_replay.rs']
    else:replay_rows=verify;checkers=[HERE/'compact_replay.py']
    sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [producer,HERE/'replay.py',Path(__file__),*checkers]}
    proof=dest/'compact-proof.json';reason=None;checked=None;budget={'seconds':a.seconds,'max_rows':a.max_rows}
    if a.engine=='python':
        C=Compact(Ring(*params),seconds=a.seconds,max_rows=a.max_rows);result=C.run(gens)
        data={'trace':C.trace}
    else:
        tokens=inp.read_text().split();tokens[5:8]=map(str,[a.seconds,a.max_rows,a.max_terms]);inp.write_text(' '.join(tokens)+'\n')
        budget.update(max_terms=a.max_terms,max_queue=a.max_queue,proof_byte_cap=256*2**20,rss_cap_bytes=2*2**30)
        binary=dest/'producer';producer_bytes=producer.read_bytes();snapshot=dest/'compact_producer.rs';snapshot.write_bytes(producer_bytes);sources[producer.name]=hashlib.sha256(producer_bytes).hexdigest();subprocess.run(['rustc','-O',str(snapshot),'-o',str(binary)],check=True)
        with (dest/'stderr.log').open('w') as err:
            start=time.monotonic();command=[str(binary),str(inp),str(proof),str(a.max_queue)]
            if a.prune_coprime:command.append('coprime')
            child=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=err,text=True,start_new_session=True)
            while True:
                try:
                    stdout,_=child.communicate(timeout=1)
                    if child.returncode:
                        reason='producer-error';result={'status':reason,'returncode':child.returncode};break
                    result=json.loads(stdout);break
                except subprocess.TimeoutExpired:
                    if time.monotonic()-start>a.seconds+30:reason='external-time-cap'
                    if proof.exists() and proof.stat().st_size>256*2**20:reason='proof-byte-cap'
                    if shutil.disk_usage(dest).free<2*2**30:reason='free-disk-floor'
                    rss=subprocess.run(['ps','-o','rss=','-p',str(child.pid)],capture_output=True,text=True).stdout.strip()
                    if rss and int(rss)>2*2**20:reason='rss-cap-2GiB'
                    if reason:
                        os.killpg(child.pid,signal.SIGTERM)
                        try:child.communicate(timeout=5)
                        except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.communicate()
                        result={'status':reason,'claim':'partial file only; no mathematical claim'};break
        if not reason:data=json.loads(proof.read_text())
    if not reason:
        data['input_sha256']=hashlib.sha256(inp.read_bytes()).hexdigest();proof.write_text(json.dumps(data,separators=(',',':'))+'\n')
        checked=replay_rows(inp,proof);assert checked['contains_one']==result['contains_one']
        assert checked['verified_rows']==result['compact_rows']
    if checked:
        result['producer_claim']=result.get('claim');result['claim']='retained polynomial membership independently replayed; compact closedness remains unproved'
        if checked['contains_one']:result['claim']='constant polynomial 1 independently derived within the recorded degree; lower bound remains separate'
    result.update(pair_strategy='coprime-pruned' if a.prune_coprime else 'all-bounded',replay_engine=a.replay_engine,engine=a.engine,budget=budget,replay=checked,source_sha256=sources,input_sha256=hashlib.sha256(inp.read_bytes()).hexdigest(),original_input_sha256=original_sha,proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest() if proof.exists() else None)
    (dest/'result.json').write_text(json.dumps(result,indent=2)+'\n');(HERE/'results'/(a.tag+'.json')).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
    if reason:raise SystemExit(1)
if __name__=='__main__':main()
