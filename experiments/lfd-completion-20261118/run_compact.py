#!/usr/bin/env python3
"""Bounded compact search; replays every retained derivation independently."""
import argparse,hashlib,json,shutil
from pathlib import Path
from replay import read_input
from compact_candidate import Ring,Compact
from compact_replay import verify
HERE=Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser();ap.add_argument('input',type=Path);ap.add_argument('tag');ap.add_argument('--seconds',type=int,default=30);ap.add_argument('--max-rows',type=int,default=500);a=ap.parse_args()
    if shutil.disk_usage(HERE).free<3*2**30:raise RuntimeError('need 3 GiB free')
    dest=HERE/'runs'/a.tag
    if dest.exists():raise FileExistsError(dest)
    dest.mkdir();inp=dest/'input.txt';inp.write_bytes(a.input.read_bytes())
    sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'compact_candidate.py',HERE/'compact_replay.py',HERE/'replay.py',Path(__file__)]}
    params,gens=read_input(inp);C=Compact(Ring(*params),seconds=a.seconds,max_rows=a.max_rows);result=C.run(gens)
    proof=dest/'compact-proof.json';proof.write_text(json.dumps({'input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest(),'trace':C.trace},separators=(',',':'))+'\n')
    checked=verify(inp,proof);assert checked['contains_one']==result['contains_one']
    result.update(budget={'seconds':a.seconds,'max_rows':a.max_rows},replay=checked,source_sha256=sources,input_sha256=hashlib.sha256(inp.read_bytes()).hexdigest(),proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest())
    (dest/'result.json').write_text(json.dumps(result,indent=2)+'\n');(HERE/'results'/(a.tag+'.json')).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':main()
