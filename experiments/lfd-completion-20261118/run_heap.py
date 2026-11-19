#!/usr/bin/env python3
"""Guarded compressed proof production on an existing, frozen input."""
import argparse,gzip,hashlib,json,shutil,subprocess,time
from pathlib import Path
HERE=Path(__file__).resolve().parent

def digest(path,decoded=False):
    h=hashlib.sha256()
    with (gzip.open(path,'rb') if decoded else path.open('rb')) as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('input',type=Path);ap.add_argument('tag');ap.add_argument('--seconds',type=int,default=1800);ap.add_argument('--engine',choices=['heap','indexed','restart'],default='heap');ap.add_argument('--resume',type=Path);args=ap.parse_args()
    resume_record=None
    if args.resume:
        if args.engine!='restart':raise ValueError('--resume requires --engine restart')
        from verify_compressed import verify
        verify(args.resume)
        old=(args.resume/'input.txt').read_text().split();new=args.input.read_text().split()
        if old[:5]!=new[:5] or old[8:]!=new[8:]:raise ValueError('resume changes the polynomial system or degree')
        resume_record={'run_dir':str(args.resume.resolve()),'result_sha256':digest(args.resume/'result.json')}
    dest=HERE/'runs'/args.tag
    if dest.exists():raise FileExistsError(dest)
    if shutil.disk_usage(HERE).free<3*2**30:raise RuntimeError('need at least 3 GiB free before launch')
    dest.mkdir();tokens=args.input.read_text().split();tokens[5]=str(args.seconds)
    (dest/'input.txt').write_text(' '.join(tokens)+'\n')
    source=HERE/(args.engine+'_producer.rs')
    source_hashes={p.name:digest(p) for p in [source,Path(__file__)]}
    binary=dest/'producer';subprocess.run(['rustc','-O',str(source),'-o',str(binary)],check=True)
    outputs=[dest/'derivations.txt.gz',dest/'basis.txt.gz'];start=time.monotonic();reason=None
    with (dest/'stderr.log').open('w') as err:
        command=[str(binary),str(dest/'input.txt'),*map(str,outputs)]
        if args.resume:command.extend(str(args.resume/name) for name in ['derivations.txt.gz','basis.txt.gz'])
        child=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=err,text=True,start_new_session=True)
        while True:
            try:
                stdout,_=child.communicate(timeout=2)
                if child.returncode:raise RuntimeError(f'producer exit {child.returncode}')
                result=json.loads(stdout);break
            except subprocess.TimeoutExpired:
                if shutil.disk_usage(dest).free<2*2**30:reason='free-disk-floor'
                if sum(p.stat().st_size for p in outputs if p.exists())>2*2**30:reason='artifact-byte-cap'
                if time.monotonic()-start>args.seconds+120:reason='external-timeout'
                if reason:
                    # Terminate the whole producer/compressor process group.
                    import os,signal
                    os.killpg(child.pid,signal.SIGTERM)
                    try:child.communicate(timeout=5)
                    except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.communicate()
                    result={'status':reason,'scope':'partial compressed files; no complete closure claim'};break
    result.update(resume_from=resume_record,engine=args.engine,source_sha256=source_hashes,input_sha256=digest(dest/'input.txt'),original_input_sha256=digest(args.input),elapsed_seconds=time.monotonic()-start)
    result['artifacts']={p.name:{'sha256':digest(p),'bytes':p.stat().st_size} for p in outputs}
    if not reason:
        for p in outputs:result['artifacts'][p.name]['decoded_sha256']=digest(p,True)
    (dest/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    (HERE/'results'/(args.tag+'.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':main()
