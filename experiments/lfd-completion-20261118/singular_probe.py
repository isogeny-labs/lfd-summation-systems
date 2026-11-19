#!/usr/bin/env python3
"""Bounded local Sage/libSingular candidate generation; not a certificate."""
import argparse,hashlib,json,os,signal,subprocess,sys,time,shutil
from pathlib import Path
HERE=Path(__file__).resolve().parent

def worker(inp,dest):
    from sage.all import GF,PolynomialRing
    from sage.version import version
    from replay import read_input
    params,gens=read_input(inp);p,m,k,a,D=params;r=m*k;n=r+a
    P=PolynomialRing(GF(p),names=tuple('z'+str(i) for i in reversed(range(n))),order='deglex');variables=list(reversed(P.gens()))
    def expand(g):
        out=P.zero()
        for code,c in g.items():
            term=P(c)
            for j in range(r):
                if code>>j&1:term*=variables[j]
            for j in range(a):term*=variables[r+j]**((code>>(r+6*j))%64)
            out+=term
        return out
    chains=[variables[b*k+j]**2-(variables[b*k+j+1] if j+1<k else 1) for b in range(m) for j in range(k)]
    start=time.monotonic();G=P.ideal(chains+[expand(g) for g in gens]).groebner_basis(algorithm='libsingular:std',deg_bound=D)
    polys=[]
    for g in G:
        normal={}
        for powers,c in g.dict().items():
            powers=list(reversed(powers));code=0
            for b in range(m):code+=(sum(int(powers[b*k+j])*(1<<j) for j in range(k))%(1<<k))<<(b*k)
            for j in range(a):
                assert powers[r+j]<64;code+=int(powers[r+j])<<(r+6*j)
            normal[code]=(normal.get(code,0)+int(c))%p
        normal={x:c for x,c in normal.items() if c}
        if normal:polys.append(sorted(normal.items()))
    result={'status':'completed','scope':'local CAS candidate only; no V_D membership or compact closedness certificate','parameters':params,'sage_version':version,'algorithm':'libsingular:std','degree_bound':D,'seconds':time.monotonic()-start,'raw_polynomials':len(G),'normal_polynomials':len(polys),'terms':sum(map(len,polys)),'contains_one_candidate':any(len(g)==1 and g[0][0]==0 for g in polys)}
    (dest/'candidate.json').write_text(json.dumps({'input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest(),'polynomials':polys},separators=(',',':'))+'\n')
    (dest/'worker-result.json').write_text(json.dumps(result,indent=2)+'\n')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('input',type=Path);ap.add_argument('tag');ap.add_argument('--worker',action='store_true');ap.add_argument('--seconds',type=int,default=120);a=ap.parse_args()
    if a.worker:worker(a.input,Path(a.tag));return
    if shutil.disk_usage(HERE).free<3*2**30:raise RuntimeError('need 3 GiB free')
    d=HERE/'runs'/a.tag
    if d.exists():raise FileExistsError(d)
    d.mkdir();inp=d/'input.txt';inp.write_bytes(a.input.read_bytes());start=time.monotonic();reason=None
    source_sha=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    with (d/'stdout.log').open('w') as out,(d/'stderr.log').open('w') as err:
        child=subprocess.Popen(['sage','-python',str(Path(__file__).resolve()),str(inp),str(d),'--worker'],stdout=out,stderr=err,start_new_session=True)
        while child.poll() is None:
            if time.monotonic()-start>a.seconds:reason='wall-time-cap'
            rss=subprocess.run(['ps','-o','rss=','-p',str(child.pid)],capture_output=True,text=True).stdout.strip()
            if rss and int(rss)>2*2**20:reason='rss-cap-2GiB'
            if shutil.disk_usage(d).free<2*2**30:reason='free-disk-floor'
            if reason:
                os.killpg(child.pid,signal.SIGTERM)
                try:child.wait(timeout=5)
                except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()
                break
            time.sleep(1)
    if reason:result={'status':'capped','reason':reason,'scope':'no completed candidate or mathematical claim'}
    elif child.returncode:result={'status':'producer-error','returncode':child.returncode,'scope':'no mathematical claim; see stderr.log'}
    else:result=json.loads((d/'worker-result.json').read_text())
    result.update(input_sha256=hashlib.sha256(inp.read_bytes()).hexdigest(),source_sha256=source_sha,elapsed_seconds=time.monotonic()-start)
    if (d/'candidate.json').exists():result['candidate_sha256']=hashlib.sha256((d/'candidate.json').read_bytes()).hexdigest()
    (d/'result.json').write_text(json.dumps(result,indent=2)+'\n');(HERE/'results'/(a.tag+'.json')).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
