#!/usr/bin/env python3
"""Count carry-free leading-monomial coverage; not a closure certificate."""
import argparse,gzip,hashlib,itertools,json
from pathlib import Path
HERE=Path(__file__).resolve().parent

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for b in iter(lambda:stream.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def count_coverage(leads,r,D):
    profile=[]
    for w in range(min(r,D)+1):
        count=0
        for bits in itertools.combinations(range(r),w):
            mon=sum(1<<i for i in bits);sub=mon
            while True:
                if sub in leads:count+=1;break
                if sub==0:break
                sub=(sub-1)&mon
        profile.append(count)
    return profile

def main():
    ap=argparse.ArgumentParser();ap.add_argument('run_dir',type=Path);a=ap.parse_args();d=a.run_dir
    rec=json.loads((d/'result.json').read_text());params=list(map(int,(d/'input.txt').read_text().split()[:5]));p,m,k,aux,D=params
    if aux:raise ValueError('chain-only coverage tool')
    basis=d/'basis.txt.gz';expected=rec['artifacts']['basis.txt.gz']['sha256']
    if digest(basis)!=expected:raise ValueError('basis digest mismatch')
    leads=set();raw_profile=[0]*(D+1);rows=0
    with gzip.open(basis,'rt') as f:
        for line in f:
            fields=line.split(maxsplit=4);ident,nt,lead,c=map(int,fields[:4])
            assert ident==rows and nt>0 and c==1 and lead.bit_count()<=D and lead<2**(m*k)
            assert lead not in leads;leads.add(lead);raw_profile[lead.bit_count()]+=1;rows+=1
    assert rows==rec['derived_rank']
    profile=count_coverage(leads,m*k,D)
    out={'scope':'carry-free leading coverage; row membership must be independently verified; no closedness claim','parameters':params,'source_sha256':digest(Path(__file__)),'basis_sha256':expected,'retained_rows':rows,'retained_leading_counts':raw_profile,'covered_leading_counts':profile,'covered_count':sum(profile),'multiplication_rule':'For each covered monomial M choose a row with leading mask L subset M; multiply by M xor L. Its total degree is at most D.'}
    (HERE/'results'/(d.name+'-leading-coverage.json')).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
