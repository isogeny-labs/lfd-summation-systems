#!/usr/bin/env python3
"""Attempt to falsify the compact candidate against the direct closure."""
import argparse,hashlib,itertools,json,random,subprocess,tempfile
from pathlib import Path
from compact_candidate import Ring,Compact
from compact_replay import verify
HERE=Path(__file__).resolve().parent

def text_input(R,gens):
    return f'{R.p} {R.m} {R.k} {R.a} {R.D} 10 100000 1000000\n{len(gens)}\n'+''.join(str(len(g))+' '+' '.join(f'{x} {c}' for x,c in sorted(g.items()))+'\n' for g in gens)

def direct(binary,dest,R,gens):
    inp=dest/'input';inp.write_text(text_input(R,gens))
    out=json.loads(subprocess.run([str(binary),str(inp),str(dest/'ledger'),str(dest/'basis')],capture_output=True,text=True,check=True).stdout)
    if out['status'] not in ['closed','unit-witness']:raise ValueError('reference cap')
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--trials',type=int,default=200);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--exhaustive-univariate',action='store_true');ap.add_argument('--compiled',action='store_true');args=ap.parse_args()
    rng=random.Random(20261119);records=[]
    sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'compact_candidate.py',Path(__file__),HERE/'sparse_producer.rs',HERE/'compact_replay.py',HERE/'replay.py',HERE/'compact_producer.rs']}
    with tempfile.TemporaryDirectory(prefix='compact-falsify-') as tmp:
        d=Path(tmp);binary=d/'reference'
        subprocess.run(['rustc','-O',str(HERE/'sparse_producer.rs'),'-o',str(binary)],check=True)
        compiled=d/'compact'
        if args.compiled:subprocess.run(['rustc','-O',str(HERE/'compact_producer.rs'),'-o',str(compiled)],check=True)
        # The harness must catch these deliberately incomplete criteria.
        R=Ring(17,1,2,0,2);g=[{2:1}];ref=direct(binary,d,R,g)
        bad=Compact(R,carry=False).run(g);assert ref['contains_one'] and not bad['contains_one']
        R=Ring(17,1,1,3,3);x,y,z=[1<<(1+6*j) for j in range(3)];g=[{x+y:1,0:1},{x+z:1,0:1}]
        ref=direct(binary,d,R,g);bad=Compact(R,pairs=False);bad.run(g)
        assert not ref['contains_one'] and bad.candidate_rank()!=ref['derived_rank']
        cases=[]
        for i in range(args.trials):
            R=Ring(17,1+i%2,1+(i//2)%2,(i//4)%2,1+(i//8)%4)
            codes=list(R.codes());gens=[]
            for _ in range(1+i%3):
                chosen=rng.sample(codes,min(len(codes),2+i%5));g={x:rng.randrange(1,R.p) for x in chosen}
                if i%2==0:
                    g[0]=(-sum(c for x,c in g.items() if x!=0))%R.p
                g={x:c for x,c in g.items() if c};gens.append(g)
            cases.append((R,gens))
        if args.exhaustive_univariate:
            for D in [1,2,3]:
                for coeffs in itertools.product(range(5),repeat=4):
                    cases.append((Ring(5,1,2,0,D),[{x:c for x,c in enumerate(coeffs) if c}]))
        for i,(R,gens) in enumerate(cases):
            ref=direct(binary,d,R,gens);candidate=Compact(R,seconds=5);out=candidate.run(gens)
            proof=d/'proof.json';trace=candidate.trace;python_out=out;trace_identity=None
            if args.compiled:
                out=json.loads(subprocess.run([str(compiled),str(d/'input'),str(proof)],capture_output=True,text=True,check=True).stdout)
                trace=json.loads(proof.read_text())['trace'];trace_identity=trace==json.loads(json.dumps(candidate.trace))
                if not trace_identity:raise AssertionError(f'compiled derivation differs in case {i}')
            proof.write_text(json.dumps({'input_sha256':hashlib.sha256((d/'input').read_bytes()).hexdigest(),'trace':trace}))
            membership=verify(d/'input',proof)
            assert membership['contains_one']==out['contains_one']
            rank=None if out['contains_one'] else candidate.candidate_rank() if not candidate.capped else None
            ok=not candidate.capped and out['contains_one']==ref['contains_one'] and (out['contains_one'] or rank==ref['derived_rank'])
            rec={'case':i,'parameters':vars(R),'input_sha256':hashlib.sha256(text_input(R,gens).encode()).hexdigest(),'candidate':out,'python_candidate':python_out,'compiled_trace_identity':trace_identity,'candidate_rank':rank,'reference':ref,'membership_replay':membership,'agrees':ok}
            records.append(rec)
            if not ok:
                rec['failure_kind']='resource-cap' if candidate.capped else 'reading-mismatch'
                rec['reproduction_input']=text_input(R,gens);break
    result={'scope':'small-instance falsification only; no large-cell exactness claim','source_sha256':sources,'negative_controls':{'missing_carry_detected':True,'missing_pair_detected':True},'cases':records,'all_agree':len(records)==len(cases) and all(x['agrees'] for x in records)}
    args.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'all_agree':result['all_agree'],'cases':len(records),'last':records[-1] if not result['all_agree'] else None},indent=2))
    if not result['all_agree']:raise SystemExit(1)
if __name__=='__main__':main()
