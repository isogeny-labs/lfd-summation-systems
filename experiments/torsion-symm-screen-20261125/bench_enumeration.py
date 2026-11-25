"""Same-target direct point-decomposition enumeration control."""
from __future__ import annotations
import json
import platform
from pathlib import Path
from statistics import median
from time import perf_counter,process_time
from screen import PrimeCurve, prime_solutions, prime_decompositions

HERE=Path(__file__).resolve().parent

def main():
    frozen=json.loads((HERE/'results.json').read_text())
    records=[]
    for cell in frozen['cells']:
        if cell['k']!=4: continue
        E=PrimeCurve(cell['p'],1,0);V=cell['x_values']
        for target in cell['rows']:
            for rep in range(50):
                start_wall=perf_counter();start_cpu=process_time()
                sols=prime_solutions(E,V,target['target'][0])
                count=prime_decompositions(E,V,tuple(target['target']),sols)
                cpu=process_time()-start_cpu;wall=perf_counter()-start_wall
                assert len(sols)==target['scalar_ordered_pairs']
                assert count==target['unordered_x_pairs_with_point_decomposition']
                records.append({'p':cell['p'],'target':target['target'],
                                'repetition':rep,'wall_seconds':wall,
                                'cpu_seconds':cpu,'decompositions':count})
    out={'status':'passed','method':'exact x-pair enumeration plus group-law validation',
         'host':platform.platform(),'protocol':'same frozen target and factor-base membership as Groebner arms; 50 repetitions per target; includes all no-relation attempts',
         'records':records}
    (HERE/'enumeration-timing.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
    for p in (193,257):
        rows=[r for r in records if r['p']==p]
        print(p,'median wall ms',round(1000*median(r['wall_seconds'] for r in rows),3),
              'median cpu ms',round(1000*median(r['cpu_seconds'] for r in rows),3))

if __name__=='__main__':main()
