"""Same-target SymPy grevlex solve timing, including polynomial construction."""
from __future__ import annotations
from collections import defaultdict
import json
import platform
from pathlib import Path
from statistics import median
from time import perf_counter, process_time
import sympy as sp
from screen import direct_generators, symm_generators, torsion_quotient_generators

HERE=Path(__file__).resolve().parent
ARMS=(('direct',direct_generators),('symmetric',symm_generators),('torsion_quotient',torsion_quotient_generators))


def main():
    frozen=json.loads((HERE/'results.json').read_text())
    x,y=sp.symbols('x y')
    records=[]
    for cell in frozen['cells']:
        if cell['k']!=4: continue
        p,N=cell['p'],1<<cell['k']
        for target in cell['rows']:
            t=target['target'][0]
            for rep in range(5):
                for offset in range(3):
                    name,fun=ARMS[(offset+rep)%3]
                    start_wall=perf_counter();start_cpu=process_time()
                    polys=[sp.Poly.from_dict(g,(x,y),modulus=p).as_expr()
                           for g in fun(p,N,t)]
                    basis=sp.groebner(polys,x,y,order='grevlex',modulus=p)
                    cpu=process_time()-start_cpu;wall=perf_counter()-start_wall
                    is_unit=len(basis.polys)==1 and basis.polys[0].total_degree()==0
                    assert is_unit == (target['scalar_ordered_pairs']==0),(p,t,name)
                    records.append({'p':p,'k':cell['k'],'target':target['target'],
                                    'scalar_ordered_pairs':target['scalar_ordered_pairs'],
                                    'arm':name,'repetition':rep,'wall_seconds':wall,
                                    'cpu_seconds':cpu,'unit_ideal':is_unit,
                                    'basis_polynomials':len(basis.polys)})
    out={'status':'passed','solver':'SymPy groebner, grevlex, GF(p)',
         'host':platform.platform(),'python':platform.python_version(),
         'sympy':sp.__version__,'protocol':'same 12 targets per cell, 5 repetitions per arm, rotated arm order; each timing includes input construction and Groebner computation; root extraction, point validation, relation-rank and linear algebra excluded',
         'records':records}
    (HERE/'groebner-timing.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
    for p in (193,257):
        for name,_ in ARMS:
            rows=[r for r in records if r['p']==p and r['arm']==name and r['unit_ideal']]
            print(p,name,len(rows),'median wall ms',round(1000*median(r['wall_seconds'] for r in rows),3),'median cpu ms',round(1000*median(r['cpu_seconds'] for r in rows),3))

if __name__=='__main__':main()
