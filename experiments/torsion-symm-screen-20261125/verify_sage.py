"""Independent Sage dense-space check of selected unit-ideal closure thresholds."""
from sage.all import GF, VectorSpace, vector
from screen import direct_generators, symm_generators, torsion_quotient_generators


def contains_one(generators, p, D):
    monoms=[(i,d-i) for d in range(D,-1,-1) for i in range(d,-1,-1)]
    index={m:j for j,m in enumerate(monoms)}
    space=VectorSpace(GF(p),len(monoms))
    def row(poly):
        v=vector(GF(p),len(monoms))
        for m,c in poly.items():
            if sum(m)<=D: v[index[m]]=c
        return v
    basis=space.subspace([row(g) for g in generators if g and max(map(sum,g))<=D])
    while True:
        new=[]
        for v in basis.basis():
            if any(v[col] and i+j==D for (i,j),col in index.items()):
                continue
            for axis in (0,1):
                w=vector(GF(p),len(monoms))
                for (i,j),col in index.items():
                    if v[col] and i+j<D:
                        next_m=(i+1,j) if axis==0 else (i,j+1)
                        w[index[next_m]]=v[col]
                new.append(w)
        nxt=space.subspace(list(basis.basis())+new)
        if nxt.dimension()==basis.dimension(): break
        basis=nxt
    one=vector(GF(p),len(monoms));one[index[(0,0)]]=1
    return one in basis,basis.dimension()


def main():
    cases=((97,8,5,[(direct_generators,9),(symm_generators,8),(torsion_quotient_generators,7)]),
           (193,16,5,[(direct_generators,17),(symm_generators,16),(torsion_quotient_generators,15)]))
    for p,N,t,arms in cases:
        for fun,D in arms:
            below,rank_below=contains_one(fun(p,N,t),p,D-1)
            at,rank_at=contains_one(fun(p,N,t),p,D)
            print(p,N,t,fun.__name__,D,below,at,rank_below,rank_at)
            assert not below and at

if __name__=='__main__': main()
