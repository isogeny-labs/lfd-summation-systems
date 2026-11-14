#!/usr/bin/env python3
"""Independent exact characteristic-zero polynomial Wronskian regressions."""
import random
import sympy as s
x=s.symbols("x")
rng=random.Random(11013)
count=0
for r in range(1,6):
    for sample in range(8):
        degrees=sorted(rng.sample(range(0,12),r))
        ps=[]
        for d in degrees:
            leading=rng.choice([-3,-2,-1,1,2,3])
            ps.append(leading*x**d+sum(rng.randint(-3,3)*x**i for i in range(d)))
        wr=s.Matrix([[s.diff(p,x,i) for p in ps] for i in range(r)])
        det=s.Poly(s.expand(wr.det(method="berkowitz")),x)
        predicted=sum(degrees)-r*(r-1)//2
        vandermonde=1
        for i in range(r):
            for j in range(i+1,r):
                vandermonde*=degrees[j]-degrees[i]
        leadprod=1
        for p in ps: leadprod*=s.LC(s.Poly(p,x))
        assert det.degree()==predicted,(r,degrees,det.degree(),predicted)
        assert det.LC()==leadprod*vandermonde,(r,degrees,det.LC(),leadprod*vandermonde)
        count+=1
print("N8 independent distinct-degree Wronskian regressions passed:",count)
