"""Pure-integer DAG checker. Imports no producer or repository algebra code.

Proves output rows belong to V_D modulo the chain equations. Each shift
requires its SOURCE polynomial to have weight < D, before normal-form
reduction. The degree-compatible chain relations justify that reduction.
Does not prove completeness, independence, or anything about unprovided
large-run basis rows. Generators must be supplied independently by caller.
"""
import hashlib
import json


def verify(cert, *, generators, parameters):
    m,k,p,D = parameters
    if tuple(cert[x] for x in ('m','k','p','D')) != parameters:
        raise ValueError('parameter mismatch')
    if m < 1 or k < 1 or D < 0 or p < 3 or any(p%d==0 for d in range(2, int(p**0.5)+1)) or (p-1) % (1<<k):
        raise ValueError('invalid field or chains')
    def parse(entries):
        out = {}
        for e,c in entries:
            if type(e) is not int or type(c) is not int or e<0 or e>=1<<(m*k) or e in out or not 0<c<p:
                raise ValueError('noncanonical polynomial')
            out[e]=c
        return out
    inputs = [{e:int(c)%p for e,c in g.items() if int(c)%p} for g in generators]
    if [parse(g) for g in cert['generators']] != inputs:
        raise ValueError('input mismatch')
    values=[]
    for i,node in enumerate(cert['nodes']):
        def ref(j):
            if type(j) is not int or not 0<=j<i:
                raise ValueError('forward or invalid reference')
            return values[j]
        op=node['op']
        if op=='input':
            j=node['index']
            if type(j) is not int or not 0<=j<len(inputs): raise ValueError('input index')
            row=dict(inputs[j])
        elif op=='lin':
            row={}
            for j,c in node['terms']:
                if type(c) is not int or not 0<c<p: raise ValueError('coefficient')
                for e,v in ref(j).items(): row[e]=(row.get(e,0)+c*v)%p
        elif op=='shift':
            source=ref(node['source']);j=node['variable']
            if type(j) is not int or not 0<=j<m*k: raise ValueError('variable')
            if max((e.bit_count() for e in source), default=0)>=D:
                raise ValueError('inadmissible unreduced degree')
            row={}
            b,bit=divmod(j,k);offset=b*k;N=1<<k
            for e,c in source.items():
                # Reconstruct the exponent blocks independently of producer.
                exps=[(e>>(s*k)) % N for s in range(m)]
                exps[b]=(exps[b]+2**bit)%N
                dest=sum(v<<(s*k) for s,v in enumerate(exps))
                row[dest]=(row.get(dest,0)+c)%p
        else: raise ValueError('unknown operation')
        row={e:c for e,c in row.items() if c}
        if max((e.bit_count() for e in row), default=0)>D: raise ValueError('degree cap')
        values.append(row)
    rows=[parse(row) for row in cert['rows']]
    if len(rows)!=len(cert['outputs']): raise ValueError('output count')
    for j,row in zip(cert['outputs'],rows):
        if type(j) is not int or not 0<=j<len(values) or values[j]!=row:
            raise ValueError('output mismatch')
    digest=hashlib.sha256(json.dumps(cert['rows'],separators=(',',':')).encode()).hexdigest()
    if digest!=cert['row_sha256']: raise ValueError('output digest')
    return rows
