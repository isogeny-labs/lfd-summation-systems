"""Provisional divisor/carry-pair closure. No large-cell exactness claim.

The finite-space regression harness is the oracle while the bounded
critical-pair criterion remains a proof obligation.
"""
from dataclasses import dataclass
import heapq,itertools,math,time

@dataclass
class Ring:
    p:int
    m:int
    k:int
    a:int
    D:int
    @property
    def r(self):return self.m*self.k
    def split(self,x):return x%2**self.r,tuple((x>>(self.r+6*j))%64 for j in range(self.a))
    def pack(self,mask,aux):return mask+sum(e<<(self.r+6*j) for j,e in enumerate(aux))
    def degree(self,x):
        mask,es=self.split(x);return mask.bit_count()+sum(es)
    def key(self,x):return self.degree(x),x
    def multiply(self,x,y):
        mask1,e1=self.split(x);mask2,e2=self.split(y);N=1<<self.k
        mask=sum((((mask1>>(i*self.k))%N+(mask2>>(i*self.k))%N)%N)<<(i*self.k) for i in range(self.m))
        es=tuple(a+b for a,b in zip(e1,e2));assert all(e<64 for e in es)
        return self.pack(mask,es)
    def quotient(self,x,y):
        """Return carry-free multiplier when y divides the normal monomial x."""
        mx,ex=self.split(x);my,ey=self.split(y)
        if my&mx!=my or any(a<b for a,b in zip(ex,ey)):return None
        return self.pack(mx^my,tuple(a-b for a,b in zip(ex,ey)))
    def lcm(self,x,y):
        mx,ex=self.split(x);my,ey=self.split(y)
        return self.pack(mx|my,tuple(max(a,b) for a,b in zip(ex,ey)))
    def scale_shift(self,g,mult,c=1):
        out={}
        for x,v in g.items():
            z=self.multiply(x,mult);out[z]=(out.get(z,0)+c*v)%self.p
        return {x:v for x,v in out.items() if v}
    def codes(self):
        def aux(left,n):
            if n==0:yield ();return
            for e in range(left+1):
                for tail in aux(left-e,n-1):yield (e,)+tail
        for w in range(min(self.r,self.D)+1):
            for bits in itertools.combinations(range(self.r),w):
                mask=sum(1<<i for i in bits)
                for es in aux(self.D-w,self.a):yield self.pack(mask,es)
    def size(self):return sum(math.comb(self.r,w)*math.comb(self.D-w+self.a,self.a) for w in range(min(self.r,self.D)+1))

class Compact:
    def __init__(self,ring,*,seconds=30,max_rows=10000,carry=True,pairs=True):
        self.R=ring;self.G=[];self.leads=[];self.trace=[];self.queue=[];self.serial=0
        self.start=time.monotonic();self.seconds=seconds;self.max_rows=max_rows
        self.carry=carry;self.pairs=pairs;self.processed=0;self.capped=False
    def budget(self):
        if time.monotonic()-self.start>self.seconds or len(self.G)>=self.max_rows:
            self.capped=True;raise TimeoutError('compact candidate resource cap')
    def reduce(self,row):
        R=self.R;row={x:c%R.p for x,c in row.items() if c%R.p};rem={};ops=[]
        while row:
            self.budget();x=max(row,key=R.key);c=row[x];found=False
            for j,y in enumerate(self.leads):
                mult=R.quotient(x,y)
                if mult is None:continue
                assert R.degree(y)+R.degree(mult)<=R.D
                product=R.scale_shift(self.G[j],mult)
                assert max(product,key=R.key)==x and product[x]==1
                ops.append([j,mult,c])
                for z,v in product.items():
                    t=(row.get(z,0)-c*v)%R.p
                    if t:row[z]=t
                    else:row.pop(z,None)
                found=True;break
            if not found:rem[x]=row.pop(x)
        return rem,ops
    def enqueue(self,degree,kind,i,j):
        heapq.heappush(self.queue,(degree,self.serial,kind,i,j));self.serial+=1
    def add(self,row,source):
        R=self.R;row,ops=self.reduce(row)
        if not row:return
        lead=max(row,key=R.key);inv=pow(row[lead],-1,R.p)
        row={x:c*inv%R.p for x,c in row.items()};i=len(self.G)
        self.G.append(row);self.leads.append(lead)
        self.trace.append({'source':source,'inverse':inv,'reductions':ops,'row':sorted(row.items())})
        if self.carry and R.degree(lead)<R.D:
            mask,_=R.split(lead)
            for bit in range(R.r):
                if mask>>bit&1:self.enqueue(R.degree(lead)+1,'carry',i,bit)
        if self.pairs:
            for j in range(i):
                lcm=R.lcm(lead,self.leads[j]);degree=R.degree(lcm)
                if degree<=R.D:self.enqueue(degree,'pair',i,j)
    def run(self,gens):
        R=self.R
        try:
            for j,g in enumerate(gens):
                if max(map(R.degree,g),default=0)<=R.D:self.add(g,['input',j])
                if 0 in self.leads:break
            while self.queue and 0 not in self.leads:
                self.budget();_,_,kind,i,j=heapq.heappop(self.queue)
                if kind=='carry':raw=R.scale_shift(self.G[i],1<<j)
                else:
                    lcm=R.lcm(self.leads[i],self.leads[j]);left=R.scale_shift(self.G[i],R.quotient(lcm,self.leads[i]));right=R.scale_shift(self.G[j],R.quotient(lcm,self.leads[j]))
                    raw=dict(left)
                    for x,c in right.items():raw[x]=(raw.get(x,0)-c)%R.p
                    raw={x:c for x,c in raw.items() if c}
                self.add(raw,[kind,i,j]);self.processed+=1
        except TimeoutError:pass
        return {'status':'capped' if self.capped else 'provisional-complete','contains_one':0 in self.leads,'compact_rows':len(self.G),'terms':sum(map(len,self.G)),'pairs_processed':self.processed,'seconds':time.monotonic()-self.start,'claim':'candidate only; compact closedness criterion unproved'}
    def candidate_rank(self):
        if self.capped:raise ValueError('no rank on capped candidate')
        if self.R.size()>100000:raise ValueError('enumeration cap')
        return sum(any(self.R.quotient(x,y) is not None for y in self.leads) for x in self.R.codes())
