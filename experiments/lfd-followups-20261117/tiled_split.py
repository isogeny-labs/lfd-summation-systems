"""Float32 tiled split closure; exact float64 arithmetic, no full basis copies.

The basis keeps all columns but only 4-byte residues. Updates touch a small
row/column tile, and the product queue references current low basis rows
instead of storing another dense copy. Shifting each low pivot once is
valid because subsequent changes use later low pivots, also enqueued.
"""
import time
from collections import deque
import numpy as np
from csd import elim
from csd.ic_bounds_m3_batched import TWO53


class TiledClosure:
    def __init__(self,space,p,batch=64,tile_rows=64,tile_cols=2048,max_basis_gib=8):
        self.sp,self.p,self.n,self.batch=space,p,len(space),batch
        if p>=2**24 or self.n*(p-1)**2>=TWO53:raise ValueError('exactness bound')
        self.tile_rows,self.tile_cols=tile_rows,tile_cols
        self.max_basis_bytes=max_basis_gib*2**30;self.stop_reason=None
        self.blocks=[];self.piv=[];self.locations=[];self.queue=deque();self.contains_one=False
    def _reduce(self,P):
        P=np.asarray(P,dtype=np.float64).copy()%self.p
        start=0
        for block in self.blocks:
            count=len(block);c=P[:,self.piv[start:start+count]].copy()
            # Snapshot coefficients before any columns are changed.
            for a in range(0,self.n,self.tile_cols):
                P[:,a:a+self.tile_cols]-=c@block[:,a:a+self.tile_cols].astype(np.float64)
                np.mod(P[:,a:a+self.tile_cols],self.p,out=P[:,a:a+self.tile_cols])
            start+=count
        return P
    def _absorb(self,R):
        p=self.p;R=R[np.any(R!=0,axis=1)]
        if not len(R):return
        out=[];rank,piv,_=elim.fp_echelon(R.astype(np.int64),p,out)
        if not rank:return
        E=out[0].astype(np.int64)
        for i in range(rank-1,0,-1):
            c=E[:i,piv[i]].copy()
            E[:i]=(E[:i]-np.outer(c,E[i]))%p
        Ef=E.astype(np.float64)
        for block in self.blocks:
            c=block[:,piv].astype(np.float64)
            for a in range(0,self.n,self.tile_cols):
                piece=block[:,a:a+self.tile_cols].astype(np.float64)-c@Ef[:,a:a+self.tile_cols]
                block[:,a:a+self.tile_cols]=np.mod(piece,p)
        for a in range(0,rank,self.tile_rows):
            block=Ef[a:a+self.tile_rows].astype(np.float32)
            bid=len(self.blocks);self.blocks.append(block)
            for j in range(len(block)):
                pivot=int(piv[a+j]);self.piv.append(pivot);self.locations.append((bid,j))
                if pivot==self.n-1:self.contains_one=True
                if self.sp.deg[pivot]<self.sp.D:
                    idx=len(self.piv)-1
                    self.queue.extend((idx,var) for var in range(len(self.sp.cols)))
    def add_inputs(self,rows):
        if rows:self._absorb(self._reduce(np.asarray(rows)))
    def _product(self,items):
        P=np.zeros((len(items),self.n))
        for i,(idx,var) in enumerate(items):
            block,row=self.locations[idx];v=self.blocks[block][row]
            sup=np.flatnonzero(v);dest=self.sp.cols[var,sup]
            if np.any(dest<0):raise AssertionError('non-low queued row')
            P[i,dest]=v[sup]
        return P
    def saturate(self,deadline=None,log=None):
        batches=0
        while self.queue and not self.contains_one:
            if deadline and time.time()>deadline:
                self.stop_reason="deadline";return False
            if (len(self.piv)+self.batch)*self.n*4>self.max_basis_bytes:
                self.stop_reason="basis-memory-cap";return False
            items=[self.queue.popleft() for _ in range(min(self.batch,len(self.queue)))]
            self._absorb(self._reduce(self._product(items)));batches+=1
            if log and batches%20==0:log(dict(batches=batches,rank=len(self.piv),queue=len(self.queue)))
        return True
    def verify_closed(self,deadline=None):
        if self.contains_one:return True
        for i,pivot in enumerate(self.piv):
            if self.sp.deg[pivot]>=self.sp.D:continue
            items=[(i,j) for j in range(len(self.sp.cols))]
            if np.any(self._reduce(self._product(items))):return False
            if deadline and time.time()>deadline:return None
        return True
    def contains(self,row):return not np.any(self._reduce(np.asarray(row)[None,:]))
    @property
    def rank(self):return self.n if self.contains_one else len(self.piv)
    def low_dimension(self,d):
        return int(np.sum(self.sp.deg<=d)) if self.contains_one else int(np.sum(self.sp.deg[self.piv]<=d))
