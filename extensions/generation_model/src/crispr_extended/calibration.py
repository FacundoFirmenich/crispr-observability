"""MIT. Exact source-fragment scores and split-conformal locus-block sets.

These are probabilistic coverage sets under exchangeable blocks, not literal
source-base certificates. A score above the computed cap is infinity.
"""
import math
import numpy as np
from numba import njit
from .core import _encode

@njit(cache=True)
def _fragment_distance(h, read, cap):
    n=len(h);m=len(read);inf=cap+1
    row=np.full(m+1,inf,np.int32);best=inf
    for i in range(n+1):
        nxt=np.full(m+1,inf,np.int32)
        for j in range(m+1):
            if i<n and j<=cap:
                if j+1<nxt[j]:nxt[j]=j+1
                if j<m:
                    c=j+int(h[i]!=read[j])
                    if c<nxt[j+1]:nxt[j+1]=c
            cost=row[j]
            if cost>cap:continue
            if j==m and cost<best:best=cost
            if j<m and cost<cap and cost+1<row[j+1]:row[j+1]=cost+1
            if i<n:
                if cost<cap and cost+1<nxt[j]:nxt[j]=cost+1
                if j<m:
                    c=cost+int(h[i]!=read[j])
                    if c<nxt[j+1]:nxt[j+1]=c
        row=nxt
    return best

def fragment_distance(sequence, read, cap=32):
    if type(cap) is not int or not 0<=cap<=1000:raise ValueError('Invalid cap')
    c=int(_fragment_distance(_encode(sequence,True),_encode(read),cap))
    return c if c<=cap else None

def allele_scores(read,catalogue,cap=32):
    scores={}
    for h in catalogue:
        allele=str(h['allele']);c=fragment_distance(h['sequence'],read,cap)
        if allele not in scores or (c is not None and (scores[allele] is None or c<scores[allele])):scores[allele]=c
    if not scores:raise ValueError('Empty catalogue')
    return scores

def block_quantile(scores,alpha_numerator,alpha_denominator):
    if not scores or not 0<alpha_numerator<alpha_denominator:raise ValueError('Invalid calibration')
    n=len(scores);rank=((n+1)*(alpha_denominator-alpha_numerator)+alpha_denominator-1)//alpha_denominator
    ordered=sorted(math.inf if s is None else s for s in scores)
    q=math.inf if rank>n else ordered[rank-1]
    return {'n_blocks':n,'rank':rank,'q':None if math.isinf(q) else int(q),'infinity':math.isinf(q),'alpha':[alpha_numerator,alpha_denominator]}

def prediction_set(scores,quantile):
    if quantile['infinity']:return sorted(scores)
    return sorted(a for a,c in scores.items() if c is not None and c<=quantile['q'])
