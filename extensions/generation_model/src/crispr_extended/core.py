"""Exact min-plus DP over finite generation histories.

Positive-length source fragments; free source start/end; unit substitutions,
insertions and deletions; all gap transitions admitted. Insertions belong to
the boundary before source coordinate i, including at fragment endpoints.
No score-optimal traceback or biological genotype is used as input.
"""
import numpy as np
from numba import njit

@njit(cache=True)
def _advance(s, i, b, consumes, lo, hi, target):
    # s=0 unavailable; 1..L+1 exact prefix; L+2 mismatch;
    # L+3 completed E; L+4 completed NE. i is next source coordinate.
    L=len(target)
    if consumes and i==lo-1:
        return 1 if b<4 else 0
    if lo<=i<hi:
        if b>=4:return 0
        if 1<=s<=L+1:
            k=s-1
            return s+1 if k<L and b==target[k] else L+2
        return s
    if consumes and i==hi:
        if b>=4:return 0
        if s==L+1:return L+3
        if 2<=s<=L+2:return L+4
        return 0
    return s

@njit(cache=True)
def _dp(h,read,lo,hi,target,cap):
    n=len(h);m=len(read);L=len(target);K=L+5;INF=cap+1
    row=np.full((m+1,K),INF,np.int32)
    answer=np.full(3,INF,np.int32)  # NE, E, U
    for i in range(n+1):
        nxt=np.full((m+1,K),INF,np.int32)
        for j in range(m+1):
            # A new positive-source-length history can begin here, preceded
            # by j insertions. It becomes active only after consuming a base.
            if i<n and j<=cap:
                if j+1<=cap and j+1<nxt[j,0]:nxt[j,0]=j+1
                if j<m:
                    cost=j+int(h[i]!=read[j]);ns=_advance(0,i,read[j],True,lo,hi,target)
                    if cost<nxt[j+1,ns]:nxt[j+1,ns]=cost
            for s in range(K):
                cost=row[j,s]
                if cost>cap:continue
                if j==m:
                    label=1 if s==L+3 else 0 if s==L+4 else 2
                    if cost<answer[label]:answer[label]=cost
                if j<m and cost<cap:
                    ns=_advance(s,i,read[j],False,lo,hi,target)
                    if cost+1<row[j+1,ns]:row[j+1,ns]=cost+1
                if i<n:
                    # Deletion of either required source flank loses it.
                    ns=0 if i==lo-1 or i==hi else s
                    if cost<cap and cost+1<nxt[j,ns]:nxt[j,ns]=cost+1
                    if j<m:
                        ns=_advance(s,i,read[j],True,lo,hi,target)
                        c=cost+int(h[i]!=read[j])
                        if c<nxt[j+1,ns]:nxt[j+1,ns]=c
        row=nxt
    return answer

def _encode(s,reference=False):
    if not isinstance(s,str) or not s or set(s)-set('ACGT' if reference else 'ACGTN'):
        raise ValueError('Nonempty uppercase ACGT source/target or ACGTN read required')
    return np.array(['ACGTN'.index(x) for x in s],dtype=np.int8)

def minimum_costs(sequence,lo,hi,read,target,max_errors):
    if type(max_errors) is not int or not 0<=max_errors<=1000:raise ValueError('Invalid error cap')
    if type(lo) is not int or type(hi) is not int or not 0<lo<hi<len(sequence):raise ValueError('Both source flanks must exist in template')
    costs=_dp(_encode(sequence,True),_encode(read),lo,hi,_encode(target,True),max_errors)
    return {label:(int(c) if c<=max_errors else None) for label,c in zip(('NE','E','U'),costs)}

def certify(read,target,haplotypes,budgets=(0,1,2,3,5,8,12)):
    budgets=tuple(budgets)
    if not budgets or any(type(b) is not int or b<0 for b in budgets) or list(budgets)!=sorted(set(budgets)):raise ValueError('Budgets must be unique, sorted nonnegative integers')
    if not haplotypes:raise ValueError('Empty candidate catalogue')
    ids=[h['id'] for h in haplotypes]
    if len(set(ids))!=len(ids):raise ValueError('Duplicate haplotype identifiers')
    domains=[]
    for h in haplotypes:
        costs=minimum_costs(h['sequence'],h['lo'],h['hi'],read,target,max(budgets))
        domains.append({'id':h['id'],'allele':h.get('allele',h['id']),'minimum_errors':costs})
    decisions=[]
    for b in budgets:
        labels=sorted({lab for d in domains for lab,c in d['minimum_errors'].items() if c is not None and c<=b})
        alleles=sorted({str(d['allele']) for d in domains if any(c is not None and c<=b for c in d['minimum_errors'].values())})
        call=labels[0] if len(labels)==1 and labels[0]!='U' else 'U'
        reason='CERTIFIED' if call!='U' else 'NO_ADMISSIBLE_HISTORY' if not labels else 'SOURCE_FLANK_OR_BASE_UNRESOLVED_POSSIBLE' if 'U' in labels else 'EVENT_DISAGREEMENT'
        decisions.append({'budget':b,'possible_observed_labels':labels,'observed_tract_certificate':call,'reason':reason,'possible_source_alleles':alleles,'source_allele_certificate':alleles[0] if len(alleles)==1 else None})
    return {'model':'bounded-edit-fragment-v1','domains':domains,'decisions':decisions,'scope':'Conditional on candidate templates, local source domain and absolute unit-edit bound; no calibrated error probability.'}
