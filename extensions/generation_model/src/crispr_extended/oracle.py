"""Independent exhaustive history enumeration for bounded small instances.
No production DP states or transition helpers are imported.
"""
def enumerate_histories(source,lo,hi,read,target,cap):
    best={'E':None,'NE':None,'U':None};histories=0
    for start in range(len(source)):
        def walk(i,j,cost,emitted):
            nonlocal histories
            if cost>cap:return
            if j==len(read) and i>start:
                histories+=1
                survivors={p for b,p,inserted in emitted if not inserted and b in 'ACGT'}
                span=''.join(b for b,p,inserted in emitted if lo<=p<hi)
                label='U' if lo-1 not in survivors or hi not in survivors or not span or set(span)-set('ACGT') else 'E' if span==target else 'NE'
                if best[label] is None or cost<best[label]:best[label]=cost
            if i<len(source):
                walk(i+1,j,cost+1,emitted)
                if j<len(read):walk(i+1,j+1,cost+(source[i]!=read[j]),emitted+[(read[j],i,False)])
            if j<len(read):walk(i,j+1,cost+1,emitted+[(read[j],i,True)])
        walk(start,0,0,[])
    return best,histories
