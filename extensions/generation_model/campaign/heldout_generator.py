"""MIT. Independent held-out simulation, not biological observations.
No production certifier, scoring or alignment imports.
"""
import hashlib,random
def generate(reference,case,allele,length,geometry,model,replicate):
    p=case['indel']['pos'];ref=case['indel']['ref'];alt=case['indel']['alts'][0];variant=alt if allele else ref
    source=reference[p-300:p]+variant+reference[p+len(ref):p+len(ref)+300]
    event=300;end=event+len(variant)
    start={'center':event-(length-len(variant))//2,'left_flank_absent':event-5,'right_flank_absent':end+5-length}[geometry];stop=start+length
    assert 0<=start<stop<=len(source)
    identity=f'EXTENDED_MP27_V1|{p}|{allele}|{length}|{geometry}|{model}|{replicate}'
    seed=int(hashlib.sha256(identity.encode()).hexdigest()[:16],16);rng=random.Random(seed)
    substitutions={};deletions=set();insertions={}
    if model=='sparse_substitutions':
        for i in range(start,stop):
            if rng.random()<0.008:substitutions[i]=rng.choice([x for x in 'ACGT' if x!=source[i]])
    elif model=='clustered_indels':
        center=rng.randrange(start,stop);width=rng.randint(2,5)
        for i in range(center,min(center+width,stop)):
            if rng.random()<0.55:deletions.add(i)
            if rng.random()<0.55:insertions[i]=''.join(rng.choice('ACGT') for _ in range(rng.randint(1,2)))
    else:assert model=='clean'
    bases=[];origins=[];inserted=[];edits=[]
    for i in range(start,stop):
        for b in insertions.get(i,''):
            bases.append(b);origins.append(i);inserted.append(True);edits.append(['I',i,b])
        if i in deletions:edits.append(['D',i,source[i]]);continue
        b=substitutions.get(i,source[i]);bases.append(b);origins.append(i);inserted.append(False)
        if b!=source[i]:edits.append(['X',i,source[i],b])
    return {'id':hashlib.sha256(identity.encode()).hexdigest(),'site':case['id'],'sequence':''.join(bases),'source_origins':origins,'inserted':inserted,'source_start':start,'source_end':stop,'source_event_start':event,'source_event_end':end,'allele':allele,'model':model,'geometry':geometry,'length':length,'replicate':replicate,'seed':seed,'edits':edits}
