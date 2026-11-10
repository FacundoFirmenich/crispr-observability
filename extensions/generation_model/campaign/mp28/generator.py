"""MIT. Independent source and edit-ledger generator; no predictor imports."""
import hashlib,random

def source_templates(reference,case):
    v=case['indel'];p=v['pos'];ref=v['ref'];assert reference[p:p+len(ref)]==ref
    out=[]
    for allele,variant in enumerate((ref,v['alts'][0])):
        seq=reference[p-300:p]+variant+reference[p+len(ref):p+len(ref)+300]
        end=300+len(variant);b=list(seq)
        for i in (275,end+25):b[i]='CGTA'['ACGT'.index(b[i])]
        nuisance=''.join(b[:end+40])+'GCA'+''.join(b[end+40:])
        out.extend([{'id':f'canonical_{allele}','allele':str(allele),'sequence':seq,'background':'canonical','event_start':300,'event_end':end},
                    {'id':f'nuisance_{allele}','allele':str(allele),'sequence':nuisance,'background':'nuisance','event_start':300,'event_end':end}])
    return out

def generate(template,case,regime,geometry,replicate,protocol):
    seq=template['sequence'];lo=template['event_start'];hi=template['event_end'];length=protocol['length']
    start={'center':lo-(length-(hi-lo))//2,'left_event_absent':lo+1,'right_event_absent':hi-1-length}[geometry]
    stop=start+length;assert 0<=start<stop<=len(seq)
    identity=f"MP28_V1|{case['id']}|{template['allele']}|{regime}|{geometry}|{replicate}"
    seed=int(hashlib.sha256(identity.encode()).hexdigest()[:16],16);rng=random.Random(seed)
    noise=protocol['noise']['shifted' if regime=='shifted' else 'baseline'];bases=[];edits=[]
    for i in range(start,stop):
        if rng.random()<noise['insertion']:
            b=rng.choice('ACGT');bases.append(b);edits.append(['I',i,b])
        if rng.random()<noise['deletion']:edits.append(['D',i,seq[i]]);continue
        b=seq[i]
        if rng.random()<noise['substitution']:
            b=rng.choice([x for x in 'ACGT' if x!=b]);edits.append(['X',i,seq[i],b])
        bases.append(b)
    return {'id':hashlib.sha256(identity.encode()).hexdigest(),'site':case['id'],'allele':template['allele'],'regime':regime,'geometry':geometry,'replicate':replicate,'seed':seed,'template_id':template['id'],'source_start':start,'source_end':stop,'sequence':''.join(bases),'edits':edits}
