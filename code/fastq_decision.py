"""Stream FASTQ into an exact-anchor DNA-event report. Python standard library only.

This declared observer is not an alignment-domain certificate or molecular truth.
SPDX-License-Identifier: GPL-3.0-only
"""
from pathlib import Path
from collections import Counter
import argparse,gzip,hashlib,json
from observability_report import analyse,render,require
ROOT=Path(__file__).resolve().parents[1]
COMP=str.maketrans('ACGTN','TGCAN')

def validate(c):
    for k in ['left_anchor','right_anchor','exact_insert']:
        require(isinstance(c.get(k),str) and bool(c[k]) and set(c[k])<=set('ACGT'),'Specify nonempty ACGT '+k)
    require(type(c.get('max_span')) is int and c['max_span']>=len(c['exact_insert']),'Invalid max_span')
    require(type(c.get('min_mean_phred')) is int and 0<=c['min_mean_phred']<=93,'Invalid Phred threshold')
    require(isinstance(c.get('reference'),str) and bool(c['reference']),'Specify reference')
    return c

def positions(s,m):
    i=s.find(m)
    while i>=0:
        yield i
        i=s.find(m,i+1)

def classify(s,q,c):
    if sum(ord(x)-33 for x in q)<c['min_mean_phred']*len(q):return 'U','quality'
    hits=[]
    for seq in [s,s.translate(COMP)[::-1]]:
        for a in positions(seq,c['left_anchor']):
            lo=a+len(c['left_anchor'])
            for b in positions(seq,c['right_anchor']):
                if lo<=b<=lo+c['max_span']:hits.append(seq[lo:b])
    if len(hits)!=1:return 'U','missing_anchor_pair' if not hits else 'multiple_anchor_pairs'
    if set(hits[0])-set('ACGT'):return 'U','ambiguous_base'
    return ('E' if hits[0]==c['exact_insert'] else 'NE'),'unique_observed_pair'

def count(path,c):
    validate(c);counts=Counter();reasons=Counter();seen=set();stream=hashlib.sha256()
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt',encoding='ascii',newline='') as f:
        while h:=f.readline():
            s=f.readline().rstrip('\r\n');plus=f.readline().rstrip('\r\n');q=f.readline().rstrip('\r\n')
            require(h.startswith('@') and plus.startswith('+') and s and len(s)==len(q),'Malformed four-line FASTQ')
            rid=h[1:].split()[0];require(rid not in seen,'Duplicate read identifier: '+rid);seen.add(rid)
            require(all(33<=ord(x)<=126 for x in q),'Invalid Phred+33 quality')
            require(set(s)<=set('ACGTN'),'Only uppercase ACGTN input is supported; normalize explicitly before analysis')
            if len(plus)>1:require(plus[1:].split()[0]==rid,'FASTQ separator identifier mismatch')
            stream.update((s+'\n'+q+'\n').encode('ascii'))
            label,reason=classify(s,q,c);counts[label]+=1;reasons[reason]+=1
    require(bool(seen),'Empty FASTQ')
    with path.open('rb') as f:
        h=hashlib.sha256()
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return {'sample_id':path.name,'universe_id':stream.hexdigest(),'source_locator':path.name,'N':len(seen),**{k:counts[k] for k in ['E','NE','U']},'compressed_or_plain_sha256':h.hexdigest(),'reasons':dict(reasons)}

def run(t,c,config,out):
    validate(config);require(not out.exists(),'Output already exists; use a new directory')
    require(t.resolve()!=c.resolve(),'Treatment and control must be different sources')
    data={'schema_version':'1.0','comparability':'same_declared_contract','contract':{
        'target':'Exact observed insert between declared anchors: '+config['exact_insert'],'observation_unit':'unique FASTQ record within each library',
        'reference':config['reference'],'coordinates':'Sequence strictly between ordered exact anchors; both orientations; maximum span '+str(config['max_span']),
        'quality_policy':'Whole-read arithmetic mean Phred >= '+str(config['min_mean_phred']),
        'candidate_domain':'All exact anchor pairs in either orientation; require exactly one pair; no alignment-domain certification',
        'accepted_label_assumption':'Accepted E/NE labels are correct under this observer; U may be either value','decision_timing':config.get('decision_timing','User-specified analysis; no prospective status inferred')},
        'observer_config':config,'treated':count(t,config),'control':count(c,config)}
    result=analyse(data);out.mkdir(parents=True)
    for name,value in [('input.json',data),('report.json',result)]: (out/name).write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')
    (out/'report.html').write_text(render(result),encoding='utf-8')
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--treated',required=True,type=Path);p.add_argument('--control',required=True,type=Path);p.add_argument('--config',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
    r=run(a.treated,a.control,json.loads(a.config.read_bytes()),a.output);print(json.dumps(r['contrast']))
if __name__=='__main__':main()
