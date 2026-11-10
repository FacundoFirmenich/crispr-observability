"""Retrospective P1 application, original MP07 engine, individual paired reads.
References: modal R1 strings of independent healthy/JEB control libraries.
Targets: two previously observed Clone10 local strings (post hoc selection).
No independent genotype truth, therapeutic efficacy or prospective validation.
"""
import sys, pathlib, json, hashlib, collections, gzip
ROOT=pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent))
from col17a1_clone10_audit import fetch, parse, rc
from recompute_panel import encode, classify_intervals, maximize, NEG
import optimal_lattice as lattice
from recompute_panel import solve

DELTAS=(0,9,18,27,36,45,54)
LEFT='GATGCCTCCCACAGTC'
RIGHT='GTCAGGCGGGGCAGCTC'
# Common target-flanking interval, within both mate coverages. A first launch
# with wider amplicon flanks was interrupted before a freeze/result existed:
# those flanks would not fit both different DNA/cDNA amplicon endpoints.
OUTER_L=LEFT
OUTER_R=RIGHT

def clip(s,l,r):
    if s.count(l)!=1 or s.count(r)!=1:return None
    a,b=s.index(l),s.index(r)
    return s[a:b+len(r)] if b>=a+len(l) else None
def sha(b):return hashlib.sha256(b).hexdigest()
def main():
    refs=[];controls=[]
    for typ in ['HKC','JEB']:
        p=ROOT/f'C17Ex52_DNA_immo{typ}_Ctl_fwd.fastq.gz.top.json'
        j=json.loads(p.read_text());s=clip(j['top'][0][0],OUTER_L,OUTER_R)
        assert s is not None
        lo=s.index(LEFT)+len(LEFT);hi=s.index(RIGHT)
        refs.append({'name':typ,'sequence':s,'lo':lo,'hi':hi,'local':s[lo:hi]})
        controls.append({'file':p.name,'sha256':sha(p.read_bytes()),'source':j})
    targets={'long40':'ACAGCTCATCGGGGTAGCAGCTCCTCCTCACACAGCTCAT',
             'short18':'GGGGTAGCAGCTCCTCAT'}
    assert [len(r['local']) for r in refs]==[32,30]
    assert [len(x) for x in targets.values()]==[40,18]
    freeze={'scope':__doc__,'references':refs,'targets':targets,'outer_left':OUTER_L,'outer_right':OUTER_R,
       'deltas':DELTAS,'controls':controls,'engine_sha256':sha((ROOT/'recompute_panel.py').read_bytes()),
       'contract':'E means exact local target string; NE means another resolved string; U includes non-admission and alignment ambiguity. Consensus requires same determinate call across both references and both mates.',
       'annotation_gate':'Published S2/S6 unavailable; no claim of author-genotype confirmation.',
       'selection':'Clone10 targets and all anchor choices retrospective; no independent accuracy estimate.'}
    (ROOT/'P1_FREEZE.json').write_text(json.dumps(freeze,indent=2)+'\n')
    meta=json.loads(fetch('https://zenodo.org/api/records/6242279'))
    files=sorted([f for f in meta['files'] if 'Clone10_' in f['key']],key=lambda f:f['key'])
    assert len(files)==4
    reads={};verified=[]
    for f in files:
        b=fetch(f['links']['self']);assert len(b)==f['size'] and 'md5:'+hashlib.md5(b).hexdigest()==f['checksum']
        reads[f['key']]=parse(b);verified.append(dict(f,sha256=sha(b)))
        print('INPUT',f['key'],len(reads[f['key']]),flush=True)
    cache={};witnesses=[]
    def evaluate(q):
        if q is None:return [['U']*len(DELTAS) for _ in targets]
        if q in cache:return cache[q]['calls']
        classes=classify_intervals(q,targets);deficits=[]
        for ref in refs:
            opt,mx,_,_,suffix=maximize(encode(ref['sequence']),encode(q),ref['lo'],ref['hi'],classes,max(DELTAS))
            assert opt==suffix[0,0,0] and all(max(x)==opt for x in mx)
            deficits.append([[None if v==NEG else int(opt-v) for v in row] for row in mx])
        calls=[]
        for ti in range(2):
            tc=[]
            for delta in DELTAS:
                sets=[{k for k,v in zip(('NE','E','U'),d[ti]) if v is not None and v<=delta} for d in deficits]
                tc.append(next(iter(sets[0])) if len(sets[0])==1 and sets[0]==sets[1] else 'U')
            # Increasing admissible alignments may revoke but cannot recover certification.
            for a,b in zip(tc,tc[1:]):assert a==b or b=='U'
            calls.append(tc)
        cache[q]={'calls':calls,'deficits':deficits}
        return calls
    stats=collections.Counter();admission=collections.Counter();totals={}
    ledger=ROOT/'P1_PAIRS.jsonl.gz'
    with gzip.open(ledger,'wt',encoding='ascii') as out:
        for typ in ['DNA','cDNA']:
            a=reads[f'C17Ex52_{typ} immoJEB_PPN_Clone10_fwd.fastq.gz']
            b=reads[f'C17Ex52_{typ} immoJEB_PPN_Clone10_rev.fastq.gz']
            assert len(a)==len(b);totals[typ]=len(a)
            for idx,((h,s),(h2,t)) in enumerate(zip(a,b)):
                assert h==h2
                qs=[clip(s,OUTER_L,OUTER_R),clip(rc(t),OUTER_L,OUTER_R)]
                cs=[evaluate(q) for q in qs]
                admission[typ,sum(q is not None for q in qs)]+=1
                pair=[]
                for ti,name in enumerate(targets):
                    cc=[]
                    for di,delta in enumerate(DELTAS):
                        x,y=cs[0][ti][di],cs[1][ti][di]
                        call=x if x==y else 'U';cc.append(call)
                        stats[typ,name,delta,call]+=1
                        if x!='U' and y!='U' and x!=y:stats[typ,name,delta,'mate_conflict']+=1
                    pair.append(cc)
                out.write(json.dumps({'library':typ,'id':h,'segment_hashes':[sha(q.encode()) if q else None for q in qs],
                                      'mate_calls':cs,'pair_calls':pair})+'\n')
                if idx%5000==0:print('PROGRESS',typ,idx,'unique_segments',len(cache),flush=True)
    # Full scored/path-projection witness checks on both targets, each reference.
    for name,target in targets.items():
        q=OUTER_L+refs[1]['sequence'][len(OUTER_L):refs[1]['lo']]+target+refs[1]['sequence'][refs[1]['hi']:]
        assert q in cache, 'Observed target segment required for witness check'
        for ref in refs:
            case={'reference_sequence':ref['sequence'],'observed_sequence':q,'lo':ref['lo'],'hi':ref['hi'],
                  'target_kind':'HEK3_DNA','programmed_local':target,'unedited_local':refs[1]['local']}
            opt,ws=solve(case,lattice,cutoff=54)
            witnesses.append({'target':name,'reference':ref['name'],'optimum':opt,'witnesses':ws})
    result={'freeze_sha256':sha((ROOT/'P1_FREEZE.json').read_bytes()),'files':verified,'totals':totals,
       'admission':[{'library':t,'admitted_mates':n,'pairs':v} for (t,n),v in sorted(admission.items())],
       'counts':[{'library':t,'target':n,'delta':d,'call':c,'pairs':v} for (t,n,d,c),v in sorted(stats.items())],
       'unique_segments':len(cache),'ledger_sha256':sha(ledger.read_bytes()),'witness_checks':len(witnesses)}
    for typ,n in totals.items():
        for target in targets:
            for d in DELTAS:assert sum(stats[typ,target,d,c] for c in ['E','NE','U'])==n
    for filename,value in [('P1_RESULT.json',result),('P1_SEGMENTS.json',cache),('P1_WITNESSES.json',witnesses)]:
        (ROOT/filename).write_text(json.dumps(value,indent=2)+'\n')
    print(json.dumps(result),flush=True)
if __name__=='__main__':main()
