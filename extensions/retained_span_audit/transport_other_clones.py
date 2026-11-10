"""Fixed computational transport to all 17 other deposited COL17A1 clones.
Same-study held-out reads, not independent biological truth or external-study validation.
Targets/references from Clone10 analysis; primary padding9/deficit36 fixed now.
Both padding0/9 and the original complete deficit grid are retained.
"""
import pathlib,json,hashlib,gzip,collections,re,datetime,time
from col17a1_clone10_audit import fetch,parse,rc
from recompute_panel import encode,classify_intervals,maximize,NEG
ROOT=pathlib.Path(__file__).resolve().parent
OUT=ROOT/'TRANSPORT_MP14'
def sha(b):return hashlib.sha256(b).hexdigest()
def clip(s,l,r):
    if s.count(l)!=1 or s.count(r)!=1:return None
    a,b=s.index(l),s.index(r)
    return s[a:b+len(r)] if b>=a+len(l) else None
def download(f):
    for attempt in range(3):
        try:
            b=fetch(f['links']['self'])
            assert len(b)==f['size'] and 'md5:'+hashlib.md5(b).hexdigest()==f['checksum']
            return b
        except Exception:
            if attempt==2:raise
            time.sleep(2)
def main():
    OUT.mkdir(exist_ok=True)
    freeze_path=OUT/'PROTOCOL.json'
    if freeze_path.exists():
        protocol=json.loads(freeze_path.read_text())
    else:
        meta=json.loads(fetch('https://zenodo.org/api/records/6242279'));groups={}
        for f in meta['files']:
            m=re.fullmatch(r'C17Ex52_(DNA|cDNA) immoJEB_PPN_Clone(\d+)_(fwd|rev).fastq.gz',f['key'])
            if m and m[2]!='10':groups.setdefault(m[2],[]).append(f)
        assert len(groups)==17 and all(len(v)==4 for v in groups.values())
        protocol={'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':__doc__,
          'primary':{'padding':9,'deficit':36},'padding':[0,9],'deltas':[0,9,18,27,36,45,54],
          'groups':groups,'parent_freeze_sha256':sha((ROOT/'P1_FREEZE.json').read_bytes()),
          'engine_sha256':sha((ROOT/'recompute_panel.py').read_bytes()),
          'script_sha256':sha(pathlib.Path(__file__).read_bytes()),
          'comparison':'R1 P1 call versus separately observed R2 exact-flank local string; disagreement is not independently adjudicated biological error.',
          'retention':'All raw pair denominators retained. All17 clones reported. No parameter changes based on new read outcomes.'}
        freeze_path.write_text(json.dumps(protocol,indent=2)+'\n')
    assert sha((ROOT/'P1_FREEZE.json').read_bytes())==protocol['parent_freeze_sha256']
    assert sha((ROOT/'recompute_panel.py').read_bytes())==protocol['engine_sha256']
    assert sha(pathlib.Path(__file__).read_bytes())==protocol['script_sha256']
    freeze=json.loads((ROOT/'P1_FREEZE.json').read_text());refs=freeze['references'];targets=freeze['targets']
    L,R=freeze['outer_left'],freeze['outer_right'];pads=protocol['padding'];ds=protocol['deltas']
    cache={};cachepath=OUT/'SEGMENTS.json'
    if cachepath.exists():cache=json.loads(cachepath.read_text())
    # Reuse deterministic baseline results only when query strings and frozen rules match.
    old=json.loads((ROOT/'BOUNDARY_SEGMENTS.json').read_text())
    def evaluate(q):
        if q is None:return {str(p):{'calls':[['U']*len(ds) for _ in targets]} for p in pads}
        h=sha(q.encode())
        if h in cache:assert cache[h]['sequence']==q;return cache[h]['results']
        if h in old:ans={str(p):old[h][str(p)] for p in pads}
        else:
            ans={}
            for p in pads:
                ts={k:(L[-p:] if p else '')+v+R[:p] for k,v in targets.items()}
                labels=classify_intervals(q,ts);deficits=[]
                for ref in refs:
                    opt,mx,_,_,suffix=maximize(encode(ref['sequence']),encode(q),ref['lo']-p,ref['hi']+p,labels,max(ds))
                    assert opt==suffix[0,0,0] and all(max(x)==opt for x in mx)
                    deficits.append([[None if v==NEG else int(opt-v) for v in row] for row in mx])
                calls=[]
                for ti in range(2):
                    cc=[]
                    for d in ds:
                        sets=[{k for k,v in zip(('NE','E','U'),z[ti]) if v is not None and v<=d} for z in deficits]
                        cc.append(next(iter(sets[0])) if len(sets[0])==1 and sets[0]==sets[1] else 'U')
                    for a,b in zip(cc,cc[1:]):assert a==b or b=='U'
                    calls.append(cc)
                ans[str(p)]={'calls':calls,'deficits':deficits}
        cache[h]={'sequence':q,'results':ans};return ans
    completed=[]
    for clone in sorted(protocol['groups'],key=int):
        resultpath=OUT/f'clone{clone}_result.json';ledgerpath=OUT/f'clone{clone}_pairs.jsonl.gz'
        if resultpath.exists():
            v=json.loads(resultpath.read_text());assert sha(ledgerpath.read_bytes())==v['ledger_sha256']
            completed.append(v);print('REUSED',clone,flush=True);continue
        reads={};verified=[]
        for f in protocol['groups'][clone]:
            b=download(f);reads[f['key']]=parse(b);verified.append(dict(f,sha256=sha(b)))
        counts=collections.Counter();cross=collections.Counter();admission=collections.Counter();totals={};literal=collections.Counter()
        with gzip.open(ledgerpath,'wt',encoding='ascii') as out:
            for typ in ['DNA','cDNA']:
                a=reads[f'C17Ex52_{typ} immoJEB_PPN_Clone{clone}_fwd.fastq.gz'];b=reads[f'C17Ex52_{typ} immoJEB_PPN_Clone{clone}_rev.fastq.gz']
                assert len(a)==len(b);totals[typ]=len(a)
                for (h,s),(h2,t) in zip(a,b):
                    assert h==h2
                    qs=[clip(s,L,R),clip(rc(t),L,R)];ev=[evaluate(q) for q in qs]
                    admission[typ,sum(q is not None for q in qs)]+=1
                    direct=[['U' if q is None or 'N' in q[len(L):-len(R)] else ('E' if q[len(L):-len(R)]==target else 'NE') for target in targets.values()] for q in qs]
                    for ti,name in enumerate(targets):literal[typ,name,direct[0][ti],direct[1][ti]]+=1
                    for p in pads:
                        for ti,name in enumerate(targets):
                            for di,d in enumerate(ds):
                                x,y=(e[str(p)]['calls'][ti][di] for e in ev)
                                call=x if x==y else 'U';counts[typ,p,name,d,call]+=1
                                if x!='U' and y!='U' and x!=y:counts[typ,p,name,d,'mate_conflict']+=1
                                cross[typ,p,name,d,x,direct[1][ti]]+=1
                    out.write(json.dumps({'id':h,'library':typ,'segment_hashes':[sha(q.encode()) if q else None for q in qs],'literal_mate_labels':direct})+'\n')
        v={'clone':clone,'protocol_sha256':sha(freeze_path.read_bytes()),'totals':totals,'files':verified,
           'ledger_sha256':sha(ledgerpath.read_bytes()),
           'admission':[{'library':t,'admitted_mates':n,'pairs':v} for (t,n),v in sorted(admission.items())],
           'counts':[{'library':t,'padding':p,'target':a,'delta':d,'call':c,'pairs':v} for (t,p,a,d,c),v in sorted(counts.items())],
           'r1_vs_r2_literal':[{'library':t,'padding':p,'target':a,'delta':d,'r1':x,'r2_literal':y,'pairs':v} for (t,p,a,d,x,y),v in sorted(cross.items())],
           'literal_pairs':[{'library':t,'target':a,'r1':x,'r2':y,'pairs':v} for (t,a,x,y),v in sorted(literal.items())]}
        for typ,n in totals.items():
            for p in pads:
                for a in targets:
                    for d in ds:
                        assert sum(counts[typ,p,a,d,c] for c in ['E','NE','U'])==n
                        assert sum(cross[typ,p,a,d,x,y] for x in ['E','NE','U'] for y in ['E','NE','U'])==n
        cachepath.write_text(json.dumps(cache,separators=(',',':')))
        resultpath.write_text(json.dumps(v,indent=2)+'\n');completed.append(v)
        print('COMPLETE',clone,totals,'PRIMARY_POSITIVE',[(t,a,counts[t,9,a,36,'E']) for t in totals for a in targets],'segments',len(cache),flush=True)
    assert len(completed)==17
    summary={'protocol_sha256':sha(freeze_path.read_bytes()),'clones':17,'raw_pairs':sum(sum(v['totals'].values()) for v in completed),
        'segments':len(cache),'status':'COMPLETE_FIXED_SAME_STUDY_TRANSPORT','result_files':[f'clone{c}_result.json' for c in sorted(protocol['groups'],key=int)]}
    (OUT/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
