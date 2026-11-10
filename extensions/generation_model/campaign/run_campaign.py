"""Frozen, resumable CPU validation. Writes only its own Kaggle output root."""
import bisect,collections,gzip,hashlib,json,os,platform,shutil,sys,time,traceback
import concurrent.futures,multiprocessing
from pathlib import Path
HERE=Path(__file__).resolve().parent;REPO=HERE.parent;sys.path.insert(0,str(REPO/'src'))
from crispr_extended import certify
from heldout_generator import generate
OUT=Path(os.environ.get('EXTENDED_OUTPUT','/kaggle/working/EXTENDED_MP27'))
P=json.loads((REPO/'protocols/MP27.json').read_text());B=json.loads((HERE/'input_binding.json').read_text())
def digest(p):
    with open(p,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
PH=hashlib.sha256(json.dumps(P,sort_keys=True,separators=(',',':')).encode()).hexdigest()
ENGINE=digest(REPO/'src/crispr_extended/core.py')
def save(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix(path.suffix+'.tmp')
    with tmp.open('w') as f:json.dump(obj,f,indent=2);f.flush();os.fsync(f.fileno())
    tmp.replace(path)
def source_roots():
    roots={}
    for p in Path('/kaggle/input').rglob('MANIFEST.json'):
        h=digest(p)
        for name in ('MP24','MP25'):
            if h==B[name]['manifest_sha256']:
                assert name not in roots;roots[name]=p.parent
    assert set(roots)=={'MP24','MP25'}
    return roots
def verified(root,relative):
    m=json.loads((root/'MANIFEST.json').read_text());row=next(r for r in m['files'] if r['path']==relative);p=root/relative
    assert p.resolve().is_relative_to(root.resolve()) and digest(p)==row['sha256'];return p
def select_heldout(root,excluded):
    import pysam
    vcf=verified(root,'inputs/benchmark.vcf.gz');bed=verified(root,'inputs/benchmark.bed')
    intervals=sorted(tuple(map(int,l.split()[1:3])) for l in bed.read_text().splitlines() if l.startswith('chr20\t'));starts=[x[0] for x in intervals]
    allrows=[]
    with pysam.VariantFile(vcf) as vs:
        for v in vs:
            if v.chrom=='chr20':allrows.append({'pos':v.start,'ref':v.ref,'alts':list(v.alts or []),'filter':list(v.filter)})
    positions=[v['pos'] for v in allrows];groups=collections.defaultdict(list)
    for v in allrows:
        if v['pos'] in excluded or len(v['alts'])!=1 or v['filter'] not in ([],['PASS']):continue
        ref=v['ref'];alt=v['alts'][0];n=abs(len(ref)-len(alt))
        if min(len(ref),len(alt))!=1 or not 2<=n<=30 or set(ref+alt)-set('ACGT'):continue
        lo=v['pos']-300;hi=v['pos']+len(ref)+300;k=bisect.bisect_right(starts,lo)-1
        if k<0 or not intervals[k][0]<=lo<hi<=intervals[k][1]:continue
        neighbors=allrows[bisect.bisect_left(positions,lo):bisect.bisect_left(positions,hi)]
        if len(neighbors)!=1:continue
        group=('ins' if len(alt)>len(ref) else 'del')+'_'+('2-5' if n<=5 else '6-10' if n<=10 else '11-30')
        case={'id':f"chr20_{v['pos']+1}",'stratum':group,'indel':v,'rank':hashlib.sha256(f"EXTENDED_MP27_HOLDOUT_V1|{v['pos']}".encode()).hexdigest()};groups[group].append(case)
    selected=[c for g in sorted(groups) for c in sorted(groups[g],key=lambda c:c['rank'])[:20]]
    assert selected
    selection={'selected':selected,'available':{g:len(x) for g,x in groups.items()},'excluded_positions':sorted(excluded),'selection_before_new_outcomes':True}
    save(OUT/'HELDOUT_SELECTION.json',selection);return selected
def templates(reference,case,pad):
    v=case['indel'];p=v['pos'];ref=v['ref'];alt=v['alts'][0];assert reference[p:p+len(ref)]==ref
    rr=reference[p-300:p+len(ref)+300];aa=rr[:300]+alt+rr[300+len(ref):]
    return reference[p-pad:p]+alt+reference[p+len(ref):p+len(ref)+pad],[{'id':'reference','allele':'0','sequence':rr,'lo':300-pad,'hi':300+len(ref)+pad},{'id':'alternate','allele':'1','sequence':aa,'lo':300-pad,'hi':300+len(alt)+pad}]
def truth(r,pad):
    lo=r['source_event_start']-pad;hi=r['source_event_end']+pad
    survivors={i for i,ins in zip(r['source_origins'],r['inserted']) if not ins}
    span=''.join(b for b,i in zip(r['sequence'],r['source_origins']) if lo<=i<hi)
    valid=lo-1 in survivors and hi in survivors and bool(span) and not set(span)-set('ACGT')
    return span if valid else None
def verify_generated(r,reference,case):
    v=case['indel'];p=v['pos'];ref=v['ref'];variant=v['alts'][0] if r['allele'] else ref;c=r['source_event_start']
    src=reference[p-c:p]+variant+reference[p+len(ref):p+len(ref)+c]
    ins=collections.defaultdict(list);sub={};deleted=set()
    for e in r['edits']:
        assert r['source_start']<=e[1]<r['source_end']
        if e[0]=='I':ins[e[1]].append(e[2])
        elif e[0]=='D':assert src[e[1]]==e[2];deleted.add(e[1])
        else:assert e[0]=='X' and src[e[1]]==e[2] and e[3]!=e[2];sub[e[1]]=e[3]
    bases=[];origins=[];inserted=[]
    for i in range(r['source_start'],r['source_end']):
        for b in ins[i]:bases.append(b);origins.append(i);inserted.append(True)
        if i not in deleted:bases.append(sub.get(i,src[i]));origins.append(i);inserted.append(False)
    assert ''.join(bases)==r['sequence'] and origins==r['source_origins'] and inserted==r['inserted']
def generated_records(reference,case):
    for allele in (0,1):
        for length in (100,150,200):
            for geometry in ('center','left_flank_absent','right_flank_absent'):
                for model in ('clean','sparse_substitutions','clustered_indels'):
                    for rep in range(1 if model=='clean' else 2):yield generate(reference,case,allele,length,geometry,model,rep)
def existing_records(root,case,real):
    path=verified(root,'sites/'+case['id']+('/READS.jsonl.gz' if real else '/RECORDS.jsonl.gz'))
    with gzip.open(path,'rt') as f:
        for line in f:yield json.loads(line)
def counter_json(c):return [{'key':list(k),'n':v} for k,v in sorted(c.items())]
def run_site(stage,case,iterator,reference):
    dest=OUT/stage/case['id'];dest.mkdir(parents=True,exist_ok=True);existing=dest/'SUMMARY.json'
    if existing.exists():
        s=json.loads(existing.read_text());assert s['protocol_sha256']==PH and s['engine_sha256']==ENGINE and digest(dest/'RECORDS.jsonl.gz')==s['records_sha256'];return s
    counts=collections.Counter();failures=[];n=0;tick=time.monotonic();real=stage=='real';temp=dest/'RECORDS.jsonl.gz.tmp'
    with gzip.open(temp,'wt') as out:
        for r in iterator:
            n+=1
            if not real:verify_generated(r,reference,case)
            results=[]
            for pad in P['paddings']:
                target,haps=templates(reference,case,pad);result=certify(r['sequence'],target,haps,P['budgets'])
                if real:observed='UNAVAILABLE';source_allele={'E':'1','NE':'0','UNKNOWN':'UNKNOWN'}[r['truth_alt']];errors=None
                else:
                    span=truth(r,pad);observed='U' if span is None else 'E' if span==target else 'NE';source_allele=str(r['allele']);errors=len(r['edits'])
                    if 'truth' in r:assert observed==r['truth'][str(pad)]['observed_label']
                for d in result['decisions']:
                    inside=not real and errors<=d['budget'];call=d['observed_tract_certificate'];ac=d['source_allele_certificate'] or 'U'
                    key=(pad,d['budget'],observed,call,source_allele,ac,str(inside),r.get('model','real'),r.get('geometry','real'))
                    counts[key]+=1
                    if inside and (observed not in d['possible_observed_labels'] or source_allele not in d['possible_source_alleles']):failures.append({'id':r['id'],'padding':pad,'budget':d['budget'],'kind':'TRUE_HISTORY_LABEL_EXCLUDED','observed':observed,'source_allele':source_allele,'decision':d})
                results.append({'padding':pad,'observed_truth':observed,'source_allele_truth':source_allele,'literal_error_count':errors,**result})
            assert [d['possible_source_alleles'] for d in results[0]['decisions']]==[d['possible_source_alleles'] for d in results[1]['decisions']]
            identity=r.get('id') or r['read_id']+':'+str(r['flag']&192)
            record={'id':identity,'site':case['id'],'results':results}
            if stage=='heldout':record['generated_observation']=r
            out.write(json.dumps(record,separators=(',',':'))+'\n')
    temp.replace(dest/'RECORDS.jsonl.gz')
    s={'protocol_sha256':PH,'engine_sha256':ENGINE,'stage':stage,'case':case,'records':n,'counts':counter_json(counts),'soundness_failures':failures,'records_sha256':digest(dest/'RECORDS.jsonl.gz'),'seconds':time.monotonic()-tick};save(existing,s)
    print('SITE_COMPLETE',stage,case['id'],n,len(failures),round(s['seconds'],2),flush=True);return s
def manifest():
    save(OUT/'MANIFEST.json',{'protocol_sha256':PH,'engine_sha256':ENGINE,'files':[{'path':p.relative_to(OUT).as_posix(),'sha256':digest(p),'MB':p.stat().st_size/1e6} for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='MANIFEST.json' and not p.name.endswith('.tmp')]})
def site_job(stage,case):
    iterator=generated_records(REFERENCE,case) if stage=='heldout' else existing_records(ROOTS['MP25' if stage=='real' else 'MP24'],case,stage=='real')
    return run_site(stage,case,iterator,REFERENCE)
def main():
    global REFERENCE,ROOTS
    started=time.monotonic();OUT.mkdir(parents=True,exist_ok=True)
    # Resume only exact protocol/engine checkpoints, never rename old results.
    for p in Path('/kaggle/input').rglob('EXTENDED_CHECKPOINT.json'):
        c=json.loads(p.read_text())
        if c.get('protocol_sha256')==PH and c.get('engine_sha256')==ENGINE:
            m=json.loads((p.parent/'MANIFEST.json').read_text())
            for row in m['files']:
                src=p.parent/row['path'];assert src.resolve().is_relative_to(p.parent.resolve()) and digest(src)==row['sha256'];dst=OUT/row['path'];dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst)
    try:
        roots=source_roots();sels={k:json.loads(verified(root,'SELECTION.json').read_text()) for k,root in roots.items()}
        reference_paths=[p for p in Path('/kaggle/input').rglob('chr20.fa.gz') if digest(p)==B['reference']['sha256']];assert reference_paths
        with gzip.open(reference_paths[0],'rt') as f:assert next(f).strip()=='>chr20';reference=''.join(l.strip() for l in f).upper()
        assert len(reference)==64444167
        excluded=set(sels['MP25']['excluded_prior_sites'])|{c['indel']['pos'] for c in sels['MP25']['selected']}
        heldout=select_heldout(roots['MP25'],excluded)
        save(OUT/'PROTOCOL.json',P);save(OUT/'INPUT_BINDING.json',B);save(OUT/'ENVIRONMENT.json',{'python':sys.version,'platform':platform.platform(),'cpu_count':os.cpu_count(),'numpy':__import__('numpy').__version__,'numba':__import__('numba').__version__,'engine_sha256':ENGINE})
        plan=[('development',c) for c in sels['MP24']['selected']]+[('real',c) for c in sels['MP25']['selected']]+[('heldout',c) for c in heldout]
        REFERENCE=reference;ROOTS=roots;summaries=[];jobs=iter(plan)
        # Warm the compiled kernel before fork; read-only reference is shared.
        certify('ACGT','CG',[{'id':'warmup','sequence':'ACGT','lo':1,'hi':3}],(0,1))
        workers=min(4,os.cpu_count() or 1)
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('fork')) as pool:
            active={}
            for stage,case in list(plan[:workers]):active[pool.submit(site_job,stage,case)]=(stage,case)
            jobs=iter(plan[workers:])
            while active:
                done,_=concurrent.futures.wait(active,return_when=concurrent.futures.FIRST_COMPLETED)
                for future in done:
                    stage,case=active.pop(future);s=future.result();summaries.append(s)
                    save(OUT/'EXTENDED_CHECKPOINT.json',{'protocol_sha256':PH,'engine_sha256':ENGINE,'completed_sites':len(summaries),'planned_sites':len(plan),'last':[stage,case['id']],'seconds':time.monotonic()-started})
                    if time.monotonic()-started<=9.5*3600:
                        job=next(jobs,None)
                        if job:active[pool.submit(site_job,*job)]=job
        counts=collections.Counter();failures=0;records=collections.Counter();nv=collections.Counter()
        for s in summaries:
            records[s['stage']]+=s['records'];failures+=len(s['soundness_failures'])
            for row in s['counts']:
                k=row['key'];counts[(s['stage'],)+tuple(k)]+=row['n']
                if s['stage']=='heldout' and k[0]==5 and k[1]==2 and k[7]=='clean' and k[8]=='center' and k[2]==k[3] and k[2] in ('E','NE'):nv[k[2]]+=row['n']
        complete=len(summaries)==len(plan)
        result={'status':'COMPLETE' if complete else 'CHECKPOINT_CONTINUABLE','protocol_sha256':PH,'engine_sha256':ENGINE,'records_by_stage':dict(records),'sites_completed':len(summaries),'sites_planned':len(plan),'soundness_failures':failures,'in_domain_soundness_gate':'PASS' if complete and failures==0 else 'FAIL' if failures else 'PENDING','nonvacuity_budget2_clean_center_padding5':dict(nv),'nonvacuity_gate':'PASS' if nv['E']>0 and nv['NE']>0 else 'FAIL' if complete else 'PENDING','seconds':time.monotonic()-started,'counts':counter_json(counts),'scope':'Conditional computational validation; no independent biological accuracy'}
        save(OUT/'RESULT.json',result);print('CAMPAIGN_RESULT',json.dumps({k:v for k,v in result.items() if k!='counts'}),flush=True)
    except Exception as e:
        save(OUT/'FAILURE.json',{'type':type(e).__name__,'message':str(e),'traceback':traceback.format_exc()});print(traceback.format_exc(),flush=True);raise
    finally:manifest()
if __name__=='__main__':main()
