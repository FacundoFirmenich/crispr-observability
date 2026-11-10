"""MIT. Frozen MP28 locus-block calibration and distribution-shift campaign."""
import bisect,collections,concurrent.futures,gzip,hashlib,json,multiprocessing,os,platform,shutil,sys,time,traceback
from pathlib import Path
HERE=Path(__file__).resolve().parent;R=HERE.parents[1];sys.path.insert(0,str(R/'src'))
from crispr_extended.calibration import allele_scores,block_quantile,prediction_set,fragment_distance
from generator import source_templates,generate
P=json.loads((R/'protocols/MP28.json').read_text());B=json.loads((HERE/'input_binding.json').read_text())
OUT=Path(os.environ.get('MP28_OUTPUT','/kaggle/working/MP28'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def object_sha(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
PH=object_sha(P);EH=object_sha({p.relative_to(R).as_posix():sha(p) for p in [HERE/'run.py',HERE/'generator.py',R/'src/crispr_extended/calibration.py',R/'src/crispr_extended/core.py']})
def save(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix(p.suffix+'.tmp')
    with t.open('w') as f:json.dump(x,f,indent=2);f.flush();os.fsync(f.fileno())
    t.replace(p)
def verified(root,name):
    row=next(x for x in json.loads((root/'MANIFEST.json').read_text())['files'] if x['path']==name)
    p=root/name;assert p.resolve().is_relative_to(root.resolve()) and sha(p)==row['sha256'];return p
def sources():
    roots={}
    for p in Path('/kaggle/input').rglob('MANIFEST.json'):
        h=sha(p)
        for name in ('MP25','MP27'):
            if h==B[name]['manifest_sha256']:roots[name]=p.parent
    assert set(roots)=={'MP25','MP27'}
    return roots
def selection(roots):
    import pysam
    old=json.loads(verified(roots['MP25'],'SELECTION.json').read_text())
    prior=json.loads(verified(roots['MP27'],'HELDOUT_SELECTION.json').read_text())
    excluded=set(old['excluded_prior_sites'])|{c['indel']['pos'] for c in old['selected']}|{c['indel']['pos'] for c in prior['selected']}
    bed=verified(roots['MP25'],'inputs/benchmark.bed');vcf=verified(roots['MP25'],'inputs/benchmark.vcf.gz')
    intervals=sorted(tuple(map(int,l.split()[1:3])) for l in bed.read_text().splitlines() if l.startswith('chr20\t'));starts=[x[0] for x in intervals]
    rows=[]
    with pysam.VariantFile(vcf) as vs:
        for v in vs:
            if v.chrom=='chr20':rows.append({'pos':v.start,'ref':v.ref,'alts':list(v.alts or []),'filter':list(v.filter)})
    positions=[v['pos'] for v in rows];pool=[]
    for v in rows:
        if v['pos'] in excluded or len(v['alts'])!=1 or v['filter'] not in ([],['PASS']):continue
        ref=v['ref'];alt=v['alts'][0];n=abs(len(ref)-len(alt))
        if min(len(ref),len(alt))!=1 or not 2<=n<=30 or set(ref+alt)-set('ACGT'):continue
        lo=v['pos']-300;hi=v['pos']+len(ref)+300;k=bisect.bisect_right(starts,lo)-1
        if k<0 or not intervals[k][0]<=lo<hi<=intervals[k][1]:continue
        if bisect.bisect_left(positions,hi)-bisect.bisect_left(positions,lo)!=1:continue
        group=('ins' if len(alt)>len(ref) else 'del')+'_'+('2-5' if n<=5 else '6-10' if n<=10 else '11-30')
        pool.append({'id':f"chr20_{v['pos']+1}",'indel':v,'stratum':group,'rank':hashlib.sha256(f"MP28_LOCUS_SPLIT_V1|{v['pos']}".encode()).hexdigest()})
    pool.sort(key=lambda c:c['rank']);n=P['n_calibration'];nt=P['n_test'];assert len(pool)>=n+nt
    chosen={'calibration':pool[:n],'test':pool[n:n+nt],'eligible_count':len(pool),'excluded_positions':sorted(excluded),'selection_before_outcomes':True}
    save(OUT/'SELECTION.json',chosen);return chosen
def validate_record(r,template):
    seq=template['sequence'];ins={};changes={};deletions=set()
    for e in r['edits']:
        i=e[1];assert r['source_start']<=i<r['source_end']
        if e[0]=='I':assert i not in ins;ins[i]=e[2]
        elif e[0]=='D':assert seq[i]==e[2];deletions.add(i)
        else:assert e[0]=='X' and seq[i]==e[2] and e[2]!=e[3];changes[i]=e[3]
    actual=''.join(ins.get(i,'')+('' if i in deletions else changes.get(i,seq[i])) for i in range(r['source_start'],r['source_end']))
    assert actual==r['sequence']
def maximum(values):return None if any(x is None for x in values) else max(values)
def run_site(stage,case):
    dest=OUT/stage/case['id'];done=dest/'SUMMARY.json'
    if done.exists():
        s=json.loads(done.read_text());assert s['protocol_sha256']==PH and s['engine_sha256']==EH and sha(dest/'RECORDS.jsonl.gz')==s['records_sha256'];return s
    templates=source_templates(REFERENCE,case);catalogues={'canonical':[t for t in templates if t['background']=='canonical'],'expanded':templates}
    regimes=['baseline','nuisance'] if stage=='calibration' else ['baseline','shifted','nuisance'];dest.mkdir(parents=True,exist_ok=True);temp=dest/'RECORDS.jsonl.gz.tmp';truthscores=collections.defaultdict(list);failures=[];n=0
    with gzip.open(temp,'wt') as f:
        for regime in regimes:
            for allele in ('0','1'):
                template=next(t for t in templates if t['allele']==allele and t['background']==('nuisance' if regime=='nuisance' else 'canonical'))
                for geometry in P['geometries']:
                    for rep in range(P['replicates']):
                        r=generate(template,case,regime,geometry,rep,P);validate_record(r,template)
                        scores={name:allele_scores(r['sequence'],cat,P['score_cap']) for name,cat in catalogues.items()};r['scores']=scores
                        for a in ('0','1'):
                            c=scores['canonical'][a];e=scores['expanded'][a]
                            assert c is None or (e is not None and e<=c)
                        truecat='expanded' if regime=='nuisance' else 'canonical';cost=scores[truecat][allele]
                        if len(r['edits'])<=P['score_cap'] and (cost is None or cost>len(r['edits'])):failures.append({'id':r['id'],'kind':'TRUE_HISTORY_EXCLUDED'})
                        for cat in catalogues:truthscores[(regime,cat)].append(scores[cat][allele])
                        f.write(json.dumps(r,separators=(',',':'))+'\n');n+=1
    temp.replace(dest/'RECORDS.jsonl.gz')
    s={'protocol_sha256':PH,'engine_sha256':EH,'stage':stage,'case':case,'records':n,'records_sha256':sha(dest/'RECORDS.jsonl.gz'),'soundness_failures':failures,'block_scores':{'|'.join(k):maximum(v) for k,v in truthscores.items()}}
    save(done,s);print('LOCUS_COMPLETE',stage,case['id'],n,flush=True);return s
def manifest():
    save(OUT/'MANIFEST.json',{'protocol_sha256':PH,'engine_sha256':EH,'files':[{'path':p.relative_to(OUT).as_posix(),'sha256':sha(p),'MB':p.stat().st_size/1e6} for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='MANIFEST.json' and not p.name.endswith('.tmp')]})
def stage_run(stage,cases,started):
    results=[];workers=min(4,os.cpu_count() or 1)
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('fork')) as pool:
        iterator=iter(cases);active={}
        for case in list(cases[:workers]):active[pool.submit(run_site,stage,case)]=case
        iterator=iter(cases[workers:])
        while active:
            done,_=concurrent.futures.wait(active,return_when=concurrent.futures.FIRST_COMPLETED)
            for future in done:
                case=active.pop(future);results.append(future.result())
                save(OUT/'CHECKPOINT.json',{'protocol_sha256':PH,'engine_sha256':EH,'stage':stage,'completed':len(results),'planned':len(cases),'last':case['id'],'seconds':time.monotonic()-started})
                if time.monotonic()-started<9.5*3600:
                    case=next(iterator,None)
                    if case:active[pool.submit(run_site,stage,case)]=case
    return sorted(results,key=lambda s:s['case']['id'])
def summarize(tests,freeze):
    rows=[];strata=collections.defaultdict(collections.Counter)
    for name,regime,cat,calibrator,status in P['comparisons']:
        for an,ad in P['alphas']:
            q=freeze['quantiles'][f'{calibrator}|{an}/{ad}'];counts=collections.Counter();blocks=0
            for site in tests:
                covered=True;seen=0
                with gzip.open(OUT/'test'/site['case']['id']/'RECORDS.jsonl.gz','rt') as f:
                    for line in f:
                        r=json.loads(line)
                        if r['regime']!=regime:continue
                        labels=prediction_set(r['scores'][cat],q);inside=r['allele'] in labels;covered=covered and inside;seen+=1
                        outcome='empty' if not labels else 'ambiguous' if len(labels)==2 else 'correct_singleton' if inside else 'wrong_singleton'
                        counts['records']+=1;counts['included']+=inside;counts[outcome]+=1
                        for group in ('geometry:'+r['geometry'],'stratum:'+site['case']['stratum']):
                            c=strata[(name,f'{an}/{ad}',group)];c['records']+=1;c['included']+=inside;c[outcome]+=1
                assert seen==12;blocks+=covered
            rows.append({'comparison':name,'alpha':[an,ad],'scope':status,'threshold':q,'blocks':len(tests),'covered_blocks':blocks,**dict(counts)})
    return rows,[{'comparison':k[0],'alpha':k[1],'group':k[2],**dict(v)} for k,v in sorted(strata.items())]
def main():
    global REFERENCE
    started=time.monotonic();OUT.mkdir(parents=True,exist_ok=True)
    try:
        for p in Path('/kaggle/input').rglob('CHECKPOINT.json'):
            checkpoint=json.loads(p.read_text())
            if checkpoint.get('protocol_sha256')==PH and checkpoint.get('engine_sha256')==EH:
                for row in json.loads((p.parent/'MANIFEST.json').read_text())['files']:
                    src=verified(p.parent,row['path']);dst=OUT/row['path'];dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst)
        roots=sources();chosen=selection(roots);ref=verified(roots['MP25'],'inputs/chr20.fa.gz');assert sha(ref)==B['reference']['sha256']
        with gzip.open(ref,'rt') as f:assert next(f).strip()=='>chr20';REFERENCE=''.join(l.strip() for l in f).upper()
        assert len(REFERENCE)==64444167
        save(OUT/'PROTOCOL.json',P);save(OUT/'INPUT_BINDING.json',B);save(OUT/'ENVIRONMENT.json',{'python':sys.version,'platform':platform.platform(),'numpy':__import__('numpy').__version__,'numba':__import__('numba').__version__,'pysam':__import__('pysam').__version__,'engine_sha256':EH})
        fragment_distance('ACGT','ACGT',32)
        calibration=stage_run('calibration',chosen['calibration'],started)
        if len(calibration)!=P['n_calibration']:save(OUT/'RESULT.json',{'status':'CHECKPOINT_CONTINUABLE','stage':'calibration'});return
        quantiles={}
        for key,scorekey in P['calibrators'].items():
            values=[s['block_scores'][scorekey.replace('/','|')] for s in calibration]
            for an,ad in P['alphas']:quantiles[f'{key}|{an}/{ad}']=block_quantile(values,an,ad)
        freeze={'protocol_sha256':PH,'engine_sha256':EH,'quantiles':quantiles,'calibration_summary_hashes':{s['case']['id']:sha(OUT/'calibration'/s['case']['id']/'SUMMARY.json') for s in calibration},'before_test_scoring':True}
        frozen=OUT/'CALIBRATION_FREEZE.json'
        if frozen.exists():assert json.loads(frozen.read_text())==freeze
        else:save(frozen,freeze)
        print('CALIBRATION_FROZEN',json.dumps(quantiles),flush=True)
        tests=stage_run('test',chosen['test'],started);rows,groups=summarize(tests,freeze);failures=sum(len(s['soundness_failures']) for s in calibration+tests)
        result={'status':'COMPLETE' if len(tests)==P['n_test'] else 'CHECKPOINT_CONTINUABLE','protocol_sha256':PH,'engine_sha256':EH,'calibration_blocks':len(calibration),'test_blocks':len(tests),'calibration_records':sum(s['records'] for s in calibration),'test_records':sum(s['records'] for s in tests),'soundness_failures':failures,'seconds':time.monotonic()-started,'comparisons':rows,'stratified':groups,'scope':P['coverage_boundary']}
        save(OUT/'RESULT.json',result);print('CAMPAIGN_RESULT',json.dumps({k:v for k,v in result.items() if k!='stratified'}),flush=True)
    except Exception as e:
        save(OUT/'FAILURE.json',{'type':type(e).__name__,'message':str(e),'traceback':traceback.format_exc()});raise
    finally:manifest()

if __name__=='__main__':main()
