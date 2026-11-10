"""MIT. Verify MP28 source, read every result, recount and preserve exact bytes."""
import collections,concurrent.futures,gzip,hashlib,io,json,math,sys,time,zipfile
from pathlib import Path
HERE=Path(__file__).resolve().parent;R=HERE.parents[1];sys.path.insert(0,str(HERE.parent))
import monitor as transport
from kaggle.api.kaggle_api_extended import KaggleApi
from kagglesdk.kernels.types.kernels_api_service import ApiGetKernelRequest,ApiListKernelSessionOutputRequest
transport.TAG='v0.2.0';KERNEL='fjfmad/crispr-mp28-locus-calibration-v1'
def sha(b):return hashlib.sha256(b).hexdigest()
def save(name,x):(R/'results'/name).write_text(json.dumps(x,indent=2)+'\n',encoding='utf8')
def load(data,name):return json.loads(data['MP28/'+name])
def labels(scores,q):return sorted(scores) if q['infinity'] else sorted(a for a,c in scores.items() if c is not None and c<=q['q'])
def templates(reference,case):
    v=case['indel'];pos=v['pos'];out={}
    assert reference[pos:pos+len(v['ref'])]==v['ref']
    for allele,event in enumerate([v['ref'],v['alts'][0]]):
        s=reference[pos-300:pos]+event+reference[pos+len(v['ref']):pos+len(v['ref'])+300];end=300+len(event);b=list(s)
        for i in [275,end+25]:b[i]={'A':'C','C':'G','G':'T','T':'A'}[b[i]]
        b[end+40:end+40]=list('GCA')
        out[f'canonical_{allele}']=s;out[f'nuisance_{allele}']=''.join(b)
    return out
def independent_distance(h,r):
    row=list(range(len(r)+1));best=len(r)+len(h)
    for a in h:
        nxt=[0]
        for j,b in enumerate(r):nxt.append(min(nxt[-1]+1,row[j+1]+1,row[j]+int(a!=b)))
        best=min(best,nxt[-1]);row=nxt
    return best
def verify(data):
    p=load(data,'PROTOCOL.json');selection=load(data,'SELECTION.json');freeze=load(data,'CALIBRATION_FREEZE.json');result=load(data,'RESULT.json')
    assert p==json.loads((R/'protocols/MP28.json').read_text())
    refbytes=transport.download('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz',21_000_000)
    assert sha(refbytes)=='271ca0ee4247c2fa30c8491764ad8acd4bdde8eb3cf314f0720235010ae85c9f'
    reference=''.join(gzip.decompress(refbytes).decode().splitlines()[1:]).upper();cases={stage:{c['id']:c for c in selection[stage]} for stage in ('calibration','test')}
    assert not set(cases['calibration'])&set(cases['test'])
    assert all(c['indel']['pos'] not in selection['excluded_positions'] for d in cases.values() for c in d.values())
    allrecords={};block_scores={};replayed=0;scorechecked=0;seen=set();failures=0
    for name,b in sorted(data.items()):
        if not name.endswith('/SUMMARY.json'):continue
        s=json.loads(b);stage=s['stage'];site=s['case']['id'];assert s['case']==cases[stage][site]
        recordname=name.rsplit('/',1)[0]+'/RECORDS.jsonl.gz';assert sha(data[recordname])==s['records_sha256']
        records=[json.loads(line) for line in gzip.decompress(data[recordname]).splitlines()];assert len(records)==s['records'];source=templates(reference,s['case']);scores=collections.defaultdict(list)
        for record in records:
            assert record['id'] not in seen;seen.add(record['id']);sequence=source[record['template_id']];edits=collections.defaultdict(list)
            for e in record['edits']:
                assert record['source_start']<=e[1]<record['source_end'];edits[e[1]].append(e)
            replay=[]
            for i in range(record['source_start'],record['source_end']):
                emitted=sequence[i]
                for e in edits[i]:
                    if e[0]=='I':replay.append(e[2])
                    elif e[0]=='D':assert emitted==e[2];emitted=''
                    else:assert e[0]=='X' and emitted==e[2] and e[2]!=e[3];emitted=e[3]
                replay.append(emitted)
            assert ''.join(replay)==record['sequence'];replayed+=1
            for cat in ('canonical','expanded'):
                scores[record['regime']+'|'+cat].append(record['scores'][cat][record['allele']])
            for a in ('0','1'):
                c=record['scores']['canonical'][a];e=record['scores']['expanded'][a];assert c is None or (e is not None and e<=c)
            truecat='expanded' if record['regime']=='nuisance' else 'canonical';truthscore=record['scores'][truecat][record['allele']]
            if len(record['edits'])<=p['score_cap'] and (truthscore is None or truthscore>len(record['edits'])):failures+=1
            # Predeclared deterministic audit subset: first 36 test records in path order.
            if stage=='test' and scorechecked<36:
                exact={key:independent_distance(seq,record['sequence']) for key,seq in source.items()}
                for a in ('0','1'):
                    for cat in ('canonical','expanded'):
                        cost=exact['canonical_'+a] if cat=='canonical' else min(exact['canonical_'+a],exact['nuisance_'+a])
                        assert record['scores'][cat][a]==(cost if cost<=p['score_cap'] else None)
                scorechecked+=1
        expected={k:None if None in v else max(v) for k,v in scores.items()};assert expected==s['block_scores']
        if stage=='calibration':assert sha(b)==freeze['calibration_summary_hashes'][site]
        block_scores[(stage,site)]=expected;allrecords[(stage,site)]=records
    assert replayed==result['calibration_records']+result['test_records'];assert failures==result['soundness_failures']
    for key,scorekey in p['calibrators'].items():
        values=sorted(math.inf if v[scorekey.replace('/','|')] is None else v[scorekey.replace('/','|')] for (stage,site),v in block_scores.items() if stage=='calibration')
        for an,ad in p['alphas']:
            rank=((len(values)+1)*(ad-an)+ad-1)//ad;q=math.inf if rank>len(values) else values[rank-1];actual=freeze['quantiles'][f'{key}|{an}/{ad}']
            assert actual['rank']==rank and actual['q']==(None if math.isinf(q) else q) and actual['infinity']==math.isinf(q)
    for row in result['comparisons']:
        spec=next(x for x in p['comparisons'] if x[0]==row['comparison']);name,regime,cat,calibrator,_=spec;an,ad=row['alpha'];q=freeze['quantiles'][f'{calibrator}|{an}/{ad}'];counts=collections.Counter();covered=blocks=0
        for (stage,site),records in allrecords.items():
            if stage!='test':continue
            block=True;n=0
            for record in records:
                if record['regime']!=regime:continue
                call=labels(record['scores'][cat],q);inside=record['allele'] in call;n+=1;block=block and inside
                outcome='empty' if not call else 'ambiguous' if len(call)==2 else 'correct_singleton' if inside else 'wrong_singleton'
                counts['records']+=1;counts['included']+=inside;counts[outcome]+=1
                if regime=='nuisance':assert set(labels(record['scores']['canonical'],q))<=set(labels(record['scores']['expanded'],q))
            assert n==12;blocks+=1;covered+=block
        assert blocks==row['blocks'] and covered==row['covered_blocks']
        for key in ('records','included','empty','ambiguous','correct_singleton','wrong_singleton'):assert counts[key]==row.get(key,0),(row['comparison'],key)
    return {'status':'PASS_FULL_RECORD_REPLAY_AND_RECOUNT','records_replayed':replayed,'independent_score_check_records':scorechecked,'quantile_recount':'PASS','catalogue_monotonicity':'PASS','soundness_failures':failures,'scope':'All generated strings and primary counts replayed; 36 records independently rescored. Controlled simulation, not biological calibration.'}
def observe():
    k=KaggleApi();k.authenticate();status=str(k.kernels_status(KERNEL).status)
    q=ApiGetKernelRequest();q.user_name,q.kernel_slug=KERNEL.split('/');q.version_label=k._kernel_version_label('1')
    with k.build_kaggle_client() as c:remote=c.kernels.kernels_api_client.get_kernel(q)
    expected=json.loads((HERE/'kaggle/SOURCE_FREEZE.json').read_text())['kernel_source_sha256'];assert remote.metadata.is_private and sha(remote.blob.source.encode())==expected
    launch={'kernel':KERNEL,'version':1,'private':True,'status':status,'full_source_readback':'PASS','source_sha256':expected};save('MP28_LAUNCH_VERIFIED.json',launch)
    if any(x in status.upper() for x in ('RUNNING','QUEUED')):print(json.dumps(launch));return False
    files={};token=None
    while True:
        q=ApiListKernelSessionOutputRequest();q.user_name,q.kernel_slug=KERNEL.split('/');q.version_label=k._kernel_version_label('1');k._set_paging(q,200,token)
        with k.build_kaggle_client() as c:page=c.kernels.kernels_api_client.list_kernel_session_output(q)
        for f in page.files or []:files[f.file_name]=f.url
        token=page.next_page_token
        if not token:break
    if 'MP28/MANIFEST.json' not in files:raise RuntimeError('Terminal manifest absent: '+status)
    raw=transport.download(files['MP28/MANIFEST.json'],5_000_000);manifest=json.loads(raw);data={'MP28/MANIFEST.json':raw};assert sum(r['MB'] for r in manifest['files'])<200
    def read(row):
        name='MP28/'+row['path'];assert '..' not in Path(row['path']).parts
        b=transport.download(files[name],round(row['MB']*1e6)+1);assert sha(b)==row['sha256'];return name,b
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        for name,b in pool.map(read,manifest['files']):data[name]=b
    failure=load(data,'FAILURE.json') if 'MP28/FAILURE.json' in data else None
    verification={'status':'EXECUTION_FAILURE','failure':failure} if failure else verify(data)
    buf=io.BytesIO()
    with zipfile.ZipFile(buf,'w',zipfile.ZIP_DEFLATED) as z:
        for name,b in sorted(data.items()):
            info=zipfile.ZipInfo(name,date_time=(2026,11,10,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;z.writestr(info,b)
    custody=transport.custody('MP28_OUTPUTS.zip',buf.getvalue());save('MP28_CUSTODY.json',custody);save('MP28_VERIFICATION.json',verification)
    if 'MP28/RESULT.json' in data:save('MP28_RESULT.json',load(data,'RESULT.json'))
    save('MP28_CALIBRATION_FREEZE.json',load(data,'CALIBRATION_FREEZE.json')) if 'MP28/CALIBRATION_FREEZE.json' in data else None
    print(json.dumps({'launch':launch,'verification':verification,'custody':custody}));return True
if __name__=='__main__':observe()
