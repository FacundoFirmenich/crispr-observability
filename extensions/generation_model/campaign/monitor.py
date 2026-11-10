"""Bounded read-only Kaggle observation and exact result custody in this repo.
No scientific reruns, cancellations or changes to other repositories.
"""
import collections,concurrent.futures,gzip,hashlib,io,json,subprocess,sys,time,urllib.request,urllib.parse,zipfile
from datetime import datetime,timezone
from pathlib import Path
from kaggle.api.kaggle_api_extended import KaggleApi
from kagglesdk.kernels.types.kernels_api_service import ApiGetKernelRequest,ApiListKernelSessionOutputRequest
R=Path(__file__).resolve().parents[1];OUT=R/'results';REPO='FacundoFirmenich/crispr-observability-extended';TAG='v0.1.0';KERNEL='fjfmad/crispr-extended-mp27-generation-v1'
def gh(*args):return subprocess.check_output(['gh',*args],text=True,encoding='utf8')
def sha(b):return hashlib.sha256(b).hexdigest()
def save(name,x):(OUT/name).write_text(json.dumps(x,indent=2)+'\n',encoding='utf8')
def download(url,limit=250_000_000):
    with urllib.request.urlopen(url,timeout=180) as f:b=f.read(limit+1)
    assert len(b)<=limit;return b
class Safe(urllib.request.HTTPRedirectHandler):
    def redirect_request(self,req,*args,**kw):
        r=super().redirect_request(req,*args,**kw)
        if r:r.remove_header('Authorization')
        return r
def asset_bytes(repo,asset_id):
    req=urllib.request.Request(f'https://api.github.com/repos/{repo}/releases/assets/{asset_id}',headers={'Authorization':'Bearer '+gh('auth','token').strip(),'Accept':'application/octet-stream'})
    with urllib.request.build_opener(Safe).open(req,timeout=180) as f:b=f.read(250_000_001)
    assert len(b)<250_000_001;return b
def custody(name,payload):
    assert json.loads(gh('api','repos/'+REPO))['private']
    release=json.loads(gh('api',f'repos/{REPO}/releases/tags/{TAG}'));asset=next((a for a in release['assets'] if a['name']==name),None)
    if asset is None:
        url=release['upload_url'].split('{')[0]+'?'+urllib.parse.urlencode({'name':name})
        req=urllib.request.Request(url,data=payload,method='POST',headers={'Authorization':'Bearer '+gh('auth','token').strip(),'Content-Type':'application/zip'})
        with urllib.request.urlopen(req,timeout=180) as f:asset=json.load(f)
    readback=asset_bytes(REPO,asset['id']);assert sha(readback)==sha(payload) and len(readback)==len(payload)
    return {'status':'REMOTE_FULL_READBACK_PASS','repository':REPO,'tag':TAG,'asset':name,'asset_id':asset['id'],'sha256':sha(payload),'MB':len(payload)/1e6}
def prior_records():
    specs=[('MP24',732743229,'79aa7b1510d814fe236b5ad1f2ae42bc5f027192c4dac5f4d29a467c489078b7'),('MP25',732724176,'e87a2863f8a9a84ffc3b297558e2157de0ad6dbd93215819c81f16375517d80c')];out={}
    for prefix,aid,expected in specs:
        b=asset_bytes('FacundoFirmenich/crispr-p1-private-custody',aid);assert sha(b)==expected
        custody(prefix+'_INPUT_EVIDENCE.zip',b)
        z=zipfile.ZipFile(io.BytesIO(b));rows={}
        for name in z.namelist():
            if name.startswith(prefix+'/sites/') and name.endswith(('/RECORDS.jsonl.gz','/READS.jsonl.gz')):
                for line in gzip.decompress(z.read(name)).splitlines():
                    r=json.loads(line);identity=r.get('id') or r['read_id']+':'+str(r['flag']&192);key=(r['site'],identity);assert key not in rows
                    rows[key]={'model':r.get('model','real'),'geometry':r.get('geometry','real'),'allele':str(r['allele']) if 'allele' in r else {'E':'1','NE':'0','UNKNOWN':'UNKNOWN'}[r['truth_alt']],'errors':len(r['edits']) if 'edits' in r else None,'truth':r.get('truth')}
        out['development' if prefix=='MP24' else 'real']=rows
    return out
def observe(k):
    status=str(k.kernels_status(KERNEL).status);snap={'utc':datetime.now(timezone.utc).isoformat(),'kernel':KERNEL,'version':1,'status':status};save('MONITOR_STATUS.json',snap)
    if any(s in status.upper() for s in ('RUNNING','QUEUED')):return False
    q=ApiGetKernelRequest();q.user_name,q.kernel_slug=KERNEL.split('/');q.version_label=k._kernel_version_label('1')
    with k.build_kaggle_client() as c:r=c.kernels.kernels_api_client.get_kernel(q)
    expected=json.loads((R/'campaign/kaggle/SOURCE_FREEZE.json').read_text())['kernel_source_sha256'];assert r.metadata.is_private and sha(r.blob.source.encode())==expected
    files={};token=None
    while True:
        q=ApiListKernelSessionOutputRequest();q.user_name,q.kernel_slug=KERNEL.split('/');q.version_label=k._kernel_version_label('1');k._set_paging(q,200,token)
        with k.build_kaggle_client() as c:r=c.kernels.kernels_api_client.list_kernel_session_output(q)
        for f in r.files or []:assert f.file_name not in files;files[f.file_name]=f.url
        token=r.next_page_token
        if not token:break
    name='EXTENDED_MP27/MANIFEST.json';assert name in files,'No terminal manifest available'
    raw=download(files[name],5_000_000);m=json.loads(raw);data={name:raw}
    assert sum(row['MB'] for row in m['files'])<200
    def read_output(row):
        name='EXTENDED_MP27/'+row['path'];assert '..' not in Path(row['path']).parts and name in files
        b=download(files[name],int(round(row['MB']*1e6))+1);assert sha(b)==row['sha256'];return name,b
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        for name,b in pool.map(read_output,m['files']):data[name]=b
    failure=json.loads(data.get('EXTENDED_MP27/FAILURE.json',b'{}'));result=json.loads(data.get('EXTENDED_MP27/RESULT.json',b'{}'))
    if failure:verification={'status':'EXECUTION_FAILURE','failure':failure,'execution':result}
    else:
        prior=prior_records();total=collections.Counter();recounted=collections.Counter();failures=0;seen=set();allrows=0
        refbytes=download('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz',21_000_000)
        assert sha(refbytes)=='271ca0ee4247c2fa30c8491764ad8acd4bdde8eb3cf314f0720235010ae85c9f'
        reference=''.join(gzip.decompress(refbytes).decode().splitlines()[1:]).upper()
        for name,b in data.items():
            if not name.endswith('/SUMMARY.json'):continue
            s=json.loads(b);stage=s['stage'];site=s['case']['id'];actual=collections.Counter();records_name=name.rsplit('/',1)[0]+'/RECORDS.jsonl.gz';assert sha(data[records_name])==s['records_sha256'];n=0
            for line in gzip.decompress(data[records_name]).splitlines():
                r=json.loads(line);key=(stage,site,r['id']);assert key not in seen;seen.add(key);n+=1
                if stage=='heldout':
                    g=r['generated_observation'];source={'model':g['model'],'geometry':g['geometry'],'allele':str(g['allele']),'errors':len(g['edits'])}
                    v=s['case']['indel'];pos=v['pos'];variant=v['alts'][0] if g['allele'] else v['ref'];template=reference[pos-300:pos]+variant+reference[pos+len(v['ref']):pos+len(v['ref'])+300]
                    ins=collections.defaultdict(list);deleted=set();changes={}
                    for e in g['edits']:
                        assert g['source_start']<=e[1]<g['source_end']
                        if e[0]=='I':ins[e[1]].append(e[2])
                        elif e[0]=='D':assert template[e[1]]==e[2] and e[1] not in deleted;deleted.add(e[1])
                        else:assert e[0]=='X' and template[e[1]]==e[2] and e[1] not in changes;changes[e[1]]=e[3]
                    bases=[];origins=[];inserted=[]
                    for i in range(g['source_start'],g['source_end']):
                        for b in ins[i]:bases.append(b);origins.append(i);inserted.append(True)
                        if i not in deleted:bases.append(changes.get(i,template[i]));origins.append(i);inserted.append(False)
                    assert ''.join(bases)==g['sequence'] and origins==g['source_origins'] and inserted==g['inserted']
                else:source=prior[stage][site,r['id']]
                for config in r['results']:
                    pad=config['padding'];truth=config['observed_truth'];allele=config['source_allele_truth'];errors=config['literal_error_count'];assert allele==source['allele'] and errors==source['errors']
                    if stage=='development':assert truth==source['truth'][str(pad)]['observed_label']
                    if stage=='real':assert truth=='UNAVAILABLE'
                    if stage=='heldout':
                        lo=300-pad;hi=300+len(variant)+pad;survivors={i for i,ins in zip(origins,inserted) if not ins}
                        span=''.join(b for b,i in zip(bases,origins) if lo<=i<hi)
                        target=reference[pos-pad:pos]+v['alts'][0]+reference[pos+len(v['ref']):pos+len(v['ref'])+pad]
                        expected='U' if lo-1 not in survivors or hi not in survivors or not span else 'E' if span==target else 'NE';assert truth==expected
                    previous=set()
                    for d in config['decisions']:
                        budget=d['budget'];labels={lab for domain in config['domains'] for lab,c in domain['minimum_errors'].items() if c is not None and c<=budget};aa={str(domain['allele']) for domain in config['domains'] if any(c is not None and c<=budget for c in domain['minimum_errors'].values())}
                        assert labels==set(d['possible_observed_labels']) and aa==set(d['possible_source_alleles']) and previous<=labels;previous=labels
                        call=next(iter(labels)) if len(labels)==1 and 'U' not in labels else 'U';ac=next(iter(aa)) if len(aa)==1 else 'U';assert call==d['observed_tract_certificate'] and ac==(d['source_allele_certificate'] or 'U')
                        inside=stage!='real' and errors<=budget
                        if inside and (truth not in labels or allele not in aa):failures+=1
                        kk=(pad,budget,truth,call,allele,ac,str(inside),source['model'],source['geometry']);actual[kk]+=1;total[(stage,)+kk]+=1;allrows+=1
            assert actual==collections.Counter({tuple(row['key']):row['n'] for row in s['counts']}) and n==s['records'];recounted[stage]+=n
        assert total==collections.Counter({tuple(row['key']):row['n'] for row in result['counts']}) and dict(recounted)==result['records_by_stage'] and failures==result['soundness_failures']
        verification={'status':'PASS_INTERNAL_RECOUNT','execution':{k:v for k,v in result.items() if k!='counts'},'records_recounted':dict(recounted),'decision_rows_recounted':allrows,'scope':'Recount from frozen class minima, not a second full DP or independent biological validation'}
        save('MP27_RESULT.json',result)
    save('MP27_VERIFICATION.json',verification);data['VERIFICATION.json']=json.dumps(verification,indent=2).encode()
    buf=io.BytesIO()
    with zipfile.ZipFile(buf,'w',zipfile.ZIP_DEFLATED) as z:
        for name,b in sorted(data.items()):
            zi=zipfile.ZipInfo(name,date_time=(2026,11,10,0,0,0));zi.compress_type=zipfile.ZIP_DEFLATED;z.writestr(zi,b)
    receipt=custody('EXTENDED_MP27_OUTPUTS.zip',buf.getvalue());save('MP27_CUSTODY.json',receipt);snap.update(verification=verification['status'],custody=receipt);save('MONITOR_STATUS.json',snap);print(json.dumps(snap),flush=True);return True
if __name__=='__main__':
    k=KaggleApi();k.authenticate();deadline=time.monotonic()+12*3600
    while True:
        try:
            if observe(k):break
        except Exception as e:save('MONITOR_ERROR.json',{'utc':datetime.now(timezone.utc).isoformat(),'type':type(e).__name__,'message':str(e)})
        if len(sys.argv)<2 or sys.argv[1]!='watch' or time.monotonic()>deadline:break
        time.sleep(120)
