"""Finish a whole engine job in staging; verify/fsync before exposing final outputs."""
import argparse, gzip, hashlib, json, os, shutil, subprocess, time, zipfile
from concurrent.futures import ThreadPoolExecutor
import run_grid as g
import run_profiles as profiles
from checked_run import verify

def execute(job,spec):
    manifest=next(x for x in json.loads((g.BASE/'INPUT_MANIFEST.json').read_text()) if x['run']==job['run'])
    inp=g.BASE/'inputs'/(job['run']+'.Q30.fastq.gz');assert g.sha(inp)==manifest['input_sha256']
    rec=g.BASE/'executions'/job['reference']/job['setting']/(job['run']+'.json')
    final=g.BASE/'runs'/job['reference']/job['setting']/('CRISPResso_on_'+job['run'])
    if rec.exists():
        saved=g.BASE/'failed_attempts/post_run_output_integrity'/job['reference']/job['setting']/job['run'];saved.mkdir(parents=True,exist_ok=False)
        for p in rec.parent.glob(job['run']+'.*'):shutil.copy2(p,saved/p.name)
        if final.exists():shutil.copytree(final,saved/final.name)
    live=__import__('pathlib').Path('/tmp/mp02_engine_jobs')/job['reference']/job['setting']/job['run'];assert not live.exists()
    cmd=profiles.command(job,spec);cmd[cmd.index('-o')+1]=str(live)
    obj={**job,'command':cmd,'final_output':str(final),'state':'RUNNING','started_utc':g.now(),'input_sha256_before':g.sha(inp)};g.dump(rec,obj);t=time.monotonic()
    with rec.with_suffix('.stdout').open('w') as o,rec.with_suffix('.stderr').open('w') as e:p=subprocess.run(cmd,stdout=o,stderr=e)
    obj.update(returncode=p.returncode,elapsed_seconds=time.monotonic()-t,finished_utc=g.now(),input_sha256_after=g.sha(inp))
    if p.returncode or obj['input_sha256_before']!=obj['input_sha256_after']:
        obj['state']='FAILED';g.dump(rec,obj);raise RuntimeError(job['id'])
    folder=live/('CRISPResso_on_'+job['run']);files={}
    expected={name:(s,q) for name,s,plus,q in g.rows(inp)};seen=set()
    for name,s,plus,q in g.rows(folder/'CRISPResso_output.fastq.gz'):
        assert name not in seen and (s,q)==expected[name];seen.add(name)
    assert len(seen)==manifest['Q30_pass']
    for f in folder.rglob('*'):
        if not f.is_file():continue
        if f.suffix=='.gz':gzip.open(f,'rb').read()
        if f.suffix=='.zip':
            with zipfile.ZipFile(f) as z:assert z.testzip() is None
        with f.open('rb') as h:os.fsync(h.fileno())
        files[str(f.relative_to(folder))]={'sha256':g.sha(f),'bytes':f.stat().st_size}
    assert g.sha(inp)==manifest['input_sha256']
    final.parent.mkdir(parents=True,exist_ok=True)
    ready=final.with_name(final.name+'.closed_staging')
    assert not ready.exists();shutil.copytree(folder,ready)
    for name,meta in files.items():
        assert g.sha(ready/name)==meta['sha256']
        with (ready/name).open('rb') as h:os.fsync(h.fileno())
    if final.exists():shutil.rmtree(final)
    os.replace(ready,final);fd=os.open(str(final.parent),os.O_DIRECTORY);os.fsync(fd);os.close(fd)
    for name,meta in files.items():assert g.sha(final/name)==meta['sha256']
    obj.update(state='COMPLETED',integrity='ALL_GZIP_ZIP_RAW_IDENTITIES_VERIFIED_AND_FSYNCED',files=files);g.dump(rec,obj)
    print(json.dumps({'id':job['id'],'state':obj['state'],'seconds':obj['elapsed_seconds']}),flush=True)
    return obj
def main(mode,workers):
    verify(mode+'_staged_before');spec=json.loads((g.BASE/'GRID.json').read_text());jobs=json.loads((g.BASE/'JOBS.json').read_text())
    for f in ['GRID_FREEZE.json','STOCK_PROFILE_FREEZE.json']:
        for p,h in json.loads((g.BASE/f).read_text())['files'].items():assert g.sha(g.BASE/p)==h
    if mode=='baseline_retry':
        failed={x['id'] for x in json.loads((g.BASE/'BASELINE_FILE_INTEGRITY_AUDIT.json').read_text()) if x['state']=='FAIL'};jobs=[j for j in jobs if j['id'] in failed]
    elif mode=='stock':jobs=[{**j,'setting':'stock_default_no_guide','id':j['id'].replace('/baseline/','/stock_default_no_guide/')} for j in jobs if j['setting']=='baseline']
    else:jobs=[j for j in jobs if j['setting']!='baseline']
    with ThreadPoolExecutor(max_workers=workers) as pool:outcomes=list(pool.map(lambda j:execute(j,spec),jobs))
    verify(mode+'_staged_after');g.dump(g.BASE/(mode.upper()+'_STATUS.json'),{'state':'COMPLETED','jobs':len(outcomes),'workers':workers,'utc':g.now(),'method':'staged output + gzip/zip/raw identity verification + fsync + final hash verification'})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['baseline_retry','stock','oat']);p.add_argument('--workers',type=int,default=3);a=p.parse_args();main(a.mode,a.workers)
