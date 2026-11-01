"""Execute the frozen source-authenticated HEK3 external panel, sequentially."""
from pathlib import Path
import argparse,json,gzip,hashlib,subprocess,time,os
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--work',type=Path,required=True);a=p.parse_args();b=a.work.resolve();freeze=json.loads((b/'ENGINE_FREEZE.json').read_text());assert sha(b/'JOBS.json')==freeze['jobs_sha256'] and sha(b/'RULE_FREEZE.json')==freeze['rule_sha256'];inputs=[]
 for r in json.loads((b/'ACQUISITION.json').read_text()):
  rid=r['run'];src=b/'inputs'/(rid+'.fastq.gz');assert sha(src)==r['sha256'];dst=b/'filtered'/(rid+'.Q30.fastq.gz');dst.parent.mkdir(exist_ok=True);n=qpass=0;logical=hashlib.sha256()
  if not dst.exists():
   with gzip.open(src,'rt') as f,dst.open('wb') as out:
    with gzip.GzipFile(filename='',mode='wb',fileobj=out,mtime=0) as g:
     while h:=f.readline():
      s,plus,q=f.readline().strip(),f.readline(),f.readline().strip();assert h.startswith('@') and plus.startswith('+') and len(s)==len(q);n+=1;logical.update((s+'\n'+q+'\n').encode())
      if sum(ord(x)-33 for x in q)>=30*len(q):qpass+=1;g.write((h.rstrip()+'\n'+s+'\n+\n'+q+'\n').encode())
    out.flush();os.fsync(out.fileno())
  else:
   old=next(x for x in json.loads((b/'FILTERED_INPUTS.json').read_text()) if x['run']==rid);n,qpass=old['N'],old['Q30']
  with gzip.open(dst,'rt') as f:assert sum(1 for _ in f)==4*qpass
  inputs.append({'run':rid,'N':n,'Q30':qpass,'raw_sha256':r['sha256'],'filtered_sha256':sha(dst)})
 (b/'FILTERED_INPUTS.json').write_text(json.dumps(inputs,indent=2)+'\n')
 for j in json.loads((b/'JOBS.json').read_text()):
  p=b/'executions'/j['mode']/(j['run']+'.json');p.parent.mkdir(parents=True,exist_ok=True)
  if p.exists():assert json.loads(p.read_text())['state']=='COMPLETED';continue
  inp=next(x for x in inputs if x['run']==j['run']);assert sha(b/'filtered'/(j['run']+'.Q30.fastq.gz'))==inp['filtered_sha256'];t=time.monotonic()
  with p.with_suffix('.stdout').open('w') as o,p.with_suffix('.stderr').open('w') as e:r=subprocess.run(j['command'],stdout=o,stderr=e)
  obj={**j,'state':'FAILED','returncode':r.returncode,'seconds':time.monotonic()-t};p.write_text(json.dumps(obj,indent=2)+'\n');assert r.returncode==0
  folder=b/'runs'/j['mode']/('CRISPResso_on_'+j['run'])
  with gzip.open(folder/'CRISPResso_output.fastq.gz','rt') as f:assert sum(1 for _ in f)==4*inp['Q30']
  obj.update(state='COMPLETED',annotation_EOF_and_Q30_count_verified=True);p.write_text(json.dumps(obj,indent=2)+'\n');print({k:obj[k] for k in ['run','mode','state','seconds']},flush=True)
 (b/'RUN_STATUS.json').write_text(json.dumps({'state':'COMPLETED','jobs':12,'source_study':'Anzalone2019 PRJNA565979 HEK3 Fig2b_10','raw_records':sum(x['N'] for x in inputs)},indent=2)+'\n')
if __name__=='__main__':main()
