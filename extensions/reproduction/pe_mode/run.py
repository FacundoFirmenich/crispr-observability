"""Execute only the commands frozen by prepare.py; retain every execution record."""
from pathlib import Path
import argparse,datetime,hashlib,json,subprocess,time
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work',type=Path,required=True);a=p.parse_args();w=a.work.resolve();f=json.loads((w/'FREEZE.json').read_text());assert sha(w/'JOBS.json')==f['jobs_sha256'] and sha(w/'INPUT_MANIFEST.json')==f['inputs_sha256'];inputs={x['run']:x for x in json.loads((w/'INPUT_MANIFEST.json').read_text())};(w/'executions').mkdir(exist_ok=True)
 for job in json.loads((w/'JOBS.json').read_text()):
  rid=job['run'];assert sha(w/'inputs'/(rid+'.Q30.fastq.gz'))==inputs[rid]['filtered_sha256'];p=w/'executions'/(rid+'.json');assert not p.exists(),'Preserve prior attempts and use a new work directory'
  record={**job,'state':'RUNNING','started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()};p.write_text(json.dumps(record,indent=2)+'\n');t=time.monotonic()
  with p.with_suffix('.stdout').open('w') as o,p.with_suffix('.stderr').open('w') as e:r=subprocess.run(job['command'],stdout=o,stderr=e)
  record.update(returncode=r.returncode,state='COMPLETED' if r.returncode==0 else 'FAILED',seconds=time.monotonic()-t);p.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({'run':rid,'state':record['state']}),flush=True);assert r.returncode==0
if __name__=='__main__':main()
