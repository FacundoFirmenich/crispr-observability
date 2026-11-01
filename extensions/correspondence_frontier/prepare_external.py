"""Prepare a fresh, separate twelve-job HEK3 engine execution from public real reads."""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess
HERE=Path(__file__).resolve().parent
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work',type=Path,required=True);p.add_argument('--engine',type=Path,required=True);a=p.parse_args();w=a.work.resolve();engine=a.engine.resolve();assert not w.exists() and not w.is_relative_to(HERE)
 assert subprocess.check_output([str(engine),'--version'],text=True).strip()=='CRISPResso 2.3.4'
 w.mkdir(parents=True)
 for name in ['RULE_FREEZE.json','ACQUISITION.json'] :shutil.copy2(HERE/'inputs/hek3'/name,w/name)
 shutil.copytree(HERE/'inputs/hek3/inputs',w/'inputs')
 data=json.loads((HERE/'provenance/HEK3_COMMANDS.json').read_text());jobs=data['jobs']
 for j in jobs:j['command']=[str(engine) if i==0 else x.replace('${HEK3_WORK}',str(w)) for i,x in enumerate(j['command'])]
 (w/'JOBS.json').write_text(json.dumps(jobs,indent=2)+'\n')
 sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
 (w/'ENGINE_FREEZE.json').write_text(json.dumps({'jobs_sha256':sha(w/'JOBS.json'),'rule_sha256':sha(w/'RULE_FREEZE.json'),'scope':'Fresh repetition of the already retrospective fixed-rule experiment; no new prospective status.'},indent=2)+'\n')
 print('Prepared',len(jobs),'jobs. Execute run_external.py --work',w)
if __name__=='__main__':main()
