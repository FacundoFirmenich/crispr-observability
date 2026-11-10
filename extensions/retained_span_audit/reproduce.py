"""Single entry point: hash checks, graph validation and authentic-data replay."""
from pathlib import Path
import argparse, json, subprocess, sys, time, platform
from verify_release import validate

ROOT=Path(__file__).resolve().parent
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mode',choices=['quick','full'],default='quick')
    ap.add_argument('--output',type=Path,required=True,help='New directory; existing directories are rejected.')
    args=ap.parse_args();out=args.output.resolve();assert not out.exists(),'Choose a new output directory'
    start=time.monotonic();integrity=validate();out.mkdir(parents=True)
    replay=out/'COL17_REPLAY.json'
    subprocess.run([sys.executable,str(ROOT/'replay_col17_mp21.py'),'--mode','full' if args.mode=='full' else 'example','--output',str(replay)],check=True,cwd=ROOT)
    result=dict(status='PASS',mode=args.mode,integrity=integrity,col17=json.loads(replay.read_text()),
                environment=dict(python=platform.python_version(),system=platform.system(),machine=platform.machine()),
                independent_biological_validation='NOT_ESTABLISHED',fresh_raw_fastq_extraction='NOT_EXECUTED')
    if args.mode=='full':
        import verify_alignment_oracle_mp21 as oracle
        oracle.OUT=out;oracle.main()
        result['independent_oracle']=json.loads((out/'INDEPENDENT_ORACLE.json').read_text())
    result['seconds']=round(time.monotonic()-start,3)
    (out/'REPRODUCTION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('col17','independent_oracle')},indent=2))
if __name__=='__main__':main()
