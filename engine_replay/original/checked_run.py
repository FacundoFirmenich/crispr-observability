"""Pre/post immutable-input checks around each frozen execution block."""
import argparse, gzip, hashlib, json, os, subprocess
import run_grid as g
def verify(label):
    results=[]
    for m in json.loads((g.BASE/'INPUT_MANIFEST.json').read_text()):
        p=g.BASE/'inputs'/(m['run']+'.Q30.fastq.gz')
        content=p.read_bytes();assert hashlib.sha256(content).hexdigest()==m['input_sha256']
        lines=gzip.decompress(content).splitlines();assert len(lines)%4==0;h=hashlib.sha256();seen=set()
        for i in range(0,len(lines),4):
            assert lines[i].startswith(b'@') and lines[i+2].startswith(b'+') and len(lines[i+1])==len(lines[i+3])
            assert lines[i] not in seen;seen.add(lines[i]);h.update(lines[i+1]+b'\n'+lines[i+3]+b'\n')
        assert len(lines)//4==m['Q30_pass'] and h.hexdigest()==m['Q30_sequence_quality_stream_sha256']
        results.append({'run':m['run'],'sha256':hashlib.sha256(content).hexdigest(),'Q30_stream_sha256':h.hexdigest(),'Q30_pass':len(lines)//4,'state':'PASS'})
    g.dump(g.BASE/'input_checks'/(label+'.json'),{'utc':g.now(),'checks':results,'state':'PASS'})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['baseline','stock','oat']);p.add_argument('--workers',default=3,type=int);a=p.parse_args()
    verify(a.mode+'_before')
    r=subprocess.run(['python3',str(g.BASE/'run_profiles.py'),a.mode,'--workers',str(a.workers)])
    verify(a.mode+'_after')
    assert r.returncode==0
