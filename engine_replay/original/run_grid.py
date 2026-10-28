"""Fixed finite grid replay on authenticated original single-end reads.

The third-party engine is fetched/installed separately and is never bundled.
No parameter is varied using a window-only surrogate.
"""
from pathlib import Path
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse, csv, datetime, gzip, hashlib, io, json, os, subprocess, time, zipfile

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
OLD = ROOT/'restored/m2/advance_m2'
RAW = ROOT/'publication_verification/P1/advance_m7/replay/project/evidence/v22/v21/sources/fastq'
ENGINE = ROOT/'_mp02_engine/env/bin/CRISPResso'
SOURCE = ROOT/'_mp02_engine/CRISPResso2'
COMMIT = '1fc4b629e5400a48d6569b0ab4b8d1766f33636f'
GRID = [('baseline',-20,-2,60),('open_m10',-10,-2,60),('open_m30',-30,-2,60),('extend_m1',-20,-1,60),('extend_m4',-20,-4,60),('score_50',-20,-2,50),('score_70',-20,-2,70)]
EVENTS = ('full_programmed','native_peptide')
CODON = dict(zip(('TTT TTC TTA TTG TCT TCC TCA TCG TAT TAC TAA TAG TGT TGC TGA TGG CTT CTC CTA CTG CCT CCC CCA CCG CAT CAC CAA CAG CGT CGC CGA CGG ATT ATC ATA ATG ACT ACC ACA ACG AAT AAC AAA AAG AGT AGC AGA AGG GTT GTC GTA GTG GCT GCC GCA GCG GAT GAC GAA GAG GGT GGC GGA GGG').split(), 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))

def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dump(path,obj):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp'); tmp.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n'); tmp.replace(path)
def rows(path):
    with gzip.open(path,'rt') as f:
        while h:=f.readline():
            s=f.readline().strip(); plus=f.readline().strip(); q=f.readline().strip()
            assert h.startswith('@') and plus.startswith('+') and len(s)==len(q)
            yield h.split()[0][1:],s,plus,q
def project(ar,aq,ref):
    assert len(ar)==len(aq) and ar.replace('-','')==ref['sequence']
    lo=ref['core_start_0']; hi=ref['coding_end_exclusive_0']; pos=-1; observed={}; coding=[]
    for a,b in zip(ar,aq):
        if a!='-':
            pos+=1; observed[pos]=b
            if lo<=pos<hi and b!='-': coding.append(b)
        elif lo<=pos+1<hi and b!='-': coding.append(b)
    assert pos+1==len(ref['sequence'])
    if any(observed.get(i) not in 'ACGT' for i in (lo-1,hi)): return ('U','U')
    coding=''.join(coding)
    if not coding or any(c not in 'ACGT' for c in coding): return ('U','U')
    native=len(coding)%3==0 and ''.join(CODON[coding[i:i+3]] for i in range(0,len(coding),3))==ref['native_peptide']
    exact=coding==ref['exact_programmed_core']+'AT'
    return ('E' if exact else 'N','E' if native else 'N')
def summarize(c,N): return {'N':N,'E':c['E'],'NE':c['N'],'U':c['U']}
def prepare():
    assert not (BASE/'GRID_FREEZE.json').exists(), 'Freeze already exists; do not overwrite.'
    design=json.loads((OLD/'DESIGN.json').read_text())
    refs={k:json.loads((OLD/p).read_text()) for k,p in [('mutant','REFERENCE.json'),('programmed','REFERENCE_PROGRAMMED.json')]}
    spec={'id':'MP02_ALIGNMENT_SENSITIVITY','created_utc':now(),'state':'FIXED_BEFORE_NEW_GRID_RESULTS','engine_commit':COMMIT,'scope':'all nine original raw runs; two frozen references; seven one-at-a-time settings; 126 full engine jobs',
      'grid':[{'name':n,'gap_open':o,'gap_extend':e,'min_score':s} for n,o,e,s in GRID],
      'references':refs,'runs':design['runs'],'input_policy':design['quality'],'projection_policy':design['adapter'],
      'fixed_options':{'gap_incentive':0,'exclude_left':0,'exclude_right':0,'n_processes':1,'all_internal_quality_options':0,'trim_sequences':False,'detailed_allele_table':True,'fastq_output':True,'suppress_report':True,'suppress_plots':True},
      'analysis':'Raw-denominator E/NE/U, per-read transitions against freshly rerun baseline, native software quantities separately; exact baseline match to frozen numerical outputs is required before interpreting sensitivity.',
      'limits':'Post hoc finite parameter sensitivity on the same biological genealogy; not error calibration or independent replication. Only open/extend/min-score change; window stays frozen.'}
    dump(BASE/'GRID.json',spec)
    dump(BASE/'GRID_FREEZE.json',{'created_utc':now(),'files':{'GRID.json':sha(BASE/'GRID.json'),'run_grid.py':sha(__file__)},'status':'FROZEN','new_grid_outputs_opened':False})
    inputs=[]
    for r in design['runs']:
        rid=r['run_accession']; src=RAW/(rid+'.fastq.gz'); dest=BASE/'inputs'/(rid+'.Q30.fastq.gz'); dest.parent.mkdir(exist_ok=True)
        stream=hashlib.sha256(); qstream=hashlib.sha256(); ids=set(); N=K=0
        with dest.open('wb') as b,gzip.GzipFile(filename='',mode='wb',fileobj=b,mtime=0) as gz:
            for name,s,plus,q in rows(src):
                assert name not in ids; ids.add(name); N+=1; stream.update((s+'\n'+q+'\n').encode())
                if sum(ord(c)-33 for c in q)>=30*len(q):
                    K+=1; qstream.update((s+'\n'+q+'\n').encode()); gz.write(('@'+name+'\n'+s+'\n+\n'+q+'\n').encode())
        assert N==r['raw_reads'] and stream.hexdigest()==r['sequence_quality_stream_sha256']
        frozen=json.loads((OLD/'results'/(rid+'.execution.json')).read_text()); assert frozen['Q30_pass']==K
        inputs.append({'run':rid,'source':str(src.relative_to(ROOT)),'source_sha256':sha(src),'N_raw':N,'raw_sequence_quality_stream_sha256':stream.hexdigest(),'Q30_pass':K,'Q30_sequence_quality_stream_sha256':qstream.hexdigest(),'input_sha256':sha(dest)})
    dump(BASE/'INPUT_MANIFEST.json',inputs)
    jobs=[]
    for ref in refs:
        for n,o,e,s in GRID:
            for r in design['runs']:
                rid=r['run_accession']; jobs.append({'id':f'{ref}/{n}/{rid}','reference':ref,'setting':n,'run':rid,'state':'UNEXECUTED'})
    dump(BASE/'JOBS.json',jobs)
    print(json.dumps({'status':'PREPARED','jobs':len(jobs),'raw_records':sum(x['N_raw'] for x in inputs),'freeze_sha256':sha(BASE/'GRID_FREEZE.json')}),flush=True)
def command(job,spec):
    ref=spec['references'][job['reference']]; g=next(x for x in spec['grid'] if x['name']==job['setting']); rid=job['run']
    return [str(ENGINE),'-r1',str(BASE/'inputs'/(rid+'.Q30.fastq.gz')),'-a',ref['sequence'],'-an','observed_mutant','-n',rid,'-o',str(BASE/'runs'/job['reference']/job['setting']),
      '--default_min_aln_score',str(g['min_score']),'--needleman_wunsch_gap_open',str(g['gap_open']),'--needleman_wunsch_gap_extend',str(g['gap_extend']),'--needleman_wunsch_gap_incentive','0',
      '--quantification_window_coordinates',f"{ref['core_start_0']}-{ref['coding_end_exclusive_0']-1}",'--exclude_bp_from_left','0','--exclude_bp_from_right','0','--min_average_read_quality','0','--min_single_bp_quality','0','--min_bp_quality_or_N','0','--n_processes','1','--write_detailed_allele_table','--fastq_output','--suppress_report','--suppress_plots']
def execute(job,spec):
    record=BASE/'executions'/job['reference']/job['setting']/(job['run']+'.json')
    if record.exists() and json.loads(record.read_text()).get('state')=='COMPLETED': return json.loads(record.read_text())
    assert not record.exists(), 'Preserve failed run before retrying.'
    cmd=command(job,spec); obj={**job,'command':cmd,'state':'RUNNING','started_utc':now()}; dump(record,obj); t=time.monotonic()
    with record.with_suffix('.stdout').open('w') as o,record.with_suffix('.stderr').open('w') as e:
        p=subprocess.run(cmd,stdout=o,stderr=e)
    obj.update(state='COMPLETED' if p.returncode==0 else 'FAILED',returncode=p.returncode,elapsed_seconds=time.monotonic()-t,finished_utc=now()); dump(record,obj)
    print(json.dumps({'id':job['id'],'state':obj['state'],'seconds':obj['elapsed_seconds']}),flush=True)
    return obj
def run(mode,workers):
    freeze=json.loads((BASE/'GRID_FREEZE.json').read_text())
    for p,h in freeze['files'].items(): assert sha(BASE/p)==h, 'Frozen file changed: '+p
    assert subprocess.check_output(['git','-C',str(SOURCE),'rev-parse','HEAD'],text=True).strip()==COMMIT
    subprocess.run(['git','-C',str(SOURCE),'diff','--exit-code'],check=True)
    spec=json.loads((BASE/'GRID.json').read_text()); jobs=json.loads((BASE/'JOBS.json').read_text())
    if mode=='pilot': jobs=[j for j in jobs if j['reference']=='mutant' and j['setting']=='baseline' and j['run']=='SRR29623272']
    with ThreadPoolExecutor(max_workers=workers) as pool:
        outcomes=list(pool.map(lambda j:execute(j,spec),jobs))
    dump(BASE/('PILOT_STATUS.json' if mode=='pilot' else 'RUN_STATUS.json'),{'state':'COMPLETED' if all(x['state']=='COMPLETED' for x in outcomes) else 'FAILED','finished_utc':now(),'jobs':len(outcomes),'workers':workers})
    assert all(x['state']=='COMPLETED' for x in outcomes)

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('mode',choices=['prepare','pilot','run']); p.add_argument('--workers',type=int,default=3); a=p.parse_args()
    if a.mode=='prepare': prepare()
    else: run(a.mode,a.workers)
