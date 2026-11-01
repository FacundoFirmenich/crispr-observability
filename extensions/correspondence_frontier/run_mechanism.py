"""Frozen retrospective mechanism ablation; same raw records, no biological simulation."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse,datetime,gzip,hashlib,json,subprocess,time
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,v):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2)+'\n')
def rc(s):return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def prepare(work,engine):
 assert not work.exists() and not work.is_relative_to(ROOT);work.mkdir(parents=True)
 assert subprocess.check_output([str(engine),'--version'],text=True).strip()=='CRISPResso 2.3.4'
 refs=[json.loads((ROOT/'raw_replay/metadata'/p).read_text()) for p in ['REFERENCE.json','REFERENCE_PROGRAMMED.json']]
 runs=json.loads((ROOT/'engine_replay/dependencies/m2/DESIGN.json').read_text())['runs'];con=json.loads((ROOT/'extensions/reproduction/pe_mode/CONSTRUCT.json').read_text());jobs=[];inputs=[]
 for r in runs:
  rid=r['run_accession'];src=ROOT/'raw_replay/inputs/fastq'/(rid+'.fastq.gz');dst=work/'inputs'/(rid+'.Q30.fastq.gz');dst.parent.mkdir(exist_ok=True);n=k=0;h=hashlib.sha256()
  with gzip.open(src,'rt') as f,dst.open('wb') as out,gzip.GzipFile(filename='',mode='wb',fileobj=out,mtime=0) as g:
   while head:=f.readline():
    s,plus,q=f.readline().strip(),f.readline(),f.readline().strip();assert head.startswith('@') and plus.startswith('+') and len(s)==len(q);n+=1;h.update((s+'\n'+q+'\n').encode())
    if sum(ord(c)-33 for c in q)>=30*len(q):k+=1;g.write((head.rstrip()+'\n'+s+'\n+\n'+q+'\n').encode())
  assert n==r['raw_reads'] and h.hexdigest()==r['sequence_quality_stream_sha256']
  with gzip.open(dst,'rt') as f:assert sum(1 for _ in f)==4*k
  inputs.append({'run':rid,'N':n,'Q30':k,'source_sha256':sha(src),'filtered_sha256':sha(dst)})
  for mode in ['M_rc','P_rc','MP_fw','MP_rc','PE_zero']:
   cmd=[str(engine),'-r1',str(dst),'-n',rid,'-o',str(work/'runs'/mode),'--n_processes','1','--write_detailed_allele_table','--fastq_output','--suppress_report','--suppress_plots']
   if mode=='PE_zero':
    cmd+=['-a',rc(refs[0]['sequence']),'-an','observed_mutant','--prime_editing_pegRNA_spacer_seq',con['spacer_RNA_5to3'],'--prime_editing_pegRNA_extension_seq',con['RTT_RNA_5to3']+con['PBS_RNA_5to3'],'--prime_editing_pegRNA_scaffold_seq',con['scaffold_RNA_5to3'],'--prime_editing_nicking_guide_seq',con['nsgRNA_spacer_RNA_5to3'],'--needleman_wunsch_gap_incentive','0']
   else:
    inds=[0] if mode=='M_rc' else [1] if mode=='P_rc' else [0,1];rev=mode.endswith('_rc');seqs=[rc(refs[i]['sequence']) if rev else refs[i]['sequence'] for i in inds];wins=[]
    for i in inds:
     ref=refs[i];lo,hi=ref['core_start_0'],ref['coding_end_exclusive_0']
     if rev:lo,hi=len(ref['sequence'])-hi,len(ref['sequence'])-lo
     wins.append(f'{lo}-{hi-1}')
    cmd+=['-a',','.join(seqs),'-an',','.join('observed_mutant' if i==0 else 'programmed_reference' for i in inds),'--default_min_aln_score','60','--needleman_wunsch_gap_open','-20','--needleman_wunsch_gap_extend','-2','--needleman_wunsch_gap_incentive','0','--quantification_window_coordinates',','.join(wins),'--exclude_bp_from_left','0','--exclude_bp_from_right','0','--min_average_read_quality','0','--min_single_bp_quality','0','--min_bp_quality_or_N','0']
   jobs.append({'mode':mode,'run':rid,'command':cmd})
 dump(work/'JOBS.json',jobs);dump(work/'INPUTS.json',inputs)
 dump(work/'FREEZE.json',{'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'timing':'Known retrospective CFTR data; MP05 results known. New ablation outputs not inspected before this command freeze. Not prospective biological validation.','source_commit':'1fc4b629e5400a48d6569b0ab4b8d1766f33636f','source_script_sha256':sha(Path(__file__)),'jobs_sha256':sha(work/'JOBS.json'),'inputs_sha256':sha(work/'INPUTS.json'),'factors':'Reverse-complement orientation, one versus competing reference assignment, and removal of dedicated PE gap incentives. This is a declared partial factorial; unsupported forward PE orientation is not invented.','primary_outputs':'Per-record coding correspondence after normalization, raw-N event hulls, same-record label transitions and attribution patterns for the 89 MP05 discordances. Native tool metrics separate.','limits':'Operational sensitivity, not comparative biological accuracy, independent truth or a new cohort.'})
 print(json.dumps({'state':'FROZEN','jobs':len(jobs),'raw_records':sum(x['N'] for x in inputs)}),flush=True)
def execute(work,j):
 p=work/'executions'/j['mode']/(j['run']+'.json')
 if p.exists():assert json.loads(p.read_text())['state']=='COMPLETED';return
 inp=next(x for x in json.loads((work/'INPUTS.json').read_text()) if x['run']==j['run']);assert sha(work/'inputs'/(j['run']+'.Q30.fastq.gz'))==inp['filtered_sha256']
 obj={**j,'state':'RUNNING','start_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()};dump(p,obj);t=time.monotonic()
 with p.with_suffix('.stdout').open('w') as o,p.with_suffix('.stderr').open('w') as e:r=subprocess.run(j['command'],stdout=o,stderr=e)
 obj.update(state='COMPLETED' if r.returncode==0 else 'FAILED',returncode=r.returncode,seconds=time.monotonic()-t);dump(p,obj);print(json.dumps({'mode':j['mode'],'run':j['run'],'state':obj['state'],'seconds':obj['seconds']}),flush=True);assert r.returncode==0
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','run']);p.add_argument('--work',type=Path,required=True);p.add_argument('--engine',type=Path);p.add_argument('--workers',type=int,default=3);a=p.parse_args();w=a.work.resolve()
 if a.action=='prepare':assert a.engine;prepare(w,a.engine.resolve());return
 f=json.loads((w/'FREEZE.json').read_text());assert sha(w/'JOBS.json')==f['jobs_sha256'] and sha(w/'INPUTS.json')==f['inputs_sha256'] and sha(Path(__file__))==f['source_script_sha256']
 with ThreadPoolExecutor(max_workers=a.workers) as pool:list(pool.map(lambda j:execute(w,j),json.loads((w/'JOBS.json').read_text())))
 dump(w/'RUN_STATUS.json',{'state':'COMPLETED','jobs':45,'biological_cohorts_added':0})
if __name__=='__main__':main()
