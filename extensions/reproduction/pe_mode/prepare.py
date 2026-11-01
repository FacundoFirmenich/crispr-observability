"""Prepare authenticated, oriented CRISPResso2 2.3.4 PE jobs before running them."""
from pathlib import Path
import argparse,datetime,gzip,hashlib,json,subprocess
ROOT=Path(__file__).resolve().parents[3]
ARCHIVE_SHA256='7e4f2ae9ce9bad27f00afd53abf195431f68e02fc06b39da30d4f7696fca6468'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work',type=Path,required=True);p.add_argument('--engine',type=Path,required=True);p.add_argument('--source-archive',type=Path,required=True);a=p.parse_args()
 w=a.work.resolve();assert not w.exists() and not w.is_relative_to(ROOT);assert sha(a.source_archive)==ARCHIVE_SHA256
 engine=a.engine.resolve();version=subprocess.check_output([str(engine),'--version'],text=True).strip();assert version=='CRISPResso 2.3.4'
 w.mkdir(parents=True);(w/'inputs').mkdir();con=json.loads((Path(__file__).parent/'CONSTRUCT.json').read_text());ref=json.loads((ROOT/'raw_replay/metadata/REFERENCE.json').read_text());s=ref['sequence'].translate(str.maketrans('ACGT','TGCA'))[::-1];design=json.loads((ROOT/'engine_replay/dependencies/m2/DESIGN.json').read_text());jobs=[];inputs=[]
 for r in design['runs']:
  rid=r['run_accession'];src=ROOT/'raw_replay/inputs/fastq'/(rid+'.fastq.gz');dst=w/'inputs'/(rid+'.Q30.fastq.gz');n=k=0;digest=hashlib.sha256()
  with gzip.open(src,'rt') as f,dst.open('wb') as out,gzip.GzipFile(filename='',fileobj=out,mode='wb',mtime=0) as g:
   while h:=f.readline():
    seq,plus,q=f.readline().strip(),f.readline(),f.readline().strip();assert h.startswith('@') and plus.startswith('+') and len(seq)==len(q);n+=1;digest.update((seq+'\n'+q+'\n').encode())
    if sum(ord(c)-33 for c in q)>=30*len(q):k+=1;g.write((h.rstrip()+'\n'+seq+'\n+\n'+q+'\n').encode())
  assert n==r['raw_reads'] and digest.hexdigest()==r['sequence_quality_stream_sha256']
  with gzip.open(dst,'rt') as check:assert sum(1 for _ in check)==4*k
  inputs.append({'run':rid,'N':n,'Q30':k,'raw_sha256':sha(src),'filtered_sha256':sha(dst)})
  cmd=[str(engine),'-r1',str(dst),'-a',s,'-an','observed_mutant','-n',rid,'-o',str(w/'runs'),'--n_processes','1','--write_detailed_allele_table','--fastq_output','--suppress_report','--suppress_plots','--prime_editing_pegRNA_spacer_seq',con['spacer_RNA_5to3'],'--prime_editing_pegRNA_extension_seq',con['RTT_RNA_5to3']+con['PBS_RNA_5to3'],'--prime_editing_pegRNA_scaffold_seq',con['scaffold_RNA_5to3'],'--prime_editing_nicking_guide_seq',con['nsgRNA_spacer_RNA_5to3']]
  jobs.append({'run':rid,'command':cmd})
 dump(w/'INPUT_MANIFEST.json',inputs);dump(w/'JOBS.json',jobs);dump(w/'FREEZE.json',{'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'engine':version,'upstream_commit':'1fc4b629e5400a48d6569b0ab4b8d1766f33636f','source_archive_sha256':ARCHIVE_SHA256,'construct_sha256':sha(Path(__file__).parent/'CONSTRUCT.json'),'jobs_sha256':sha(w/'JOBS.json'),'inputs_sha256':sha(w/'INPUT_MANIFEST.json'),'timing':'New computational preparation of a known retrospective comparison. Not prospective biological validation.','no_sequence_check_or_edited_reference_override':True})
 print(json.dumps({'status':'PREPARED_NO_JOBS_EXECUTED','jobs':len(jobs),'raw_records':sum(x['N'] for x in inputs),'Q30_records':sum(x['Q30'] for x in inputs)}))
if __name__=='__main__':main()
