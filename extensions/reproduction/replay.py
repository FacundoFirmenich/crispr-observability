"""Clean, offline replay of recovered original geometry, catalogue and contract code.

This is an internal computational reproduction, not molecular truth calibration.
Every scientific output is recomputed before comparison with historical outputs.
"""
from pathlib import Path
import argparse, json, hashlib, shutil, subprocess, sys, os, gzip
BASE=Path(__file__).resolve().parent;REPO=BASE.parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def expected_files(folder):
 for p in sorted(folder.glob('*.json*')):
  data=gzip.decompress(p.read_bytes()) if p.name.endswith('.gz') else p.read_bytes()
  yield p.name.removesuffix('.gz'),hashlib.sha256(data).hexdigest()
def run(cmd,log,pythonpath=None):
 env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
 if pythonpath is not None:env['PYTHONPATH']=str(pythonpath)
 with log.with_suffix('.stdout').open('w') as o,log.with_suffix('.stderr').open('w') as e:p=subprocess.run(cmd,stdout=o,stderr=e,env=env)
 if p.returncode:raise RuntimeError('Execution failed; see '+str(log))
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work',type=Path,required=True);p.add_argument('--only',choices=['all','geometry','external','cftr','contract'],default='all');a=p.parse_args()
 w=a.work.resolve();assert not w.exists() and not w.is_relative_to(REPO),'Use a new directory outside the repository';w.mkdir(parents=True)
 source=json.loads((BASE/'SOURCE_MANIFEST.json').read_text())['files']
 for name,rec in source.items():assert sha(BASE/name)==rec['sha256'] and (BASE/name).stat().st_size==rec['bytes'],name
 report={'status':'RUNNING','original_source_files_verified':len(source),'scope':'Clean internal computational reproduction. No new biological records, accuracy or independent human reproduction.'}
 if a.only in ['all','geometry']:
  g=BASE/'geometry';root=w/'v22';shutil.copytree(g/'source',root);raw=root/'v21/sources/fastq';raw.mkdir(parents=True)
  runs=json.loads((root/'v21/results/READ_PHASE_RESULTS.json').read_text());seen=set()
  for r in runs:
   h=r['sequence_quality_stream_sha256']
   if h in seen:continue
   seen.add(h);name=r['run_accession']+'.fastq.gz';src=REPO/'raw_replay/inputs/fastq'/name
   if not src.exists():src=g/'additional_fastq'/name
   shutil.copy2(src,raw/name)
  out=w/'geometry';out.mkdir()
  run([sys.executable,str(g/'original/separate_quality_geometry.py'),'--root',str(root),'--out',str(out),'--design',str(g/'original/DESIGN.json')],w/'geometry_main')
  run([sys.executable,str(g/'original/veto_common_quality.py'),'--root',str(root),'--out',str(out)],w/'geometry_veto')
  run([sys.executable,str(g/'original/verify_m1.py'),str(out)],w/'geometry_independent')
  matched={}
  for name,expected_hash in expected_files(g/'expected'):
   # The veto embeds the exact unchanged classifier hash. Pure scientific files are byte compared.
   assert sha(out/name)==expected_hash,name;matched[name]=expected_hash
  report['geometry']={'status':'PASS','unique_streams':len(seen),'raw_records':sum(r['N'] for r in json.loads((out/'VALIDATION.json').read_text())['inputs']),'scientific_files_byte_identical':matched}
 if a.only in ['all','external']:
  src=BASE/'external';dst=w/'external';shutil.copytree(src,dst);(dst/'results').mkdir()
  for name in ['run_external.py','verify_external.py','explain_witnesses.py']:run([sys.executable,str(dst/'code'/name)],w/name)
  matched={}
  for name,expected_hash in expected_files(src/'expected'):assert sha(dst/'results'/name)==expected_hash,name;matched[name]=expected_hash
  run([sys.executable,'-m','unittest','discover','-s',str(dst/'tests')],w/'external_tests')
  report['external']={'status':'PASS','scientific_files_byte_identical':matched,'engine_merging_or_alignment_rerun':False}
 if a.only in ['all','cftr']:
  src=BASE/'cftr';dst=w/'cftr';shutil.copytree(src,dst);(dst/'results').mkdir();shutil.copytree(BASE/'external/software/src',dst/'src')
  for rec in json.loads((dst/'INPUT_MANIFEST.json').read_text()):
   rel=Path(rec['path'])
   if (dst/rel).exists():continue
   name=rel.name;kind=rel.parts[1]
   source=REPO/'raw_replay/inputs/fastq'/name if kind=='raw' else REPO/'raw_replay/inputs/alleles'/kind.removesuffix('_alleles')/name
   (dst/rel).parent.mkdir(exist_ok=True);shutil.copy2(source,dst/rel)
  for name in ['reconstruct_and_audit.py','compose_targets.py','verify_independent.py']:run([sys.executable,str(dst/'code'/name)],w/name)
  matched={}
  for name,expected_hash in expected_files(src/'expected'):assert sha(dst/'results'/name)==expected_hash,name;matched[name]=expected_hash
  run([sys.executable,'-m','unittest','discover','-s',str(dst/'tests')],w/'cftr_tests',dst/'src')
  report['cftr']=json.loads((dst/'results/INDEPENDENT_VERIFICATION.json').read_text())
  report['cftr']['scientific_files_byte_identical']=matched
 if a.only in ['all','contract']:
  src=BASE/'contract';dst=w/'contract';shutil.copytree(src,dst);(dst/'results').mkdir();(dst/'delivery').mkdir()
  run([sys.executable,str(dst/'code/event_contract.py'),'build'],w/'contract_build')
  name='OBSERVED_EVENT_CONTRACT_1_0.json';assert sha(dst/'delivery'/name)==sha(src/'expected'/name)
  report['contract']=json.loads((dst/'results/CONTRACT_VERIFICATION.json').read_text())
 report['status']='PASS';(w/'REPRODUCTION_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
