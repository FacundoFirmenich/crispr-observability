"""Exercise the exact alignment envelope on real CFTR adverses and hash comparators."""
from pathlib import Path
from collections import Counter
import argparse,csv,gzip,hashlib,importlib.util,json,sys
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(Path(__file__).parent))
from optimal_lattice import audit,rescore,project
def module(name,p):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
raw=module('raw',ROOT/'raw_replay/run_raw_replay.py')
def rc(s):return s.translate(raw.COMP)[::-1]
def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--input',type=Path);p.add_argument('--focus-new',type=Path);p.add_argument('--focus-old',type=Path);a=p.parse_args();o=a.output.resolve();assert not o.exists();o.mkdir(parents=True)
 if a.input:cases=json.loads(a.input.read_text())['cases']
 else:
  from CRISPResso2 import CRISPResso2Align as c2
  import numpy as np
  refs=[raw.load(ROOT/'raw_replay/metadata'/n) for n in ['REFERENCE.json','REFERENCE_PROGRAMMED.json']]
  wanted={r['read_id'] for r in csv.DictReader(a.focus_new.open())};wanted.update(r['read_id'] for r in json.loads(a.focus_old.read_text()))
  design=raw.load(ROOT/'engine_replay/dependencies/m2/DESIGN.json');chosen=[]
  for run in design['runs']:
   rid=run['run_accession'];records=list(raw.fastq(ROOT/'raw_replay/inputs/fastq'/(rid+'.fastq.gz')))
   quality=[r for r in records if sum(ord(x)-33 for x in r[2])>=30*len(r[2])]
   controls=sorted(quality,key=lambda r:hashlib.sha256(r[0].encode()).digest())[:3]
   chosen.extend((rid,r,'known_adverse' if r[0] in wanted else 'hash_comparator') for r in quality if r[0] in wanted or r in controls)
  cases=[];matrix=c2.make_matrix()
  for rid,(name,seq,q),why in chosen:
   for ri,ref in enumerate(refs):
    f=[c2.global_align(s,ref['sequence'],matrix,np.zeros(len(ref['sequence'])+1,dtype=np.int64),-20,-2) for s in [seq,rc(seq)]];f=max(f,key=lambda x:x[2]);oriented=f[0].replace('-','')
    rev=[c2.global_align(s,rc(ref['sequence']),matrix,np.zeros(len(ref['sequence'])+1,dtype=np.int64),-20,-2) for s in [seq,rc(seq)]];rev=max(rev,key=lambda x:x[2]);back=(rc(rev[0]),rc(rev[1]),rev[2]);assert back[0].replace('-','')==oriented
    cases.append({'run':rid,'read_id':name,'selection':why,'reference':'mutant' if ri==0 else 'programmed','reference_sequence':ref['sequence'],'observed_sequence':oriented,'lo':ref['core_start_0'],'hi':ref['coding_end_exclusive_0'],'native_peptide':ref['native_peptide'],'exact_programmed_coding':ref['exact_programmed_core']+'AT','C2_forward':{'observed':f[0],'reference':f[1],'reported_identity_percent':f[2]},'C2_reverse_normalized':{'observed':back[0],'reference':back[1],'reported_identity_percent':back[2]}})
  (o/'REAL_ALIGNMENT_CASES.json').write_text(json.dumps({'timing':'Known CFTR adverses plus three lowest SHA256(Q30 read identifier) comparators per stream; no synthetic data.','engine':'CRISPResso2 2.3.4 unmodified extension; zero incentive; -20/-2; make_matrix default EDNAFULL-compatible ACGTN entries','cases':cases},indent=2)+'\n')
 results=[]
 for c in cases:
  r=audit(c['reference_sequence'],c['observed_sequence'],c['lo'],c['hi'])
  n=len(c['reference_sequence']);reverse=audit(rc(c['reference_sequence']),rc(c['observed_sequence']),n-c['hi'],n-c['lo'],boundary_convention='right_insertions')
  r['complete_set_orientation_equivariant']=(reverse['optimal_integer_objective']==r['optimal_integer_objective'] and reverse['optimal_path_count']==r['optimal_path_count'] and sorted(rc(s) for s in reverse['valid_coding_strings'])==r['valid_coding_strings'] and reverse['some_optimal_projection_unresolved']==r['some_optimal_projection_unresolved'])
  assert r['complete_set_orientation_equivariant'],'Full optimal-path set fails orientation check'
  for k in ['C2_forward','C2_reverse_normalized']:
   t=c[k];s=rescore(t['reference'],t['observed']);assert s==r['optimal_integer_objective'],(c['read_id'],c['reference'],k,s,r['optimal_integer_objective'])
   coding=project(t['reference'],t['observed'],c['lo'],c['hi']);assert coding is None and r['some_optimal_projection_unresolved'] or coding in r['valid_coding_strings']
  fw=project(c['C2_forward']['reference'],c['C2_forward']['observed'],c['lo'],c['hi']);rv=project(c['C2_reverse_normalized']['reference'],c['C2_reverse_normalized']['observed'],c['lo'],c['hi'])
  def native(s):return s==c['programmed_local'] if c.get('target_kind')=='HEK3_DNA' else len(s)%3==0 and ''.join(raw.CODONS[s[i:i+3]] for i in range(0,len(s),3))==c['native_peptide']
  def second(s):return s==c['unedited_local'] if c.get('target_kind')=='HEK3_DNA' else s==c['exact_programmed_coding']
  values=r['valid_coding_strings'];r.update(run=c['run'],read_id=c['read_id'],reference=c['reference'],selection=c['selection'],C2_forward_coding=fw,C2_reverse_coding=rv,orientation_changes_coding=fw!=rv,native_values=sorted({native(s) for s in values}),exact_values=sorted({second(s) for s in values}))
  r['native_invariant']=not r['some_optimal_projection_unresolved'] and len(r['native_values'])==1;r['exact_invariant']=not r['some_optimal_projection_unresolved'] and len(r['exact_values'])==1
  results.append(r)
 if cases[0].get('target_kind')=='HEK3_DNA':
  for r in results:
   for old,new in [('native_values','programmed_local_values'),('exact_values','unedited_local_values'),('native_invariant','programmed_local_invariant'),('exact_invariant','unedited_local_invariant')]:r[new]=r.pop(old)
 (o/'OPTIMAL_ALIGNMENT_ENVELOPES.json').write_text(json.dumps({'cases':results},indent=2)+'\n')
 first,second_id=('programmed_local_invariant','unedited_local_invariant') if cases[0].get('target_kind')=='HEK3_DNA' else ('native_invariant','exact_invariant')
 summary={'status':'PASS','actual_read_reference_cases':len(results),'distinct_actual_read_ids':len({r['read_id'] for r in results}),'C2_tracebacks_verified_optimal':2*len(results),'actual_cases_with_full_path_enumeration':sum(r['all_paths_independently_enumerated'] for r in results),'complete_set_orientation_equivariant_cases':sum(r['complete_set_orientation_equivariant'] for r in results),'orientation_changes_coding_cases':sum(r['orientation_changes_coding'] for r in results),'multiple_valid_coding_cases':sum(len(r['valid_coding_strings'])>1 for r in results),first+'_cases':sum(r[first] for r in results),second_id+'_cases':sum(r[second_id] for r in results),'whole_string_invariant_cases':sum(r['all_optimal_strings_identical'] for r in results),'some_optimal_projection_unresolved_cases':sum(r['some_optimal_projection_unresolved'] for r in results),'largest_optimal_path_count':max(r['optimal_path_count'] for r in results),'not_claimed':'Classification accuracy, independent biological truth, all scoring objectives or random representative sampling.'}
 (o/'LATTICE_VERIFICATION.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
