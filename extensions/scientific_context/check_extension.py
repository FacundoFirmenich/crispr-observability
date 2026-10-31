"""Rebuild added reference/QC outputs and verify actual saved carrier results.

No network, new alignments, molecular simulations or educational study.
"""
from pathlib import Path
import hashlib,json,subprocess,sys,tempfile
ROOT=Path(__file__).resolve().parent
REPOSITORY=ROOT.parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def main():
 manifest=json.loads((ROOT/'SOURCE_MANIFEST.json').read_text())
 for rel,e in manifest['files'].items():
  p=ROOT/rel;assert p.is_file()and p.stat().st_size==e['bytes']and sha(p)==e['sha256'],rel
 with tempfile.TemporaryDirectory(prefix='crispr-reference-audit-')as td:
  t=Path(td);m=t/'molecular';q=t/'raw_qc';c=t/'chemical'
  subprocess.run([sys.executable,'-B',str(ROOT/'analyze_molecular_references.py'),'--output',str(m)],check=True,stdout=subprocess.DEVNULL)
  for p in m.iterdir():assert p.read_bytes()==(ROOT/'derived/molecular'/p.name).read_bytes(),p.name
  subprocess.run([sys.executable,'-B',str(ROOT/'audit_raw_qc.py'),'--raw-root',str(REPOSITORY/'raw_replay'),'--output',str(q)],check=True,stdout=subprocess.DEVNULL)
  for p in q.iterdir():assert p.read_bytes()==(ROOT/'derived/raw_qc'/p.name).read_bytes(),p.name
  subprocess.run([sys.executable,'-B',str(ROOT/'analyze_chemical_context.py'),'--output',str(c)],check=True,stdout=subprocess.DEVNULL)
  for p in c.iterdir():assert p.read_bytes()==(ROOT/'derived/chemical'/p.name).read_bytes(),p.name
 counts=0
 for name in ['canonical_local_DNA','programmed_local_DNA','declared_native_local_peptide']:
  d=ROOT/'accessible'/name;expected=''.join((ROOT/'derived/molecular'/(name+'.fasta')).read_text().splitlines()[1:])
  manifests=list(d.rglob('*.run.json'));assert len(manifests)==1
  m=json.loads(manifests[0].read_text())
  for a in m['artifacts']:
   p=manifests[0].parent/a['path'];assert p.stat().st_size==a['bytes']and sha(p)==a['sha256']
  for carrier,status in [('decode','exact-manifest-matched'),('audio-decode','exact-wav-embedded'),('score-decode','exact-score-embedded')]:
   r=json.loads((d/carrier/'decode.report.json').read_text());actual=''.join((d/carrier/'sequence.fasta').read_text().splitlines()[1:])
   assert actual==expected and r['exact'] is True and r['status']==status
   assert r['normalized_sequence_sha256']==hashlib.sha256(expected.encode()).hexdigest();counts+=1
 raw=json.loads((ROOT/'derived/raw_qc/RAW_QC_SUMMARY.json').read_text());mol=json.loads((ROOT/'derived/molecular/MOLECULAR_REFERENCE_AUDIT.json').read_text())
 assert raw['raw_records']==177006 and raw['quality_pass_records']==175010 and raw['quality_fail_records']==1996
 assert raw['reads_trimmed_or_removed']==0 and raw['FastQC_or_MultiQC_executed'] is False
 assert mol['canonical_and_programmed_DNA_differ_at_bases']==5
 assert mol['canonical_local_peptide_coordinates_1_inclusive']==[500,513]
 assert sorted(x['sequence_length']for x in mol['reference_set_local_peptide_counterexample'])==[604,1419,1480]
 print(json.dumps({'status':'PASS','source_files_verified':len(manifest['files']),'reference_raw_and_chemical_results_rebuilt_byte_identically':True,'matched_chemical_comparators':3,'raw_records':177006,'historical_quality_failures_reproduced':1996,'stored_actual_carrier_roundtrip_reports_verified':counts,'carrier_decoders_rerun_by_this_CI':False,'new_biological_samples_or_accuracy_calibration':False}))

if __name__=='__main__':main()
