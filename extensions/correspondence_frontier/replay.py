"""Offline actual-data replay of MP06 certificates, external transport and optimal DAGs."""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys
HERE=Path(__file__).resolve().parent
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work',type=Path,required=True);p.add_argument('--part',choices=['all','mechanism','external','cftr_lattice','hek3_lattice'],default='all');a=p.parse_args();w=a.work.resolve();assert not w.exists() and not w.is_relative_to(HERE);w.mkdir(parents=True)
 for name,item in json.loads((HERE/'INPUT_MANIFEST.json').read_text())['files'].items():
  data=(HERE/name).read_bytes();assert len(data)==item['bytes'] and hashlib.sha256(data).hexdigest()==item['sha256'],name
 parts=['mechanism','external','cftr_lattice','hek3_lattice'] if a.part=='all' else [a.part]
 for part in parts:
  out=w/part
  if part=='mechanism':cmd=['analyze_mechanism.py','--work',str(HERE/'inputs/cftr'),'--PE-runs',str(HERE/'inputs/PE_default'),'--focus',str(HERE/'inputs/FOCUS_89.csv')]
  elif part=='external':cmd=['analyze_external.py','--work',str(HERE/'inputs/hek3')]
  else:cmd=['audit_lattice.py','--input',str(HERE/'inputs'/(part.split('_')[0]+'_REAL_ALIGNMENT_CASES.json'))]
  subprocess.run([sys.executable,'-B',str(HERE/cmd[0]),*cmd[1:],'--output',str(out)],check=True)
  expected=HERE/'expected'/part
  assert {p.name for p in out.iterdir()}=={p.name for p in expected.iterdir()},part
  for file in out.iterdir():assert file.read_bytes()==(expected/file.name).read_bytes(),str(file)
  if part=='external':
   from verify_marginals import verify
   verify(out/'PUBLISHED_HUMAN_MARGINALS_IDENTIFICATION.json')
  print('EXACT_SCIENTIFIC_REPLAY',part,flush=True)
 (w/'VERIFICATION.json').write_text(json.dumps({'status':'PASS','parts':parts,'scope':'Exact reconstruction on real frozen inputs, including whole original raw denominators and complete optimal-path projection enumeration. Not new engine execution, independent human reproduction or biological accuracy.'},indent=2)+'\n')
if __name__=='__main__':main()
