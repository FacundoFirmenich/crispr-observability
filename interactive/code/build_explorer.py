"""Build the self-contained public explorer from a verified numerical model."""
from pathlib import Path
import argparse,json,hashlib

root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser()
p.add_argument('--model',type=Path,default=root/'evidence_model.json')
p.add_argument('--template',type=Path,default=root/'explorer.template.html')
p.add_argument('--output',type=Path,default=root/'epistemic_graph.html')
a=p.parse_args()
model=json.loads(a.model.read_text())
assert model['schemaVersion']=='1.0'
for row in model['observations']+model['sensitivityObservations']:
 assert row['N']>0 and all(isinstance(row[k],int) and row[k]>=0 for k in ['N','E','NE','U'])
 assert row['E']+row['NE']+row['U']==row['N']
 assert abs(row['L']-row['E']/row['N'])<1e-12
 assert abs(row['H']-(row['E']+row['U'])/row['N'])<1e-12
encoded=json.dumps(model,ensure_ascii=False,separators=(',',':')).replace('<','\\u003c').replace('>','\\u003e').replace('&','\\u0026')
fingerprint=hashlib.sha256(encoded.encode()).hexdigest()
template=a.template.read_text()
assert template.count('__EVIDENCE_MODEL__')==1 and template.count('__EVIDENCE_SHA256__')==1
result=template.replace('__EVIDENCE_MODEL__',encoded).replace('__EVIDENCE_SHA256__',fingerprint)
a.output.parent.mkdir(parents=True,exist_ok=True)
a.output.write_text(result)
receipt={'html':str(a.output.name),'html_sha256':hashlib.sha256(result.encode()).hexdigest(),'embedded_json_sha256':fingerprint,'input_json_sha256':hashlib.sha256(a.model.read_bytes()).hexdigest(),'observations':len(model['observations']),'contrasts':len(model['contrasts']),'sensitivity_event_rows':len(model['sensitivityObservations']),'sensitivity_contrasts':len(model['sensitivityContrasts']),'graph_nodes':len(model['graph']['nodes']),'AND_hyperedges':len(model['graph']['hyperedges']),'dependencies':'none; CSS, JavaScript and frozen data are embedded','public_source_scope':'Concise data/claim model only; no full manuscript or supplement text'}
(root/'BUILD_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
