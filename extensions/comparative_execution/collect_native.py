from pathlib import Path
import urllib.request,tarfile,io,json,csv,hashlib,argparse
base='https://github.com/FacundoFirmenich/crispr-observability/releases/download/comparison-evidence-43493655364/'
def read_asset(name):
 b=urllib.request.urlopen(base+name).read();t=tarfile.open(fileobj=io.BytesIO(b));return b,t
def main(output):
 assert not output.exists(), "Choose a new result path"
 out={'source_release':base,'rows':[],'execution':[],'sources':[]};families={}
 for offset in [0,4,8]:
  name='amplican-evidence.tar.gz' if offset==0 else f'amplican-anchor-{offset}.tar.gz';b,t=read_asset(name);out['sources'].append({'url':base+name,'sha256':hashlib.sha256(b).hexdigest()})
  parse=lambda n:json.loads(t.extractfile('evidence/'+n).read())
  manifest={r['run_accession']:r for r in parse('INPUT_MANIFEST.json')};ex=parse('EXECUTION.json');out['execution'].append({'offset':offset,'native_seconds':ex['steps'][-1]['seconds'],'total_seconds':ex['elapsed_seconds']})
  rows=list(csv.DictReader(io.TextIOWrapper(t.extractfile('evidence/native/config_summary.csv'))))
  unresolved={};families[offset]={}
  for x in csv.DictReader(io.TextIOWrapper(t.extractfile('evidence/native/alignments/unassigned_reads.csv'))):
   unresolved[x['Barcode']]=unresolved.get(x['Barcode'],0)+int(x['Total']);key=(x['Barcode'],x['Forward'],x['Reverse']);assert key not in families[offset];families[offset][key]=int(x['Total'])
  for x in rows:
   rid=x['ID'];m=manifest[rid];assigned=int(x['Reads']);assert assigned+unresolved.get(rid,0)==m['Q30']
   out['rows'].append({'run':rid,'condition':m['condition'],'preparation':m['source_replica'],'anchor_offset':offset,'raw_N':m['raw_reads'],'Q30':m['Q30'],'assigned':assigned,'unassigned_Q30':unresolved.get(rid,0),'native_filtered':int(x['Reads_Filtered']),'native_edited':int(x['Reads_Edited']),'native_HDR':int(x['HDR'])})
 ref=json.loads(t.extractfile('evidence/REFERENCE.json').read())['sequence'];out['admission_transitions']=[]
 for a,b in [(0,4),(4,8),(0,8)]:
  for rid in sorted(manifest):
   ua={k:v for k,v in families[a].items() if k[0]==rid};ub={k:v for k,v in families[b].items() if k[0]==rid}
   for k in ua.keys()&ub.keys():assert ua[k]==ub[k]
   admitted={k:v for k,v in ua.items() if k not in ub};lost={k:v for k,v in ub.items() if k not in ua};gain=sum(admitted.values());loss=sum(lost.values())
   old=next(x for x in out['rows'] if x['run']==rid and x['anchor_offset']==a);new=next(x for x in out['rows'] if x['run']==rid and x['anchor_offset']==b);assert new['assigned']-old['assigned']==gain-loss
   out['admission_transitions'].append({'run':rid,'from_offset':a,'to_offset':b,'newly_assigned_weighted_records':gain,'newly_unassigned_weighted_records':loss,'newly_assigned_matching_reference_4_24_at_read_4_24':sum(v for k,v in admitted.items() if k[1][4:24]==ref[4:24]),'linkage':'Exact full observed forward/reverse sequence family and unchanged multiplicity within run; not molecular identity or retained original read IDs'})
 output.write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps([x for x in out['rows'] if x['condition']=='PE treated'],indent=2));print(out['execution'])
 print('total assigned',[(a,sum(r['assigned'] for r in out['rows'] if r['anchor_offset']==a)) for a in [0,4,8]])
 print('weighted transitions',[(a,b,sum(x['newly_assigned_weighted_records'] for x in out['admission_transitions'] if x['from_offset']==a and x['to_offset']==b),sum(x['newly_unassigned_weighted_records'] for x in out['admission_transitions'] if x['from_offset']==a and x['to_offset']==b)) for a,b in [(0,4),(4,8),(0,8)]])
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
