"""Retrospective exact-sequence P1 stress test on DNA-modal positive-supported spans."""
import pathlib,json,hashlib,datetime
from recompute_panel import encode,classify_intervals,maximize,NEG
R=pathlib.Path(__file__).resolve().parent;O=R/'ADMISSION_MP15B'
sha=lambda b:hashlib.sha256(b).hexdigest()
f=json.loads((R/'P1_FREEZE.json').read_text());support=json.loads((O/'SEQUENCE_SUPPORT.json').read_text())
refs=[]
for c in f['controls']:
 full=c['source']['top'][0][0];a=full.index('GGGGACAGCCGCCTCCTGTC');b=full.index('ATGAGCACAGGAGGAGGTGG')+20
 q=full[a:b];refs.append({'name':'HKC' if 'HKC' in c['file'] else 'JEB','sequence':q,'lo':10,'hi':len(q)-10})
protocol={'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':__doc__,'window':'Trim ten bases at each distal span end; target is DNA-modal span trimmed identically. Original scoring, both reference backgrounds, delta grid unchanged. No selection by P1 outcome. Modal sequence selection is retrospective, not independent validation.', 'deltas':[0,9,18,27,36,45,54],'references':refs,'targets':[{'clone':x['clone'],'event':x['event'],'DNA_modal':x['modal']['DNA'],'cDNA_modal':x['modal']['cDNA']} for x in support['events']], 'engine_sha256':sha((R/'recompute_panel.py').read_bytes()),'script_sha256':sha(pathlib.Path(__file__).read_bytes())}
pp=O/'P1_POSITIVE_PROTOCOL.json';assert not pp.exists();pp.write_text(json.dumps(protocol,indent=2)+'\n')
results=[]
for row in protocol['targets']:
 if not row['DNA_modal']:
  results.append(dict(row,status='NO_OBSERVED_POSITIVE_TARGET'));continue
 q=row['DNA_modal'][0]['sequence'];target=q[10:-10]
 targets={'observed_modal':target,'unedited_JEB':next(r['sequence'][10:-10] for r in refs if r['name']=='JEB')};labels=classify_intervals(q,targets);deficits=[]
 for r in refs:
  opt,mx,_,_,suffix=maximize(encode(r['sequence']),encode(q),r['lo'],r['hi'],labels,54)
  assert opt==suffix[0,0,0] and all(max(v)==opt for v in mx)
  deficits.append({'reference':r['name'],'optimum':int(opt),'event_deficits':[None if z==NEG else int(opt-z) for z in mx[0]]})
 calls=[]
 for d in protocol['deltas']:
  sets=[{label for label,z in zip(('NE','E','U'),r['event_deficits']) if z is not None and z<=d} for r in deficits]
  calls.append(next(iter(sets[0])) if len(sets[0])==1 and sets[0]==sets[1] else 'U')
 same=bool(row['cDNA_modal']) and row['cDNA_modal'][0]['sequence']==q
 results.append(dict(row,status='COMPUTED',reference_deficits=deficits,calls=calls,same_DNA_cDNA_modal=same))
 print(row['clone'],row['event'],calls,'same_modal',same,flush=True)
out={'scope':__doc__,'protocol_sha256':sha(pp.read_bytes()),'results':results,'deltas':protocol['deltas'],'certified_modal_targets_by_delta':{str(d):sum(x.get('calls',[None]*7)[i]=='E' for x in results) for i,d in enumerate(protocol['deltas'])},'limitations':'Selected modal sequence stress test, not all-read P1 validation, independent biological truth, or sensitivity/specificity. Read support counts represent repeated observed strings, not independent genotypes.'}
(O/'P1_POSITIVE_RESULT.json').write_text(json.dumps(out,indent=2)+'\n');print(out['certified_modal_targets_by_delta'])
