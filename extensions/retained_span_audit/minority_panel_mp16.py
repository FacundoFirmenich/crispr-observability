"""Retrospective DNA-discovered recurrent nonmodal sequence events, separately corroborated in cDNA."""
import pathlib,json,hashlib,gzip,collections,datetime
from recompute_panel import encode,classify_intervals,maximize,NEG,solve
import optimal_lattice as lattice
R=pathlib.Path(__file__).resolve().parent;S=R/'ADMISSION_MP15B';O=R/'NONMODAL_MP16'
sha=lambda b:hashlib.sha256(b).hexdigest()
p=json.loads((O/'PROTOCOL.json').read_text());candidates=[];sequences={}
for clone,ts in sorted(p['targets'].items(),key=lambda z:int(z[0])):
 seqs=json.loads((S/f'clone{clone}_segments.json').read_text());sequences[clone]=seqs;ct=collections.Counter()
 with gzip.open(S/f'clone{clone}_pairs.jsonl.gz','rt') as f:
  for line in f:
   j=json.loads(line)
   if j['library']!='DNA' or j['event_relative_JEB'] is None:continue
   hs=set(filter(None,j['segments']+[j['merged_segment']]))
   if len(hs)==1:ct[next(iter(hs))]+=1
 for h,n in sorted(ct.items()):
  if n>=5 and seqs[h] not in [t['sequence'] for t in ts]:candidates.append({'clone':clone,'hash':h,'sequence':seqs[h],'DNA_pairs':n,'net_length_JEB':len(seqs[h])-120})
protocol={'scope':__doc__,'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'threshold_DNA_pairs':5,'selection':'All DNA nonmodal spans with >=5 sequence-unambiguous length-determinate supporting pairs; no cDNA or P1 selection. Prior source reads have been seen, so retrospective, not blind. These are sequence events, not adjudicated biological alleles.','candidates':candidates,'references':p['references'],'deltas':p['deltas'],'script_sha256':sha(pathlib.Path(__file__).read_bytes())}
pp=O/'MINORITY_PROTOCOL.json';assert not pp.exists();pp.write_text(json.dumps(protocol,indent=2)+'\n')
rna={}
for clone in sequences:
 ct=collections.Counter()
 with gzip.open(S/f'clone{clone}_pairs.jsonl.gz','rt') as f:
  for line in f:
   j=json.loads(line)
   if j['library']!='cDNA' or j['event_relative_JEB'] is None:continue
   hs=set(filter(None,j['segments']+[j['merged_segment']]))
   if len(hs)==1:ct[next(iter(hs))]+=1
 rna[clone]=ct
results=[];witnesses=[]
for row in candidates:
 q=row['sequence'];target=q[10:-10];defs=[]
 for ref in p['references']:
  case={'target_kind':'HEK3_DNA','reference_sequence':ref['sequence'],'observed_sequence':q,'lo':ref['lo'],'hi':ref['hi'],'programmed_local':target,'unedited_local':next(x['sequence'][10:-10] for x in p['references'] if x['name']=='JEB')}
  opt,w=solve(case,lattice,54);deficit=[opt-w['programmed_local'][k]['score'] if k in w['programmed_local'] else None for k in range(3)]
  defs.append({'reference':ref['name'],'optimum':opt,'deficits':deficit});witnesses.append({'clone':row['clone'],'hash':row['hash'],'reference':ref['name'],'optimum':opt,'witnesses':w})
 calls=[]
 for d in p['deltas']:
  sets=[{lab for lab,z in zip(('NE','E','U'),x['deficits']) if z is not None and z<=d} for x in defs];calls.append(next(iter(sets[0])) if len(sets[0])==1 and sets[0]==sets[1] else 'U')
 results.append(dict(row,cDNA_pairs=rna[row['clone']][row['hash']],reference_scores=defs,calls=calls))
out={'scope':__doc__,'candidates':len(results),'distinct_sequences':len({x['hash'] for x in results}),'cDNA_observed_ge1':sum(x['cDNA_pairs']>=1 for x in results),'cDNA_supported_ge5':sum(x['cDNA_pairs']>=5 for x in results),'DNA_support_total':sum(x['DNA_pairs'] for x in results),'cDNA_support_total':sum(x['cDNA_pairs'] for x in results),'certified_by_delta':{str(d):sum(x['calls'][i]=='E' for x in results) for i,d in enumerate(p['deltas'])},'results':results,'reference_witness_cases':len(witnesses)}
(O/'MINORITY_RESULT.json').write_text(json.dumps(out,indent=2)+'\n');(O/'MINORITY_WITNESSES.json').write_text(json.dumps(witnesses,separators=(',',':')));print(json.dumps({k:v for k,v in out.items() if k!='results'}))
