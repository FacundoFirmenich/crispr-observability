"""Exhaustive fixed-target P1 scoring of every retained MP15B distal span, including nonmodal observations."""
import pathlib,json,hashlib,datetime,gzip,collections
from recompute_panel import encode,classify_intervals,maximize,NEG
R=pathlib.Path(__file__).resolve().parent;S=R/'ADMISSION_MP15B';O=R/'NONMODAL_MP16'
sha=lambda b:hashlib.sha256(b).hexdigest()
def main():
 O.mkdir(exist_ok=True);pp=O/'PROTOCOL.json'
 if not pp.exists():
  prior=json.loads((S/'P1_POSITIVE_PROTOCOL.json').read_text());support=json.loads((S/'SEQUENCE_SUPPORT.json').read_text());targets={}
  for row in prior['targets']:
   if row['DNA_modal']:q=row['DNA_modal'][0]['sequence'];label=f"S2_{row['event']:+d}"
   else:
    assert row['clone']=='3' and row['event']==6
    q=support['clone3_minus6']['modal_observations']['DNA'][0]['sequence'];label='observed_minus6_S2_plus6_discordant'
   targets.setdefault(row['clone'],[]).append({'label':label,'sequence':q,'event_relative_JEB':len(q)-120,'source_published_event':row['event']})
  p={'scope':__doc__,'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'references':prior['references'],'deltas':prior['deltas'],'targets':targets,
   'selection':'Every retained span in every clone, no frequency filtering. Frozen MP15 DNA-modal targets; Clone3 observed -6 tracked separately from unobserved published +6. Never redefine each query as its own positive target.',
   'pair_rule':'Score every available spanning mate and exact-merged span. Require all available representations to give the same determinate label; otherwise U. No observed span yields U. All original pairs retained. Report nonmodal if any available full span differs from all frozen target spans.',
   'script_sha256':sha(pathlib.Path(__file__).read_bytes()),'engine_sha256':sha((R/'recompute_panel.py').read_bytes()),'source_protocol_sha256':sha((S/'PROTOCOL.json').read_bytes()),
   'inputs':{p.name:sha(p.read_bytes()) for p in S.glob('clone*.json*') if '_segments' in p.name or '_pairs' in p.name}}
  pp.write_text(json.dumps(p,indent=2)+'\n')
 p=json.loads(pp.read_text());assert sha(pathlib.Path(__file__).read_bytes())==p['script_sha256'];refs=p['references'];ds=p['deltas'];cache={}
 for clone,targets in sorted(p['targets'].items(),key=lambda x:int(x[0])):
  dest=O/f'clone{clone}_scores.json';summary=O/f'clone{clone}_counts.json'
  if summary.exists():print('REUSED',clone,flush=True);continue
  seqpath=S/f'clone{clone}_segments.json';ledger=S/f'clone{clone}_pairs.jsonl.gz'
  assert sha(seqpath.read_bytes())==p['inputs'][seqpath.name] and sha(ledger.read_bytes())==p['inputs'][ledger.name]
  seqs=json.loads(seqpath.read_text());scores={}
  ts={t['label']:t['sequence'][10:-10] for t in targets}
  if len(ts)==1:ts['duplicate_engine_slot']=next(iter(ts.values()))
  key=tuple(ts.values())
  for h,q in seqs.items():
   assert sha(q.encode())==h
   if (h,key) in cache:scores[h]=cache[h,key];continue
   labels=classify_intervals(q,ts);deficits=[]
   for r in refs:
    opt,mx,_,_,suffix=maximize(encode(r['sequence']),encode(q),r['lo'],r['hi'],labels,54)
    assert opt==suffix[0,0,0] and all(max(x)==opt for x in mx)
    deficits.append({'reference':r['name'],'optimum':int(opt),'deficits':[[None if z==NEG else int(opt-z) for z in row] for row in mx]})
   calls=[]
   for ti in range(len(targets)):
    cc=[]
    for d in ds:
     sets=[{label for label,z in zip(('NE','E','U'),ref['deficits'][ti]) if z is not None and z<=d} for ref in deficits]
     cc.append(next(iter(sets[0])) if len(sets[0])==1 and sets[0]==sets[1] else 'U')
    assert all(a==b or b=='U' for a,b in zip(cc,cc[1:]));calls.append(cc)
   scores[h]={'reference_scores':deficits,'calls':calls};cache[h,key]=scores[h]
  dest.write_text(json.dumps(scores,separators=(',',':')))
  counts=collections.Counter();support=collections.Counter();n=collections.Counter();examples={}
  with gzip.open(ledger,'rt') as f:
   for line in f:
    j=json.loads(line);lib=j['library'];hs=set(h for h in j['segments']+[j['merged_segment']] if h);n[lib]+=1
    nonmodal=bool(hs) and any(seqs[h] not in {t['sequence'] for t in targets} for h in hs)
    group='nonmodal' if nonmodal else 'modal_only' if hs else 'no_span';support[lib,group]+=1
    for ti,t in enumerate(targets):
     for di,d in enumerate(ds):
      labels={scores[h]['calls'][ti][di] for h in hs};call=next(iter(labels)) if len(labels)==1 and 'U' not in labels else 'U'
      counts[lib,t['label'],d,group,call]+=1
      if group=='nonmodal' and d in (0,36,45):
       k=(ti,d,call)
       if k not in examples:examples[k]={'id':j['id'],'library':lib,'target_index':ti,'delta':d,'call':call,'hashes':sorted(hs)}
  rows=[{'library':lib,'target':t['label'],'delta':d,'group':g,'call':call,'pairs':counts[lib,t['label'],d,g,call]} for lib in n for t in targets for d in ds for g in ('modal_only','nonmodal','no_span') for call in ('E','NE','U')]
  for lib in n:
   for t in targets:
    for d in ds:assert sum(x['pairs'] for x in rows if x['library']==lib and x['target']==t['label'] and x['delta']==d)==n[lib]
  summary.write_text(json.dumps({'clone':clone,'raw':dict(n),'unique_spans':len(seqs),'scores_sha256':sha(dest.read_bytes()),'counts':rows,'support':[{'library':lib,'group':g,'pairs':v} for (lib,g),v in support.items()],'witness_candidates':list(examples.values())},indent=2)+'\n')
  print('COMPLETE',clone,len(seqs),dict(n),flush=True)
 print('ALL_CLONES_COMPLETE',flush=True)
if __name__=='__main__':main()
