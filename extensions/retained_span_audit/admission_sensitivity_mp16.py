import pathlib,json,gzip,collections
R=pathlib.Path(__file__).resolve().parent;O=R/'NONMODAL_MP16';S=R/'ADMISSION_MP15B'
p=json.loads((O/'PROTOCOL.json').read_text());ct=collections.Counter();support=collections.Counter()
for c,ts in p['targets'].items():
 scores=json.loads((O/f'clone{c}_scores.json').read_text());seq=json.loads((S/f'clone{c}_segments.json').read_text())
 with gzip.open(S/f'clone{c}_pairs.jsonl.gz','rt') as f:
  for line in f:
   j=json.loads(line);hs=set(filter(None,j['segments']+[j['merged_segment']]));lib=j['library'];strict=j['event_relative_JEB'] is not None
   group='no_span' if not hs else 'modal_only' if all(seq[h] in [t['sequence'] for t in ts] for h in hs) else 'nonmodal'
   support[lib,'MP15B_admitted' if strict else 'MP15B_unresolved']+=1
   for di,d in enumerate(p['deltas']):
    cc=[]
    for ti in range(len(ts)):
     labels={scores[h]['calls'][ti][di] for h in hs};cc.append(next(iter(labels)) if len(labels)==1 and 'U' not in labels else 'U')
    call='any_E' if 'E' in cc else 'all_NE' if cc and all(x=='NE' for x in cc) else 'unresolved'
    for policy in ('any_available_span','MP15B_gate'):
     actual=call if policy=='any_available_span' or strict else 'unresolved';ct[lib,policy,d,group,actual]+=1
rows=[dict(zip(('library','policy','delta','group','call'),k),pairs=v) for k,v in sorted(ct.items())]
(O/'ADMISSION_SENSITIVITY.json').write_text(json.dumps({'scope':'Separate scoring availability from the stricter MP15B paired admission contract; no gate substitution hidden.','counts':rows},indent=2)+'\n')
for lib in ('DNA','cDNA'):
 for policy in ('any_available_span','MP15B_gate'):
  a=collections.Counter()
  for x in rows:
   if x['library']==lib and x['policy']==policy and x['delta']==36:a[x['call']]+=x['pairs']
  print(lib,policy,dict(a))
