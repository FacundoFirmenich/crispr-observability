"""Same-clone DNA/cDNA exact distal-span evidence, without biological truth labels."""
import pathlib,json,gzip,collections,difflib
R=pathlib.Path(__file__).resolve().parent;O=R/'ADMISSION_MP15B'
p=json.loads((O/'PROTOCOL.json').read_text());out=[]
for c in sorted(p['groups'],key=int):
 seq=json.loads((O/f'clone{c}_segments.json').read_text());counts={t:collections.Counter() for t in ('DNA','cDNA')}
 with gzip.open(O/f'clone{c}_pairs.jsonl.gz','rt') as f:
  for line in f:
   j=json.loads(line)
   if j['event_relative_JEB'] is None:continue
   hs=set(h for h in j['segments']+[j['merged_segment']] if h is not None)
   if len(hs)==1:counts[j['library']][next(iter(hs))]+=1
 for e in p['published'][c]:
  k=e['indel_relative_JEB'];ct={t:collections.Counter({h:n for h,n in counts[t].items() if len(seq[h])-120==k}) for t in counts};shared=set(ct['DNA'])&set(ct['cDNA'])
  out.append({'clone':c,'event':k,'shared_exact_sequences':len(shared),'DNA_pairs_shared':sum(ct['DNA'][h] for h in shared),'cDNA_pairs_shared':sum(ct['cDNA'][h] for h in shared),'modal':{t:[{'hash':h,'sequence':seq[h],'pairs':n} for h,n in ct[t].most_common(1)] for t in ct}})
result={'scope':__doc__,'events':out,'events_shared_ge100_each_library':sum(min(x['DNA_pairs_shared'],x['cDNA_pairs_shared'])>=100 for x in out)}
# Preserve unexpected Clone3 -6 evidence separately from the published +6 target.
seq=json.loads((O/'clone3_segments.json').read_text());cc={t:collections.Counter() for t in ('DNA','cDNA')}
with gzip.open(O/'clone3_pairs.jsonl.gz','rt') as f:
 for line in f:
  j=json.loads(line)
  if j['event_relative_JEB']==-6:
   for h in set(x for x in j['segments']+[j['merged_segment']] if x):cc[j['library']][h]+=1
control=json.loads((R/'P1_FREEZE.json').read_text())['controls'][1]['source']['top'][0][0][90:210]
result['clone3_minus6']={'published_plus6_kept_unchanged':True,'reference_JEB':control,'modal_observations':{t:[{'sequence':seq[h],'pairs_with_sequence':n,'descriptive_alignment_not_truth':list(difflib.SequenceMatcher(None,control,seq[h],autojunk=False).get_opcodes())} for h,n in v.most_common(1)] for t,v in cc.items()}}
(O/'SEQUENCE_SUPPORT.json').write_text(json.dumps(result,indent=2)+'\n');print('SHARED_SEQUENCE_POSITIVE_EVENTS_GE100_BOTH',result['events_shared_ge100_each_library'],'/',len(out))
