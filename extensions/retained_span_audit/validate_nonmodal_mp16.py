"""Exhaustive alignment-witness rescore/reprojection and separate ledger aggregation."""
import pathlib,json,hashlib,gzip,collections,csv
from recompute_panel import solve
import optimal_lattice as lattice
R=pathlib.Path(__file__).resolve().parent;O=R/'NONMODAL_MP16';S=R/'ADMISSION_MP15B'
sha=lambda b:hashlib.sha256(b).hexdigest()
def main():
 p=json.loads((O/'PROTOCOL.json').read_text());aggregate=collections.Counter();support=collections.Counter();reports=[];witnesses=[];allrows=[];distinct=set();physical=collections.Counter()
 for clone,ts in sorted(p['targets'].items(),key=lambda z:int(z[0])):
  seqs=json.loads((S/f'clone{clone}_segments.json').read_text());scores=json.loads((O/f'clone{clone}_scores.json').read_text());result=json.loads((O/f'clone{clone}_counts.json').read_text());assert sha((O/f'clone{clone}_scores.json').read_bytes())==result['scores_sha256'];distinct.update(seqs)
  selected={h for x in result['witness_candidates'] for h in x['hashes']};count=0;cases=0
  for h,q in seqs.items():
   for ref,expected in zip(p['references'],scores[h]['reference_scores']):
    case={'target_kind':'HEK3_DNA','reference_sequence':ref['sequence'],'observed_sequence':q,'lo':ref['lo'],'hi':ref['hi'],'programmed_local':ts[0]['sequence'][10:-10],'unedited_local':ts[-1]['sequence'][10:-10]}
    opt,w=solve(case,lattice,54);assert opt==expected['optimum'];cases+=1
    for i,name in enumerate(('programmed_local','unedited_local')):
     assert [opt-w[name][k]['score'] if k in w[name] else None for k in range(3)]==expected['deficits'][i]
     count+=len(w[name])
    if h in selected:witnesses.append({'clone':clone,'segment_hash':h,'reference':ref['name'],'optimum':opt,'targets':w})
  ct=collections.Counter();ss=collections.Counter();n=collections.Counter()
  with gzip.open(S/f'clone{clone}_pairs.jsonl.gz','rt') as f:
   for line in f:
    j=json.loads(line);lib=j['library'];n[lib]+=1;hs=set(filter(None,j['segments']+[j['merged_segment']]))
    group='no_span' if not hs else 'modal_only' if all(seqs[h] in [t['sequence'] for t in ts] for h in hs) else 'nonmodal';ss[lib,group]+=1
    for di,d in enumerate(p['deltas']):
     pc=[]
     for ti,t in enumerate(ts):
      calls=[]
      for h in hs:
       sets=[]
       for ref in scores[h]['reference_scores']:
        sets.append([lab for lab,z in zip(('NE','E','U'),ref['deficits'][ti]) if z is not None and z<=d])
       calls.append(sets[0][0] if len(sets[0])==1 and sets[0]==sets[1] else 'U')
      call=calls[0] if calls and calls[0]!='U' and all(x==calls[0] for x in calls) else 'U';ct[lib,t['label'],d,group,call]+=1;pc.append(call)
     physical[lib,d,group,'any_E' if 'E' in pc else 'all_NE' if pc and all(x=='NE' for x in pc) else 'unresolved']+=1
  assert dict(n)==result['raw']
  for x in result['counts']:
   key=x['library'],x['target'],x['delta'],x['group'],x['call'];assert ct[key]==x['pairs'];aggregate[x['library'],x['delta'],x['group'],x['call']]+=x['pairs'];allrows.append(dict(clone=clone,**x))
  support.update(ss);reports.append({'clone':clone,'unique_spans':len(seqs),'reference_cases':cases,'attaining_witnesses_rescored_reprojected':count,'raw_pairs':sum(n.values()),'status':'PASS'})
  print('VALIDATED',clone,len(seqs),cases,count,flush=True)
 out={'status':'ALL_RETAINED_SPANS_AND_ALL_PAIR_AGGREGATES_VALIDATED','clones':reports,'raw_pairs':sum(x['raw_pairs'] for x in reports),'clone_span_cases':sum(x['unique_spans'] for x in reports),'distinct_sequence_hashes':len(distinct),'witnesses_checked':sum(x['attaining_witnesses_rescored_reprojected'] for x in reports),'support':[{'library':k[0],'group':k[1],'pairs':v} for k,v in sorted(support.items())], 'aggregate_event_comparisons':[dict(zip(('library','delta','group','call'),k),comparisons=v) for k,v in sorted(aggregate.items())], 'physical_pair_counts':[dict(zip(('library','delta','group','call'),k),pairs=v) for k,v in sorted(physical.items())],'validator_sha256':sha(pathlib.Path(__file__).read_bytes())}
 (O/'VALIDATION.json').write_text(json.dumps(out,indent=2)+'\n');(O/'SELECTED_WITNESSES.json').write_text(json.dumps(witnesses,separators=(',',':')))
 with (O/'COMPLETE_COUNTS.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(allrows[0]));w.writeheader();w.writerows(allrows)
 print('ALL_VALIDATED',out['raw_pairs'],out['clone_span_cases'],out['witnesses_checked'],flush=True)
if __name__=='__main__':main()
