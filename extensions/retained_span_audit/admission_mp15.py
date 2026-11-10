"""Retrospective distal-anchor length-event validation, separate from exact-string P1."""
import pathlib,json,hashlib,datetime,collections,gzip,re
from transport_other_clones import download,clip
from col17a1_clone10_audit import fetch,parse,rc
R=pathlib.Path(__file__).resolve().parent
O=R/'ADMISSION_MP15'
def sha(b):return hashlib.sha256(b).hexdigest()
def main():
 O.mkdir(exist_ok=True); pp=O/'PROTOCOL.json'
 if not pp.exists():
  f=json.loads((R/'P1_FREEZE.json').read_text()); pub=json.loads((R/'TRANSPORT_MP14/PUBLISHED_RECONCILIATION.json').read_text())
  s=f['controls'][1]['source']['top'][0][0]; a,b=s[90:110],s[190:210]
  assert s.count(a)==s.count(b)==1
  meta=json.loads(fetch('https://zenodo.org/api/records/6242279')); groups={}
  for x in meta['files']:
   m=re.fullmatch(r'C17Ex52_(DNA|cDNA) immoJEB_PPN_Clone(\d+)_(fwd|rev).fastq.gz',x['key'])
   if m and m[2] in pub['clones']:groups.setdefault(m[2],[]).append(x)
  assert len(groups)==18 and all(len(v)==4 for v in groups.values())
  p={'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':__doc__,
   'left':a,'right':b,'JEB_span':120,'HKC_span':122,'JEB_control_positions_zero_based':[90,110,190,210],
   'old_left':f['outer_left'],'old_right':f['outer_right'],'groups':groups,'published':pub['clones'],
   'adequate_observed_positive_support':100,
   'rule':'Unique ordered exact 20-base distal anchors in each oriented mate; no N in clipped span. Both mate span lengths must agree. Event is span length minus JEB120. Preserve sequence-discordant equal-length pairs separately. No merging, imputation or quality correction.',
   'limits':'Retrospective redesign after MP14, not unseen validation. Length equivalence is coarser than allele sequence identity. S2 is same-study aggregate annotation, not per-read truth. Positive support >=100 is an operational coverage threshold, not a powered accuracy study.',
   'script_sha256':sha(pathlib.Path(__file__).read_bytes()),'source_annotation_sha256':sha((R/'TRANSPORT_MP14/PUBLISHED_RECONCILIATION.json').read_bytes())}
  pp.write_text(json.dumps(p,indent=2)+'\n')
 p=json.loads(pp.read_text()); assert p['script_sha256']==sha(pathlib.Path(__file__).read_bytes())
 for clone,files in sorted(p['groups'].items(),key=lambda x:int(x[0])):
  result=O/f'clone{clone}.json'; ledger=O/f'clone{clone}_pairs.jsonl.gz'; segments=O/f'clone{clone}_segments.json'
  if result.exists():
   v=json.loads(result.read_text()); assert sha(ledger.read_bytes())==v['ledger_sha256'] and sha(segments.read_bytes())==v['segments_sha256']; print('REUSED',clone,flush=True);continue
  reads={}; verified=[]; seqs={}; rows=[]
  for x in files:
   z=download(x);reads[x['key']]=parse(z);verified.append(dict(x,sha256=sha(z)))
  with gzip.open(ledger,'wt',encoding='ascii') as out:
   for typ in ('DNA','cDNA'):
    x=reads[f'C17Ex52_{typ} immoJEB_PPN_Clone{clone}_fwd.fastq.gz'];y=reads[f'C17Ex52_{typ} immoJEB_PPN_Clone{clone}_rev.fastq.gz'];assert len(x)==len(y)
    ct=collections.Counter(); dist=collections.Counter(); cross=collections.Counter(); discord=collections.Counter()
    for (id,s),(id2,t) in zip(x,y):
     assert id==id2; t=rc(t); old=all(clip(q,p['old_left'],p['old_right']) is not None for q in (s,t))
     qs=[clip(q,p['left'],p['right']) for q in (s,t)];qs=[q if q is not None and 'N' not in q else None for q in qs]
     hs=[]
     for q in qs:
      h=sha(q.encode()) if q is not None else None;hs.append(h)
      if h:seqs[h]=q
     both=all(q is not None for q in qs); ct['raw']+=1;ct['old_both']+=int(old);ct['new_both']+=int(both);cross[int(old),int(both)]+=1
     event=None
     if both:
      if len(qs[0])==len(qs[1]):
       event=len(qs[0])-p['JEB_span'];dist[event]+=1;ct['determinate_length']+=1
       if qs[0]==qs[1]:ct['identical_sequence']+=1
       else:ct['same_length_different_sequence']+=1;discord[event]+=1
      else:ct['length_conflict']+=1
     ct['unresolved']+=int(event is None)
     out.write(json.dumps({'id':id,'library':typ,'old_both':old,'segments':hs,'event_relative_JEB':event},separators=(',',':'))+'\n')
    assert ct['raw']==ct['determinate_length']+ct['unresolved']
    events=[]
    for e in p['published'][clone]:
     k=e['indel_relative_JEB']; n=dist[k]
     events.append(dict(e,E=n,NE=ct['determinate_length']-n,U=ct['unresolved'],same_length_sequence_discordant=discord[k],support_at_least_100=n>=p['adequate_observed_positive_support']))
    rows.append({'library':typ,'counts':dict(ct),'admission_cross':[{'old':a,'new':b,'pairs':v} for (a,b),v in sorted(cross.items())], 'length_distribution':dict(sorted(dist.items())),'events':events})
  segments.write_text(json.dumps(seqs,separators=(',',':')))
  result.write_text(json.dumps({'clone':clone,'protocol_sha256':sha(pp.read_bytes()),'files':verified,'libraries':rows,'ledger_sha256':sha(ledger.read_bytes()),'segments_sha256':sha(segments.read_bytes())},indent=2)+'\n')
  print('COMPLETE',clone,[(x['library'],x['counts'],[(e['indel_relative_JEB'],e['E']) for e in x['events']]) for x in rows],flush=True)
 print('MP15_ALL_18_COMPLETE',flush=True)
if __name__=='__main__':main()
