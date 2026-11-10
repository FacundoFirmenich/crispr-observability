"""Independent replay from source FASTQ; alternative exact-overlap implementation."""
import pathlib,json,gzip,hashlib,collections,csv
from transport_other_clones import download
R=pathlib.Path(__file__).resolve().parent;O=R/'ADMISSION_MP15B'
def sha(b):return hashlib.sha256(b).hexdigest()
def extract(s,p):
 a=p['left'];b=p['right']
 if s.count(a)!=1 or s.count(b)!=1:return None
 i=s.find(a);j=s.find(b)
 if j<i+len(a):return None
 q=s[i:j+len(b)]
 return q if 'N' not in q else None
def records(b):
 lines=gzip.decompress(b).decode('ascii').splitlines();assert len(lines)%4==0
 for i in range(0,len(lines),4):
  h,s,plus,quality=lines[i:i+4];assert h.startswith('@') and plus.startswith('+') and len(s)==len(quality)
  yield h.split()[0].removesuffix('/1').removesuffix('/2'),s
def main():
 p=json.loads((O/'PROTOCOL.json').read_text());report=[];events=[];agg={t:collections.Counter() for t in ('DNA','cDNA')}
 for c,files in sorted(p['groups'].items(),key=lambda z:int(z[0])):
  result=json.loads((O/f'clone{c}.json').read_text());sq=json.loads((O/f'clone{c}_segments.json').read_text())
  assert all(sha(v.encode())==k for k,v in sq.items())
  ledger=O/f'clone{c}_pairs.jsonl.gz';assert sha(ledger.read_bytes())==result['ledger_sha256']
  observed={t:[] for t in ('DNA','cDNA')}
  with gzip.open(ledger,'rt') as f:
   for line in f:
    j=json.loads(line);observed[j['library']].append(j)
  for row in result['libraries']:
   typ=row['library'];raw=[]
   for side in ('fwd','rev'):
    source=next(f for f in files if f['key']==f'C17Ex52_{typ} immoJEB_PPN_Clone{c}_{side}.fastq.gz')
    b=download(source);assert sha(b)==next(f['sha256'] for f in result['files'] if f['key']==source['key']);raw.append(list(records(b)))
   assert len(raw[0])==len(raw[1])==len(observed[typ]);dist=collections.Counter();ct=collections.Counter();transition=collections.Counter()
   for (id,s),(id2,t),j in zip(*raw,observed[typ]):
    assert id==id2==j['id'];t=t.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
    q1,q2=extract(s,p),extract(t,p);d=None;m=None;merged=None
    if q1 is not None and q2 is not None and len(q1)==len(q2):d=len(q1)-120
    candidates=[]
    if set(s+t)<=set('ACGT'):
     start=0
     while True:
      i=s.find(t[:30],start)
      if i<0:break
      k=len(s)-i
      if 30<=k<=len(t) and s[i:]==t[:k]:candidates.append(k)
      start=i+1
     if len(candidates)==1:
      merged=extract(s+t[candidates[0]:],p)
      if merged is not None:m=len(merged)-120
    choices=set(x for x in (d,m) if x is not None);event=next(iter(choices)) if len(choices)==1 else None
    if q1 is not None and q2 is not None and len(q1)!=len(q2):event=None
    assert event==j['event_relative_JEB'] and d==j['direct_event'] and m==j['merged_event']
    for q,h in zip((q1,q2),j['segments']):assert (q is None and h is None) or sq[h]==q
    assert (merged is None and j['merged_segment'] is None) or sq[j['merged_segment']]==merged
    ct['raw']+=1;ct['determinate_length']+=int(event is not None);ct['unresolved']+=int(event is None)
    transition[int(j['old_both']),int(event is not None)]+=1
    if event is not None:dist[event]+=1
   for k,v in ct.items():assert v==row['counts'][k]
   assert dict(dist)=={int(k):v for k,v in row['length_distribution'].items() if v}
   agg[typ].update(row['counts']);agg[typ]['newly_determinate_vs_old_admission']+=transition[0,1];agg[typ]['lost_vs_old_admission']+=transition[1,0]
   for e in row['events']:
    assert e['E']==dist[e['indel_relative_JEB']] and e['E']+e['NE']+e['U']==ct['raw']
    events.append(dict(clone=c,library=typ,raw=ct['raw'],**e))
   report.append({'clone':c,'library':typ,'raw_source_replayed':ct['raw'],'status':'PASS'})
  print('INDEPENDENT_RAW_REPLAY_PASS',c,flush=True)
 out={'status':'ALL_RAW_PAIRS_REPLAYED','pairs':sum(x['raw_source_replayed'] for x in report),'source_files':72,'libraries':report,'aggregate':{k:dict(v) for k,v in agg.items()},'event_strata':len(events),'support_ge100':sum(e['E']>=100 for e in events),'unsupported':[e for e in events if e['E']<100],'events':events,'validator_sha256':sha(pathlib.Path(__file__).read_bytes())}
 (O/'VALIDATION.json').write_text(json.dumps(out,indent=2)+'\n')
 with (O/'EVENT_COUNTS.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(events[0]));w.writeheader();w.writerows(events)
 print(json.dumps({k:v for k,v in out.items() if k not in ('events','libraries')}),flush=True)
if __name__=='__main__':main()
