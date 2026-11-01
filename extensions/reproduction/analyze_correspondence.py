"""Reconstruct certificate ablations and local synonymous variants from real records.

The stronger certificate verifies string correspondence and usually costs
resolution. It is not an accuracy comparison or a full-allele/function assay.
"""
from pathlib import Path
from collections import Counter,defaultdict
from fractions import Fraction as F
import argparse,csv,gzip,hashlib,importlib.util,itertools,json
ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('raw',ROOT/'raw_replay/run_raw_replay.py');raw=importlib.util.module_from_spec(spec);spec.loader.exec_module(raw)
OPS=('target','mutant_projection','programmed_projection')
def peptide(s):return ''.join(raw.CODONS[s[i:i+3]] for i in range(0,len(s),3)) if s and len(s)%3==0 else None
def records():
 raw.verify_inputs();design=json.loads((ROOT/'engine_replay/dependencies/m2/DESIGN.json').read_text());refs=[raw.load(ROOT/'raw_replay/metadata'/n) for n in ['REFERENCE.json','REFERENCE_PROGRAMMED.json']]
 for r in design['runs']:
  rid=r['run_accession'];maps=[raw.allele_map(ROOT/'raw_replay/inputs/alleles'/kind/(rid+'.zip'),ref)[0] for kind,ref in zip(['mutant','programmed'],refs)]
  for name,s,q in raw.fastq(ROOT/'raw_replay/inputs/fastq'/(rid+'.fastq.gz')):
   qual=sum(ord(c)-33 for c in q)>=30*len(q);codes=[None]*3
   if qual:
    rc=s.translate(raw.COMP)[::-1];codes[0]=raw.target_observer([(s,q,'F'),(rc,q[::-1],'R')],refs[0]['exact_programmed_core'])['coding']
    for i,mp in enumerate(maps):
     keys={x for x in [s,rc] if x in mp};assert len(keys)<=1
     if keys:codes[i+1]=mp[next(iter(keys))][1][0]
   yield r,name,s,qual,dict(zip(OPS,codes))
def write_csv(p,rows):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=a.output.resolve();assert not out.exists();out.mkdir(parents=True)
 exact=raw.load(ROOT/'raw_replay/metadata/REFERENCE.json')['exact_programmed_core']+'AT';subsets=[s for n in [1,2,3] for s in itertools.combinations(OPS,n)]
 counts=defaultdict(Counter);vcounts=Counter();examples={};meta={};nraw=Counter();same_label_distinct=Counter();rawh={}
 for r,name,s,q,codes in records():
  rid=r['run_accession'];meta[rid]=r;nraw[rid]+=1
  vals=list(codes.values())
  if all(x is not None for x in vals) and len(set(vals))>1 and len({peptide(x)=='GTIKENIIFGVSYD' for x in vals})==1:same_label_distinct[rid]+=1
  for subset in subsets:
   cs=[codes[x] for x in subset];st=cs[0] if all(x is not None for x in cs) and len(set(cs))==1 else None
   for event in ['native_local','exact_programmed']:
    label='U' if st is None else ('E' if (peptide(st)=='GTIKENIIFGVSYD' if event=='native_local' else st==exact) else 'NE')
    counts[rid,'+'.join(subset),event][label]+=1
  if all(x is not None for x in vals) and len(set(vals))==1 and peptide(vals[0])=='GTIKENIIFGVSYD' and vals[0]!=exact:
   vcounts[rid,vals[0]]+=1;examples.setdefault((rid,vals[0]),name)
 rows=[]
 for (rid,sub,event),c in sorted(counts.items()):
  n=nraw[rid];assert sum(c.values())==n
  rows.append({'run':rid,'preparation':meta[rid]['source_replica'],'condition':meta[rid]['condition'],'certificate':sub,'event':event,'N':n,'E':c['E'],'NE':c['NE'],'U':c['U'],'L_exact':str(F(c['E'],n)),'H_exact':str(F(c['E']+c['U'],n))})
 write_csv(out/'CORRESPONDENCE_ABLATION_COUNTS.csv',rows)
 contrasts=[]
 for t in rows:
  if t['condition']!='PE treated':continue
  for baseline in ['mock','untreated']:
   c=next(x for x in rows if x['preparation']==t['preparation'] and x['condition']==baseline and x['certificate']==t['certificate'] and x['event']==t['event'])
   lo=F(t['L_exact'])-F(c['H_exact']);hi=F(t['H_exact'])-F(c['L_exact']);contrasts.append({k:t[k] for k in ['preparation','certificate','event']}|{'control':baseline,'L_exact':str(lo),'H_exact':str(hi),'L_pp':float(100*lo),'H_pp':float(100*hi),'width_pp':float(100*(hi-lo))})
 write_csv(out/'CORRESPONDENCE_ABLATION_CONTRASTS.csv',contrasts)
 variants=[]
 for (rid,s),n in sorted(vcounts.items()):
  diffs=[{'local_codon_1based':i//3+1,'canonical_residue':500+i//3,'programmed_codon':exact[i:i+3],'observed_codon':s[i:i+3],'amino_acid':raw.CODONS[s[i:i+3]]} for i in range(0,len(s),3) if s[i:i+3]!=exact[i:i+3]]
  assert len(s)==len(exact)==42 and all(raw.CODONS[x['programmed_codon']]==x['amino_acid'] for x in diffs)
  variants.append({'run':rid,'preparation':meta[rid]['source_replica'],'condition':meta[rid]['condition'],'observed_coding':s,'coding_sha256':hashlib.sha256(s.encode()).hexdigest(),'records':n,'example_read_id':examples[rid,s],'synonymous_codon_differences':diffs})
 report={'status':'PASS','raw_records':sum(nraw.values()),'certificate_subsets':len(subsets),'native_label_agreement_with_unequal_strings':sum(same_label_distinct.values()),'native_nonexact_programmed_witness_records':sum(vcounts.values()),'native_nonexact_distinct_strings':len({s for rid,s in vcounts}),'native_nonexact_by_condition':dict(Counter({condition:sum(n for (rid,s),n in vcounts.items() if meta[rid]['condition']==condition) for condition in ['PE treated','mock','untreated']})),'variants':variants,'limits':'Observed local coding strings and synonymous differences only. A sequencing error, mixed product or true allele origin is not adjudicated. No full-allele, expression, function, independent truth or optimal-observer inference.'}
 assert report['raw_records']==177006 and report['native_label_agreement_with_unequal_strings']==18 and report['native_nonexact_programmed_witness_records']==128
 (out/'NATIVE_NONEXACT_LOCAL_STRINGS.json').write_text(json.dumps(report,indent=2)+'\n')
 w='+'.join(OPS);assert sum(x['U'] for x in rows if x['certificate']==w and x['event']=='native_local')==6281
 ratios=[]
 for prep in [1,2,3]:
  q=[x for x in contrasts if x['preparation']==prep and x['event']=='native_local' and x['control']=='mock'];m=next(x for x in q if x['certificate']=='mutant_projection');wc=next(x for x in q if x['certificate']==w);ratios.append({'preparation':prep,'witness_to_mutant_width_ratio':wc['width_pp']/m['width_pp'],'witness_width_pp':wc['width_pp'],'mutant_width_pp':m['width_pp']})
 (out/'ABLATION_SUMMARY.json').write_text(json.dumps({'status':'PASS','raw_records':177006,'count_rows':len(rows),'contrast_rows':len(contrasts),'width_ratios':ratios,'subset_certificate_semantics':'Every required operator resolves and returns exactly the same observed coding string. Fewer required operators give a weaker correspondence certificate. These certificates do not measure accuracy.','native_nonexact_witness_records':128,'native_nonexact_distinct_strings':report['native_nonexact_distinct_strings']},indent=2)+'\n')
 print(json.dumps({'status':'PASS','records':177006,'native_nonexact_records':128,'distinct_strings':report['native_nonexact_distinct_strings'],'width_ratios':ratios},indent=2))
if __name__=='__main__':main()
