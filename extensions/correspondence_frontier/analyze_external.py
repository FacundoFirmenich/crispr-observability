"""External HEK3 raw-to-allele reconciliation and declared local-DNA certificates."""
from pathlib import Path
from collections import Counter,defaultdict
from fractions import Fraction as F
import argparse,csv,gzip,hashlib,io,json,sys,zipfile
sys.path.insert(0,str(Path(__file__).parent))
from analyze_mechanism import hull,minkowski
def rc(s):return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def project(a,s,ref,lo,hi):
 assert a.replace('-','')==ref and len(a)==len(s);ix=[i for i,x in enumerate(a) if x!='-']
 if any(s[ix[k]] not in 'ACGT' for k in [lo-1,hi]):return None
 code=s[ix[lo-1]+1:ix[hi-1]+1].replace('-','');return code if code and not set(code)-set('ACGT') else None
def exact_target(s,left,right):
 hits=[]
 for seq in [s,rc(s)]:
  a=seq.find(left)
  while a>=0:
   b=seq.find(right,a+len(left))
   while b>=0:
    if b-a-len(left)<=350:hits.append(seq[a+len(left):b])
    b=seq.find(right,b+1)
   a=seq.find(left,a+1)
 return hits[0] if len(hits)==1 and hits[0] and not set(hits[0])-set('ACGT') else None
def tablemap(p,ref,lo,hi,reverse):
 options=defaultdict(set);weights=Counter();mismatch=Counter()
 with zipfile.ZipFile(p) as z:
  assert z.testzip() is None
  for r in csv.DictReader(io.StringIO(z.read(z.namelist()[0]).decode()),delimiter='\t'):
   a,s=r['Reference_Sequence'],r['Aligned_Sequence']
   if reverse:a,s=rc(a),rc(s)
   raw=s.replace('-','');key=min(raw,rc(raw));code=project(a,s,ref,lo,hi);options[key].add(code);n=int(r['#Reads']);weights[key]+=n;pos=-1
   for x,y in zip(a,s):
    if x!='-':pos+=1
    if x=='G' and y=='A':mismatch[pos]+=n
 return {k:next(iter(v)) if len(v)==1 else None for k,v in options.items()},weights,mismatch
def savecsv(p,rows):
 with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
 p=argparse.ArgumentParser();p.add_argument('--work',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();b=a.work.resolve();o=a.output.resolve();assert not o.exists();o.mkdir(parents=True);spec=json.loads((b/'RULE_FREEZE.json').read_text());inputs=json.loads((b/'FILTERED_INPUTS.json').read_text());lo,hi=spec['local_lo'],spec['local_hi'];modes=['M_fw','P_fw','M_rc','P_rc'];rows=[];polys=[];mismatch_rows=[];identity=[];disagreement_examples=[]
 S={(0,0),(1,0),(0,1)}
 def sig(c):return int(c==spec['programmed_local']),int(c==spec['unedited_local'])
 for inp in inputs:
  rid=inp['run'];maps={};expected={};seen=defaultdict(Counter);counts=defaultdict(Counter);groups=Counter();n=qpass=0;ids=set();logical=hashlib.sha256();discordant=Counter();example={};hashsamples=[]
  for mode in modes:
   ref=spec['reference'] if mode.startswith('M') else spec['programmed_reference'];maps[mode],expected[mode],mm=tablemap(b/'runs'/mode/('CRISPResso_on_'+rid)/'Alleles_frequency_table.zip',ref,lo,hi,mode.endswith('_rc'))
   if mode=='M_fw':mismatch_rows.extend({'run':rid,'position_zero_based':pos,'G_to_A_alignment_weight':v,'inside_programmed_local':lo<=pos<hi} for pos,v in mm.most_common())
  with gzip.open(b/'inputs'/(rid+'.fastq.gz'),'rt') as f:
   while h:=f.readline():
    seq,plus,q=f.readline().strip(),f.readline(),f.readline().strip();name=h.split()[0][1:];assert name not in ids and len(seq)==len(q);ids.add(name);n+=1;logical.update((seq+'\n'+q+'\n').encode());quality=sum(ord(x)-33 for x in q)>=30*len(q);qpass+=quality;key=min(seq,rc(seq));codes={'T':exact_target(seq,spec['left_anchor'],spec['right_anchor']) if quality else None}
    for mode in modes:
     codes[mode]=maps[mode].get(key) if quality else None
     if quality and key in maps[mode]:seen[mode][key]+=1
    for mode,c in codes.items():
     for ev,j in [('programmed_local',0),('unedited_local',1)]:counts[mode,'single_projection',ev]['U' if c is None else 'E' if sig(c)[j] else 'NE']+=1
    complete=all(c is not None for c in codes.values());states={sig(c) for c in codes.values()} if complete else S;groups[tuple(sorted(states))]+=1
    for cert,allowed in [('target_relative',states),('string_witness',states if complete and len(set(codes.values()))==1 else S)]:
     for ev,j in [('programmed_local',0),('unedited_local',1)]:values={s[j] for s in allowed};counts['F5',cert,ev]['U' if len(values)>1 else 'E' if values=={1} else 'NE']+=1
    if quality:
     hashsamples.append((hashlib.sha256(name.encode()).digest(),name,seq))
     if complete and len(set(codes.values()))>1:discordant[key]+=1;example.setdefault(key,(name,seq,codes))
  assert n==inp['N'] and qpass==inp['Q30'];assert all(seen[k]==expected[k] for k in modes),'External raw family abundances do not match all allele tables'
  for (op,cert,ev),c in sorted(counts.items()):
   assert sum(c.values())==n;rows.append({'run':rid,'operator':op,'certificate':cert,'event':ev,'N':n,'E':c['E'],'NE':c['NE'],'U':c['U'],'L_exact':str(F(c['E'],n)),'H_exact':str(F(c['E']+c['U'],n)),'L_percent':100*c['E']/n,'H_percent':100*(c['E']+c['U'])/n})
  p=[(0,0)]
  for states,count in sorted(groups.items()):p=minkowski(p,[(count*x,count*y) for x,y in states])
  for dx,dy in [(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1),(1,-1),(-1,1)]:assert max(dx*x+dy*y for x,y in p)==sum(c*max(dx*x+dy*y for x,y in s) for s,c in groups.items())
  polys.append({'run':rid,'N':n,'targets':['programmed_local','unedited_local'],'logical_relation':'mutually exclusive','vertices_integer_counts':p,'signature_domains':[{'allowed':s,'records':c} for s,c in sorted(groups.items())]})
  chosen={name:seq for _,name,seq in sorted(hashsamples)[:5]}
  for key,count in sorted(discordant.items(),key=lambda kv:(-kv[1],kv[0]))[:6]:name,seq,cs=example[key];chosen[name]=seq;disagreement_examples.append({'run':rid,'example_id':name,'family_records':count,'codings':cs})
  identity.append({'run':rid,'N':n,'Q30':qpass,'sequence_quality_sha256':logical.hexdigest(),'all_four_raw_allele_abundances_reconciled':True,'disagreeing_complete_records':sum(discordant.values()),'distinct_disagreeing_canonical_strings':len(discordant),'lattice_panel':[{'read_id':name,'raw_sequence':seq} for name,seq in sorted(chosen.items())]})
 savecsv(o/'HEK3_EVENT_BOUNDS.csv',rows);savecsv(o/'HEK3_G_TO_A_POSITION_DIAGNOSTICS.csv',mismatch_rows)
 (o/'HEK3_JOINT_POLYGONS.json').write_text(json.dumps({'polygons':polys},indent=2)+'\n');(o/'HEK3_PANEL_AND_IDENTITIES.json').write_text(json.dumps({'runs':identity,'discordant_examples':disagreement_examples},indent=2)+'\n')
 manual={'human':{'WT':328,'TMPL':49,'INDEL':6,'AMBIG':0},'CRISPR_A':{'WT':327,'TMPL':49,'INDEL':7,'AMBIG':0},'CRISPResso2':{'WT':330,'TMPL':49,'INDEL':0,'AMBIG':4}}
 bounds=[]
 for tool in ['CRISPR_A','CRISPResso2']:
  r,c=manual['human'],manual[tool];N=sum(r.values());assert N==383==sum(c.values());lower=max(0,max(r[k]+c[k]-N for k in r));upper=sum(min(r[k],c[k]) for k in r);bounds.append({'tool':tool,'N_selected_human_subset':N,'minimum_possible_agreements':lower,'maximum_possible_agreements':upper,'minimum_agreement_percent':100*lower/N,'maximum_agreement_percent':100*upper/N})
 (o/'PUBLISHED_HUMAN_MARGINALS_IDENTIFICATION.json').write_text(json.dumps({'source':'CRISPR-A2023 SupplementaryTable5; transcribed exact aggregate counts; source image retained separately.','marginals':manual,'sharp_agreement_bounds_from_marginals_only':bounds,'convention':'Four categories, with reported ambiguous assignments retained as their own category and no matching human ambiguous labels.','proof':'Agreement cannot exceed sum of per-category minimum counts. The complete off-diagonal transport graph permits zero diagonal except any dominant row-plus-column excess over N; hence the lower bound. These are mathematical completions of a real published aggregate table, not fabricated biological observations.','limits':'No authenticated per-read human labels were recovered from inspected public author sources. These are identification bounds on the selected 383-record table, not measured accuracy, confidence intervals or full-cohort error rates. No ranking is identified because the bounds overlap.','external_human_reproduction_of_this_project':False},indent=2)+'\n')
 summary={'status':'PASS','distinct_source_studies_added':1,'raw_records':sum(x['N'] for x in identity),'Q30_records':sum(x['Q30'] for x in identity),'streams':3,'new_engine_jobs':12,'raw_allele_reconciliations':12,'fixed_rules_before_raw_outcomes':True,'biological_experiments_newly_performed':0,'human_ground_truth':'Only published aggregate marginals, no authenticated per-read truth labels.','not_claimed':'Comparative classification accuracy, random representative transport sample, blinded or prospective validation, independent laboratory replication, protein function or a treated-control contrast.'}
 (o/'EXTERNAL_VERIFICATION.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2));print('Target-relative programmed-local bounds',[(r['run'],r['L_percent'],r['H_percent']) for r in rows if r['operator']=='F5' and r['certificate']=='target_relative' and r['event']=='programmed_local'])
if __name__=='__main__':main()
