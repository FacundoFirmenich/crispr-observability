"""Per-record mechanism decomposition and task-relative joint identification hulls."""
from pathlib import Path
from collections import Counter,defaultdict
from fractions import Fraction as F
import argparse,csv,gzip,hashlib,importlib.util,itertools,json,re,sys,zipfile,io
ROOT=Path(__file__).resolve().parents[2]
def module(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
raw=module('raw',ROOT/'raw_replay/run_raw_replay.py');cor=module('cor',ROOT/'extensions/reproduction/analyze_correspondence.py')
def rc(s):return s.translate(raw.COMP)[::-1]
def signature(c,exact):
 n=cor.peptide(c)=='GTIKENIIFGVSYD';e=c==exact;assert not e or n;return (int(n),int(e))
STATES={(0,0),(1,0),(1,1)}
def hull(points):
 p=sorted(set(points))
 if len(p)<=1:return p
 def cross(o,a,b):return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
 low=[];high=[]
 for q in p:
  while len(low)>=2 and cross(low[-2],low[-1],q)<=0:low.pop()
  low.append(q)
 for q in reversed(p):
  while len(high)>=2 and cross(high[-2],high[-1],q)<=0:high.pop()
  high.append(q)
 return low[:-1]+high[:-1]
def minkowski(p,q):return hull([(a[0]+b[0],a[1]+b[1]) for a in p for b in q])
def polytope(groups):
 p=[(0,0)]
 for states,n in sorted(groups.items()):p=minkowski(p,[(n*x,n*y) for x,y in states])
 return p
def read_mode(folder,reverse,refs):
 values={};stats=Counter();abundances=Counter();choices=defaultdict(set)
 with zipfile.ZipFile(folder/'Alleles_frequency_table.zip') as z:
  assert z.testzip() is None;assert len(z.namelist())==1
  table=csv.DictReader(io.StringIO(z.read(z.namelist()[0]).decode()),delimiter='\t')
  for row in table:
   seq=row['Aligned_Sequence'].replace('-','');key=min(seq,rc(seq));n=int(row['#Reads']);rn=row['Reference_Name'];code=None
   if rn in ['observed_mutant','programmed_reference','Prime-edited','Scaffold-incorporated']:
    idx=0 if rn=='observed_mutant' else 1;aref,aseq=row['Reference_Sequence'],row['Aligned_Sequence']
    if reverse:aref,aseq=rc(aref),rc(aseq)
    code,_=raw.project_coding(aref,aseq,refs[idx])
   choices[key].add((code,rn));abundances[key]+=n
 for key,options in choices.items():
  code,rn=next(iter(options)) if len(options)==1 else (None,'ambiguous');values[key]=(sys.intern(code) if code is not None else None,rn);stats['resolved_weight' if code is not None else 'unresolved_weight']+=abundances[key]
 return values,dict(stats),abundances
def savecsv(p,rows):
 with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
 p=argparse.ArgumentParser();p.add_argument('--work',type=Path,required=True);p.add_argument('--PE-runs',type=Path,required=True);p.add_argument('--focus',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=a.output.resolve();assert not out.exists();out.mkdir(parents=True);assert json.loads((a.work/'RUN_STATUS.json').read_text())['state']=='COMPLETED'
 refs=[raw.load(ROOT/'raw_replay/metadata'/n) for n in ['REFERENCE.json','REFERENCE_PROGRAMMED.json']];exact=refs[0]['exact_programmed_core']+'AT'
 families={'F3':('T','M_fw','P_fw'),'F4':('T','M_fw','P_fw','PE_default'),'F9':('T','M_fw','P_fw','M_rc','P_rc','MP_fw','MP_rc','PE_default','PE_zero')};events={'native':lambda s:s[0],'exact_programmed':lambda s:s[1],'native_nonexact':lambda s:int(s==(1,0))}
 focus={r['read_id'] for r in csv.DictReader(a.focus.open())};focus_rows=[];counts=defaultdict(Counter);groups=defaultdict(Counter);meta={};N=Counter();diagnostics={};effects=Counter();profile_counts=defaultdict(Counter)
 for rid,records in itertools.groupby(cor.records(),key=lambda x:x[0]['run_accession']):
  maps={};expected={};seen=defaultdict(Counter);diagnostics[rid]={}
  for mode in ['M_rc','P_rc','MP_fw','MP_rc','PE_zero','PE_default']:
   folder=(a.PE_runs/('CRISPResso_on_'+rid)) if mode=='PE_default' else (a.work/'runs'/mode/('CRISPResso_on_'+rid));maps[mode],diagnostics[rid][mode],expected[mode]=read_mode(folder,not mode.endswith('_fw'),refs)
  for r,name,seq,q,old in records:
   N[rid]+=1;meta[rid]=r;codes={'T':old['target'],'M_fw':old['mutant_projection'],'P_fw':old['programmed_projection']};assignment={}
   for mode,mp in maps.items():
    key=min(seq,rc(seq));value=mp.get(key) if q else None
    if value:codes[mode]=value[0];assignment[mode]=value[1];seen[mode][key]+=1
    else:codes[mode]=None;assignment[mode]='quality_U'
   for mode,c in codes.items():
    for event,fun in events.items():profile_counts[rid,mode,event]['U' if c is None else 'E' if fun(signature(c,exact)) else 'NE']+=1
   for mode in ['MP_fw','MP_rc']:
    source=('M_' if assignment[mode]=='observed_mutant' else 'P_')+mode.split('_')[1]
    if assignment[mode] in ['observed_mutant','programmed_reference'] and codes[mode]!=codes[source]:effects[mode+'_versus_same_selected_single_reference']+=1
   for family,ops in families.items():
    cs=[codes[k] for k in ops];complete=all(c is not None for c in cs);states={signature(c,exact) for c in cs} if complete else STATES;key=tuple(sorted(states));groups[rid,family][key]+=1
    for cert,allowed in [('target_relative',states),('string_witness',states if complete and len(set(cs))==1 else STATES)]:
     for event,fun in events.items():vals={fun(s) for s in allowed};label='U' if len(vals)>1 else 'E' if vals=={1} else 'NE';counts[rid,family,cert,event][label]+=1
   if name in focus:
    fw=old['target'];pe=codes['PE_default'];row={'run':rid,'read_id':name,**{k:'' if c is None else c for k,c in codes.items()},'orientation_reproduces_PE_suffix':codes['M_rc']==pe or codes['P_rc']==pe,'PE_zero_equals_PE_default':codes['PE_zero']==pe,'PE_zero_equals_MP_rc':codes['PE_zero']==codes['MP_rc'],'MP_fw_equals_old_target':codes['MP_fw']==fw};focus_rows.append(row)
  for mode in maps:assert seen[mode]==expected[mode],('Raw-to-allele abundance mismatch',rid,mode)
 rows=[]
 for (rid,fam,cert,event),c in sorted(counts.items()):
  n=N[rid];assert sum(c.values())==n;rows.append({'run':rid,'preparation':meta[rid]['source_replica'],'condition':meta[rid]['condition'],'family':fam,'certificate':cert,'event':event,'N':n,'E':c['E'],'NE':c['NE'],'U':c['U'],'L_exact':str(F(c['E'],n)),'H_exact':str(F(c['E']+c['U'],n))})
 contrasts=[]
 for t in rows:
  if t['condition']!='PE treated':continue
  for control in ['mock','untreated']:
   c=next(x for x in rows if x['preparation']==t['preparation'] and x['condition']==control and all(x[k]==t[k] for k in ['family','certificate','event']));lo=F(t['L_exact'])-F(c['H_exact']);hi=F(t['H_exact'])-F(c['L_exact']);contrasts.append({k:t[k] for k in ['preparation','family','certificate','event']}|{'control':control,'L_exact':str(lo),'H_exact':str(hi),'L_pp':float(100*lo),'H_pp':float(100*hi)})
 polygons=[]
 for (rid,fam),g in sorted(groups.items()):
  vertices=polytope(g);n=N[rid]
  # Independent support-function checks over every candidate signature group.
  for dx,dy in [(1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1),(1,1),(-1,-1)]:assert max(dx*x+dy*y for x,y in vertices)==sum(count*max(dx*x+dy*y for x,y in states) for states,count in g.items())
  polygons.append({'run':rid,'preparation':meta[rid]['source_replica'],'condition':meta[rid]['condition'],'family':fam,'N':n,'integer_count_hull_vertices':vertices,'signature_domains':[{'allowed_signatures':states,'records':count} for states,count in sorted(g.items())],'scope':'Sharp convex hull of joint native/exact record counts under the declared finite extraction-family candidate domains; not a continuous finite identified set or biological truth region.'})
 savecsv(out/'FAMILY_CERTIFICATE_COUNTS.csv',rows);savecsv(out/'FAMILY_CERTIFICATE_CONTRASTS.csv',contrasts);savecsv(out/'FOCUS_89_MECHANISM.csv',focus_rows)
 (out/'JOINT_IDENTIFICATION_POLYGONS.json').write_text(json.dumps({'polygons':polygons},indent=2)+'\n')
 savecsv(out/'PROFILE_EVENT_COUNTS.csv',[{'run':rid,'mode':mode,'event':event,'N':N[rid],**{k:c[k] for k in ['E','NE','U']}} for (rid,mode,event),c in sorted(profile_counts.items())])
 assert len(focus_rows)==89 and sum(N.values())==177006
 def unknown(fam,cert,event):return sum(r['U'] for r in rows if r['family']==fam and r['certificate']==cert and r['event']==event)
 assert unknown('F3','string_witness','native')==6281 and unknown('F4','string_witness','native')==6374
 summary={'status':'PASS','new_engine_jobs':45,'records':177006,'focus_records':89,'focus_orientation_reproduces_PE_suffix':sum(r['orientation_reproduces_PE_suffix'] for r in focus_rows),'focus_PE_zero_equals_PE_default':sum(r['PE_zero_equals_PE_default'] for r in focus_rows),'focus_PE_zero_equals_MP_rc':sum(r['PE_zero_equals_MP_rc'] for r in focus_rows),'focus_MP_fw_equals_old_target':sum(r['MP_fw_equals_old_target'] for r in focus_rows),'same_selected_single_reference_mismatch_counts':dict(effects),'unknown_counts':{f:{c:{e:unknown(f,c,e) for e in events} for c in ['string_witness','target_relative']} for f in families},'joint_polygons':len(polygons),'joint_support_directions_checked':8*len(polygons),'finite_target_domains':'Required-operator U is wildcard over all three logically possible native/exact signatures; not dropped. Resolved finite candidate signatures need not have identical DNA.','not_claimed':'New cohorts, biological accuracy, independent human replication, universal model coverage or measured analyst utility.'}
 (out/'MECHANISM_AND_CERTIFICATE_VERIFICATION.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
