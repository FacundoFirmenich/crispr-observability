"""Analyze authenticated dedicated PE mode without equating software classes to events.

Unique selected reference assignments are projected on the frozen local coding
interval. Multi-reference ties, missing assignments and unresolved flanks stay U.
Original raw denominators and complete event/non-event/unknown partitions remain.
"""
from pathlib import Path
from collections import Counter,defaultdict
from fractions import Fraction as F
import argparse,csv,gzip,hashlib,importlib.util,json,re
ROOT=Path(__file__).resolve().parents[3]
def module(name,p):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
raw=module('raw',ROOT/'raw_replay/run_raw_replay.py');cor=module('cor',ROOT/'extensions/reproduction/analyze_correspondence.py')
def rc(s):return s.translate(raw.COMP)[::-1]
def unwrap(x):return x['value'] if isinstance(x,dict) and 'value' in x else x
def read_annotations(path,refs):
 m={};diagnostic=Counter();classes=Counter()
 with gzip.open(path/'CRISPResso_output.fastq.gz','rt') as f:
  while h:=f.readline():
   s,plus,q=f.readline().strip(),f.readline().rstrip(),f.readline().strip();name=h.split()[0][1:];assert name not in m
   tags=dict(re.findall(r'([A-Z_]+)=([^ ]*)',plus));aln=tags['ALN'].split('&');coding=None;reason='unmapped' if aln==['NA'] else 'ambiguous_reference' if len(aln)!=1 else 'unresolved_projection'
   if len(aln)==1 and aln[0]!='NA':
    rn=aln[0];assert rn in ['observed_mutant','Prime-edited','Scaffold-incorporated']
    aref,aseq=rc(tags['ALN_REF']),rc(tags['ALN_SEQ']);assert aseq.replace('-','') in {s,rc(s)},'Alignment imputes or removes observed sequence'
    ref=refs[0 if rn=='observed_mutant' else 1];coding,_=raw.project_coding(aref,aseq,ref)
    if coding is not None:reason='resolved'
   diagnostic[reason]+=1;classes[tags.get('CLASS','UNMAPPED')]+=1;m[name]={'coding':coding,'reason':reason,'assignment':'&'.join(aln),'class':tags.get('CLASS','UNMAPPED'),'sequence':s}
 return m,diagnostic,classes
def csvout(p,rows):
 with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--runs',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=a.output.resolve();assert not out.exists();out.mkdir(parents=True)
 refs=[raw.load(ROOT/'raw_replay/metadata'/n) for n in ['REFERENCE.json','REFERENCE_PROGRAMMED.json']];exact=refs[0]['exact_programmed_core']+'AT';design=raw.load(ROOT/'engine_replay/dependencies/m2/DESIGN.json');metadata={r['run_accession']:r for r in design['runs']}
 maps={};native=[];sources=[];diagnostics={};class_diagnostics={}
 for rid,r in metadata.items():
  path=a.runs/('CRISPResso_on_'+rid);info=raw.load(path/'CRISPResso2_info.json');assert info['running_info']['version']=='2.3.4'
  assert rc(unwrap(info['results']['refs']['observed_mutant']['sequence']))==refs[0]['sequence']
  assert rc(unwrap(info['results']['refs']['Prime-edited']['sequence']))==refs[1]['sequence']
  quant=list(csv.DictReader((path/'CRISPResso_quantification_of_editing_frequency.txt').open(),delimiter='\t'))
  for q in quant:native.append({'run':rid,'preparation':r['source_replica'],'condition':r['condition'],**q})
  maps[rid],diagnostics[rid],class_diagnostics[rid]=read_annotations(path,refs)
  sources.append({'run':rid,'info_sha256':raw.sha(path/'CRISPResso2_info.json'),'annotated_fastq_sha256':raw.sha(path/'CRISPResso_output.fastq.gz'),'quantification_sha256':raw.sha(path/'CRISPResso_quantification_of_editing_frequency.txt')})
 counts=defaultdict(Counter);N=Counter();Q=Counter();joint=Counter();changed=Counter();examples=[];pe_witness_rejections=Counter();disagreement_labels=Counter();prefixes=Counter()
 for r,name,s,quality,old in cor.records():
  rid=r['run_accession'];N[rid]+=1;v=maps[rid].pop(name,None);Q[rid]+=quality
  assert (v is not None)==quality,'Annotated output identifiers must exactly match whole-read Q30 input'
  if v is not None:assert v['sequence'] in {s,rc(s)},'Raw-record sequence identity changed'
  coding=v['coding'] if v else None;old_w=next(iter(old.values())) if all(x is not None for x in old.values()) and len(set(old.values()))==1 else None
  augmented=old_w if old_w is not None and coding==old_w else None
  if old_w is not None and augmented is None:pe_witness_rejections[rid]+=1
  for op,st in [('PE_selected_reference',coding),('old_three_operator_witness',old_w),('PE_corroborated_witness',augmented)]:
   for event in ['native_local','exact_programmed']:
    label='U' if st is None else 'E' if (cor.peptide(st)=='GTIKENIIFGVSYD' if event=='native_local' else st==exact) else 'NE';counts[rid,op,event][label]+=1
  if coding is not None and old_w is not None and coding!=old_w:
   changed[rid]+=1
   old_native=cor.peptide(old_w)=='GTIKENIIFGVSYD';new_native=cor.peptide(coding)=='GTIKENIIFGVSYD';old_exact=old_w==exact;new_exact=coding==exact
   disagreement_labels[old_native,new_native,old_exact,new_exact]+=1
   prefix=old_w[:-len(coding)] if old_w.endswith(coding) else None
   prefixes[prefix]+=1
   examples.append({'run':rid,'read_id':name,'old_witness_coding':old_w,'PE_selected_coding':coding,'native_class':v['class'],'old_native':old_native,'PE_native':new_native,'old_exact':old_exact,'PE_exact':new_exact,'old_length':len(old_w),'PE_length':len(coding),'extra_old_prefix_if_PE_is_suffix':prefix})
  if v:
   local='U' if coding is None else 'native' if cor.peptide(coding)=='GTIKENIIFGVSYD' else 'non_native';ex='U' if coding is None else 'exact' if coding==exact else 'nonexact';joint[rid,v['class'],local,ex]+=1
 assert not any(maps.values())
 rows=[]
 for (rid,op,event),c in sorted(counts.items()):
  n=N[rid];assert sum(c.values())==n;rows.append({'run':rid,'preparation':metadata[rid]['source_replica'],'condition':metadata[rid]['condition'],'operator':op,'event':event,'N':n,'E':c['E'],'NE':c['NE'],'U':c['U'],'L_exact':str(F(c['E'],n)),'H_exact':str(F(c['E']+c['U'],n)),'L_percent':100*c['E']/n,'H_percent':100*(c['E']+c['U'])/n})
 contrasts=[]
 for t in rows:
  if t['condition']!='PE treated':continue
  for control in ['mock','untreated']:
   c=next(x for x in rows if x['preparation']==t['preparation'] and x['condition']==control and x['operator']==t['operator'] and x['event']==t['event']);lo=F(t['L_exact'])-F(c['H_exact']);hi=F(t['H_exact'])-F(c['L_exact']);contrasts.append({'preparation':t['preparation'],'control':control,'operator':t['operator'],'event':t['event'],'L_exact':str(lo),'H_exact':str(hi),'L_pp':float(100*lo),'H_pp':float(100*hi)})
 comparison=[]
 for rid,r in metadata.items():
  qs=[x for x in native if x['run']==rid];aligned=int(qs[0]['Reads_aligned_all_amplicons']);prime=next(x for x in qs if x['Amplicon']=='Prime-edited');scaf=next(x for x in qs if x['Amplicon']=='Scaffold-incorporated');bounds=next(x for x in rows if x['run']==rid and x['operator']=='PE_selected_reference' and x['event']=='native_local');w=next(x for x in rows if x['run']==rid and x['operator']=='old_three_operator_witness' and x['event']=='native_local')
  comparison.append({'run':rid,'preparation':r['source_replica'],'condition':r['condition'],'raw_N':N[rid],'Q30_input':Q[rid],'all_amplicons_aligned':aligned,'Prime_edited_assigned':int(prime['Reads_aligned']),'Prime_edited_assignment_percent_of_all_aligned':100*int(prime['Reads_aligned'])/aligned,'Prime_edited_UNMODIFIED_percent_within_Prime_edited':prime['Unmodified%'],'Scaffold_incorporated_assigned':int(scaf['Reads_aligned']),'multi_reference_ambiguous':diagnostics[rid]['ambiguous_reference'],'unmapped_Q30':diagnostics[rid]['unmapped'],'PE_projected_native_L_percent_raw_N':bounds['L_percent'],'PE_projected_native_H_percent_raw_N':bounds['H_percent'],'old_witness_native_L_percent_raw_N':w['L_percent'],'old_witness_native_H_percent_raw_N':w['H_percent']})
 csvout(out/'PE_NATIVE_TOOL_METRICS.csv',native);csvout(out/'PE_EVENT_BOUNDS.csv',rows);csvout(out/'PE_EVENT_CONTRASTS.csv',contrasts);csvout(out/'PE_TOOL_VS_EVENT_COMPARISON.csv',comparison);csvout(out/'PE_CLASS_VS_LOCAL_EVENTS.csv',[dict(run=rid,native_class=cl,projected_native=nl,projected_exact=ex,records=n) for (rid,cl,nl,ex),n in sorted(joint.items())])
 report={'status':'PASS','raw_records':sum(N.values()),'Q30_records':sum(Q.values()),'successful_dedicated_PE_jobs':len(metadata),'engine_version':'2.3.4','authenticates_reference_orientation_and_automatic_programmed_reference':True,'adapter':'Unique selected reference only; RC aligned strings to frozen forward coordinates; immediate observed flanks must be ACGT; no sequence imputation; multi-reference ties remain U.','PE_projection_different_from_old_resolved_witness':dict(changed),'old_witness_records_not_PE_corroborated':dict(pe_witness_rejections),'different_coding_examples':examples,'disagreement_label_transitions':[{'old_native':a,'PE_native':b,'old_exact':c,'PE_exact':d,'records':n} for (a,b,c,d),n in sorted(disagreement_labels.items())],'old_prefix_signatures':[{'prefix':k,'records':v} for k,v in sorted(prefixes.items(),key=lambda kv:str(kv[0]))],'annotation_diagnostics':{k:dict(v) for k,v in diagnostics.items()},'source_outputs':sources,'all_three_PE_vs_mock_native_lower_bounds_positive':all(x['L_pp']>0 for x in contrasts if x['operator']=='PE_selected_reference' and x['event']=='native_local' and x['control']=='mock'),'limits':'Tool native assignments and window-relative modification rates have their stated tool denominators. They are not native local peptide or exact-DNA rates. Projected bounds assume correct accepted labels. No comparative accuracy, independent cohort, wet-lab truth, full allele or function demonstrated.'}
 if examples:csvout(out/'PE_STRING_DISAGREEMENTS.csv',examples)
 (out/'PE_ANALYSIS_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ['source_outputs','different_coding_examples','annotation_diagnostics']},indent=2))
if __name__=='__main__':main()
