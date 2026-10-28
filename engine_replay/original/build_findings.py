"""Condense verified numerical outcomes without elevating their evidential scope."""
import json, subprocess
from collections import defaultdict
import run_grid as g

s=json.loads((g.BASE/'SUMMARIES.json').read_text());t=json.loads((g.BASE/'TRANSITIONS.json').read_text());c=json.loads((g.BASE/'CONTRASTS.json').read_text());v=json.loads((g.BASE/'VERIFICATION.json').read_text())
joints=json.loads((g.BASE/'JOINT_STABILITY.json').read_text())
assert v['state']=='PASS' and v['stage']=='all' and len(s)==144
subprocess.run(['git','-C',str(g.SOURCE),'diff','--exit-code'],check=True)
changes=[]
for x in t:
    tr=x['transitions'];off=sum(n for k,n in tr.items() if k[0]!=k[2]);conf=tr['E>N']+tr['N>E']
    if off:changes.append({'reference':x['reference'],'setting':x['setting'],'run':x['run'],'event':x['event'],'labels_changed_vs_baseline':off,'E_NE_conflicts_vs_baseline':conf,'transitions':tr})
totals=defaultdict(lambda:{'labels_changed_vs_baseline':0,'E_NE_conflicts_vs_baseline':0})
for x in changes:
    key=(x['reference'],x['setting'],x['event'])
    for k in ['labels_changed_vs_baseline','E_NE_conflicts_vs_baseline']:totals[key][k]+=x[k]
aggregated=[{'reference':k[0],'setting':k[1],'event':k[2],**val} for k,val in totals.items()]
native='native_local_peptide_all_observed_sequences';exact='CTT_plus_full_SE2';primary=[]
for prep in [1,2,3]:
    for e in [exact,native]:
        xs=[x for x in c if x['preparation']==prep and x['control']=='mock' and x['event']==e and x['setting']!='stock_default_no_guide']
        lo=min(xs,key=lambda x:x['L_pp']);hi=max(xs,key=lambda x:x['H_pp'])
        primary.append({'preparation':prep,'event':e,'OAT_reference_settings':len(xs),'L_min_pp':lo['L_pp'],'L_min_at':{k:lo[k] for k in ['reference','setting']},'H_max_pp':hi['H_pp'],'H_max_at':{k:hi[k] for k in ['reference','setting']},'all_lower_endpoints_positive':all(x['L_pp']>0 for x in xs)})
jt=defaultdict(lambda:defaultdict(int))
for x in joints:
    for k,n in x['counts'].items():jt[x['scope'],x['event']][k]+=n
joint_totals=[{'scope':k[0],'event':k[1],'N_raw':sum(z.values()),'counts':dict(z)} for k,z in jt.items()]
stock=[]
for x in s:
    if x['setting']!='stock_default_no_guide':continue
    b=next(z for z in s if z['run']==x['run'] and z['reference']==x['reference'] and z['setting']=='baseline');eff=x['effective_configuration'];idx=eff['include_idxs']
    stock.append({'run':x['run'],'reference':x['reference'],'projected_counts_equal_baseline':x['projected']==b['projected'],'modified_percent_aligned_baseline':float(b['native_software_quantification']['Modified%']),'modified_percent_aligned_stock':float(x['native_software_quantification']['Modified%']),'QW_first_last_count':[min(idx),max(idx),len(idx)],'effective_gap_incentive_nonzero_positions':sum(z!=0 for z in eff['effective_gap_incentive_vector'])})
out={'state':'VERIFIED_FINITE_PARAMETER_SENSITIVITY','jobs':144,'custom_baseline_jobs':18,'one_at_a_time_perturbation_jobs':108,'stock_default_no_guide_jobs':18,'biological_raw_records':177006,'baseline_exact_match_to_frozen':'18/18; native quantification tables and both projected molecular events','third_party_source_modified':False,
  'changes_aggregated_by_profile':aggregated,'changed_run_event_rows':changes,'joint_stability_totals':joint_totals,'primary_PE_mock_extrema':primary,'stock_profile':stock,
  'stock_separate_estimand':'MODIFIED% is the software classification over its effective quantification window, not native CFTR peptide frequency.',
  'limits':['Post hoc grid with fixed choices before new replay results, not a prospective holdout.','Same nine primary biological streams; parameter jobs are not independent replicates.','No molecular truth calibration; E/NE disagreements are computational sensitivity.','Ranges cover only the two declared references and finite grid.','Observed local file truncations were preserved and excluded until successful verified reconstruction/replay.']}
g.dump(g.BASE/'FINDINGS.json',out)
statuses=[]
for p in sorted((g.BASE/'executions').rglob('*.json')):
    x=json.loads(p.read_text());assert x['state']=='COMPLETED';statuses.append({k:x[k] for k in ['id','reference','setting','run','state','returncode','elapsed_seconds']})
assert len(statuses)==144;g.dump(g.BASE/'JOB_STATUS.json',statuses)
print(json.dumps({'jobs':144,'nonzero_transition_profile_rows':len(aggregated),'primary_extrema':primary,'stock_projected_equal_baseline_all':all(x['projected_counts_equal_baseline'] for x in stock)}),flush=True)
